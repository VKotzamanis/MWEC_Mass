function st = body_checks(geo, body, rho, name, opts)
%BODY_CHECKS  Shared assertions of the F5, F6 tests on one body (helper, not a test).
%
%   st = body_checks(geo, body, rho, name, opts)
%
%   geo: the S1b the body was built from; body: S4 (build_body); rho: densities by region.
%   opts.step (default true): write the bodies with write_step and import them with step_check.py;
%   opts.mesh_size: element size passed to step_check.py for its mesh volumes (default its own);
%   opts.n_gauss: F6 order (default F6's).
%   Asserted:
%   - S4 fields; validate_brep on body.brep and on the closed shells of body.shells (precast: the
%     fused solid all_outer with the cavities all_voids; thin shell: layer and void as solids).
%   - I2: every edge in the loops of a B-spline face is the same curve (contract S1) as one boundary
%     row of that face, control points bitwise; each row that is not collapsed is exactly one edge
%     of the face's loops and a collapsed row none (a loop of 3, or 2, edges).
%   - I7: the outer faces of every module enclose a positive volume, and no region volume is
%     negative beyond the rounding bound below.
%   - I1 (rule 11): per module the region volumes of F6 sum to the module volume enclosed by the
%     module's outer faces alone; the module volumes sum to the hull volume integrated on the
%     unsplit outer patches of geo with the cut parameters u* as extra breakpoints.
%   - STEP (opts.step): METRE, closed volumes, the solid count of the mode.
%   st: bp (F6), V_module_outer [N x 1], V_hull, err_module [N x 1], bound_module, err_hull,
%   bound_hull, abs_moments [1 x 10] (sum over the lateral faces of the absolute F6 face integrals
%   [V, int x, int y, int z, int x^2, int y^2, int z^2, int xy, int xz, int yz]), n_faces,
%   area (of the lateral faces), n_nodes (largest number of Gauss nodes of one outer patch), occ,
%   mesh, ref (F6 volumes of the imported solids, same order).

if nargin < 5
    opts = struct();
end
if ~isfield(opts, 'step')
    opts.step = true;
end
fopts = struct();
if isfield(opts, 'n_gauss')
    fopts.n_gauss = opts.n_gauss;
end
check(isequal(sort(fieldnames(body))', sort({'design', 'inner_t', 'planes', 'brep', 'shells', 'voids'})), ...
    '%s: S4 fields', name);
B = body.brep;
mwecmass.output.step.validate_brep(B);
precast = strcmp(body.design.mode, 'modular_precast');
if precast
    names = {'uhpc', 'air'};
    mass_brep = B;
    mass_brep.bodies = struct('name', [geo.hull_name '_UHPC_all'], 'kind', 'solid', ...
        'shells', {[{body.shells.all_outer}, body.shells.all_voids]});
else
    names = {'ballast', 'shell', 'air'};
    mass_brep = B;
    mass_brep.bodies = struct('name', {'layer', 'void'}, 'kind', 'solid', ...
        'shells', {{body.shells.layer}, {body.shells.void}});
    mass_brep.bodies = mass_brep.bodies(~cellfun(@isempty, {body.shells.layer, body.shells.void}));
end
mwecmass.output.step.validate_brep(mass_brep);
shared_edges(B, name);

bp = mwecmass.solid.body_properties(body, rho, fopts);
ng = bp.quad;
N = numel(body.design.edges) - 1;

% I1 per module: the outer faces alone (same quadrature code and nodes as F6, so the module volume
% and the region volumes differ only by the order of summation of the same face integrals)
Vout = zeros(N, 1);
absf = zeros(N, 1);
nf = zeros(N, 1);
absmom = zeros(1, 10);
area = 0;
for k = 1:numel(B.faces)
    f = B.faces(k);
    s = B.surfaces{f.surface};
    if ~strcmp(s.type, 'bspline') || all(reshape(s.ctrl(:, :, 3), [], 1) == s.ctrl(1, 1, 3))
        continue
    end
    [m, ~, ~, mom, a] = face_volume(s, ng);
    if ~f.same_sense
        m = -m;
    end
    absmom = absmom + abs(mom);
    area = area + a;
    i = f.module(1);
    absf(i) = absf(i) + abs(m);
    nf(i) = nf(i) + 1;
    if strcmp(f.role, 'outer')
        if strcmp(f.outside, 'exterior')
            Vout(i) = Vout(i) + m;
        else
            Vout(i) = Vout(i) - m;
        end
    end
end
err_m = zeros(N, 1);
% recursive summation of nf face integrals into at most three region sums and their sum: each
% partial sum is rounded once per addition, so the two totals differ by at most (nf + 4) eps times
% the sum of the absolute face integrals
bound_m = (nf + 4) * eps .* absf;
for i = 1:N
    Vr = cellfun(@(r) bp.modules(i).(['V_' r]), names);
    for r = 1:numel(names)
        check(Vr(r) >= -bound_m(i), '%s: module %d region %s volume %.3e < 0', name, i, names{r}, Vr(r));
    end
    check(Vout(i) > 0, '%s: module %d outer faces enclose %.3e m3 (outer normals not out of the hull)', name, i, Vout(i));
    err_m(i) = abs(sum(Vr) - Vout(i));
    check(err_m(i) <= bound_m(i), '%s: module %d volume closure %.3e > %.3e', name, i, err_m(i), bound_m(i));
    check(abs(bp.modules(i).V - sum(Vr)) <= bound_m(i), '%s: module %d V differs from its regions', name, i);
end

% I1 hull: unsplit outer patches, breakpoints at their knots and at the cuts u* of every plane
planes = body.planes;
Vh = 0;
T = 0;
nmax = 0;
for p = 1:numel(geo.outer)
    P = geo.outer(p);
    if P.z_range(1) == P.z_range(2)
        continue
    end
    br = unique(P.surf.knots{1});
    for z = planes(planes > min(P.z_range) & planes < max(P.z_range))
        [~, ~, us] = mwecmass.solid.split_bspline_surface(P, z);
        br = unique([br(:); us]);
    end
    [m, scale, nn] = face_volume(P.surf, ng, br);
    if ~P.outward
        m = -m;
    end
    Vh = Vh + m;
    T = T + scale;
    nmax = max(nmax, nn);
end
% per node the split pieces and the unsplit patch are evaluated from control points that differ by
% the rounding of knot insertion; with the evaluation, the cross product and the product with x
% this is at most 64 ulp of |x| |S_u| |S_v| per node (each step a few convex combinations or one
% product), and the face sums add the recursive summation bound of their node count
bound_h = (64 + nmax + numel(B.faces)) * eps * T;
err_h = abs(sum(Vout) - Vh);
check(err_h <= bound_h, '%s: sum of module volumes %.17g differs from the hull %.17g by %.3e > %.3e', ...
    name, sum(Vout), Vh, err_h, bound_h);
check(abs(sum([bp.modules.V]) - sum(Vout)) <= sum(bound_m), '%s: sum of F6 module volumes', name);

fprintf('%s: %d faces, %d edges; hull %.12f m3; closure: max per module %.2e (bound %.2e), hull %.2e (bound %.2e)\n', ...
    name, numel(B.faces), numel(B.edges), Vh, max(err_m), max(bound_m), err_h, bound_h);
for i = 1:N
    md = bp.modules(i);
    fprintf('  module %d [%8.4f, %8.4f]: V %.9f', i, body.design.edges(i), body.design.edges(i + 1), md.V);
    for r = 1:numel(names)
        fprintf(', V_%s %.9f', names{r}, md.(['V_' names{r}]));
    end
    fprintf(', mass %.4f kg, CG_body [%.2e %.2e %.9f]\n', md.mass, md.CG_body);
end
fprintf('  total mass %.6f kg, CG_body [%.3e %.3e %.12f] m, Ixx %.6f, Iyy %.6f, Izz %.6f kg m2 (about the CG)\n', ...
    bp.total.mass, bp.total.CG_body, diag(bp.total.I_cg));

st = struct('bp', bp, 'V_module_outer', Vout, 'V_hull', Vh, 'err_module', err_m, 'bound_module', bound_m, ...
    'err_hull', err_h, 'bound_hull', bound_h, 'abs_moments', absmom, 'n_faces', sum(nf), 'area', area, 'n_nodes', nmax, 'occ', [], 'mesh', [], 'ref', []);
if ~opts.step
    return
end
args = {};
if isfield(opts, 'mesh_size')
    args = {'--mesh-size', sprintf('%.17g', opts.mesh_size)};
end
file = [tempname() '.step'];
cleanup = onCleanup(@() delete_if(file));
mwecmass.output.step.write_step(B, file);
res = stp_check(file, args{:});
check(strcmp(res.declared_length_unit, 'METRE'), '%s: length unit %s', name, res.declared_length_unit);
check(res.open_edges_volumes == 0, '%s: %d open mesh edges on the solids', name, res.open_edges_volumes);
if precast
    check(res.n_volumes == N, '%s: %d solids imported, expected %d', name, res.n_volumes, N);
    ref = arrayfun(@(m) m.V_uhpc, bp.modules(:)');
else
    nb = any(strcmp({B.bodies.name}, [geo.hull_name '_STEEL_ballast']));
    check(res.n_volumes == nb, '%s: %d solids imported, expected %d', name, res.n_volumes, nb);
    check(strcmp(B.bodies(end).name, [geo.hull_name '_STEEL_shell']) && strcmp(B.bodies(end).kind, 'sheet'), ...
        '%s: the shell sheet is missing', name);
    ref = sum(bp.regions.ballast.V) * ones(1, nb);
end
occ = res.occ_volumes(:)';
mesh = res.mesh_volumes(:)';
mwecmass.output.step.write_step(mass_brep, file);
r2 = stp_check(file, args{:});
check(r2.open_edges_volumes == 0 && r2.n_volumes == numel(mass_brep.bodies), '%s: %d closed solids from the fused shells', ...
    name, r2.n_volumes);
if precast
    ref = [ref, sum(ref)];
else
    rv = [sum(bp.regions.shell.V), sum(bp.regions.air.V)];
    ref = [ref, rv(~cellfun(@isempty, {body.shells.layer, body.shells.void}))];
end
occ = [occ, r2.occ_volumes(:)'];
mesh = [mesh, r2.mesh_volumes(:)'];
fprintf('  STEP import (%d + %d solids, METRE, closed): F6 %s\n', res.n_volumes, r2.n_volumes, sprintf('%.9f ', ref));
fprintf('    OCC getMass - F6 (relative) %s\n', sprintf('%9.2e ', (occ - ref) ./ ref));
fprintf('    mesh volume - F6 (relative) %s\n', sprintf('%9.2e ', (mesh - ref) ./ ref));
st.occ = occ;
st.mesh = mesh;
st.ref = ref;
end

function [m, scale, nn, mom, area] = face_volume(s, ng, br)
% volume integral of x n_x over the face with F6's Gauss rule (spans between br, default the
% distinct u knots); scale = sum over the nodes of |w| |x| |S_u| |S_v|; mom: the ten F6 integrals
if nargin < 3
    br = unique(s.knots{1});
end
[xg, wg] = gauss_legendre(ng);
kv = unique(s.knots{2});
au = br(1:end - 1);
bu = br(2:end);
av = kv(1:end - 1);
bv = kv(2:end);
uu = (au(:)' + bu(:)') / 2 + (bu(:)' - au(:)') / 2 .* xg;
vv = (av(:)' + bv(:)') / 2 + (bv(:)' - av(:)') / 2 .* xg;
wu = (bu(:)' - au(:)') / 2 .* wg;
wv = (bv(:)' - av(:)') / 2 .* wg;
[U, V] = ndgrid(uu(:), vv(:));
[Wu, Wv] = ndgrid(wu(:), wv(:));
[P, Su, Sv] = mwecmass.solid.eval_bspline_surface(s, U(:), V(:));
nx = Su(:, 2) .* Sv(:, 3) - Su(:, 3) .* Sv(:, 2);
w = Wu(:) .* Wv(:) .* nx;
m = w' * P(:, 1);
scale = (abs(Wu(:) .* Wv(:)) .* abs(P(:, 1)))' * (sqrt(sum(Su.^2, 2)) .* sqrt(sum(Sv.^2, 2)));
nn = numel(U);
x = P(:, 1);
y = P(:, 2);
z = P(:, 3);
area = (Wu(:) .* Wv(:))' * sqrt(sum(cross(Su, Sv, 2).^2, 2));
mom = w' * [x, x.^2 / 2, x .* y, x .* z, x.^3 / 3, x .* y.^2, x .* z.^2, x.^2 .* y / 2, x.^2 .* z / 2, x .* y .* z];
end

function [x, w] = gauss_legendre(n)
b = (1:n - 1) ./ sqrt(4 * (1:n - 1).^2 - 1);
[Vv, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D));
w = 2 * Vv(1, i)'.^2;
end

function shared_edges(B, name)
% I2 on every B-spline face: loop edges are the face's boundary rows, one edge per row that is not
% collapsed
for k = 1:numel(B.faces)
    s = B.surfaces{B.faces(k).surface};
    if ~strcmp(s.type, 'bspline')
        continue
    end
    [nu, nv, ~] = size(s.ctrl);
    W = s.weights;
    rows = {row_curve(s.ctrl(:, 1, :), s.knots{1}, wsel(W, ':', 1)), ...
            row_curve(s.ctrl(nu, :, :), s.knots{2}, wsel(W, nu, ':')), ...
            row_curve(s.ctrl(:, nv, :), s.knots{1}, wsel(W, ':', nv)), ...
            row_curve(s.ctrl(1, :, :), s.knots{2}, wsel(W, 1, ':'))};
    live = cellfun(@(r) ~all(all(r.ctrl == r.ctrl(1, :))), rows);
    E = abs([B.faces(k).loops{:}]);
    check(numel(E) == sum(live), '%s: face %d has %d loop edges for %d rows that are not collapsed', name, k, numel(E), sum(live));
    hit = zeros(1, 4);
    for e = E
        c = B.curves(B.edges(e).curve);
        found = 0;
        for b = find(live)
            if same_curve(c, rows{b})
                found = b;
                break
            end
        end
        check(found > 0, '%s: face %d edge %d is not the same curve as a boundary row of the face', name, k, e);
        hit(found) = hit(found) + 1;
    end
    check(isequal(hit(live), ones(1, sum(live))), '%s: face %d rows not matched one to one by its edges', name, k);
end
end

function r = row_curve(C, kn, w)
r = struct('ctrl', reshape(C, [], 3), 'knots', kn(:)', 'weights', w(:));
end

function w = wsel(W, i, j)
if isempty(W)
    w = [];
else
    w = W(i, j);
end
end

function tf = same_curve(a, b)
% contract S1: control points and weights bitwise (empty = ones), knots equal after the map to
% [0, 1] within 4 * 2^-52 (three correctly rounded operations on values in [0, 1]), either direction
n = size(a.ctrl, 1);
tf = false;
if size(b.ctrl, 1) ~= n || numel(a.knots) ~= numel(b.knots)
    return
end
wa = a.weights(:);
wb = b.weights(:);
if isempty(wa), wa = ones(n, 1); end
if isempty(wb), wb = ones(n, 1); end
ka = (a.knots(:)' - a.knots(1)) / (a.knots(end) - a.knots(1));
kb = (b.knots(:)' - b.knots(1)) / (b.knots(end) - b.knots(1));
kr = (b.knots(end) - b.knots(end:-1:1)) / (b.knots(end) - b.knots(1));
tf = (isequal(a.ctrl, b.ctrl) && isequal(wa, wb) && all(abs(ka - kb) <= 4 * 2^-52)) || ...
     (isequal(a.ctrl, flipud(b.ctrl)) && isequal(wa, flipud(wb)) && all(abs(ka - kr(:)') <= 4 * 2^-52));
end

function delete_if(f)
if exist(f, 'file')
    delete(f);
end
end

function check(cond, varargin)
if ~cond
    error('body_checks:fail', varargin{:});
end
end
