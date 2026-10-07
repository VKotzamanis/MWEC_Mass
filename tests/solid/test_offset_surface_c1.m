function test_offset_surface_c1()
%TEST_OFFSET_SURFACE_C1  F2 offset_surface and F2b void_closing_distance on the C1 hull (contract section 3, C1 oracles).
%   Independent references: the contract's C1 oracles (neck d_close 0.10 m; smallest convex
%   principal radius 0.100 m, so no fold at 25.4 or 76.2 mm and the offset folds and closes from
%   d = 0.100 m on), the Python evaluator tests/reference/c1_reference.py (offset half-widths
%   from the deck text) and the v1.0 erosion result for module 4 (void 3.548 m3, at t = 76.2 mm).
%   t_local on the dense grid is measured by this test's own projection on the outer patches.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
deck = fullfile(root, 'Input', 'C1.ms2');
evalc('model = mwecmass.geometry.MS2Parser.parse(deck);');
geo = mwecmass.solid.outer_nurbs(model);
zr = geo.z_range;
prim = [1 2];

% oracle: smallest convex principal radius of the written outer surface
[rmin, zat] = convex_radius(geo.outer(prim));
fprintf('C1: smallest convex principal radius on a 41-point-per-span grid %.15f m at z = %.4f (oracle 0.100 m)\n', rmin, zat);

% F2b: the neck closes at d = 0.10 m; the precast range ends at the v1.0 wall bottom
dc = mwecmass.solid.void_closing_distance(model, [], geo, zr);
% two evaluations of the chord ends on the flat neck sides and one norm: within 8 ulp of 0.1
check(abs(dc - 0.1) <= 8 * eps(0.1), 'C1: d_close %.17g, oracle 0.10', dc);
z_wall = -0.7;
dcp = mwecmass.solid.void_closing_distance(model, [], geo, [zr(1) z_wall]);
fprintf('C1: d_close %.17g m over the hull, %.17g m over [%g, %g]\n', dc, dcp, zr(1), z_wall);

% modular precast: t = t_min = 76.2 mm over the full hull and over the precast void range
t_min = 0.0762;
tic;
[full, rep] = mwecmass.solid.offset_surface(model, [], geo, t_min, zr, struct('t_min', t_min));
fprintf('C1 precast full range: %.1f s\n', toc);
check(numel(full.patches) == numel(geo.outer), 'C1: one inner patch per outer patch at 76.2 mm (no fold piece)');
check(isequal({full.patches.name}, strcat({geo.outer.name}, '_inner1')), 'C1: inner patch names');
dist = @(X) patch_distance(geo.outer, X);
inner_set_checks(geo, full, t_min, dist, 'C1 t = 0.0762', struct('dense', prim));
check(all(full.patches(1).pole) && isequal(full.patches(1).c0_u, geo.outer(1).c0_u), ...
    'C1: inner surface1 keeps the poles and the c0_u row of the outer');
check(rep.ok, 'C1: report');
tic;
part = mwecmass.solid.offset_surface(model, [], geo, t_min, [zr(1) z_wall], struct('t_min', t_min));
fprintf('C1 precast range [%g, %g]: %.1f s, inner up to z = %.17g\n', zr(1), z_wall, toc, part.z_range(2));
check(part.z_range(2) >= z_wall, 'C1: the precast set does not reach the wall bottom');
inner_set_checks(geo, part, t_min, dist, 'C1 precast range', struct('dense', prim, 'z_clip', [zr(1) z_wall]));

% the Python reference: half-widths of the exact offset profile at d
zs = [-3.15 -3.0 -2.7 -2.2 -1.6 -1.2 -1.05 -1.0 -0.8 -0.3 0.5 0.99];
ref = python_offset(root, deck, full.d, zs);
fprintf('C1 t = 0.0762, d = %.17g: inner half-width (F4, flat side) vs the Python exact offset\n', full.d);
dmax = 0;
for i = 1:numel(zs)
    L = slice_lateral(full.patches, zs(i));
    hw = max(L.pts(:, 1));
    fprintf('  z = %6.2f: F2 %.12f m, reference %.12f m, difference %9.2e m\n', zs(i), hw, ref(i), hw - ref(i));
    dmax = max(dmax, abs(hw - ref(i)));
end
fprintf('  largest |difference| %.3g m (eps_fit/2 = %.3g m)\n', dmax, full.eps_fit / 2);

% module 4 of v1.0 (body [-1.3375, -0.7], above the v1.0 ballast level): void volume
edges = [-1.3375 -0.7];
[V8, V12] = void_volume(part.patches, edges);
fprintf('C1 module 4 void [%g, %g] at d = %.6f: %.9f m3 (Gauss 8 per row interval), %.9f m3 (12), v1.0 erosion at t = 0.0762: 3.548 m3\n', ...
    edges, part.d, V8, V12);

% the offset nodes agree with MS2Parser points moved by d along T1 surface_normals
cache = mwecmass.geometry.precompute_boundary_cache(model, 100);
uu = [0.05; 0.2; 0.4; 0.6; 0.7; 0.8; 0.9; 0.95];
vv = [0.1; 0.5; 0.9];
[UU, VV] = ndgrid(uu, vv);
gap = 0;
for p = prim
    name = model.visible_surfs{p};
    Pp = zeros(numel(UU), 3);
    for i = 1:numel(UU)
        Pp(i, :) = model.eval_surface(name, UU(i), VV(i));
    end
    n = mwecmass.solid.surface_normals(model, cache, name, UU(:), VV(:));
    X = Pp - full.d * n;
    gap = max(gap, max(patch_distance(full.patches, X)));
end
fprintf('C1: MS2Parser points moved by d along -surface_normals lie within %.3g m of the inner set (eps_fit/2 = %.3g m)\n', ...
    gap, full.eps_fit / 2);

% thin shell: t = t_min = 25.4 mm, slender neck
t_min = 0.0254;
tic;
thin = mwecmass.solid.offset_surface(model, [], geo, t_min, zr, struct('t_min', t_min));
fprintf('C1 thin shell: %.1f s\n', toc);
check(numel(thin.patches) == numel(geo.outer), 'C1: one inner patch per outer patch at 25.4 mm (no fold piece)');
inner_set_checks(geo, thin, t_min, dist, 'C1 t = 0.0254', struct('dense', prim));
hw = zeros(0, 1);
for z = linspace(-0.45, 0.95, 15)
    L = slice_lateral(thin.patches, z);
    hw(end + 1) = max(L.pts(:, 1)); %#ok<AGROW>
end
fprintf('C1 neck at 25.4 mm: inner half-width %.9f .. %.9f m on z in [-0.45, 0.95] (0.1 - d = %.9f m)\n', ...
    min(hw), max(hw), 0.1 - thin.d);

% knots_from: the same knots at the same t give the same faces; at nearby t a smooth property
same = mwecmass.solid.offset_surface(model, [], geo, t_min, zr, struct('t_min', t_min, 'knots_from', thin));
check(same.refit && isequal({same.patches.surf}, {thin.patches.surf}), 'C1: knots_from at the same t changes the faces');
fprintf('C1: knots_from at the same t reproduces every face bitwise\n');
dt = 2e-4;
A = zeros(1, 5);
okr = false(1, 5);
for j = 0:4
    s = mwecmass.solid.offset_surface(model, [], geo, t_min + j * dt, zr, struct('t_min', t_min, 'knots_from', thin));
    check(isequal(cellfun(@(c) numel(c.knots{1}), {s.patches.surf}), cellfun(@(c) numel(c.knots{1}), {thin.patches.surf})), ...
        'C1: knots_from changed the knot vectors');
    L = slice_lateral(s.patches, -2.0);
    A(j + 1) = L.area;
    okr(j + 1) = s.report.ok;
end
d1 = diff(A);
fprintf('C1 refit on fixed knots, t = 0.0254 + k 2e-4 (k = 0..4): inner area at z = -2 %s m2; differences %s; second differences %s; M1-M3 %s\n', ...
    mat2str(A, 12), mat2str(d1, 6), mat2str(diff(d1), 3), mat2str(okr));

% the offset folds and the neck closes from d = 0.100 m on: t_max = d_close - eps_fit/2
expect_error(@() mwecmass.solid.offset_surface(model, [], geo, dc - 0.01 * t_min / 2, zr, struct('t_min', t_min)), ...
    'mwecmass:solid:VoidClosed');
fprintf('C1: VoidClosed at t = t_max = %.17g m (d = d_close)\n', dc - 0.01 * t_min / 2);
end

function L = slice_lateral(P, z)
lat = arrayfun(@(e) e.z_range(1) ~= e.z_range(2), P);
L = mwecmass.solid.slice_bspline_surface(P(lat), z);
end

function [V8, V12] = void_volume(P, edges)
% integral of the F4 section area over z, Gauss-Legendre between the heights of the knot rows of
% the inner faces (the area is smooth between them)
rows = zeros(0, 1);
for k = 1:numel(P)
    ku = unique(P(k).surf.knots{1});
    Q = mwecmass.solid.eval_bspline_surface(P(k).surf, ku(:), repmat(P(k).surf.knots{2}(1), numel(ku), 1));
    rows = [rows; Q(:, 3)]; %#ok<AGROW>
end
br = unique([edges(1); rows(rows > edges(1) & rows < edges(2)); edges(2)]);
V8 = 0;
V12 = 0;
for n = [8 12]
    [x, w] = gauss(n);
    V = 0;
    for j = 1:numel(br) - 1
        zz = br(j) + (br(j + 1) - br(j)) * (x + 1) / 2;
        for i = 1:n
            L = slice_lateral(P, zz(i));
            V = V + (br(j + 1) - br(j)) / 2 * w(i) * L.area;
        end
    end
    if n == 8
        V8 = V;
    else
        V12 = V;
    end
end
end

function [rmin, zat] = convex_radius(E)
% principal curvatures from the first and second fundamental forms with the inward normal; a
% positive curvature bends toward the hull's inside (convex)
rmin = Inf;
zat = NaN;
for k = 1:numel(E)
    s = E(k).surf;
    us = grid_in(s.knots{1}, 41);
    vs = grid_in(s.knots{2}, 41);
    [UU, VV] = ndgrid(us, vs);
    [S, Su, Sv, Suu, Suv, Svv] = mwecmass.solid.eval_bspline_surface(s, UU(:), VV(:));
    n = cross(Su, Sv, 2);
    ln = sqrt(sum(n.^2, 2));
    n = n ./ ln;
    if E(k).outward
        n = -n;
    end
    e = sum(Su .* Su, 2); f = sum(Su .* Sv, 2); g = sum(Sv .* Sv, 2);
    l = sum(Suu .* n, 2); m = sum(Suv .* n, 2); nn = sum(Svv .* n, 2);
    det1 = e .* g - f.^2;
    H = (e .* nn - 2 * f .* m + g .* l) ./ (2 * det1);
    K = (l .* nn - m.^2) ./ det1;
    kmax = H + sqrt(max(H.^2 - K, 0));
    ok = ln > sqrt(eps) * sqrt(e .* g);
    kmax(~ok) = -Inf;
    [kk, i] = max(kmax);
    if 1 / kk < rmin
        rmin = 1 / kk;
        zat = S(i, 3);
    end
end
end

function s = grid_in(k, m)
ku = unique(k);
s = zeros(0, 1);
for j = 1:numel(ku) - 1
    s = [s; ku(j) + (1:m - 1)' / m * (ku(j + 1) - ku(j))]; %#ok<AGROW>
end
end

function dist = patch_distance(E, X)
% distance from the points X to the patches E: nearest point of a grid on every knot-span cell,
% then projected Gauss-Newton steps on (u, v) of that patch, clamped to its parameter box
dist = Inf(size(X, 1), 1);
for k = 1:numel(E)
    s = E(k).surf;
    us = [grid_in(s.knots{1}, 8); unique(s.knots{1})'];
    vs = [grid_in(s.knots{2}, 8); unique(s.knots{2})'];
    [UU, VV] = ndgrid(us, vs);
    G = mwecmass.solid.eval_bspline_surface(s, UU(:), VV(:));
    u = zeros(size(X, 1), 1);
    v = u;
    for c = 1:500:size(X, 1)
        i = c:min(c + 499, size(X, 1));
        D = (X(i, 1) - G(:, 1)').^2 + (X(i, 2) - G(:, 2)').^2 + (X(i, 3) - G(:, 3)').^2;
        [~, j] = min(D, [], 2);
        u(i) = UU(j);
        v(i) = VV(j);
    end
    lo = [s.knots{1}(1) s.knots{2}(1)];
    hi = [s.knots{1}(end) s.knots{2}(end)];
    for it = 1:60
        [S, Su, Sv] = mwecmass.solid.eval_bspline_surface(s, u, v);
        r = X - S;
        a = sum(Su .* Su, 2); b = sum(Su .* Sv, 2); c2 = sum(Sv .* Sv, 2);
        g1 = sum(Su .* r, 2); g2 = sum(Sv .* r, 2);
        lam = 1e-14 * (a + c2) + realmin;
        det2 = (a + lam) .* (c2 + lam) - b.^2;
        du = ((c2 + lam) .* g1 - b .* g2) ./ det2;
        dv = ((a + lam) .* g2 - b .* g1) ./ det2;
        un = min(max(u + du, lo(1)), hi(1));
        vn = min(max(v + dv, lo(2)), hi(2));
        if isequal(un, u) && isequal(vn, v)
            break
        end
        u = un;
        v = vn;
    end
    S = mwecmass.solid.eval_bspline_surface(s, u, v);
    dist = min(dist, sqrt(sum((X - S).^2, 2)));
end
end

function hw = python_offset(root, deck, d, zs)
cmd = sprintf('python3 -I "%s" --deck "%s" --offset %.17g%s', fullfile(root, 'tests', 'reference', 'c1_reference.py'), ...
    deck, d, sprintf(' %.17g', zs));
[status, txt] = system(cmd);
if status ~= 0
    error('test_offset_surface_c1:python', 'reference evaluator failed:\n%s', txt);
end
r = jsondecode(txt);
hw = zeros(size(zs));
for i = 1:numel(zs)
    if iscell(r)
        hw(i) = r{i}.x_half_width;
    else
        hw(i) = r(i).x_half_width;
    end
end
end

function [x, w] = gauss(n)
b = (1:n - 1) ./ sqrt(4 * (1:n - 1).^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D));
w = 2 * V(1, i)'.^2;
end

function expect_error(f, id)
try
    f();
catch err
    if ~strcmp(err.identifier, id)
        error('test_offset_surface_c1:fail', 'expected %s, got %s (%s)', id, err.identifier, err.message);
    end
    return
end
error('test_offset_surface_c1:fail', 'expected %s, got no error', id);
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'solid'));
end

function check(cond, varargin)
if ~cond
    error('test_offset_surface_c1:fail', varargin{:});
end
end
