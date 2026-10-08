function test_outer_nurbs_c1()
%TEST_OUTER_NURBS_C1  F1 on Input/C1.ms2: S1/S1b fields, C1 oracles, seams, orientation, and I9
%   (F4 sections of the NURBS against T1 outer_rows on MS2Parser).
%   The C1 oracles of contract section 3 (read from the deck by construction) serve as the
%   independent reference.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
evalc('model = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''Input'', ''C1.ms2''));');
tic;
geo = mwecmass.solid.outer_nurbs(model);
fprintf('outer_nurbs(C1): %.2f s\n', toc);

s1 = {'name', 'source', 'type', 'flips', 'surf', 'outward', 'exact', 'z_of_u', 'u_range', 'z_range', ...
    'offset_kind', 'pole', 'c0_u', 'c0_v', 'seam_u0', 'seam_u1', 'seam_v0', 'seam_v1', 'visible', 'fit', 'swap_uv'};
check(isequal(sort(fieldnames(geo))', sort({'hull_name', 'outer', 'z_range', 'analytic', 'flat'})), 'S1b fields');
check(isequal(sort(fieldnames(geo.outer))', sort(s1)), 'S1 fields');
check(strcmp(geo.hull_name, 'C1') && isempty(geo.analytic) && isempty(geo.flat), 'S1b values');
P = geo.outer;
check(numel(P) == 8 && isequal({P.name}, model.visible_surfs) && isequal([P.visible], 1:8), 'one entry per visible surface');
check(all([P.exact]) && all([P.z_of_u]) && ~any([P.swap_uv]) && all(cellfun(@isempty, {P.fit})), 'all 8 exact');

% C1 oracles
rev = strncmp({P.name}, 'surface1', 8);
for k = 1:8
    p = P(k);
    Z = p.surf.ctrl(:, :, 3);
    check(all(all(Z == Z(:, 1))), '%s: rows have one z', p.name);
    check(p.z_range(1) == 1.1, '%s: top at z = 1.1', p.name);
    % the keel row is the parser's bead BeadBottom, which the parser evaluates to -3.25 exactly
    check(p.z_range(2) == -3.25, '%s: keel at z = -3.25', p.name);
    check(numel(p.c0_u) == 1 && isempty(p.c0_v), '%s: one c0_u row, no c0_v', p.name);
    i = find(p.surf.knots{1} == p.c0_u, 1) - 1;
    check(all(p.surf.ctrl(i, :, 3) == 1), '%s: the c0_u row is the arc-to-curve1 joint z = 1.0', p.name);
    if rev(k)
        check(strcmp(p.offset_kind, 'rev_z') && all(p.pole), '%s: rev_z with poles at keel and neck top', p.name);
    else
        check(strcmp(p.offset_kind, 'ruled_parallel') && ~any(p.pole), '%s: ruled_parallel', p.name);
        check(all(p.surf.ctrl(1, :, 1) == 0) && all(p.surf.ctrl(end, :, 1) == 0), ...
            '%s: ridge row and keel row at x = 0', p.name);
        mx = find(strcmp({P.name}, mirror_name(p.name)));
        check(isequal(p.seam_u0, [mx 4]) && isequal(p.seam_u1, [mx 2]), '%s: ridge and keel shared with its X-mirror', p.name);
    end
end
check(isequal(geo.z_range, [-3.25 1.1]), 'geo.z_range = [-3.25 1.1] bitwise');
fprintf('C1: geo.z_range [%.17g %.17g], c0_u %.17g, weights of the arc row %.17g\n', geo.z_range, P(1).c0_u, ...
    P(1).surf.weights(2, 1));

% seams: mutual, bitwise equal rows (C1: every seam from one array or a mirror plane)
n_seam = check_seams(P);
fprintf('C1: %d seam boundaries, all mutual and bitwise equal\n', n_seam);

% same geo with the general-path inputs
cache = mwecmass.geometry.precompute_boundary_cache(model, 100);
geo2 = mwecmass.solid.outer_nurbs(model, cache, struct('t_min', 0.0254));
check(isequal(geo, geo2), 'outer_nurbs(model, cache, opts) differs from outer_nurbs(model)');

% I9: F4 sections of the NURBS against T1 outer_rows (points on MS2Parser's surfaces)
zs = [-3.2; -3.0; -2.706; -2.2; -1.619; -1.2; -1.05; -1.0; -0.95; -0.7; -0.531; -0.3; 0.2; 0.556; 0.99; 1.0; 1.01; 1.08];
rows = mwecmass.solid.outer_rows(model, cache, zs);
[~, orient_t1] = mwecmass.solid.surface_normals(model, cache, 1, 0.5, 0.5);
fprintf('I9   z [m]    points  max |NURBS - T1| [m]  max ratio to T1 accuracy\n');
worst_ratio = 0;
for q = 1:numel(zs)
    L = mwecmass.solid.slice_bspline_surface(P, zs(q));
    check(L.simple && numel(L.pieces) == 8, 'C1 loop at z = %g', zs(q));
    r = rows(q);
    [dev, bnd, nrm] = distance_to_loop(P, L, r.pts);
    i = 1:7:r.n;
    n1 = mwecmass.solid.surface_normals(model, cache, r.patch(i), r.u(i), r.v(i), orient_t1);
    check(all(sum(n1 .* nrm(i, :), 2) > 0), 'z = %g: outward flag disagrees with T1 surface_normals', zs(q));
    ratio = max(dev ./ bnd);
    worst_ratio = max(worst_ratio, ratio);
    fprintf('I9 %7.3f  %6d  %.3e             %.3f\n', zs(q), r.n, max(dev), ratio);
end
check(worst_ratio <= 1, 'I9: NURBS sections differ from T1 outer_rows beyond T1''s in-plane accuracy');
end

function [dmin, bound, nrm] = distance_to_loop(P, L, X)
% distance from each point X(i, :) to the loop's exact row curves (Newton in v from the nearest of
% 201 samples), the T1 accuracy bound there and the outward unit normal of the NURBS. T1 places its
% points within 16 ulp of z_k on the parser surface (or at a parameter fixed to 2 ulp), so their
% in-plane error is at most min(16 ulp |n_z|/|n_h|, sqrt(2 * 16 ulp / |z_ss|)) (first- and
% second-order root location, z_ss the second derivative of z across the row), plus the rounding of
% both evaluations (64 ulp of the coordinates).
m = size(X, 1);
dmin = Inf(m, 1);
kbest = zeros(m, 1);
vbest = zeros(m, 1);
for k = 1:numel(L.pieces)
    c = L.pieces(k).curve;
    s = linspace(c.knots(1), c.knots(end), 201)';
    Q = mwecmass.solid.eval_bspline_curve(c, s);
    D = (X(:, 1) - Q(:, 1)').^2 + (X(:, 2) - Q(:, 2)').^2;
    [~, j] = min(D, [], 2);
    x = s(j);
    for it = 1:40
        [C, C1, C2] = mwecmass.solid.eval_bspline_curve(c, x);
        dC = C(:, 1:2) - X(:, 1:2);
        g = sum(dC .* C1(:, 1:2), 2);
        H = sum(C1(:, 1:2).^2, 2) + sum(dC .* C2(:, 1:2), 2);
        x = min(max(x - g ./ H, c.knots(1)), c.knots(end));
    end
    C = mwecmass.solid.eval_bspline_curve(c, x);
    d = sqrt(sum((C(:, 1:2) - X(:, 1:2)).^2, 2));
    better = d < dmin;
    dmin(better) = d(better);
    kbest(better) = k;
    vbest(better) = x(better);
end
nrm = zeros(m, 3);
bound = zeros(m, 1);
for k = unique(kbest)'
    i = kbest == k;
    pc = L.pieces(k);
    p = P(pc.patch);
    [~, Su, Sv, Suu, Suv, Svv] = mwecmass.solid.eval_bspline_surface(p.surf, repmat(pc.u, nnz(i), 1), vbest(i));
    n = cross(Su, Sv, 2);
    n = n ./ sqrt(sum(n.^2, 2));
    if ~p.outward
        n = -n;
    end
    nrm(i, :) = n;
    nh = sqrt(sum(n(:, 1:2).^2, 2));
    ulp16 = 16 * eps(max(1, abs(X(i, 3))));
    b1 = ulp16 .* abs(n(:, 3)) ./ max(nh, realmin);
    % direction across the row on the surface: w = Su - (Su.Sv / Sv.Sv) Sv, i.e. (du, dv) = (1, -F/G)
    E = sum(Su .* Su, 2);
    F = sum(Su .* Sv, 2);
    G = sum(Sv .* Sv, 2);
    dv = -F ./ G;
    w2 = E + 2 * F .* dv + G .* dv.^2;
    zss = (Suu(:, 3) + 2 * Suv(:, 3) .* dv + Svv(:, 3) .* dv.^2) ./ w2;
    b2 = sqrt(2 * ulp16 ./ max(abs(zss), realmin));
    bound(i) = min(b1, b2) + 64 * eps * max(1, max(abs(X(i, :)), [], 2));
end
end

function n_seam = check_seams(P)
fields = {'seam_v0', 'seam_u1', 'seam_v1', 'seam_u0'};
n_seam = 0;
for k = 1:numel(P)
    for b = 1:4
        nb = P(k).(fields{b});
        if isempty(nb)
            check(b == 2 && P(k).pole(2) || b == 4 && P(k).pole(1), '%s: boundary %d open and not a pole', P(k).name, b);
            continue
        end
        check(isequal(P(nb(1)).(fields{nb(2)}), [k b]), '%s: seam %d not mutual', P(k).name, b);
        [ra, wa, ka] = boundary(P(k).surf, b);
        [rb, wb, kb] = boundary(P(nb(1)).surf, nb(2));
        same = isequal(ra, rb) && isequal(wa, wb);
        rev = isequal(ra, flipud(rb)) && isequal(wa, flipud(wb));
        check((same || rev) && isequal(ka, kb), '%s: seam %d rows or knots differ', P(k).name, b);
        n_seam = n_seam + 1;
    end
end
end

function [r, w, kn] = boundary(s, b)
W = s.weights;
if isempty(W)
    W = ones(size(s.ctrl, 1), size(s.ctrl, 2));
end
switch b
    case 1, r = squeeze(s.ctrl(:, 1, :)); w = W(:, 1); kn = s.knots{1};
    case 2, r = squeeze(s.ctrl(end, :, :)); w = W(end, :)'; kn = s.knots{2};
    case 3, r = squeeze(s.ctrl(:, end, :)); w = W(:, end); kn = s.knots{1};
    case 4, r = squeeze(s.ctrl(1, :, :)); w = W(1, :)'; kn = s.knots{2};
end
end

function m = mirror_name(name)
if ~isempty(strfind(name, '_mirrX')) %#ok<STREMP>
    m = strrep(name, '_mirrX', '');
elseif ~isempty(strfind(name, '_mirrY')) %#ok<STREMP>
    m = strrep(name, '_mirrY', '_mirrX_mirrY');
else
    m = [name '_mirrX'];
end
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
end

function check(cond, varargin)
if ~cond
    error('test_outer_nurbs_c1:fail', varargin{:});
end
end
