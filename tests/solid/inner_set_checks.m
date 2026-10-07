function st = inner_set_checks(geo, inner, t_min, dist_fn, name, opts)
%INNER_SET_CHECKS  Shared assertions of the F2 tests on one S2 inner set (helper, not a test).
%
%   st = inner_set_checks(geo, inner, t_min, dist_fn, name, opts)
%
%   dist_fn(X) -> [n x 1]: distance of the points X [n x 3] to the outer surface as written, computed
%   by the caller independently of offset_surface (closed form or its own projection).
%   opts.dense: indices of the inner patches checked on the dense grid (default all; mirrors are
%   asserted to be flipped copies of their primary, so a test may check primaries only).
%   opts.z_clip: [z_lo z_hi] of the requested z_range; a boundary row at or beyond an end of it is
%   the open end of a partial set (default: none).
%   Asserted: the S2 and S2r field sets, d = t + eps_fit/2 and eps_fit = 0.01 t_min (bitwise), the
%   fit report (ok, no cap), every patch z_of_u with monotone row heights, every boundary that is
%   neither a pole nor a clipped end shared with one neighbour as the same curve (contract S1),
%   normals into the void, mirrors as flipped copies, M1 and M2 (AGENTS section 5 item 9.4) on a
%   dense grid of points that are neither fitting nodes nor S2r check points (fractions j/7 of
%   every u span and j/9 of every v span), and M3 on horizontal sections at heights that are not
%   knot rows: one simple closed loop inside the outer section.

if nargin < 6
    opts = struct();
end
P = inner.patches;
want = {'t', 'd', 'eps_fit', 'z_range', 'z_lo', 'refit', 'patches', 'flat', 'report'};
check(isequal(sort(fieldnames(inner)'), sort(want)), '%s: S2 fields', name);
check(inner.eps_fit == 0.01 * t_min && inner.d == inner.t + inner.eps_fit / 2, '%s: d = t + eps_fit/2', name);
rf = {'n_nodes', 'n_knots', 'n_passes', 'n_removed', 'n_check', 't_local_min', 't_local_max', 'M1', 'M2', 'M3', 'M3_reason'};
check(isequal(sort(fieldnames(inner.report.patches)'), sort(rf)) && numel(inner.report.patches) == numel(P), ...
    '%s: S2r fields', name);
check(inner.report.ok && ~inner.report.cap_reached, '%s: fit report not ok', name);
zc = [-Inf Inf];
if isfield(opts, 'z_clip')
    zc = opts.z_clip;
end
dense = 1:numel(P);
if isfield(opts, 'dense')
    dense = opts.dense;
end

fprintf('%s: %d inner patches, z_lo %.17g, z_range %s\n', name, numel(P), inner.z_lo, mat2str(inner.z_range, 17));
fprintf('  %-28s %-8s %-8s %6s %7s %7s %13s %13s\n', 'patch', 'knots u', 'knots v', 'passes', 'removed', 'n_check', ...
    'tl_min [m]', 'tl_max [m]');
for k = 1:numel(P)
    r = inner.report.patches(k);
    fprintf('  %-28s %-8d %-8d %6d %7d %7d %13.10f %13.10f\n', P(k).name, r.n_knots(1), r.n_knots(2), r.n_passes, ...
        r.n_removed, r.n_check, r.t_local_min, r.t_local_max);
end

% rows of one z, monotone row heights; seams; mirrors
n_seam = 0;
for k = 1:numel(P)
    Z = P(k).surf.ctrl(:, :, 3);
    check(all(all(Z == Z(:, 1))), '%s: %s rows of more than one z', name, P(k).name);
    dz = diff(Z(:, 1));
    check(~(any(dz > 0) && any(dz < 0)), '%s: %s row heights not monotone', name, P(k).name);
    check(isequal(P(k).z_range, [Z(1, 1) Z(end, 1)]), '%s: %s z_range', name, P(k).name);
    for b = 1:4
        nb = P(k).(seam_field(b));
        c = boundary_curve(P(k).surf, b);
        if isempty(nb)
            collapsed = all(all(c.ctrl == c.ctrl(1, :)));
            clipped = any(b == [2 4]) && (c.ctrl(1, 3) >= zc(2) || c.ctrl(1, 3) <= zc(1));
            check(collapsed || clipped, '%s: %s boundary %d open', name, P(k).name, b);
            continue
        end
        check(isequal(P(nb(1)).(seam_field(nb(2))), [k b]), '%s: %s seam %d not mutual', name, P(k).name, b);
        check(same_curve(c, boundary_curve(P(nb(1)).surf, nb(2))), '%s: %s seam %d not the same curve', ...
            name, P(k).name, b);
        n_seam = n_seam + 1;
    end
    prim = find([P.visible] == min([P(strcmp({P.source}, P(k).source)).visible]));
    j = prim(strcmp(regexprep({P(prim).name}, '^.*_inner', ''), regexprep(P(k).name, '^.*_inner', '')));
    if ~isempty(j) && j ~= k
        S = P(j).surf.ctrl;
        rel = setxor(P(j).flips, P(k).flips);
        for f = 1:numel(rel)
            col = find(strcmp(rel{f}, {'X', 'Y'}));
            S(:, :, col) = -S(:, :, col);
        end
        check(isequal(S, P(k).surf.ctrl) && isequal(P(j).surf.knots, P(k).surf.knots) && ...
            isequal(P(j).surf.weights, P(k).surf.weights), '%s: %s is not the flipped copy of %s', name, P(k).name, P(j).name);
    end
end

% dense grid: M1, M2 and the normal direction
eps_fit = inner.eps_fit;
t = inner.t;
tl_all = zeros(0, 1);
n_into = 0;
for k = dense
    s = P(k).surf;
    us = fractions(s.knots{1}, 7);
    vs = fractions(s.knots{2}, 9);
    [UU, VV] = ndgrid(us, vs);
    [X, Su, Sv] = mwecmass.solid.eval_bspline_surface(s, UU(:), VV(:));
    tl = dist_fn(X);
    tl_all = [tl_all; tl]; %#ok<AGROW>
    n = cross(Su, Sv, 2);
    ln = sqrt(sum(n.^2, 2));
    ok = ln > 0;
    n = n(ok, :) ./ ln(ok);
    if ~P(k).outward
        n = -n;
    end
    % a step of eps_fit into the void raises t_local (the face is smooth on that scale away from a crease)
    i = find(ok);
    i = i(ceil(numel(i) / 2));
    j = find(find(ok) == i);
    into = dist_fn(X(i, :) + eps_fit * n(j, :)) > tl(i);
    check(into, '%s: %s normal does not point into the void', name, P(k).name);
    n_into = n_into + 1;
    m1 = all(tl >= t_min);
    m2 = all(tl >= t & tl <= t + eps_fit);
    check(m1 && m2, '%s: %s dense grid t_local [%.10f %.10f] outside [t, t + eps_fit] = [%.10f %.10f] (t_min %.10f)', ...
        name, P(k).name, min(tl), max(tl), t, t + eps_fit, t_min);
end
fprintf('  dense grid (%d points on %d patches, no fitting node, no S2r check point): t_local %.12f .. %.12f m, band [%.12f, %.12f]: M1, M2 pass\n', ...
    numel(tl_all), numel(dense), min(tl_all), max(tl_all), t, t + eps_fit);

% dense sections: M3
lat = arrayfun(@(e) e.z_range(1) ~= e.z_range(2), P);
zz = reshape([P(lat).z_range], 2, [])';
zlo = max(min(zz(:)), zc(1));
zhi = min(max(zz(:)), zc(2));
rows = zeros(0, 1);
for k = find(lat)
    ku = unique(P(k).surf.knots{1});
    Q = mwecmass.solid.eval_bspline_surface(P(k).surf, ku(:), repmat(P(k).surf.knots{2}(1), numel(ku), 1));
    rows = [rows; Q(:, 3)]; %#ok<AGROW>
end
hz = zlo + (1:79)' / 80 * (zhi - zlo);
hz = hz(~ismember(hz, rows));
lato = arrayfun(@(e) e.z_range(1) ~= e.z_range(2), geo.outer);
amin = Inf;
for h = hz'
    Li = mwecmass.solid.slice_bspline_surface(P(lat), h);
    Lo = mwecmass.solid.slice_bspline_surface(geo.outer(lato), h);
    check(Li.simple, '%s: section at z = %.17g not simple', name, h);
    check(all(inpolygon(Li.pts(:, 1), Li.pts(:, 2), Lo.pts(:, 1), Lo.pts(:, 2))), '%s: section at z = %.17g leaves the outer section', name, h);
    check(Li.area > 0, '%s: section at z = %.17g has no area', name, h);
    amin = min(amin, Li.area / Lo.area);
end
fprintf('  M3: %d seam boundaries mutual and the same curve; %d sections (z %.4f .. %.4f) one simple loop inside the outer section, smallest inner/outer area %.4f; %d normals into the void\n', ...
    n_seam, numel(hz), hz(1), hz(end), amin, n_into);
st = struct('tl_min', min(tl_all), 'tl_max', max(tl_all), 'n_dense', numel(tl_all), 'n_sections', numel(hz), 'n_seam', n_seam);
end

function s = fractions(k, m)
ku = unique(k);
s = zeros(0, 1);
for j = 1:numel(ku) - 1
    s = [s; ku(j) + (1:m - 1)' / m * (ku(j + 1) - ku(j))]; %#ok<AGROW>
end
end

function f = seam_field(b)
fields = {'seam_v0', 'seam_u1', 'seam_v1', 'seam_u0'};
f = fields{b};
end

function tf = same_curve(a, b)
% contract S1: degree, bitwise control points and weights (empty = ones), knots equal after the
% affine map to [0, 1] (4 * 2^-52: three correctly rounded operations); either direction
wa = a.weights;
wb = b.weights;
ka = (a.knots - a.knots(1)) / (a.knots(end) - a.knots(1));
kb = (b.knots - b.knots(1)) / (b.knots(end) - b.knots(1));
kr = fliplr(1 - kb);
tf = a.degree == b.degree && numel(ka) == numel(kb) && ...
    ((isequal(a.ctrl, b.ctrl) && isequal(wa, wb) && all(abs(ka - kb) <= 4 * 2^-52)) || ...
    (isequal(a.ctrl, flipud(b.ctrl)) && isequal(wa, flipud(wb)) && all(abs(ka - kr) <= 4 * 2^-52)));
end

function c = boundary_curve(s, b)
W = s.weights;
if isempty(W)
    W = ones(size(s.ctrl, 1), size(s.ctrl, 2));
end
switch b
    case 1, c = struct('degree', s.degree(1), 'ctrl', reshape(s.ctrl(:, 1, :), [], 3), 'knots', s.knots{1}, 'weights', W(:, 1));
    case 2, c = struct('degree', s.degree(2), 'ctrl', reshape(s.ctrl(end, :, :), [], 3), 'knots', s.knots{2}, 'weights', W(end, :)');
    case 3, c = struct('degree', s.degree(1), 'ctrl', reshape(s.ctrl(:, end, :), [], 3), 'knots', s.knots{1}, 'weights', W(:, end));
    case 4, c = struct('degree', s.degree(2), 'ctrl', reshape(s.ctrl(1, :, :), [], 3), 'knots', s.knots{2}, 'weights', W(1, :)');
end
end

function check(cond, varargin)
if ~cond
    error('inner_set_checks:fail', varargin{:});
end
end
