function hs = hydrostatics_at_draft(geo, vs, opts)
%HYDROSTATICS_AT_DRAFT  Hydrostatics of the outer surface at one vertical shift (contract F7, S7).
%
%   hs = mwecmass.solid.hydrostatics_at_draft(geo, vs, opts)
%
%   geo: S1b (outer_nurbs); vs [m]: vertical shift, the waterline is z_body = -vs. opts.n_gauss:
%   Gauss-Legendre points per knot span (default as body_properties).
%   The outer pieces below the waterline (lateral pieces cut there by split_bspline_surface) are
%   integrated with body_properties (divergence theorem; the waterplane has n_x = 0 and adds 0):
%   V_sub and CB. S_wet integrates |S_u x S_v| over the same pieces and over the constant-z pieces
%   strictly below the waterline, with the same Gauss rule, and adds the area of every flat region
%   (geo.flat, general path) strictly below it, by Green's theorem on the rows that bound it. The waterplane is
%   slice_bspline_surface (F4) of the lateral pieces reaching the waterline from below: Aw, centre
%   of flotation (xF, yF), I_wp_xx = int y^2 dA - Aw yF^2 and I_wp_yy = int x^2 dA - Aw xF^2
%   (about the centre of flotation), I_wp_yy_origin = int x^2 dA, BM_L = I_wp_yy / V_sub, KM =
%   CB(3) + BM_L (world z). CB_body is body frame, CB world (CB_body + [0 0 vs]).
%   Waterline at or above z_max: submersion 'full', V_sub = V_hull, CB = hull centroid, Aw,
%   I_wp_*, BM_L = 0, KM = CB(3), xF = yF = NaN. At or below z_min: 'none', V_sub = S_wet = Aw = 0,
%   CB, xF, yF, BM_L, KM = NaN. waterline is empty in both cases.

if nargin < 3 || isempty(opts)
    opts = struct();
end
zw = -vs;
zmin = geo.z_range(1);
zmax = geo.z_range(2);
hs = struct('vs', vs, 'draft', -(zmin + vs), 'submersion', 'partial', 'V_sub', 0, ...
    'CB_body', NaN(1, 3), 'CB', NaN(1, 3), 'S_wet', 0, 'Aw', 0, 'xF', NaN, 'yF', NaN, ...
    'I_wp_xx', 0, 'I_wp_yy', 0, 'I_wp_yy_origin', 0, 'BM_L', NaN, 'KM', NaN, 'waterline', []);
if zw <= zmin
    hs.submersion = 'none';
    return
end
full = zw >= zmax;
lat = cell(1, 0);
sgn = zeros(1, 0);
flat = cell(1, 0);
lateral = false(1, numel(geo.outer));
for k = 1:numel(geo.outer)
    p = geo.outer(k);
    zr = p.z_range;
    if zr(1) == zr(2)
        if zr(1) < zw || full
            flat{end + 1} = p.surf; %#ok<AGROW>
        end
        continue
    end
    lateral(k) = true;
    if min(zr) >= zw
        continue
    end
    if full || max(zr) <= zw
        q = p;
    else
        [lo, hi] = mwecmass.solid.split_bspline_surface(p, zw);
        if zr(1) < zr(2)
            q = lo;
        else
            q = hi;
        end
    end
    lat{end + 1} = q.surf; %#ok<AGROW>
    sgn(end + 1) = q.outward; %#ok<AGROW>
end

faces = struct('surface', num2cell(1:numel(lat)), 'same_sense', num2cell(logical(sgn)), ...
    'loops', {{}}, 'role', 'outer', 'module', 1, 'inside', 'submerged', 'outside', 'exterior');
sub = struct('design', struct('mode', 'displacement', 'edges', [zmin; zmax]), ...
    'brep', struct('surfaces', {lat}, 'faces', faces));
bp = mwecmass.solid.body_properties(sub, struct('submerged', 1), opts);
R = bp.regions.submerged;
ng = bp.quad;
hs.V_sub = R.V(1);
hs.CB_body = R.S(1, :) / R.V(1);
hs.CB = hs.CB_body + [0 0 vs];
area = 0;
for k = 1:numel(lat)
    area = area + face_area(lat{k}, ng);
end
for k = 1:numel(flat)
    area = area + face_area(flat{k}, ng);
end
if isfield(geo, 'flat')
    for j = 1:numel(geo.flat)
        if geo.flat(j).z < zw || full
            area = area + flat_area(geo, j, ng);
        end
    end
end
hs.S_wet = area;
if full
    hs.submersion = 'full';
    hs.BM_L = 0;
    hs.KM = hs.CB(3);
    return
end
reach = lateral & arrayfun(@(p) min(p.z_range) < zw, geo.outer);
L = mwecmass.solid.slice_bspline_surface(geo.outer(reach), zw);
hs.waterline = L;
hs.Aw = L.area;
hs.xF = L.centroid(1);
hs.yF = L.centroid(2);
hs.I_wp_xx = L.I(1) - L.area * hs.yF^2;
hs.I_wp_yy = L.I(2) - L.area * hs.xF^2;
hs.I_wp_yy_origin = L.I(2);
hs.BM_L = hs.I_wp_yy / hs.V_sub;
hs.KM = hs.CB(3) + hs.BM_L;
end

function a = flat_area(geo, j, ng)
% area of flat region j from the rows that name it (S1 seam [0 j]) by Green's theorem: a row
% bounds the flat face in the direction opposite to its use in the loop of its lateral face, so
% with outward normals both ways the u0 row runs along +v when S_u x S_v points out of the hull and
% the u1 row along -v; the loop is then counter-clockwise seen from the flat's normal [0 0 normal_z]
[xg, wg] = gauss_legendre(ng);
a = 0;
for k = 1:numel(geo.outer)
    p = geo.outer(k);
    s = p.surf;
    rows = {'seam_u0', 1, 1; 'seam_u1', size(s.ctrl, 1), -1};
    for r = 1:2
        sm = p.(rows{r, 1});
        if numel(sm) ~= 2 || sm(1) ~= 0 || sm(2) ~= j
            continue
        end
        i = rows{r, 2};
        w = [];
        if ~isempty(s.weights)
            w = s.weights(i, :)';
        end
        c = struct('degree', s.degree(2), 'ctrl', reshape(s.ctrl(i, :, :), [], 3), 'knots', s.knots{2}, 'weights', w);
        kv = unique(c.knots);
        g = 0;
        for q = 1:numel(kv) - 1
            t = (kv(q) + kv(q + 1)) / 2 + (kv(q + 1) - kv(q)) / 2 * xg;
            [C, Cs] = mwecmass.solid.eval_bspline_curve(c, t);
            g = g + (kv(q + 1) - kv(q)) / 2 * (wg' * (C(:, 1) .* Cs(:, 2)));
        end
        a = a + rows{r, 3} * (2 * p.outward - 1) * g;
    end
end
a = a * geo.flat(j).normal_z;
end

function a = face_area(s, ng)
% integral of |S_u x S_v| over the face (Gauss-Legendre, ng points per knot span in u and v)
[xg, wg] = gauss_legendre(ng);
ku = unique(s.knots{1});
kv = unique(s.knots{2});
au = ku(1:end - 1);
bu = ku(2:end);
av = kv(1:end - 1);
bv = kv(2:end);
uu = (au(:)' + bu(:)') / 2 + (bu(:)' - au(:)') / 2 .* xg;
vv = (av(:)' + bv(:)') / 2 + (bv(:)' - av(:)') / 2 .* xg;
wu = (bu(:)' - au(:)') / 2 .* wg;
wv = (bv(:)' - av(:)') / 2 .* wg;
[U, V] = ndgrid(uu(:), vv(:));
[Wu, Wv] = ndgrid(wu(:), wv(:));
[~, Su, Sv] = mwecmass.solid.eval_bspline_surface(s, U(:), V(:));
n = cross(Su, Sv, 2);
a = (Wu(:) .* Wv(:))' * sqrt(sum(n.^2, 2));
end

function [x, w] = gauss_legendre(n)
% Golub-Welsch: nodes and weights of n-point Gauss-Legendre on [-1, 1] (columns)
b = (1:n - 1) ./ sqrt(4 * (1:n - 1).^2 - 1);
[Vv, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D));
w = 2 * Vv(1, i)'.^2;
end
