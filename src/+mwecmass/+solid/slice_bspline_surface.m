function loop = slice_bspline_surface(patches, z)
%SLICE_BSPLINE_SURFACE  Horizontal section of z_of_u patches at height z as one closed loop (contract F4, S5).
%
%   loop = mwecmass.solid.slice_bspline_surface(patches, z)
%
%   patches: S1 array. Constant-z patches (z_range(1) = z_range(2)) and patches whose z range does
%   not reach z are skipped. Each other patch contributes its iso-u row at z: the cut row of
%   split_bspline_surface (F3b), or its end row when z is an end of its z range; a row whose
%   control points all coincide (a pole) contributes nothing. The rows are chained end to end by
%   mutually nearest end points into one loop (bitwise-equal ends on shared seams), oriented
%   counter-clockwise seen from +z; otherwise mwecmass:solid:SectionNotClosed.
%   loop: z, pieces(k) (patch, u, curve, dir), pts [n x 3] (8 points per knot span of every piece,
%   counter-clockwise, no repeated point), area, centroid [1x2], I = [Ixx Iyy Ixy] =
%   [int y^2, int x^2, int x y] dA about the origin, simple (pts has no crossing segments).
%   Moments by Green's theorem: A = loop int x dy, int x dA = 1/2 loop int x^2 dy, int y dA =
%   -1/2 loop int y^2 dx, int x^2 dA = 1/3 loop int x^3 dy, int y^2 dA = -1/3 loop int y^3 dx,
%   int x y dA = 1/2 loop int x^2 y dy, each by 16-point Gauss-Legendre on every knot span
%   (exact for polynomial rows up to degree 8, convergent for rational ones).

pieces = struct('patch', {}, 'u', {}, 'curve', {}, 'dir', {});
for k = 1:numel(patches)
    pk = patches(k);
    zr = pk.z_range;
    if zr(1) == zr(2) || z < min(zr) || z > max(zr)
        continue
    end
    surf = pk.surf;
    if z == zr(1)
        u = pk.u_range(1);
        row = surf.ctrl(1, :, :);
        w = row_weights(surf, 1);
    elseif z == zr(2)
        u = pk.u_range(2);
        row = surf.ctrl(end, :, :);
        w = row_weights(surf, size(surf.ctrl, 1));
    else
        [~, hi, u] = mwecmass.solid.split_bspline_surface(pk, z);
        row = hi.surf.ctrl(1, :, :);
        w = row_weights(hi.surf, 1);
    end
    row = reshape(row, [], 3);
    if all(all(row == row(1, :)))
        continue
    end
    curve = struct('degree', surf.degree(2), 'ctrl', row, 'knots', surf.knots{2}, 'weights', w);
    pieces(end + 1) = struct('patch', k, 'u', u, 'curve', curve, 'dir', 1); %#ok<AGROW>
end
if isempty(pieces)
    error('mwecmass:solid:SectionNotClosed', 'slice_bspline_surface: no patch crosses z = %.17g', z);
end
pieces = chain(pieces, z);

g = zeros(1, 6);
for k = 1:numel(pieces)
    g = g + pieces(k).dir * green(pieces(k).curve);
end
if g(1) < 0
    pieces = pieces(end:-1:1);
    for k = 1:numel(pieces)
        pieces(k).dir = -pieces(k).dir;
    end
    g = -g;
end
pts = zeros(0, 3);
for k = 1:numel(pieces)
    c = pieces(k).curve;
    ku = unique(c.knots);
    s = zeros(0, 1);
    for j = 1:numel(ku) - 1
        s = [s; ku(j) + (0:7)' / 8 * (ku(j + 1) - ku(j))]; %#ok<AGROW>
    end
    if pieces(k).dir < 0
        s = flipud([s(2:end); ku(end)]);
    end
    q = mwecmass.solid.eval_bspline_curve(c, s);
    q(:, 3) = z;
    pts = [pts; q]; %#ok<AGROW>
end
keep = [true; any(diff(pts, 1, 1) ~= 0, 2)];
pts = pts(keep, :);
if size(pts, 1) > 1 && isequal(pts(end, :), pts(1, :))
    pts = pts(1:end - 1, :);
end
loop = struct('z', z, 'pieces', pieces, 'pts', pts, 'area', g(1), 'centroid', g(2:3) / g(1), ...
    'I', g([5 4 6]), 'simple', is_simple(pts));
end

function w = row_weights(surf, i)
if isempty(surf.weights)
    w = [];
else
    w = surf.weights(i, :)';
end
end

function pieces = chain(pieces, z)
% Order the pieces into one loop by mutually nearest end points (as T1 outer_rows joins its pieces).
m = numel(pieces);
E = zeros(2 * m, 3);
for k = 1:m
    E(2 * k - 1, :) = pieces(k).curve.ctrl(1, :);
    E(2 * k, :) = pieces(k).curve.ctrl(end, :);
end
if m == 1
    if ~isequal(E(1, :), E(2, :))
        error('mwecmass:solid:SectionNotClosed', 'slice_bspline_surface: the section at z = %.17g is open', z);
    end
    return
end
D = sqrt(sum((permute(E, [1 3 2]) - permute(E, [3 1 2])).^2, 3));
D(1:2 * m + 1:end) = Inf;
[~, partner] = min(D, [], 2);
if any(partner(partner) ~= (1:2 * m)')
    error('mwecmass:solid:SectionNotClosed', 'slice_bspline_surface: end points at z = %.17g do not pair up', z);
end
order = zeros(1, m);
dirs = zeros(1, m);
e = 1;
for it = 1:m
    k = ceil(e / 2);
    order(it) = k;
    if mod(e, 2) == 1
        dirs(it) = 1;
        ex = e + 1;
    else
        dirs(it) = -1;
        ex = e - 1;
    end
    e = partner(ex);
    if ceil(e / 2) == order(1)
        break
    end
end
if e ~= 1 || numel(unique(order(order > 0))) ~= m
    error('mwecmass:solid:SectionNotClosed', 'slice_bspline_surface: the section at z = %.17g is not one loop', z);
end
pieces = pieces(order);
for k = 1:m
    pieces(k).dir = dirs(k);
end
end

function g = green(c)
% [A, int x dA, int y dA, int x^2 dA, int y^2 dA, int x y dA] of one curve in its own direction.
[xg, wg] = gauss_legendre(16);
ku = unique(c.knots);
g = zeros(1, 6);
for j = 1:numel(ku) - 1
    a = ku(j);
    b = ku(j + 1);
    s = (a + b) / 2 + (b - a) / 2 * xg;
    [C, Cs] = mwecmass.solid.eval_bspline_curve(c, s);
    x = C(:, 1);
    y = C(:, 2);
    dx = Cs(:, 1);
    dy = Cs(:, 2);
    f = [x .* dy, x.^2 / 2 .* dy, -y.^2 / 2 .* dx, x.^3 / 3 .* dy, -y.^3 / 3 .* dx, x.^2 .* y / 2 .* dy];
    g = g + (b - a) / 2 * (wg' * f);
end
end

function [x, w] = gauss_legendre(n)
% Golub-Welsch: nodes and weights of n-point Gauss-Legendre on [-1, 1].
persistent cache
if isempty(cache)
    cache = cell(1, 64);
end
if ~isempty(cache{n})
    x = cache{n}(:, 1);
    w = cache{n}(:, 2);
    return
end
b = (1:n - 1) ./ sqrt(4 * (1:n - 1).^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D));
w = 2 * V(1, i)'.^2;
cache{n} = [x w];
end

function tf = is_simple(P)
% No two non-adjacent segments of the closed polygon P cross or touch.
n = size(P, 1);
A = P(:, 1:2);
B = A([2:n, 1], :);
tf = true;
for i = 1:n - 2
    j = i + 2:n;
    if i == 1
        j = j(j ~= n);
    end
    if isempty(j)
        continue
    end
    Ai = repmat(A(i, :), numel(j), 1);
    Bi = repmat(B(i, :), numel(j), 1);
    d1 = orient(A(j, :), B(j, :), Ai);
    d2 = orient(A(j, :), B(j, :), Bi);
    d3 = orient(Ai, Bi, A(j, :));
    d4 = orient(Ai, Bi, B(j, :));
    cross_ = d1 .* d2 < 0 & d3 .* d4 < 0;
    touch = (d1 == 0 & on_box(A(j, :), B(j, :), Ai)) | (d2 == 0 & on_box(A(j, :), B(j, :), Bi)) | ...
        (d3 == 0 & on_box(Ai, Bi, A(j, :))) | (d4 == 0 & on_box(Ai, Bi, B(j, :)));
    if any(cross_ | touch)
        tf = false;
        return
    end
end
end

function tf = on_box(a, b, c)
% c (collinear with segment a-b) lies on the segment
tf = c(:, 1) >= min(a(:, 1), b(:, 1)) & c(:, 1) <= max(a(:, 1), b(:, 1)) & ...
    c(:, 2) >= min(a(:, 2), b(:, 2)) & c(:, 2) <= max(a(:, 2), b(:, 2));
end

function d = orient(a, b, c)
d = (b(:, 1) - a(:, 1)) .* (c(:, 2) - a(:, 2)) - (b(:, 2) - a(:, 2)) .* (c(:, 1) - a(:, 1));
end
