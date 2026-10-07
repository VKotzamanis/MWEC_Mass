function loop = slice_bspline_surface(patches, z)
%SLICE_BSPLINE_SURFACE  Stand-in of contract F4: the S5 section loop of u-degree-1 patches at height z.
%
%   loop = mwecmass.solid.slice_bspline_surface(patches, z)
%
%   patches: S1 array (only .surf is read). The iso-u row at z is the control row interpolated as in
%   split_bspline_surface (closed form for u-degree 1); constant-z patches are skipped. Rows are
%   chained by bitwise-equal end points into one closed loop, oriented counter-clockwise seen from
%   +z, else mwecmass:solid:SectionNotClosed.
%   Area, centroid and I = [Ixx Iyy Ixy] = [int y^2, int x^2, int x y] dA about the origin by
%   Green's theorem: A = loop int x dy, int x dA = 1/2 loop int x^2 dy, int y dA = -1/2 loop int
%   y^2 dx, int x^2 dA = 1/3 loop int x^3 dy, int y^2 dA = -1/3 loop int y^3 dx, int x y dA =
%   1/2 loop int x^2 y dy. Closed forms per piece: a degree-1 segment is polynomial (Gauss-Legendre
%   with 3 nodes, exact to degree 5); a rational quadratic with weights (1, cos a, 1) and equal legs
%   is a circular arc of half-angle a (Piegl & Tiller, The NURBS Book, 2nd ed., section 7.5), whose
%   integrands are trigonometric polynomials integrated exactly term by term. Any other curve, a
%   patch of u-degree other than 1, or one that is not z_of_u errors mwecmass:standin:NotAnalytic.

pieces = struct('patch', {}, 'u', {}, 'curve', {}, 'dir', {});
for k = 1:numel(patches)
    surf = patches(k).surf;
    if surf.degree(1) ~= 1
        error('mwecmass:standin:NotAnalytic', 'slice_bspline_surface stand-in: patch %d has u-degree %d', k, surf.degree(1));
    end
    [nu, nv, ~] = size(surf.ctrl);
    Z = surf.ctrl(:, :, 3);
    if any(any(Z ~= Z(:, 1)))
        error('mwecmass:standin:NotAnalytic', 'slice_bspline_surface stand-in: patch %d is not z_of_u', k);
    end
    Z = Z(:, 1);
    if all(Z == Z(1))
        continue
    end
    W = surf.weights;
    if isempty(W)
        W = ones(nu, nv);
    end
    t = surf.knots{1}(:)';
    i = find(Z(1:end - 1) ~= Z(2:end) & min(Z(1:end - 1), Z(2:end)) <= z & max(Z(1:end - 1), Z(2:end)) >= z, 1);
    if isempty(i)
        continue
    end
    a0 = W(i, 1);
    a1 = W(i + 1, 1);
    s = a0 * (z - Z(i)) / (a0 * (z - Z(i)) + a1 * (Z(i + 1) - z));
    if s == 0
        row = surf.ctrl(i, :, :);
        wrow = W(i, :);
    elseif s == 1
        row = surf.ctrl(i + 1, :, :);
        wrow = W(i + 1, :);
    elseif isequal(W(i, :), W(i + 1, :))
        row = (1 - s) * surf.ctrl(i, :, :) + s * surf.ctrl(i + 1, :, :);
        wrow = W(i, :);
    else
        wrow = (1 - s) * W(i, :) + s * W(i + 1, :);
        row = ((1 - s) * surf.ctrl(i, :, :) .* W(i, :) + s * surf.ctrl(i + 1, :, :) .* W(i + 1, :)) ./ wrow;
    end
    curve = struct('degree', surf.degree(2), 'ctrl', reshape(row, nv, 3), 'knots', surf.knots{2}, ...
        'weights', []);
    if ~isempty(surf.weights)
        curve.weights = wrow(:);
    end
    pieces(end + 1) = struct('patch', k, 'u', t(i + 1) + s * (t(i + 2) - t(i + 1)), 'curve', curve, 'dir', 1); %#ok<AGROW>
end
if isempty(pieces)
    error('mwecmass:solid:SectionNotClosed', 'slice_bspline_surface: no patch crosses z = %g', z);
end

% chain by bitwise-equal end points
n = numel(pieces);
order = 1;
dirs = 1;
used = false(1, n);
used(1) = true;
cur = pieces(1).curve.ctrl(end, :);
first = pieces(1).curve.ctrl(1, :);
while ~isequal(cur, first) || numel(order) < n
    found = false;
    for k = find(~used)
        c = pieces(k).curve.ctrl;
        if isequal(c(1, :), cur)
            order(end + 1) = k; dirs(end + 1) = 1; cur = c(end, :); %#ok<AGROW>
        elseif isequal(c(end, :), cur)
            order(end + 1) = k; dirs(end + 1) = -1; cur = c(1, :); %#ok<AGROW>
        else
            continue
        end
        used(k) = true;
        found = true;
        break
    end
    if ~found
        error('mwecmass:solid:SectionNotClosed', 'slice_bspline_surface: section at z = %g is not one closed loop', z);
    end
end
pieces = pieces(order);
for k = 1:n
    pieces(k).dir = dirs(k);
end

g = zeros(1, 6);
for k = 1:n
    g = g + pieces(k).dir * green(pieces(k).curve);
end
if g(1) < 0
    pieces = pieces(end:-1:1);
    for k = 1:n
        pieces(k).dir = -pieces(k).dir;
    end
    g = -g;
end
m = 8;
pts = zeros(0, 3);
for k = 1:n
    c = pieces(k).curve;
    s = c.knots(1) + (0:m - 1)' / m * (c.knots(end) - c.knots(1));
    if pieces(k).dir < 0
        s = c.knots(end) - (0:m - 1)' / m * (c.knots(end) - c.knots(1));
    end
    pts = [pts; mwecmass.solid.eval_bspline_curve(c, s)]; %#ok<AGROW>
end
loop = struct('z', z, 'pieces', pieces, 'pts', pts, 'area', g(1), 'centroid', g(2:3) / g(1), ...
    'I', g([5 4 6]), 'simple', true);
end

function g = green(c)
% [A, Sx, Sy, Iyy, Ixx, Ixy] contributions of one curve in its own direction.
P = c.ctrl(:, 1:2);
if c.degree == 1 && isempty(c.weights) || c.degree == 1 && all(c.weights == c.weights(1))
    g = zeros(1, 6);
    xg = [-sqrt(3/5) 0 sqrt(3/5)];
    wg = [5 8 5] / 9;
    for j = 1:size(P, 1) - 1
        p = P(j, :);
        q = P(j + 1, :);
        s = (xg + 1) / 2;
        x = p(1) + s * (q(1) - p(1));
        y = p(2) + s * (q(2) - p(2));
        dx = q(1) - p(1);
        dy = q(2) - p(2);
        f = [x * dy; x.^2 / 2 * dy; -y.^2 / 2 * dx; x.^3 / 3 * dy; -y.^3 / 3 * dx; x.^2 .* y / 2 * dy];
        g = g + (f * wg' / 2)';
    end
    return
end
w = c.weights;
if c.degree ~= 2 || size(P, 1) ~= 3 || isempty(w) || w(1) ~= w(3) || ...
        sum((P(2, :) - P(1, :)).^2) ~= sum((P(2, :) - P(3, :)).^2) || any(c.ctrl(:, 3) ~= c.ctrl(1, 3))
    error('mwecmass:standin:NotAnalytic', 'slice_bspline_surface stand-in: curve is neither a line nor a circular arc');
end
ca = w(2) / w(1);
if ca == sqrt(2) / 2
    cen = P(1, :) + P(3, :) - P(2, :);
else
    M = (P(1, :) + P(3, :)) / 2;
    L = norm(P(2, :) - P(1, :));
    cen = P(2, :) + (M - P(2, :)) * (L / (sqrt(1 - ca^2) * norm(M - P(2, :))));
end
rho = norm(P(1, :) - cen);
a = P(1, :) - cen;
b = P(3, :) - cen;
th0 = atan2(a(2), a(1));
sweep = 2 * acos(ca) * sign(a(1) * b(2) - a(2) * b(1));
th1 = th0 + sweep;
% trigonometric polynomials as coefficient vectors over e^{i k theta}, k = -1..1
X = [rho / 2, cen(1), rho / 2];
Y = [1i * rho / 2, cen(2), -1i * rho / 2];
DX = [-1i * rho / 2, 0, 1i * rho / 2];
DY = [rho / 2, 0, rho / 2];
I = @(F) integ(F, th0, th1);
g = real([I(conv(X, DY)), I(conv(conv(X, X), DY)) / 2, -I(conv(conv(Y, Y), DX)) / 2, ...
    I(conv(conv(conv(X, X), X), DY)) / 3, -I(conv(conv(conv(Y, Y), Y), DX)) / 3, ...
    I(conv(conv(conv(X, X), Y), DY)) / 2]);
end

function v = integ(F, a, b)
K = (numel(F) - 1) / 2;
v = 0;
for k = -K:K
    if k == 0
        v = v + F(K + 1) * (b - a);
    else
        v = v + F(k + K + 1) * (exp(1i * k * b) - exp(1i * k * a)) / (1i * k);
    end
end
end
