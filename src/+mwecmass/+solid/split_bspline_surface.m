function [lo, hi, u_star] = split_bspline_surface(patch, z)
%SPLIT_BSPLINE_SURFACE  Cut a z_of_u patch at height z by knot insertion (contract F3b).
%
%   [lo, hi, u_star] = mwecmass.solid.split_bspline_surface(patch, z)
%
%   patch: S1 entry with z_of_u true (every control row has one z, weights separable), z strictly
%   between z(u0) and z(u1), in either order. u* solves z(u*) = z, where z(u) = sum N_i a_i z_i /
%   sum N_i a_i (rows z_i, row weights a_i); u* is bracketed between knot values and refined by
%   regula falsi (Illinois) to adjacent doubles. u* is inserted to multiplicity du (Piegl & Tiller, The NURBS
%   Book, 2nd ed., algorithm A5.1, on homogeneous rows); a coordinate or weight that is equal in
%   the two rows being combined is copied, every new row takes one z (computed from column 1) and
%   the cut row takes z itself, so lo and hi keep z_of_u and share the cut row bitwise. lo is the
%   piece [u0, u*], hi the piece [u*, u1]; both keep the parent parameter.
%   Errors: mwecmass:solid:ZNotOneParameter (patch not z_of_u), mwecmass:solid:ZOutside (z not
%   strictly inside the z range), mwecmass:solid:ZNotMonotonic (z(u) = z has more than one root).

surf = patch.surf;
nu = size(surf.ctrl, 1);
Z = surf.ctrl(:, :, 3);
if ~patch.z_of_u || any(any(Z ~= Z(:, 1)))
    error('mwecmass:solid:ZNotOneParameter', 'split_bspline_surface: %s is not z_of_u', patch.name);
end
Zr = Z(:, 1);
if ~((z > Zr(1) && z < Zr(end)) || (z < Zr(1) && z > Zr(end)))
    error('mwecmass:solid:ZOutside', 'split_bspline_surface: z = %.17g not strictly inside (%.17g, %.17g)', ...
        z, Zr(1), Zr(end));
end
W = surf.weights;
if isempty(W)
    a = ones(nu, 1);
else
    a = W(:, 1);
end
p = surf.degree(1);
knots = surf.knots{1}(:)';
u_star = z_root(Zr, a, knots, p, z, patch.name);

[ctrl, Wn, knots_new] = insert_rows(surf.ctrl, W, knots, p, u_star, p - sum(knots == u_star));
f = find(knots_new == u_star, 1);
ctrl(f - 1, :, 3) = z;
lo = patch;
hi = patch;
lo.surf.ctrl = ctrl(1:f - 1, :, :);
hi.surf.ctrl = ctrl(f - 1:end, :, :);
lo.surf.knots{1} = [knots_new(1:f - 1), repmat(u_star, 1, p + 1)];
hi.surf.knots{1} = [repmat(u_star, 1, p + 1), knots_new(f + p:end)];
if isempty(W)
    lo.surf.weights = [];
    hi.surf.weights = [];
else
    lo.surf.weights = Wn(1:f - 1, :);
    hi.surf.weights = Wn(f - 1:end, :);
end
lo.u_range = [patch.u_range(1) u_star];
hi.u_range = [u_star patch.u_range(2)];
lo.z_range = [patch.z_range(1) z];
hi.z_range = [z patch.z_range(2)];
lo.pole = [patch.pole(1) collapsed(lo.surf.ctrl(end, :, :))];
hi.pole = [collapsed(hi.surf.ctrl(1, :, :)) patch.pole(2)];
lo.c0_u = patch.c0_u(patch.c0_u > lo.u_range(1) & patch.c0_u < u_star);
hi.c0_u = patch.c0_u(patch.c0_u > u_star & patch.c0_u < hi.u_range(2));
lo.seam_u1 = [];
hi.seam_u0 = [];
end

function u_star = z_root(Zr, a, knots, p, z, name)
% Unique root of z(u) = z, else ZNotMonotonic. c_i = a_i (z_i - z) are the coefficients of the
% numerator of z(u) - z; it has at most as many sign changes as they have (variation diminishing),
% and it vanishes on a whole span only if every coefficient active there is zero.
c = Zr - z;
sg = sign(c);
nzs = sg(sg ~= 0);
if sum(diff(nzs) ~= 0) ~= 1
    error('mwecmass:solid:ZNotMonotonic', 'split_bspline_surface: z(u) = %.17g has more than one root on %s', z, name);
end
ku = unique(knots);
for j = 1:numel(ku) - 1
    span = find(knots <= ku(j), 1, 'last');
    if all(c(span - p:span) == 0)
        error('mwecmass:solid:ZNotMonotonic', 'split_bspline_surface: z(u) = %.17g on an interval of %s', z, name);
    end
end
zc = struct('degree', p, 'ctrl', c, 'knots', knots, 'weights', a);
f = mwecmass.solid.eval_bspline_curve(zc, ku(:));
f(1) = c(1);
f(end) = c(end);
hit = find(f == 0);
chg = find(sign(f(1:end - 1)) .* sign(f(2:end)) < 0);
if numel(hit) + numel(chg) ~= 1
    error('mwecmass:solid:ZNotMonotonic', 'split_bspline_surface: z(u) = %.17g has more than one root on %s', z, name);
end
if ~isempty(hit)
    u_star = ku(hit);
    return
end
lo = ku(chg);
hi = ku(chg + 1);
flo = f(chg);
fhi = f(chg + 1);
% Illinois regula falsi, kept bracketing, until the bracket is two adjacent doubles
side = 0;
for it = 1:200
    x = (flo * hi - fhi * lo) / (flo - fhi);
    if ~(x > lo && x < hi)
        x = lo + (hi - lo) / 2;
    end
    if x <= lo || x >= hi
        break
    end
    fx = mwecmass.solid.eval_bspline_curve(zc, x);
    if fx == 0
        u_star = x;
        return
    end
    if sign(fx) == sign(flo)
        lo = x;
        flo = fx;
        if side == 1
            fhi = fhi / 2;
        end
        side = 1;
    else
        hi = x;
        fhi = fx;
        if side == -1
            flo = flo / 2;
        end
        side = -1;
    end
end
if abs(fhi) < abs(flo)
    u_star = hi;
else
    u_star = lo;
end
end

function tf = collapsed(row)
row = reshape(row, [], 3);
tf = all(all(row == row(1, :)));
end

function [ctrl, W, knots] = insert_rows(ctrl, W, knots, p, u, r)
% Insert u r times into the u knot vector (Piegl & Tiller A5.1 on homogeneous rows). Where the two
% rows being combined have an equal coordinate (or weight) the result copies it; each new row
% takes the z of column 1, so rows keep one z.
if r <= 0
    return
end
rational = ~isempty(W);
[np, nv, ~] = size(ctrl);
if ~rational
    W = ones(np, nv);
end
k = find(knots <= u, 1, 'last');
if k > np
    k = np;
    while knots(k) == knots(k + 1)
        k = k - 1;
    end
end
s = sum(knots == u);
UQ = [knots(1:k), repmat(u, 1, r), knots(k + 1:end)];
Q = zeros(np + r, nv, 3);
QW = zeros(np + r, nv);
Q(1:k - p, :, :) = ctrl(1:k - p, :, :);
QW(1:k - p, :) = W(1:k - p, :);
Q(k - s + r:np + r, :, :) = ctrl(k - s:np, :, :);
QW(k - s + r:np + r, :) = W(k - s:np, :);
R = ctrl(k - p:k - s, :, :);
RW = W(k - p:k - s, :);
L = k - p;
for j = 1:r
    L = k - p + j;
    for i = 0:p - j - s
        alpha = (u - knots(L + i)) / (knots(i + k + 1) - knots(L + i));
        [R(i + 1, :, :), RW(i + 1, :)] = combine(R(i + 2, :, :), RW(i + 2, :), R(i + 1, :, :), RW(i + 1, :), alpha);
    end
    Q(L, :, :) = R(1, :, :);
    QW(L, :) = RW(1, :);
    Q(k + r - j - s, :, :) = R(p - j - s + 1, :, :);
    QW(k + r - j - s, :) = RW(p - j - s + 1, :);
end
for i = L + 1:k - s - 1
    Q(i, :, :) = R(i - L + 1, :, :);
    QW(i, :) = RW(i - L + 1, :);
end
ctrl = Q;
knots = UQ;
if rational
    W = QW;
else
    W = [];
end
end

function [P, w] = combine(P1, w1, P0, w0, alpha)
% alpha * (w1 P1, w1) + (1 - alpha) * (w0 P0, w0) in homogeneous form, back to Euclidean.
w = alpha * w1 + (1 - alpha) * w0;
same_w = w1 == w0;
w(same_w) = w1(same_w);
P = (alpha * P1 .* w1 + (1 - alpha) * P0 .* w0) ./ w;
same = P1 == P0;
P(same) = P1(same);
P(1, :, 3) = P(1, 1, 3);
end
