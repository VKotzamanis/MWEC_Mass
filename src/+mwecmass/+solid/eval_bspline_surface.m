function [S, Su, Sv, Suu, Suv, Svv] = eval_bspline_surface(surf, u, v)
%EVAL_BSPLINE_SURFACE  Point and partial derivatives of a (rational) tensor-product B-spline surface.
%
%   [S, Su, Sv] = mwecmass.solid.eval_bspline_surface(surf, u, v)
%   [S, Su, Sv, Suu, Suv, Svv] = mwecmass.solid.eval_bspline_surface(surf, u, v)
%
%   surf: T9 'bspline' surface (degree [du dv], ctrl [nu x nv x 3], knots {ku, kv}, weights [] or
%   [nu x nv]); u, v: parameters of equal length. Outputs are [numel(u) x 3]. Parameters outside
%   the knot ranges are evaluated on the first or last span.
%   Tensor-product basis from eval_bspline_curve (Piegl & Tiller, The NURBS Book, 2nd ed., A2.3);
%   rational surfaces by the quotient rule on the homogeneous surface (ibid., eq. (4.19)).

u = u(:);
v = v(:);
[nu, nv, ~] = size(surf.ctrl);
nd = 0;
if nargout > 1, nd = 1; end
if nargout > 3, nd = 2; end
Bu = basis(surf.degree(1), surf.knots{1}, nu, u, nd);
Bv = basis(surf.degree(2), surf.knots{2}, nv, v, nd);
W = surf.weights;
rational = ~isempty(W);
if ~rational
    W = ones(nu, nv);
end
m = numel(u);
H = cell(1, 6);
for k = 1:6
    H{k} = zeros(m, 4);
end
for c = 1:4
    if c <= 3
        Pc = surf.ctrl(:, :, c) .* W;
    else
        Pc = W;
    end
    A0 = Bu{1} * Pc;
    H{1}(:, c) = sum(A0 .* Bv{1}, 2);
    if nd >= 1
        A1 = Bu{2} * Pc;
        H{2}(:, c) = sum(A1 .* Bv{1}, 2);
        H{3}(:, c) = sum(A0 .* Bv{2}, 2);
    end
    if nd >= 2
        H{4}(:, c) = sum((Bu{3} * Pc) .* Bv{1}, 2);
        H{5}(:, c) = sum(A1 .* Bv{2}, 2);
        H{6}(:, c) = sum(A0 .* Bv{3}, 2);
    end
end
if ~rational
    S = H{1}(:, 1:3);
    if nd >= 1
        Su = H{2}(:, 1:3);
        Sv = H{3}(:, 1:3);
    end
    if nd >= 2
        Suu = H{4}(:, 1:3);
        Suv = H{5}(:, 1:3);
        Svv = H{6}(:, 1:3);
    end
    return
end
w0 = H{1}(:, 4);
S = H{1}(:, 1:3) ./ w0;
if nd >= 1
    wu = H{2}(:, 4);
    wv = H{3}(:, 4);
    Su = (H{2}(:, 1:3) - S .* wu) ./ w0;
    Sv = (H{3}(:, 1:3) - S .* wv) ./ w0;
end
if nd >= 2
    Suu = (H{4}(:, 1:3) - 2 * Su .* wu - S .* H{4}(:, 4)) ./ w0;
    Suv = (H{5}(:, 1:3) - Su .* wv - Sv .* wu - S .* H{5}(:, 4)) ./ w0;
    Svv = (H{6}(:, 1:3) - 2 * Sv .* wv - S .* H{6}(:, 4)) ./ w0;
end
end

function B = basis(p, knots, n, s, nd)
basis_curve = struct('degree', p, 'ctrl', eye(n), 'knots', knots, 'weights', []);
B = cell(1, 3);
switch nd
    case 0
        B{1} = mwecmass.solid.eval_bspline_curve(basis_curve, s);
    case 1
        [B{1}, B{2}] = mwecmass.solid.eval_bspline_curve(basis_curve, s);
    otherwise
        [B{1}, B{2}, B{3}] = mwecmass.solid.eval_bspline_curve(basis_curve, s);
end
end
