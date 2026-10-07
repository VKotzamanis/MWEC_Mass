function [S, Su, Sv] = eval_bspline_surface(surf, u, v)
%EVAL_BSPLINE_SURFACE  Stand-in of contract F3 (surface): point and partial derivatives of a T9 'bspline' surface.
%
%   [S, Su, Sv] = mwecmass.solid.eval_bspline_surface(surf, u, v)
%
%   surf: degree [du dv], ctrl [nu x nv x 3], knots {ku, kv}, weights [] or [nu x nv]; u, v column
%   vectors of equal length. S, Su, Sv are [numel(u) x 3].
%   Source: tensor-product rational B-spline, Piegl & Tiller, The NURBS Book, 2nd ed., eq. (4.15)
%   with the basis of algorithms A2.2-A2.3 and the quotient rule of eq. (4.19).
%   Exact up to rounding for every NURBS surface, so this stand-in has no analytic requirement
%   (its signature carries neither a geo nor a patch).

u = u(:);
v = v(:);
[nu, nv, ~] = size(surf.ctrl);
W = surf.weights;
if isempty(W)
    W = ones(nu, nv);
end
cu = struct('degree', surf.degree(1), 'ctrl', zeros(nu, 3), 'knots', surf.knots{1}, 'weights', []);
cv = struct('degree', surf.degree(2), 'ctrl', zeros(nv, 3), 'knots', surf.knots{2}, 'weights', []);
[Nu, dNu] = basis(cu, nu, u);
[Nv, dNv] = basis(cv, nv, v);
m = numel(u);
S = zeros(m, 3);
Su = zeros(m, 3);
Sv = zeros(m, 3);
for k = 1:m
    a = Nu(k, :);
    da = dNu(k, :);
    b = Nv(k, :)';
    db = dNv(k, :)';
    w = a * W * b;
    wu = da * W * b;
    wv = a * W * db;
    H = zeros(1, 3);
    Hu = zeros(1, 3);
    Hv = zeros(1, 3);
    for c = 1:3
        Pw = surf.ctrl(:, :, c) .* W;
        H(c) = a * Pw * b;
        Hu(c) = da * Pw * b;
        Hv(c) = a * Pw * db;
    end
    S(k, :) = H / w;
    Su(k, :) = (Hu - S(k, :) * wu) / w;
    Sv(k, :) = (Hv - S(k, :) * wv) / w;
end
end

function [N, dN] = basis(curve, n, s)
% Basis values and derivatives from the curve evaluator: a curve with identity control points
% returns the basis functions themselves, one column of ctrl per basis function.
m = numel(s);
N = zeros(m, n);
dN = zeros(m, n);
for i = 1:n
    e = zeros(n, 3);
    e(i, 1) = 1;
    curve.ctrl = e;
    [C, Cs] = mwecmass.solid.eval_bspline_curve(curve, s);
    N(:, i) = C(:, 1);
    dN(:, i) = Cs(:, 1);
end
end
