function [C, Cs] = eval_bspline_curve(curve, s)
%EVAL_BSPLINE_CURVE  Stand-in of contract F3 (curve): point and first derivative of a (rational) B-spline.
%
%   [C, Cs] = mwecmass.solid.eval_bspline_curve(curve, s)
%
%   curve is a T9 curve struct (degree, ctrl [n x 3], knots, weights [] or [n x 1]); s is a column
%   vector of parameters inside the knot range. C and Cs are [numel(s) x 3].
%   Source: basis functions and derivatives by the Cox-de Boor recursion, Piegl & Tiller, The NURBS
%   Book, 2nd ed., algorithms A2.1-A2.3; rational quotient rule, eq. (4.8).
%   The evaluation is exact up to rounding for every NURBS, so this stand-in takes no geo and has
%   no analytic requirement (its signature carries neither a geo nor a patch).

s = s(:);
P = curve.ctrl;
n = size(P, 1);
w = curve.weights(:);
if isempty(w)
    w = ones(n, 1);
end
[N, dN] = mwecmass_standin_basis(curve.knots(:)', curve.degree, n, s);
H = N * (P .* w);
W = N * w;
dH = dN * (P .* w);
dW = dN * w;
C = H ./ W;
Cs = (dH - C .* dW) ./ W;
end

function [N, dN] = mwecmass_standin_basis(t, p, n, s)
% Rows: parameters, columns: the n basis functions of degree p and their first derivatives.
m = numel(s);
N = zeros(m, n);
dN = zeros(m, n);
for k = 1:m
    x = s(k);
    span = find(t(1:n) <= x, 1, 'last');
    if x >= t(n + 1)
        span = n;
        while t(span) == t(span + 1)
            span = span - 1;
        end
    end
    % Piegl & Tiller A2.3 with n = 1 derivative, 1-based: span i has t(i) <= x < t(i+1)
    ndu = zeros(p + 1, p + 1);
    left = zeros(1, p + 1);
    right = zeros(1, p + 1);
    ndu(1, 1) = 1;
    for j = 1:p
        left(j + 1) = x - t(span + 1 - j);
        right(j + 1) = t(span + j) - x;
        saved = 0;
        for r = 0:j - 1
            ndu(j + 1, r + 1) = right(r + 2) + left(j - r + 1);
            tmp = ndu(r + 1, j) / ndu(j + 1, r + 1);
            ndu(r + 1, j + 1) = saved + right(r + 2) * tmp;
            saved = left(j - r + 1) * tmp;
        end
        ndu(j + 1, j + 1) = saved;
    end
    vals = ndu(:, p + 1)';
    ders = zeros(1, p + 1);
    if p >= 1
        for r = 0:p
            d = 0;
            if r >= 1
                d = d + ndu(r, p) / ndu(p + 1, r);
            end
            if r <= p - 1
                d = d - ndu(r + 1, p) / ndu(p + 1, r + 1);
            end
            ders(r + 1) = p * d;
        end
    end
    cols = span - p:span;
    N(k, cols) = vals;
    dN(k, cols) = ders;
end
end
