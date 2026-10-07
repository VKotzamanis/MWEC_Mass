function [C, Cs, Css] = eval_bspline_curve(curve, s)
%EVAL_BSPLINE_CURVE  Point and parameter derivatives of a (rational) B-spline curve.
%
%   [C, Cs, Css] = mwecmass.solid.eval_bspline_curve(curve, s)
%
%   curve: T9 curve struct, degree p, ctrl [n x k] (k = 3 for a point curve; any k is accepted,
%   so ctrl = eye(n) returns the basis functions themselves), knots (clamped, length n+p+1),
%   weights [] or [n x 1]. s: parameters (column). C, Cs, Css: [numel(s) x k], the point and its
%   first and second derivative. A parameter outside the knot range is evaluated on the first or
%   last span (polynomial continuation).
%   Basis functions and derivatives: Piegl & Tiller, The NURBS Book, 2nd ed., algorithm A2.3;
%   rational quotient rule, eq. (4.8).

s = s(:);
P = curve.ctrl;
n = size(P, 1);
p = curve.degree;
nd = max(nargout - 1, 0);
B = basis_dense(curve.knots(:)', p, n, s, nd);
w = curve.weights(:);
if isempty(w)
    C = B{1} * P;
    if nd >= 1, Cs = B{2} * P; end
    if nd >= 2, Css = B{3} * P; end
    return
end
Pw = P .* w;
W0 = B{1} * w;
C = (B{1} * Pw) ./ W0;
if nd >= 1
    W1 = B{2} * w;
    Cs = (B{2} * Pw - C .* W1) ./ W0;
end
if nd >= 2
    W2 = B{3} * w;
    Css = (B{3} * Pw - 2 * Cs .* W1 - C .* W2) ./ W0;
end
end

function B = basis_dense(t, p, n, s, nd)
% B{k+1}: [m x n] matrix of the k-th derivatives of the n basis functions at the parameters s.
m = numel(s);
if n > p + 1
    span = sum(s >= t(p + 2:n), 2) + p + 1;
else
    span = repmat(p + 1, m, 1);
end
% ndu(:, a, b): Piegl & Tiller A2.3 table ndu[a-1][b-1], one row per parameter
ndu = zeros(m, p + 1, p + 1);
ndu(:, 1, 1) = 1;
left = zeros(m, p + 1);
right = zeros(m, p + 1);
for j = 1:p
    left(:, j + 1) = s - t(span + 1 - j)';
    right(:, j + 1) = t(span + j)' - s;
    saved = zeros(m, 1);
    for r = 0:j - 1
        ndu(:, j + 1, r + 1) = right(:, r + 2) + left(:, j - r + 1);
        tmp = ndu(:, r + 1, j) ./ ndu(:, j + 1, r + 1);
        ndu(:, r + 1, j + 1) = saved + right(:, r + 2) .* tmp;
        saved = left(:, j - r + 1) .* tmp;
    end
    ndu(:, j + 1, j + 1) = saved;
end
ders = zeros(m, nd + 1, p + 1);
for r = 0:p
    ders(:, 1, r + 1) = ndu(:, r + 1, p + 1);
end
nk = min(nd, p);
for r = 0:p
    a = zeros(m, 2, p + 1);
    s1 = 1;
    s2 = 2;
    a(:, 1, 1) = 1;
    for k = 1:nk
        d = zeros(m, 1);
        rk = r - k;
        pk = p - k;
        if r >= k
            a(:, s2, 1) = a(:, s1, 1) ./ ndu(:, pk + 2, rk + 1);
            d = a(:, s2, 1) .* ndu(:, rk + 1, pk + 1);
        end
        if rk >= -1
            j1 = 1;
        else
            j1 = -rk;
        end
        if r - 1 <= pk
            j2 = k - 1;
        else
            j2 = p - r;
        end
        for j = j1:j2
            a(:, s2, j + 1) = (a(:, s1, j + 1) - a(:, s1, j)) ./ ndu(:, pk + 2, rk + j + 1);
            d = d + a(:, s2, j + 1) .* ndu(:, rk + j + 1, pk + 1);
        end
        if r <= pk
            a(:, s2, k + 1) = -a(:, s1, k) ./ ndu(:, pk + 2, r + 1);
            d = d + a(:, s2, k + 1) .* ndu(:, r + 1, pk + 1);
        end
        ders(:, k + 1, r + 1) = d;
        tmp = s1;
        s1 = s2;
        s2 = tmp;
    end
end
fac = p;
for k = 1:nk
    ders(:, k + 1, :) = ders(:, k + 1, :) * fac;
    fac = fac * (p - k);
end
B = cell(1, nd + 1);
rows = (1:m)';
for k = 0:nd
    M = zeros(m, n);
    for r = 0:p
        M(sub2ind([m n], rows, span - p + r)) = ders(:, k + 1, r + 1);
    end
    B{k + 1} = M;
end
end
