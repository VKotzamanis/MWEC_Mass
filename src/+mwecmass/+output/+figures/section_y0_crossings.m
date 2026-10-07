function x = section_y0_crossings(loop)
%SECTION_Y0_CROSSINGS x of the points where a closed section loop crosses y = 0.
%
%   x = mwecmass.output.figures.section_y0_crossings(loop)
%
%   loop: S5 section loop (pieces(k).curve, the exact curve of each piece; clamped, positive
%   weights). The crossings are found on the exact curves, by bracketing that cannot miss a root:
%   each piece is cut into its Bezier spans (knot insertion on the homogeneous control points
%   [x*w, y*w, w]), and every span is subdivided by de Casteljau until its control values y*w are
%   either all > 0, all <= 0 (by the convex hull property of positive weights no transition of
%   [y > 0] lies inside), or change sign once and monotonically (by the variation diminishing
%   property exactly one root lies inside). That root is solved by the Illinois form of regula
%   falsi on the span, to the last floating-point step. A span whose parameter interval shrinks to
%   adjacent floats without resolving is a touch: it counts one crossing when [y > 0] differs at
%   its two ends and none otherwise (no floating-point point of the span is on the other side).
%   The half-open rule [y > 0] counts a point exactly on y = 0 once, whichever side the loop
%   leaves it to, and a touch of a vertex from above twice, so a closed loop has an even number of
%   crossings; end points are the end control points (clamped curves), so two pieces that share a
%   vertex see the same sign of y there.
%
%   Raises mwecmass:figures:SectionTopology when a weight is not positive, a curve is not clamped,
%   the subdivision budget is exhausted, or the loop does not have an even number of crossings.
%   Returns the sorted x of the crossings as a column vector.

max_spans = 1e5;
x = zeros(0, 1);
for k = 1:numel(loop.pieces)
    c = loop.pieces(k).curve;
    w = c.weights(:);
    if isempty(w)
        w = ones(size(c.ctrl, 1), 1);
    end
    if ~all(w > 0)
        error('mwecmass:figures:SectionTopology', ...
            'section_y0_crossings: piece %d has a weight that is not positive', k);
    end
    spans = bezier_spans([c.ctrl(:, 1:2) .* w, w], c.knots(:)', c.degree, k);
    budget = max_spans;
    for j = 1:numel(spans)
        stack = {spans{j}};
        ta = 0;
        tb = 1;
        while ~isempty(stack)
            budget = budget - 1;
            if budget < 0
                error('mwecmass:figures:SectionTopology', ...
                    'section_y0_crossings: subdivision budget exhausted on piece %d', k);
            end
            P = stack{end};
            a = ta(end);
            b = tb(end);
            stack(end) = [];
            ta(end) = [];
            tb(end) = [];
            above = P(:, 2) > 0;
            n_change = nnz(above(1:end - 1) ~= above(2:end));
            if n_change == 0
                continue
            end
            d = diff(P(:, 2));
            mid = a + (b - a) / 2;
            if n_change == 1 && (all(d >= 0) || all(d <= 0))
                x(end + 1, 1) = crossing_x(P); %#ok<AGROW>
            elseif mid <= a || mid >= b
                if above(1) ~= above(end)
                    x(end + 1, 1) = crossing_x(P); %#ok<AGROW>
                end
            else
                [left, right] = split_half(P);
                stack(end + 1:end + 2) = {right, left};
                ta(end + 1:end + 2) = [mid, a];
                tb(end + 1:end + 2) = [b, mid];
            end
        end
    end
end
if mod(numel(x), 2) ~= 0
    error('mwecmass:figures:SectionTopology', ...
        'section_y0_crossings: %d crossings of y = 0 on a closed loop (an even number is required)', numel(x));
end
x = sort(x);
end

function spans = bezier_spans(Pw, U, p, k)
% Bezier control polygons (homogeneous rows [x*w y*w w]) of the non-empty knot spans of a clamped
% B-spline: every interior knot is raised to multiplicity p by Boehm knot insertion.
if numel(U) ~= size(Pw, 1) + p + 1 || any(U(1:p + 1) ~= U(1)) || any(U(end - p:end) ~= U(end))
    error('mwecmass:figures:SectionTopology', 'section_y0_crossings: piece %d is not a clamped curve', k);
end
interior = U(p + 2:end - p - 1);
for u = unique(interior)
    s = nnz(U == u);
    for r = 1:p - s
        [Pw, U] = insert_knot(Pw, U, p, u, nnz(U == u));
    end
end
breaks = unique(U);
spans = cell(1, numel(breaks) - 1);
for j = 1:numel(spans)
    first = (j - 1) * p + 1;
    spans{j} = Pw(first:first + p, :);
end
end

function [Q, V] = insert_knot(P, U, p, u, s)
% One Boehm insertion of u (current multiplicity s < p) into the homogeneous polygon P.
kk = find(U <= u, 1, 'last');
n = size(P, 1);
Q = zeros(n + 1, size(P, 2));
Q(1:kk - p, :) = P(1:kk - p, :);
for i = kk - p + 1:kk - s
    alpha = (u - U(i)) / (U(i + p) - U(i));
    Q(i, :) = alpha * P(i, :) + (1 - alpha) * P(i - 1, :);
end
Q(kk - s + 1:n + 1, :) = P(kk - s:n, :);
V = [U(1:kk), u, U(kk + 1:end)];
end

function [L, R] = split_half(P)
% de Casteljau split at t = 1/2.
m = size(P, 1);
L = zeros(m, size(P, 2));
R = L;
Q = P;
for i = 1:m
    L(i, :) = Q(1, :);
    R(m - i + 1, :) = Q(end, :);
    Q = (Q(1:end - 1, :) + Q(2:end, :)) / 2;
end
end

function q = bez_point(P, t)
% de Casteljau evaluation of a homogeneous polygon; t = 0 and 1 return the end rows exactly.
Q = P;
while size(Q, 1) > 1
    Q = (1 - t) * Q(1:end - 1, :) + t * Q(2:end, :);
end
q = Q;
end

function xc = crossing_x(P)
% x of the root of y between the ends of a span whose end values have different [y > 0]: the
% Illinois form of regula falsi on the span, ended when y is exactly 0 or the bracket collapses to
% adjacent floats. A root at an end (y = 0 exactly) is returned without iterating.
if P(end, 2) == 0
    xc = P(end, 1) / P(end, 3);
    return
end
if P(1, 2) == 0
    xc = P(1, 1) / P(1, 3);
    return
end
ta = 0;
tb = 1;
ya = P(1, 2);
yb = P(end, 2);
side = 0;
for it = 1:200
    t = (ta * yb - tb * ya) / (yb - ya);
    if ~(t > ta && t < tb)
        t = ta + (tb - ta) / 2;
    end
    q = bez_point(P, t);
    xc = q(1) / q(3);
    if q(2) == 0
        return
    end
    if (q(2) > 0) == (ya > 0)
        ta = t;
        ya = q(2);
        if side == 1
            yb = yb / 2;
        end
        side = 1;
    else
        tb = t;
        yb = q(2);
        if side == -1
            ya = ya / 2;
        end
        side = -1;
    end
    if ta + (tb - ta) / 2 == ta || ta + (tb - ta) / 2 == tb
        return
    end
end
end
