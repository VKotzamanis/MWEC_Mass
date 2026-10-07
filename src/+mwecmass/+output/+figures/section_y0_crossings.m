function x = section_y0_crossings(loop, n_s)
%SECTION_Y0_CROSSINGS x of the points where a closed section loop crosses y = 0.
%
%   x = mwecmass.output.figures.section_y0_crossings(loop, n_s)
%
%   loop: S5 section loop (pieces(k).curve, the exact curve of each piece). The crossings are
%   found on the exact curves: each piece is sampled at n_s parameters (default 33) and a change
%   of [y > 0] between neighbouring samples is solved for y = 0 on the curve by the Illinois form
%   of regula falsi, to the last floating-point step. The half-open rule [y > 0] counts a point
%   exactly on y = 0 once, whichever side the loop leaves it to, and a touch twice, so the number
%   of crossings of a closed loop is even. End samples are the end control points (clamped
%   curves), so two pieces that share a vertex see the same y there. A piece whose control
%   polygon lies on one side of y = 0 cannot cross it (convex hull of positive weights) and is
%   not sampled. Returns the sorted x of the crossings as a column vector.

if nargin < 2 || isempty(n_s)
    n_s = 33;
end
x = zeros(0, 1);
for k = 1:numel(loop.pieces)
    c = loop.pieces(k).curve;
    w = c.weights;
    positive_weights = isempty(w) || all(w > 0);
    y = c.ctrl(:, 2);
    if positive_weights && (all(y > 0) || all(y <= 0))
        continue
    end
    s = linspace(c.knots(1), c.knots(end), n_s)';
    p = point_at(c, s);
    above = p(:, 2) > 0;
    for j = find(above(1:end - 1) ~= above(2:end))'
        x(end + 1, 1) = crossing_x(c, s(j), s(j + 1)); %#ok<AGROW>
    end
end
if mod(numel(x), 2) ~= 0
    error('mwecmass:figures:SectionTopology', ...
        'section_y0_crossings: %d crossings of y = 0 on a closed loop (an even number is required)', numel(x));
end
x = sort(x);
end

function xc = crossing_x(c, s_a, s_b)
% x of the root of y(s) between two samples with different [y > 0]: the Illinois form of regula
% falsi on the exact curve, ended when y is exactly 0, the iterate stops moving, or the bracket
% collapses to adjacent floats. A root at an end sample (y = 0 exactly) is returned without iterating.
p_a = point_at(c, s_a);
p_b = point_at(c, s_b);
if p_b(2) == 0
    xc = p_b(1);
    return
end
if p_a(2) == 0
    xc = p_a(1);
    return
end
ya = p_a(2);
yb = p_b(2);
side = 0;
s_old = NaN;
xc = p_b(1);
for it = 1:200
    s = (s_a * yb - s_b * ya) / (yb - ya);
    if ~(s > min(s_a, s_b) && s < max(s_a, s_b))
        s = s_a + (s_b - s_a) / 2;
    end
    p = point_at(c, s);
    xc = p(1);
    if p(2) == 0 || s == s_old
        return
    end
    s_old = s;
    if (p(2) > 0) == (ya > 0)
        s_a = s;
        ya = p(2);
        if side == 1
            yb = yb / 2;
        end
        side = 1;
    else
        s_b = s;
        yb = p(2);
        if side == -1
            ya = ya / 2;
        end
        side = -1;
    end
    if s_a + (s_b - s_a) / 2 == s_a || s_a + (s_b - s_a) / 2 == s_b
        return
    end
end
end

function p = point_at(c, s)
% Curve points; the end parameters return the end control points exactly.
p = mwecmass.solid.eval_bspline_curve(c, s(:));
at_first = s(:) == c.knots(1);
at_last = s(:) == c.knots(end);
p(at_first, :) = repmat(c.ctrl(1, :), nnz(at_first), 1);
p(at_last, :) = repmat(c.ctrl(end, :), nnz(at_last), 1);
end
