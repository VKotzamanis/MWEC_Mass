function [c, rep] = fit_bspline_surface(prob, opts)
%FIT_BSPLINE_SURFACE  Adaptive, error-bounded least-squares B-spline fit with judged metrics.
%
%   [c, rep] = mwecmass.solid.fit_bspline_surface(prob, opts)
%   n = mwecmass.solid.fit_bspline_surface('max_passes')   default refinement cap
%
%   Fits the curve that an inner surface is built from (the offset profile of a rev_z piece; the
%   surface follows from it by the outer surface's own v-construction), after Piegl & Tiller, The
%   NURBS Book, 2nd ed., sections 9.4 (least squares
%   with end constraints) and 5.4 (knot removal):
%     1. fit a non-rational B-spline of degree prob.degree through nodes prob.sample(s) placed at
%        every knot and opts.nodes_per_span points inside every span;
%     2. judge it with prob.judge(c), which returns one pass flag per knot span (the metrics are
%        evaluated between the faces as written, at check points that are not nodes);
%     3. insert a knot at the middle of every failing span and refit, up to opts.max_passes;
%     4. remove every interior knot whose removal keeps all spans passing.
%   prob fields: degree; knots (initial interior knots, single, in (0, 1)); breaks (interior
%   knots of full multiplicity, never removed); sample (handle s -> [numel(s) x m] points, s in
%   [0, 1]); fix0, fix1 ([1 x m] end values, NaN where free); tie0, tie1 (logical [1 x m]: the
%   second, second-to-last, control value equals the first, last, in that coordinate, a tangent
%   in the coordinate plane); monotone (coordinate whose control values are kept monotone between
%   the end values by pooling adjacent violators; 0 for none); judge (handle c -> [span_ok, info]).
%   A coordinate whose node values are all equal, and equal to its fixed end values, is reproduced
%   exactly (partition of unity); otherwise the fixed end values take precedence.
%   opts: max_passes (default: the value returned by the 'max_passes' query, a limit set at the
%   T3 checkpoint, contract section 9 item 2), knots_fixed (true: fit once on prob.knots and
%   prob.breaks, no insertion or removal), nodes_per_span (default 8).
%   rep: n_nodes, n_knots, n_passes, n_removed, cap_reached, info (judge output of c).

if ischar(prob)
    if strcmp(prob, 'max_passes')
        c = 12;
        return
    end
    error('mwecmass:solid:BadQuery', 'fit_bspline_surface: unknown query %s', prob);
end
if nargin < 2
    opts = struct();
end
max_passes = get_opt(opts, 'max_passes', mwecmass.solid.fit_bspline_surface('max_passes'));
fixed = get_opt(opts, 'knots_fixed', false);
nps = get_opt(opts, 'nodes_per_span', 8);
p = prob.degree;
kn = sort(prob.knots(:)');
br = sort(prob.breaks(:)');

[c, nn] = fit_once(prob, kn, br, p, nps);
[ok, info] = prob.judge(c);
passes = 1;
cap = false;
if ~fixed
    while ~all(ok)
        if passes >= max_passes
            cap = true;
            break
        end
        ku = unique(c.knots);
        bad = find(~ok(:)');
        kn = sort([kn, (ku(bad) + ku(bad + 1)) / 2]);
        [c, nn] = fit_once(prob, kn, br, p, nps);
        [ok, info] = prob.judge(c);
        passes = passes + 1;
    end
end
removed = 0;
if ~fixed && ~cap
    i = 1;
    while i <= numel(kn)
        trial = kn([1:i - 1, i + 1:end]);
        [ct, nt] = fit_once(prob, trial, br, p, nps);
        [okt, it] = prob.judge(ct);
        if all(okt)
            kn = trial;
            c = ct;
            nn = nt;
            info = it;
            removed = removed + 1;
        else
            i = i + 1;
        end
    end
end
rep = struct('n_nodes', nn, 'n_knots', numel(c.knots), 'n_passes', passes, 'n_removed', removed, ...
    'cap_reached', cap, 'info', info);
end

function v = get_opt(opts, name, default)
v = default;
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
end
end

function [c, nn] = fit_once(prob, kn, br, p, nps)
knots = [zeros(1, p + 1), sort([kn, reshape(repmat(br, p, 1), 1, [])]), ones(1, p + 1)];
ku = unique(knots);
s = ku(:);
for j = 1:numel(ku) - 1
    s = [s; ku(j) + ((1:nps)' - 0.5) / nps * (ku(j + 1) - ku(j))]; %#ok<AGROW>
end
s = sort(s);
Q = prob.sample(s);
nn = numel(s);
n = numel(knots) - p - 1;
N = mwecmass.solid.eval_bspline_curve(struct('degree', p, 'ctrl', eye(n), 'knots', knots, 'weights', []), s);
m = size(Q, 2);
P = zeros(n, m);
for col = 1:m
    val = NaN(n, 1);
    if isfield(prob, 'fix0') && ~isnan(prob.fix0(col))
        val(1) = prob.fix0(col);
    end
    if isfield(prob, 'fix1') && ~isnan(prob.fix1(col))
        val(n) = prob.fix1(col);
    end
    % equal node values are reproduced exactly when the fixed ends agree with them; a fixed end
    % (a point shared with a neighbour) takes precedence
    if all(Q(:, col) == Q(1, col)) && all(val(~isnan(val)) == Q(1, col))
        P(:, col) = Q(1, col);
        continue
    end
    tie2 = isfield(prob, 'tie0') && prob.tie0(col) && ~isnan(val(1));
    tie3 = isfield(prob, 'tie1') && prob.tie1(col) && ~isnan(val(n));
    if tie2
        val(2) = val(1);
    end
    if tie3
        val(n - 1) = val(n);
    end
    free = isnan(val);
    rhs = Q(:, col) - N(:, ~free) * val(~free);
    val(free) = N(:, free) \ rhs;
    if isfield(prob, 'monotone') && prob.monotone == col
        val = monotone(val, ~free);
    end
    P(:, col) = val;
end
c = struct('degree', p, 'ctrl', P, 'knots', knots, 'weights', []);
end

function v = monotone(v, keep)
% pool adjacent violators toward the direction of the end values; values in keep stay
dirn = sign(v(end) - v(1));
if dirn == 0
    return
end
w = dirn * v;
blocks = num2cell(w(:)');
idx = num2cell(1:numel(w));
k = 1;
while k < numel(blocks)
    if mean(blocks{k}) > mean(blocks{k + 1})
        blocks{k} = [blocks{k}, blocks{k + 1}];
        idx{k} = [idx{k}, idx{k + 1}];
        blocks(k + 1) = [];
        idx(k + 1) = [];
        k = max(k - 1, 1);
    else
        k = k + 1;
    end
end
out = w;
for k = 1:numel(blocks)
    out(idx{k}) = mean(blocks{k});
end
lo = min(w([1 end]));
hi = max(w([1 end]));
out = min(max(out, lo), hi);
out(keep) = w(keep);
v = dirn * out;
end
