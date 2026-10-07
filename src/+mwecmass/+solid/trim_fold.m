function keep = trim_fold(F, ranges, n_samples)
%TRIM_FOLD  Trim an ordered chain of plane offset curves where it overlaps itself.
%
%   keep = mwecmass.solid.trim_fold(F, ranges)
%   keep = mwecmass.solid.trim_fold(F, ranges, n_samples)
%
%   F: cell of function handles, F{k}(s) -> [Q, dQ] with Q, dQ [numel(s) x 2]: the exact offset
%   curve k of the chain (a normal offset of a profile segment, or a fan arc closing a concave
%   crease) and its derivative, on the parameter interval ranges(k, :). The curves follow each
%   other along the chain; consecutive curves either join (end of k = start of k+1) or leave a gap
%   that crosses itself (a convex crease, where the two offsets overlap).
%   Where the chain crosses itself, the part between the two crossing points is a fold or the
%   overlap of a convex crease and is removed (the trimmed points lie closer than the offset
%   distance to the profile). Crossings are found on polylines of n_samples points per curve
%   (default 400) and refined by Newton's method on the exact curves to adjacent doubles.
%   keep: struct array, in chain order, of the kept parts: k (curve index), s0, s1 (parameter
%   interval), X0, X1 [1x2] (end points; the point of a crossing is computed once and given to
%   both parts that meet there), trim0, trim1 (true where the part ends at a crossing).

if nargin < 3 || isempty(n_samples)
    n_samples = 400;
end
K = numel(F);
% polyline of every curve
S = cell(K, 1);
Q = cell(K, 1);
for k = 1:K
    S{k} = linspace(ranges(k, 1), ranges(k, 2), n_samples)';
    Q{k} = F{k}(S{k});
end
% chain positions: curve k, segment j (between samples j and j+1)
segs = zeros(0, 2);
A = zeros(0, 2);
B = zeros(0, 2);
for k = 1:K
    n = size(Q{k}, 1);
    segs = [segs; [repmat(k, n - 1, 1), (1:n - 1)']]; %#ok<AGROW>
    A = [A; Q{k}(1:n - 1, :)]; %#ok<AGROW>
    B = [B; Q{k}(2:n, :)]; %#ok<AGROW>
end
m = size(segs, 1);
hits = zeros(0, 2);
for i = 1:m - 1
    j = (i + 1:m)';
    % adjacent segments of one curve share a point; consecutive curves that join share one too
    same_curve_next = segs(j, 1) == segs(i, 1) & segs(j, 2) == segs(i, 2) + 1;
    joint = segs(j, 1) == segs(i, 1) + 1 & segs(j, 2) == 1 & segs(i, 2) == size(Q{segs(i, 1)}, 1) - 1 & ...
        all(B(i, :) == A(j, :), 2);
    j = j(~same_curve_next & ~joint);
    if isempty(j)
        continue
    end
    Ai = repmat(A(i, :), numel(j), 1);
    Bi = repmat(B(i, :), numel(j), 1);
    d1 = orient(A(j, :), B(j, :), Ai);
    d2 = orient(A(j, :), B(j, :), Bi);
    d3 = orient(Ai, Bi, A(j, :));
    d4 = orient(Ai, Bi, B(j, :));
    x = j(d1 .* d2 < 0 & d3 .* d4 < 0);
    hits = [hits; [repmat(i, numel(x), 1), x]]; %#ok<AGROW>
end
% outermost loops first: a crossing inside a removed loop is gone with it
cuts = zeros(0, 2);
pos = 1;
for q = 1:size(hits, 1)
    if hits(q, 1) < pos
        continue
    end
    c = hits(q, :);
    inner = hits(hits(:, 1) == c(1), :);
    c = inner(end, :);
    cuts(end + 1, :) = c; %#ok<AGROW>
    pos = c(2) + 1;
end
keep = struct('k', {}, 's0', {}, 's1', {}, 'X0', {}, 'X1', {}, 'trim0', {}, 'trim1', {});
start = struct('k', 1, 's', ranges(1, 1), 'X', [], 'trim', false);
done = false(K, 1);
for q = 1:size(cuts, 1)
    i = cuts(q, 1);
    j = cuts(q, 2);
    [si, sj, X] = refine(F, segs(i, :), segs(j, :), S, Q);
    ki = segs(i, 1);
    kj = segs(j, 1);
    for k = start.k:ki
        p = struct('k', k, 's0', ranges(k, 1), 's1', ranges(k, 2), 'X0', [], 'X1', [], 'trim0', false, 'trim1', false);
        if k == start.k
            p.s0 = start.s;
            p.X0 = start.X;
            p.trim0 = start.trim;
        end
        if k == ki
            p.s1 = si;
            p.X1 = X;
            p.trim1 = true;
        end
        keep(end + 1) = p; %#ok<AGROW>
        done(k) = true;
    end
    start = struct('k', kj, 's', sj, 'X', X, 'trim', true);
end
for k = start.k:K
    p = struct('k', k, 's0', ranges(k, 1), 's1', ranges(k, 2), 'X0', [], 'X1', [], 'trim0', false, 'trim1', false);
    if k == start.k
        p.s0 = start.s;
        p.X0 = start.X;
        p.trim0 = start.trim;
    end
    keep(end + 1) = p; %#ok<AGROW>
end
for q = 1:numel(keep)
    if isempty(keep(q).X0)
        keep(q).X0 = F{keep(q).k}(keep(q).s0);
    end
    if isempty(keep(q).X1)
        keep(q).X1 = F{keep(q).k}(keep(q).s1);
    end
end
end

function [si, sj, X] = refine(F, a, b, S, Q)
% Newton on F{ka}(si) = F{kb}(sj), started at the polyline crossing
ka = a(1);
kb = b(1);
P1 = Q{ka}(a(2), :);
P2 = Q{ka}(a(2) + 1, :);
R1 = Q{kb}(b(2), :);
R2 = Q{kb}(b(2) + 1, :);
M = [P2 - P1; -(R2 - R1)]';
w = M \ (R1 - P1)';
si = S{ka}(a(2)) + w(1) * (S{ka}(a(2) + 1) - S{ka}(a(2)));
sj = S{kb}(b(2)) + w(2) * (S{kb}(b(2) + 1) - S{kb}(b(2)));
for it = 1:60
    [Pa, Da] = F{ka}(si);
    [Pb, Db] = F{kb}(sj);
    r = (Pa - Pb)';
    J = [Da', -Db'];
    step = J \ r;
    si_new = si - step(1);
    sj_new = sj - step(2);
    if si_new == si && sj_new == sj
        break
    end
    si = si_new;
    sj = sj_new;
end
Pa = F{ka}(si);
Pb = F{kb}(sj);
X = (Pa + Pb) / 2;
end

function d = orient(a, b, c)
d = (b(:, 1) - a(:, 1)) .* (c(:, 2) - a(:, 2)) - (b(:, 2) - a(:, 2)) .* (c(:, 1) - a(:, 1));
end
