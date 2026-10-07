function [inner, rep] = offset_surface(model, cache, geo, t, z_range, opts)
%OFFSET_SURFACE  Inner surface of a shell of design thickness t: the normal offset of the outer surface (contract F2, S2, S2r).
%
%   [inner, rep] = mwecmass.solid.offset_surface(model, cache, geo, t, z_range, opts)
%
%   geo: S1b of outer_nurbs; t [m]; z_range [m, body]: the inner surface must cover these heights;
%   opts.t_min [m] (required): eps_fit = 0.01 t_min, and the surface is built at the normal
%   distance d = t + eps_fit/2 (AGENTS section 5 item 9.2). opts.max_passes: refinement cap
%   (default fit_bspline_surface('max_passes')); opts.n_gauss: Gauss points per knot span of the
%   check points (default 6); opts.knots_from: an S2 set whose pieces and knot vectors are reused,
%   only the control points refitted (no refinement, M1-M3 reported).
%   d >= d_close (void_closing_distance over z_range) raises mwecmass:solid:VoidClosed.
%
%   Construction, per visible source patch (mirrors are flipped copies, as in F1):
%   - rev_z: the profile in the meridian plane of the first column is split at its C0 knots; each
%     segment is offset by d along its inward normal; at a convex crease or a fold the offsets
%     overlap and are trimmed at their crossing (trim_fold); a concave crease gets a face of its
%     own, the crease point offset by d along the fan of normals from one side's to the other's
%     (a circular arc of radius d, exact). Straight segments are translated exactly; curved ones
%     are fitted (fit_bspline_surface) with end points shared with their neighbours, an axis end
%     kept on the axis with a horizontal tangent. Each smooth run between creases becomes one
%     patch: the profile revolved with the outer surface's own rational arc (rows of one z).
%   - ruled_parallel, planar (every boundary a convex crease with a planar neighbour): the face
%     moved by d along its inward normal and trimmed by its neighbours' offset planes (corners at
%     the crossing of three planes).
%   - ruled_parallel with equal rulings and a tangent-continuous seam with a rev_z patch: ruled
%     between that patch's inner boundary curve (taken bitwise, so the seam is one curve) and its
%     translate by the ruling.
%   Every other patch takes the general path (mwecmass.solid.fit_z_faces).
%   Only the outer part whose points lie within d of z_range is offset, so a set over a partial
%   range is open at its clipped end (F5 cuts it there).
%   M1-M3 are judged between the faces as written (inner as fitted, outer as in geo) at the check
%   points (Gauss points and the midpoint of every knot span): t_local = distance to the outer
%   surface; M1: t_local >= t_min; M2: t <= t_local <= t + eps_fit; M3: rows of one z with
%   monotone heights, every non-pole boundary shared with one neighbour (same curve), and every
%   horizontal section of the inner set one simple closed loop inside the outer section.
%   A fitted profile that reaches the refinement cap raises mwecmass:solid:FitNotConverged.
%   Results are cached per (hull, t bitwise, z_range, t_min, options).

if nargin < 6 || isempty(opts)
    opts = struct();
end
if ~isfield(opts, 't_min') || isempty(opts.t_min)
    error('mwecmass:solid:FitInputMissing', 'offset_surface: opts.t_min is required');
end
persistent store
if isempty(store)
    store = containers.Map();
end
key = cache_key(model, geo, t, z_range, opts);
if store.isKey(key)
    v = store(key);
    inner = v{1};
    rep = v{2};
    return
end

eps_fit = 0.01 * opts.t_min;
d = t + eps_fit / 2;
zr = sort(z_range(:)');
d_close = mwecmass.solid.void_closing_distance(model, cache, geo, zr);
if d >= d_close
    error('mwecmass:solid:VoidClosed', 'offset_surface: d = %.17g m >= d_close = %.17g m (t_max = %.17g m)', ...
        d, d_close, d_close - eps_fit / 2);
end
C = struct('d', d, 't', t, 't_min', opts.t_min, 'eps_fit', eps_fit, 'zr', zr, ...
    'max_passes', get_opt(opts, 'max_passes', mwecmass.solid.fit_bspline_surface('max_passes')), ...
    'n_gauss', get_opt(opts, 'n_gauss', 6), 'knots_from', get_opt(opts, 'knots_from', []), ...
    'outer', geo.outer);

P = geo.outer;
nv = max([P.visible]);
prim = primaries(P, nv);
built = cell(1, nv);
fitrep = cell(1, nv);
general_vis = false(1, nv);
% rev_z patches first (ruled neighbours take their boundary curves)
for v = 1:nv
    G = P([P.visible] == v);
    if isempty(G) || prim(v) ~= v || ~strcmp(G(1).offset_kind, 'rev_z')
        continue
    end
    [E, fr, ok] = rev_inner(G, C);
    if ok
        built{v} = E;
        fitrep{v} = fr;
    else
        general_vis(v) = true;
    end
end
built = mirror_all(P, prim, built, nv);
for v = 1:nv
    G = P([P.visible] == v);
    if isempty(G) || prim(v) ~= v || ~isempty(built{v}) || general_vis(v)
        continue
    end
    E = [];
    if strcmp(G(1).offset_kind, 'ruled_parallel')
        E = planar_inner(G, P, C);
        if isempty(E)
            E = ruled_inner(G, P, built);
        end
    end
    if isempty(E)
        general_vis(v) = true;
    else
        built{v} = E;
    end
end
built = mirror_all(P, prim, built, nv);
for v = 1:nv
    if isempty(built{v}) && ~general_vis(prim(v))
        general_vis(v) = true;
    end
end

patches = [];
for v = 1:nv
    patches = [patches, built{v}]; %#ok<AGROW>
end
flat = struct('z', {}, 'normal_z', {}, 'visible', {});
if any(general_vis)
    general = [false(1, numel(patches)), true(1, nnz(general_vis))];
    stubs = P([]);
    for v = find(general_vis)
        G = P([P.visible] == v);
        s = G(1);
        s.surf = [];
        s.exact = false;
        stubs(end + 1) = s; %#ok<AGROW>
    end
    fopts = struct('t_min', opts.t_min, 'max_passes', C.max_passes, 'kind', 'inner', 'd', d, 'geo', geo);
    [patches, flat] = mwecmass.solid.fit_z_faces(model, cache, [patches, stubs], general, fopts);
end
patches = find_seams(patches);
rep = judge_set(patches, C, fitrep, prim);
zz = reshape([patches.z_range], 2, [])';
inner = struct('t', t, 'd', d, 'eps_fit', eps_fit, 'z_range', [min(zz(:)) max(zz(:))], ...
    'z_lo', min(zz(:)), 'refit', ~isempty(C.knots_from), 'patches', patches, 'flat', flat, 'report', rep);
if rep.cap_reached
    bad = find(~([rep.patches.M1] & [rep.patches.M2]));
    error('mwecmass:solid:FitNotConverged', 'offset_surface: refinement cap %d reached on %s (t_local %s m)', ...
        C.max_passes, strjoin({patches(bad).name}, ', '), mat2str([[rep.patches(bad).t_local_min]; ...
        [rep.patches(bad).t_local_max]], 6));
end
store(key) = {inner, rep};
end

% =============================================================== helpers

function v = get_opt(opts, name, default)
v = default;
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
end
end

function key = cache_key(model, geo, t, z_range, opts)
% the hull (deck and a fingerprint of the control nets), t bitwise and the options
f = 0;
g = 0;
for k = 1:numel(geo.outer)
    s = geo.outer(k).surf;
    if isempty(s)
        continue
    end
    x = s.ctrl(:);
    f = f + sum(x .* (1:numel(x))') + k * sum(s.knots{1});
    g = g + sum(x.^2) + numel(x);
end
kf = '';
if isfield(opts, 'knots_from') && ~isempty(opts.knots_from)
    kf = [num2hex(opts.knots_from.t), sprintf('%d', numel(opts.knots_from.patches))];
end
o = [get_opt(opts, 't_min', NaN), get_opt(opts, 'max_passes', NaN), get_opt(opts, 'n_gauss', NaN)];
key = [model.filename, '|', num2hex(f), num2hex(g), '|', num2hex(t), '|', reshape(num2hex(z_range(:))', 1, []), ...
    '|', reshape(num2hex(o(:))', 1, []), '|', kf];
end

function prim = primaries(P, nv)
% the visible patch whose inner set is built for each visible patch: per ultimate source, the one
% with the fewest flips (mirrors are flipped copies of it)
prim = zeros(1, nv);
for v = 1:nv
    G = P([P.visible] == v);
    if isempty(G)
        continue
    end
    cand = unique([P(strcmp({P.source}, G(1).source)).visible]);
    nf = arrayfun(@(c) numel(P(find([P.visible] == c, 1)).flips), cand);
    [~, i] = min(nf);
    prim(v) = cand(i);
end
end

function built = mirror_all(P, prim, built, nv)
for v = 1:nv
    if prim(v) == v || ~isempty(built{v}) || isempty(built{prim(v)})
        continue
    end
    G = P([P.visible] == v);
    F0 = P(find([P.visible] == prim(v), 1)).flips;
    rel = setxor(G(1).flips, F0);
    E = built{prim(v)};
    for k = 1:numel(E)
        for f = 1:numel(rel)
            col = find(strcmp(rel{f}, {'X', 'Y'}), 1);
            E(k).surf.ctrl(:, :, col) = -E(k).surf.ctrl(:, :, col);
            E(k).outward = ~E(k).outward;
        end
        E(k).name = strrep(E(k).name, P(find([P.visible] == prim(v), 1)).name, G(1).name);
        E(k).flips = G(1).flips;
        E(k).visible = v;
    end
    built{v} = E;
end
end

function e = new_entry(src, name, surf, outward, kind)
e = src;
e.name = name;
e.surf = surf;
e.outward = outward;
e.exact = false;
e.z_of_u = true;
e.offset_kind = kind;
e.fit = [];
for f = {'seam_u0', 'seam_u1', 'seam_v0', 'seam_v1'}
    e.(f{1}) = [];
end
e.u_range = [surf.knots{1}(1) surf.knots{1}(end)];
e.z_range = [surf.ctrl(1, 1, 3) surf.ctrl(end, 1, 3)];
e.pole = [collapsed(surf.ctrl(1, :, :)) collapsed(surf.ctrl(end, :, :))];
e.c0_u = full_knots(surf.knots{1}, surf.degree(1));
e.c0_v = full_knots(surf.knots{2}, surf.degree(2));
end

function c = full_knots(k, p)
u = unique(k);
u = u(u > k(1) & u < k(end));
c = u(arrayfun(@(x) sum(k == x), u) >= p);
end

function tf = collapsed(row)
row = reshape(row, [], 3);
tf = all(all(row == row(1, :)));
end

% =============================================================== rev_z

function [E, fr, ok] = rev_inner(G, C)
E = [];
fr = [];
ok = false;
[~, o] = sort(arrayfun(@(e) e.u_range(1), G));
G = G(o);
% axis and the arc pattern of the outer revolution
ax = [];
for k = 1:numel(G)
    s = G(k).surf;
    for r = [1 size(s.ctrl, 1)]
        if collapsed(s.ctrl(r, :, :))
            ax = reshape(s.ctrl(r, 1, 1:2), 1, 2);
        end
    end
end
s = G(1).surf;
W = s.weights;
if isempty(W)
    W = ones(size(s.ctrl, 1), size(s.ctrl, 2));
end
if isempty(ax)
    i = find(any(any(s.ctrl(:, :, 1:2) ~= s.ctrl(:, 1, 1:2), 2), 3), 1);
    if isempty(i) || size(s.ctrl, 2) < 3
        return
    end
    A = reshape(s.ctrl(i, 1, 1:2), 1, 2);
    Cc = reshape(s.ctrl(i, 2, 1:2), 1, 2);
    B = reshape(s.ctrl(i, 3, 1:2), 1, 2);
    M = [Cc - A; Cc - B];
    ax = (M \ [A * (Cc - A)'; B * (Cc - B)'])';
end
R = sqrt((s.ctrl(:, 1, 1) - ax(1)).^2 + (s.ctrl(:, 1, 2) - ax(2)).^2);
iref = find(R > 0 & W(:, 1) == 1, 1);
if isempty(iref)
    iref = find(R > 0, 1);
end
if isempty(iref)
    return
end
U = (reshape(s.ctrl(iref, :, 1:2), [], 2) - ax) / R(iref);
b = W(iref, :) / W(iref, 1);
e1 = U(1, :);

% profile segments in the meridian plane of the first column, split at every C0 knot
segs = struct('ctrl', {}, 'w', {}, 'knots', {}, 'ua', {}, 'ub', {}, 'entry', {}, 'linear', {}, 'flatz', {});
for k = 1:numel(G)
    s = G(k).surf;
    Wk = s.weights;
    if isempty(Wk)
        Wk = ones(size(s.ctrl, 1), size(s.ctrl, 2));
    end
    Pc = reshape(s.ctrl(:, 1, :), [], 3);
    r = (Pc(:, 1) - ax(1)) * e1(1) + (Pc(:, 2) - ax(2)) * e1(2);
    r(Pc(:, 1) == ax(1) & Pc(:, 2) == ax(2)) = 0;
    c2 = struct('degree', s.degree(1), 'ctrl', [r, Pc(:, 3)], 'knots', s.knots{1}, 'weights', Wk(:, 1));
    if all(Wk(:, 1) == 1)
        c2.weights = [];
    end
    cuts = [s.knots{1}(1), full_knots(s.knots{1}, s.degree(1)), s.knots{1}(end)];
    for j = 1:numel(cuts) - 1
        sg = sub_curve(c2, cuts(j), cuts(j + 1));
        segs(end + 1) = struct('ctrl', sg.ctrl, 'w', sg.weights, 'knots', sg.knots, 'ua', cuts(j), ...
            'ub', cuts(j + 1), 'entry', k, 'linear', sg.degree == 1, ...
            'flatz', all(sg.ctrl(:, 2) == sg.ctrl(1, 2))); %#ok<AGROW>
    end
end
for k = 1:numel(segs)
    segs(k).degree = G(segs(k).entry).surf.degree(1);
end

% inward normal sign in (r, z): compare with the 3-D normal at the middle of the first non-degenerate segment
sn = 0;
for k = 1:numel(segs)
    ent = G(segs(k).entry);
    um = (segs(k).ua + segs(k).ub) / 2;
    vm = ent.surf.knots{2}(1);
    [~, Su, Sv] = mwecmass.solid.eval_bspline_surface(ent.surf, um, vm);
    n3 = cross(Su, Sv);
    if norm(n3) == 0
        continue
    end
    if ent.outward
        n3 = -n3;
    end
    [~, T] = seg_eval(segs(k), um);
    n2 = [n3(1:2) * e1', n3(3)];
    sn = sign(rot(T) * n2');
    if sn ~= 0
        break
    end
end
if sn == 0
    return
end

% clip to the heights whose offsets can reach z_range
lo = C.zr(1) - C.d;
hi = C.zr(2) + C.d;
keepseg = true(1, numel(segs));
open0 = false;
open1 = false;
for k = 1:numel(segs)
    z = segs(k).ctrl(:, 2);
    if max(z) < lo || min(z) > hi
        keepseg(k) = false;
    end
end
first = find(keepseg, 1);
last = find(keepseg, 1, 'last');
if first > 1
    open0 = true;
end
if last < numel(segs)
    open1 = true;
end
segs = segs(first:last);
for k = [1 numel(segs)]
    z = segs(k).ctrl(:, 2);
    if ~segs(k).flatz && (min(z) < lo || max(z) > hi)
        bounds = [lo hi];
        bounds = bounds(bounds > min(z) & bounds < max(z));
        for zb = bounds
            ub = seg_root(segs(k), zb);
            zs = seg_eval(segs(k), segs(k).ua);
            inside_a = zs(2) >= lo && zs(2) <= hi;
            if inside_a
                segs(k) = restrict(segs(k), segs(k).ua, ub);
                if k == numel(segs), open1 = true; end
            else
                segs(k) = restrict(segs(k), ub, segs(k).ub);
                if k == 1, open0 = true; end
            end
        end
    end
end

% joints: smooth, convex crease or concave crease
nseg = numel(segs);
joint = cell(1, nseg - 1);
for k = 1:nseg - 1
    a = segs(k).ctrl(end, :) - segs(k).ctrl(end - 1, :);
    bb = segs(k + 1).ctrl(2, :) - segs(k + 1).ctrl(1, :);
    if all(a == 0) || all(bb == 0)
        joint{k} = 'smooth';
        continue
    end
    cr = a(1) * bb(2) - a(2) * bb(1);
    % parallel up to the rounding of the control points (16 ulp of the leg lengths)
    if abs(cr) <= 16 * eps * norm(a) * norm(bb) && a * bb' > 0
        joint{k} = 'smooth';
    else
        nA = sn * rot(a) / norm(a);
        if bb * nA' > 0
            joint{k} = 'convex';
        else
            joint{k} = 'concave';
        end
    end
end

% exact offset curves and fans, in chain order
F = {};
rng = zeros(0, 2);
what = zeros(0, 2);
for k = 1:nseg
    F{end + 1} = @(u) offset_point(segs(k), u, C.d, sn); %#ok<AGROW>
    rng(end + 1, :) = [segs(k).ua segs(k).ub]; %#ok<AGROW>
    what(end + 1, :) = [1 k]; %#ok<AGROW>
    if k < nseg && strcmp(joint{k}, 'concave')
        X = segs(k).ctrl(end, :);
        [~, Ta] = seg_eval(segs(k), segs(k).ub);
        [~, Tb] = seg_eval(segs(k + 1), segs(k + 1).ua);
        na = sn * rot(Ta) / norm(Ta);
        nb = sn * rot(Tb) / norm(Tb);
        F{end + 1} = @(th) fan_point(X, na, nb, C.d, th); %#ok<AGROW>
        rng(end + 1, :) = [0 1]; %#ok<AGROW>
        what(end + 1, :) = [2 k]; %#ok<AGROW>
    end
end
keep = mwecmass.solid.trim_fold(F, rng);
keep = clip_axis(keep, F);
% a crossing on the offset of a straight segment with a constant coordinate (a horizontal or a
% vertical line in the meridian plane) takes that coordinate of the offset line exactly
for q = 1:numel(keep) - 1
    if ~keep(q).trim1
        continue
    end
    X = keep(q).X1;
    for w = [keep(q).k, keep(q + 1).k]
        if what(w, 1) ~= 1 || ~segs(what(w, 2)).linear
            continue
        end
        sg = segs(what(w, 2));
        [~, T] = seg_eval(sg, sg.ua);
        n = sn * rot(T) / norm(T);
        for c = 1:2
            if all(sg.ctrl(:, c) == sg.ctrl(1, c))
                X(c) = sg.ctrl(1, c) + C.d * n(c);
            end
        end
    end
    keep(q).X1 = X;
    keep(q + 1).X0 = X;
end

% runs: maximal sequences of kept parts joined smoothly without a trim
runs = {};
cur = 1;
for q = 2:numel(keep)
    a = keep(q - 1);
    bq = keep(q);
    wa = what(a.k, :);
    wb = what(bq.k, :);
    smooth = wa(1) == 1 && wb(1) == 1 && wb(2) == wa(2) + 1 && strcmp(joint{wa(2)}, 'smooth') && ...
        ~a.trim1 && ~bq.trim0 && segs(wa(2)).flatz == segs(wb(2)).flatz;
    if smooth
        cur(end + 1) = q; %#ok<AGROW>
    else
        runs{end + 1} = cur; %#ok<AGROW>
        cur = q;
    end
end
runs{end + 1} = cur;

fr = struct('name', {}, 'n_nodes', {}, 'n_knots', {}, 'n_passes', {}, 'n_removed', {}, 'cap_reached', {});
names = G(1).name;
kf = C.knots_from;
for q = 1:numel(runs)
    parts = keep(runs{q});
    w1 = what(parts(1).k, :);
    w2 = what(parts(end).k, :);
    name = sprintf('%s_inner%d', names, q);
    if w1(1) == 2
        % fan face: the crease point offset along the fan of normals, an exact circular arc
        X = segs(w1(2)).ctrl(end, :);
        c2 = arc2(X, parts(1).X0, parts(1).X1, C.d);
        nominal = [0 1];
        frq = struct('name', name, 'n_nodes', 0, 'n_knots', numel(c2.knots), 'n_passes', 0, 'n_removed', 0, 'cap_reached', false);
    elseif all(arrayfun(@(p) segs(what(p.k, 2)).linear, parts))
        % straight segments: translated exactly, ends at the trim points
        pts = parts(1).X0;
        br = 0;
        L = [];
        for p = parts
            pts = [pts; p.X1]; %#ok<AGROW>
            L(end + 1) = segs(what(p.k, 2)).ub - segs(what(p.k, 2)).ua; %#ok<AGROW>
        end
        sig = [0 cumsum(L) / sum(L)];
        sig(end) = 1;
        pts = snap_axis(pts, parts, segs, what, open0, open1, nseg);
        c2 = struct('degree', 1, 'ctrl', pts, 'knots', [0 sig 1], 'weights', []);
        nominal = [segs(w1(2)).ua segs(w2(2)).ub];
        frq = struct('name', name, 'n_nodes', 0, 'n_knots', numel(c2.knots), 'n_passes', 0, 'n_removed', 0, 'cap_reached', false);
    else
        [c2, frq] = fit_run(parts, segs, what, C, sn, open0, open1, nseg, name, kf);
        nominal = [segs(w1(2)).ua segs(w2(2)).ub];
    end
    surf = revolve(c2, nominal, ax, U, b, G(1).surf);
    ent = new_entry(G(1), name, surf, ~G(1).outward, 'rev_z');
    E = [E, ent]; %#ok<AGROW>
    fr(end + 1) = frq; %#ok<AGROW>
end
ok = true;
end

function keep = clip_axis(keep, F)
% where the offset of a chain end crosses the axis (the hoop direction folds where the hoop radius
% is below d, e.g. at a pole whose curvature radius is below d), the part beyond the axis meets
% its rotated copies: the chain then starts (ends) where it comes back to the axis, a cone tip
for pass = 1:2
    if isempty(keep)
        return
    end
    q = 1;
    if pass == 2
        q = numel(keep);
    end
    p = keep(q);
    s = linspace(p.s0, p.s1, 401)';
    Q = F{p.k}(s);
    % the end samples are the chain ends themselves (on the axis at a pole, up to rounding)
    neg = [false; Q(2:end - 1, 1) < 0; false];
    if ~any(neg)
        continue
    end
    if pass == 1
        i = find(neg, 1, 'last');
        if i == numel(s)
            keep(q) = [];
            continue
        end
        sc = root_r(F{p.k}, s(i), s(i + 1));
        keep(q).s0 = sc;
        X = F{p.k}(sc);
        keep(q).X0 = [0 X(2)];
        keep(q).trim0 = true;
    else
        i = find(neg, 1);
        if i == 1
            keep(q) = [];
            continue
        end
        sc = root_r(F{p.k}, s(i - 1), s(i));
        keep(q).s1 = sc;
        X = F{p.k}(sc);
        keep(q).X1 = [0 X(2)];
        keep(q).trim1 = true;
    end
end
end

function x = root_r(f, a, b)
% bisection of the radial coordinate of f to adjacent doubles
fa = f(a);
fa = fa(1);
while true
    x = a + (b - a) / 2;
    if x <= a || x >= b
        return
    end
    fx = f(x);
    if sign(fx(1)) == sign(fa)
        a = x;
        fa = fx(1);
    else
        b = x;
    end
end
end

function pts = snap_axis(pts, parts, segs, what, open0, open1, nseg)
% an end of the chain on the axis (a pole of the outer surface) stays on it
if ~open0 && what(parts(1).k, 1) == 1 && what(parts(1).k, 2) == 1 && segs(1).ctrl(1, 1) == 0 && ~parts(1).trim0
    pts(1, 1) = 0;
end
if ~open1 && what(parts(end).k, 1) == 1 && what(parts(end).k, 2) == nseg && segs(nseg).ctrl(end, 1) == 0 && ...
        ~parts(end).trim1
    pts(end, 1) = 0;
end
end

function [c2, frq] = fit_run(parts, segs, what, C, sn, open0, open1, nseg, name, kf)
% least-squares cubic through the exact offset points of a smooth run, judged in the meridian
% plane against the outer profile (t_local = distance to the outer segments)
L = arrayfun(@(p) p.s1 - p.s0, parts);
sig = [0 cumsum(L) / sum(L)];
sig(end) = 1;
prob = struct();
prob.degree = 3;
prob.breaks = sig(2:end - 1);
prob.sample = @(s) run_sample(s, parts, sig, segs, what, C.d, sn);
X0 = parts(1).X0;
X1 = parts(end).X1;
X0 = snap_axis([X0; X1], parts, segs, what, open0, open1, nseg);
X1 = X0(2, :);
X0 = X0(1, :);
prob.fix0 = X0;
prob.fix1 = X1;
s1 = segs(what(parts(1).k, 2));
s2 = segs(what(parts(end).k, 2));
prob.tie0 = [false, ~parts(1).trim0 && X0(1) == 0 && s1.ctrl(1, 2) == s1.ctrl(2, 2)];
prob.tie1 = [false, ~parts(end).trim1 && X1(1) == 0 && s2.ctrl(end, 2) == s2.ctrl(end - 1, 2)];
prob.monotone = 2;
% initial knots: the outer knots inside the run and enough to keep every span's turn of the
% outer tangent within pi/4 (dense where the outer surface curves sharply)
kn = [];
for j = 1:numel(parts)
    sg = segs(what(parts(j).k, 2));
    ku = unique(sg.knots);
    ku = ku(ku > parts(j).s0 & ku < parts(j).s1);
    grid = [parts(j).s0, ku, parts(j).s1];
    for g = 1:numel(grid) - 1
        [~, Ta] = seg_eval(sg, grid(g));
        [~, Tb] = seg_eval(sg, grid(g + 1));
        ang = acos(max(-1, min(1, (Ta * Tb') / (norm(Ta) * norm(Tb)))));
        m = max(1, ceil(ang / (pi / 4)));
        uu = grid(g) + (1:m - 1) / m * (grid(g + 1) - grid(g));
        kn = [kn, sig(j) + ([uu, grid(g + 1)] - parts(j).s0) / (parts(j).s1 - parts(j).s0) * (sig(j + 1) - sig(j))]; %#ok<AGROW>
    end
end
kn = unique(kn);
kn = kn(kn > 0 & kn < 1 & ~ismember(kn, prob.breaks));
prob.knots = kn;
prob.judge = @(c) judge_run(c, segs, C);
fopts = struct('max_passes', C.max_passes);
if ~isempty(kf)
    i = find(strcmp({kf.patches.name}, name), 1);
    if isempty(i)
        error('mwecmass:solid:KnotsFromMismatch', 'offset_surface: knots_from has no patch %s', name);
    end
    k = kf.patches(i).surf.knots{1};
    k = (k - k(1)) / (k(end) - k(1));
    ki = k(4:end - 4);
    [u, ~, j] = unique(ki);
    mult = accumarray(j(:), 1)';
    prob.breaks = u(mult >= 3);
    prob.knots = u(mult < 3);
    fopts.knots_fixed = true;
end
[c2, r] = mwecmass.solid.fit_bspline_surface(prob, fopts);
frq = struct('name', name, 'n_nodes', r.n_nodes, 'n_knots', r.n_knots, 'n_passes', r.n_passes, ...
    'n_removed', r.n_removed, 'cap_reached', r.cap_reached);
end

function Q = run_sample(s, parts, sig, segs, what, d, sn)
Q = zeros(numel(s), 2);
for j = 1:numel(parts)
    in = s >= sig(j) & s <= sig(j + 1);
    if ~any(in)
        continue
    end
    u = parts(j).s0 + (s(in) - sig(j)) / (sig(j + 1) - sig(j)) * (parts(j).s1 - parts(j).s0);
    Q(in, :) = offset_point(segs(what(parts(j).k, 2)), u, d, sn);
end
end

function [ok, info] = judge_run(c, segs, C)
% M1/M2 in the meridian plane: distance of check points of the fitted profile (Gauss points and
% the midpoint of every span) to the whole outer profile chain
ku = unique(c.knots);
[xg, ~] = gauss_legendre(C.n_gauss);
nsp = numel(ku) - 1;
s = zeros(0, 1);
sp = zeros(0, 1);
for j = 1:nsp
    sj = [ku(j) + (ku(j + 1) - ku(j)) * (xg + 1) / 2; (ku(j) + ku(j + 1)) / 2];
    s = [s; sj]; %#ok<AGROW>
    sp = [sp; repmat(j, numel(sj), 1)]; %#ok<AGROW>
end
X = mwecmass.solid.eval_bspline_curve(c, s);
dist = profile_distance(segs, X);
ok = true(1, nsp);
for j = 1:nsp
    dj = dist(sp == j);
    ok(j) = all(dj >= C.t & dj <= C.t + C.eps_fit & dj >= C.t_min);
end
info = struct('t_local_min', min(dist), 't_local_max', max(dist));
end

function dist = profile_distance(segs, X)
% distance from the points X to the outer profile segments: Newton with descent safeguard on each
% segment from the three nearest of 81 samples (several local minima occur behind a fold),
% skipping segments whose sample box is farther than the best sample distance so far
dist = Inf(size(X, 1), 1);
G = cell(1, numel(segs));
for k = 1:numel(segs)
    s = linspace(segs(k).ua, segs(k).ub, 81)';
    G{k} = {s, seg_eval(segs(k), s)};
    D = min((X(:, 1) - G{k}{2}(:, 1)').^2 + (X(:, 2) - G{k}{2}(:, 2)').^2, [], 2);
    dist = min(dist, sqrt(D));
end
for k = 1:numel(segs)
    sg = segs(k);
    s = G{k}{1};
    Q = G{k}{2};
    ext = max(max(Q, [], 1) - min(Q, [], 1));
    lb = sqrt(sum(max(0, max(min(Q, [], 1) - X, X - max(Q, [], 1))).^2, 2));
    idx = find(lb <= dist + ext);
    if isempty(idx)
        continue
    end
    Xi = X(idx, :);
    D = (Xi(:, 1) - Q(:, 1)').^2 + (Xi(:, 2) - Q(:, 2)').^2;
    [~, order] = sort(D, 2);
    for start = 1:3
    u = s(order(:, start));
    dcur = sqrt(sum((seg_eval(sg, u) - Xi).^2, 2));
    for it = 1:30
        [P0, P1, P2] = seg_eval(sg, u);
        r = P0 - Xi;
        g = sum(r .* P1, 2);
        H = sum(P1 .* P1, 2) + sum(r .* P2, 2);
        H(H <= 0) = sum(P1(H <= 0, :).^2, 2);
        du = -g ./ max(H, realmin);
        % descent: halve a step that does not decrease the distance
        lam = ones(size(u));
        for h = 1:12
            un = min(max(u + lam .* du, sg.ua), sg.ub);
            dn = sqrt(sum((seg_eval(sg, un) - Xi).^2, 2));
            worse = dn > dcur;
            if ~any(worse)
                break
            end
            lam(worse) = lam(worse) / 2;
        end
        un(worse) = u(worse);
        dn(worse) = dcur(worse);
        if isequal(un, u)
            break
        end
        u = un;
        dcur = dn;
    end
    dist(idx) = min(dist(idx), dcur);
    end
end
end

function [P, T, A] = seg_eval(sg, u)
c = struct('degree', sg.degree, 'ctrl', sg.ctrl, 'knots', sg.knots, 'weights', sg.w);
switch nargout
    case {0, 1}
        P = mwecmass.solid.eval_bspline_curve(c, u(:));
    case 2
        [P, T] = mwecmass.solid.eval_bspline_curve(c, u(:));
    otherwise
        [P, T, A] = mwecmass.solid.eval_bspline_curve(c, u(:));
end
end

function [Q, dQ] = offset_point(sg, u, d, sn)
% exact offset P + d n of a profile segment, n = sn * rot(T) / |T|, and its derivative
[P, T, A] = seg_eval(sg, u);
nt = sqrt(sum(T.^2, 2));
n = sn * [T(:, 2), -T(:, 1)] ./ nt;
Q = P + d * n;
if nargout > 1
    dn = sn * ([A(:, 2), -A(:, 1)] ./ nt - [T(:, 2), -T(:, 1)] .* (sum(T .* A, 2) ./ nt.^3));
    dQ = T + d * dn;
end
end

function [Q, dQ] = fan_point(X, na, nb, d, th)
a0 = atan2(na(2), na(1));
a1 = atan2(nb(2), nb(1));
da = mod(a1 - a0 + pi, 2 * pi) - pi;
ang = a0 + th(:) * da;
Q = X + d * [cos(ang), sin(ang)];
dQ = d * da * [-sin(ang), cos(ang)];
end

function c = arc2(X, A, B, d)
% circular arc of radius d about X from A to B (angle below pi): rational quadratic
va = A - X;
vb = B - X;
cth = (va * vb') / d^2;
M = X + (va + vb) / (1 + cth);
c = struct('degree', 2, 'ctrl', [A; M; B], 'knots', [0 0 0 1 1 1], 'weights', [1; sqrt((1 + cth) / 2); 1]);
end

function v = rot(T)
v = [T(:, 2), -T(:, 1)];
end

function u = seg_root(sg, z)
% parameter where the segment's height is z (monotone segment): bisection to adjacent doubles
a = sg.ua;
b = sg.ub;
fa = seg_eval(sg, a);
fa = fa(2) - z;
while true
    m = a + (b - a) / 2;
    if m <= a || m >= b
        break
    end
    fm = seg_eval(sg, m);
    fm = fm(2) - z;
    if sign(fm) == sign(fa)
        a = m;
        fa = fm;
    else
        b = m;
    end
end
u = a;
end

function sg = restrict(sg, a, b)
c = sub_curve(struct('degree', sg.degree, 'ctrl', sg.ctrl, 'knots', sg.knots, 'weights', sg.w), a, b);
sg.ctrl = c.ctrl;
sg.w = c.weights;
sg.knots = c.knots;
sg.ua = a;
sg.ub = b;
sg.flatz = all(c.ctrl(:, 2) == c.ctrl(1, 2));
end

function c = sub_curve(c, a, b)
% the part [a, b] of c by knot insertion (Piegl & Tiller A5.1), keeping its parameter
p = c.degree;
for x = [a b]
    if x > c.knots(1) && x < c.knots(end)
        r = p - sum(c.knots == x);
        if r > 0
            c = insert(c, x, r);
        end
    end
end
i0 = 1;
if a > c.knots(1)
    i0 = find(c.knots == a, 1) - 1;
end
i1 = size(c.ctrl, 1);
if b < c.knots(end)
    i1 = find(c.knots == b, 1) - 1;
end
k = c.knots(c.knots > a & c.knots < b);
c.ctrl = c.ctrl(i0:i1, :);
if ~isempty(c.weights)
    c.weights = c.weights(i0:i1);
end
c.knots = [repmat(a, 1, p + 1), k, repmat(b, 1, p + 1)];
end

function c = insert(c, u, r)
p = c.degree;
knots = c.knots;
n = size(c.ctrl, 1);
w = c.weights;
rational = ~isempty(w);
if ~rational
    w = ones(n, 1);
end
k = find(knots <= u, 1, 'last');
s = sum(knots == u);
Q = zeros(n + r, size(c.ctrl, 2));
QW = zeros(n + r, 1);
Q(1:k - p, :) = c.ctrl(1:k - p, :);
QW(1:k - p) = w(1:k - p);
Q(k - s + r:n + r, :) = c.ctrl(k - s:n, :);
QW(k - s + r:n + r) = w(k - s:n);
Rp = c.ctrl(k - p:k - s, :);
Rw = w(k - p:k - s);
L = k - p;
for j = 1:r
    L = k - p + j;
    for i = 0:p - j - s
        al = (u - knots(L + i)) / (knots(i + k + 1) - knots(L + i));
        wn = al * Rw(i + 2) + (1 - al) * Rw(i + 1);
        if Rw(i + 2) == Rw(i + 1)
            wn = Rw(i + 1);
        end
        Pn = (al * Rw(i + 2) * Rp(i + 2, :) + (1 - al) * Rw(i + 1) * Rp(i + 1, :)) / wn;
        same = Rp(i + 2, :) == Rp(i + 1, :);
        Pn(same) = Rp(i + 1, same);
        Rp(i + 1, :) = Pn;
        Rw(i + 1) = wn;
    end
    Q(L, :) = Rp(1, :);
    QW(L) = Rw(1);
    Q(k + r - j - s, :) = Rp(p - j - s + 1, :);
    QW(k + r - j - s) = Rw(p - j - s + 1);
end
for i = L + 1:k - s - 1
    Q(i, :) = Rp(i - L + 1, :);
    QW(i) = Rw(i - L + 1);
end
c.ctrl = Q;
c.knots = [knots(1:k), repmat(u, 1, r), knots(k + 1:end)];
if rational
    c.weights = QW;
end
end

function surf = revolve(c2, nominal, ax, U, b, outer)
% the inner profile (r, z) on the outer surface's arc pattern: row i = axis + r_i U_j at height z_i
n = size(c2.ctrl, 1);
m = size(U, 1);
ctrl = zeros(n, m, 3);
for i = 1:n
    r = c2.ctrl(i, 1);
    ctrl(i, :, 1) = ax(1) + r * U(:, 1)';
    ctrl(i, :, 2) = ax(2) + r * U(:, 2)';
    ctrl(i, :, 3) = c2.ctrl(i, 2);
end
a = c2.weights;
if isempty(a)
    a = ones(n, 1);
end
k = nominal(1) + (c2.knots - c2.knots(1)) / (c2.knots(end) - c2.knots(1)) * (nominal(2) - nominal(1));
k(1:c2.degree + 1) = nominal(1);
k(end - c2.degree:end) = nominal(2);
surf = struct('type', 'bspline', 'degree', [c2.degree outer.degree(2)], 'ctrl', ctrl, ...
    'knots', {{k, outer.knots{2}}}, 'weights', a(:) * b(:)');
end

% =============================================================== ruled_parallel

function E = planar_inner(G, P, C)
% a planar face whose four boundaries are convex creases with planar neighbours: translated by d
% and trimmed by the neighbours' offset planes
E = [];
if numel(G) ~= 1
    return
end
e = G(1);
s = e.surf;
if ~isequal(s.degree, [1 1]) || ~isempty(s.weights) && any(s.weights(:) ~= s.weights(1))
    return
end
[n, c0] = plane_of(e);
if isempty(n)
    return
end
fields = {'seam_v0', 'seam_u1', 'seam_v1', 'seam_u0'};
NB = zeros(4, 4);
for bnd = 1:4
    nb = e.(fields{bnd});
    if isempty(nb)
        return
    end
    q = P(nb(1));
    [n2, c2] = plane_of(q);
    if isempty(n2) || ~isequal(q.surf.degree, [1 1])
        return
    end
    % convex: the neighbour lies on the inner side of this face's plane
    ctr = mean(reshape(q.surf.ctrl, [], 3), 1);
    if n * ctr' - c0 >= 0
        return
    end
    NB(bnd, :) = [n2, c2];
end
% offset planes n.x = c - d (outward unit normals)
own = [n, c0 - C.d];
NB(:, 4) = NB(:, 4) - C.d;
corner_b = [1 4; 1 2; 3 4; 3 2];
ij = [1 1; 2 1; 1 2; 2 2];
ctrl = zeros(2, 2, 3);
for q = 1:4
    Pl = [own; NB(corner_b(q, 1), :); NB(corner_b(q, 2), :)];
    x = (Pl(:, 1:3) \ Pl(:, 4))';
    for ax3 = 1:3
        for r = 1:3
            nr = Pl(r, 1:3);
            if nr(ax3) ~= 0 && all(nr([1:ax3 - 1, ax3 + 1:3]) == 0)
                x(ax3) = Pl(r, 4) / nr(ax3);
            end
        end
    end
    ctrl(ij(q, 1), ij(q, 2), :) = reshape(x, 1, 1, 3);
end
ctrl(:, 2, 3) = ctrl(:, 1, 3);
surf = s;
surf.ctrl = ctrl;
E = new_entry(e, [e.name '_inner'], surf, ~e.outward, e.offset_kind);
E.u_range = e.u_range;
end

function [n, c] = plane_of(e)
% outward unit normal and offset of a planar bilinear face (its four corners coplanar, exactly)
n = [];
c = [];
if isempty(e.surf) || ~isequal(e.surf.degree, [1 1])
    return
end
Q = reshape(e.surf.ctrl, [], 3);
nn = cross(Q(2, :) - Q(1, :), Q(3, :) - Q(1, :));
if norm(nn) == 0
    return
end
nn = nn / norm(nn);
for k = 1:3
    if sum(nn ~= 0) == 1 && nn(k) ~= 0
        nn(k) = sign(nn(k));
    end
end
dev = abs((Q - Q(1, :)) * nn');
if any(dev > 16 * eps * max(1, max(abs(Q(:)))))
    return
end
[~, Su, Sv] = mwecmass.solid.eval_bspline_surface(e.surf, mean(e.surf.knots{1}([1 end])), mean(e.surf.knots{2}([1 end])));
if (cross(Su, Sv) * nn' > 0) ~= e.outward
    nn = -nn;
end
n = nn;
c = n * Q(1, :)';
end

function E = ruled_inner(G, P, built)
% ruled between the inner boundary of the rev_z neighbour across a tangent-continuous v-seam and
% its translate by the ruling (every control ruling of the patch equal, bitwise)
E = [];
R0 = [];
for k = 1:numel(G)
    s = G(k).surf;
    if size(s.ctrl, 2) ~= 2
        return
    end
    R = reshape(s.ctrl(:, 2, :) - s.ctrl(:, 1, :), [], 3);
    if any(any(R ~= R(1, :))) || (~isempty(R0) && ~isequal(R(1, :), R0))
        return
    end
    R0 = R(1, :);
end
e = G(1);
side = 0;
for bnd = [1 3]
    nb = e.(seam_name(bnd));
    if isempty(nb) || ~any(nb(2) == [1 3])
        continue
    end
    q = P(nb(1));
    if ~strcmp(q.offset_kind, 'rev_z') || isempty(built{q.visible}) || ~smooth_seam(e, bnd, q, nb(2))
        continue
    end
    if ~isequal(boundary_curve(e.surf, bnd).ctrl, boundary_curve(q.surf, nb(2)).ctrl)
        continue
    end
    side = bnd;
    break
end
if side == 0
    return
end
nbe = built{q.visible};
for j = 1:numel(nbe)
    f = nbe(j).surf;
    col = 1;
    if nb(2) == 3
        col = size(f.ctrl, 2);
    end
    c1 = f.ctrl(:, col, :);
    shift = reshape(R0, 1, 1, 3);
    if side == 1
        ctrl = cat(2, c1, c1 + shift);
    else
        ctrl = cat(2, c1 - shift, c1);
    end
    surf = struct('type', 'bspline', 'degree', [f.degree(1) 1], 'ctrl', ctrl, 'knots', {{f.knots{1}, [0 0 1 1]}}, ...
        'weights', []);
    if ~isempty(f.weights)
        surf.weights = [f.weights(:, col) f.weights(:, col)];
    end
    E = [E, new_entry(e, sprintf('%s_inner%d', e.name, j), surf, ~e.outward, e.offset_kind)]; %#ok<AGROW>
end
end

function f = seam_name(b)
fields = {'seam_v0', 'seam_u1', 'seam_v1', 'seam_u0'};
f = fields{b};
end

function tf = smooth_seam(e, be, q, bq)
% the cross-boundary control legs of the two patches are parallel and point the same way (rounding
% of the control points: 16 ulp of the leg lengths) at every row with nonzero legs
s = e.surf;
f = q.surf;
if be == 1
    la = s.ctrl(:, 1, :) - s.ctrl(:, 2, :);
else
    la = s.ctrl(:, end, :) - s.ctrl(:, end - 1, :);
end
if bq == 1
    lb = f.ctrl(:, 2, :) - f.ctrl(:, 1, :);
else
    lb = f.ctrl(:, end - 1, :) - f.ctrl(:, end, :);
end
la = reshape(la, [], 3);
lb = reshape(lb, [], 3);
if size(la, 1) ~= size(lb, 1)
    tf = false;
    return
end
tf = true;
for i = 1:size(la, 1)
    if all(la(i, :) == 0) || all(lb(i, :) == 0)
        continue
    end
    cr = norm(cross(la(i, :), lb(i, :)));
    if cr > 16 * eps * norm(la(i, :)) * norm(lb(i, :)) || la(i, :) * lb(i, :)' <= 0
        tf = false;
        return
    end
end
end

% =============================================================== seams

function P = find_seams(P)
fields = {'seam_v0', 'seam_u1', 'seam_v1', 'seam_u0'};
n = numel(P);
B = cell(n, 4);
for k = 1:n
    for b = 1:4
        B{k, b} = boundary_curve(P(k).surf, b);
        P(k).(fields{b}) = [];
    end
end
for k = 1:n
    for b = 1:4
        ck = B{k, b};
        if ~isempty(P(k).(fields{b})) || all(all(ck.ctrl == ck.ctrl(1, :)))
            continue
        end
        for j = 1:n
            hit = false;
            for c = 1:4
                if (j == k && c == b) || ~isempty(P(j).(fields{c}))
                    continue
                end
                if same_curve(ck, B{j, c})
                    P(k).(fields{b}) = [j c];
                    P(j).(fields{c}) = [k b];
                    hit = true;
                    break
                end
            end
            if hit
                break
            end
        end
    end
end
end

function c = boundary_curve(s, b)
W = s.weights;
switch b
    case 1
        c = struct('degree', s.degree(1), 'ctrl', reshape(s.ctrl(:, 1, :), [], 3), 'knots', s.knots{1}, 'weights', []);
        if ~isempty(W), c.weights = W(:, 1); end
    case 2
        c = struct('degree', s.degree(2), 'ctrl', reshape(s.ctrl(end, :, :), [], 3), 'knots', s.knots{2}, 'weights', []);
        if ~isempty(W), c.weights = W(end, :)'; end
    case 3
        c = struct('degree', s.degree(1), 'ctrl', reshape(s.ctrl(:, end, :), [], 3), 'knots', s.knots{1}, 'weights', []);
        if ~isempty(W), c.weights = W(:, end); end
    case 4
        c = struct('degree', s.degree(2), 'ctrl', reshape(s.ctrl(1, :, :), [], 3), 'knots', s.knots{2}, 'weights', []);
        if ~isempty(W), c.weights = W(1, :)'; end
end
end

function tf = same_curve(a, b)
% contract S1: same degree and number of knots, bitwise control points and weights (empty =
% ones), knots equal after the affine map to [0, 1] (4 * 2^-52); either direction
tf = false;
if a.degree ~= b.degree || numel(a.knots) ~= numel(b.knots) || size(a.ctrl, 1) ~= size(b.ctrl, 1)
    return
end
wa = a.weights;
wb = b.weights;
if isempty(wa), wa = ones(size(a.ctrl, 1), 1); end
if isempty(wb), wb = ones(size(b.ctrl, 1), 1); end
ka = (a.knots - a.knots(1)) / (a.knots(end) - a.knots(1));
kb = (b.knots - b.knots(1)) / (b.knots(end) - b.knots(1));
if isequal(a.ctrl, b.ctrl) && isequal(wa(:), wb(:)) && all(abs(ka - kb) <= 4 * 2^-52)
    tf = true;
    return
end
kr = fliplr(1 - kb);
tf = isequal(a.ctrl, flipud(b.ctrl)) && isequal(wa(:), flipud(wb(:))) && all(abs(ka - kr) <= 4 * 2^-52);
end

% =============================================================== metrics

function rep = judge_set(E, C, fitrep, prim)
% S2r: M1, M2 at the check points of every patch built on its own (mirrors are exact copies and
% take the numbers of their primary), M3 per patch and on horizontal sections of the whole set
fields = {'seam_v0', 'seam_u1', 'seam_v1', 'seam_u0'};
pr = struct('n_nodes', {}, 'n_knots', {}, 'n_passes', {}, 'n_removed', {}, 'n_check', {}, 't_local_min', {}, ...
    't_local_max', {}, 'M1', {}, 'M2', {}, 'M3', {}, 'M3_reason', {});
[xg, ~] = gauss_legendre(C.n_gauss);
cap = false;
zcheck = zeros(0, 1);
for k = 1:numel(E)
    e = E(k);
    q = struct('n_nodes', 0, 'n_knots', [numel(e.surf.knots{1}), numel(e.surf.knots{2})], 'n_passes', 0, ...
        'n_removed', 0, 'n_check', 0, 't_local_min', NaN, 't_local_max', NaN, 'M1', true, 'M2', true, ...
        'M3', true, 'M3_reason', '');
    v = e.visible;
    pv = prim(v);
    fr = fitrep{pv};
    base = regexprep(e.name, '^.*_inner', '_inner');
    if ~isempty(fr)
        i = find(cellfun(@(x) ~isempty(regexp(x, [base '$'], 'once')), {fr.name}), 1);
        if ~isempty(i)
            q.n_nodes = fr(i).n_nodes;
            q.n_passes = fr(i).n_passes;
            q.n_removed = fr(i).n_removed;
            cap = cap || fr(i).cap_reached;
        end
    end
    j = find([E.visible] == pv & cellfun(@(x) ~isempty(regexp(x, [base '$'], 'once')), {E.name}), 1);
    if pv ~= v && ~isempty(j) && j < k
        same = pr(j);
        same.n_knots = q.n_knots;
        pr(k) = same;
        continue
    end
    us = check_params(e.surf.knots{1}, xg);
    vs = check_params(e.surf.knots{2}, xg);
    [UU, VV] = ndgrid(us, vs);
    X = mwecmass.solid.eval_bspline_surface(e.surf, UU(:), VV(:));
    tl = surface_distance(C.outer, X);
    q.n_check = numel(tl);
    q.t_local_min = min(tl);
    q.t_local_max = max(tl);
    q.M1 = all(tl >= C.t_min);
    q.M2 = all(tl >= C.t & tl <= C.t + C.eps_fit);
    Z = e.surf.ctrl(:, :, 3);
    dz = diff(Z(:, 1));
    reason = {};
    if any(any(Z ~= Z(:, 1)))
        reason{end + 1} = 'rows of more than one z'; %#ok<AGROW>
    end
    if any(dz > 0) && any(dz < 0)
        reason{end + 1} = 'row heights not monotone'; %#ok<AGROW>
    end
    for b = 1:4
        c = boundary_curve(e.surf, b);
        if isempty(e.(fields{b})) && ~all(all(c.ctrl == c.ctrl(1, :))) && ~at_clip(c, C)
            reason{end + 1} = sprintf('boundary %d has no neighbour', b); %#ok<AGROW>
        end
    end
    if ~isempty(reason)
        q.M3 = false;
        q.M3_reason = strjoin(reason, '; ');
    end
    pr(k) = q;
    if e.z_range(1) ~= e.z_range(2)
        zz = mwecmass.solid.eval_bspline_surface(e.surf, us, repmat(vs(1), numel(us), 1));
        zcheck = [zcheck; zz(:, 3)]; %#ok<AGROW>
    end
end
% horizontal sections of the inner set: one simple closed loop inside the outer section
zz = reshape([E.z_range], 2, [])';
zlo = min(zz(:));
zhi = max(zz(:));
flatz = zz(zz(:, 1) == zz(:, 2), 1);
zcheck = unique(zcheck(zcheck > zlo & zcheck < zhi & ~ismember(zcheck, flatz)));
if numel(zcheck) > 40
    zcheck = zcheck(round(linspace(1, numel(zcheck), 40)));
end
sec_reason = '';
for h = zcheck'
    try
        Li = mwecmass.solid.slice_bspline_surface(E, h);
        Lo = mwecmass.solid.slice_bspline_surface(C.outer, h);
        if ~Li.simple
            sec_reason = sprintf('section at z = %.6g is not simple', h);
        elseif ~all(inpolygon(Li.pts(:, 1), Li.pts(:, 2), Lo.pts(:, 1), Lo.pts(:, 2)))
            sec_reason = sprintf('section at z = %.6g leaves the outer section', h);
        end
    catch err
        sec_reason = sprintf('section at z = %.6g: %s', h, err.message);
    end
    if ~isempty(sec_reason)
        break
    end
end
if ~isempty(sec_reason)
    for k = 1:numel(pr)
        pr(k).M3 = false;
        pr(k).M3_reason = strjoin([{pr(k).M3_reason}, {sec_reason}], '; ');
    end
end
rep = struct('patches', pr, 'ok', all([pr.M1]) && all([pr.M2]) && all([pr.M3]), 'cap_reached', cap);
end

function tf = at_clip(c, C)
% a boundary row at the clipped end of a set built over part of the hull
tf = all(c.ctrl(:, 3) == c.ctrl(1, 3)) && (c.ctrl(1, 3) >= C.zr(2) || c.ctrl(1, 3) <= C.zr(1));
end

function s = check_params(k, xg)
ku = unique(k);
s = zeros(0, 1);
for j = 1:numel(ku) - 1
    s = [s; ku(j) + (ku(j + 1) - ku(j)) * (xg + 1) / 2; (ku(j) + ku(j + 1)) / 2]; %#ok<AGROW>
end
end

function dmin = surface_distance(E, X)
% distance from every point of X to the faces E: on every knot-span cell of every face, Newton on
% (u, v) clamped to the cell, with descent safeguard and an active bound, from the three nearest of
% its 9 x 9 grid points (the distance may have several local minima in a cell, e.g. behind a fold;
% a cell is polynomial or rational without interior knots, so the iteration does not cross a
% crease). The grid points of all cells give an upper bound first; a cell is then searched only
% for the points whose distance to its grid's bounding box, grown by the box's size, is below it.
m = size(X, 1);
cells = struct('k', {}, 'lo', {}, 'hi', {}, 'U', {}, 'V', {}, 'G', {});
dmin = Inf(m, 1);
for k = 1:numel(E)
    s = E(k).surf;
    if isempty(s)
        continue
    end
    ku = unique(s.knots{1});
    kv = unique(s.knots{2});
    for a = 1:numel(ku) - 1
        for b = 1:numel(kv) - 1
            % at an interior knot the evaluator takes the next span: stay just below it
            hi = [ku(a + 1) - (a + 1 < numel(ku)) * eps(ku(a + 1)), kv(b + 1) - (b + 1 < numel(kv)) * eps(kv(b + 1))];
            us = ku(a) + (0:8)' / 8 * (hi(1) - ku(a));
            vs = kv(b) + (0:8)' / 8 * (hi(2) - kv(b));
            [UU, VV] = ndgrid(us, vs);
            G = mwecmass.solid.eval_bspline_surface(s, UU(:), VV(:));
            cells(end + 1) = struct('k', k, 'lo', [ku(a) kv(b)], 'hi', hi, 'U', UU(:), 'V', VV(:), 'G', G); %#ok<AGROW>
            D = (X(:, 1) - G(:, 1)').^2 + (X(:, 2) - G(:, 2)').^2 + (X(:, 3) - G(:, 3)').^2;
            dmin = min(dmin, sqrt(min(D, [], 2)));
        end
    end
end
for c = cells
    s = E(c.k).surf;
    G = c.G;
    ext = max(max(G, [], 1) - min(G, [], 1));
    lb = sqrt(sum(max(0, max(min(G, [], 1) - X, X - max(G, [], 1))).^2, 2));
    idx = find(lb <= dmin + ext);
    if isempty(idx)
        continue
    end
    Xi = X(idx, :);
    D = (Xi(:, 1) - G(:, 1)').^2 + (Xi(:, 2) - G(:, 2)').^2 + (Xi(:, 3) - G(:, 3)').^2;
    [~, order] = sort(D, 2);
    for start = 1:3
        j = order(:, start);
        u = c.U(j);
        v = c.V(j);
        dcur = sqrt(sum((mwecmass.solid.eval_bspline_surface(s, u, v) - Xi).^2, 2));
        for it = 1:25
            [S, Su, Sv, Suu, Suv, Svv] = mwecmass.solid.eval_bspline_surface(s, u, v);
            r = S - Xi;
            g1 = sum(r .* Su, 2);
            g2 = sum(r .* Sv, 2);
            h11 = sum(Su .* Su, 2) + sum(r .* Suu, 2);
            h12 = sum(Su .* Sv, 2) + sum(r .* Suv, 2);
            h22 = sum(Sv .* Sv, 2) + sum(r .* Svv, 2);
            bad = ~(h11 .* h22 - h12.^2 > 0 & h11 > 0);
            h11(bad) = sum(Su(bad, :).^2, 2) + realmin;
            h22(bad) = sum(Sv(bad, :).^2, 2) + realmin;
            h12(bad) = 0;
            dt = h11 .* h22 - h12.^2;
            du = -(h22 .* g1 - h12 .* g2) ./ dt;
            dv = -(h11 .* g2 - h12 .* g1) ./ dt;
            % a parameter held at a bound of the cell: Newton in the other one alone
            hu = (u <= c.lo(1) & du < 0) | (u >= c.hi(1) & du > 0);
            hv = (v <= c.lo(2) & dv < 0) | (v >= c.hi(2) & dv > 0);
            du(hv & ~hu) = -g1(hv & ~hu) ./ h11(hv & ~hu);
            dv(hv) = 0;
            dv(hu & ~hv) = -g2(hu & ~hv) ./ h22(hu & ~hv);
            du(hu) = 0;
            du(~isfinite(du)) = 0;
            dv(~isfinite(dv)) = 0;
            % descent: a step toward a stationary point that is not a minimum is halved until the
            % distance decreases, or dropped
            lam = ones(size(u));
            for h = 1:12
                un = min(max(u + lam .* du, c.lo(1)), c.hi(1));
                vn = min(max(v + lam .* dv, c.lo(2)), c.hi(2));
                dn = sqrt(sum((mwecmass.solid.eval_bspline_surface(s, un, vn) - Xi).^2, 2));
                worse = dn > dcur;
                if ~any(worse)
                    break
                end
                lam(worse) = lam(worse) / 2;
            end
            un(worse) = u(worse);
            vn(worse) = v(worse);
            dn(worse) = dcur(worse);
            if isequal(un, u) && isequal(vn, v)
                break
            end
            u = un;
            v = vn;
            dcur = dn;
        end
        dmin(idx) = min(dmin(idx), dcur);
    end
end
end

function [x, w] = gauss_legendre(n)
b = (1:n - 1) ./ sqrt(4 * (1:n - 1).^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D));
w = 2 * V(1, i)'.^2;
end
