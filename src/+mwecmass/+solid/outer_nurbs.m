function geo = outer_nurbs(model, cache, opts)
%OUTER_NURBS  Hull outer surface as NURBS patches from the .ms2 entity tree (contract F1, S1, S1b).
%
%   geo = mwecmass.solid.outer_nurbs(model)
%   geo = mwecmass.solid.outer_nurbs(model, cache, opts)
%
%   model: MS2Parser.parse output. cache: T1 boundary cache (precompute_boundary_cache) and
%   opts.t_min [m]: needed only when a patch takes the general path (error
%   mwecmass:solid:FitInputMissing otherwise). opts.max_passes: refinement cap of fitted faces;
%   opts.force_general (default false): every patch takes the general path.
%
%   Deck check first (the deck is re-read with the parser's line rules): an entity that the
%   visible surfaces depend on (named on the line of a visible surface, of a mirror's source, or
%   recursively of an entity already named) but that the parser did not read raises
%   mwecmass:solid:UnsupportedEntity; a Symmetry entry or MirrSurf plane other than x = 0 or y = 0
%   raises mwecmass:solid:UnsupportedMirror.
%
%   Exact path. Points: FramePoint, MirrPoint, AbsBead (on the converted curve at the parser's
%   parameter). Curves: BCurve as is (the parser's clamped uniform knots), Line degree 1, Arc as
%   a rational quadratic (weights 1, cos(theta/2), 1), BSubCurve by knot insertion at the beads,
%   PolyCurve2 joined with C0 knots at k/n after exact degree elevation to the highest degree,
%   ProjCurve by zeroing one coordinate of the control points, EdgeSnake as the boundary row or
%   column of the converted parent surface (edges 1 = v0, 2 = u1, 3 = v1, 4 = u0). Surfaces:
%   RevSurf about a vertical axis as profile x rational arc (segments of at most 90 degrees,
%   multiples of 90 degrees with exact 0 and +-1), RuledSurf degree 1 in v between its two curves
%   after one common degree and knot vector (both curves must share the parser's parameter map
%   and their weights), mirrors by negating one coordinate of the control points. The parameter
%   maps of the parser (uniform angle on arcs and revolutions, k/n pieces of a PolyCurve2) are
%   tracked so that beads, sub-curves and rulings land where the parser puts them.
%   A patch whose z depends on v alone is stored with u and v exchanged (swap_uv), so that z
%   depends on u; it takes the exact path when its row heights are monotone. Every exact patch is
%   split at the ends of constant-z intervals strictly inside its z range (existing knots, no root
%   search) and then, with split_bspline_surface, at every patch-corner height strictly inside its
%   z range. Then, repeated until nothing changes: a constant-z patch with a vertex (a corner of
%   another entry) inside one of its boundaries, or with a corner inside a row of another entry,
%   joins a flat region (general path) with every constant-z patch touching it at that height; a
%   patch with a vertex inside one of its rows, not counting the end points of seams between merged
%   patches, takes the general path; seams are boundaries that share both end points and are the
%   same curve (same degree, bitwise control points and weights, knots equal after the affine map
%   to [0, 1]), and a boundary shared end to end that is not the same curve sends the later patch
%   in visible order to the general path. A boundary with no neighbour that is not a pole raises
%   mwecmass:solid:HullNotClosed. Patches of the general path are handed to
%   mwecmass.solid.fit_z_faces. outward: S_u x S_v points out of the hull, decided on the sections
%   of the patches (F4).

if nargin < 2
    cache = [];
end
if nargin < 3 || isempty(opts)
    opts = struct();
end
deck = read_deck(model.filename);
check_deck(model, deck);
[~, stem] = fileparts(model.filename);

names = model.visible_surfs(:)';
topo = model.classify_visible_surfaces();
nvis = numel(names);
ctx = struct('model', model, 'memo', containers.Map());
force = isfield(opts, 'force_general') && ~isempty(opts.force_general) && opts.force_general;

P = repmat(empty_entry(), 1, nvis);
general = false(1, nvis);
for k = 1:nvis
    e = empty_entry();
    e.name = names{k};
    ent = model.entities(names{k});
    e.source = names{k};
    if strcmp(ent.type, 'MirrSurf')
        i = find(strcmp({topo.mirrors.name}, names{k}), 1);
        e.source = topo.mirrors(i).ultimate_source;
        e.flips = topo.mirrors(i).effective_flips;
    end
    e.type = model.entities(e.source).type;
    e.visible = k;
    c = convert_patch(ctx, e.source);
    if c.ok
        surf = c.surf;
        for f = 1:numel(e.flips)
            surf = flip_surface(surf, e.flips{f});
        end
        e.surf = surf;
        e.exact = true;
        e.z_of_u = true;
        e.swap_uv = c.swap;
        e.offset_kind = c.offset_kind;
        e = refresh(e);
    else
        general(k) = true;
    end
    P(k) = e;
end
if force
    general(:) = true;
end
for k = find(general)
    P(k).exact = false;
end

[P, general] = split_constant_z(P, general);
[P, general] = split_corner_heights(P, general);
[P, general] = find_seams(P, general, unparsed_lines(model, deck));

flat = struct('z', {}, 'normal_z', {}, 'visible', {});
if any(general)
    if isempty(cache) || ~isfield(opts, 't_min') || isempty(opts.t_min)
        error('mwecmass:solid:FitInputMissing', ['outer_nurbs: %s takes the general path, which needs ' ...
            'the boundary cache and opts.t_min'], strjoin(unique({P(general).name}), ', '));
    end
    fopts = struct('t_min', opts.t_min, 'kind', 'outer');
    if isfield(opts, 'max_passes')
        fopts.max_passes = opts.max_passes;
    end
    for k = find(general)
        P(k).exact = false;
    end
    [P, flat] = mwecmass.solid.fit_z_faces(model, cache, P, general, fopts);
end
P = set_outward(P);
zr = [P.z_range];
geo = struct('hull_name', stem, 'outer', P, 'z_range', [min(zr) max(zr)], 'analytic', [], 'flat', flat);
end

% =============================================================== entries

function e = empty_entry()
e = struct('name', '', 'source', '', 'type', '', 'flips', {{}}, 'surf', [], 'outward', false, ...
    'exact', false, 'z_of_u', false, 'u_range', [], 'z_range', [], 'offset_kind', '', ...
    'pole', [false false], 'c0_u', [], 'c0_v', [], 'seam_u0', [], 'seam_u1', [], ...
    'seam_v0', [], 'seam_v1', [], 'visible', [], 'fit', [], 'swap_uv', false);
end

function e = refresh(e)
s = e.surf;
e.u_range = [s.knots{1}(1) s.knots{1}(end)];
e.z_range = [s.ctrl(1, 1, 3) s.ctrl(end, 1, 3)];
e.pole = [collapsed(s.ctrl(1, :, :)) collapsed(s.ctrl(end, :, :))];
e.c0_u = full_knots(s.knots{1}, s.degree(1));
e.c0_v = full_knots(s.knots{2}, s.degree(2));
end

function c = full_knots(k, p)
% interior knot values of multiplicity >= p
u = unique(k);
u = u(u > k(1) & u < k(end));
c = u(arrayfun(@(x) sum(k == x), u) >= p);
end

function tf = collapsed(row)
row = reshape(row, [], 3);
tf = all(all(row == row(1, :)));
end

function surf = flip_surface(surf, plane)
switch plane
    case 'X'
        surf.ctrl(:, :, 1) = -surf.ctrl(:, :, 1);
    case 'Y'
        surf.ctrl(:, :, 2) = -surf.ctrl(:, :, 2);
    otherwise
        error('mwecmass:solid:UnsupportedMirror', 'outer_nurbs: mirror plane %s', plane);
end
end

% =============================================================== deck check

function deck = read_deck(filename)
% Entity lines of the deck, assembled with the line rules of MS2Parser.parse.
fid = fopen(filename, 'r');
raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
fclose(fid);
lines = raw{1};
deck = struct('symmetry', {{}}, 'type', {{}}, 'name', {{}}, 'tokens', {{}}, 'line', {{}});
for i = 1:numel(lines)
    L = strtrim(lines{i});
    if startsWith(L, 'Symmetry:')
        deck.symmetry = lower(strsplit(strtrim(L(10:end))));
    elseif strcmp(L, 'BeginModel;')
        break
    end
end
in_model = false;
buf = '';
for i = 1:numel(lines)
    L = strtrim(lines{i});
    if strcmp(L, 'BeginModel;')
        in_model = true;
        continue
    end
    if strcmp(L, 'EndModel;')
        break
    end
    if ~in_model || isempty(L) || startsWith(L, 'Attribute:')
        continue
    end
    if ~isempty(buf)
        if startsWith(lines{i}, ' ') || startsWith(lines{i}, char(9)) || startsWith(L, '{') || startsWith(L, 'A:')
            buf = [buf ' ' L]; %#ok<AGROW>
            if endsWith(L, ';')
                deck = add_line(deck, buf);
                buf = '';
            end
            continue
        end
        deck = add_line(deck, buf);
        buf = '';
    end
    if endsWith(L, ';')
        deck = add_line(deck, L);
    else
        buf = L;
    end
end
if ~isempty(buf)
    deck = add_line(deck, buf);
end
end

function deck = add_line(deck, line)
line = strtrim(regexprep(line, ';$', ''));
if isempty(line)
    return
end
tok = strsplit(line);
if numel(tok) < 2
    return
end
deck.type{end + 1} = tok{1};
deck.name{end + 1} = tok{2};
tok = strrep(strrep(tok(3:end), '{', ''), '}', '');
deck.tokens{end + 1} = tok(~cellfun(@isempty, tok));
deck.line{end + 1} = line;
end

function check_deck(model, deck)
bad = deck.symmetry(~ismember(deck.symmetry, {'x', 'y', ''}));
if ~isempty(bad)
    error('mwecmass:solid:UnsupportedMirror', 'outer_nurbs: Symmetry plane %s of %s (only x = 0 and y = 0)', ...
        strjoin(bad, ' '), model.filename);
end
todo = {};
for k = 1:numel(model.visible_surfs)
    nm = model.visible_surfs{k};
    if any(strcmp(deck.name, nm))
        todo{end + 1} = nm; %#ok<AGROW>
    else
        [src, ~] = model.resolve_mirror_chain(nm);
        todo{end + 1} = src; %#ok<AGROW>
    end
end
seen = containers.Map();
while ~isempty(todo)
    nm = todo{end};
    todo(end) = [];
    if seen.isKey(nm)
        continue
    end
    seen(nm) = true;
    i = find(strcmp(deck.name, nm), 1, 'last');
    if isempty(i)
        continue
    end
    if ~model.entities.isKey(nm)
        error('mwecmass:solid:UnsupportedEntity', ['outer_nurbs: %s (%s) is referenced by the hull but ' ...
            'MS2Parser does not read it: %s'], nm, deck.type{i}, deck.line{i});
    end
    if strcmp(deck.type{i}, 'MirrSurf')
        sl = find(strcmp(deck.tokens{i}, '/'), 1);
        plane = deck.tokens{i}{sl + 2};
        if isempty(strfind(plane, 'X=0')) && isempty(strfind(plane, 'Y=0')) %#ok<STREMP>
            error('mwecmass:solid:UnsupportedMirror', 'outer_nurbs: MirrSurf %s mirrors in %s (only x = 0 and y = 0)', nm, plane);
        end
    end
    tok = deck.tokens{i};
    hit = tok(ismember(tok, deck.name));
    todo = [todo, hit(~strcmp(hit, nm))]; %#ok<AGROW>
end
end

function lines = unparsed_lines(model, deck)
lines = {};
for i = 1:numel(deck.name)
    if ~model.entities.isKey(deck.name{i})
        lines{end + 1} = deck.line{i}; %#ok<AGROW>
    end
end
end

% =============================================================== conversion

function c = convert_patch(ctx, name)
c = struct('ok', false, 'reason', '', 'surf', [], 'swap', false, 'offset_kind', '');
sf = conv_surface(ctx, name);
if ~sf.ok
    c.reason = sf.reason;
    return
end
surf = sf.surf;
Z = surf.ctrl(:, :, 3);
if all(all(Z == Z(:, 1)))
    c.swap = false;
elseif all(all(Z == Z(1, :)))
    surf = swap_surface(surf);
    c.swap = true;
else
    c.reason = 'z depends on both parameters';
    return
end
Zr = surf.ctrl(:, 1, 3);
dz = diff(Zr);
if any(dz > 0) && any(dz < 0)
    c.reason = 'the row heights are not monotone';
    return
end
c.ok = true;
c.surf = surf;
c.offset_kind = sf.kind;
end

function surf = swap_surface(surf)
surf.ctrl = permute(surf.ctrl, [2 1 3]);
surf.degree = surf.degree([2 1]);
surf.knots = surf.knots([2 1]);
if ~isempty(surf.weights)
    surf.weights = surf.weights.';
end
end

function pt = conv_point(ctx, name)
% pt.ok, pt.p [1x3]
key = ['P:' name];
if ctx.memo.isKey(key)
    pt = ctx.memo(key);
    return
end
pt = struct('ok', false, 'reason', '', 'p', [NaN NaN NaN]);
m = ctx.model;
if ~m.entities.isKey(name)
    pt.reason = sprintf('point %s is not defined', name);
    ctx.memo(key) = pt;
    return
end
e = m.entities(name);
switch e.type
    case 'FramePoint'
        base = [0 0 0];
        par = '';
        if ~isempty(e.params.parent1)
            par = e.params.parent1;
        elseif ~isempty(e.params.parent2)
            par = e.params.parent2;
        end
        if ~isempty(par)
            b = conv_point(ctx, par);
            if ~b.ok
                pt.reason = b.reason;
                ctx.memo(key) = pt;
                return
            end
            base = b.p;
        end
        pt.p = base + e.params.offset;
        pt.ok = true;
    case 'MirrPoint'
        b = conv_point(ctx, e.params.source);
        if b.ok
            p = b.p;
            col = find(strcmp(e.params.plane, {'X', 'Y', 'Z'}), 1);
            p(col) = -p(col);
            pt.p = p;
            pt.ok = true;
        else
            pt.reason = b.reason;
        end
    case 'AbsBead'
        cv = conv_curve(ctx, e.params.parent_curve);
        if cv.ok
            s = pm_eval(cv.pmap, e.params.parameter);
            pt.p = mwecmass.solid.eval_bspline_curve(cv.curve, s);
            pt.ok = true;
        else
            pt.reason = cv.reason;
        end
    otherwise
        pt.reason = sprintf('%s is a %s, not converted exactly', name, e.type);
end
ctx.memo(key) = pt;
end

function cv = conv_curve(ctx, name)
% cv.ok, cv.curve (T9 curve struct, knots in [0, 1]), cv.pmap (parser parameter -> curve
% parameter), cv.ends {start entity, end entity} ('' when the end is no named point)
key = ['C:' name];
if ctx.memo.isKey(key)
    cv = ctx.memo(key);
    return
end
cv = struct('ok', false, 'reason', '', 'curve', [], 'pmap', [], 'ends', {{'', ''}});
m = ctx.model;
if ~m.entities.isKey(name)
    cv.reason = sprintf('curve %s is not defined', name);
    ctx.memo(key) = cv;
    return
end
e = m.entities(name);
p = e.params;
switch e.type
    case 'BCurve'
        n = numel(p.ctrl_pt_names);
        if n <= p.degree || p.degree < 1
            cv.reason = sprintf('BCurve %s has %d points for degree %d', name, n, p.degree);
        else
            [Q, ok, why] = points(ctx, p.ctrl_pt_names);
            if ok
                cv = ok_curve(cv, struct('degree', p.degree, 'ctrl', Q, ...
                    'knots', mwecmass.geometry.MS2Parser.make_clamped_knots(n, p.degree), 'weights', []), ...
                    pm_identity(), p.ctrl_pt_names([1 end]));
            else
                cv.reason = why;
            end
        end
    case 'Line'
        [Q, ok, why] = points(ctx, {p.pt_start, p.pt_end});
        if ok
            cv = ok_curve(cv, line_curve(Q(1, :), Q(2, :)), pm_identity(), {p.pt_start, p.pt_end});
        else
            cv.reason = why;
        end
    case 'Arc'
        [Q, ok, why] = points(ctx, {p.pt_start, p.pt_center, p.pt_end});
        if ok
            [crv, pm] = arc_curve(Q(1, :), Q(2, :), Q(3, :));
            if isempty(crv)
                cv.reason = sprintf('Arc %s: end point not at the start radius', name);
                ctx.memo(key) = cv;
                return
            end
            ends = {p.pt_start, p.pt_end};
            if crv.degree == 1 && isequal(crv.ctrl(1, :), crv.ctrl(2, :))
                ends = {'', ''};
            end
            cv = ok_curve(cv, crv, pm, ends);
        else
            cv.reason = why;
        end
    case 'PolyCurve2'
        cv = poly_curve(ctx, cv, p.curve_names);
    case 'BSubCurve'
        cv = sub_curve(ctx, cv, p.bead_names);
    case 'ProjCurve'
        sc = conv_curve(ctx, p.source);
        if sc.ok
            crv = sc.curve;
            col = find(strcmp(p.proj_plane, {'X', 'Y', 'Z'}), 1);
            crv.ctrl(:, col) = 0;
            cv = ok_curve(cv, crv, sc.pmap, {'', ''});
        else
            cv.reason = sc.reason;
        end
    case 'EdgeSnake'
        sf = conv_surface(ctx, p.surface_name);
        if sf.ok
            s = sf.surf;
            W = s.weights;
            switch p.edge_index
                case 1
                    crv = edge_curve(s, s.ctrl(:, 1, :), W, 1, 1);
                    pm = sf.umap;
                case 2
                    crv = edge_curve(s, s.ctrl(end, :, :), W, 2, size(s.ctrl, 1));
                    pm = sf.vmap;
                case 3
                    crv = edge_curve(s, s.ctrl(:, end, :), W, 1, size(s.ctrl, 2));
                    pm = sf.umap;
                case 4
                    crv = edge_curve(s, s.ctrl(1, :, :), W, 2, 1);
                    pm = sf.vmap;
                otherwise
                    crv = [];
                    pm = [];
            end
            if isempty(crv)
                cv.reason = sprintf('EdgeSnake %s: edge %g', name, p.edge_index);
            else
                cv = ok_curve(cv, crv, pm, {'', ''});
            end
        else
            cv.reason = sf.reason;
        end
    otherwise
        cv.reason = sprintf('%s is a %s, not converted exactly', name, e.type);
end
ctx.memo(key) = cv;
end

function cv = ok_curve(cv, crv, pm, ends)
cv.ok = true;
cv.curve = crv;
cv.pmap = pm;
cv.ends = ends;
end

function crv = edge_curve(s, rows, W, dir, idx)
crv = struct('degree', s.degree(dir), 'ctrl', reshape(rows, [], 3), 'knots', s.knots{dir}, 'weights', []);
if ~isempty(W)
    if dir == 1
        crv.weights = W(:, idx);
    else
        crv.weights = W(idx, :)';
    end
end
end

function [Q, ok, why] = points(ctx, list)
Q = zeros(numel(list), 3);
ok = true;
why = '';
for i = 1:numel(list)
    pt = conv_point(ctx, list{i});
    if ~pt.ok
        ok = false;
        why = pt.reason;
        return
    end
    Q(i, :) = pt.p;
end
end

function crv = line_curve(A, B)
crv = struct('degree', 1, 'ctrl', [A; B], 'knots', [0 0 1 1], 'weights', []);
end

function [crv, pm] = arc_curve(S, Cn, E)
% MS2Parser arc_evaluate: radius |S - Cn|, angle theta in (0, pi) from S toward E, uniform in t
vs = S - Cn;
ve = E - Cn;
r = norm(vs);
pm = pm_identity();
if r < 1e-14
    crv = line_curve(Cn, Cn);
    return
end
cp = cross(vs, ve);
if norm(cp) < 1e-14
    crv = line_curve(S, E);
    return
end
er = vs / r;
et = cross(cp / norm(cp), er);
th = atan2(dot(ve, et), dot(ve, er));
% The parser's arc ends at radius r; the named end point is taken (so that the arc joins the next
% curve by entity) when |E - Cn| equals r up to the rounding of the data. With m = max|S, Cn, E|,
% each coordinate of S - Cn and E - Cn is off by at most 2 eps(m) (both operands rounded to
% doubles, eps(m)/2 each, and the subtraction, eps(2m)/2), so each norm moves by at most
% 2 sqrt(3) eps(m), and each norm is computed to 2 eps of its value (squares, sums, square root).
m = max(abs([S Cn E]));
if abs(norm(ve) - r) > 4 * sqrt(3) * eps(m) + 4 * eps(max(r, norm(ve)))
    crv = [];
    return
end
cth = dot(vs, ve) / (r * norm(ve));
P1 = Cn + (vs + ve) / (1 + cth);
crv = struct('degree', 2, 'ctrl', [S; P1; E], 'knots', [0 0 0 1 1 1], 'weights', [1; sqrt((1 + cth) / 2); 1]);
pm = pm_arc(th / 2);
end

function cv = poly_curve(ctx, cv, list)
nc = numel(list);
parts = cell(1, nc);
for k = 1:nc
    parts{k} = conv_curve(ctx, list{k});
    if ~parts{k}.ok
        cv.reason = parts{k}.reason;
        return
    end
end
p = max(cellfun(@(x) x.curve.degree, parts));
rational = any(cellfun(@(x) ~isempty(x.curve.weights), parts));
ctrl = zeros(0, 3);
w = zeros(0, 1);
knots = zeros(1, p + 1);
pms = cell(1, nc);
for k = 1:nc
    c = crv_elevate(parts{k}.curve, p);
    wk = c.weights;
    if isempty(wk)
        wk = ones(size(c.ctrl, 1), 1);
    end
    if k > 1
        if ~isequal(c.ctrl(1, :), ctrl(end, :))
            cv.reason = sprintf('PolyCurve2 pieces %s and %s do not join', list{k - 1}, list{k});
            return
        end
        if wk(1) ~= w(end)
            wk = wk * (w(end) / wk(1));
            wk(1) = w(end);
        end
        ctrl = [ctrl; c.ctrl(2:end, :)]; %#ok<AGROW>
        w = [w; wk(2:end)]; %#ok<AGROW>
    else
        ctrl = c.ctrl;
        w = wk;
    end
    kk = c.knots(p + 2:end - p - 1);
    knots = [knots, (k - 1 + kk) / nc]; %#ok<AGROW>
    if k < nc
        knots = [knots, repmat(k / nc, 1, p)]; %#ok<AGROW>
    end
    pms{k} = parts{k}.pmap;
end
knots = [knots, ones(1, p + 1)];
crv = struct('degree', p, 'ctrl', ctrl, 'knots', knots, 'weights', []);
if rational
    crv.weights = w;
end
cv = ok_curve(cv, crv, pm_concat(pms), {parts{1}.ends{1}, parts{end}.ends{2}});
end

function cv = sub_curve(ctx, cv, beads)
m = ctx.model;
if isempty(beads) || ~all(cellfun(@(b) m.entities.isKey(b), beads))
    cv.reason = 'BSubCurve beads are not defined';
    return
end
b1 = m.entities(beads{1});
b2 = m.entities(beads{end});
if ~strcmp(b1.type, 'AbsBead') || ~strcmp(b2.type, 'AbsBead')
    cv.reason = 'BSubCurve on beads that are not AbsBead';
    return
end
par = conv_curve(ctx, b1.params.parent_curve);
if ~par.ok
    cv.reason = par.reason;
    return
end
t1 = b1.params.parameter;
t2 = b2.params.parameter;
s1 = pm_eval(par.pmap, t1);
s2 = pm_eval(par.pmap, t2);
if s1 == s2
    cv.reason = 'BSubCurve of zero length';
    return
end
crv = crv_extract(par.curve, min(s1, s2), max(s1, s2));
if s1 > s2
    crv = crv_reverse(crv);
end
% the end points are the beads themselves (the same curve point, evaluated once)
e1 = conv_point(ctx, beads{1});
e2 = conv_point(ctx, beads{end});
crv.ctrl(1, :) = e1.p;
crv.ctrl(end, :) = e2.p;
cv = ok_curve(cv, crv, pm_restrict(par.pmap, t1, t2, s1, s2), {beads{1}, beads{end}});
end

function sf = conv_surface(ctx, name)
% sf.ok, sf.surf (parser orientation), sf.umap, sf.vmap, sf.kind ('rev_z' | 'ruled_parallel' | '')
key = ['S:' name];
if ctx.memo.isKey(key)
    sf = ctx.memo(key);
    return
end
sf = struct('ok', false, 'reason', '', 'surf', [], 'umap', [], 'vmap', [], 'kind', '');
m = ctx.model;
if ~m.entities.isKey(name)
    sf.reason = sprintf('surface %s is not defined', name);
    ctx.memo(key) = sf;
    return
end
e = m.entities(name);
switch e.type
    case 'RevSurf'
        sf = rev_surface(ctx, sf, e);
    case 'RuledSurf'
        sf = ruled_surface(ctx, sf, e);
    case 'MirrSurf'
        src = conv_surface(ctx, e.params.source);
        sf = src;
        if src.ok
            sf.surf = flip_surface(src.surf, e.params.mirror_plane);
        end
    otherwise
        sf.reason = sprintf('%s is a %s, not converted exactly', name, e.type);
end
ctx.memo(key) = sf;
end

function sf = rev_surface(ctx, sf, e)
m = ctx.model;
pr = conv_curve(ctx, e.params.profile);
if ~pr.ok
    sf.reason = pr.reason;
    return
end
if ~m.entities.isKey(e.params.axis) || ~strcmp(m.entities(e.params.axis).type, 'Line')
    sf.reason = 'RevSurf axis is not a Line';
    return
end
ax_e = m.entities(e.params.axis);
[A, ok, why] = points(ctx, {ax_e.params.pt_start, ax_e.params.pt_end});
if ~ok
    sf.reason = why;
    return
end
B = A(2, :);
A = A(1, :);
% the 16-ulp bound of T1's root search on composed parser evaluations
vertical = abs(A(1) - B(1)) <= 16 * eps(max(1, abs(B(1)))) && abs(A(2) - B(2)) <= 16 * eps(max(1, abs(B(2))));
if ~vertical || A(3) == B(3)
    sf.reason = 'RevSurf about an axis that is not vertical';
    return
end
ax = B(1);
ay = B(2);
sgn = sign(B(3) - A(3));
Q = pr.curve.ctrl;
axis_ends = {ax_e.params.pt_start, ax_e.params.pt_end};
if any(strcmp(pr.ends{1}, axis_ends))
    Q(1, 1:2) = [ax ay];
end
if any(strcmp(pr.ends{2}, axis_ends))
    Q(end, 1:2) = [ax ay];
end
phi0 = e.params.angle_start;
Th = e.params.angle_end - phi0;
if Th == 0
    sf.reason = 'RevSurf of zero angle';
    return
end
nseg = ceil(abs(Th) / 90);
D = Th / nseg;
if abs(D) == 90
    wc = sqrt(2) / 2;
    den = 1;
else
    wc = cos(D * pi / 360);
    den = 1 + cos(D * pi / 180);
end
n = size(Q, 1);
nv = 2 * nseg + 1;
ctrl = zeros(n, nv, 3);
cs = zeros(nseg + 1, 2);
for j = 0:nseg
    [cs(j + 1, 1), cs(j + 1, 2)] = cos_sin_deg(sgn * (phi0 + j * D));
end
for i = 1:n
    rx = Q(i, 1) - ax;
    ry = Q(i, 2) - ay;
    q = [rx * cs(:, 1) - ry * cs(:, 2), rx * cs(:, 2) + ry * cs(:, 1)];
    for j = 0:nseg
        ctrl(i, 2 * j + 1, :) = [ax + q(j + 1, 1), ay + q(j + 1, 2), Q(i, 3)];
        if j < nseg
            mid = (q(j + 1, :) + q(j + 2, :)) / den;
            ctrl(i, 2 * j + 2, :) = [ax + mid(1), ay + mid(2), Q(i, 3)];
        end
    end
end
wv = ones(1, nv);
wv(2:2:end) = wc;
a = pr.curve.weights;
if isempty(a)
    a = ones(n, 1);
end
kv = [0 0 0, reshape(repmat((1:nseg - 1) / nseg, 2, 1), 1, []), 1 1 1];
sf.surf = struct('type', 'bspline', 'degree', [pr.curve.degree 2], 'ctrl', ctrl, ...
    'knots', {{pr.curve.knots, kv}}, 'weights', a * wv);
sf.umap = pr.pmap;
pms = cell(1, nseg);
for j = 1:nseg
    pms{j} = pm_arc(abs(D) * pi / 360);
end
sf.vmap = pm_concat(pms);
sf.kind = 'rev_z';
sf.ok = true;
end

function [c, s] = cos_sin_deg(phi)
k = phi / 90;
if k == round(k)
    i = mod(round(k), 4) + 1;
    cc = [1 0 -1 0];
    ss = [0 1 0 -1];
    c = cc(i);
    s = ss(i);
else
    c = cos(phi * pi / 180);
    s = sin(phi * pi / 180);
end
end

function sf = ruled_surface(ctx, sf, e)
c1 = conv_curve(ctx, e.params.curve1);
c2 = conv_curve(ctx, e.params.curve2);
if ~c1.ok || ~c2.ok
    sf.reason = [c1.reason c2.reason];
    return
end
if ~pm_equal(c1.pmap, c2.pmap)
    sf.reason = 'RuledSurf curves with different parameter maps';
    return
end
[a, b] = compatible(c1.curve, c2.curve);
if isempty(a.weights) ~= isempty(b.weights) || ~isequal(a.weights, b.weights)
    wa = a.weights;
    wb = b.weights;
    if isempty(wa), wa = ones(size(a.ctrl, 1), 1); end
    if isempty(wb), wb = ones(size(b.ctrl, 1), 1); end
    if ~isequal(wa, wb)
        sf.reason = 'RuledSurf curves with different weights';
        return
    end
    a.weights = wa;
end
n = size(a.ctrl, 1);
ctrl = zeros(n, 2, 3);
ctrl(:, 1, :) = reshape(a.ctrl, n, 1, 3);
ctrl(:, 2, :) = reshape(b.ctrl, n, 1, 3);
W = [];
if ~isempty(a.weights)
    W = [a.weights a.weights];
end
sf.surf = struct('type', 'bspline', 'degree', [a.degree 1], 'ctrl', ctrl, 'knots', {{a.knots, [0 0 1 1]}}, ...
    'weights', W);
sf.umap = c1.pmap;
sf.vmap = pm_identity();
R = b.ctrl - a.ctrl;
nz = find(any(R ~= 0, 2), 1);
parallel = ~isempty(nz) && all(all(cross(R, repmat(R(nz, :), n, 1), 2) == 0));
sf.kind = '';
if parallel
    sf.kind = 'ruled_parallel';
end
sf.ok = true;
end

function [a, b] = compatible(a, b)
% one degree and one knot vector for two curves (exact degree elevation and knot insertion)
p = max(a.degree, b.degree);
a = crv_elevate(a, p);
b = crv_elevate(b, p);
u = unique([a.knots b.knots]);
for x = u
    ma = sum(a.knots == x);
    mb = sum(b.knots == x);
    if ma < mb
        a = crv_insert(a, x, mb - ma);
    elseif mb < ma
        b = crv_insert(b, x, ma - mb);
    end
end
end

% =============================================================== curve operations

function c = crv_insert(c, x, r)
n = size(c.ctrl, 1);
W = c.weights;
[Q, QW, k] = insert_rows(reshape(c.ctrl, n, 1, 3), W, c.knots, c.degree, x, r);
c.ctrl = reshape(Q, [], 3);
c.weights = QW;
c.knots = k;
end

function c = crv_extract(c, a, b)
% the piece of c between parameters a < b, knots renormalised to [0, 1]
p = c.degree;
for x = [a b]
    if x > c.knots(1) && x < c.knots(end)
        c = crv_insert(c, x, p - sum(c.knots == x));
    end
end
i0 = find(c.knots == a, 1) - 1;
if a == c.knots(1)
    i0 = 1;
end
i1 = find(c.knots == b, 1) - 1;
if b == c.knots(end)
    i1 = size(c.ctrl, 1);
end
k = c.knots(c.knots > a & c.knots < b);
c.ctrl = c.ctrl(i0:i1, :);
if ~isempty(c.weights)
    c.weights = c.weights(i0:i1);
end
c.knots = [zeros(1, p + 1), (k - a) / (b - a), ones(1, p + 1)];
end

function c = crv_reverse(c)
c.ctrl = flipud(c.ctrl);
if ~isempty(c.weights)
    c.weights = flipud(c.weights);
end
c.knots = 1 - fliplr(c.knots);
c.knots(1:c.degree + 1) = 0;
c.knots(end - c.degree:end) = 1;
end

function c = crv_elevate(c, p)
% exact degree elevation: the homogeneous curve is a spline of degree p on the knots with every
% distinct value repeated (p - degree) more times; its control points solve the interpolation at
% the Greville abscissae (Schoenberg-Whitney, so the system is regular). A coordinate or weight
% that is equal for every input control point is copied, and the end points are kept.
q = c.degree;
if p == q
    return
end
u = unique(c.knots);
k = [];
for x = u
    k = [k, repmat(x, 1, sum(c.knots == x) + p - q)]; %#ok<AGROW>
end
n = numel(k) - p - 1;
g = zeros(n, 1);
for i = 1:n
    g(i) = mean(k(i + 1:i + p));
end
w = c.weights;
rational = ~isempty(w);
if ~rational
    w = ones(size(c.ctrl, 1), 1);
end
Bo = mwecmass.solid.eval_bspline_curve(struct('degree', q, 'ctrl', eye(size(c.ctrl, 1)), 'knots', c.knots, ...
    'weights', []), g);
Bn = mwecmass.solid.eval_bspline_curve(struct('degree', p, 'ctrl', eye(n), 'knots', k, 'weights', []), g);
H = Bn \ (Bo * [c.ctrl .* w, w]);
wn = H(:, 4);
Q = H(:, 1:3) ./ wn;
for col = 1:3
    if all(c.ctrl(:, col) == c.ctrl(1, col))
        Q(:, col) = c.ctrl(1, col);
    end
end
Q(1, :) = c.ctrl(1, :);
Q(end, :) = c.ctrl(end, :);
wn(1) = w(1);
wn(end) = w(end);
if all(w == w(1))
    wn(:) = w(1);
end
c.degree = p;
c.ctrl = Q;
c.knots = k;
if rational
    c.weights = wn;
else
    c.weights = [];
end
end

function [ctrl, W, knots] = insert_rows(ctrl, W, knots, p, u, r)
% Insert u r times into the first-parameter knot vector (Piegl & Tiller, The NURBS Book, 2nd ed.,
% A5.1, on homogeneous rows). A coordinate or weight equal in the two rows being combined is
% copied; each new row takes the z of its first column, so rows keep one z.
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
w = alpha * w1 + (1 - alpha) * w0;
same_w = w1 == w0;
w(same_w) = w1(same_w);
P = (alpha * P1 .* w1 + (1 - alpha) * P0 .* w0) ./ w;
same = P1 == P0;
P(same) = P1(same);
P(1, :, 3) = P(1, 1, 3);
end

% =============================================================== parameter maps
% A map is a struct array of pieces, increasing in t: parser interval [t0, t1] -> curve interval
% [s0, s1]; kind 'lin' (affine) or 'arc': the parser's uniform angle on a rational quadratic arc of
% half-angle ha, restricted to the part [ta, tb] of the full arc's parser interval [0, 1].

function pm = pm_identity()
pm = struct('t0', 0, 't1', 1, 's0', 0, 's1', 1, 'kind', 'lin', 'ha', 0, 'ta', 0, 'tb', 1);
end

function pm = pm_arc(ha)
pm = struct('t0', 0, 't1', 1, 's0', 0, 's1', 1, 'kind', 'arc', 'ha', ha, 'ta', 0, 'tb', 1);
end

function s = arc_sigma(tau, ha)
% rational quadratic arc of half-angle ha, weights (1, cos ha, 1): the point at angle phi from the
% bisector has curve parameter sigma with tan(phi/2) = tan(ha/2) (2 sigma - 1); the parser's angle
% is (2 tau - 1) ha from the bisector
s = 0.5 * (1 + tan((2 * tau - 1) * ha / 2) / tan(ha / 2));
s(tau == 0) = 0;
s(tau == 1) = 1;
end

function s = pm_eval(pm, t)
i = find([pm.t0] <= t, 1, 'last');
if isempty(i)
    i = 1;
end
q = pm(i);
if t == q.t0
    s = q.s0;
    return
end
if t == q.t1
    s = q.s1;
    return
end
if q.t0 == q.s0 && q.t1 == q.s1 && strcmp(q.kind, 'lin')
    s = t;
    return
end
tau = (t - q.t0) / (q.t1 - q.t0);
if strcmp(q.kind, 'arc')
    sa = arc_sigma(q.ta, q.ha);
    sb = arc_sigma(q.tb, q.ha);
    tau = (arc_sigma(q.ta + tau * (q.tb - q.ta), q.ha) - sa) / (sb - sa);
end
s = q.s0 + tau * (q.s1 - q.s0);
end

function out = pm_concat(pms)
nc = numel(pms);
out = pm_identity();
out(:) = [];
for k = 1:nc
    for q = pms{k}
        q.t0 = (k - 1 + q.t0) / nc;
        q.t1 = (k - 1 + q.t1) / nc;
        q.s0 = (k - 1 + q.s0) / nc;
        q.s1 = (k - 1 + q.s1) / nc;
        out(end + 1) = q; %#ok<AGROW>
    end
end
out = pm_merge(out);
end

function out = pm_restrict(pm, t1, t2, s1, s2)
% map of the sub-curve whose parser parameter is t1 + t (t2 - t1) and whose curve parameter is
% (s - s1) / (s2 - s1)
lo = min(t1, t2);
hi = max(t1, t2);
out = pm_identity();
out(:) = [];
for q = pm
    a = max(q.t0, lo);
    b = min(q.t1, hi);
    if b <= a
        continue
    end
    r = q;
    if strcmp(q.kind, 'arc')
        r.ta = q.ta + (a - q.t0) / (q.t1 - q.t0) * (q.tb - q.ta);
        r.tb = q.ta + (b - q.t0) / (q.t1 - q.t0) * (q.tb - q.ta);
    end
    sa = (pm_eval(pm, a) - s1) / (s2 - s1);
    sb = (pm_eval(pm, b) - s1) / (s2 - s1);
    ta = (a - t1) / (t2 - t1);
    tb = (b - t1) / (t2 - t1);
    if ta > tb
        [ta, tb] = deal(tb, ta);
        [sa, sb] = deal(sb, sa);
        [r.ta, r.tb] = deal(r.tb, r.ta);
    end
    r.t0 = ta;
    r.t1 = tb;
    r.s0 = sa;
    r.s1 = sb;
    out(end + 1) = r; %#ok<AGROW>
end
[~, o] = sort([out.t0]);
out = pm_merge(out(o));
end

function pm = pm_merge(pm)
% join neighbouring identity pieces
k = 1;
while k < numel(pm)
    a = pm(k);
    b = pm(k + 1);
    if strcmp(a.kind, 'lin') && strcmp(b.kind, 'lin') && a.t0 == a.s0 && a.t1 == a.s1 && ...
            b.t0 == b.s0 && b.t1 == b.s1 && a.t1 == b.t0
        pm(k).t1 = b.t1;
        pm(k).s1 = b.s1;
        pm(k + 1) = [];
    else
        k = k + 1;
    end
end
end

function tf = pm_equal(a, b)
tf = numel(a) == numel(b);
for i = 1:numel(a)
    if ~tf
        return
    end
    tf = a(i).t0 == b(i).t0 && a(i).t1 == b(i).t1 && a(i).s0 == b(i).s0 && a(i).s1 == b(i).s1 && ...
        strcmp(a(i).kind, b(i).kind) && a(i).ha == b(i).ha && a(i).ta == b(i).ta && a(i).tb == b(i).tb;
end
end

% =============================================================== splits

function [P, general] = split_constant_z(P, general)
% every exact patch at the ends of each constant-z interval strictly inside its z range
out = P([]);
g = false(1, 0);
for k = 1:numel(P)
    e = P(k);
    if general(k) || e.z_range(1) == e.z_range(2)
        out(end + 1) = e; %#ok<AGROW>
        g(end + 1) = general(k); %#ok<AGROW>
        continue
    end
    s = e.surf;
    p = s.degree(1);
    Zr = s.ctrl(:, 1, 3);
    ku = unique(s.knots{1});
    cut = [];
    for j = 1:numel(ku) - 1
        span = find(s.knots{1} <= ku(j), 1, 'last');
        zz = Zr(span - p:span);
        if all(zz == zz(1)) && zz(1) > min(Zr) && zz(1) < max(Zr)
            cut = [cut, ku(j), ku(j + 1)]; %#ok<AGROW>
        end
    end
    cut = unique(cut);
    cut = cut(cut > ku(1) & cut < ku(end));
    pieces = split_at_knots(e, cut);
    out = [out, pieces]; %#ok<AGROW>
    g = [g, false(1, numel(pieces))]; %#ok<AGROW>
end
P = out;
general = g;
end

function pieces = split_at_knots(e, cut)
pieces = e([]);
rest = e;
for x = cut
    s = rest.surf;
    p = s.degree(1);
    [ctrl, W, knots] = insert_rows(s.ctrl, s.weights, s.knots{1}, p, x, p - sum(s.knots{1} == x));
    f = find(knots == x, 1);
    lo = rest;
    hi = rest;
    lo.surf.ctrl = ctrl(1:f - 1, :, :);
    hi.surf.ctrl = ctrl(f - 1:end, :, :);
    lo.surf.knots{1} = [knots(1:f - 1), repmat(x, 1, p + 1)];
    hi.surf.knots{1} = [repmat(x, 1, p + 1), knots(f + p:end)];
    if ~isempty(W)
        lo.surf.weights = W(1:f - 1, :);
        hi.surf.weights = W(f - 1:end, :);
    end
    pieces(end + 1) = refresh(lo); %#ok<AGROW>
    rest = refresh(hi);
end
pieces(end + 1) = rest;
end

function [P, general] = split_corner_heights(P, general)
% every exact patch at every patch-corner height strictly inside its z range (F3b); new row end
% points are identified with the corners and with each other at that height
has = arrayfun(@(e) ~isempty(e.surf), P);
H = zeros(0, 3);
for k = find(has)
    H = [H; corners(P(k).surf)]; %#ok<AGROW>
end
if isempty(H)
    return
end
hz = unique(H(:, 3))';
out = P([]);
g = false(1, 0);
% new row ends: entry, column, height, |C_u| eps(u*) of the column at the cut, column degree
newends = zeros(0, 5);
for k = 1:numel(P)
    e = P(k);
    if general(k) || e.z_range(1) == e.z_range(2)
        out(end + 1) = e; %#ok<AGROW>
        g(end + 1) = general(k); %#ok<AGROW>
        continue
    end
    zlo = min(e.z_range);
    zhi = max(e.z_range);
    hs = hz(hz > zlo & hz < zhi);
    if e.z_range(1) > e.z_range(2)
        hs = fliplr(hs);
    end
    rest = e;
    for h = hs
        [lo, hi, us] = mwecmass.solid.split_bspline_surface(rest, h);
        out(end + 1) = refresh(lo); %#ok<AGROW>
        g(end + 1) = false; %#ok<AGROW>
        for b = [1 3]
            c = boundary_curve(rest.surf, b);
            [~, Cu] = mwecmass.solid.eval_bspline_curve(c, us);
            j = 1;
            if b == 3
                j = size(lo.surf.ctrl, 2);
            end
            newends = [newends; numel(out), j, h, norm(Cu(1:2)) * eps(us), c.degree]; %#ok<AGROW>
        end
        rest = refresh(hi);
    end
    out(end + 1) = rest; %#ok<AGROW>
    g(end + 1) = false; %#ok<AGROW>
end
P = out;
general = g;
% identification: an end point at a corner of another patch takes that corner; end points of two
% new rows at one point take the first one's value
for q = 1:size(newends, 1)
    k = newends(q, 1);
    j = newends(q, 2);
    h = newends(q, 3);
    Eq = reshape(P(k).surf.ctrl(end, j, :), 1, 3);
    C = H(H(:, 3) == h, :);
    target = [];
    if ~isempty(C)
        d = sqrt(sum((C - Eq).^2, 2));
        [dm, i] = min(d);
        if dm <= rounding_bound(max(abs([Eq C(i, :)])), newends(q, 5), newends(q, 4))
            target = C(i, :);
        end
    end
    if isempty(target)
        for r = 1:q - 1
            if newends(r, 3) ~= h
                continue
            end
            E2 = reshape(P(newends(r, 1)).surf.ctrl(end, newends(r, 2), :), 1, 3);
            if norm(E2 - Eq) <= rounding_bound(max(abs([Eq E2])), max(newends([q r], 5)), ...
                    newends(q, 4) + newends(r, 4))
                target = E2;
                break
            end
        end
    end
    if ~isempty(target)
        P(k).surf.ctrl(end, j, :) = reshape(target, 1, 1, 3);
        % the next piece of the same patch starts with this row
        P(k + 1).surf.ctrl(1, j, :) = reshape(target, 1, 1, 3);
    end
end
end

function b = rounding_bound(m, p, moved)
% Largest distance between two computations of one exact point of a converted curve of degree p
% with control points of magnitude <= m: each is a corner (deck data, at most 4 roundings: the
% decimal, the revolution offset and two products), a curve evaluation or a cut by p knot
% insertions (Piegl & Tiller A5.1: per coordinate p + 1 weighted terms, each at most 6 rounded
% operations: the factor, two products, the weight, the sum and the division), so each coordinate
% of each point is off by at most 6 (p + 1) eps(m), the distance by 12 sqrt(3) (p + 1) eps(m) for
% the two; moved: the shift of the point by its parameter (|C_u| eps(u) for a parameter found to
% adjacent doubles).
b = 12 * sqrt(3) * (p + 1) * eps(m) + moved;
end

function C = corners(s)
C = [reshape(s.ctrl(1, 1, :), 1, 3); reshape(s.ctrl(1, end, :), 1, 3); ...
    reshape(s.ctrl(end, 1, :), 1, 3); reshape(s.ctrl(end, end, :), 1, 3)];
end

% =============================================================== seams

function [P, general] = find_seams(P, general, unparsed)
% Contract F1 order of the rules, repeated until nothing changes: (1a) constant-z patches with a
% vertex inside one of their rows (every boundary of a constant-z patch is a row), or with a corner
% inside a row of another patch, join a flat region with every constant-z patch that touches it at
% that height; (1b) vertices counted with the end points of seams between merged patches dropped:
% a patch with a vertex inside one of its rows takes the general path; (2) a boundary shared end to
% end (the same two end points) that is not the same curve sends the later patch in visible order
% to the general path. Vertices here are the corners of the entries.
fields = {'seam_v0', 'seam_u1', 'seam_v1', 'seam_u0'};
n = numel(P);
has = arrayfun(@(e) ~isempty(e.surf), P);
flatz = false(1, n);
for k = find(has)
    flatz(k) = P(k).z_range(1) == P(k).z_range(2);
end
B = cell(n, 4);
for k = find(has)
    for b = 1:4
        B{k, b} = boundary_curve(P(k).surf, b);
    end
end
C = cell(1, n);
for k = find(has)
    C{k} = corners(P(k).surf);
end
changed = true;
while changed
    changed = false;
    % (1a) flat regions
    merged = general & flatz;
    for k = find(~general & flatz)
        if any_inside(C, B, k, row_set(flatz, k), setdiff(find(has), k), []) || ...
                corner_in_other_row(C, B, flatz, has, k)
            merged(k) = true;
        end
    end
    grow = true;
    while grow
        grow = false;
        for k = find(~merged & flatz & ~general)
            for j = find(merged)
                if P(j).z_range(1) == P(k).z_range(1) && touches(C, B, k, j)
                    merged(k) = true;
                    grow = true;
                    break
                end
            end
        end
    end
    for k = find(merged & ~general)
        general = mark_piece_general(P, general, k);
        changed = true;
    end
    if changed
        continue
    end
    % (1b) vertices inside rows of the exact entries, seam ends between merged patches dropped
    drop = dropped_corners(C, B, P, merged);
    for k = find(~general & has)
        others = setdiff(find(has), k);
        if any_inside(C, B, k, row_set(flatz, k), others, drop)
            general = mark_general(P, general, k);
            changed = true;
        end
    end
    if changed
        continue
    end
    % (2) seams among the exact entries; tie-break for boundaries shared end to end
    ex = find(~general & has);
    for k = 1:n
        for b = 1:4
            P(k).(fields{b}) = [];
        end
    end
    for k = ex
        for b = 1:4
            ck = B{k, b};
            if collapsed(ck.ctrl) || ~isempty(P(k).(fields{b}))
                continue
            end
            for j = ex
                hit = false;
                for c = 1:4
                    if (j == k && c == b) || ~isempty(P(j).(fields{c}))
                        continue
                    end
                    if same_ends(ck, B{j, c}) && same_curve(ck, B{j, c})
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
    open = zeros(0, 2);
    for k = ex
        for b = 1:4
            if isempty(P(k).(fields{b})) && ~collapsed(B{k, b}.ctrl)
                open(end + 1, :) = [k b]; %#ok<AGROW>
            end
        end
    end
    for q = 1:size(open, 1)
        k = open(q, 1);
        ck = B{k, open(q, 2)};
        cand = open(open(:, 1) ~= k, :);
        keep = false(size(cand, 1), 1);
        for i = 1:size(cand, 1)
            keep(i) = same_ends(ck, B{cand(i, 1), cand(i, 2)});
        end
        cand = cand(keep, :);
        if isempty(cand)
            continue
        end
        % the partner among boundaries with the same two end points: the one whose midpoint is nearest
        mk = curve_mid(ck);
        dist = arrayfun(@(i) norm(curve_mid(B{cand(i, 1), cand(i, 2)}) - mk), 1:size(cand, 1));
        [~, i] = min(dist);
        j = cand(i, 1);
        later = j;
        if P(k).visible > P(j).visible
            later = k;
        end
        general = mark_general(P, general, later);
        changed = true;
        break
    end
    if changed
        continue
    end
    if ~isempty(open) && ~any(general)
        k = open(1, 1);
        if isempty(unparsed)
            unparsed = {'(none)'};
        end
        error('mwecmass:solid:HullNotClosed', ['outer_nurbs: boundary %d of %s has no neighbour and is not ' ...
            'a pole. Deck lines MS2Parser did not read:\n%s'], open(1, 2), P(k).name, ...
            strjoin(unparsed, sprintf('\n')));
    end
end
end

function rows = row_set(flatz, k)
% the row boundaries of entry k: u0 and u1, and every boundary of a constant-z entry
rows = [2 4];
if flatz(k)
    rows = 1:4;
end
end

function tf = any_inside(C, B, k, rows, others, drop)
% a corner of one of the entries others (not listed in drop, rows [entry corner]) lies strictly
% inside one of the boundaries rows of entry k
tf = false;
for b = rows
    c = B{k, b};
    if collapsed(c.ctrl)
        continue
    end
    for j = others
        for i = 1:4
            if ~isempty(drop) && any(drop(:, 1) == j & drop(:, 2) == i)
                continue
            end
            if strictly_inside(c, C{j}(i, :))
                tf = true;
                return
            end
        end
    end
end
end

function tf = corner_in_other_row(C, B, flatz, has, k)
% a corner of entry k lies strictly inside a row of another entry
tf = false;
for j = setdiff(find(has), k)
    for b = row_set(flatz, j)
        c = B{j, b};
        if collapsed(c.ctrl)
            continue
        end
        for i = 1:4
            if strictly_inside(c, C{k}(i, :))
                tf = true;
                return
            end
        end
    end
end
end

function tf = touches(C, B, k, j)
% two constant-z entries at one height touch: a corner of one lies on a boundary of the other
tf = false;
for pair = [k j; j k]'
    for i = 1:4
        X = C{pair(1)}(i, :);
        for b = 1:4
            if on_curve(B{pair(2), b}, X)
                tf = true;
                return
            end
        end
    end
end
end

function drop = dropped_corners(C, B, P, merged)
% corners of merged patches that are end points of seams between merged patches at one height:
% they lie on a boundary of another merged patch of the same height
drop = zeros(0, 2);
idx = find(merged);
for k = idx
    for i = 1:4
        X = C{k}(i, :);
        for j = idx
            if j == k || P(j).z_range(1) ~= P(k).z_range(1)
                continue
            end
            if any(arrayfun(@(b) on_curve(B{j, b}, X), 1:4))
                drop(end + 1, :) = [k i]; %#ok<AGROW>
                break
            end
        end
    end
end
end

function general = mark_general(P, general, k)
% a mirror takes its source's path: every entry of the same ultimate source goes along
src = P(k).source;
for j = 1:numel(P)
    if strcmp(P(j).source, src)
        general(j) = true;
    end
end
end

function general = mark_piece_general(P, general, k)
% one piece of a patch joins a flat region; the same piece of every mirror of its source goes along
vis = [P.visible];
grp = find(vis == P(k).visible);
i = find(grp == k);
for v = unique(vis(strcmp({P.source}, P(k).source)))
    g = find(vis == v);
    if numel(g) == numel(grp)
        general(g(i)) = true;
    end
end
end

function [d, x, c1] = project(c, X)
% distance from X to the curve c: Newton from the nearest of 401 samples, to adjacent doubles
s = linspace(c.knots(1), c.knots(end), 401)';
Q = mwecmass.solid.eval_bspline_curve(c, s);
[~, m] = min(sum((Q - X).^2, 2));
x = s(m);
for it = 1:60
    [Cx, C1, C2] = mwecmass.solid.eval_bspline_curve(c, x);
    H = C1 * C1' + (Cx - X) * C2';
    if H <= 0
        H = C1 * C1';
    end
    xn = min(max(x - ((Cx - X) * C1') / H, c.knots(1)), c.knots(end));
    if xn == x || ~isfinite(xn)
        break
    end
    x = xn;
end
[Cx, c1] = mwecmass.solid.eval_bspline_curve(c, x);
d = norm(Cx - X);
end

function tf = on_curve(c, X)
% X lies on the curve c (its end points included) up to the rounding of two computations of one
% point (rounding_bound)
if isequal(X, c.ctrl(1, :)) || isequal(X, c.ctrl(end, :))
    tf = true;
    return
end
m = max(abs([c.ctrl(:); X(:)]));
tf = false;
if outside_hull_box(c, X, rounding_bound(m, c.degree, 0))
    return
end
[d, x, c1] = project(c, X);
tf = d <= rounding_bound(m, c.degree, norm(c1) * eps(x));
end

function tf = strictly_inside(c, X)
% X lies on the curve c and is neither of its end points (rounding_bound for both tests)
tf = false;
if X(3) ~= c.ctrl(1, 3) && all(c.ctrl(:, 3) == c.ctrl(1, 3))
    return
end
if isequal(X, c.ctrl(1, :)) || isequal(X, c.ctrl(end, :))
    return
end
m = max(abs([c.ctrl(:); X(:)]));
bnd = rounding_bound(m, c.degree, 0);
if norm(X - c.ctrl(1, :)) <= bnd || norm(X - c.ctrl(end, :)) <= bnd || outside_hull_box(c, X, bnd)
    return
end
[d, x, c1] = project(c, X);
tf = x > c.knots(1) && x < c.knots(end) && d <= rounding_bound(m, c.degree, norm(c1) * eps(x));
end

function tf = outside_hull_box(c, X, bnd)
% the curve lies in the convex hull of its control points (positive weights): X is off the curve
% when it lies farther than bnd (the rounding of the two computations, rounding_bound) outside
% their bounding box
tf = any(X < min(c.ctrl, [], 1) - bnd | X > max(c.ctrl, [], 1) + bnd);
end

function X = curve_mid(c)
X = mwecmass.solid.eval_bspline_curve(c, (c.knots(1) + c.knots(end)) / 2);
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

function tf = same_ends(a, b)
tf = (isequal(a.ctrl(1, :), b.ctrl(1, :)) && isequal(a.ctrl(end, :), b.ctrl(end, :))) || ...
    (isequal(a.ctrl(1, :), b.ctrl(end, :)) && isequal(a.ctrl(end, :), b.ctrl(1, :)));
end

function tf = same_curve(a, b)
% contract S1: same degree and number of knots, bitwise control points and weights (empty = ones),
% knots equal after the affine map to [0, 1] up to 4 * 2^-52 (three correctly rounded operations
% on a value in [0, 1]); either direction
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
kr = fliplr((b.knots(end) - b.knots) / (b.knots(end) - b.knots(1)));
tf = isequal(a.ctrl, flipud(b.ctrl)) && isequal(wa(:), flipud(wb(:))) && all(abs(ka - kr) <= 4 * 2^-52);
end

% =============================================================== orientation

function P = set_outward(P)
% outward = S_u x S_v points out of the hull. Lateral patch: at a height inside its z range, the
% horizontal part of the normal against the outward normal of the counter-clockwise section loop
% (F4) at the patch's row. Constant-z patch: the hull lies below it when its centre is inside the
% section of the faces that reach its height from below and not inside the one from above.
lat = arrayfun(@(e) ~isempty(e.surf) && e.z_range(1) ~= e.z_range(2), P);
for k = 1:numel(P)
    e = P(k);
    if isempty(e.surf)
        continue
    end
    if lat(k)
        P(k).outward = lateral_outward(P(lat), find(find(lat) == k), e);
    else
        P(k).outward = flat_outward(P(lat), e);
    end
end
end

function out = lateral_outward(L, kk, e)
for frac = [0.5 1/3 2/3 0.25 0.75]
    h = e.z_range(1) + frac * (e.z_range(2) - e.z_range(1));
    loop = mwecmass.solid.slice_bspline_surface(L, h);
    i = find([loop.pieces.patch] == kk, 1);
    if isempty(i)
        continue
    end
    pc = loop.pieces(i);
    kv = pc.curve.knots;
    for v = kv(1) + [0.5 0.3 0.7] * (kv(end) - kv(1))
        [~, T] = mwecmass.solid.eval_bspline_curve(pc.curve, v);
        T = pc.dir * T;
        o = [T(2), -T(1)];
        [~, Su, Sv] = mwecmass.solid.eval_bspline_surface(e.surf, pc.u, v);
        n = cross(Su, Sv);
        d = n(1:2) * o';
        if d ~= 0
            out = d > 0;
            return
        end
    end
end
error('mwecmass:solid:OrientationUndecided', 'outer_nurbs: orientation of %s is undecided', e.name);
end

function out = flat_outward(L, e)
h = e.z_range(1);
ku = e.surf.knots{1};
kv = e.surf.knots{2};
[S, Su, Sv] = mwecmass.solid.eval_bspline_surface(e.surf, mean(ku([1 end])), mean(kv([1 end])));
n = cross(Su, Sv);
zr = reshape([L.z_range], 2, [])';
below = L(min(zr, [], 2) < h & max(zr, [], 2) >= h);
above = L(min(zr, [], 2) <= h & max(zr, [], 2) > h);
in_b = inside(below, h, S);
in_a = inside(above, h, S);
if in_b == in_a || n(3) == 0
    error('mwecmass:solid:OrientationUndecided', 'outer_nurbs: orientation of the constant-z patch %s is undecided', e.name);
end
out = (n(3) > 0) == in_b;
end

function tf = inside(L, h, S)
tf = false;
if isempty(L)
    return
end
loop = mwecmass.solid.slice_bspline_surface(L, h);
tf = inpolygon(S(1), S(2), loop.pts(:, 1), loop.pts(:, 2));
end
