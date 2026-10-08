function body = build_body(geo, design, inner)
%BUILD_BODY  B-rep of a realised design: modules, shells, voids and ballast (contract F5, S3, S4).
%
%   body = mwecmass.solid.build_body(geo, design, inner)
%
%   geo: S1b (outer_nurbs). design: S3 (mode, edges, vs, t, z_ballast, solid_modules). inner: S2
%   sets (offset_surface), one per distinct t of the hollow modules; module i uses the set whose t
%   equals design.t(i) bitwise.
%
%   Layout (rule 11): in module i the full outer section is solid below z_ballast; above it a
%   hollow module is a shell around the void of its set; solid_modules are solid throughout.
%   Regions: uhpc and air (modular_precast), ballast, shell and air (thin_shell).
%
%   Lateral faces are the outer and inner pieces split at their c0_u and c0_v rows (existing
%   knots of full multiplicity: no control point moves) and cut by split_bspline_surface at the
%   planes (interior module edges and z_ballast) strictly inside their z range. Inner pieces are
%   kept only inside the void of a module that uses their set. Constant-z pieces are faces as
%   they are: outer ones always; an inner one unless a plane face takes its place (at z_ballast,
%   and in modular precast at a module edge). A collapsed boundary row (pole) is left out of the
%   face's loop.
%
%   Plane faces: at each plane the loops of the lateral faces ending there from below and from
%   above (outer and, where the module on that side has a void there, inner) are chains of the
%   rows already built. Every edge of these loops is classified on both sides by the loops it
%   belongs to and, for the others, by the winding number of an interior point of the edge on the
%   exact curves (Bezier subdivision until the point lies outside the control-point box). A plane
%   face covers the area where the regions below and above differ as S4 requires (cap,
%   joint_step, ballast_top; a cap also between the uhpc of two precast modules); its loops are
%   the edges where that face type begins or ends, one face per outer loop with its holes. Where
%   a face of one side is outside the hull on the other side, the flat part is already an outer
%   (or inner) constant-z face. At a joint_step the loop of the larger t must lie inside the loop
%   of the smaller t: no crossing of the two exact loops (Bezier pieces separated by their
%   control-point boxes or fat lines) and one point of the larger-t loop inside the other, else
%   mwecmass:solid:JointNotNested.
%
%   Edges: one edge per shared boundary, taken from the first face that builds it; a later face
%   boundary that is the same curve (contract S1: degree, control points and weights bitwise, an
%   empty weight vector equal to all ones, knots equal after the map to [0, 1]; either direction)
%   references it. Vertices are bitwise-equal points (-0 = +0).
%
%   body: design, inner_t, planes, brep (T9 writer struct; extra face fields role, module, inside,
%   outside), shells (precast all_outer, all_voids; thin shell layer, void) and voids(i): the void
%   of module i, z_lo and z_hi [m, body] (NaN: none) and open_lo, open_hi (true where the void
%   ends in a pole or line of its set, so its section has no area there).
%   Flat regions (geo.flat, inner(k).flat; general path): one plane face each at its height, normal
%   [0 0 normal_z], whose loops are the chains of the rows of lateral faces that name it (S1 seam
%   [0 j]); an outer flat region is an outer face of the module on its material side, an inner one
%   a face of the void kept by the rules of constant-z inner pieces.
%   Errors: mwecmass:solid:BadEdges, MissingInnerSet, SectionNotClosed, JointNotNested.

e = design.edges(:);
N = numel(e) - 1;
zmin = geo.z_range(1);
zmax = geo.z_range(2);
if N < 1 || e(1) ~= zmin || e(end) ~= zmax || any(diff(e) <= 0)
    error('mwecmass:solid:BadEdges', 'build_body: module edges must rise from z_min to z_max of the hull');
end
precast = strcmp(design.mode, 'modular_precast');
if precast
    solid_name = 'uhpc';
else
    solid_name = 'shell';
end
zb = design.z_ballast;
t = design.t(:);
set_of = zeros(N, 1);
for i = 1:N
    if any(design.solid_modules == i) || ~isfinite(t(i))
        continue
    end
    k = find(arrayfun(@(s) isequal(s.t, t(i)), inner), 1);
    if isempty(k)
        error('mwecmass:solid:MissingInnerSet', 'build_body: no inner set for t = %.17g (module %d)', t(i), i);
    end
    set_of(i) = k;
end
planes = e(2:end - 1);
if zb > zmin && zb < zmax
    planes = [planes; zb];
end
planes = unique(planes);

% extent of each set and whether it ends in a constant-z area
nset = numel(inner);
sinfo = struct('z_lo', cell(1, nset), 'z_hi', [], 'flat_lo', [], 'flat_hi', []);
for k = unique(set_of(set_of > 0))'
    zr = reshape([inner(k).patches.z_range], 2, [])';
    lo = min(zr(:));
    hi = max(zr(:));
    cz = zr(:, 1) == zr(:, 2);
    sinfo(k).z_lo = lo;
    sinfo(k).z_hi = hi;
    fz = [];
    if isfield(inner(k), 'flat') && ~isempty(inner(k).flat)
        fz = [inner(k).flat.z];
    end
    sinfo(k).flat_lo = any(cz & zr(:, 1) == lo) || any(fz == lo);
    sinfo(k).flat_hi = any(cz & zr(:, 1) == hi) || any(fz == hi);
end
% void of each module: (lo, hi) where its set leaves air above z_ballast
voids = struct('z_lo', num2cell(NaN(N, 1)), 'z_hi', NaN, 'open_lo', false, 'open_hi', false);
for i = 1:N
    k = set_of(i);
    if k == 0
        continue
    end
    lo = max([e(i), zb, sinfo(k).z_lo]);
    hi = min(e(i + 1), sinfo(k).z_hi);
    if lo < hi
        voids(i).z_lo = lo;
        voids(i).z_hi = hi;
        voids(i).open_lo = lo == sinfo(k).z_lo && ~sinfo(k).flat_lo;
        voids(i).open_hi = hi == sinfo(k).z_hi && ~sinfo(k).flat_hi;
    end
end

brep = struct('vertices', zeros(0, 3), ...
    'curves', struct('degree', {}, 'ctrl', {}, 'knots', {}, 'weights', {}), ...
    'edges', struct('vertices', {}, 'curve', {}), 'surfaces', {{}}, ...
    'faces', struct('surface', {}, 'same_sense', {}, 'loops', {}, 'role', {}, 'module', {}, ...
    'inside', {}, 'outside', {}), ...
    'bodies', struct('name', {}, 'kind', {}, 'shells', {}));
vmap = containers.Map('KeyType', 'char', 'ValueType', 'double');
emap = containers.Map('KeyType', 'char', 'ValueType', 'any');
% rows of lateral faces at their lowest and highest height: set (0 outer), module, z, edge (0: collapsed)
rows = struct('set', {}, 'module', {}, 'z', {}, 'edge', {}, 'top', {});
% rows that bound a flat region j (S1 seam [0 j]) of the outer surface (set 0) or of inner set k
frows = struct('set', {}, 'flat', {}, 'edge', {});

% outer faces
for p = 1:numel(geo.outer)
    for piece = c0_pieces(geo.outer(p))
        zr = piece.z_range;
        if zr(1) == zr(2)
            h = zr(1);
            up = flat_up(piece);
            if up
                m = module_below(h);
                ballast = h <= zb;
            else
                m = module_above(h);
                ballast = h < zb;
            end
            add_face(piece.surf, piece.outward, 'outer', m, region(ballast), 'exterior', 0);
            continue
        end
        for slab = cut_at(piece, planes)
            zm = mean(slab.z_range);
            add_face(slab.surf, slab.outward, 'outer', module_at(zm), region(zm < zb), 'exterior', 0, ...
                flat_ids(geo.outer(p), slab));
        end
    end
end

% inner faces: each piece of a used set cut once at every plane (so that the rows two modules
% share at a joint come from one knot insertion), kept inside the void of a module using that set
for k = unique(set_of(set_of > 0))'
    for q = 1:numel(inner(k).patches)
        for piece = c0_pieces(inner(k).patches(q))
            zr = piece.z_range;
            if zr(1) == zr(2)
                [keep, i] = keep_inner_flat(zr(1), void_above(piece), k);
                if keep
                    add_inner(piece, i, k);
                end
                continue
            end
            for slab = cut_at(piece, planes)
                zm = mean(slab.z_range);
                i = module_at(zm);
                if set_of(i) == k && zm > voids(i).z_lo && zm < voids(i).z_hi
                    add_inner(slab, i, k, flat_ids(inner(k).patches(q), slab));
                end
            end
        end
    end
end

% flat regions (general path): one plane face each, bounded by the rows that name it
gflat = [];
if isfield(geo, 'flat')
    gflat = geo.flat;
end
for j = 1:numel(gflat)
    h = gflat(j).z;
    if gflat(j).normal_z > 0
        add_flat(0, j, h, 1, 'outer', module_below(h), region(h <= zb), 'exterior');
    else
        add_flat(0, j, h, -1, 'outer', module_above(h), region(h < zb), 'exterior');
    end
end
for k = unique(set_of(set_of > 0))'
    if ~isfield(inner(k), 'flat')
        continue
    end
    for j = 1:numel(inner(k).flat)
        h = inner(k).flat(j).z;
        [keep, i] = keep_inner_flat(h, inner(k).flat(j).normal_z > 0, k);
        if keep
            add_flat(k, j, h, sign(inner(k).flat(j).normal_z), 'inner', i, solid_name, 'air');
        end
    end
end

% plane faces
for zp = planes'
    plane_faces(zp);
end

% bodies and shells
hull = geo.hull_name;
if precast
    for i = 1:N
        brep.bodies(end + 1) = struct('name', sprintf('%s_UHPC_module_%d', hull, i), 'kind', 'solid', ...
            'shells', {solid_shells(bound('uhpc', i))});
    end
    shells = struct('all_outer', -bound('exterior', []), 'all_voids', {components(bound('air', []))});
else
    if ~isempty(bound('ballast', []))
        brep.bodies(end + 1) = struct('name', sprintf('%s_STEEL_ballast', hull), 'kind', 'solid', ...
            'shells', {solid_shells(bound('ballast', []))});
    end
    out = -bound('exterior', []);
    sheet = out(arrayfun(@(f) any(strcmp({brep.faces(abs(f)).inside, brep.faces(abs(f)).outside}, 'shell')), out));
    if ~isempty(sheet)
        brep.bodies(end + 1) = struct('name', sprintf('%s_STEEL_shell', hull), 'kind', 'sheet', 'shells', {{sheet}});
    end
    shells = struct('layer', bound('shell', []), 'void', bound('air', []));
end
used = unique(set_of(set_of > 0 & isfinite([voids.z_lo]')));
inner_t = [];
if ~isempty(used)
    inner_t = sort([inner(used).t]);
end
body = struct('design', design, 'inner_t', inner_t, 'planes', planes(:)', 'brep', brep, ...
    'shells', shells, 'voids', voids);

% ------------------------------------------------------------ nested helpers

    function m = module_at(z)
        m = find(e(1:end - 1) <= z & z < e(2:end), 1);
        if isempty(m)
            m = N;
        end
    end

    function m = module_below(z)
        m = find(e(1:end - 1) < z & z <= e(2:end), 1);
        if isempty(m)
            m = 1;
        end
    end

    function m = module_above(z)
        m = find(e(1:end - 1) <= z & z < e(2:end), 1);
        if isempty(m)
            m = N;
        end
    end

    function r = region(ballast)
        if ~precast && ballast
            r = 'ballast';
        else
            r = solid_name;
        end
    end

    function [keep, i] = keep_inner_flat(h, up, k)
        % a constant-z part of inner set k at height h (up: the void lies above it) is a face of the
        % void of module i unless a plane face takes its place (z_ballast; precast module edges)
        if up
            i = module_above(h);
            keep = set_of(i) == k && h >= voids(i).z_lo && h < voids(i).z_hi && ...
                ~(h == voids(i).z_lo && (h == zb || (precast && h == e(i))));
        else
            i = module_below(h);
            keep = set_of(i) == k && h > voids(i).z_lo && h <= voids(i).z_hi && ...
                ~(precast && h == e(i + 1));
        end
    end

    function add_inner(piece, i, k, fl)
        if nargin < 4
            fl = [0 0];
        end
        if piece.outward
            add_face(piece.surf, true, 'inner', i, solid_name, 'air', k, fl);
        else
            add_face(piece.surf, true, 'inner', i, 'air', solid_name, k, fl);
        end
    end

    function add_flat(set, j, h, nz, role, m, inside, outside)
        % plane face of flat region j at height h, normal [0 0 nz]: its loops are the chains of the
        % rows that name it; a loop inside an even number of the others is an outer loop
        % (counter-clockwise seen from the normal), the others are holes of the smallest outer loop
        % around them
        ed = unique([frows([frows.set] == set & [frows.flat] == j).edge]);
        if isempty(ed)
            error('mwecmass:solid:SectionNotClosed', 'build_body: no rows bound flat region %d at z = %.17g', j, h);
        end
        loops = chain_all(ed, h, true);
        n = numel(loops);
        A = cellfun(@(Lp) loop_area(Lp), loops);
        depth = zeros(1, n);
        pts = cell(1, n);
        for a = 1:n
            pts{a} = curve_midpoint(brep.curves(brep.edges(abs(loops{a}(1))).curve));
        end
        for a = 1:n
            for b = 1:n
                if a ~= b && winding(loops{b}, pts{a}) ~= 0
                    depth(a) = depth(a) + 1;
                end
            end
        end
        outer = find(mod(depth, 2) == 0);
        for a = 1:n
            want = nz * (2 * (mod(depth(a), 2) == 0) - 1);
            if sign(A(a)) ~= want
                loops{a} = -loops{a}(end:-1:1);
            end
        end
        for o = outer
            holes = find(depth == depth(o) + 1);
            holes = holes(arrayfun(@(b) winding(loops{o}, pts{b}) ~= 0, holes));
            brep.surfaces{end + 1} = struct('type', 'plane', 'origin', [0 0 h], 'normal', [0 0 nz]);
            brep.faces(end + 1) = struct('surface', numel(brep.surfaces), 'same_sense', true, ...
                'loops', {[loops(o), loops(holes)]}, 'role', role, 'module', m, 'inside', inside, 'outside', outside);
        end
    end

    function add_face(s, outward, role, m, inside, other, set, fl)
        % lateral or constant-z face of a z_of_u piece; outward: S_u x S_v points from inside to other;
        % fl: the flat regions its u0 and u1 rows bound (0: none)
        [nu, nv, ~] = size(s.ctrl);
        W = s.weights;
        bnd = {{s.ctrl(:, 1, :), s.knots{1}, col(W, 1), s.degree(1), 1}, ...
               {s.ctrl(nu, :, :), s.knots{2}, row(W, nu), s.degree(2), 1}, ...
               {s.ctrl(:, nv, :), s.knots{1}, col(W, nv), s.degree(1), -1}, ...
               {s.ctrl(1, :, :), s.knots{2}, row(W, 1), s.degree(2), -1}};
        loop = zeros(1, 0);
        bedge = zeros(1, 4);
        for b = 1:4
            C = reshape(bnd{b}{1}, [], 3);
            if all(all(C == C(1, :)))
                continue
            end
            crv = struct('degree', bnd{b}{4}, 'ctrl', C, 'knots', bnd{b}{2}(:)', 'weights', bnd{b}{3});
            bedge(b) = edge_for(crv);
            loop(end + 1) = bnd{b}{5} * bedge(b); %#ok<AGROW>
        end
        if outward
            ins = inside; outs = other;
        else
            ins = other; outs = inside;
        end
        brep.surfaces{end + 1} = struct('type', 'bspline', 'degree', s.degree, 'ctrl', s.ctrl, ...
            'knots', {s.knots}, 'weights', s.weights);
        brep.faces(end + 1) = struct('surface', numel(brep.surfaces), 'same_sense', true, ...
            'loops', {{loop}}, 'role', role, 'module', m, 'inside', ins, 'outside', outs);
        z0 = s.ctrl(1, 1, 3);
        z1 = s.ctrl(nu, 1, 3);
        if z0 ~= z1
            % row 4 (u0) and row 2 (u1): the face's rows at its two end heights
            rows(end + 1) = struct('set', set, 'module', m, 'z', z0, 'edge', bedge(4), 'top', z0 > z1);
            rows(end + 1) = struct('set', set, 'module', m, 'z', z1, 'edge', bedge(2), 'top', z1 > z0);
            if nargin >= 8
                ub = [4 2];
                for r = find(fl ~= 0 & bedge(ub) ~= 0)
                    frows(end + 1) = struct('set', set, 'flat', fl(r), 'edge', bedge(ub(r))); %#ok<AGROW>
                end
            end
        end
    end

    function idx = edge_for(crv)
        a = vertex(crv.ctrl(1, :));
        b = vertex(crv.ctrl(end, :));
        key = sprintf('%d_%d', min(a, b), max(a, b));
        if emap.isKey(key)
            for kk = emap(key)
                s = same_curve(brep.curves(brep.edges(kk).curve), crv);
                if s ~= 0
                    idx = s * kk;
                    return
                end
            end
        else
            emap(key) = zeros(1, 0);
        end
        brep.curves(end + 1) = crv;
        brep.edges(end + 1) = struct('vertices', [a b], 'curve', numel(brep.curves));
        idx = numel(brep.edges);
        emap(key) = [emap(key), idx];
    end

    function i = vertex(xyz)
        key = sprintf('%s', num2hex(xyz + 0)');
        if vmap.isKey(key)
            i = vmap(key);
        else
            brep.vertices(end + 1, :) = xyz + 0;
            i = size(brep.vertices, 1);
            vmap(key) = i;
        end
    end

    function plane_faces(zp)
        mb = module_below(zp);
        ma = module_above(zp);
        Ob = section_loop(0, [], zp, true);
        Oa = section_loop(0, [], zp, false);
        Ib = [];
        Ia = [];
        if set_of(mb) > 0 && voids(mb).z_hi == zp && ~voids(mb).open_hi
            Ib = section_loop(set_of(mb), mb, zp, true);
        end
        if set_of(ma) > 0 && voids(ma).z_lo == zp && ~voids(ma).open_lo
            Ia = section_loop(set_of(ma), ma, zp, false);
        end
        if ~isempty(Ib) && ~isempty(Ia) && set_of(mb) ~= set_of(ma)
            if t(mb) < t(ma)
                check_nested(Ia, Ib, zp);
            else
                check_nested(Ib, Ia, zp);
            end
        end
        L = {Ob, Oa, Ib, Ia};
        present = ~cellfun(@isempty, L);
        E = unique(abs([L{:}]));
        keyL = cell(1, numel(E));
        keyR = cell(1, numel(E));
        for j = 1:numel(E)
            inL = false(1, 4);
            inR = false(1, 4);
            pt = [];
            for q = find(present)
                s = sign(L{q}(abs(L{q}) == E(j)));
                if ~isempty(s)
                    inL(q) = s(1) > 0;
                    inR(q) = s(1) < 0;
                else
                    if isempty(pt)
                        pt = curve_midpoint(brep.curves(brep.edges(E(j)).curve));
                    end
                    w = winding(L{q}, pt) ~= 0;
                    inL(q) = w;
                    inR(q) = w;
                end
            end
            keyL{j} = classify(inL, zp, mb, ma, present);
            keyR{j} = classify(inR, zp, mb, ma, present);
        end
        keys = unique([keyL(~cellfun(@isempty, keyL)), keyR(~cellfun(@isempty, keyR))]);
        for q = 1:numel(keys)
            K = keys{q};
            sel = zeros(1, 0);
            for j = 1:numel(E)
                l = strcmp(keyL{j}, K);
                r = strcmp(keyR{j}, K);
                if l && ~r
                    sel(end + 1) = E(j); %#ok<AGROW>
                elseif r && ~l
                    sel(end + 1) = -E(j); %#ok<AGROW>
                end
            end
            add_plane_faces(zp, K, sel);
        end
    end

    function key = classify(in, zp, mb, ma, present)
        % face type between the region just below and just above a point of the plane ('' none)
        key = '';
        if ~in(1) || ~in(2)
            return
        end
        if present(3) && in(3)
            rb = 'air';
        elseif precast
            rb = 'uhpc';
        else
            rb = region(zp <= zb);
        end
        if present(4) && in(4)
            ra = 'air';
        elseif precast
            ra = 'uhpc';
        else
            ra = region(zp < zb);
        end
        if strcmp(rb, ra)
            if precast && strcmp(rb, 'uhpc') && mb ~= ma
                role = 'cap';
            else
                return
            end
        elseif ~precast
            if zp ~= zb
                return
            end
            role = 'ballast_top';
        elseif strcmp(ra, 'air')
            if present(3)
                role = 'joint_step';
            elseif zp == zb
                role = 'ballast_top';
            else
                role = 'cap';
            end
        else
            if present(4)
                role = 'joint_step';
            else
                role = 'cap';
            end
        end
        key = sprintf('%s|%s|%s|%d|%d', role, rb, ra, mb, ma);
    end

    function add_plane_faces(zp, K, sel)
        parts = strsplit(K, '|');
        mb = str2double(parts{4});
        ma = str2double(parts{5});
        m = mb;
        if ma ~= mb
            m = [mb ma];
        end
        loops = chain_all(sel, zp);
        A = cellfun(@(Lp) loop_area(Lp), loops);
        outer = find(A > 0);
        holes = find(A < 0);
        owner = zeros(size(holes));
        for h = 1:numel(holes)
            pt = curve_midpoint(brep.curves(brep.edges(abs(loops{holes(h)}(1))).curve));
            best = Inf;
            for o = outer
                if A(o) < best && winding(loops{o}, pt) ~= 0
                    best = A(o);
                    owner(h) = o;
                end
            end
            if owner(h) == 0
                error('mwecmass:solid:SectionNotClosed', 'build_body: a hole loop at z = %.17g lies in no outer loop', zp);
            end
        end
        for o = outer
            lps = [loops(o), loops(holes(owner == o))];
            brep.surfaces{end + 1} = struct('type', 'plane', 'origin', [0 0 zp], 'normal', [0 0 1]);
            brep.faces(end + 1) = struct('surface', numel(brep.surfaces), 'same_sense', true, ...
                'loops', {lps}, 'role', parts{1}, 'module', m, 'inside', parts{2}, 'outside', parts{3});
        end
    end

    function L = section_loop(set, m, zp, below)
        % signed edges of the closed loop of the rows at zp of the lateral faces on one side,
        % counter-clockwise seen from +z; [] when every row there is collapsed
        pick = arrayfun(@(r) r.set == set && r.z == zp && r.top == below && ...
            (set == 0 || r.module == m), rows);
        ed = unique([rows(pick).edge]);
        ed = ed(ed ~= 0);
        if isempty(ed)
            L = [];
            return
        end
        loops = chain_all(abs(ed), zp, true);
        if numel(loops) ~= 1
            error('mwecmass:solid:SectionNotClosed', 'build_body: the rows at z = %.17g form %d loops', zp, numel(loops));
        end
        L = loops{1};
        if loop_area(L) < 0
            L = -L(end:-1:1);
        end
    end

    function loops = chain_all(sel, zp, free_direction)
        % closed chains of the signed edges sel (free_direction: edges may be reversed)
        if nargin < 3
            free_direction = false;
        end
        loops = {};
        rest = sel(:)';
        while ~isempty(rest)
            Lp = rest(1);
            rest(1) = [];
            [first, cur] = edge_ends(Lp);
            while cur ~= first
                hit = 0;
                for kk = 1:numel(rest)
                    [a, b] = edge_ends(rest(kk));
                    if a == cur
                        hit = kk;
                        s = rest(kk);
                        break
                    elseif free_direction && b == cur
                        hit = kk;
                        s = -rest(kk);
                        break
                    end
                end
                if hit == 0
                    error('mwecmass:solid:SectionNotClosed', 'build_body: the edges at z = %.17g do not close', zp);
                end
                rest(hit) = [];
                Lp(end + 1) = s; %#ok<AGROW>
                [~, cur] = edge_ends(s);
            end
            loops{end + 1} = Lp; %#ok<AGROW>
        end
    end

    function [a, b] = edge_ends(s)
        v = brep.edges(abs(s)).vertices;
        if s > 0
            a = v(1); b = v(2);
        else
            a = v(2); b = v(1);
        end
    end

    function A = loop_area(Lp)
        A = 0;
        for kk = 1:numel(Lp)
            A = A + sign(Lp(kk)) * curve_area(brep.curves(brep.edges(abs(Lp(kk))).curve));
        end
    end

    function w = winding(Lp, pt)
        ang = 0;
        for kk = 1:numel(Lp)
            c = brep.curves(brep.edges(abs(Lp(kk))).curve);
            ang = ang + sign(Lp(kk)) * curve_angle(c, pt);
        end
        w = round(ang / (2 * pi));
    end

    function check_nested(Lin, Lout, zp)
        % Lin (larger t) must lie inside Lout (smaller t): no crossing, one point inside
        A = loop_pieces(Lin);
        B = loop_pieces(Lout);
        if ~pieces_disjoint(A, B)
            error('mwecmass:solid:JointNotNested', ...
                'build_body: the inner loops at the joint z = %.17g cross or touch (no proof of separation)', zp);
        end
        pt = curve_midpoint(brep.curves(brep.edges(abs(Lin(1))).curve));
        if winding(Lout, pt) == 0
            error('mwecmass:solid:JointNotNested', ...
                'build_body: the void of the larger t at z = %.17g is not inside the void of the smaller t', zp);
        end
    end

    function P = loop_pieces(Lp)
        P = {};
        for kk = 1:numel(Lp)
            P = [P, bezier_pieces(brep.curves(brep.edges(abs(Lp(kk))).curve))]; %#ok<AGROW>
        end
    end

    function L = bound(rg, m)
        L = zeros(1, 0);
        for kk = 1:numel(brep.faces)
            f = brep.faces(kk);
            if strcmp(f.inside, rg) && (isempty(m) || f.module(1) == m)
                L(end + 1) = kk; %#ok<AGROW>
            elseif strcmp(f.outside, rg) && (isempty(m) || f.module(end) == m)
                L(end + 1) = -kk; %#ok<AGROW>
            end
        end
    end

    function C = components(L)
        % face lists connected through shared edges
        n = numel(L);
        lab = 1:n;
        Ef = cell(1, n);
        for kk = 1:n
            lp = brep.faces(abs(L(kk))).loops;
            Ef{kk} = abs([lp{:}]);
        end
        changed = true;
        while changed
            changed = false;
            for a = 1:n
                for b = a + 1:n
                    if lab(a) ~= lab(b) && ~isempty(intersect(Ef{a}, Ef{b}))
                        lab(lab == max(lab(a), lab(b))) = min(lab(a), lab(b));
                        changed = true;
                    end
                end
            end
        end
        C = {};
        for l = unique(lab)
            C{end + 1} = L(lab == l); %#ok<AGROW>
        end
    end

    function sh = solid_shells(L)
        % the component holding exterior faces first; further components are cavities, oriented
        % as standalone solids (normals out of the cavity)
        C = components(L);
        isout = false(1, numel(C));
        for c = 1:numel(C)
            for f = abs(C{c})
                isout(c) = isout(c) || strcmp(brep.faces(f).inside, 'exterior') || ...
                    strcmp(brep.faces(f).outside, 'exterior');
            end
        end
        sh = C(isout);
        for c = find(~isout)
            sh{end + 1} = -C{c}; %#ok<AGROW>
        end
    end
end

% ------------------------------------------------------------ pieces of z_of_u patches

function P = c0_pieces(patch)
% the patch split at its c0_u and c0_v rows (existing knots of full multiplicity)
P = patch;
for u = patch.c0_u(:)'
    Q = P([]);
    for k = 1:numel(P)
        s = P(k).surf;
        if u > s.knots{1}(1) && u < s.knots{1}(end)
            [a, b] = split_u_at_knot(P(k), u);
            Q = [Q, a, b]; %#ok<AGROW>
        else
            Q = [Q, P(k)]; %#ok<AGROW>
        end
    end
    P = Q;
end
for v = patch.c0_v(:)'
    Q = P([]);
    for k = 1:numel(P)
        s = P(k).surf;
        if v > s.knots{2}(1) && v < s.knots{2}(end)
            [a, b] = split_v_at_knot(P(k), v);
            Q = [Q, a, b]; %#ok<AGROW>
        else
            Q = [Q, P(k)]; %#ok<AGROW>
        end
    end
    P = Q;
end
end

function [lo, hi] = split_u_at_knot(patch, u)
s = patch.surf;
p = s.degree(1);
kn = s.knots{1}(:)';
f = find(kn == u, 1);
lo = patch;
hi = patch;
lo.surf.ctrl = s.ctrl(1:f - 1, :, :);
hi.surf.ctrl = s.ctrl(f - 1:end, :, :);
lo.surf.knots{1} = [kn(1:f - 1), repmat(u, 1, p + 1)];
hi.surf.knots{1} = [repmat(u, 1, p + 1), kn(f + p:end)];
if ~isempty(s.weights)
    lo.surf.weights = s.weights(1:f - 1, :);
    hi.surf.weights = s.weights(f - 1:end, :);
end
z = s.ctrl(f - 1, 1, 3);
lo.u_range = [patch.u_range(1) u];
hi.u_range = [u patch.u_range(2)];
lo.z_range = [patch.z_range(1) z];
hi.z_range = [z patch.z_range(2)];
lo.c0_u = patch.c0_u(patch.c0_u > lo.u_range(1) & patch.c0_u < u);
hi.c0_u = patch.c0_u(patch.c0_u > u & patch.c0_u < hi.u_range(2));
lo.pole = [patch.pole(1) collapsed(lo.surf.ctrl(end, :, :))];
hi.pole = [collapsed(hi.surf.ctrl(1, :, :)) patch.pole(2)];
end

function [lo, hi] = split_v_at_knot(patch, v)
s = patch.surf;
p = s.degree(2);
kn = s.knots{2}(:)';
f = find(kn == v, 1);
lo = patch;
hi = patch;
lo.surf.ctrl = s.ctrl(:, 1:f - 1, :);
hi.surf.ctrl = s.ctrl(:, f - 1:end, :);
lo.surf.knots{2} = [kn(1:f - 1), repmat(v, 1, p + 1)];
hi.surf.knots{2} = [repmat(v, 1, p + 1), kn(f + p:end)];
if ~isempty(s.weights)
    lo.surf.weights = s.weights(:, 1:f - 1);
    hi.surf.weights = s.weights(:, f - 1:end);
end
lo.c0_v = patch.c0_v(patch.c0_v < v);
hi.c0_v = patch.c0_v(patch.c0_v > v);
lo.pole = [collapsed(lo.surf.ctrl(1, :, :)) collapsed(lo.surf.ctrl(end, :, :))];
hi.pole = [collapsed(hi.surf.ctrl(1, :, :)) collapsed(hi.surf.ctrl(end, :, :))];
end

function P = cut_at(piece, zs)
% the lateral piece cut by split_bspline_surface at every height of zs strictly inside its z range,
% pieces in the order of u
zr = piece.z_range;
zs = zs(zs > min(zr) & zs < max(zr));
zs = unique(zs);
if zr(1) > zr(2)
    zs = flipud(zs(:));
end
P = piece([]);
rest = piece;
for z = zs(:)'
    [lo, hi] = mwecmass.solid.split_bspline_surface(rest, z);
    P(end + 1) = lo; %#ok<AGROW>
    rest = hi;
end
P(end + 1) = rest;
end

function fl = flat_ids(parent, piece)
% flat regions bounded by the u0 and u1 rows of piece, a piece of parent (S1 seam [0 j] on a row
% the piece keeps)
fl = [0 0];
f = {'seam_u0', 'seam_u1'};
for q = 1:2
    if ~isfield(parent, f{q})
        continue
    end
    sm = parent.(f{q});
    if numel(sm) == 2 && sm(1) == 0 && piece.z_range(q) == parent.z_range(q)
        fl(q) = sm(2);
    end
end
end

function up = flat_up(piece)
% outward normal of a constant-z outer piece points up (+z)
up = normal_z(piece.surf) > 0;
if ~piece.outward
    up = ~up;
end
end

function up = void_above(piece)
% the void lies above a constant-z inner piece (its normal into the void points up)
up = normal_z(piece.surf) > 0;
if ~piece.outward
    up = ~up;
end
end

function nz = normal_z(s)
% sign of the z component of S_u x S_v, summed over interior Gauss points of the knot spans
ku = unique(s.knots{1});
kv = unique(s.knots{2});
x = [0.2113248654051871; 0.7886751345948129];
uu = ku(1:end - 1)' + x' .* diff(ku)';
vv = kv(1:end - 1)' + x' .* diff(kv)';
[U, V] = ndgrid(uu(:), vv(:));
[~, Su, Sv] = mwecmass.solid.eval_bspline_surface(s, U(:), V(:));
nz = sum(Su(:, 1) .* Sv(:, 2) - Su(:, 2) .* Sv(:, 1));
if nz == 0
    error('mwecmass:solid:DegenerateFace', 'build_body: a constant-z face has no normal direction');
end
end

function tf = collapsed(r)
r = reshape(r, [], 3);
tf = all(all(r == r(1, :)));
end

function w = col(W, j)
if isempty(W)
    w = [];
else
    w = W(:, j);
end
end

function w = row(W, i)
if isempty(W)
    w = [];
else
    w = W(i, :)';
end
end

% ------------------------------------------------------------ curves

function s = same_curve(a, b)
% contract S1 same curve: 1 same direction, -1 reversed, 0 not the same curve
s = 0;
n = size(a.ctrl, 1);
if a.degree ~= b.degree || numel(a.knots) ~= numel(b.knots) || size(b.ctrl, 1) ~= n
    return
end
wa = a.weights(:);
wb = b.weights(:);
if isempty(wa), wa = ones(n, 1); end
if isempty(wb), wb = ones(n, 1); end
ka = (a.knots - a.knots(1)) / (a.knots(end) - a.knots(1));
kb = (b.knots - b.knots(1)) / (b.knots(end) - b.knots(1));
% four rounded operations at most on values in [0, 1] (contract S1)
tol = 4 * 2^-52;
if isequal(a.ctrl, b.ctrl) && isequal(wa, wb) && all(abs(ka - kb) <= tol)
    s = 1;
    return
end
kr = (b.knots(end) - b.knots(end:-1:1)) / (b.knots(end) - b.knots(1));
if isequal(a.ctrl, flipud(b.ctrl)) && isequal(wa, flipud(wb)) && all(abs(ka - kr) <= tol)
    s = -1;
end
end

function p = curve_midpoint(c)
C = mwecmass.solid.eval_bspline_curve(c, (c.knots(1) + c.knots(end)) / 2);
p = C(1, 1:2);
end

function A = curve_area(c)
% loop integral of x dy along the curve (Gauss-Legendre, 16 points per knot span)
[xg, wg] = gauss_legendre16();
ku = unique(c.knots);
A = 0;
for j = 1:numel(ku) - 1
    a = ku(j);
    b = ku(j + 1);
    s = (a + b) / 2 + (b - a) / 2 * xg;
    [C, Cs] = mwecmass.solid.eval_bspline_curve(c, s);
    A = A + (b - a) / 2 * (wg' * (C(:, 1) .* Cs(:, 2)));
end
end

function ang = curve_angle(c, pt)
% angle swept by the curve about pt: Bezier pieces subdivided until pt lies outside the box of
% their control points, where the curve and its chord sweep the same angle
ang = 0;
P = bezier_pieces(c);
for k = 1:numel(P)
    ang = ang + piece_angle(P{k}, pt, 0);
end
end

function a = piece_angle(H, pt, depth)
X = H(:, 1:2) ./ H(:, 3);
lo = min(X, [], 1);
hi = max(X, [], 1);
if any(pt < lo) || any(pt > hi)
    d0 = X(1, :) - pt;
    d1 = X(end, :) - pt;
    a = atan2(d0(1) * d1(2) - d0(2) * d1(1), d0 * d1');
    return
end
if depth >= 60
    error('mwecmass:solid:PointOnLoop', 'build_body: a test point lies on a section loop');
end
[H1, H2] = bezier_halves(H);
a = piece_angle(H1, pt, depth + 1) + piece_angle(H2, pt, depth + 1);
end

function P = bezier_pieces(c)
% homogeneous 2-D Bezier pieces [wx wy w] of a (rational) B-spline curve (knot insertion to
% multiplicity p at every interior knot)
p = c.degree;
n = size(c.ctrl, 1);
w = c.weights(:);
if isempty(w)
    w = ones(n, 1);
end
Q = [c.ctrl(:, 1:2) .* w, w];
kn = c.knots(:)';
for u = unique(kn(p + 2:end - p - 1))
    r = p - sum(kn == u);
    for it = 1:r
        k = find(kn <= u, 1, 'last');
        R = Q;
        Q = zeros(size(R, 1) + 1, 3);
        for i = 1:size(Q, 1)
            if i <= k - p
                Q(i, :) = R(i, :);
            elseif i > k
                Q(i, :) = R(i - 1, :);
            else
                al = (u - kn(i)) / (kn(i + p) - kn(i));
                Q(i, :) = al * R(i, :) + (1 - al) * R(i - 1, :);
            end
        end
        kn = [kn(1:k), u, kn(k + 1:end)];
    end
end
m = (size(Q, 1) - 1) / p;
P = cell(1, m);
for j = 1:m
    P{j} = Q((j - 1) * p + 1:j * p + 1, :);
end
end

function [A, B] = bezier_halves(H)
% de Casteljau at 1/2 on homogeneous control points
n = size(H, 1);
A = zeros(n, 3);
B = zeros(n, 3);
T = H;
A(1, :) = T(1, :);
B(n, :) = T(n, :);
for r = 1:n - 1
    T = (T(1:end - 1, :) + T(2:end, :)) / 2;
    A(r + 1, :) = T(1, :);
    B(n - r, :) = T(end, :);
end
end

function tf = pieces_disjoint(A, B)
% no point of the Bezier pieces A lies on a piece of B: separated by the boxes of their control
% points or by the fat line of one of them, after subdivision
tf = true;
for i = 1:numel(A)
    for j = 1:numel(B)
        if ~pair_disjoint(A{i}, B{j}, 0)
            tf = false;
            return
        end
    end
end
end

function tf = pair_disjoint(H, G, depth)
X = H(:, 1:2) ./ H(:, 3);
Y = G(:, 1:2) ./ G(:, 3);
if any(max(X, [], 1) < min(Y, [], 1)) || any(max(Y, [], 1) < min(X, [], 1))
    tf = true;
    return
end
if fat_separated(X, Y) || fat_separated(Y, X)
    tf = true;
    return
end
if depth >= 60
    tf = false;
    return
end
if max(max(X, [], 1) - min(X, [], 1)) >= max(max(Y, [], 1) - min(Y, [], 1))
    [H1, H2] = bezier_halves(H);
    tf = pair_disjoint(H1, G, depth + 1) && pair_disjoint(H2, G, depth + 1);
else
    [G1, G2] = bezier_halves(G);
    tf = pair_disjoint(H, G1, depth + 1) && pair_disjoint(H, G2, depth + 1);
end
end

function tf = fat_separated(X, Y)
% Y's control points all lie outside the band of X's control points about X's chord
d = X(end, :) - X(1, :);
L = norm(d);
tf = false;
if L == 0
    return
end
nrm = [-d(2) d(1)] / L;
dx = (X - X(1, :)) * nrm';
dy = (Y - X(1, :)) * nrm';
tf = min(dy) > max(dx) || max(dy) < min(dx);
end

function [x, w] = gauss_legendre16()
persistent xs ws
if isempty(xs)
    n = 16;
    b = (1:n - 1) ./ sqrt(4 * (1:n - 1).^2 - 1);
    [V, D] = eig(diag(b, 1) + diag(b, -1));
    [xs, i] = sort(diag(D));
    ws = 2 * V(1, i)'.^2;
end
x = xs;
w = ws;
end
