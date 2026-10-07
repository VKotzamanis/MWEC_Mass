function body = build_body(geo, design, inner)
%BUILD_BODY  Stand-in of contract F5: the S4 body of an S3 design on a stand-in fixture.
%
%   body = mwecmass.solid.build_body(geo, design, inner)
%
%   geo: S1b with geo.analytic set (stand-in F1), else mwecmass:standin:NotAnalytic. inner: S2 sets
%   (the F2 stand-in, which returns sti_inner_box for the box, or sti_inner_box directly), one
%   per distinct t of design.t.
%   Layout of contract S3 (module air intervals from sti_closed_form 'layout'); faces per mode as
%   in S4. Every face is built once, edges and vertices are shared by index. Lateral faces are the
%   u-degree-1 patches cut at their knot rows (c0_u: every interior knot of degree 1) and at the
%   planes by knot insertion with alpha from the closed form z(u) of split_bspline_surface; a
%   collapsed row (pole) is omitted from the loop. A row cut at a plane takes that height bitwise
%   (knot insertion moves z by rounding only), so the faces' control rows, their loops and the
%   F6b comparisons with z agree exactly. Plane faces (caps, joint_step, ballast_top)
%   take their loops from the rows of the lateral faces at that height. At z_ballast equal to the
%   inner z_lo of a module the layout gives its void bottom as ballast_top, so no inner end disk is
%   written there; ballast_top is then the void-bottom disk (precast: uhpc/air) or, in thin shell,
%   an annulus (ballast/shell) and a disk (ballast/air) (contract S3, S4). Joint loops of a fixture
%   are concentric, so the joint_step nesting of F5 holds by construction (JointNotNested cannot
%   occur here).
%
%   Extra face fields: role, inside, outside (region against / along the face normal, 'exterior'
%   outside the hull) and module: the module of the face, or [module of inside, module of
%   outside] for a face on a module joint. body.analytic (stand-in only, read by the F6 and F6b
%   stand-ins): fixture description, offset distance d per module, eps_fit.

if ~isstruct(geo) || ~isfield(geo, 'analytic') || isempty(geo.analytic)
    error('mwecmass:standin:NotAnalytic', 'build_body stand-in: geo.analytic is empty (not a stand-in fixture)');
end
if exist('sti_closed_form', 'file') ~= 2
    addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'fixtures'), '-end');
end
fx = geo.analytic;
e = design.edges(:);
N = numel(e) - 1;
zmin = geo.z_range(1);
zmax = geo.z_range(2);
if e(1) ~= zmin || e(end) ~= zmax || any(diff(e) <= 0)
    error('mwecmass:solid:BadEdges', 'build_body: module edges must rise from z_min to z_max');
end
precast = strcmp(design.mode, 'modular_precast');
t = design.t(:);
d = NaN(N, 1);
set_of = zeros(N, 1);
for i = 1:N
    if any(design.solid_modules == i) || ~isfinite(t(i))
        continue
    end
    k = find(arrayfun(@(s) isequal(s.t, t(i)), inner), 1);
    if isempty(k)
        error('mwecmass:solid:MissingInnerSet', 'build_body: no inner set for t = %g (module %d)', t(i), i);
    end
    d(i) = inner(k).d;
    set_of(i) = k;
end
zb = design.z_ballast;
lay = sti_closed_form('layout', fx, design, d);
planes = e(2:end - 1);
if zb > zmin && zb < zmax
    planes = [planes; zb];
end
planes = unique(planes);

brep = struct('vertices', zeros(0, 3), ...
    'curves', struct('degree', {}, 'ctrl', {}, 'knots', {}, 'weights', {}), ...
    'edges', struct('vertices', {}, 'curve', {}), 'surfaces', {{}}, ...
    'faces', struct('surface', {}, 'same_sense', {}, 'loops', {}, 'role', {}, 'module', {}, ...
    'inside', {}, 'outside', {}), ...
    'bodies', struct('name', {}, 'kind', {}, 'shells', {}));

% outer lateral faces
for p = 1:numel(geo.outer)
    P = geo.outer(p);
    s = P.surf;
    tk = unique(s.knots{1});
    for j = 1:numel(tk) - 1
        Za = row_z(s, tk(j));
        Zb = row_z(s, tk(j + 1));
        if Za == Zb
            if Za == zmin
                m = 1;
            else
                m = N;
            end
            add_patch_face(piece(s, tk(j), tk(j + 1)), 'outer', m, P.outward, material(Za, Za == zmin), 'exterior');
            continue
        end
        cuts = planes(planes > min(Za, Zb) & planes < max(Za, Zb));
        us = sort([tk(j); tk(j + 1); arrayfun(@(z) u_of_z(s, z), cuts)]);
        if Za < Zb
            hz = [Za; sort(cuts(:)); Zb];
        else
            hz = [Za; sort(cuts(:), 'descend'); Zb];
        end
        for q = 1:numel(us) - 1
            zm = (hz(q) + hz(q + 1)) / 2;
            add_patch_face(piece(s, us(q), us(q + 1), hz(q), hz(q + 1)), 'outer', module_of(zm), P.outward, material(zm, false), 'exterior');
        end
    end
end

% inner lateral faces, module by module
if precast
    solid_name = 'uhpc';
else
    solid_name = 'shell';
end
for i = 1:N
    if ~lay(i).air
        continue
    end
    S = inner(set_of(i));
    for q = 1:numel(S.patches)
        Q = S.patches(q);
        s = Q.surf;
        tk = unique(s.knots{1});
        Za = row_z(s, tk(1));
        Zb = row_z(s, tk(end));
        if Za == Zb
            if (Za == lay(i).a && strcmp(lay(i).bottom, 'inner')) || (Za == lay(i).b && strcmp(lay(i).top, 'inner'))
                add_patch_face(piece(s, tk(1), tk(end)), 'inner', i, Q.outward, solid_name, 'air');
            end
            continue
        end
        lo = max(lay(i).a, min(Za, Zb));
        hi = min(lay(i).b, max(Za, Zb));
        if lo >= hi
            continue
        end
        us = sort([u_of_z(s, lo), u_of_z(s, hi)]);
        hz = [lo hi];
        if Za > Zb
            hz = [hi lo];
        end
        add_patch_face(piece(s, us(1), us(2), hz(1), hz(2)), 'inner', i, Q.outward, solid_name, 'air');
    end
end

% plane faces
if precast
    for j = 2:N
        z = e(j);
        i = j - 1;
        below = lay(i).air && lay(i).b == z;
        above = lay(j).air && lay(j).a == z;
        mods = [i j];
        if ~below && ~above
            plane_face(z, 0, 0, 'cap', mods, 'uhpc', 'uhpc');
        elseif below && ~above
            plane_face(z, 0, set_of(i), 'cap', mods, 'uhpc', 'uhpc');
            if strcmp(lay(i).top, 'joint')
                plane_face(z, set_of(i), 0, 'cap', mods, 'air', 'uhpc');
            end
        elseif above && ~below
            plane_face(z, 0, set_of(j), 'cap', mods, 'uhpc', 'uhpc');
            if strcmp(lay(j).bottom, 'ballast_top')
                plane_face(z, set_of(j), 0, 'ballast_top', mods, 'uhpc', 'air');
            elseif strcmp(lay(j).bottom, 'joint')
                plane_face(z, set_of(j), 0, 'cap', mods, 'uhpc', 'air');
            end
        else
            if d(i) <= d(j)
                ks = set_of(i); kl = set_of(j);
            else
                ks = set_of(j); kl = set_of(i);
            end
            plane_face(z, 0, ks, 'cap', mods, 'uhpc', 'uhpc');
            if d(i) > d(j)
                plane_face(z, ks, kl, 'joint_step', mods, 'uhpc', 'air');
            elseif d(i) < d(j)
                plane_face(z, ks, kl, 'joint_step', mods, 'air', 'uhpc');
            end
        end
    end
    for i = 1:N
        if lay(i).air && strcmp(lay(i).bottom, 'ballast_top') && lay(i).a > e(i)
            plane_face(zb, set_of(i), 0, 'ballast_top', i, 'uhpc', 'air');
        end
    end
elseif zb > zmin && zb < zmax
    m_in = module_of(zb);
    m_out = m_in;
    jj = find(e == zb);
    if ~isempty(jj)
        m_in = jj - 1;
        m_out = jj;
    end
    mods = unique([m_in m_out], 'stable');
    cut = find(arrayfun(@(L) L.air && L.a == zb && strcmp(L.bottom, 'ballast_top'), lay), 1);
    if isempty(cut)
        plane_face(zb, 0, 0, 'ballast_top', mods, 'ballast', 'shell');
    else
        plane_face(zb, 0, set_of(cut), 'ballast_top', mods, 'ballast', 'shell');
        plane_face(zb, set_of(cut), 0, 'ballast_top', mods, 'ballast', 'air');
    end
end

% bodies and shells
hull = geo.hull_name;
if precast
    for i = 1:N
        brep.bodies(end + 1) = struct('name', sprintf('%s_UHPC_module_%d', hull, i), 'kind', 'solid', ...
            'shells', {solid_shells(bound('uhpc', i))});
    end
    voids = components(bound('air', []));
    shells = struct('all_outer', -bound('exterior', []), 'all_voids', {voids});
else
    if any(bound('ballast', []))
        brep.bodies(end + 1) = struct('name', sprintf('%s_STEEL_ballast', hull), 'kind', 'solid', ...
            'shells', {solid_shells(bound('ballast', []))});
    end
    out = -bound('exterior', []);
    sheet = out(arrayfun(@(f) strcmp(brep.faces(abs(f)).inside, 'shell') || strcmp(brep.faces(abs(f)).outside, 'shell'), out));
    if ~isempty(sheet)
        brep.bodies(end + 1) = struct('name', sprintf('%s_STEEL_shell', hull), 'kind', 'sheet', 'shells', {{sheet}});
    end
    shells = struct('layer', bound('shell', []), 'void', bound('air', []));
end
used = set_of(set_of > 0);
eps_fit = [];
inner_t = [];
if ~isempty(used)
    eps_fit = inner(used(1)).eps_fit;
    inner_t = unique([inner(unique(used)).t]);
end
body = struct('design', design, 'inner_t', inner_t, 'planes', planes(:)', ...
    'brep', brep, 'shells', shells, 'analytic', struct('fixture', fx, 'd', d, 'eps_fit', eps_fit));

% ------------------------------------------------------------ nested helpers

    function m = module_of(z)
        m = find(e(1:end - 1) <= z & z < e(2:end), 1);
        if isempty(m)
            m = N;
        end
    end

    function r = material(z, at_bottom)
        if precast
            r = 'uhpc';
        elseif z < zb || (at_bottom && zb > zmin) || zb >= zmax
            r = 'ballast';
        else
            r = 'shell';
        end
    end

    function add_patch_face(s, role, m, outward, solid, other)
        loop = [];
        [nu, nv, ~] = size(s.ctrl);
        W = s.weights;
        bnd = {{squeeze(s.ctrl(:, 1, :)), s.knots{1}, col(W, 'c', 1), 1}, ...
               {squeeze(s.ctrl(nu, :, :)), s.knots{2}, col(W, 'r', nu), 1}, ...
               {squeeze(s.ctrl(:, nv, :)), s.knots{1}, col(W, 'c', nv), -1}, ...
               {squeeze(s.ctrl(1, :, :)), s.knots{2}, col(W, 'r', 1), -1}};
        degs = [s.degree(1), s.degree(2), s.degree(1), s.degree(2)];
        for b = 1:4
            C = bnd{b}{1};
            if all(all(C == C(1, :)))
                continue
            end
            crv = struct('degree', degs(b), 'ctrl', C, 'knots', bnd{b}{2}, 'weights', bnd{b}{3});
            loop(end + 1) = bnd{b}{4} * edge_for(crv); %#ok<AGROW>
        end
        if outward
            inside = solid; outside = other;
        else
            inside = other; outside = solid;
        end
        brep.surfaces{end + 1} = s;
        brep.faces(end + 1) = struct('surface', numel(brep.surfaces), 'same_sense', true, ...
            'loops', {{loop}}, 'role', role, 'module', m, 'inside', inside, 'outside', outside);
    end

    function plane_face(z, outer_set, hole_set, role, m, inside, outside)
        loops = {section_edges(z, outer_set)};
        if hole_set > 0
            h = section_edges(z, hole_set);
            loops{2} = -h(end:-1:1);
        end
        brep.surfaces{end + 1} = struct('type', 'plane', 'origin', [0 0 z], 'normal', [0 0 1]);
        brep.faces(end + 1) = struct('surface', numel(brep.surfaces), 'same_sense', true, ...
            'loops', {loops}, 'role', role, 'module', m, 'inside', inside, 'outside', outside);
    end

    function L = section_edges(z, set)
        % signed edges of the section of the outer patches (set 0) or inner set at z, CCW from +z
        if set == 0
            pats = geo.outer;
        else
            pats = inner(set).patches;
        end
        segs = [];
        for kk = 1:numel(pats)
            sf = pats(kk).surf;
            Z = sf.ctrl(:, 1, 3);
            if all(Z == Z(1)) || z < min(Z) || z > max(Z)
                continue
            end
            r = piece(sf, u_of_z(sf, z), u_of_z(sf, z), z, z);
            crv = struct('degree', sf.degree(2), 'ctrl', squeeze(r.ctrl(1, :, :)), 'knots', sf.knots{2}, ...
                'weights', col(r.weights, 'r', 1));
            segs(end + 1) = edge_for(crv); %#ok<AGROW>
        end
        L = chain(segs);
    end

    function L = chain(segs)
        L = segs(1);
        [first, cur] = edge_ends(L);
        rest = segs(2:end);
        while ~isempty(rest)
            hit = 0;
            for kk = 1:numel(rest)
                [a, b] = edge_ends(rest(kk));
                if a == cur
                    hit = rest(kk);
                elseif b == cur
                    hit = -rest(kk);
                end
                if hit
                    rest(kk) = [];
                    break
                end
            end
            if ~hit
                error('mwecmass:solid:SectionNotClosed', 'build_body stand-in: section loop does not close');
            end
            L(end + 1) = hit; %#ok<AGROW>
            [~, cur] = edge_ends(hit);
        end
        if cur ~= first
            error('mwecmass:solid:SectionNotClosed', 'build_body stand-in: section loop does not close');
        end
        V = zeros(numel(L), 2);
        for kk = 1:numel(L)
            [a, ~] = edge_ends(L(kk));
            V(kk, :) = brep.vertices(a, 1:2);
        end
        area = sum(V(:, 1) .* V([2:end 1], 2) - V([2:end 1], 1) .* V(:, 2));
        if area < 0
            L = -L(end:-1:1);
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

    function idx = edge_for(crv)
        for kk = 1:numel(brep.edges)
            c = brep.curves(brep.edges(kk).curve);
            if c.degree ~= crv.degree
                continue
            end
            if isequal(c.ctrl, crv.ctrl) && isequal(c.knots, crv.knots) && isequal(c.weights, crv.weights)
                idx = kk;
                return
            end
            kr = crv.knots(1) + crv.knots(end) - crv.knots(end:-1:1);
            if isequal(c.ctrl, flipud(crv.ctrl)) && isequal(c.knots, kr) && isequal(c.weights, flipud(crv.weights))
                idx = -kk;
                return
            end
        end
        brep.curves(end + 1) = crv;
        brep.edges(end + 1) = struct('vertices', [vertex(crv.ctrl(1, :)), vertex(crv.ctrl(end, :))], ...
            'curve', numel(brep.curves));
        idx = numel(brep.edges);
    end

    function i = vertex(xyz)
        i = find(all(brep.vertices == xyz, 2), 1);
        if isempty(i)
            brep.vertices(end + 1, :) = xyz;
            i = size(brep.vertices, 1);
        end
    end

    function L = bound(region, m)
        L = zeros(1, 0);
        for kk = 1:numel(brep.faces)
            f = brep.faces(kk);
            if strcmp(f.inside, region) && (isempty(m) || f.module(1) == m)
                L(end + 1) = kk; %#ok<AGROW>
            elseif strcmp(f.outside, region) && (isempty(m) || f.module(end) == m)
                L(end + 1) = -kk; %#ok<AGROW>
            end
        end
    end

    function C = components(L)
        % face lists connected through shared edges
        n = numel(L);
        lab = 1:n;
        E = cell(1, n);
        for kk = 1:n
            lp = brep.faces(abs(L(kk))).loops;
            E{kk} = [];
            for qq = 1:numel(lp)
                E{kk} = [E{kk}, abs(lp{qq}(:)')];
            end
        end
        changed = true;
        while changed
            changed = false;
            for a = 1:n
                for b = a + 1:n
                    if lab(a) ~= lab(b) && ~isempty(intersect(E{a}, E{b}))
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
        % outer shell first (it holds the exterior faces); further components are cavities,
        % oriented as standalone solids (normals out of the cavity)
        C = components(L);
        outer = false(1, numel(C));
        for c = 1:numel(C)
            for f = abs(C{c})
                outer(c) = outer(c) || strcmp(brep.faces(f).inside, 'exterior') || strcmp(brep.faces(f).outside, 'exterior');
            end
        end
        sh = C(outer);
        for c = find(~outer)
            sh{end + 1} = -C{c}; %#ok<AGROW>
        end
    end
end

% ------------------------------------------------------------ patch arithmetic (u-degree 1)

function w = col(W, kind, i)
if isempty(W)
    w = [];
elseif kind == 'c'
    w = W(:, i);
else
    w = W(i, :)';
end
end

function z = row_z(s, u)
r = piece(s, u, u);
z = r.ctrl(1, 1, 3);
end

function u = u_of_z(s, z)
% parameter of the iso-u row at height z (closed form of split_bspline_surface)
t = s.knots{1}(:)';
Z = s.ctrl(:, 1, 3);
W = s.weights;
if isempty(W)
    W = ones(size(s.ctrl, 1), size(s.ctrl, 2));
end
hit = find(Z == z & [true; Z(2:end) ~= Z(1:end - 1)] & [Z(1:end - 1) ~= Z(2:end); true], 1);
if ~isempty(hit)
    u = t(hit + 1);
    return
end
i = find((Z(1:end - 1) < z & Z(2:end) > z) | (Z(1:end - 1) > z & Z(2:end) < z), 1);
a0 = W(i, 1);
a1 = W(i + 1, 1);
s_ = a0 * (z - Z(i)) / (a0 * (z - Z(i)) + a1 * (Z(i + 1) - z));
u = t(i + 1) + s_ * (t(i + 2) - t(i + 1));
end

function r = piece(s, ua, ub, za, zb)
% the patch between u = ua and u = ub inside one knot span, by knot insertion (degree 1 in u);
% with za, zb the height of each row is set to the plane height it was cut at
[ra, wa] = row_at(s, ua);
[rb, wb] = row_at(s, ub);
if nargin > 3
    ra(:, :, 3) = za;
    rb(:, :, 3) = zb;
end
r = s;
r.ctrl = cat(1, ra, rb);
r.knots{1} = [ua ua ub ub];
if isempty(s.weights)
    r.weights = [];
else
    r.weights = [wa; wb];
end
end

function [row, w] = row_at(s, u)
t = s.knots{1}(:)';
W = s.weights;
if isempty(W)
    W = ones(size(s.ctrl, 1), size(s.ctrl, 2));
end
n = size(s.ctrl, 1);
j = find(t(2:n + 1) == u, 1);
if ~isempty(j)
    row = s.ctrl(j, :, :);
    w = W(j, :);
    return
end
i = find(t(2:n) < u & t(3:n + 1) > u, 1);
a = (u - t(i + 1)) / (t(i + 2) - t(i + 1));
if isequal(W(i, :), W(i + 1, :))
    row = (1 - a) * s.ctrl(i, :, :) + a * s.ctrl(i + 1, :, :);
    w = W(i, :);
else
    w = (1 - a) * W(i, :) + a * W(i + 1, :);
    row = ((1 - a) * s.ctrl(i, :, :) .* W(i, :) + a * s.ctrl(i + 1, :, :) .* W(i + 1, :)) ./ w;
end
end
