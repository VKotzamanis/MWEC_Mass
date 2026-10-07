function out = sti_closed_form(what, varargin)
%STI_CLOSED_FORM  Closed forms of the stand-in fixtures cylinder.ms2 and box.ms2.
%
%   fx  = sti_closed_form('fixture', src)          fixture description; src is a fixture name
%                                                  ('cylinder' | 'box'), an MS2Parser model of a
%                                                  fixture deck, or a geo (S1b) of one
%   P   = sti_closed_form('patches', fx, delta)    S1 patches of the hull shrunk inward by delta
%                                                  (delta = 0: the exact NURBS of the deck, contract F1)
%   dc  = sti_closed_form('d_close', fx)           offset distance at which the void closes (F2b)
%   sec = sti_closed_form('section', fx, delta)    moments of the horizontal section at offset delta
%   lay = sti_closed_form('layout', fx, design, d) air interval of every module (S3 layout)
%   reg = sti_closed_form('regions', fx, design, d)        V, S, J per region and module (S6 data)
%   rs  = sti_closed_form('region_sections', fx, design, d, i, z)  section moments per region at z
%   hs  = sti_closed_form('hydrostatics', fx, vs)  S7 numbers at vertical shift vs (no loop)
%   hl  = sti_closed_form('hull', fx)              V, S, J of the whole hull
%
%   d is the offset distance per module [N x 1] (d = t + eps_fit/2, contract section 0; NaN for a
%   module without void). Body frame, metres. A source that is not a fixture deck errors with
%   mwecmass:standin:NotAnalytic.
%
%   Sources of the closed forms (every section is constant in z, so every region is a union of
%   prisms):
%   - cylinder of radius R about the z axis between z0 and z1, shrunk by d: radius R - d between
%     z0 + d and z1 - d (the sharp rims are convex creases, so the offsets of the side and of the
%     end disks meet without a fold); box [x0,x1] x [y0,y1] x [z0,z1] shrunk by d: each face moved
%     inward by d (contract section 3, Stand-in kit SK).
%   - disk of radius r centred at (cx, cy): A = pi r^2, int x^2 dA = pi r^4/4 + A cx^2 (parallel
%     axis); rectangle: products of int_x0^x1 x^k dx = (x1^(k+1) - x0^(k+1))/(k+1).
%   - prism of section A between za and zb: int dV = A (zb - za), int z dV = A (zb^2 - za^2)/2,
%     int z^2 dV = A (zb^3 - za^3)/3, int x z dV = (int x dA)(zb^2 - za^2)/2.
%   - hydrostatics of a vertical prism: V_sub = A (zw - z0), CB at mid-depth, waterplane = section,
%     I_wp about the centre of flotation by the parallel-axis theorem, S_wet = A + perimeter * depth.

switch what
    case 'fixture'
        out = fixture(varargin{1});
    case 'patches'
        out = patches(varargin{:});
    case 'd_close'
        fx = varargin{1};
        if strcmp(fx.kind, 'cylinder')
            out = min(fx.R, diff(fx.z) / 2);
        else
            out = min([diff(fx.x), diff(fx.y), diff(fx.z)]) / 2;
        end
    case 'section'
        out = section(varargin{:});
    case 'layout'
        out = layout(varargin{:});
    case 'regions'
        out = regions(varargin{:});
    case 'region_sections'
        out = region_sections(varargin{:});
    case 'hydrostatics'
        out = hydrostatics(varargin{:});
    case 'hull'
        out = hull(varargin{1});
    otherwise
        error('sti_closed_form: unknown request %s', what);
end
end

% ---------------------------------------------------------------- fixture

function fx = fixture(src)
if isstruct(src)
    if isfield(src, 'analytic') && ~isempty(src.analytic)
        fx = src.analytic;
        return
    end
    if isfield(src, 'hull_name') && any(strcmp(src.hull_name, {'cylinder', 'box'}))
        fx = fixture(src.hull_name);
        return
    end
    not_analytic('the geo is not one of the stand-in fixture decks');
end
if ischar(src)
    if ~any(strcmp(src, {'cylinder', 'box'}))
        not_analytic(sprintf('%s is not a stand-in fixture name', src));
    end
    deck = fullfile(fileparts(mfilename('fullpath')), [src '.ms2']);
    evalc('model = mwecmass.geometry.MS2Parser.parse(deck);');
    fx = fixture(model);
    return
end
if ~isa(src, 'mwecmass.geometry.MS2Parser')
    not_analytic('expected a fixture name, an MS2Parser model or a geo');
end
model = src;
[~, stem] = fileparts(model.filename);
has = @(names) all(cellfun(@(n) model.entities.isKey(n), names));
topo = model.classify_visible_surfaces();
fx = struct('kind', '', 'name', stem, 'deck', model.filename, 'z', [], 'R', [], 'x', [], 'y', [], ...
    'surfs', struct('name', model.visible_surfs, 'flips', {{}}));
for k = 1:numel(fx.surfs)
    i = find(strcmp({topo.mirrors.name}, fx.surfs(k).name));
    if ~isempty(i)
        fx.surfs(k).flips = topo.mirrors(i).effective_flips;
    end
end
if strcmp(stem, 'cylinder') && has({'cylinder', 'profile', 'axis', 'K', 'P1', 'P2', 'T'})
    K = model.eval_point('K');
    P1 = model.eval_point('P1');
    P2 = model.eval_point('P2');
    T = model.eval_point('T');
    e = model.entities('cylinder');
    ok = isequal(K(1:2), [0 0]) && isequal(T(1:2), [0 0]) && P1(2) == 0 && P2(2) == 0 && ...
        P1(1) == P2(1) && P1(1) > 0 && P1(3) == K(3) && P2(3) == T(3) && T(3) > K(3) && ...
        e.params.angle_start == 0 && e.params.angle_end == 90 && numel(model.visible_surfs) == 4;
    if ~ok
        not_analytic('cylinder.ms2 does not have the fixture layout');
    end
    fx.kind = 'cylinder';
    fx.R = P1(1);
    fx.z = [K(3) T(3)];
elseif strcmp(stem, 'box') && has({'box_side1', 'box_side2', 'box_side3', 'box_side4', 'box_bottom', ...
        'box_top', 'B1', 'B3', 'T1', 'T3'})
    B1 = model.eval_point('B1');
    B3 = model.eval_point('B3');
    T1 = model.eval_point('T1');
    if numel(model.visible_surfs) ~= 6 || ~isequal(T1(1:2), B1(1:2))
        not_analytic('box.ms2 does not have the fixture layout');
    end
    fx.kind = 'box';
    fx.x = [B3(1) B1(1)];
    fx.y = [B1(2) B3(2)];
    fx.z = [B1(3) T1(3)];
else
    not_analytic(sprintf('%s is not a stand-in fixture deck', model.filename));
end
end

function not_analytic(msg)
error('mwecmass:standin:NotAnalytic', 'stand-in: %s (closed forms exist only for tests/standins/fixtures)', msg);
end

% ---------------------------------------------------------------- patches (S1)

function P = patches(fx, delta)
if nargin < 2
    delta = 0;
end
z0 = fx.z(1) + delta;
z1 = fx.z(2) - delta;
names = {fx.surfs.name};
P = repmat(empty_patch(), 1, numel(names));
for k = 1:numel(names)
    p = empty_patch();
    p.name = names{k};
    p.flips = fx.surfs(k).flips;
    p.exact = delta == 0;
    p.z_of_u = true;
    p.u_range = [0 1];
    if strcmp(fx.kind, 'cylinder')
        r = fx.R - delta;
        w = sqrt(2) / 2;
        ctrl = zeros(4, 3, 3);
        ctrl(1, :, :) = arc_row(0, z0);
        ctrl(2, :, :) = arc_row(r, z0);
        ctrl(3, :, :) = arc_row(r, z1);
        ctrl(4, :, :) = arc_row(0, z1);
        fl = flip_mask(p.flips);
        for c = 1:2
            if fl(c)
                ctrl(:, :, c) = -ctrl(:, :, c);
            end
        end
        p.source = 'cylinder';
        p.type = 'RevSurf';
        p.surf = struct('type', 'bspline', 'degree', [1 2], 'ctrl', ctrl, ...
            'knots', {{[0 0 1/3 2/3 1 1], [0 0 0 1 1 1]}}, 'weights', repmat([1 w 1], 4, 1));
        % S_u x S_v of the unmirrored quarter points to the axis; every mirror reverses it.
        p.outward = mod(sum(fl), 2) == 1;
        p.z_range = [z0 z1];
        p.offset_kind = 'rev_z';
        p.pole = [true true];
        p.c0_u = [1/3 2/3];
        p.c0_v = [];
        p.seam_v0 = [find_flip(fx, [fl(1) ~fl(2)]) 1];
        p.seam_v1 = [find_flip(fx, [~fl(1) fl(2)]) 3];
    else
        x0 = fx.x(1) + delta; x1 = fx.x(2) - delta;
        y0 = fx.y(1) + delta; y1 = fx.y(2) - delta;
        C = [x1 y0; x1 y1; x0 y1; x0 y0];
        side = @(j) find(strcmp(names, sprintf('box_side%d', j)));
        bottom = find(strcmp(names, 'box_bottom'));
        top = find(strcmp(names, 'box_top'));
        bnd = [2 3 4 1];   % boundary of box_bottom / box_top that side j meets
        p.source = p.name;
        p.type = 'RuledSurf';
        p.offset_kind = 'ruled_parallel';
        p.pole = [false false];
        p.c0_u = [];
        p.c0_v = [];
        ctrl = zeros(2, 2, 3);
        j = sscanf(p.name, 'box_side%d');
        if ~isempty(j)
            jn = mod(j, 4) + 1;
            ctrl(:, 1, :) = [C(j, :) z0; C(j, :) z1];
            ctrl(:, 2, :) = [C(jn, :) z0; C(jn, :) z1];
            p.outward = false;
            p.z_range = [z0 z1];
            p.seam_v0 = [side(mod(j - 2, 4) + 1) 3];
            p.seam_v1 = [side(jn) 1];
            p.seam_u0 = [bottom bnd(j)];
            p.seam_u1 = [top bnd(j)];
        else
            zc = z0;
            if strcmp(p.name, 'box_top')
                zc = z1;
            end
            ctrl(:, 1, :) = [x0 y0 zc; x1 y0 zc];
            ctrl(:, 2, :) = [x0 y1 zc; x1 y1 zc];
            p.outward = strcmp(p.name, 'box_top');
            p.z_range = [zc zc];
            b = 4 - strcmp(p.name, 'box_top') * 2;   % sides meet the bottom with u0 (4), the top with u1 (2)
            p.seam_v0 = [side(4) b];
            p.seam_u1 = [side(1) b];
            p.seam_v1 = [side(2) b];
            p.seam_u0 = [side(3) b];
        end
        p.surf = struct('type', 'bspline', 'degree', [1 1], 'ctrl', ctrl, ...
            'knots', {{[0 0 1 1], [0 0 1 1]}}, 'weights', []);
    end
    P(k) = p;
end

    function row = arc_row(rr, zz)
        row = reshape([rr 0 zz; rr rr zz; 0 rr zz], 1, 3, 3);
    end
end

function p = empty_patch()
p = struct('name', '', 'source', '', 'type', '', 'flips', {{}}, 'surf', [], 'outward', false, ...
    'exact', false, 'z_of_u', false, 'u_range', [], 'z_range', [], 'offset_kind', '', ...
    'pole', [false false], 'c0_u', [], 'c0_v', [], 'seam_u0', [], 'seam_u1', [], ...
    'seam_v0', [], 'seam_v1', []);
end

function fl = flip_mask(flips)
fl = [any(strcmp(flips, 'X')) any(strcmp(flips, 'Y'))];
end

function k = find_flip(fx, want)
k = find(arrayfun(@(s) isequal(flip_mask(s.flips), want), fx.surfs));
end

% ---------------------------------------------------------------- sections and regions

function s = section(fx, delta)
% A, Sx = int x dA, Sy = int y dA, Ixx = int y^2 dA, Iyy = int x^2 dA, Ixy = int x y dA, perimeter
if strcmp(fx.kind, 'cylinder')
    r = fx.R - delta;
    A = pi * r^2;
    s = struct('A', A, 'Sx', 0, 'Sy', 0, 'Ixx', pi * r^4 / 4, 'Iyy', pi * r^4 / 4, 'Ixy', 0, ...
        'perimeter', 2 * pi * r);
else
    x = fx.x + [delta -delta];
    y = fx.y + [delta -delta];
    px = @(k) (x(2)^(k + 1) - x(1)^(k + 1)) / (k + 1);
    py = @(k) (y(2)^(k + 1) - y(1)^(k + 1)) / (k + 1);
    s = struct('A', px(0) * py(0), 'Sx', px(1) * py(0), 'Sy', px(0) * py(1), ...
        'Ixx', px(0) * py(2), 'Iyy', px(2) * py(0), 'Ixy', px(1) * py(1), ...
        'perimeter', 2 * (diff(x) + diff(y)));
end
end

function lay = layout(fx, design, d)
% Air interval [a, b] of every module (empty: no air) and the kind of its bottom and top.
e = design.edges(:);
N = numel(e) - 1;
d = d(:);
lay = struct('a', cell(N, 1), 'b', [], 'air', false, 'bottom', '', 'top', '');
zb = design.z_ballast;
for i = 1:N
    lay(i).a = NaN;
    lay(i).b = NaN;
    if any(design.solid_modules == i) || ~isfinite(d(i))
        continue
    end
    z_lo = fx.z(1) + d(i);
    z_hi = fx.z(2) - d(i);
    a = max([e(i), zb, z_lo]);
    b = min(e(i + 1), z_hi);
    if a < b
        lay(i).a = a;
        lay(i).b = b;
        lay(i).air = true;
        if a == z_lo
            lay(i).bottom = 'inner';
        elseif a == zb
            lay(i).bottom = 'ballast_top';
        else
            lay(i).bottom = 'joint';
        end
        if b == z_hi
            lay(i).top = 'inner';
        else
            lay(i).top = 'joint';
        end
    end
end
end

function names = region_names(mode)
if strcmp(mode, 'modular_precast')
    names = {'uhpc', 'air'};
else
    names = {'ballast', 'shell', 'air'};
end
end

function reg = regions(fx, design, d)
e = design.edges(:);
N = numel(e) - 1;
names = region_names(design.mode);
lay = layout(fx, design, d);
for r = 1:numel(names)
    reg.(names{r}) = struct('V', zeros(N, 1), 'S', zeros(N, 3), 'J', zeros(3, 3, N));
end
so = section(fx, 0);
reg.V_module = so.A * (e(2:end) - e(1:end - 1));
for i = 1:N
    br = [e(i); e(i + 1)];
    for z = [design.z_ballast, lay(i).a, lay(i).b]
        if z > e(i) && z < e(i + 1)
            br(end + 1) = z; %#ok<AGROW>
        end
    end
    br = unique(br);
    for k = 1:numel(br) - 1
        za = br(k);
        zz = br(k + 1);
        rs = region_sections(fx, design, d, i, (za + zz) / 2, lay);
        for r = 1:numel(names)
            s = rs.(names{r});
            if s.A == 0
                continue
            end
            [V, S, J] = prism(s, za, zz);
            reg.(names{r}).V(i) = reg.(names{r}).V(i) + V;
            reg.(names{r}).S(i, :) = reg.(names{r}).S(i, :) + S;
            reg.(names{r}).J(:, :, i) = reg.(names{r}).J(:, :, i) + J;
        end
    end
end
end

function rs = region_sections(fx, design, d, i, z, lay)
% Section moments of every region of module i at height z (z inside the module).
if nargin < 6
    lay = layout(fx, design, d);
end
names = region_names(design.mode);
zero = struct('A', 0, 'Sx', 0, 'Sy', 0, 'Ixx', 0, 'Iyy', 0, 'Ixy', 0, 'perimeter', 0);
for r = 1:numel(names)
    rs.(names{r}) = zero;
end
so = section(fx, 0);
precast = strcmp(design.mode, 'modular_precast');
if precast
    solid = 'uhpc';
else
    solid = 'shell';
end
if z < design.z_ballast
    if precast
        rs.uhpc = so;
    else
        rs.ballast = so;
    end
elseif lay(i).air && z > lay(i).a && z < lay(i).b
    si = section(fx, d(i));
    rs.air = si;
    rs.(solid) = minus(so, si);
else
    rs.(solid) = so;
end
end

function s = minus(a, b)
f = {'A', 'Sx', 'Sy', 'Ixx', 'Iyy', 'Ixy'};
s = a;
for k = 1:numel(f)
    s.(f{k}) = a.(f{k}) - b.(f{k});
end
s.perimeter = a.perimeter + b.perimeter;
end

function [V, S, J] = prism(s, za, zb)
h = zb - za;
z1 = (zb^2 - za^2) / 2;
z2 = (zb^3 - za^3) / 3;
V = s.A * h;
S = [s.Sx * h, s.Sy * h, s.A * z1];
J = [s.Iyy * h, s.Ixy * h, s.Sx * z1; ...
     s.Ixy * h, s.Ixx * h, s.Sy * z1; ...
     s.Sx * z1, s.Sy * z1, s.A * z2];
end

function hl = hull(fx)
s = section(fx, 0);
[hl.V, hl.S, hl.J] = prism(s, fx.z(1), fx.z(2));
hl.centroid = hl.S / hl.V;
hl.area = 2 * s.A + s.perimeter * diff(fx.z);
end

% ---------------------------------------------------------------- hydrostatics (S7 numbers)

function hs = hydrostatics(fx, vs)
s = section(fx, 0);
z0 = fx.z(1);
z1 = fx.z(2);
zw = -vs;
hs = struct('vs', vs, 'draft', -(z0 + vs), 'submersion', '', 'V_sub', 0, 'CB_body', NaN(1, 3), ...
    'CB', NaN(1, 3), 'S_wet', 0, 'Aw', 0, 'xF', NaN, 'yF', NaN, 'I_wp_xx', 0, 'I_wp_yy', 0, ...
    'I_wp_yy_origin', 0, 'BM_L', NaN, 'KM', NaN);
if zw <= z0
    hs.submersion = 'none';
    return
end
if zw >= z1
    hl = hull(fx);
    hs.submersion = 'full';
    hs.V_sub = hl.V;
    hs.CB_body = hl.centroid;
    hs.S_wet = hl.area;
    hs.BM_L = 0;
else
    hs.submersion = 'partial';
    hs.V_sub = s.A * (zw - z0);
    hs.CB_body = [s.Sx / s.A, s.Sy / s.A, (z0 + zw) / 2];
    hs.S_wet = s.A + s.perimeter * (zw - z0);
    hs.Aw = s.A;
    hs.xF = s.Sx / s.A;
    hs.yF = s.Sy / s.A;
    hs.I_wp_xx = s.Ixx - s.A * hs.yF^2;
    hs.I_wp_yy = s.Iyy - s.A * hs.xF^2;
    hs.I_wp_yy_origin = s.Iyy;
    hs.BM_L = hs.I_wp_yy / hs.V_sub;
end
hs.CB = hs.CB_body + [0 0 vs];
hs.KM = hs.CB(3) + hs.BM_L;
end
