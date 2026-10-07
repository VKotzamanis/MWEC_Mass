function test_sk_body()
%TEST_SK_BODY  Stand-in F5, F6, F6b on both fixtures and both modes: B-rep validity, shared edges,
%   STEP import (closed solids, solid count), volume closure, sections.
%   The closed forms of sti_closed_form are the independent oracle of the OCC volumes (printed).
%   Reads body.analytic (stand-in marker of F5): J1, which merges the real F5, F6, F6b, deletes or
%   rewrites this test.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
evalc('mc = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''tests'', ''standins'', ''fixtures'', ''cylinder.ms2''));');
evalc('mb = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''tests'', ''standins'', ''fixtures'', ''box.ms2''));');
gc = mwecmass.solid.outer_nurbs(mc);
gb = mwecmass.solid.outer_nurbs(mb);
t_min = 0.1;
opts = struct('t_min', t_min);
ec = [-3; -2; -1; 0; 1];
eb = [-2.5; -1.5; -0.5; 0.5];
ic = [mwecmass.solid.offset_surface(mc, [], gc, 0.1, [-3 1], opts), ...
      mwecmass.solid.offset_surface(mc, [], gc, 0.2, [-3 1], opts)];
ib = [sti_inner_box(gb, 0.1, t_min), sti_inner_box(gb, 0.15, t_min)];
d1 = ic(1).d;
d2 = ic(2).d;

% {geo, inner sets, mode, edges, t, z_ballast, solid_modules, label}
cases = {
    gc, ic, 'modular_precast', ec, [0.1; 0.1; 0.2; NaN], -2.5, 4, 'cylinder precast, ballast inside module 1, joint_step, wall on top'
    gc, ic, 'modular_precast', ec, [0.1; 0.2; 0.2; NaN], -2, 4, 'cylinder precast, ballast at an edge'
    gc, ic, 'modular_precast', ec, [0.1; 0.1; 0.1; NaN], -1.5, 4, 'cylinder precast, ballast spilled into module 2'
    gc, ic, 'modular_precast', ec, [0.2; 0.1; 0.1; 0.1], -2.95, [], 'cylinder precast, ballast below the inner z_lo, closed void'
    gc, ic, 'thin_shell', ec, 0.1 * ones(4, 1), -2.5, [], 'cylinder thin shell, ballast inside module 1'
    gc, ic, 'thin_shell', ec, 0.1 * ones(4, 1), -1, [], 'cylinder thin shell, ballast at an edge'
    gc, ic, 'thin_shell', ec, 0.1 * ones(4, 1), -2.95, [], 'cylinder thin shell, ballast below the inner z_lo'
    gc, ic, 'thin_shell', ec, 0.2 * ones(4, 1), -3, [], 'cylinder thin shell, no ballast'
    gb, ib, 'modular_precast', eb, [0.1; 0.15; NaN], -2, 3, 'box precast, ballast inside module 1'
    gb, ib, 'modular_precast', eb, [0.15; 0.1; 0.1], -1.5, [], 'box precast, ballast at an edge'
    gb, ib, 'thin_shell', eb, 0.1 * ones(3, 1), -1, [], 'box thin shell, ballast in module 2'
    gb, ib, 'thin_shell', eb, 0.15 * ones(3, 1), -2.45, [], 'box thin shell, ballast below the inner z_lo'
    gb, [], 'modular_precast', eb, NaN(3, 1), -2.5, 1:3, 'box precast, every module solid (no inner set)'
    gc, [], 'modular_precast', ec, NaN(4, 1), -2, 1:4, 'cylinder precast, every module solid (no inner set)'
    gc, ic, 'modular_precast', [-3; -3 + d1; -1; 0; 1], [0.1; 0.1; 0.1; NaN], -3, 4, 'cylinder precast, inner z_lo on a module edge'
    gc, ic, 'modular_precast', [-3; -2; 1 - d1; 1], [0.1; 0.1; 0.1], -2.5, [], 'cylinder precast, inner z_hi on a module edge'
    gc, ic, 'modular_precast', [-3; -2; 1 - d2; 1], [0.1; 0.2; 0.1], -2.5, [], 'cylinder precast, thicker lower z_hi on a joint'
    gc, ic, 'modular_precast', ec, [0.1; 0.1; 0.2; NaN], ic(1).z_lo, 4, 'cylinder precast, z_ballast = inner z_lo of module 1'
    gc, ic, 'thin_shell', ec, 0.1 * ones(4, 1), ic(1).z_lo, [], 'cylinder thin shell, z_ballast = inner z_lo of module 1'
    gb, ib, 'modular_precast', eb, [0.1; 0.15; NaN], ib(1).z_lo, 3, 'box precast, z_ballast = inner z_lo of module 1'
    gb, ib, 'thin_shell', eb, 0.1 * ones(3, 1), ib(1).z_lo, [], 'box thin shell, z_ballast = inner z_lo of module 1'
    };
n_at_zlo = 0;
worst = struct('closure', 0, 'hull', 0, 'occ', 0, 'section', 0);
for c = 1:size(cases, 1)
    [geo, sets, mode, edges, t, zb, solid, label] = cases{c, :};
    design = struct('mode', mode, 'edges', edges, 'vs', 0.5, 't', t, 'z_ballast', zb, 'solid_modules', solid);
    body = mwecmass.solid.build_body(geo, design, sets);
    check(isequal(sort(fieldnames(body))', sort({'design', 'inner_t', 'planes', 'brep', 'shells', 'analytic'})), 'S4 fields');
    brep = body.brep;
    mwecmass.output.step.validate_brep(brep);
    n_shared = shared_edges(brep);
    fx = geo.analytic;
    N = numel(edges) - 1;
    precast = strcmp(mode, 'modular_precast');
    if precast
        rho = struct('uhpc', 2500, 'air', 1.2);
        names = {'uhpc', 'air'};
    else
        rho = struct('ballast', 7500, 'shell', 7850, 'air', 1.2);
        names = {'ballast', 'shell', 'air'};
    end
    bp = mwecmass.solid.body_properties(body, rho, []);
    reg = sti_closed_form('regions', fx, design, body.analytic.d);
    so = sti_closed_form('section', fx, 0);
    for i = 1:N
        Vsum = 0;
        for r = 1:numel(names)
            check(bp.modules(i).(['V_' names{r}]) >= 0, 'negative region volume');
            Vsum = Vsum + bp.modules(i).(['V_' names{r}]);
        end
        % each region is a sum of at most four prisms A * dz whose breakpoints enter through one
        % subtraction each; 16 eps on A * (sum of |breakpoints|) bounds their rounding
        bound = 16 * eps * so.A * (abs(edges(i)) + abs(edges(i + 1)) + sum(abs(fx.z)));
        err = abs(Vsum - reg.V_module(i));
        check(err <= bound, '%s: module %d volume closure %.3e > %.3e', label, i, err, bound);
        worst.closure = max(worst.closure, err / reg.V_module(i));
    end
    if isfinite(body.analytic.d(1)) && zb == fx.z(1) + body.analytic.d(1)
        check_ballast_at_z_lo(body, label, reg, bp, so);
        n_at_zlo = n_at_zlo + 1;
    end
    hl = sti_closed_form('hull', fx);
    err = abs(sum(reg.V_module) - hl.V);
    check(err <= 16 * eps * so.A * sum(abs(edges)), '%s: sum of module volumes differs from the hull', label);
    worst.hull = max(worst.hull, err / hl.V);

    % STEP import of the bodies
    file = [tempname() '.step'];
    cleanup = onCleanup(@() delete_if(file));
    mwecmass.output.step.write_step(brep, file);
    res = stp_check(file);
    check(strcmp(res.declared_length_unit, 'METRE') && res.open_edges_volumes == 0, '%s: STEP not closed', label);
    occ = res.occ_volumes(:)';
    if precast
        check(res.n_volumes == N, '%s: %d solids imported, expected %d', label, res.n_volumes, N);
        ref = arrayfun(@(m) m.V_uhpc, bp.modules);
        all_brep = brep;
        all_brep.bodies = struct('name', [geo.hull_name '_UHPC_all'], 'kind', 'solid', ...
            'shells', {[{body.shells.all_outer}, body.shells.all_voids]});
        mwecmass.output.step.validate_brep(all_brep);
        mwecmass.output.step.write_step(all_brep, file);
        r2 = stp_check(file);
        check(r2.n_volumes == 1 && r2.open_edges_volumes == 0, '%s: fused solid not one closed volume', label);
        occ = [occ, r2.occ_volumes]; %#ok<AGROW>
        ref = [ref, sum(ref)]; %#ok<AGROW>
    else
        has_ballast = zb > fx.z(1);
        check(res.n_volumes == has_ballast, '%s: %d solids imported', label, res.n_volumes);
        check(any(strcmp({brep.bodies.name}, [geo.hull_name '_STEEL_shell'])) && ...
            strcmp(brep.bodies(end).kind, 'sheet'), '%s: shell sheet missing', label);
        ref = sum(bp.regions.ballast.V) * ones(1, has_ballast);
        mass_brep = brep;
        mass_brep.bodies = struct('name', {'layer', 'void'}, 'kind', 'solid', ...
            'shells', {{body.shells.layer}, {body.shells.void}});
        mwecmass.output.step.validate_brep(mass_brep);
        mwecmass.output.step.write_step(mass_brep, file);
        r2 = stp_check(file);
        check(r2.n_volumes == 2 && r2.open_edges_volumes == 0, '%s: layer and void not closed', label);
        occ = [occ, r2.occ_volumes(:)']; %#ok<AGROW>
        ref = [ref, sum(bp.regions.shell.V), sum(bp.regions.air.V)]; %#ok<AGROW>
    end
    worst.occ = max([worst.occ, abs(occ - ref) ./ ref]);
    fprintf('%-70s faces %3d, shared edges %3d, OCC - closed form (relative): %s\n', label, ...
        numel(brep.faces), n_shared, sprintf('%.1e ', (occ - ref) ./ ref));

    % sections (F6b) inside every module, away from its planes
    for i = 1:N
        z = edges(i) + 0.6 * (edges(i + 1) - edges(i));
        sec = mwecmass.solid.body_section(body, z);
        check(sec.module == i && abs(sec.outer.area - so.A) <= 64 * eps * so.A, '%s: outer section', label);
        rs = sti_closed_form('region_sections', fx, design, body.analytic.d, i, z);
        check(sec.solid == (rs.air.A == 0), '%s: module %d solid flag', label, i);
        if ~sec.solid
            err = abs(sec.inner.area - rs.air.A) / rs.air.A;
            check(err <= 64 * eps, '%s: void section area', label);
            worst.section = max(worst.section, err);
        end
    end
end
check(n_at_zlo == 4, 'expected four cases with z_ballast = inner z_lo, ran %d', n_at_zlo);
fprintf('largest relative volume closure error per module %.3e, sum of modules vs hull %.3e\n', worst.closure, worst.hull);
fprintf('largest |OCC getMass - closed form| / closed form %.3e (printed, not asserted)\n', worst.occ);
fprintf('largest relative difference of F6b void areas from the closed form %.3e\n', worst.section);

body = mwecmass.solid.build_body(gc, struct('mode', 'thin_shell', 'edges', ec, 'vs', 0.5, ...
    't', 0.1 * ones(4, 1), 'z_ballast', -2.5, 'solid_modules', []), ic);
sec = mwecmass.solid.body_section(body, -2.7);
check(sec.solid && sec.module == 1, 'section below the ballast top is solid');

% F6b side: 'above' is the default and the half-open rule, 'below' the faces and module below
% a module edge or z_ballast; heights: z_min, z_ballast, module edges, z_max
designs = {
    gc, ic, struct('mode', 'modular_precast', 'edges', ec, 'vs', 0.5, 't', [0.1; 0.2; 0.1; NaN], 'z_ballast', -2.5, 'solid_modules', 4)
    gb, ib, struct('mode', 'modular_precast', 'edges', eb, 'vs', 0.5, 't', [0.1; 0.15; NaN], 'z_ballast', -2, 'solid_modules', 3)
    gc, ic, struct('mode', 'thin_shell', 'edges', ec, 'vs', 0.5, 't', 0.1 * ones(4, 1), 'z_ballast', -1, 'solid_modules', [])
    gb, ib, struct('mode', 'thin_shell', 'edges', eb, 'vs', 0.5, 't', 0.15 * ones(3, 1), 'z_ballast', -1, 'solid_modules', [])
    };
worst_side = 0;
n_side = 0;
for c = 1:size(designs, 1)
    body = mwecmass.solid.build_body(designs{c, 1}, designs{c, 3}, designs{c, 2});
    [w, n] = check_sides(body);
    worst_side = max(worst_side, w);
    n_side = n_side + n;
end
fprintf('F6b side: %d side/height checks, largest |outer above - outer below| %.3e (relative to the section area)\n', ...
    n_side, worst_side);
expect_error(@() mwecmass.solid.body_section(body, -0.5, 'left'), 'mwecmass:solid:BadSide');
end

function check_ballast_at_z_lo(body, label, reg, bp, so)
% contract S3, S4: at z_ballast = inner z_lo no inner face lies at that height; the ballast_top
% faces are the only faces in that plane
zb = body.design.z_ballast;
faces = body.brep.faces;
precast = strcmp(body.design.mode, 'modular_precast');
at = [];
inner_at = false;
for k = 1:numel(faces)
    s = body.brep.surfaces{faces(k).surface};
    if strcmp(s.type, 'plane')
        const = s.origin(3) == zb;
    else
        const = all(all(s.ctrl(:, :, 3) == zb));
    end
    if const
        at(end + 1) = k; %#ok<AGROW>
        inner_at = inner_at || strcmp(faces(k).role, 'inner');
    end
end
check(~inner_at, '%s: an inner face lies at z_ballast', label);
roles = {faces(at).role};
check(all(strcmp(roles, 'ballast_top')), '%s: faces at z_ballast: %s', label, strjoin(roles, ' '));
loops = arrayfun(@(k) numel(faces(k).loops), at);
pair = arrayfun(@(k) [faces(k).inside '/' faces(k).outside], at, 'UniformOutput', false);
if precast
    check(numel(at) == 1 && strcmp(pair{1}, 'uhpc/air') && loops == 1, '%s: expected one uhpc/air disk, got %s', label, strjoin(pair, ' '));
else
    [~, o] = sort(pair);
    check(numel(at) == 2 && isequal(pair(o), {'ballast/air', 'ballast/shell'}) && isequal(loops(o), [1 2]), ...
        '%s: expected a ballast/shell annulus and a ballast/air disk, got %s', label, strjoin(pair, ' '));
end
% the lateral inner faces of module 1 start at z_ballast and are whole (one face per patch)
inn = find(arrayfun(@(f) strcmp(f.role, 'inner') && f.module(1) == 1, faces));
lo = arrayfun(@(k) min(body.brep.surfaces{faces(k).surface}.ctrl(:, 1, 3)), inn);
hi = arrayfun(@(k) max(body.brep.surfaces{faces(k).surface}.ctrl(:, 1, 3)), inn);
check(numel(inn) == 4 && all(lo == zb) && all(hi == body.design.edges(2)), ...
    '%s: inner faces of module 1 (%d) do not start at z_ballast (%.17g, lo %s)', label, numel(inn), zb, mat2str(lo, 17));
% region volumes of module 1: full section below z_ballast, wall of section A - A_in above
e = body.design.edges;
si = sti_closed_form('section', body.analytic.fixture, body.analytic.d(1));
bound = 16 * eps * so.A * (abs(zb) + abs(e(1)) + abs(e(2)));
below = so.A * (zb - e(1));
if precast
    check(abs(bp.modules(1).V_uhpc - (below + (so.A - si.A) * (e(2) - zb))) <= bound, '%s: uhpc volume of module 1', label);
    check(abs(bp.modules(1).V_air - si.A * (e(2) - zb)) <= bound, '%s: air volume of module 1', label);
else
    check(abs(bp.modules(1).V_ballast - below) <= bound, '%s: ballast volume of module 1', label);
    check(abs(bp.modules(1).V_shell - (so.A - si.A) * (e(2) - zb)) <= bound, '%s: shell volume of module 1', label);
    check(abs(bp.modules(1).V_air - si.A * (e(2) - zb)) <= bound, '%s: air volume of module 1', label);
end
fprintf('%-70s faces at z_ballast: %s\n', label, strjoin(pair, ', '));
end

function [worst, n] = check_sides(body)
% F6b at the heights where the side matters, against the layout of the closed form
fx = body.analytic.fixture;
e = body.design.edges(:);
N = numel(e) - 1;
d = body.analytic.d;
lay = sti_closed_form('layout', fx, body.design, d);
zs = unique([e; body.design.z_ballast]);
zs = zs(zs >= e(1) & zs <= e(end))';
worst = 0;
n = 0;
for z = zs
    for side = {'above', 'below'}
        s = side{1};
        sec = mwecmass.solid.body_section(body, z, s);
        if strcmp(s, 'above')
            m = find(e(1:N) <= z & z < e(2:N + 1), 1);
            if isempty(m), m = N; end
        else
            m = find(e(1:N) < z & z <= e(2:N + 1), 1);
            if isempty(m), m = 1; end
        end
        check(sec.module == m, 'z = %g %s: module %d, expected %d', z, s, sec.module, m);
        if strcmp(s, 'above')
            air = lay(m).air && lay(m).a <= z && z < lay(m).b;
        else
            air = lay(m).air && lay(m).a < z && z <= lay(m).b;
        end
        check(sec.solid == ~air, 'z = %g %s: solid flag', z, s);
        so = sti_closed_form('section', fx, 0);
        check(abs(sec.outer.area - so.A) <= 64 * eps * so.A, 'z = %g %s: outer area', z, s);
        if air
            si = sti_closed_form('section', fx, d(m));
            check(abs(sec.inner.area - si.A) <= 64 * eps * si.A, 'z = %g %s: inner area', z, s);
        else
            check(isempty(sec.inner), 'z = %g %s: inner loop of a solid section', z, s);
        end
        n = n + 1;
    end
    a = mwecmass.solid.body_section(body, z);
    b = mwecmass.solid.body_section(body, z, 'above');
    check(isequal(a, b), 'z = %g: the default side is not above', z);
    below = mwecmass.solid.body_section(body, z, 'below');
    worst = max(worst, abs(a.outer.area - below.outer.area) / a.outer.area);
    pa = sortrows(a.outer.pts);
    pb = sortrows(below.outer.pts);
    check(isequal(pa, pb), 'z = %g: outer loops of the two sides differ', z);
end
end

function n = shared_edges(brep)
% I2: the curve of every edge in a B-spline face's loop equals one boundary row of that face,
% bitwise (control points, knots, weights), in either direction
n = 0;
uses = zeros(1, numel(brep.edges));
for k = 1:numel(brep.faces)
    s = brep.surfaces{brep.faces(k).surface};
    for L = brep.faces(k).loops
        uses(abs(L{1})) = uses(abs(L{1})) + 1;
    end
    if ~strcmp(s.type, 'bspline')
        continue
    end
    W = s.weights;
    rows = {struct('ctrl', squeeze(s.ctrl(:, 1, :)), 'knots', s.knots{1}, 'weights', wcol(W, 1, 'c')), ...
            struct('ctrl', squeeze(s.ctrl(end, :, :)), 'knots', s.knots{2}, 'weights', wcol(W, size(W, 1), 'r')), ...
            struct('ctrl', squeeze(s.ctrl(:, end, :)), 'knots', s.knots{1}, 'weights', wcol(W, size(W, 2), 'c')), ...
            struct('ctrl', squeeze(s.ctrl(1, :, :)), 'knots', s.knots{2}, 'weights', wcol(W, 1, 'r'))};
    for e = abs(brep.faces(k).loops{1})
        c = brep.curves(brep.edges(e).curve);
        ok = false;
        for b = 1:4
            r = rows{b};
            ok = ok || (isequal(c.ctrl, r.ctrl) && isequal(c.knots, r.knots) && isequal(c.weights, r.weights)) || ...
                (isequal(c.ctrl, flipud(r.ctrl)) && isequal(c.weights, flipud(r.weights)) && ...
                 isequal(c.knots, r.knots(1) + r.knots(end) - r.knots(end:-1:1)));
        end
        check(ok, 'face %d edge %d: curve is not a boundary row of the face', k, e);
    end
end
n = sum(uses >= 2);
end

function w = wcol(W, i, kind)
if isempty(W)
    w = [];
elseif kind == 'c'
    w = W(:, i);
else
    w = W(i, :)';
end
end

function delete_if(f)
if exist(f, 'file')
    delete(f);
end
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'step'));
addpath(fullfile(root, 'tests', 'standins'), '-end');
addpath(fullfile(root, 'tests', 'standins', 'fixtures'), '-end');
end

function expect_error(f, id)
try
    f();
catch err
    check(strcmp(err.identifier, id), 'expected %s, got %s', id, err.identifier);
    return
end
error('test_sk_body:fail', 'expected error %s', id);
end

function check(cond, varargin)
if ~cond
    error('test_sk_body:fail', varargin{:});
end
end
