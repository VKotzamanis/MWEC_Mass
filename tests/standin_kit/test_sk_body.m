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
    };
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
fprintf('largest relative volume closure error per module %.3e, sum of modules vs hull %.3e\n', worst.closure, worst.hull);
fprintf('largest |OCC getMass - closed form| / closed form %.3e (printed, not asserted)\n', worst.occ);
fprintf('largest relative difference of F6b void areas from the closed form %.3e\n', worst.section);

body = mwecmass.solid.build_body(gc, struct('mode', 'thin_shell', 'edges', ec, 'vs', 0.5, ...
    't', 0.1 * ones(4, 1), 'z_ballast', -2.5, 'solid_modules', []), ic);
sec = mwecmass.solid.body_section(body, -2.7);
check(sec.solid && sec.module == 1, 'section below the ballast top is solid');
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

function check(cond, varargin)
if ~cond
    error('test_sk_body:fail', varargin{:});
end
end
