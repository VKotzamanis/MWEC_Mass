function test_thin_shell_stage3()
%TEST_THIN_SHELL_STAGE3  Thin-shell Stage 3 (contract F14) on the cylinder fixture.
%   Runs mwecmass.realise.thin_shell.run on the SK stand-ins of F1-F7, F9, F10 (the real functions
%   win once merged; this is the join test) for these Stage-2 solutions: one built from a
%   thin-shell design (draft kept; run with two acceptance bands), one with an unreachable GM
%   (closest fail at the Stage-2 draft), one whose displaced mass at the Stage-2 draft is below
%   the lightest buildable design (draft released), one whose ballast spills past module 1, one
%   that ends the flotation curve with z_ballast at the hull bottom (below the inner z_lo), and one
%   that no draft within the bounds floats (closest fail with the draft released). The figure and
%   STEP calls run against test doubles that record them.
%   The closed forms of sti_closed_form (prisms) are the independent oracle of the hull volume;
%   the prism formulas of the cylinder (written out below) are the oracle of the no-float case.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
stub_dir = make_stubs();
cleanup = onCleanup(@() remove_stubs(stub_dir));
global T7_STUB_CALLS %#ok<GVMIS>

config = sti_config('cylinder');
config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
config.rho_shell = 7500;
config.rho_ballast = 7500;
config.rho_air = 1.2;
config.steel_t_min = 0.0254;
config.steel_t_init = 0.0254;
config.vertical_shift_bounds = [-0.9 2.9];
config.output = struct('console_echo', true, 'save', struct('stage3', ...
    struct('steel_solve_log', false, 'steel_solve', true, 'step', true)));
rho = struct('ballast', config.rho_ballast, 'shell', config.rho_shell, 'air', config.rho_air);
geo = config.hull_solid;
fx = geo.analytic;
N = numel(config.strip_edges) - 1;
eps_fit = 0.01 * config.steel_t_min;
t_max = mwecmass.solid.void_closing_distance(config.ms2_model, [], geo, geo.z_range) - eps_fit / 2;
V_hull = sti_closed_form('hull', fx).V;
fprintf('t_max = d_close - eps_fit/2 = %.6f m\n', t_max);

% Stage 2 from a thin-shell design at vs = 0 (t = 30 mm, z_ballast floating it): its module
% densities as strip densities; the strip model puts each module's mass at uniform density, so
% Z_CG and Iyy of Stage 2 differ from the design's and Stage 3 has to move t and z_ballast.
vs = 0;
t_design = 0.03;
zb_design = fzero(@(z) design_flotation(config, rho, vs, t_design, z), [fx.z(1), fx.z(2)]);
[~, st] = sti_stage2(config, vs, [NaN; 100; 100; 100]);
r_design = sti_realised(config, thin_design(config, vs, t_design, zb_design), rho, st, ...
    struct('t_min', config.steel_t_min));
rho2 = [r_design.modules.rho_eff]';
% Case 1 runs twice: with the author's mass_acceptable_pct (10 %) and with 20 %, a test input that
% exercises the accepted branch (the strip model's T_pitch differs from the design's by more
% than 10 % here).
% Cases 5 and 6 use a light shell and a denser ballast (rho_air < rho_shell <= rho_ballast), so
% flotation needs ballast above module 1 (case 5) or no ballast at a thick shell (case 6).
steel = [7500, 7500];
light = [1100, 1500];
wide = [-0.9 2.9];
cases = struct('label', {}, 'vs', {}, 'rho2', {}, 'gm_scale', {}, 'pct', {}, 'escalation', {}, ...
    'dens', {}, 'bounds', {});
cases(1) = struct('label', 'Stage 2 from a thin-shell design', 'vs', vs, 'rho2', rho2, ...
    'gm_scale', 1, 'pct', 10, 'escalation', 'fixed_draft', 'dens', steel, 'bounds', wide);
cases(2) = struct('label', 'Stage 2 from a thin-shell design, 20 % band', 'vs', vs, 'rho2', rho2, ...
    'gm_scale', 1, 'pct', 20, 'escalation', 'fixed_draft', 'dens', steel, 'bounds', wide);
cases(3) = struct('label', 'unreachable GM (5 x GM_Stage2)', 'vs', vs, 'rho2', rho2, ...
    'gm_scale', 5, 'pct', 10, 'escalation', 'fixed_draft', 'dens', steel, 'bounds', wide);
cases(4) = struct('label', 'mass unreachable at the Stage-2 draft', 'vs', 2.4, ...
    'rho2', [NaN; 50; 50; 50], 'gm_scale', 1, 'pct', 10, 'escalation', 'draft_free', ...
    'dens', steel, 'bounds', wide);
cases(5) = struct('label', 'ballast spilled past module 1', 'vs', 0, 'rho2', [NaN; 900; 100; 50], ...
    'gm_scale', 1, 'pct', 10, 'escalation', 'fixed_draft', 'dens', light, 'bounds', wide);
cases(6) = struct('label', 'end of the flotation curve (z_ballast at the hull bottom)', 'vs', -0.9, ...
    'rho2', [NaN; 50; 50; 50], 'gm_scale', 0.05, 'pct', 10, 'escalation', 'fixed_draft', ...
    'dens', light, 'bounds', wide);
cases(7) = struct('label', 'no draft within the bounds floats', 'vs', 2.4, 'rho2', [NaN; 50; 50; 50], ...
    'gm_scale', 1, 'pct', 10, 'escalation', 'draft_free', 'dens', steel, 'bounds', [2.0 2.9]);

for c = 1:numel(cases)
    cs = cases(c);
    fprintf('\n== %s\n', cs.label);
    config.rho_shell = cs.dens(1);
    config.rho_ballast = cs.dens(2);
    config.vertical_shift_bounds = cs.bounds;
    rho = struct('ballast', config.rho_ballast, 'shell', config.rho_shell, 'air', config.rho_air);
    [f3, ~] = sti_stage2(config, cs.vs, cs.rho2);
    f3.GM_L = cs.gm_scale * f3.GM_L;
    rho_x = f3.densities_at_nodes;
    config.mass_acceptable_pct = cs.pct;
    T7_STUB_CALLS = {};
    opt_results = struct('Final3D', f3, 'stage2_3d', struct('properties', f3));
    [results, final_props] = mwecmass.realise.thin_shell.run(config, [cs.vs; rho_x], opt_results);
    r = results.stage3;

    % S8 and final_props describe the realised body, never Stage 2
    check(isequaln(results.Final3D, f3), 'results.Final3D changed');
    check(isequaln(results.stage2_3d.properties, final_props), 'stage2_3d.properties is not final_props');
    check(isequaln(rmfield(final_props, {'stage3_status', 'stage3_check'}), r.props), 'final_props is not S8 props');
    check(strcmp(final_props.stage3_status, r.status) && isequaln(final_props.stage3_check, r.check), ...
        'final_props status fields');
    check(strcmp(r.mode, 'thin_shell') && strcmp(r.hull_name, 'cylinder') && isempty(r.k_star) && ...
        isempty(r.V_uhpc_target), 'S8 header');
    check(strcmp(r.escalation, cs.escalation), 'escalation %s, expected %s', r.escalation, cs.escalation);
    check(strcmp(r.status, 'accepted') == r.check.pass, 'status does not follow the check');
    check(strcmp(r.status, 'accepted') == isempty(r.reason), 'reason must be empty exactly when accepted');
    check(isequal(r.stage2.vs, cs.vs) && isequal(r.stage2.rho, rho_x) && r.stage2.GM == f3.GM_L && ...
        r.stage2.Z_CG == f3.CG_total(3) && r.stage2.T_pitch == f3.periods.pitch, 'S8.stage2');
    if strcmp(cs.escalation, 'fixed_draft')
        check(r.vs == cs.vs && r.design.vs == cs.vs, 'the draft moved although mass balance was reachable');
    end
    check(numel(r.solver) == 1 && strcmp(r.solver.step, cs.escalation), 'solver entries');
    check(r.solver.max_eq_violation >= 0, 'solver max_eq_violation');

    % one uniform shell thickness within [t_min, t_max]; z_ballast inside the hull
    t = r.design.t;
    check(numel(t) == N && all(t == t(1)) && t(1) >= config.steel_t_min && t(1) <= t_max, 't = %g', t(1));
    zb = r.design.z_ballast;
    check(zb >= fx.z(1) && zb <= fx.z(2) && isempty(r.design.solid_modules), 'z_ballast %g', zb);

    % identity: the reported design re-evaluated through the stand-in kit's S8 builder
    ref = sti_realised(config, r.design, rho, r.stage2, struct('t_min', config.steel_t_min, ...
        'escalation', r.escalation));
    check(isequaln(ref.props, r.props) && isequaln(ref.check, r.check), 'props differ from sti_realised');
    check(isequaln(ref.modules, r.modules), 'modules differ from sti_realised');

    % volume closure (rule 11): each module V is the sum of its three region volumes, rounding of
    % three terms; the hull is four closed-form prisms of the same section, rounding of four terms
    V = [r.modules.V];
    for i = 1:N
        m = r.modules(i);
        check(abs(m.V_ballast + m.V_shell + m.V_air - m.V) <= 4 * eps * m.V, 'module %d closure', i);
        check(m.rho_stage2 == rho_x(i), 'module %d rho_stage2', i);
    end
    check(abs(sum(V) - V_hull) <= 8 * eps * V_hull, 'sum of module volumes %.17g vs hull %.17g', sum(V), V_hull);

    % figure and STEP called for every status, with the realised design
    check(numel(T7_STUB_CALLS) == 2 && strcmp(T7_STUB_CALLS{1}{1}, 'plot_steel_solve') && ...
        strcmp(T7_STUB_CALLS{2}{1}, 'export_stage3'), 'F13 and F11 calls');
    check(isequaln(T7_STUB_CALLS{1}{2}.props, r.props) && isequaln(T7_STUB_CALLS{2}{2}.props, r.props), ...
        'F13/F11 got another design');
    check(numel(r.step_files) == 1 && strcmp(r.step_files.name, 'stub'), 'step_files not stored');

    fl = r.check.equalities(1);
    gm = r.check.equalities(2);
    fprintf(['%s: %s; t %.6f m, z_ballast %.6f m, vs %.4f m; flotation residual %.3e, GM residual %.3e; ' ...
        'rel. dev. Z_CG %.3g, GM %.3g, T_heave %.3g, T_pitch %.3g; solver exitflag %d, %d iterations\n'], ...
        cs.label, r.status, t(1), zb, r.vs, fl.residual, gm.residual, r.check.metrics.rel_dev, ...
        r.solver.exitflag, r.solver.iterations);
    fprintf('module rho_eff %s kg/m^3; Stage 2 %s\n', mat2str([r.modules.rho_eff], 6), mat2str(rho_x', 6));

    switch c
        case {1, 2}
            % two unknowns, two equalities at a reachable draft: the step solves both
            check(fl.pass && gm.pass, 'equalities not met at the Stage-2 draft');
            fprintf('design t %.6f m, z_ballast %.6f m\n', t_design, zb_design);
            if c == 2
                check(strcmp(r.status, 'accepted'), 'not accepted within the 20 %% band');
                d1 = results_1.stage3.design;
                check(isequal(r.design, d1), 'the band changed the design');
            else
                results_1 = results;
            end
        case 3
            check(strcmp(r.status, 'failed') && any(strcmp(r.check.failed, 'GM')) && ...
                ~isempty(strfind(r.reason, 'GM')), 'unreachable GM must end as a closest fail naming GM');
            check(fl.pass, 'the closest design must keep mass balance when it can');
        case 4
            check(r.vs ~= cs.vs, 'draft not released');
        case 5
            e = config.strip_edges;
            m1 = r.modules(1);
            check(fl.pass, 'flotation not held with ballast spilled');
            check(zb > e(2) && zb < e(3), 'z_ballast %.6f not in module 2', zb);
            check(isnan(m1.t) && m1.h_ballast == e(2) - e(1) && m1.V_ballast == m1.V && ...
                m1.V_shell == 0 && m1.V_air == 0, 'module 1 below z_ballast is not all ballast');
            check(r.modules(2).h_ballast == zb - e(2) && r.modules(3).h_ballast == 0, 'h_ballast');
        case 6
            inner = mwecmass.solid.offset_surface(config.ms2_model, [], geo, t(1), geo.z_range, ...
                struct('t_min', config.steel_t_min));
            fprintf('inner z_lo %.6f m\n', inner.z_lo);
            check(zb == fx.z(1) && zb <= inner.z_lo, 'z_ballast %.6f is not the hull bottom', zb);
            check(r.solver.exitflag == 0 && fl.pass, 'curve end: exitflag %d, flotation %g', ...
                r.solver.exitflag, fl.residual);
            check(all([r.modules.V_ballast] == 0), 'ballast volume below the hull bottom');
        case 7
            % too heavy at every draft in the bounds: the smallest flotation violation is the
            % lightest design (t_min, no ballast) at the deepest allowed draft (vs = 2.0)
            d = config.steel_t_min + eps_fit / 2;
            R = fx.R;
            H = fx.z(2) - fx.z(1);
            M_lo = rho.shell * pi * R^2 * H - (rho.shell - rho.air) * pi * (R - d)^2 * (H - 2 * d);
            V_sub = pi * R^2 * (-cs.bounds(1) - fx.z(1));
            r_cf = M_lo / (config.RHO_WATER * V_sub) - 1;
            % M_lo is the difference of two terms of size rho_shell V_hull; their rounding,
            % relative to M_lo, bounds the error of the residual (a few roundings each)
            bound = 16 * eps * rho.shell * V_hull / M_lo * (1 + r_cf);
            fprintf('flotation residual %.15g, closed form %.15g, difference %.3g (bound %.3g)\n', ...
                fl.residual, r_cf, fl.residual - r_cf, bound);
            check(strcmp(r.status, 'failed') && ~fl.pass, 'no-float case must fail on flotation');
            check(r.vs == cs.bounds(1) && t(1) == config.steel_t_min && zb == fx.z(1), ...
                'closest design (%g, %g, %g) is not the deepest draft with the lightest build', r.vs, t(1), zb);
            check(abs(fl.residual - r_cf) <= bound, 'flotation residual differs from the closed form');
    end
end
end

function r = design_flotation(config, rho, vs, t, zb)
inner = mwecmass.solid.offset_surface(config.ms2_model, [], config.hull_solid, t, config.hull_solid.z_range, ...
    struct('t_min', config.steel_t_min));
body = mwecmass.solid.build_body(config.hull_solid, thin_design(config, vs, t, zb), inner);
bp = mwecmass.solid.body_properties(body, rho, struct());
hs = mwecmass.solid.hydrostatics_at_draft(config.hull_solid, vs, struct());
r = bp.total.mass / (config.RHO_WATER * hs.V_sub) - 1;
end

function d = thin_design(config, vs, t, zb)
N = numel(config.strip_edges) - 1;
d = struct('mode', 'thin_shell', 'edges', config.strip_edges, 'vs', vs, 't', t * ones(N, 1), ...
    'z_ballast', zb, 'solid_modules', []);
end

function dir_ = make_stubs()
% Test doubles of F13 and F11 (owned by T8 and T10) placed before src on the path; they record
% their calls in a global.
dir_ = tempname();
mkdir(fullfile(dir_, '+mwecmass', '+output', '+figures'));
mkdir(fullfile(dir_, '+mwecmass', '+output', '+step'));
write_file(fullfile(dir_, '+mwecmass', '+output', '+figures', 'plot_steel_solve.m'), { ...
    'function fig = plot_steel_solve(realised, config) %#ok<INUSD>', ...
    'global T7_STUB_CALLS', ...
    'T7_STUB_CALLS{end + 1} = {''plot_steel_solve'', realised};', ...
    'fig = [];', 'end'});
write_file(fullfile(dir_, '+mwecmass', '+output', '+step', 'export_stage3.m'), { ...
    'function files = export_stage3(realised, out_dir)', ...
    'global T7_STUB_CALLS', ...
    'T7_STUB_CALLS{end + 1} = {''export_stage3'', realised, out_dir};', ...
    'files = struct(''name'', ''stub'', ''path'', out_dir, ''bodies'', {{}});', 'end'});
addpath(dir_);
end

function remove_stubs(dir_)
rmpath(dir_);
if exist('OCTAVE_VERSION', 'builtin')
    confirm_recursive_rmdir(false);
end
rmdir(dir_, 's');
end

function write_file(path, lines)
fid = fopen(path, 'w');
fprintf(fid, '%s\n', lines{:});
fclose(fid);
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'standins'), '-end');
addpath(fullfile(root, 'tests', 'standins', 'fixtures'), '-end');
end

function check(cond, varargin)
if ~cond
    error('test_thin_shell_stage3:fail', varargin{:});
end
end
