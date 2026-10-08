function test_precast_stage3_known_optimum()
%TEST_PRECAST_STAGE3_KNOWN_OPTIMUM  Stage 3 reaches a known optimum on the equality manifold.
%   Kernel F1-F7 are the SK stand-ins until J1 (box inner sets through sti_inner_box). The Stage-2
%   record is the final_props of a buildable design x* (t = [0.14; 0.135; 0.0762] m, z_ballast at
%   flotation) at the same draft, so x* holds both equalities and its objective is 0 (independent
%   oracle: x* itself). The Stage-2 densities are those of case C of test_precast_stage3_solve, so
%   the split differs from x* and the optimisation runs.
%   Asserted: the split check fails and the optimisation runs at the Stage-2 draft only; the stored
%   design holds both equalities (AGENTS OD10, tolerance 1e-6) and passes the check (accepted);
%   phase 2 ran and the stored objective is no larger than the objective where it started (it
%   only keeps points that improve it); F7 is called once while the draft is fixed (contract
%   section 7 item 3). The objective and the residuals are printed.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
tol_eq = 1e-6;
config = sti_config('box');
config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
config.wall_strip_index = [];
config.constructability_rho_hull = 2500;
config.constructability_rho_air = 1.2;
config.constructability_t_min = 0.0762;
config.output.save.stage3 = struct('precast_midplane', false, 'precast_strips', false, 'step', false);
geo = config.hull_solid;
t_min = config.constructability_t_min;
e = config.strip_edges(:);
vs = 0;
[f3, s2] = sti_stage2(config, vs, [NaN; 600; 600]);

% x*: shells fixed, ballast level in module 1 at flotation
ctx = struct('geo', geo, 'rho', struct('uhpc', 2500, 'air', 1.2), 'z_hollow', geo.z_range, ...
    't_max_hollow', Inf, 'knots_from', [], 'sets', [], 'set_ranges', zeros(0, 2), ...
    'inner_fn', @(t, knots_from, zr) sti_inner_box(geo, t, t_min));
hs = mwecmass.solid.hydrostatics_at_draft(geo, vs, struct());
xs = struct('mode', 'modular_precast', 'edges', e, 'vs', vs, 't', [0.14; 0.135; 0.0762], ...
    'z_ballast', e(1), 'solid_modules', zeros(1, 0));
gap = @(z) body_mass(ctx, setfield(xs, 'z_ballast', z)) - config.RHO_WATER * hs.V_sub; %#ok<SFLD>
xs.z_ballast = fzero(gap, [e(1), e(2)], optimset('Display', 'off'));
ev = mwecmass.realise.modular_precast.realise_modules(ctx, xs);
p = mwecmass.realise.evaluate_realised(ev.bp, hs, xs, config);
fprintf('x*: z_ballast %.9f m, flotation residual %.3g, Z_CG %.6f m, GM %.6f m, T_pitch %.6f s\n', ...
    xs.z_ballast, p.mass_total / (config.RHO_WATER * hs.V_sub) - 1, p.CG_total(3), p.GM_L, ...
    p.periods.pitch);
f3.mass_total = p.mass_total;
f3.CG_total = p.CG_total;
f3.GM_L = p.GM_L;
f3.periods = p.periods;

opts = struct('inner_fn', @(t, knots_from, z_range) sti_inner_box(geo, t, t_min), 'd_close_fn', []);
profile clear
profile on
r = mwecmass.realise.modular_precast.solve_and_extract(config, [vs; s2.rho], f3, opts);
profile off
info = profile('info');
names = {info.FunctionTable.FunctionName};
calls = [info.FunctionTable.NumCalls];
n_hs = sum(calls(~cellfun(@isempty, regexp(names, 'hydrostatics_at_draft$', 'once'))));

steps = {r.solver.step};
fprintf('steps %s; status %s; z_ballast %.9f m; t %s mm; vs %.6f m\n', strjoin(steps, ' > '), ...
    r.status, r.design.z_ballast, mat2str(1000 * r.design.t', 9), r.vs);
for s = r.solver
    fprintf('step %-11s exitflag %g, iterations %d, objective %.6g (phase-2 start %.6g), max equality violation %.3g\n', ...
        s.step, s.exitflag, s.iterations, s.fval, s.fval_phase2_start, s.max_eq_violation);
end
fprintf('residuals flotation %.3g, GM %.3g; rel. deviations Z_CG %.4g, GM %.4g, T_heave %.4g, T_pitch %.4g\n', ...
    r.check.equalities.residual, r.check.metrics.rel_dev);
fprintf('hydrostatics_at_draft calls: %d\n', n_hs);
check(isequal(steps, {'split', 'fixed_draft'}), 'the split failed the check; only the fixed-draft step ran');
check(all(abs([r.check.equalities.residual]) <= tol_eq), 'both equalities hold');
check(strcmp(r.status, 'accepted') && r.check.pass, 'accepted');
check(isfinite(r.solver(2).fval_phase2_start) && r.solver(2).fval <= r.solver(2).fval_phase2_start, ...
    'phase 2 ran and kept only improving points');
check(n_hs == 1, 'F7 called once at the fixed draft');
end

function M = body_mass(ctx, design)
ev = mwecmass.realise.modular_precast.realise_modules(ctx, design);
M = ev.bp.total.mass;
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

function check(cond, msg)
if ~cond
    error('test_precast_stage3_known_optimum:fail', '%s', msg);
end
end
