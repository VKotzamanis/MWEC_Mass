function test_precast_stage3_solve()
%TEST_PRECAST_STAGE3_SOLVE  Modular-precast Stage-3 optimisation, escalation and closest fail.
%   Kernel F1-F7 are the SK stand-ins until J1 (box inner sets through sti_inner_box, cylinder
%   inner sets through F2); Stage 2 is sti_stage2 (closed-form prisms), not an optimiser run.
%   Cases: (C) a split that fails the check and is optimised; (S) a ballast module too light even
%   when full with the thickest shells: the ballast spills into k*+1; (D) a hull too heavy at the
%   Stage-2 draft for every design: the draft is released last; (G) an unreachable GM
%   (forced infeasible) through F14 run.m on the cylinder: the closest fail is stored in
%   results.stage3 and final_props, flagged and reported (figures and STEP are switched off in
%   tests, contract section 3). The out-of-reach conditions of S and D are checked on mass bounds
%   built here: the mass rises with z_ballast and with every t (rho_uhpc > rho_air).
%   Asserted (exact by construction): the escalation order and the conditions that start each
%   step, the status rule (accepted exactly when F10 passes), bounds of every variable, the draft
%   kept unless the draft step ran, props and final_props from the stored body, the stored solver
%   record of the last step equal to the objective and residuals of the stored props (solve returns
%   the props and check of the evaluation that set solver(end)), volume closure (rule 11).
%   Outcomes (AGENTS section 3 items 27, 31, 32; OD10 tolerance 1e-6): C accepted with both
%   equalities and an objective no worse than a known feasible design; S and D hold flotation (D
%   at the released draft: the lightest design floats inside the vs bounds); G holds flotation
%   with GM out of reach. Residuals and deviations are printed.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
tol_eq = 1e-6;
order = {'split', 'fixed_draft', 'spill', 'draft_free'};

% C: box, three modules, no wall; the split fails the GM check, the optimisation follows
config = fixture_config('box', [], 2500, []);
[f3, s2] = sti_stage2(config, 0, [NaN; 600; 600]);
r = stage3(config, [0; s2.rho], f3, []);
common(r, [], config, s2, order, tol_eq, 'C');
check(numel(r.solver) >= 2, 'C: the split check failed and the optimisation ran');
% a design within the bounds that holds both equalities and passes the check (grader, round 1:
% t = [0.139902 0.135073 0.0762] m, z_ballast = -2.218129 m) has objective 0.00839728
check(strcmp(r.status, 'accepted') && all(abs([r.check.equalities.residual]) <= tol_eq), ...
    'C: accepted with both equalities held');
check(r.solver(end).fval <= 0.00839728, 'C: objective no worse than the known feasible design');

% S: two modules, shells capped at t_max = 0.1 m - eps_fit/2 (d_close provider of the test): k* is
% too light even when full of ballast, so flotation needs the spill into module 2
config = fixture_config('box', [], 1200, [-2.5; -1; 0.5]);
d_close_fn = @(zr) 0.1;
t_max = 0.1 - config.constructability_t_min / 200;
[f3, s2] = sti_stage2(config, 0, [1100; NaN]);
M_max = mass_of(config, struct('t', [NaN; t_max], 'z_ballast', -1));
fprintf('S: heaviest design with the ballast inside k* %.3f kg, displaced mass %.3f kg\n', M_max, ...
    f3.mass_buoyant_force);
check(M_max < f3.mass_buoyant_force, 'S: flotation out of reach inside k* (mass rises with z_ballast and t)');
r = stage3(config, [0; s2.rho], f3, d_close_fn);
common(r, [], config, s2, order, tol_eq, 'S');
check(any(strcmp({r.solver.step}, 'spill')), 'S: spill ran');
check(all(r.design.t(isfinite(r.design.t)) <= t_max), 'S: shells within t_max');
check(r.check.equalities(1).pass, 'S: flotation holds after the spill');

% D: two modules, too heavy at the Stage-2 draft for every design: the draft is released last
config = fixture_config('box', [], 2500, [-2.5; -1; 0.5]);
[f3, s2] = sti_stage2(config, 2.4, [100; 100]);
M_min = mass_of(config, struct('t', config.constructability_t_min * [1; 1], 'z_ballast', -2.5));
fprintf('D: lightest design %.3f kg, displaced mass %.3f kg\n', M_min, f3.mass_buoyant_force);
check(M_min > f3.mass_buoyant_force, 'D: flotation out of reach at the Stage-2 draft');
hull = sti_closed_form('hull', sti_closed_form('fixture', config.hull_solid));
fprintf('D: rho_w V_hull %.3f kg\n', config.RHO_WATER * hull.V);
check(M_min < config.RHO_WATER * hull.V, 'D: some draft floats the lightest design');
r = stage3(config, [2.4; s2.rho], f3, []);
common(r, [], config, s2, order, tol_eq, 'D');
check(strcmp(r.escalation, 'draft_free') && r.vs ~= s2.vs, 'D: draft released');
fprintf('D: vs %.6f m (Stage 2 %.6f m), flotation residual %.3g, GM residual %.3g\n', r.vs, s2.vs, ...
    r.check.equalities.residual);
check(r.check.equalities(1).pass, 'D: flotation holds at the released draft');

% G: unreachable GM (forced infeasible, GM_Stage2 = 100 m) through F14 run.m on the cylinder
config = fixture_config('cylinder', [], 2500, [-3; -1; 1]);
vs = 0.5;
[f3, s2] = sti_stage2(config, vs, [NaN; 300]);
f3.GM_L = 100;
s2.GM = 100;
opt = struct('Final3D', f3, 'stage2_3d', struct('properties', f3), 'constructability', []);
[results, fp] = mwecmass.realise.modular_precast.run(config, [vs; s2.rho], opt);
r = results.stage3;
check(isequaln(results.Final3D, f3) && isequaln(results.stage2_3d.properties, fp), 'G: F14 results layout');
common(r, fp, config, s2, order, tol_eq, 'G');
check(strcmp(r.status, 'failed') && any(strcmp(r.check.failed, 'GM')) && ...
    ~isempty(strfind(r.reason, 'closest fail')) && strcmp(fp.stage3_status, 'failed'), ...
    'G: closest fail flagged and reported');
check(numel(r.solver) >= 3 && r.solver(2).max_eq_violation > tol_eq, ...
    'G: the GM equality failed at the Stage-2 draft, so the spill ran');
check(r.check.equalities(1).pass && ~any(strcmp({r.solver.step}, 'draft_free')), ...
    'G: flotation held at the Stage-2 draft, so the draft is kept');
end

function M = mass_of(config, d)
% Mass of a box design at module edges config.strip_edges (stand-in kernel, sti_inner_box sets).
geo = config.hull_solid;
t_min = config.constructability_t_min;
ctx = struct('geo', geo, 'rho', struct('uhpc', config.constructability_rho_hull, 'air', 1.2), ...
    'z_hollow', geo.z_range, 't_max_hollow', Inf, 'knots_from', [], 'sets', [], ...
    'set_ranges', zeros(0, 2), 'inner_fn', @(t, knots_from, zr) sti_inner_box(geo, t, t_min));
design = struct('mode', 'modular_precast', 'edges', config.strip_edges, 'vs', 0, 't', d.t, ...
    'z_ballast', d.z_ballast, 'solid_modules', zeros(1, 0));
ev = mwecmass.realise.modular_precast.realise_modules(ctx, design);
M = ev.bp.total.mass;
end

function common(r, fp, config, s2, order, tol_eq, label)
steps = {r.solver.step};
check(isequal(steps, order(1:numel(steps))) && strcmp(r.escalation, steps{end}), ...
    [label ': steps in escalation order']);
bp = mwecmass.solid.body_properties(r.body, r.rho, struct());
check(isequaln(r.props.mass_total, bp.total.mass) && r.props.CG_total(3) == bp.total.CG_body(3) + r.vs, ...
    [label ': props from the stored body']);
if ~isempty(fp)
    check(isequal(rmfield(fp, {'stage3_status', 'stage3_check'}), r.props) && ...
        strcmp(fp.stage3_status, r.status) && isequal(fp.stage3_check, r.check), [label ': final_props']);
end
check(strcmp(r.status, 'accepted') == r.check.pass && isempty(r.reason) == r.check.pass, ...
    [label ': status rule']);
X = [r.props.CG_total(3), r.props.periods.heave, r.props.periods.pitch];
X2 = [s2.Z_CG, s2.T_heave, s2.T_pitch];
check(r.solver(end).fval == sum(((X - X2) ./ X2).^2), [label ': stored objective']);
check(r.solver(end).max_eq_violation == max(abs([r.check.equalities.residual])), ...
    [label ': stored equality violation']);
if ~any(strcmp(steps, 'draft_free'))
    check(r.vs == s2.vs && r.design.vs == s2.vs, [label ': draft kept']);
end
if numel(steps) >= 3
    check(r.solver(2).max_eq_violation > tol_eq, [label ': escalation only after failed equalities']);
end
if any(strcmp(steps, 'draft_free'))
    prev = r.solver(end - 1);
    check(prev.max_eq_violation > tol_eq, [label ': draft released only after failed equalities']);
end
e = r.design.edges;
k = r.k_star;
z_hi = e(k + 1);
if any(strcmp(steps, 'spill'))
    z_hi = e(k + 2);
end
check(r.design.z_ballast >= e(k) && r.design.z_ballast <= z_hi, [label ': ballast within its bound']);
t_min = config.constructability_t_min;
t = r.design.t;
check(all(t(isfinite(t)) >= t_min), [label ': shells >= t_min']);
check(all(isnan(t(e(2:end) <= r.design.z_ballast))), [label ': modules below the ballast level solid']);
closure(r, config, label);
fprintf(['%s: steps %s; status %s; z_ballast %.6f m; t %s mm; vs %.6f m; flotation residual %.3g; ' ...
    'GM residual %.3g\n'], label, strjoin(steps, ' > '), r.status, r.design.z_ballast, ...
    mat2str(1000 * t', 6), r.vs, r.check.equalities(1).residual, r.check.equalities(2).residual);
fprintf('%s: rel. deviations Z_CG %.4g, GM %.4g, T_heave %.4g, T_pitch %.4g\n', label, ...
    r.check.metrics.rel_dev);
for s = r.solver
    fprintf('%s: step %-11s exitflag %g, iterations %d, objective %.6g, max equality violation %.3g\n', ...
        label, s.step, s.exitflag, s.iterations, s.fval, s.max_eq_violation);
end
if r.check.equalities(1).pass
    return
end
fprintf('%s: flotation not met: %s\n', label, r.reason);
end

function closure(r, config, label)
Vsum = 0;
for i = 1:numel(r.modules)
    m = r.modules(i);
    % V is the sum of the two region volumes: equal to one rounding of that sum
    check(abs(m.V_uhpc + m.V_air - m.V) <= 2 * eps * m.V, sprintf('%s: module %d closure', label, i));
    Vsum = Vsum + m.V;
end
hull = sti_closed_form('hull', sti_closed_form('fixture', config.hull_solid));
fprintf('%s: sum of module volumes - closed-form hull volume = %.3g m^3\n', label, Vsum - hull.V);
end

function config = fixture_config(name, wall, rho_uhpc, edges)
config = sti_config(name);
if ~isempty(edges)
    config.strip_edges = edges;
end
config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
config.wall_strip_index = wall;
config.constructability_rho_hull = rho_uhpc;
config.constructability_rho_air = 1.2;
config.constructability_t_min = 0.0762;
config.output.save.stage3 = struct('precast_midplane', false, 'precast_strips', false, 'step', false);
end

function r = stage3(config, x, f3, d_close_fn)
t_min = config.constructability_t_min;
geo = config.hull_solid;
opts = struct('inner_fn', @(t, knots_from, z_range) sti_inner_box(geo, t, t_min), ...
    'd_close_fn', d_close_fn);
r = mwecmass.realise.modular_precast.solve_and_extract(config, x, f3, opts);
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
    error('test_precast_stage3_solve:fail', '%s', msg);
end
end
