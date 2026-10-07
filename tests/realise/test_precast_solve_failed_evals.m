function test_precast_solve_failed_evals()
%TEST_PRECAST_SOLVE_FAILED_EVALS  Kernel errors inside the Stage-3 optimisation are failed
%evaluations, never a stop of Stage 3 (contract F5 and section 8).
%   Box fixture with two modules on the SK stand-ins, module 2 solid (no spill: k*+1 is solid);
%   solve is called directly on a start design. The inner-set provider refits on fixed knots for
%   every t, but its adaptive fit raises mwecmass:solid:FitNotConverged for every t other than the
%   start's, so the rebuild of every optimum away from the start fails. Asserted: solve returns,
%   every step it runs is recorded in escalation order, the spill is skipped with a note, failed
%   evaluations are noted, and the stored design uses only a thickness whose adaptive fit exists.
%   Numbers are printed.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
config = sti_config('box');
config.strip_edges = [-2.5; -1; 0.5];
config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
geo = config.hull_solid;
e = config.strip_edges;
t_min = 0.0762;
t0 = 0.1;
ctx = struct('geo', geo, 'rho', struct('uhpc', 2500, 'air', 1.2), 'z_hollow', [e(1) e(2)], ...
    't_max_hollow', 0.75 - t_min / 200, 'knots_from', [], 'sets', [], 'set_ranges', zeros(0, 2), ...
    'inner_fn', @(t, knots_from, zr) picky_box(geo, t, t_min, t0, knots_from));
start = struct('mode', 'modular_precast', 'edges', e, 'vs', 0, 't', [t0; NaN], ...
    'z_ballast', -2.0, 'solid_modules', 2);
[~, s2] = sti_stage2(config, 0, [1500; NaN]);
P = struct('k', 1, 'hollow', zeros(1, 0), 't_min', t_min, 't_max', [0.75 - t_min / 200; NaN], ...
    'vs_bounds', [-0.5 2.5], 'stage2', s2, 'config', config, 'pct', 10, 'tol_eq', 1e-6, ...
    'hs_fn', @(v) mwecmass.solid.hydrostatics_at_draft(geo, v, struct()));

[sol, ctx] = mwecmass.realise.modular_precast.solve(ctx, start, P);
steps = {sol.solver.step};
fprintf('steps run: %s; escalation %s; closest %s\n', strjoin(steps, ', '), sol.escalation, sol.closest);
for k = 1:numel(sol.notes)
    fprintf('note: %s\n', sol.notes{k});
end
for s = sol.solver
    fprintf('step %-11s exitflag %g, iterations %d, objective %.6g, max equality violation %.3g\n', ...
        s.step, s.exitflag, s.iterations, s.fval, s.max_eq_violation);
end
fprintf('stored design: z_ballast %.9f m, t %s mm, vs %.6f m\n', sol.design.z_ballast, ...
    mat2str(1000 * sol.design.t', 9), sol.design.vs);
order = {'fixed_draft', 'draft_free'};
check(isequal(steps, order(1:numel(steps))) && strcmp(sol.escalation, steps{end}), ...
    'steps recorded in escalation order, spill skipped');
check(any(~cellfun(@isempty, strfind(sol.notes, 'spill not possible'))), 'skipped spill noted');
check(any(~cellfun(@isempty, strfind(sol.notes, 'failed evaluations'))), 'failed evaluations noted');
t_used = sol.design.t(isfinite(sol.design.t));
check(all(t_used == t0), 'the stored shell has an adaptive fit');
adaptive = ctx.sets(~[ctx.sets.refit]);
check(all([adaptive.t] == t0), 'only buildable adaptive sets kept');
end

function inner = picky_box(geo, t, t_min, t0, knots_from)
if isempty(knots_from) && t ~= t0
    error('mwecmass:solid:FitNotConverged', 'picky_box: no adaptive fit at t = %.17g in this test', t);
end
inner = sti_inner_box(geo, t, t_min);
inner.refit = ~isempty(knots_from);
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
    error('test_precast_solve_failed_evals:fail', '%s', msg);
end
end
