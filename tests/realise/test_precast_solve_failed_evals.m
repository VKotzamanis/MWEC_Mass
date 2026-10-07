function test_precast_solve_failed_evals()
%TEST_PRECAST_SOLVE_FAILED_EVALS  Kernel errors inside the Stage-3 optimisation are failed
%evaluations, never a stop of Stage 3 (contract F5 and section 8).
%   Box fixture on the SK stand-ins; solve is called directly on a start design. The inner-set
%   provider raises mwecmass:solid:JointNotNested for every t but the two of the start design, so
%   every trial point that moves a shell fails while ballast moves still evaluate. Asserted: solve
%   returns, every step it runs is recorded, failed evaluations are noted, and the stored design
%   uses only thicknesses that could be built. Numbers are printed.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
config = sti_config('box');
config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
geo = config.hull_solid;
e = config.strip_edges;
t_min = 0.0762;
allowed = [t_min, 0.1];
ctx = struct('geo', geo, 'rho', struct('uhpc', 2500, 'air', 1.2), 'z_hollow', [e(1) e(end)], ...
    't_max_hollow', 0.75 - t_min / 200, 'knots_from', [], 'sets', [], 'set_ranges', zeros(0, 2), ...
    'inner_fn', @(t, knots_from, zr) picky_box(geo, t, t_min, allowed));
start = struct('mode', 'modular_precast', 'edges', e, 'vs', 0, 't', [t_min; 0.1; 0.1], ...
    'z_ballast', -2.0, 'solid_modules', zeros(1, 0));
[~, s2] = sti_stage2(config, 0, [NaN; 600; 600]);
P = struct('k', 1, 'hollow', [2 3], 't_min', t_min, 't_max', (0.75 - t_min / 200) * ones(3, 1), ...
    'vs_bounds', [-0.5 2.5], 'stage2', s2, 'config', config, 'pct', 10, 'tol_eq', 1e-6, ...
    'hs_fn', @(v) mwecmass.solid.hydrostatics_at_draft(geo, v, struct()));

[sol, ctx] = mwecmass.realise.modular_precast.solve(ctx, start, P);
steps = {sol.solver.step};
fprintf('steps run: %s; escalation %s; closest %s\n', strjoin(steps, ', '), sol.escalation, sol.closest);
for k = 1:numel(sol.notes)
    fprintf('note: %s\n', sol.notes{k});
end
fprintf('stored design: z_ballast %.9f m, t %s mm, vs %.6f m\n', sol.design.z_ballast, ...
    mat2str(1000 * sol.design.t', 9), sol.design.vs);
order = {'fixed_draft', 'spill', 'draft_free'};
check(isequal(steps, order(1:numel(steps))) && strcmp(sol.escalation, steps{end}), ...
    'steps recorded in escalation order');
check(any(~cellfun(@isempty, strfind(sol.notes, 'failed evaluations'))), 'failed evaluations noted');
t_used = sol.design.t(isfinite(sol.design.t));
check(all(ismember(t_used, allowed)), 'stored shells were buildable');
check(all(ismember([ctx.sets.t], allowed)), 'only buildable sets kept');
end

function inner = picky_box(geo, t, t_min, allowed)
if ~any(t == allowed)
    error('mwecmass:solid:JointNotNested', 'picky_box: t = %.17g is not buildable in this test', t);
end
inner = sti_inner_box(geo, t, t_min);
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
