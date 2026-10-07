function [sol, ctx] = solve(ctx, start, P)
%SOLVE  Modular-precast Stage-3 optimisation from the split design: escalation and closest fail.
%
%   [sol, ctx] = mwecmass.realise.modular_precast.solve(ctx, start, P)
%
%   ctx: kernel context of solve_and_extract (realise_modules). start: S3 design of the split
%   (solve_and_extract), the start point. P: k (ballast module k*), hollow (indices of the hollow
%   modules above k*), t_min, t_max [N x 1] (t_max,i = d_close over module i - eps_fit/2, contract
%   F2b), vs_bounds [1 x 2] (world), stage2 (S8 stage2), config (read by evaluate_realised),
%   pct (mass_acceptable_pct), tol_eq (the fmincon ConstraintTolerance), hs_fn (@(vs) S7, F7).
%
%   Variables x = [z_ballast; t_k*; t_i of P.hollow] and, in the last step, vs. Objective (AGENTS
%   section 3 item 27): sum over X = Z_CG (CG_total(3), world), coupled T_heave and coupled T_pitch
%   of ((X - X2) / X2)^2. Equalities: flotation M / (rho_w V_sub) - 1 = 0 and GM / GM2 - 1 = 0.
%   Bounds: t_min <= t <= t_max,i; z_ballast from the bottom of k* to the top of k* (or of k*+1).
%   Modules wholly below z_ballast are solid (t = NaN).
%   Escalation (AGENTS section 3 items 6, 20, 31; OD6 option ii):
%     fixed_draft  vs = Stage-2 value, ballast within k*;
%     spill        only when the equalities are not met: the same ballast variable, its upper bound
%                  widened to the top of k*+1 (k*+1 must exist and be hollow);
%     draft_free   only when flotation (mass balance) is still not met: vs is a variable too.
%   Every step ends when the F10 check passes (accepted) or the equalities hold (the optimum is the
%   closest fail: it minimises the deviation, and a failed check does not release the draft).
%   When no step meets the equalities, the closest fail is the evaluated design with the smallest
%   equality violation max(|flotation|, |GM|) over every step; the solver only exposes its
%   evaluations portably, so every point evaluated (iterates, their finite-difference neighbours
%   within the bounds, and the adaptive rebuilds of each step's start and optimum) is a candidate,
%   rebuilt on its adaptive fit before it is stored.
%   A kernel error named in contract F5 and section 8 (VoidClosed, JointNotNested,
%   FitNotConverged) is a failed evaluation: objective and equality residuals Inf, which the SQP
%   line search rejects; when one leaves the solver without a gradient, the step ends at its start
%   point. It never stops Stage 3.
%   Inner surfaces keep their knot vectors during a solve (F2 opts.knots_from, the adaptive sets
%   of the step's start point); the step's optimum is rebuilt on its adaptive fit, and when that
%   fit has another knot structure the step restarts there (contract section 7 item 4) until the
%   structure no longer changes or returns to one already searched (noted).
%
%   sol: design (S3 of the stored design, built on adaptive sets), escalation (last step run),
%   solver (struct array: step, exitflag, iterations, fval, max_eq_violation of the step's stored
%   point), closest ('' accepted | 'optimum' | 'least_violation'), notes (cellstr).

store = containers.Map();
store('ctx') = ctx;
store('hs') = struct('vs', {}, 'hs', {});
store('hist') = struct('step', {}, 'x', {}, 'free_vs', {}, 'viol', {});
store('seen') = containers.Map();
store('failed') = {};

e = start.edges(:);
k = P.k;
N = numel(e) - 1;
P.start = start;
P.nh = numel(P.hollow);
sol = struct('design', start, 'escalation', '', 'closest', '', 'notes', {{}}, ...
    'solver', struct('step', {}, 'exitflag', {}, 'iterations', {}, 'fval', {}, 'max_eq_violation', {}));

x = [start.z_ballast; start.t(k); start.t(P.hollow(:))];
spill_ok = k < N && any(P.hollow == k + 1);
plan = {'fixed_draft', e(k + 1), false};
if spill_ok
    plan(end + 1, :) = {'spill', e(k + 2), false};
end
plan(end + 1, :) = {'draft_free', plan{end, 2}, true};
if ~spill_ok
    sol.notes{end + 1} = sprintf(['spill not possible: module k*+1 = %d does not exist or is solid, ' ...
        'so the ballast stays inside module %d'], k + 1, k);
end

r = [];
for s = 1:size(plan, 1)
    [name, z_hi, free_vs] = plan{s, :};
    if strcmp(name, 'draft_free')
        if ~r.failed && r.check.equalities(1).pass
            break
        end
        x = [x; start.vs]; %#ok<AGROW>
    end
    [x, r, step, notes] = run_step(store, P, name, x, z_hi, free_vs);
    sol.solver(end + 1) = step;
    sol.notes = [sol.notes, notes];
    sol.escalation = name;
    sol.design = r.design;
    if ~r.failed && r.check.pass
        ctx = store('ctx');
        return
    end
    if ~r.failed && all([r.check.equalities.pass])
        sol.closest = 'optimum';
        ctx = store('ctx');
        return
    end
end

% No step met the equalities: the evaluated design with the smallest equality violation.
hist = store('hist');
viol = [hist.viol];
[~, order] = sort(viol);
for j = order(isfinite(viol(order)))
    h = hist(j);
    Pj = P;
    Pj.free_vs = h.free_vs;
    Pj.step = h.step;
    rj = evaluate(store, Pj, h.x, true);
    if ~rj.failed
        sol.design = rj.design;
        sol.closest = 'least_violation';
        sol.notes{end + 1} = sprintf(['closest fail: no step met the equalities; stored the design ' ...
            'evaluated in step %s with the smallest equality violation (%.3g on its own adaptive ' ...
            'fit, %.3g during the solve)'], h.step, rj.viol, h.viol);
        ctx = store('ctx');
        return
    end
end
sol.design = start;
sol.closest = 'least_violation';
sol.notes{end + 1} = 'closest fail: no evaluated design could be rebuilt on its adaptive fit; the split design is stored';
ctx = store('ctx');
end

function [x, r, step, notes] = run_step(store, P, name, x0, z_hi, free_vs)
% One escalation step: fmincon on fixed knot vectors, restarted while the adaptive fit at the
% optimum changes the knot structure.
e = P.start.edges(:);
k = P.k;
P.free_vs = free_vs;
P.step = name;
lb = [e(k); P.t_min; P.t_min * ones(P.nh, 1)];
ub = [z_hi; P.t_max(k); P.t_max(P.hollow(:))];
if free_vs
    lb(end + 1) = P.vs_bounds(1);
    ub(end + 1) = P.vs_bounds(2);
end
notes = {};
opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
    'ConstraintTolerance', P.tol_eq, 'OptimalityTolerance', 1e-6, 'StepTolerance', 1e-8, ...
    'MaxIterations', 300, 'MaxFunctionEvaluations', 2000);
fprintf('      Stage-3 step %s: %d variables, z_ballast in [%.6f, %.6f] m (body)%s\n', name, ...
    numel(x0), lb(1), ub(1), mwecmass.internal.ternary(free_vs, ', vs free', ', vs fixed'));

store('failed') = {};
r0 = evaluate(store, P, x0, true);
base = [];
if ~r0.failed && ~isempty(r0.inner)
    base = r0.inner(1);
end
searched = {base};
iterations = 0;
x = x0;
while true
    ctx = store('ctx');
    ctx.knots_from = base;
    store('ctx') = ctx;
    store('seen') = containers.Map();
    try
        [x, ~, flag, out] = fmincon(@(xx) objective(store, P, xx), x, [], [], [], [], lb, ub, ...
            @(xx) constraints(store, P, xx), opts);
        iterations = iterations + out.iterations;
    catch err
        % A failed evaluation inside a finite-difference gradient leaves the solver no direction;
        % the step ends at its last start point. Any other solver error is a defect.
        if isempty(store('failed'))
            rethrow(err);
        end
        flag = NaN;
        notes{end + 1} = sprintf('step %s: the solver stopped after failed evaluations (%s)', ...
            name, err.message); %#ok<AGROW>
    end
    ctx = store('ctx');
    ctx.knots_from = [];
    store('ctx') = ctx;
    r = evaluate(store, P, x, true);
    if r.failed || isempty(base)
        break
    end
    changed = find(~arrayfun(@(s) mwecmass.realise.modular_precast.same_knots(s, base), r.inner), 1);
    if isempty(changed)
        break
    end
    adaptive = r.inner(changed);
    if any(cellfun(@(s) mwecmass.realise.modular_precast.same_knots(adaptive, s), searched))
        notes{end + 1} = sprintf(['step %s: the adaptive fit at the optimum returns to a knot ' ...
            'structure already searched; stored there'], name); %#ok<AGROW>
        break
    end
    searched{end + 1} = adaptive; %#ok<AGROW>
    base = adaptive;
end
failed = store('failed');
if ~isempty(failed)
    notes{end + 1} = sprintf('step %s: %d failed evaluations, counted as rejected points (first: %s)', ...
        name, numel(failed), failed{1}); %#ok<AGROW>
end
step = struct('step', name, 'exitflag', flag, 'iterations', iterations, 'fval', r.f, ...
    'max_eq_violation', r.viol);
if r.failed
    notes{end + 1} = sprintf('step %s: the optimum could not be rebuilt on its adaptive fit (%s)', ...
        name, r.message);
end
fprintf('      step %s: exitflag %d, %d iterations, %d knot restarts, objective %.6g, max equality violation %.3g\n', ...
    name, flag, iterations, numel(searched) - 1, r.f, r.viol);
end

function f = objective(store, P, x)
v = cached(store, P, x);
f = v(1);
end

function [c, ceq] = constraints(store, P, x)
v = cached(store, P, x);
c = [];
ceq = v(2:end);
end

function v = cached(store, P, x)
% [objective; equality residuals] of x; the solver asks for the objective and the constraints at
% the same points, so each point is built once (the key prints every double exactly).
seen = store('seen');
key = sprintf('%.17g,', x);
if isKey(seen, key)
    v = seen(key);
    return
end
r = evaluate(store, P, x, false);
v = [r.f; r.res(:)];
seen(key) = v;
end

function r = evaluate(store, P, x, adaptive)
% Realised body, final_props and F10 residuals of the design x. adaptive: built on adaptive inner
% sets (knots_from empty), else on the knots of the solve.
design = design_of_x(P, x);
r = struct('design', design, 'failed', false, 'message', '', 'f', Inf, 'res', [Inf; Inf], ...
    'viol', Inf, 'check', [], 'inner', []);
ctx = store('ctx');
keep = ctx.knots_from;
if adaptive
    ctx.knots_from = [];
end
hs = hydrostatics(store, P, design.vs);
try
    [ev, ctx] = mwecmass.realise.modular_precast.realise_modules(ctx, design);
catch err
    if ~any(strcmp(err.identifier, {'mwecmass:solid:VoidClosed', 'mwecmass:solid:JointNotNested', ...
            'mwecmass:solid:FitNotConverged'}))
        rethrow(err);
    end
    r.failed = true;
    r.message = err.message;
    store('failed') = [store('failed'), {err.message}];
    return
end
ctx.knots_from = keep;
store('ctx') = ctx;
props = mwecmass.realise.evaluate_realised(ev.bp, hs, design, P.config);
check = mwecmass.realise.check_against_stage2(props, P.stage2, P.pct, P.tol_eq, P.config.RHO_WATER);
X = [props.CG_total(3), props.periods.heave, props.periods.pitch];
X2 = [P.stage2.Z_CG, P.stage2.T_heave, P.stage2.T_pitch];
r.f = sum(((X - X2) ./ X2).^2);
r.res = [check.equalities.residual]';
r.viol = max(abs(r.res));
r.check = check;
r.inner = ev.inner;
hist = store('hist');
hist(end + 1) = struct('step', P.step, 'x', x, 'free_vs', P.free_vs, 'viol', r.viol);
store('hist') = hist;
end

function d = design_of_x(P, x)
d = P.start;
k = P.k;
d.z_ballast = x(1);
d.t(k) = x(2);
d.t(P.hollow) = x(3:2 + P.nh);
if P.free_vs
    d.vs = x(end);
end
e = d.edges(:);
d.t(e(2:end) <= d.z_ballast) = NaN;
end

function hs = hydrostatics(store, P, vs)
% F7 is called once per draft (contract section 7 item 3).
list = store('hs');
for j = 1:numel(list)
    if isequal(list(j).vs, vs)
        hs = list(j).hs;
        return
    end
end
hs = P.hs_fn(vs);
list(end + 1) = struct('vs', vs, 'hs', hs);
store('hs') = list;
end
