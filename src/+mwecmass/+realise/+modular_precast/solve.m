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
%
%   Each step solves in two phases with fmincon (sqp): phase 1 minimises (GM / GM2 - 1)^2 with
%   flotation as its equality, so it reaches both equalities when it can and otherwise ends at the
%   design that holds flotation with GM closest to GM2; phase 2 runs from a phase-1 point that
%   meets both equalities and minimises the objective with both. The step's point is the best
%   design it evaluated, ranked: both equalities held (smallest objective), then flotation held
%   (smallest |GM residual|), then neither (smallest largest residual); ties by the objective.
%   Flotation ranks first because it must hold in every reported state (AGENTS section 3 item 32).
%   Escalation (AGENTS section 3 items 6, 20, 31; OD6 option ii):
%     fixed_draft  vs = Stage-2 value, ballast within k*;
%     spill        only when the two equalities are not both met: the same ballast variable, its
%                  upper bound widened to the top of k*+1 (k*+1 must exist and be hollow);
%     draft_free   only when flotation (mass balance) is still not met: vs is a variable too.
%   The escalation ends when the F10 check passes (accepted) or both equalities hold (the step's
%   optimum is the closest fail: it minimises the deviation, and a failed check does not release
%   the draft). Otherwise the closest fail is the best design over every step, by the same rank
%   (the smallest equality violation, with flotation held where any design held it). The solver
%   exposes its evaluations portably, not its iterates, so every evaluated point (iterates, their
%   finite-difference neighbours within the bounds, the adaptive rebuilds) is a candidate; the
%   stored design is rebuilt on its adaptive fit.
%   A kernel error named in contract F5 and section 8 (VoidClosed, JointNotNested,
%   FitNotConverged), and a non-finite point proposed by the solver, is a failed evaluation:
%   objective and residuals Inf, which the SQP line search rejects; when one leaves the solver
%   without a gradient (fmincon stops with an error), the step keeps its best point. Failed
%   evaluations never stop Stage 3; an error of any other kind raised while a design is
%   evaluated is a defect and is raised.
%   Inner surfaces keep their knot vectors during a solve (F2 opts.knots_from, the adaptive set of
%   the step's start point); the step's point is rebuilt on its adaptive fit, and when that fit has
%   another knot structure the step restarts there (contract section 7 item 4) until the
%   structure no longer changes or returns to one already searched (noted).
%
%   sol: design (S3 of the stored design, built on adaptive sets), escalation (last step run),
%   solver (struct array: step, exitflag (phase 2, or phase 1 when phase 2 did not run; NaN when
%   fmincon stopped with an error), iterations, fval, max_eq_violation of the step's stored point),
%   closest ('' accepted | 'optimum' | 'least_violation'), notes (cellstr).

store = containers.Map();
store('ctx') = ctx;
store('hs') = struct('vs', {}, 'hs', {});
store('hist') = struct('step', {}, 'x', {}, 'free_vs', {}, 'res', {}, 'f', {});
store('defect') = [];
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

% No step met both equalities: the best design evaluated over every step.
hist = store('hist');
order = zeros(0, 1);
if ~isempty(hist)
    keys = cell2mat(arrayfun(@(h) rank_key(h.res, h.f, P.tol_eq), hist(:), 'UniformOutput', false));
    [~, order] = sortrows(keys);
    order = order(isfinite(keys(order, 2)));
end
for j = order'
    h = hist(j);
    Pj = P;
    Pj.free_vs = h.free_vs;
    Pj.step = h.step;
    rj = evaluate(store, Pj, h.x, true);
    if ~rj.failed
        sol.design = rj.design;
        sol.closest = 'least_violation';
        sol.notes{end + 1} = sprintf(['closest fail: no step met both equalities; stored the best ' ...
            'design evaluated (step %s): flotation residual %.3g, GM residual %.3g on its own ' ...
            'adaptive fit (%.3g, %.3g during the solve)'], h.step, rj.res(1), rj.res(2), h.res(1), ...
            h.res(2));
        ctx = store('ctx');
        return
    end
end
sol.design = start;
sol.closest = 'least_violation';
sol.notes{end + 1} = 'closest fail: no evaluated design could be rebuilt on its adaptive fit; the split design is stored';
ctx = store('ctx');
end

function key = rank_key(res, f, tol)
% [class, measure, objective]: class 0 both equalities held (measure 0), 1 flotation held
% (measure |GM residual|), 2 neither (measure the largest residual).
res = res(:)';
if all(abs(res) <= tol)
    key = [0, 0, f];
elseif abs(res(1)) <= tol
    key = [1, abs(res(2)), f];
else
    key = [2, max(abs(res)), f];
end
end

function [x, r, step, notes] = run_step(store, P, name, x0, z_hi, free_vs)
% One escalation step: the two phases on fixed knot vectors, restarted while the adaptive fit at
% the step's point changes the knot structure.
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
% The solver works on q = (x - lb) ./ w in [0, 1]: metres of ballast level and of shell
% thickness then weigh alike in its relative step and stopping tests. x is clamped to the bounds
% so that rounding never puts a thickness past t_max.
w = ub - lb;
to_x = @(q) min(max(lb + w .* q, lb), ub);
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
flag = NaN;
x = x0;
while true
    ctx = store('ctx');
    ctx.knots_from = base;
    store('ctx') = ctx;
    store('seen') = containers.Map();
    gm_sq = @(q) pick(cached(store, P, to_x(q)), 3)^2;
    flot = @(q) deal([], pick(cached(store, P, to_x(q)), 2));
    obj = @(q) pick(cached(store, P, to_x(q)), 1);
    both = @(q) deal([], pick(cached(store, P, to_x(q)), [2; 3]));
    q01 = {zeros(size(x)), ones(size(x))};
    try
        v = cached(store, P, x);
        if any(abs(v(2:3)) > P.tol_eq)
            [qs, ~, flag, out] = fmincon(gm_sq, (x - lb) ./ w, [], [], [], [], q01{:}, flot, opts);
            x = to_x(qs);
            iterations = iterations + out.iterations;
            v = cached(store, P, x);
        end
        if all(abs(v(2:3)) <= P.tol_eq)
            [qs, ~, flag, out] = fmincon(obj, (x - lb) ./ w, [], [], [], [], q01{:}, both, opts);
            x = to_x(qs);
            iterations = iterations + out.iterations;
        end
    catch err
        if ~isempty(store('defect'))
            rethrow(store('defect'));
        end
        flag = NaN;
        notes{end + 1} = sprintf('step %s: the solver stopped (%s)', name, err.message); %#ok<AGROW>
    end
    x = best_of_step(store, P, name, x, x0);
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
        notes{end + 1} = sprintf(['step %s: the adaptive fit at the step''s point returns to a knot ' ...
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
    notes{end + 1} = sprintf('step %s: the step''s point could not be rebuilt on its adaptive fit (%s)', ...
        name, r.message);
end
fprintf(['      step %s: exitflag %g, %d iterations, %d knot restarts, objective %.6g, ' ...
    'residuals flotation %.3g, GM %.3g\n'], name, flag, iterations, numel(searched) - 1, r.f, r.res);
end

function y = pick(v, idx)
y = v(idx);
end

function v = cached(store, P, x)
% [objective; flotation residual; GM residual] of x; the solver asks for the objective and the
% constraints at the same points, so each point is built once (the key prints every double
% exactly).
seen = store('seen');
key = sprintf('%.17g,', x);
if isKey(seen, key)
    v = seen(key);
    return
end
if ~all(isfinite(x))
    % a solver that lost its gradient can propose a non-finite point: a failed evaluation
    v = Inf(3, 1);
    store('failed') = [store('failed'), {'the solver proposed a non-finite point'}];
    return
end
try
    r = evaluate(store, P, x, false);
catch err
    store('defect') = err;
    rethrow(err);
end
v = [r.f; r.res(:)];
seen(key) = v;
end

function x = best_of_step(store, P, name, x, x0)
% The best design this step evaluated (rank_key); x0 when x is not finite and nothing ranks.
hist = store('hist');
if ~isempty(hist)
    hist = hist(strcmp({hist.step}, name) & arrayfun(@(h) numel(h.x) == numel(x), hist));
end
best = [];
if all(isfinite(x))
    v = cached(store, P, x);
    best = rank_key(v(2:3), v(1), P.tol_eq);
else
    x = x0;
end
for h = hist
    key = rank_key(h.res, h.f, P.tol_eq);
    if isempty(best) || is_less(key, best)
        x = h.x;
        best = key;
    end
end
end

function less = is_less(a, b)
d = find(a ~= b, 1);
less = ~isempty(d) && a(d) < b(d);
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
hist(end + 1) = struct('step', P.step, 'x', x, 'free_vs', P.free_vs, 'res', r.res, 'f', r.f);
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
