function [sol, ctx] = solve(ctx, start, P)
%SOLVE  Modular-precast Stage-3 optimisation from the split design: escalation and closest fail.
%
%   [sol, ctx] = mwecmass.realise.modular_precast.solve(ctx, start, P)
%
%   ctx: kernel context of solve_and_extract (realise_modules). start: S3 design of the split
%   (solve_and_extract), the start point. P: k (ballast module k*), hollow (indices of the hollow
%   modules above k*), t_min, t_max [N x 1] (t_max,i = d_close over module i - eps_fit/2, contract
%   F2b), vs_bounds [1 x 2] (world), stage2 (S8 stage2), config (read by evaluate_realised;
%   RHO_WATER), pct (mass_acceptable_pct), tol_eq (the fmincon ConstraintTolerance), hs_fn (@(vs) S7,
%   F7), hs_cache (optional): struct array (vs, hs) of F7 results already computed.
%
%   Variables x = [z_ballast; t_k*; t_i of P.hollow] and, in the last step, vs. Objective (AGENTS
%   section 3 item 27): sum over X = Z_CG (CG_total(3), world), coupled T_heave and coupled T_pitch
%   of ((X - X2) / X2)^2. Equalities: flotation M / (rho_w V_sub) - 1 = 0 and GM / GM2 - 1 = 0.
%   Bounds: t_min <= t <= t_max,i; z_ballast from the bottom of k* to the top of k* (or of k*+1).
%   Modules wholly below z_ballast are solid (t = NaN).
%
%   Each step (AGENTS section 3 items 27, 31, 32):
%   1. Flotation first, by a monotone 1-D root: the mass rises with z_ballast and with every t
%      (rho_uhpc > rho_air), the displaced mass with the draft. With the draft fixed: z_ballast
%      within its bounds; when no level floats the design, z_ballast at that bound and every t
%      moved together towards t_min (too heavy) or t_max,i (too light; bisection that never
%      evaluates t_max,i, where the void closes, down to the StepTolerance). With the draft free:
%      vs for the design's own mass within vs_bounds, then the same path at the nearer bound when
%      no draft floats it. No root: flotation is out of reach in the step (noted).
%   2. Both equalities by Newton steps: the minimum-norm step on the residuals over the variables
%      free of their bounds (forward differences, step sqrt(eps) as sqp's), halved until the
%      residual norm falls, until tol_eq or a step below the StepTolerance. When they do not reach
%      tol_eq, fmincon (sqp) phase 1 minimises (GM / GM2 - 1)^2 with flotation as its equality,
%      and Newton steps follow from the best point evaluated.
%   3. From the best point that holds both equalities, fmincon (sqp) phase 2 minimises the
%      objective with both. Newton steps bring its end point back onto the equalities, and phase 2
%      runs again from the best point while that objective falls, within MaxIterations sqp
%      iterations per step.
%   The step's point is the best design it evaluated that rebuilds on its adaptive fit, ranked: both
%   equalities held (smallest objective), then flotation held (smallest |GM residual|), then
%   neither (smallest largest residual); ties by the objective. Flotation ranks first because it
%   must hold in every reported state (item 32). Each step starts from the previous step's point
%   and evaluates it, so its point ranks no worse.
%   Escalation (AGENTS section 3 items 6, 20, 31; OD6 option ii):
%     fixed_draft  vs = Stage-2 value, ballast within k*;
%     spill        only when the two equalities are not both met: the same ballast variable, its
%                  upper bound widened to the top of k*+1 (k*+1 must exist and be hollow);
%     draft_free   only when the step's point does not hold flotation, so no design evaluated and
%                  built at the Stage-2 draft floats: vs is a variable too.
%   The escalation ends when the F10 check passes (accepted) or both equalities hold (the step's
%   optimum is the closest fail: it minimises the deviation, and a failed check does not release
%   the draft). Otherwise the last step's point is the closest fail (least violation). The stored
%   design is always the point of the last step run, so sol.solver(end) describes it.
%   A kernel error named in contract F5 and section 8 (VoidClosed, JointNotNested,
%   FitNotConverged), and a non-finite point proposed by the solver, is a failed evaluation:
%   objective and residuals Inf, which the SQP line search, the Newton steps and the ranking
%   reject; when fmincon stops with an error, the step goes on from its best point. Failed
%   evaluations never stop Stage 3; an error of any other kind raised while a design is
%   evaluated is a defect and is raised.
%   Inner surfaces keep their knot vectors during a solve (F2 opts.knots_from, the adaptive set of
%   the step's start point); the step's point is rebuilt on its adaptive fit, and when that fit has
%   another knot structure the step restarts there (contract section 7 item 4) until the
%   structure no longer changes or returns to one already searched (noted). F7 is called once per
%   distinct vs (contract section 7 item 3).
%
%   sol: design (S3 of the stored design, built on adaptive sets), escalation (last step run),
%   solver (struct array: step; exitflag of the last fmincon call, NaN when none ran or it stopped
%   with an error; iterations (sqp, Newton and root iterations); fval and max_eq_violation of the
%   step's point; fval_phase2_start, the objective where phase 2 first started, NaN when it did
%   not run), closest ('' accepted | 'optimum' | 'least_violation'), notes (cellstr), hs_cache
%   (P.hs_cache and every F7 result of the solve).

store = containers.Map();
store('ctx') = ctx;
hs_cache = struct('vs', {}, 'hs', {});
if isfield(P, 'hs_cache')
    hs_cache = P.hs_cache;
end
store('hs') = hs_cache;
store('hist') = struct('step', {}, 'round', {}, 'x', {}, 'free_vs', {}, 'res', {}, 'f', {});
store('defect') = [];
store('seen') = containers.Map();
store('round') = 0;

e = start.edges(:);
k = P.k;
N = numel(e) - 1;
P.start = start;
P.nh = numel(P.hollow);
P.max_iter = 300;
P.step_tol = 1e-8;
P.opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'iter', ...
    'ConstraintTolerance', P.tol_eq, 'OptimalityTolerance', 1e-6, 'StepTolerance', P.step_tol, ...
    'MaxIterations', P.max_iter, 'MaxFunctionEvaluations', 2000);
sol = struct('design', start, 'escalation', '', 'closest', '', 'notes', {{}}, ...
    'solver', struct('step', {}, 'exitflag', {}, 'iterations', {}, 'fval', {}, ...
    'max_eq_violation', {}, 'fval_phase2_start', {}), 'hs_cache', []);

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
        if r.check.equalities(1).pass
            break
        end
        x = [x; start.vs]; %#ok<AGROW>
    end
    [x, r, step, notes] = run_step(store, P, name, x, z_hi, free_vs);
    sol.solver(end + 1) = step;
    sol.notes = [sol.notes, notes];
    sol.escalation = name;
    sol.design = r.design;
    if r.check.pass
        [sol, ctx] = finish(sol, store);
        return
    end
    if all([r.check.equalities.pass])
        sol.closest = 'optimum';
        [sol, ctx] = finish(sol, store);
        return
    end
end
sol.closest = 'least_violation';
sol.notes{end + 1} = sprintf(['closest fail: no step met both equalities; stored the point of step ' ...
    '%s: flotation residual %.3g, GM residual %.3g'], sol.escalation, r.res(1), r.res(2));
[sol, ctx] = finish(sol, store);
end

function [sol, ctx] = finish(sol, store)
ctx = store('ctx');
sol.hs_cache = store('hs');
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

function c = rank_class(v, tol)
key = rank_key(v(2:3), v(1), tol);
c = key(1);
end

function [x, r, step, notes] = run_step(store, P, name, x0, z_hi, free_vs)
% One escalation step on fixed knot vectors, restarted while the adaptive fit at the step's point
% changes the knot structure.
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
% The solvers work on q = (x - lb) ./ w in [0, 1]: metres of ballast level and of shell
% thickness then weigh alike in their relative steps and stopping tests. x is clamped to the
% bounds so that rounding never puts a thickness past t_max.
S = struct('lb', lb, 'ub', ub, 'w', ub - lb);
fprintf('      Stage-3 step %s: %d variables, z_ballast in [%.6f, %.6f] m (body)%s\n', name, ...
    numel(x0), lb(1), ub(1), mwecmass.internal.ternary(free_vs, ', vs free', ', vs fixed'));

store('failed') = {};
store('notes') = {};
store('iters') = 0;
store('sqp_iters') = 0;
store('flag') = NaN;
store('f2') = NaN;
store('round') = 0;
r0 = evaluate(store, P, x0, true);
base = [];
if ~r0.failed && ~isempty(r0.inner)
    base = r0.inner(1);
end
searched = {base};
x = x0;
while true
    ctx = store('ctx');
    ctx.knots_from = base;
    store('ctx') = ctx;
    store('seen') = containers.Map();
    store('round') = store('round') + 1;
    x = phases(store, P, S, x);
    [x, r] = best_buildable(store, P, numel(x0));
    if isempty(base)
        break
    end
    changed = find(~arrayfun(@(s) mwecmass.realise.modular_precast.same_knots(s, base), r.inner), 1);
    if isempty(changed)
        break
    end
    adaptive = r.inner(changed);
    if any(cellfun(@(s) mwecmass.realise.modular_precast.same_knots(adaptive, s), searched))
        add_note(store, sprintf(['step %s: the adaptive fit at the step''s point returns to a knot ' ...
            'structure already searched; stored there'], name));
        break
    end
    searched{end + 1} = adaptive; %#ok<AGROW>
    base = adaptive;
end
ctx = store('ctx');
ctx.knots_from = [];
store('ctx') = ctx;
failed = store('failed');
if ~isempty(failed)
    add_note(store, sprintf('step %s: %d failed evaluations, counted as rejected points (first: %s)', ...
        name, numel(failed), failed{1}));
end
notes = store('notes');
step = struct('step', name, 'exitflag', store('flag'), 'iterations', store('iters'), 'fval', r.f, ...
    'max_eq_violation', r.viol, 'fval_phase2_start', store('f2'));
fprintf(['      step %s: exitflag %g, %d iterations, %d knot restarts, objective %.6g, ' ...
    'residuals flotation %.3g, GM %.3g\n'], name, step.exitflag, step.iterations, ...
    numel(searched) - 1, r.f, r.res);
end

function x = phases(store, P, S, x)
% Steps 1 to 3 of the doc on the knots of the current round.
tol = P.tol_eq;
v = cached(store, P, x);
if isfinite(v(2)) && abs(v(2)) > tol
    x = flotation_root(store, P, S, x);
end
x = best_of_step(store, P, x);
if rank_class(cached(store, P, x), tol) > 0
    x = restore(store, P, S, x, [2 3]);
end
if rank_class(cached(store, P, x), tol) > 0
    run_fmincon(store, P, S, 'gm', x);
    x = restore(store, P, S, best_of_step(store, P, x), [2 3]);
    if rank_class(cached(store, P, x), tol) == 2
        x = restore(store, P, S, x, 2);
    end
    x = best_of_step(store, P, x);
end
while store('sqp_iters') < P.max_iter
    v = cached(store, P, x);
    if rank_class(v, tol) > 0
        break
    end
    if isnan(store('f2'))
        store('f2') = v(1);
    end
    x2 = run_fmincon(store, P, S, 'objective', x);
    if rank_class(cached(store, P, x2), tol) > 0
        restore(store, P, S, x2, [2 3]);
    end
    xn = best_of_step(store, P, x);
    if ~(pick(cached(store, P, xn), 1) < v(1))
        break
    end
    x = xn;
end
x = best_of_step(store, P, x);
end

function x = flotation_root(store, P, S, x)
% Monotone 1-D roots of the flotation residual (step 1 of the doc).
if P.free_vs
    v = cached(store, P, x);
    if ~isfinite(v(4))
        return
    end
    h = @(vs) P.config.RHO_WATER * pick(hydrostatics(store, P, vs), 'V_sub') - v(4);
    ends = P.vs_bounds;
    he = [h(ends(1)), h(ends(2))];
    if he(1) * he(2) <= 0
        x(end) = root(store, h, ends);
        cached(store, P, x);
        return
    end
    [~, j] = min(abs(he));
    x(end) = ends(j);
    add_note(store, sprintf(['step %s: no draft within the vs bounds [%.6g, %.6g] m floats the ' ...
        'design (%.6g kg); vs set to %.6g m, then the design path'], P.step, ends, v(4), ends(j)));
end
ti = 2:numel(x) - P.free_vs;
v = cached(store, P, x);
if ~isfinite(v(2))
    return
end
heavy = v(2) > 0;
xe = x;
if heavy
    xe(1) = S.lb(1);
else
    xe(1) = S.ub(1);
end
ve = cached(store, P, xe);
if (heavy && ve(2) <= 0) || (~heavy && ve(2) >= 0)
    path = @(p) [p; x(2:end)];
    x = path(root(store, @(p) pick(cached(store, P, path(p)), 2), sort([x(1), xe(1)])));
    cached(store, P, x);
    return
end
xt = xe;
if heavy
    xt(ti) = S.lb(ti);
    vt = cached(store, P, xt);
    if ~(vt(2) <= 0)
        add_note(store, sprintf(['step %s: flotation out of reach: the lightest design (every t = ' ...
            't_min, z_ballast %.6g m) has flotation residual %.3g'], P.step, xt(1), vt(2)));
        return
    end
    bracket = [0, 1];
else
    % bisection towards t_max,i, which is never evaluated: the void closes there
    xt(ti) = S.ub(ti);
    bracket = [0, 0.5];
    while true
        vt = cached(store, P, xe + bracket(2) * (xt - xe));
        if vt(2) >= 0
            break
        end
        bracket = [bracket(2), (1 + bracket(2)) / 2];
        if 1 - bracket(2) < P.step_tol
            add_note(store, sprintf(['step %s: flotation out of reach: with z_ballast %.6g m and ' ...
                'every t towards t_max,i the flotation residual stays at or below %.3g'], P.step, ...
                xe(1), vt(2)));
            return
        end
    end
end
path = @(p) xe + p * (xt - xe);
x = path(root(store, @(p) pick(cached(store, P, path(p)), 2), bracket));
cached(store, P, x);
end

function p = root(store, g, ends)
% fzero on a sign-changing bracket; when it stops with an error (a failed evaluation inside the
% bracket) the end with the smaller residual is kept.
try
    [p, ~, ~, out] = fzero(g, ends, optimset('Display', 'off'));
    store('iters') = store('iters') + out.iterations;
catch err
    if ~isempty(store('defect'))
        rethrow(store('defect'));
    end
    [~, j] = min(abs([g(ends(1)), g(ends(2))]));
    p = ends(j);
    add_note(store, sprintf('the flotation root stopped (%s)', err.message));
end
end

function x = restore(store, P, S, x, idx)
% Newton steps on the residuals v(idx) (2 flotation, 3 GM) with the bounds held.
tol = P.tol_eq;
n = numel(x);
to_x = @(q) min(max(S.lb + S.w .* q, S.lb), S.ub);
q = (x - S.lb) ./ S.w;
v = cached(store, P, x);
g = v(idx);
for it = 1:P.max_iter
    if all(abs(g) <= tol) || ~all(isfinite(g))
        break
    end
    J = zeros(numel(idx), n);
    for j = 1:n
        h = sqrt(eps);
        if q(j) + h > 1
            h = -h;
        end
        qj = q;
        qj(j) = q(j) + h;
        vj = cached(store, P, to_x(qj));
        J(:, j) = (vj(idx) - g) / h;
    end
    if ~all(isfinite(J(:)))
        break
    end
    free = true(n, 1);
    while true
        d = zeros(n, 1);
        d(free) = -pinv(J(:, free)) * g;
        out = free & ((q <= 0 & d < 0) | (q >= 1 & d > 0));
        if ~any(out)
            break
        end
        free(out) = false;
    end
    store('iters') = store('iters') + 1;
    a = 1;
    accepted = false;
    while a * max(abs(d)) >= P.step_tol
        qn = min(max(q + a * d, 0), 1);
        vn = cached(store, P, to_x(qn));
        if norm(vn(idx)) < norm(g)
            q = qn;
            g = vn(idx);
            x = to_x(qn);
            accepted = true;
            break
        end
        a = a / 2;
    end
    if ~accepted
        break
    end
end
end

function x = run_fmincon(store, P, S, kind, x)
% Phase 1 ('gm') or phase 2 ('objective') with fmincon on q; returns its end point (x when it
% stops with an error or returns a non-finite point).
to_x = @(q) min(max(S.lb + S.w .* q, S.lb), S.ub);
if strcmp(kind, 'gm')
    fun = @(q) pick(cached(store, P, to_x(q)), 3)^2;
    con = @(q) deal([], pick(cached(store, P, to_x(q)), 2));
else
    fun = @(q) pick(cached(store, P, to_x(q)), 1);
    con = @(q) deal([], pick(cached(store, P, to_x(q)), [2; 3]));
end
n = numel(x);
try
    [qs, ~, flag, out] = fmincon(fun, (x - S.lb) ./ S.w, [], [], [], [], zeros(n, 1), ones(n, 1), ...
        con, P.opts);
    store('iters') = store('iters') + out.iterations;
    if strcmp(kind, 'objective')
        store('sqp_iters') = store('sqp_iters') + max(out.iterations, 1);
    end
    if all(isfinite(qs))
        x = to_x(qs);
    end
catch err
    if ~isempty(store('defect'))
        rethrow(store('defect'));
    end
    flag = NaN;
    if strcmp(kind, 'objective')
        store('sqp_iters') = P.max_iter;
    end
    add_note(store, sprintf('step %s: fmincon (%s) stopped (%s)', P.step, kind, err.message));
end
store('flag') = flag;
end

function y = pick(v, idx)
if ischar(idx)
    y = v.(idx);
else
    y = v(idx);
end
end

function add_note(store, note)
store('notes') = [store('notes'), {note}];
end

function v = cached(store, P, x)
% [objective; flotation residual; GM residual; mass] of x; the solvers ask for the objective and
% the constraints at the same points, so each point is built once per knot round (the key prints
% every double exactly).
seen = store('seen');
key = sprintf('%.17g,', x);
if isKey(seen, key)
    v = seen(key);
    return
end
if ~all(isfinite(x))
    % a solver that lost its gradient can propose a non-finite point: a failed evaluation
    v = Inf(4, 1);
    store('failed') = [store('failed'), {'the solver proposed a non-finite point'}];
    return
end
try
    r = evaluate(store, P, x, false);
catch err
    store('defect') = err;
    rethrow(err);
end
v = [r.f; r.res(:); r.M];
seen(key) = v;
end

function x = best_of_step(store, P, x)
% The best design of this step and knot round with as many variables as x (rank_key); x when
% none was built.
n = numel(x);
hist = store('hist');
hist = hist(strcmp({hist.step}, P.step) & [hist.round] == store('round') & ...
    arrayfun(@(h) numel(h.x) == n, hist));
best = [];
for h = hist
    key = rank_key(h.res, h.f, P.tol_eq);
    if isempty(best) || is_less(key, best)
        x = h.x;
        best = key;
    end
end
end

function [x, r] = best_buildable(store, P, n)
% The best design of this step (every knot round) that rebuilds on its adaptive fit.
hist = store('hist');
hist = hist(strcmp({hist.step}, P.step) & arrayfun(@(h) numel(h.x) == n, hist));
keys = cell2mat(arrayfun(@(h) rank_key(h.res, h.f, P.tol_eq), hist(:), 'UniformOutput', false));
[~, order] = sortrows(keys);
tried = containers.Map();
for j = order'
    x = hist(j).x;
    key = sprintf('%.17g,', x);
    if isKey(tried, key)
        continue
    end
    tried(key) = true;
    r = evaluate(store, P, x, true);
    if ~r.failed
        return
    end
end
error('mwecmass:realise:NoBuildablePoint', ['step %s: no evaluated design rebuilds on its adaptive ' ...
    'fit, not even the step''s start, which the previous step built'], P.step);
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
    'viol', Inf, 'M', Inf, 'check', [], 'inner', []);
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
r.M = props.mass_total;
r.check = check;
r.inner = ev.inner;
hist = store('hist');
hist(end + 1) = struct('step', P.step, 'round', store('round'), 'x', x, 'free_vs', P.free_vs, ...
    'res', r.res, 'f', r.f);
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
