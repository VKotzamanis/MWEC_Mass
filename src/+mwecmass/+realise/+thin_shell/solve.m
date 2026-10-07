function realised = solve(config, x_opt, final3d, fids)
%SOLVE  Thin-shell Stage 3: one shell thickness t and z_ballast that realise the Stage-2 design.
%
%   realised = mwecmass.realise.thin_shell.solve(config, x_opt, final3d, fids)
%
%   config: build_config output with hull_solid (S1b), ms2_model, boundary_cache, strip_edges,
%   rho_shell, rho_ballast, rho_air, steel_t_min, steel_t_init, vertical_shift_bounds,
%   mass_acceptable_pct, RHO_WATER, G and the hydrodynamic cache fields. x_opt = [vs; rho_1..N]
%   and final3d (results.Final3D) are the Stage-2 solution. fids: report destinations (default 1).
%   realised: S8 (results.stage3) of the realised design, or of the closest design when Stage 3
%   fails; never the Stage-2 properties.
%
%   Body (contract S3, thin shell): below z_ballast the full outer section is steel ballast
%   (rho_ballast); above it a steel shell of one normal thickness t (rho_shell) for every module
%   encloses air (rho_air). z_ballast may lie in any module. Variables t in [t_min, t_max] with
%   t_max = d_close - eps_fit/2 (F2b: at t_max the void closes), z_ballast in [z_min, z_max], and
%   the draft last. Equalities, held to tol_eq (the fmincon ConstraintTolerance): flotation
%   M = rho_w V_sub and GM = GM_Stage2. Objective: sum over Z_CG, coupled T_heave and coupled
%   T_pitch of ((X3 - X2)/X2)^2 (AGENTS section 3 item 27).
%
%   Escalation (AGENTS section 3 items 20, 31): (t, z_ballast) at the Stage-2 draft; the draft
%   is released only when mass balance cannot be met there, i.e. when rho_w V_sub at the Stage-2
%   draft lies outside the mass range [M(t_min, z_min), M(t_min, z_max)] (with
%   rho_air < rho_shell <= rho_ballast the mass rises with t and with z_ballast, and at z_max it
%   no longer depends on t). At the Stage-2 draft the two equalities fix t and z_ballast, so the
%   step is a solve (solve_fixed_draft), not an optimisation. With the draft released, a design
%   meeting both equalities is searched first and fmincon (SQP) then minimises the objective
%   under both equalities (solve_draft_free). The mass_acceptable_pct check on Z_CG, GM, T_heave
%   and T_pitch then decides accepted or failed. Closest fail: when the equalities hold at the
%   step's result it is the design; otherwise the evaluated design with the smallest equality
%   violation (violation_rank): designs that hold flotation within tol_eq first, ranked by the GM
%   residual; then the others ranked by the flotation residual and, at equal flotation residual,
%   by the GM residual. The draft is released to restore mass balance, which therefore comes
%   first (AGENTS section 3 items 31, 32). solver (S8) describes the step's own result; when the
%   closest design replaces it, a report line names the replacement.
%
%   Within a step the inner set is refitted on fixed knot vectors; the reported design uses the
%   adaptive fit at its own t, and a changed knot structure restarts the step there (contract
%   section 7 item 4). solver.iterations counts the design evaluations of the step.

if nargin < 4 || isempty(fids)
    fids = 1;
end
t_start = tic;
required = {'hull_solid', 'ms2_model', 'boundary_cache', 'strip_edges', 'rho_shell', 'rho_ballast', ...
    'rho_air', 'steel_t_min', 'steel_t_init', 'vertical_shift_bounds', 'mass_acceptable_pct', ...
    'RHO_WATER', 'G'};
for k = 1:numel(required)
    if ~isfield(config, required{k}) || isempty(config.(required{k}))
        error('mwecmass:realise:MissingConfig', 'thin_shell.solve: config.%s is missing or empty.', required{k});
    end
end
geo = config.hull_solid;
edges = config.strip_edges(:);
N = numel(edges) - 1;
x_opt = x_opt(:);
if numel(x_opt) ~= N + 1
    error('mwecmass:realise:BadOptVector', ...
        'thin_shell.solve: x_opt must be [vs; rho_1..rho_%d] (got %d entries).', N, numel(x_opt));
end
rho = struct('ballast', config.rho_ballast, 'shell', config.rho_shell, 'air', config.rho_air);
if ~(rho.air < rho.shell && rho.shell <= rho.ballast)
    error('mwecmass:realise:DensityOrder', ...
        'thin_shell.solve needs rho_air < rho_shell <= rho_ballast (got %g, %g, %g kg/m^3).', ...
        rho.air, rho.shell, rho.ballast);
end
t_min = config.steel_t_min;
tol_eq = 1e-6;
z_min = geo.z_range(1);
z_max = geo.z_range(2);
stage2 = struct('vs', x_opt(1), 'rho', x_opt(2:end), 'mass', final3d.mass_total, ...
    'Z_CG', final3d.CG_total(3), 'GM', final3d.GM_L, 'T_heave', final3d.periods.heave, ...
    'T_pitch', final3d.periods.pitch);

d_close = mwecmass.solid.void_closing_distance(config.ms2_model, config.boundary_cache, geo, geo.z_range);
eps_fit = 0.01 * t_min;
t_max = d_close - eps_fit / 2;
if t_min >= t_max
    error('mwecmass:realise:TMinExceedsTMax', ...
        'thin_shell.solve: t_min = %.5f m >= t_max = %.5f m (void closing distance %.5f m).', ...
        t_min, t_max, d_close);
end

ctx = struct('config', config, 'geo', geo, 'edges', edges, 'rho', rho, 'stage2', stage2, ...
    't_min', t_min, 't_max', t_max, 'tol_eq', tol_eq, ...
    'sets', containers.Map('KeyType', 'char', 'ValueType', 'any'), ...
    'adaptive_sets', containers.Map('KeyType', 'char', 'ValueType', 'any'), ...
    'hydro', containers.Map('KeyType', 'char', 'ValueType', 'any'), ...
    'state', containers.Map('KeyType', 'char', 'ValueType', 'any'));
ctx.state('history') = zeros(0, 6);
ctx.state('points') = containers.Map('KeyType', 'char', 'ValueType', 'any');

emit = @(varargin) mwecmass.output.emit(fids, varargin{:});
emit('\n    mwecmass.realise.thin_shell.solve (exact geometry):\n');
emit('      Hull z = [%.4f, %.4f] m (body); %d modules; edges %s m\n', z_min, z_max, N, mat2str(edges', 6));
emit('      t in [%.5f, %.5f] m: t_min (input), t_max = d_close - eps_fit/2, d_close = %.5f m, eps_fit = %.3g m\n', ...
    t_min, t_max, d_close, eps_fit);
emit('      Stage 2: vs %.4f m, M %.1f kg, Z_CG %.4f m, GM %.4f m, T_heave %.4f s, T_pitch %.4f s (coupled)\n', ...
    stage2.vs, stage2.mass, stage2.Z_CG, stage2.GM, stage2.T_heave, stage2.T_pitch);

t0 = min(max(config.steel_t_init, t_min), t_max);
ev0 = mwecmass.realise.thin_shell.evaluate_design_point(ctx, stage2.vs, t0, z_min, true);
ctx.state('ref') = ev0.inner;
ctx.sets(num2hex(t0)) = ev0.inner;

% mass range at the Stage-2 draft
ev_lo = point(ctx, [stage2.vs, t_min, z_min]);
ev_hi = point(ctx, [stage2.vs, t_min, z_max]);
reachable = ev_lo.ceq(1) <= tol_eq && ev_hi.ceq(1) >= -tol_eq;
emit('      Mass at the Stage-2 draft: target rho_w V_sub = %.1f kg; buildable [%.1f, %.1f] kg -> %s\n', ...
    config.RHO_WATER * ev_lo.V_sub, ev_lo.mass, ev_hi.mass, ...
    ternary(reachable, 'reachable, draft fixed', 'unreachable, draft released'));

if reachable
    step = 'fixed_draft';
    solve_once = @(p) solve_fixed_draft(ctx, stage2.vs);
else
    step = 'draft_free';
    solve_once = @(p) solve_draft_free(ctx, t0, p);
end
[ev, solver] = run_step(ctx, step, solve_once, emit);

if any(~(abs(ev.ceq) <= tol_eq))
    p = closest(ctx.state('history'), tol_eq);
    if ~isempty(p)
        ev_best = mwecmass.realise.thin_shell.evaluate_design_point(ctx, p(1), p(2), p(3), true);
        if ranks_before(violation_rank(ev_best.ceq, tol_eq), violation_rank(ev.ceq, tol_eq))
            ev = ev_best;
            emit(['      Closest evaluated design replaces the step''s result: vs %.6f m, t %.6f m, ' ...
                'z_ballast %.6f m; flotation residual %.6g, GM residual %.6g\n'], ev.design.vs, ...
                ev.design.t(1), ev.design.z_ballast, ev.ceq(1), ev.ceq(2));
        end
    end
end

realised = package(ctx, ev, step, solver, config);
emit_report(emit, realised, toc(t_start));
end

% ---------------------------------------------------------------- escalation steps

function [ev, info] = run_step(ctx, step, solve_once, emit)
% Runs one escalation step, then rebuilds its result with the adaptive inner set at the final t;
% a changed knot structure restarts the step on the new knots.
seen = {};
n0 = size(ctx.state('history'), 1);
p = [];
while true
    [p, exitflag] = solve_once(p);
    ev = mwecmass.realise.thin_shell.evaluate_design_point(ctx, p(1), p(2), p(3), true);
    ref = ctx.state('ref');
    if same_knots(ev.inner, ref) || any(cellfun(@(s) same_knots(ev.inner, s), seen))
        break
    end
    emit('      Knot structure of the inner set changed at t = %.6f m: step restarted there\n', p(2));
    seen{end + 1} = ref; %#ok<AGROW>
    ctx.state('ref') = ev.inner;
    remove(ctx.sets, keys(ctx.sets));
    remove(ctx.state('points'), keys(ctx.state('points')));
    ctx.sets(num2hex(p(2))) = ev.inner;
end
% solver fields of S8: the step's own result, also when the closest design replaces it
info = struct('step', step, 'exitflag', exitflag, 'iterations', size(ctx.state('history'), 1) - n0, ...
    'fval', ev.objective, 'max_eq_violation', max(abs(ev.ceq)));
emit('      Step %s: exitflag %d, %d design evaluations; t %.6f m, z_ballast %.6f m, vs %.6f m\n', ...
    step, exitflag, info.iterations, p(2), p(3), p(1));
end

function [p, exitflag] = solve_fixed_draft(ctx, vs)
% At a fixed draft the two equalities fix t and z_ballast (AGENTS section 5): along the
% flotation curve z_ballast(t) the mass moves from the ballast top up into the thicker shell, so
% Z_CG rises and the GM residual changes monotonically with t. Root in t of the GM residual on
% the curve, from t_min to the t where the curve reaches z_min; without a sign change the end
% with the smaller residual is the closest design (flotation held). exitflag 1: root, 0: end.
t_min = ctx.t_min;
t_c = curve_end(ctx, vs);
g = @(t) gm_on_curve(ctx, vs, t);
g_lo = g(t_min);
g_hi = g(t_c);
exitflag = 0;
if g_lo == 0
    t = t_min;
    exitflag = 1;
elseif g_hi == 0
    t = t_c;
    exitflag = 1;
elseif sign(g_lo) ~= sign(g_hi)
    t = fzero(g, [t_min, t_c]);
    exitflag = 1;
elseif abs(g_lo) <= abs(g_hi)
    t = t_min;
else
    t = t_c;
end
p = [vs, t, ballast_for_flotation(ctx, vs, t)];
end

function [p, exitflag] = solve_draft_free(ctx, t0, p_prev)
% Draft released: first a design that floats and meets GM = GM_Stage2 at t0, searched along the
% drafts where flotation can be met at t0 (rho_w V_sub falls as vs rises, so these drafts form
% one interval); then fmincon over (vs, t, z_ballast) with the objective and both equalities
% from there. exitflag: fmincon's, or that of the first search (1: root, 0: closest end) when
% fmincon cannot start because the objective is not finite there.
zr = ctx.geo.z_range;
cfg = ctx.config;
b = cfg.vertical_shift_bounds;
if isempty(p_prev)
    M_lo = mass_of(ctx, [ctx.stage2.vs, t0, zr(1)]);
    M_hi = mass_of(ctx, [ctx.stage2.vs, t0, zr(2)]);
    f_hi = @(v) cfg.RHO_WATER * displaced_volume(ctx, v) - M_hi;
    f_lo = @(v) cfg.RHO_WATER * displaced_volume(ctx, v) - M_lo;
    % smallest vs (deepest draft) at which rho_w V_sub <= M_hi; largest vs (shallowest draft) at
    % which rho_w V_sub >= M_lo
    if f_hi(b(1)) <= 0
        vs_a = b(1);
    elseif f_hi(b(2)) > 0
        vs_a = Inf;
    else
        vs_a = fzero(f_hi, b);
    end
    if f_lo(b(2)) >= 0
        vs_b = b(2);
    elseif f_lo(b(1)) < 0
        vs_b = -Inf;
    else
        vs_b = fzero(f_lo, b);
    end
    exitflag = 0;
    if vs_a > vs_b
        % no draft within the bounds floats a design at t0. The corners nearest to flotation: the
        % deepest draft with the lightest design (t_min, z_min) when the hull is too heavy
        % everywhere, the shallowest draft with the heaviest design (z_max, where the mass no
        % longer depends on t) when it is too light everywhere
        corners = [b(1), ctx.t_min, zr(1); b(2), t0, zr(2)];
        r = [flotation_of(ctx, corners(1, :)), flotation_of(ctx, corners(2, :))];
        [~, k] = min(abs(r));
        p = corners(k, :);
        return
    end
    G = @(vs) gm_of(ctx, [vs, t0, ballast_for_flotation(ctx, vs, t0)]);
    G_a = G(vs_a);
    G_b = G(vs_b);
    if sign(G_a) ~= sign(G_b)
        vs = fzero(G, [vs_a, vs_b]);
        exitflag = 1;
    elseif abs(G_a) <= abs(G_b)
        vs = vs_a;
    else
        vs = vs_b;
    end
    p = [vs, t0, ballast_for_flotation(ctx, vs, t0)];
else
    p = p_prev;
    exitflag = 0;
end
if ~isfinite(objective_of(ctx, p))
    return
end
% SQP starts from an identity Hessian, so the variables enter in units of their own size: t in
% t_min, z_ballast and vs in hull heights.
h = zr(2) - zr(1);
scale = [h; ctx.t_min; h];
map = @(x) scale' .* x(:)';
lb = [b(1); ctx.t_min; zr(1)];
ub = [b(2); ctx.t_max; zr(2)];
opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', 'StepTolerance', 1e-8, ...
    'OptimalityTolerance', 1e-6, 'ConstraintTolerance', ctx.tol_eq, 'MaxIterations', 200, ...
    'MaxFunctionEvaluations', 1000);
[x, ~, exitflag] = fmincon(@(x) objective_of(ctx, map(x)), p(:) ./ scale, [], [], [], [], ...
    lb ./ scale, ub ./ scale, @(x) constraints(ctx, map(x)), opts);
p = map(x);
end

function [c, ceq] = constraints(ctx, p)
c = [];
ev = point(ctx, p);
ceq = ev.ceq;
end

function f = objective_of(ctx, p)
ev = point(ctx, p);
f = ev.objective;
end

function r = flotation_of(ctx, p)
ev = point(ctx, p);
r = ev.ceq(1);
end

function r = gm_of(ctx, p)
ev = point(ctx, p);
r = ev.ceq(2);
end

function m = mass_of(ctx, p)
ev = point(ctx, p);
m = ev.mass;
end

function V = displaced_volume(ctx, vs)
key = num2hex(vs);
if ~isKey(ctx.hydro, key)
    ctx.hydro(key) = mwecmass.solid.hydrostatics_at_draft(ctx.geo, vs, struct());
end
hs = ctx.hydro(key);
V = hs.V_sub;
end

function t_c = curve_end(ctx, vs)
% Largest t of the flotation curve at this draft: where flotation needs no ballast. The mass
% rises with t; its upper bracket is searched by halving the distance to t_max (where the void
% closes and the kernel cannot evaluate), so no step size is chosen. When even the last
% representable t below t_max floats with no ballast, that t ends the curve.
r = @(t) flotation_of(ctx, [vs, t, ctx.geo.z_range(1)]);
t = ctx.t_min;
if r(t) >= 0
    t_c = t;
    return
end
while true
    t_next = (t + ctx.t_max) / 2;
    if t_next <= t || t_next >= ctx.t_max
        t_c = t;
        return
    end
    r_next = r(t_next);
    if ~isfinite(r_next)
        t_c = t;
        return
    end
    if r_next >= 0
        t_c = fzero(r, [t, t_next]);
        return
    end
    t = t_next;
end
end

function g = gm_on_curve(ctx, vs, t)
g = gm_of(ctx, [vs, t, ballast_for_flotation(ctx, vs, t)]);
end

function zb = ballast_for_flotation(ctx, vs, t)
% z_ballast that floats (vs, t): the mass rises with z_ballast, so a sign change brackets one
% root; without one, the end nearer to flotation.
z = ctx.geo.z_range;
r = @(zz) flotation_of(ctx, [vs, t, zz]);
r_lo = r(z(1));
r_hi = r(z(2));
if r_lo >= 0
    zb = z(1);
elseif r_hi <= 0
    zb = z(2);
else
    zb = fzero(r, z);
end
end

function ev = point(ctx, p)
% One kernel evaluation per design point, shared by every search; every evaluated design is
% recorded for the closest-fail choice. A kernel error (e.g. VoidClosed) marks the point as not
% realisable.
points = ctx.state('points');
key = [num2hex(p(1)), num2hex(p(2)), num2hex(p(3))];
if isKey(points, key)
    ev = points(key);
    return
end
try
    ev_full = mwecmass.realise.thin_shell.evaluate_design_point(ctx, p(1), p(2), p(3), false);
    ev = struct('objective', ev_full.objective, 'ceq', ev_full.ceq, 'mass', ev_full.props.mass_total, ...
        'V_sub', ev_full.props.V_sub);
catch err
    if ~strncmp(err.identifier, 'mwecmass:solid:', numel('mwecmass:solid:'))
        rethrow(err);
    end
    ev = struct('objective', NaN, 'ceq', [NaN; NaN], 'mass', NaN, 'V_sub', NaN);
end
points(key) = ev; %#ok<NASGU> containers.Map is a handle
ctx.state('history') = [ctx.state('history'); p, ev.objective, abs(ev.ceq')];
end

function p = closest(hist, tol_eq)
% The evaluated design with the smallest equality violation (violation_rank); ties: the smaller
% objective.
p = [];
ok = all(isfinite(hist(:, 5:6)), 2);
hist = hist(ok, :);
if isempty(hist)
    return
end
rank = cell2mat(arrayfun(@(k) violation_rank(hist(k, 5:6)', tol_eq), (1:size(hist, 1))', ...
    'UniformOutput', false));
[~, order] = sortrows([rank, hist(:, 4)]);
p = hist(order(1), 1:3);
end

function r = violation_rank(ceq, tol_eq)
% [tier, first, second], ordered lexicographically: designs holding flotation within tol_eq
% (tier 0, ranked by |GM residual|) come before the others (tier 1, ranked by |flotation
% residual| first and |GM residual| second). Mass balance comes first because the draft is
% released only to restore it (AGENTS section 3 items 31, 32).
if ~all(isfinite(ceq))
    r = [Inf, Inf, Inf];
elseif abs(ceq(1)) <= tol_eq
    r = [0, abs(ceq(2)), 0];
else
    r = [1, abs(ceq(1)), abs(ceq(2))];
end
end

function yes = ranks_before(a, b)
k = find(a ~= b, 1);
yes = ~isempty(k) && a(k) < b(k);
end

function same = same_knots(a, b)
same = numel(a.patches) == numel(b.patches);
k = 0;
while same && k < numel(a.patches)
    k = k + 1;
    same = isequal(a.patches(k).surf.knots, b.patches(k).surf.knots) && ...
        isequal(a.patches(k).surf.degree, b.patches(k).surf.degree);
end
end

% ---------------------------------------------------------------- S8

function realised = package(ctx, ev, step, solver, config)
e = ctx.edges;
N = numel(e) - 1;
design = ev.design;
zb = design.z_ballast;
rho_floor = NaN(N, 1);
if isfield(config, 'per_strip_density_lb') && numel(config.per_strip_density_lb) == N
    rho_floor = config.per_strip_density_lb(:);
end
names = {'ballast', 'shell', 'air'};
modules = struct('z_lo', num2cell(e(1:end - 1)), 'z_hi', num2cell(e(2:end)));
for i = 1:N
    m = ev.bp.modules(i);
    modules(i).t = NaN;
    if m.V_air > 0
        modules(i).t = design.t(i);
    end
    modules(i).h_ballast = min(max(zb - e(i), 0), e(i + 1) - e(i));
    modules(i).V = m.V;
    for r = 1:numel(names)
        modules(i).(['V_' names{r}]) = m.(['V_' names{r}]);
    end
    modules(i).mass = m.mass;
    modules(i).rho_eff = m.rho_eff;
    modules(i).rho_stage2 = ctx.stage2.rho(i);
    modules(i).rho_floor = rho_floor(i);
    modules(i).CG_world = m.CG_body + [0 0 design.vs];
end
status = 'failed';
if ev.check.pass
    status = 'accepted';
end
realised = struct('mode', 'thin_shell', 'hull_name', ctx.geo.hull_name, 'status', status, ...
    'reason', ev.check.reason, 'escalation', step, 'vs', design.vs, 'draft', ev.hs.draft, ...
    'stage2', ctx.stage2, 'rho', ctx.rho, 'design', design, 'k_star', [], 'V_uhpc_target', [], ...
    'modules', {modules}, 'props', ev.props, 'check', ev.check, 'solver', solver, ...
    'fit', ev.inner.report, 'body', ev.body, ...
    'step_files', {struct('name', {}, 'path', {}, 'bodies', {})});
end

function emit_report(emit, r, elapsed)
emit('      ---------------- Thin-shell Stage 3: %s (%s) ----------------\n', upper(r.status), r.escalation);
emit('      t %.5f m (%.3f in), z_ballast %.4f m (body), vs %.4f m, draft %.4f m\n', ...
    r.design.t(1), r.design.t(1) / 0.0254, r.design.z_ballast, r.vs, r.draft);
emit('      %-9s %12s %12s %10s %8s %s\n', 'metric', 'Stage 3', 'Stage 2', 'deviation', 'limit', 'pass');
for m = r.check.metrics
    emit('      %-9s %12.5f %12.5f %9.3f%% %7.1f%% %d\n', m.name, m.value, m.stage2, 100 * m.rel_dev, ...
        100 * m.limit, m.pass);
end
for q = r.check.equalities
    emit('      equality %-9s residual %10.3e  tol %.1e  pass %d\n', q.name, q.residual, q.tol, q.pass);
end
emit('      %-6s %9s %9s %10s %10s %10s %10s %11s %11s %11s\n', 'module', 'z_lo', 'z_hi', 'V_ballast', ...
    'V_shell', 'V_air', 'mass', 'rho_eff', 'rho_floor', 'rho_stage2');
for i = 1:numel(r.modules)
    m = r.modules(i);
    emit('      %-6d %9.4f %9.4f %10.4f %10.4f %10.4f %10.1f %11.1f %11.1f %11.1f\n', i, m.z_lo, m.z_hi, ...
        m.V_ballast, m.V_shell, m.V_air, m.mass, m.rho_eff, m.rho_floor, m.rho_stage2);
end
emit('      M %.1f kg (Stage 2 %.1f kg); elapsed %.1f s\n', r.props.mass_total, r.stage2.mass, elapsed);
if ~isempty(r.reason)
    emit('      reason: %s\n', r.reason);
end
end

function s = ternary(cond, a, b)
if cond
    s = a;
else
    s = b;
end
end
