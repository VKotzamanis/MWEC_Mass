function [realised, final_props] = solve_and_extract(config, x_opt, final3d, opts)
%SOLVE_AND_EXTRACT  Modular-precast Stage 3 from the whole Stage-2 solution: split, build, check,
%optimise, closest fail.
%
%   [realised, final_props] = mwecmass.realise.modular_precast.solve_and_extract(config, x_opt, final3d, opts)
%
%   x_opt = [vs; rho_1 ... rho_N] and final3d = results.Final3D are the Stage-2 solution. config:
%   hull_solid (S1b), ms2_model, boundary_cache, strip_edges, wall_strip_index,
%   constructability_rho_hull (UHPC density), constructability_rho_air, constructability_t_min,
%   per_strip_density_lb (Stage-2 floors, stored per module; optional), mass_acceptable_pct,
%   vertical_shift_bounds (Stage-2 bounds of vs, used when the draft is released), and the fields
%   evaluate_realised reads. opts.inner_fn (optional): provider of an S2 inner set,
%   @(t, knots_from, z_range); default mwecmass.solid.offset_surface. opts.d_close_fn (optional):
%   @(z_range) offset distance at which the void closes in z_range; default
%   mwecmass.solid.void_closing_distance (F2b). opts.escalate (optional, default true): false
%   stores the split design without the optimisation (tests of the split).
%
%   Process (AGENTS section 3 items 4 to 7, 20, 27, 28, 31, 32; contract sections 7 and 8):
%   1. Module volumes V_i of the exact outer body and the split V_uhpc,i = V_i (rho_i - rho_air) /
%      (rho_uhpc - rho_air) (split_from_stage2). The wall module, modules at rho_uhpc and every
%      module below the ballast module k* are full solid sections. Shell bounds per module:
%      t_max,i = d_close over the module's z range - eps_fit/2 (F2b; the void closes there).
%   2. Each hollow module above k* gets the shell t_i >= t_min that holds V_uhpc,i: the root of
%      V_uhpc,i(t) - V_uhpc,i on [t_min, t_max,i), bracketed by bisection towards t_max,i (never
%      evaluated, the void closes there), then fzero. Trial shells are refitted on fixed knots
%      (F2 knots_from), first those of the adaptive t_min set; the stored shell is the adaptive fit
%      at the root. When that fit has another knot structure the search restarts on it (contract
%      section 7 item 4), per module, until the structure no longer changes; a structure already
%      searched ends the loop with a note. A module whose t_min shell already holds more UHPC than
%      V_uhpc,i is built at t_min; one whose void closes first keeps the thickest shell evaluated.
%      Both are reported.
%   3. k* has the shell t_k* = t_min above its ballast (OD13). The ballast level in k* is the root
%      of M(z_ballast) = rho_w V_sub at the Stage-2 draft (flotation, OD10). When no level inside
%      k* reaches it, the module end with the smaller residual is kept and flotation fails.
%   4. The body (F5, F6) and the hydrostatics (F7, called once per draft and shared with solve)
%      at the Stage-2 draft give final_props (F9); F10 checks Z_CG, GM, coupled T_heave and
%      T_pitch against Stage 2 with mass_acceptable_pct and flotation against TOL_EQ (escalation
%      'split').
%   5. When that check fails, solve optimises from the split (fixed_draft, then spill, then
%      draft_free; see solve) and returns the accepted design or the closest fail, with the props,
%      check and body of the evaluation that gave its solver record.
%   status is 'accepted' when F10 passes on the stored design, else 'failed' with the failing
%   metrics, the escalation notes and the closest-fail rule in reason; the realised design is
%   stored either way and Stage-2 properties are never returned. stage3_report prints it.
%   realised: S8; final_props = realised.props + stage3_status, stage3_check.

% Stage-3 equality tolerance: the ConstraintTolerance of the Stage-3 fmincon (AGENTS OD10).
TOL_EQ = 1e-6;
if nargin < 4
    opts = struct();
end
geo = config.hull_solid;
e = config.strip_edges(:);
N = numel(e) - 1;
vs = x_opt(1);
rho_stage2 = x_opt(2:end);
rho_stage2 = rho_stage2(:);
rho = struct('uhpc', config.constructability_rho_hull, 'air', config.constructability_rho_air);
t_min = config.constructability_t_min;
eps_fit = 0.01 * t_min;
wall = config.wall_strip_index;
stage2 = struct('vs', vs, 'rho', rho_stage2, 'mass', final3d.mass_total, ...
    'Z_CG', final3d.CG_total(3), 'GM', final3d.GM_L, 'T_heave', final3d.periods.heave, ...
    'T_pitch', final3d.periods.pitch);
rho_floor = NaN(N, 1);
if isfield(config, 'per_strip_density_lb') && numel(config.per_strip_density_lb) == N
    rho_floor = config.per_strip_density_lb(:);
end

z_hollow = [e(1) e(end)];
if isequal(wall, N)
    z_hollow(2) = e(N);
elseif isequal(wall, 1)
    z_hollow(1) = e(2);
end
if isfield(opts, 'd_close_fn') && ~isempty(opts.d_close_fn)
    d_close_fn = opts.d_close_fn;
else
    d_close_fn = @(zr) mwecmass.solid.void_closing_distance(config.ms2_model, ...
        config.boundary_cache, geo, zr);
end
ctx = struct('geo', geo, 'rho', rho, 'z_hollow', z_hollow, ...
    't_max_hollow', d_close_fn(z_hollow) - eps_fit / 2, 'knots_from', [], 'sets', [], ...
    'set_ranges', zeros(0, 2));
if isfield(opts, 'inner_fn') && ~isempty(opts.inner_fn)
    ctx.inner_fn = opts.inner_fn;
else
    ctx.inner_fn = @(t, knots_from, zr) mwecmass.solid.offset_surface(config.ms2_model, ...
        config.boundary_cache, geo, t, zr, struct('t_min', t_min, 'knots_from', knots_from));
end

fprintf('\n    mwecmass.realise.modular_precast.solve_and_extract: split, build, check (Stage-2 draft)\n');
[ev, ctx] = mwecmass.realise.modular_precast.realise_modules(ctx, ...
    design_of(e, vs, NaN(N, 1), e(1), 1:N));
split = mwecmass.realise.modular_precast.split_from_stage2(rho_stage2, [ev.bp.modules.V]', ...
    rho.uhpc, rho.air, wall);
k = split.k_star;
hs = mwecmass.solid.hydrostatics_at_draft(geo, vs, struct());
M_target = config.RHO_WATER * hs.V_sub;
print_split(split, e, rho, t_min, M_target, stage2.mass);

t_max = NaN(N, 1);
for i = [k, find(split.hollow)']
    t_max(i) = d_close_fn([e(i), e(i + 1)]) - eps_fit / 2;
end

notes = {};
t = NaN(N, 1);
t(k) = t_min;
z_ballast = e(1);
if ~isempty(k)
    z_ballast = e(k);
end
design = design_of(e, vs, t, z_ballast, find(split.solid)');
if any(split.hollow)
    [design.t, ctx, shell_notes] = shells_for_targets(ctx, design, split, t_min, t_max);
    notes = [notes, shell_notes];
end

solver = struct('step', 'split', 'exitflag', NaN, 'iterations', 0, 'fval', NaN, ...
    'max_eq_violation', NaN, 'fval_phase2_start', NaN);
if isempty(k)
    notes{end + 1} = 'every module is solid UHPC, so no ballast level can restore flotation';
else
    [design.z_ballast, ctx, solver, note] = ballast_for_flotation(ctx, design, k, M_target, solver);
    if ~isempty(note)
        notes{end + 1} = note;
    end
end

[props, check, ev] = evaluate(ctx, design, hs, config, stage2, TOL_EQ);
X = [props.CG_total(3), props.periods.heave, props.periods.pitch];
X2 = [stage2.Z_CG, stage2.T_heave, stage2.T_pitch];
solver.fval = sum(((X - X2) ./ X2).^2);
solver.max_eq_violation = max(abs([check.equalities.residual]));
escalation = 'split';

if ~check.pass && ~isempty(k) && (~isfield(opts, 'escalate') || opts.escalate)
    vsb = [-geo.z_range(2), -geo.z_range(1)];
    if isfield(config, 'vertical_shift_bounds') && ~isempty(config.vertical_shift_bounds)
        vsb = [max(vsb(1), config.vertical_shift_bounds(1)), min(vsb(2), config.vertical_shift_bounds(2))];
    end
    P = struct('k', k, 'hollow', find(split.hollow)', 't_min', t_min, 't_max', t_max, ...
        'vs_bounds', vsb, 'stage2', stage2, 'config', config, 'pct', config.mass_acceptable_pct, ...
        'tol_eq', TOL_EQ, 'hs_fn', @(v) mwecmass.solid.hydrostatics_at_draft(geo, v, struct()), ...
        'hs_cache', struct('vs', vs, 'hs', hs));
    sol = mwecmass.realise.modular_precast.solve(ctx, design, P);
    design = sol.design;
    solver = [solver, sol.solver];
    escalation = sol.escalation;
    notes = [notes, sol.notes];
    props = sol.props;
    check = sol.check;
    ev = sol.ev;
    hs = sol.hs;
    if ~check.pass && strcmp(sol.closest, 'optimum')
        notes{end + 1} = sprintf(['closest fail: the %s result meets the equalities and fails the ' ...
            'mass_acceptable_pct check'], escalation);
    end
end

status = 'accepted';
reason = '';
if ~check.pass
    status = 'failed';
    reason = strjoin([{check.reason}, notes], '; ');
end
fit = [];
if ~isempty(ev.inner)
    fit = [ev.inner.report];
end
modules = mwecmass.realise.modular_precast.extract_strip_geometry(ev.bp, design, rho_stage2, rho_floor);
realised = struct('mode', 'modular_precast', 'hull_name', geo.hull_name, 'status', status, ...
    'reason', reason, 'escalation', escalation, 'vs', design.vs, 'draft', hs.draft, 'stage2', stage2, ...
    'rho', rho, 'design', design, 'k_star', k, 'V_uhpc_target', split.V_uhpc_target, ...
    'modules', {modules}, 'props', props, 'check', check, 'solver', solver, 'fit', {fit}, ...
    'body', ev.body, 'step_files', {struct('name', {}, 'path', {}, 'bodies', {})});
final_props = props;
final_props.stage3_status = status;
final_props.stage3_check = check;
mwecmass.realise.modular_precast.stage3_report(realised, notes);
end

function [props, check, ev] = evaluate(ctx, design, hs, config, stage2, tol_eq)
% The split design on its adaptive inner sets.
ctx.knots_from = [];
ev = mwecmass.realise.modular_precast.realise_modules(ctx, design);
props = mwecmass.realise.evaluate_realised(ev.bp, hs, design, config);
check = mwecmass.realise.check_against_stage2(props, stage2, config.mass_acceptable_pct, tol_eq, ...
    config.RHO_WATER);
end

function d = design_of(e, vs, t, z_ballast, solid_modules)
d = struct('mode', 'modular_precast', 'edges', e, 'vs', vs, 't', t, 'z_ballast', z_ballast, ...
    'solid_modules', solid_modules);
end

function [t, ctx, notes] = shells_for_targets(ctx, design, split, t_min, t_max)
% Shell thickness of every hollow module above k* that holds its V_uhpc target.
notes = {};
t = design.t;
hollow = find(split.hollow)';
probe = design;
probe.t(hollow) = t_min;
[ev, ctx] = mwecmass.realise.modular_precast.realise_modules(ctx, probe);
base = ev.inner(1);
vu_min = [ev.bp.modules.V_uhpc]';
for i = hollow
    target = split.V_uhpc_target(i);
    if vu_min(i) >= target
        t(i) = t_min;
        if vu_min(i) > target
            notes{end + 1} = sprintf(['module %d: the t_min shell holds %.6g m^3 UHPC, more than the ' ...
                'Stage-2 split %.6g m^3; built at t_min'], i, vu_min(i), target); %#ok<AGROW>
        end
        continue
    end
    searched = {base};
    restarts = 0;
    while true
        ctx.knots_from = searched{end};
        [t(i), closed, ctx] = shell_root(ctx, design, hollow, i, target, t_min, t_max(i));
        ctx.knots_from = [];
        probe.t(hollow) = t(i);
        [ev, ctx] = mwecmass.realise.modular_precast.realise_modules(ctx, probe);
        adaptive = ev.inner([ev.inner.t] == t(i));
        if mwecmass.realise.modular_precast.same_knots(adaptive, searched{end})
            break
        end
        if any(cellfun(@(s) mwecmass.realise.modular_precast.same_knots(adaptive, s), searched))
            notes{end + 1} = sprintf(['module %d: the adaptive fit at t = %.6g m returns to a knot ' ...
                'structure already searched; stored at that t'], i, t(i)); %#ok<AGROW>
            break
        end
        searched{end + 1} = adaptive; %#ok<AGROW>
        restarts = restarts + 1;
    end
    fprintf('      module %d: t = %.9f m, V_uhpc - target %.3g m^3, knot restarts %d\n', i, t(i), ...
        ev.bp.modules(i).V_uhpc - target, restarts);
    if closed
        notes{end + 1} = sprintf(['module %d: the void closes (t_max = %.6g m) before the shell ' ...
            'holds the Stage-2 split %.6g m^3; built at t = %.6g m'], i, t_max(i), target, t(i)); %#ok<AGROW>
    end
end
end

function [ti, closed, ctx] = shell_root(ctx, design, hollow, i, target, t_min, t_max)
% Root of V_uhpc,i(t) - target on [t_min, t_max) on the knots of ctx.knots_from; bisection towards
% t_max brackets it without evaluating t_max, where the void closes.
lo = t_min;
hi = t_max;
mid = lo + (hi - lo) / 2;
g_mid = -1;
while mid > lo && mid < hi
    [vu, ctx] = module_uhpc(ctx, design, hollow, mid, i);
    g_mid = vu - target;
    if g_mid >= 0
        break
    end
    lo = mid;
    mid = lo + (hi - lo) / 2;
end
closed = false;
if g_mid == 0
    ti = mid;
elseif g_mid > 0
    ti = fzero(@(tt) module_uhpc(ctx, design, hollow, tt, i) - target, [lo mid], ...
        optimset('Display', 'off'));
else
    ti = lo;
    closed = true;
end
end

function [vu, ctx] = module_uhpc(ctx, design, hollow, t, i)
d = design;
d.t(hollow) = t;
[ev, ctx] = mwecmass.realise.modular_precast.realise_modules(ctx, d);
vu = ev.bp.modules(i).V_uhpc;
end

function [z, ctx, solver, note] = ballast_for_flotation(ctx, design, k, M_target, solver)
% Ballast level in module k at which the realised mass equals the displaced mass.
note = '';
e = design.edges(:);
ends = [e(k), e(k + 1)];
f = zeros(1, 2);
for j = 1:2
    d = design;
    d.z_ballast = ends(j);
    [ev, ctx] = mwecmass.realise.modular_precast.realise_modules(ctx, d);
    f(j) = ev.bp.total.mass - M_target;
end
if f(1) <= 0 && f(2) >= 0
    if f(1) == 0
        z = ends(1);
    elseif f(2) == 0
        z = ends(2);
    else
        mass_gap = @(zz) total_mass(ctx, design, zz) - M_target;
        [z, ~, flag, out] = fzero(mass_gap, ends, optimset('Display', 'off'));
        solver.exitflag = flag;
        solver.iterations = out.iterations;
    end
    if isnan(solver.exitflag)
        solver.exitflag = 1;
    end
    return
end
[~, j] = min(abs(f));
z = ends(j);
side = {'bottom', 'top'};
note = sprintf(['flotation: with the ballast inside module %d the mass ranges over [%.6g, %.6g] kg, ' ...
    'displaced mass %.6g kg; ballast kept at the module %s'], k, f(1) + M_target, f(2) + M_target, ...
    M_target, side{j});
solver.exitflag = -2;
end

function M = total_mass(ctx, design, z)
d = design;
d.z_ballast = z;
ev = mwecmass.realise.modular_precast.realise_modules(ctx, d);
M = ev.bp.total.mass;
end

function print_split(split, e, rho, t_min, M_target, M_stage2)
fprintf('      rho_uhpc = %.1f kg/m^3, rho_air = %.2f kg/m^3, t_min = %.4f m\n', rho.uhpc, rho.air, t_min);
fprintf('      displaced mass at the Stage-2 draft %.3f kg (Stage-2 mass %.3f kg)\n', M_target, M_stage2);
fprintf('      %-3s %9s %9s %10s %10s %12s %12s  %s\n', 'mod', 'z_lo[m]', 'z_hi[m]', 'rho2', 'V[m3]', ...
    'V_uhpc[m3]', 'V_air[m3]', 'role');
for i = 1:numel(split.V)
    role = 'hollow';
    if isequal(split.wall, i)
        role = 'wall (solid)';
    elseif split.solid(i)
        role = 'solid';
    elseif isequal(split.k_star, i)
        role = 'ballast module k*';
    end
    fprintf('      %-3d %9.4f %9.4f %10.3f %10.5f %12.6f %12.6f  %s\n', i, e(i), e(i + 1), ...
        split.rho_stage2(i), split.V(i), split.V_uhpc_target(i), split.V_air_target(i), role);
end
end
