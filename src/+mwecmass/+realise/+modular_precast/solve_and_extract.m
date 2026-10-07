function [realised, final_props] = solve_and_extract(config, x_opt, final3d, opts)
%SOLVE_AND_EXTRACT  Modular-precast Stage 3 from the whole Stage-2 solution: split, build, check.
%
%   [realised, final_props] = mwecmass.realise.modular_precast.solve_and_extract(config, x_opt, final3d, opts)
%
%   x_opt = [vs; rho_1 ... rho_N] and final3d = results.Final3D are the Stage-2 solution. config:
%   hull_solid (S1b), ms2_model, boundary_cache, strip_edges, wall_strip_index,
%   constructability_rho_hull (UHPC density), constructability_rho_air, constructability_t_min,
%   per_strip_density_lb (Stage-2 floors, stored per module; optional), mass_acceptable_pct, and
%   the fields evaluate_realised reads. opts.inner_fn (optional): provider of an S2 inner set,
%   @(t, knots_from); default mwecmass.solid.offset_surface over the hollow range.
%
%   Process (AGENTS section 3 items 4, 7, 28, 32; contract section 7):
%   1. Module volumes V_i of the exact outer body and the split V_uhpc,i = V_i (rho_i - rho_air) /
%      (rho_uhpc - rho_air) (split_from_stage2). The wall module, modules at rho_uhpc and every
%      module below the ballast module k* are full solid sections.
%   2. Each hollow module above k* gets the shell t_i >= t_min that holds V_uhpc,i: the root of
%      V_uhpc,i(t) - V_uhpc,i on [t_min, t_max), t_max = d_close - eps_fit/2 over the hollow range
%      (F2b), bracketed by bisection towards t_max (never evaluated, the void closes there), then
%      fzero. Trial shells are refitted on fixed knots (F2 knots_from), first those of the adaptive
%      t_min set; the stored shell is the adaptive fit at the root. When that fit has another knot
%      structure the search restarts on it (contract section 7 item 4), per module, until the
%      structure no longer changes; a structure already searched ends the loop with a note. A
%      module whose t_min shell already holds more UHPC than V_uhpc,i is built at t_min; one whose
%      void closes first keeps the thickest shell evaluated. Both are reported.
%   3. k* has the shell t_k* = t_min above its ballast (OD13). The ballast level in k* is the root
%      of M(z_ballast) = rho_w V_sub at the Stage-2 draft (flotation, OD10). When no level inside
%      k* reaches it, the module end with the smaller residual is kept and flotation fails.
%   4. The body (F5, F6) and the hydrostatics (F7) at the Stage-2 draft give final_props (F9); F10
%      checks Z_CG, GM, coupled T_heave and T_pitch against Stage 2 with mass_acceptable_pct and
%      flotation against TOL_EQ. status is 'accepted' when F10 passes, else 'failed'; the realised
%      design is stored either way and Stage-2 properties are never returned.
%   realised: S8 (escalation 'split'); final_props = realised.props + stage3_status, stage3_check.

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
ctx = struct('geo', geo, 'rho', rho, 'z_hollow', z_hollow, 'knots_from', [], 'sets', []);
if isfield(opts, 'inner_fn') && ~isempty(opts.inner_fn)
    ctx.inner_fn = opts.inner_fn;
else
    ctx.inner_fn = @(t, knots_from) mwecmass.solid.offset_surface(config.ms2_model, ...
        config.boundary_cache, geo, t, z_hollow, struct('t_min', t_min, 'knots_from', knots_from));
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

notes = {};
t = NaN(N, 1);
t(k) = t_min;
z_ballast = e(1);
if ~isempty(k)
    z_ballast = e(k);
end
design = design_of(e, vs, t, z_ballast, find(split.solid)');
if any(split.hollow)
    [design.t, ctx, shell_notes] = shells_for_targets(ctx, design, split, config, geo, t_min, eps_fit);
    notes = [notes, shell_notes];
end

solver = struct('step', 'split', 'exitflag', NaN, 'iterations', 0, 'fval', NaN, ...
    'max_eq_violation', NaN);
if isempty(k)
    notes{end + 1} = 'every module is solid UHPC, so no ballast level can restore flotation';
else
    [design.z_ballast, ctx, solver, note] = ballast_for_flotation(ctx, design, k, M_target, solver);
    if ~isempty(note)
        notes{end + 1} = note;
    end
end

[ev, ctx] = mwecmass.realise.modular_precast.realise_modules(ctx, design); %#ok<ASGLU>
props = mwecmass.realise.evaluate_realised(ev.bp, hs, design, config);
check = mwecmass.realise.check_against_stage2(props, stage2, config.mass_acceptable_pct, TOL_EQ, ...
    config.RHO_WATER);
X = [props.CG_total(3), props.periods.heave, props.periods.pitch];
X2 = [stage2.Z_CG, stage2.T_heave, stage2.T_pitch];
solver.fval = sum(((X - X2) ./ X2).^2);
solver.max_eq_violation = max(abs([check.equalities.residual]));

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
    'reason', reason, 'escalation', 'split', 'vs', vs, 'draft', hs.draft, 'stage2', stage2, ...
    'rho', rho, 'design', design, 'k_star', k, 'V_uhpc_target', split.V_uhpc_target, ...
    'modules', {modules}, 'props', props, 'check', check, 'solver', solver, 'fit', {fit}, ...
    'body', ev.body, 'step_files', {struct('name', {}, 'path', {}, 'bodies', {})});
final_props = props;
final_props.stage3_status = status;
final_props.stage3_check = check;
print_realised(realised, split, notes);
end

function d = design_of(e, vs, t, z_ballast, solid_modules)
d = struct('mode', 'modular_precast', 'edges', e, 'vs', vs, 't', t, 'z_ballast', z_ballast, ...
    'solid_modules', solid_modules);
end

function [t, ctx, notes] = shells_for_targets(ctx, design, split, config, geo, t_min, eps_fit)
% Shell thickness of every hollow module above k* that holds its V_uhpc target.
notes = {};
t = design.t;
hollow = find(split.hollow)';
t_max = mwecmass.solid.void_closing_distance(config.ms2_model, config.boundary_cache, geo, ...
    ctx.z_hollow) - eps_fit / 2;
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
        [t(i), closed, ctx] = shell_root(ctx, design, hollow, i, target, t_min, t_max);
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
            'holds the Stage-2 split %.6g m^3; built at t = %.6g m'], i, t_max, target, t(i)); %#ok<AGROW>
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

function print_realised(r, split, notes)
fprintf('      realised modules (body frame):\n');
fprintf('      %-3s %9s %10s %12s %12s %12s %10s %10s\n', 'mod', 't[mm]', 'h_ball[m]', 'V_uhpc[m3]', ...
    'target[m3]', 'V_air[m3]', 'rho_eff', 'rho2');
for i = 1:numel(r.modules)
    m = r.modules(i);
    fprintf('      %-3d %9.3f %10.5f %12.6f %12.6f %12.6f %10.3f %10.3f\n', i, 1000 * m.t, m.h_ballast, ...
        m.V_uhpc, split.V_uhpc_target(i), m.V_air, m.rho_eff, m.rho_stage2);
end
fprintf('      z_ballast = %.6f m (body), mass %.3f kg, flotation residual %.3g\n', r.design.z_ballast, ...
    r.props.mass_total, r.check.equalities(1).residual);
fprintf('      %-8s %12s %12s %10s %8s  %s\n', 'metric', 'Stage 2', 'realised', 'dev[%]', 'limit', 'pass');
for m = r.check.metrics
    fprintf('      %-8s %12.6f %12.6f %10.4f %8.2f  %d\n', m.name, m.stage2, m.value, 100 * m.rel_dev, ...
        100 * m.limit, m.pass);
end
fprintf('      GM equality residual %.3g (Stage-3 solver constraint)\n', r.check.equalities(2).residual);
for k = 1:numel(notes)
    fprintf('      note: %s\n', notes{k});
end
fprintf('      Stage 3 status: %s%s\n', r.status, ...
    mwecmass.internal.ternary(isempty(r.reason), '', [' (' r.reason ')']));
end
