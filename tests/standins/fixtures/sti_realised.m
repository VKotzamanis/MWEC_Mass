function [realised, final_props] = sti_realised(config, design, rho, stage2, opts)
%STI_REALISED  S8 realised design (results.stage3) of an S3 design on a stand-in fixture.
%
%   [realised, final_props] = sti_realised(config, design, rho, stage2, opts)
%
%   config: sti_config, with config.hull_solid (S1b; mwecmass.solid.outer_nurbs of
%   config.ms2_model when absent). design: S3. rho: densities by region (uhpc, air | ballast,
%   shell, air) [kg/m^3]. stage2: the S8.stage2 struct (second output of sti_stage2).
%   opts.t_min [m] (required: eps_fit = 0.01 t_min), opts.escalation (default 'split'),
%   opts.tol_eq (default 1e-6, the Stage-3 fmincon ConstraintTolerance), opts.rho_uhpc and
%   opts.rho_floor [N x 1] (default NaN).
%   The design is evaluated with the public kernel names: inner sets by
%   mwecmass.solid.offset_surface (cylinder) or sti_inner_box (box), then build_body,
%   body_properties, hydrostatics_at_draft, mwecmass.realise.evaluate_realised and
%   check_against_stage2. The F9 stand-in identifies a fixture also by hull_name, so it accepts the
%   geo of the real F1 (analytic = []) between J1 and J2. V_uhpc_target is the precast split of AGENTS section 3 item 4.1,
%   V_i (rho_i - rho_air) / (rho_uhpc - rho_air). No solver runs: solver holds one entry for
%   opts.escalation with exitflag NaN. final_props = props + stage3_status, stage3_check.

if nargin < 5
    opts = struct();
end
escalation = get_opt(opts, 'escalation', 'split');
tol_eq = get_opt(opts, 'tol_eq', 1e-6);
if ~isfield(config, 'hull_solid') || isempty(config.hull_solid)
    config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
end
geo = config.hull_solid;
fx = sti_closed_form('fixture', geo);
e = design.edges(:);
N = numel(e) - 1;
precast = strcmp(design.mode, 'modular_precast');
rho_floor = get_opt(opts, 'rho_floor', NaN(N, 1));

hollow = setdiff(1:N, design.solid_modules);
ts = unique(design.t(hollow));
ts = ts(isfinite(ts));
z_top = e(end);
if ~isempty(design.solid_modules)
    z_top = e(min(design.solid_modules));
end
inner = [];
for k = 1:numel(ts)
    if strcmp(fx.kind, 'box')
        set = sti_inner_box(geo, ts(k), opts.t_min);
    else
        set = mwecmass.solid.offset_surface(config.ms2_model, config.boundary_cache, geo, ts(k), ...
            [e(1) z_top], struct('t_min', opts.t_min));
    end
    inner = [inner, set]; %#ok<AGROW>
end
body = mwecmass.solid.build_body(geo, design, inner);
bp = mwecmass.solid.body_properties(body, rho, struct());
hs = mwecmass.solid.hydrostatics_at_draft(geo, design.vs, struct());
props = mwecmass.realise.evaluate_realised(bp, hs, design, config);
check = mwecmass.realise.check_against_stage2(props, stage2, config.mass_acceptable_pct, tol_eq, config.RHO_WATER);

names = fieldnames(bp.regions)';
zb = design.z_ballast;
modules = struct('z_lo', num2cell(e(1:end - 1)), 'z_hi', num2cell(e(2:end)));
for i = 1:N
    m = bp.modules(i);
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
    modules(i).rho_stage2 = stage2.rho(i);
    modules(i).rho_floor = rho_floor(i);
    modules(i).CG_world = m.CG_body + [0 0 design.vs];
end
k_star = [];
V_uhpc_target = [];
if precast
    k_star = find(e(1:end - 1) < zb & zb <= e(2:end), 1);
    rho_uhpc = get_opt(opts, 'rho_uhpc', rho.uhpc);
    V_uhpc_target = [modules.V]' .* (stage2.rho(:) - rho.air) / (rho_uhpc - rho.air);
end
status = 'accepted';
if ~check.pass
    status = 'failed';
end
fit = [];
if ~isempty(inner)
    fit = [inner.report];
end
realised = struct('mode', design.mode, 'hull_name', geo.hull_name, 'status', status, ...
    'reason', check.reason, 'escalation', escalation, 'vs', design.vs, 'draft', hs.draft, ...
    'stage2', stage2, 'rho', rho, 'design', design, 'k_star', k_star, 'V_uhpc_target', V_uhpc_target, ...
    'modules', {modules}, 'props', props, 'check', check, ...
    'solver', struct('step', escalation, 'exitflag', NaN, 'iterations', 0, 'fval', NaN, ...
    'max_eq_violation', max(abs([check.equalities.residual]))), ...
    'fit', {fit}, 'body', body, 'step_files', {struct('name', {}, 'path', {}, 'bodies', {})});
final_props = props;
final_props.stage3_status = status;
final_props.stage3_check = check;
end

function v = get_opt(opts, name, default)
if isfield(opts, name) && ~isempty(opts.(name))
    v = opts.(name);
else
    v = default;
end
end
