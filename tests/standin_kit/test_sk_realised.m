function test_sk_realised()
%TEST_SK_REALISED  sti_config, sti_stage2, sti_realised and F9, F10 on both fixtures and modes.
%   Reads no stand-in marker of F1-F7. F9 and F10 are the real src functions (T5), which shadow
%   their stand-ins; final_props adds stage3_status and stage3_check to the F9 fields.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
schema = mwecmass.output.export_schema();
fp_names = setdiff({schema.final_props_fields.name}, {'stage3_status', 'stage3_check'});
s8 = {'mode', 'hull_name', 'status', 'reason', 'escalation', 'vs', 'draft', 'stage2', 'rho', 'design', ...
    'k_star', 'V_uhpc_target', 'modules', 'props', 'check', 'solver', 'fit', 'body', 'step_files'};
% {fixture, mode, Stage-2 strip densities (NaN: solved for flotation), t, z_ballast, solid_modules, vs}
cases = {
    'cylinder', 'modular_precast', [NaN; 250; 200; 150], [0.1; 0.1; 0.1; 0.1], -2.5, [], 0.5
    'cylinder', 'thin_shell', [NaN; 400; 250; 200], 0.0254 * ones(4, 1), -2.85, [], 0.5
    'box', 'modular_precast', [NaN; 600; 500], [0.08; 0.08; 0.08], -2.3, [], 0
    'box', 'thin_shell', [NaN; 500; 300], 0.0254 * ones(3, 1), -2.3, [], 0.5
    };
for c = 1:size(cases, 1)
    [name, mode, rho2, t, zb, solid, vs] = cases{c, :};
    config = sti_config(name);
    for f = {'RHO_WATER', 'G', 'hydro_drafts', 'hydro_z_cg', 'added_mass_diagonal', 'added_mass_full', ...
            'radiation_damping_full', 'profile', 'strip_edges', 'mass_acceptable_pct'}
        check(isfield(config, f{1}), 'sti_config lacks %s', f{1});
    end
    config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
    [f3, stage2] = sti_stage2(config, vs, rho2);
    % one density solved for flotation: M = rho_w V_sub up to the rounding of an N-term dot product
    flo = f3.mass_total / (config.RHO_WATER * f3.V_sub) - 1;
    check(abs(flo) <= 8 * eps * numel(rho2), '%s: sti_stage2 does not float (%.3e)', name, flo);
    check(f3.GM_L > 0 && isfinite(stage2.T_pitch), '%s: Stage-2 design is not stable', name);
    check(stage2.Z_CG == f3.CG_total(3) && stage2.GM == f3.GM_L && stage2.T_pitch == f3.periods.pitch, 'S8.stage2');
    design = struct('mode', mode, 'edges', config.strip_edges, 'vs', vs, 't', t, 'z_ballast', zb, ...
        'solid_modules', solid);
    if strcmp(mode, 'modular_precast')
        rho = struct('uhpc', 2500, 'air', 1.2);
    else
        rho = struct('ballast', 7500, 'shell', 7850, 'air', 1.2);
    end
    [r, fp] = sti_realised(config, design, rho, stage2, struct('t_min', t(1)));
    check(isequal(sort(fieldnames(r))', sort(s8)), '%s %s: S8 fields', name, mode);
    missing = setdiff(fp_names, fieldnames(r.props));
    check(isempty(missing), 'F9 stand-in lacks %s', strjoin(missing, ', '));
    check(isfield(fp, 'stage3_status') && isfield(fp, 'stage3_check'), 'final_props stage3 fields');
    p = r.props;
    bp = mwecmass.solid.body_properties(r.body, rho, []);
    check(p.CG_total(3) == bp.total.CG_body(3) + vs, 'I6 CG frame');
    check(p.GM_L == p.KM - p.CG_total(3) && isequal(p.cross_section, config.profile), 'F9 GM and profile');
    check(isequal(p.densities_at_nodes, [r.modules.rho_eff]'), 'F9 densities');
    mods = r.modules;
    for i = 1:numel(mods)
        V = mods(i).V_air;
        if strcmp(mode, 'modular_precast')
            V = V + mods(i).V_uhpc;
        else
            V = V + mods(i).V_ballast + mods(i).V_shell;
        end
        % the region volumes are the summands of the module volume: equal to rounding of 3 terms
        check(abs(V - mods(i).V) <= 4 * eps * mods(i).V, 'module %d region sum', i);
        check(isequal(mods(i).CG_world, bp.modules(i).CG_body + [0 0 vs]), 'module %d CG_world', i);
    end
    if strcmp(mode, 'modular_precast')
        check(~isempty(r.k_star) && numel(r.V_uhpc_target) == numel(mods), 'precast split fields');
    end
    m = r.check.metrics;
    fprintf(['%-8s %-15s Stage 2: M %.1f kg, Z_CG %.4f m, GM %.4f m, T %.3f/%.3f s | realised: M %.1f kg, ' ...
        'Z_CG %.4f m, GM %.4f m, T %.3f/%.3f s | %s\n'], name, mode, stage2.mass, stage2.Z_CG, stage2.GM, ...
        stage2.T_heave, stage2.T_pitch, p.mass_total, p.CG_total(3), p.GM_L, p.periods.heave, ...
        p.periods.pitch, r.status);
    fprintf('         relative deviations Z_CG %.3g, GM %.3g, T_heave %.3g, T_pitch %.3g; flotation residual %.3g\n', ...
        m.rel_dev, r.check.equalities(1).residual);

    % F10 on a Stage 2 equal to the realised values passes the metrics and the GM equality
    same = struct('vs', vs, 'rho', rho2, 'mass', p.mass_total, 'Z_CG', p.CG_total(3), 'GM', p.GM_L, ...
        'T_heave', p.periods.heave, 'T_pitch', p.periods.pitch);
    ck = mwecmass.realise.check_against_stage2(p, same, 10, 1e-6, config.RHO_WATER);
    check(all([ck.metrics.pass]) && ck.equalities(2).pass && ck.equalities(2).residual == 0, 'F10 identity');
    check(ck.pass == ck.equalities(1).pass && ...
        ck.equalities(1).residual == p.mass_total / (config.RHO_WATER * p.V_sub) - 1, 'F10 flotation residual');
    off = same;
    off.Z_CG = same.Z_CG * 1.25;
    ck = mwecmass.realise.check_against_stage2(p, off, 10, Inf, config.RHO_WATER);
    check(~ck.pass && isequal(ck.failed, {'Z_CG'}) && ~isempty(strfind(ck.reason, 'Z_CG')), 'F10 failure report');

    % F9 on the geo of the real F1 (analytic = []) identifies the fixture by hull_name
    hs = mwecmass.solid.hydrostatics_at_draft(config.hull_solid, vs, struct());
    real_cfg = config;
    real_cfg.hull_solid.analytic = [];
    p2 = mwecmass.realise.evaluate_realised(bp, hs, design, real_cfg);
    check(isequaln(p2, p), '%s %s: F9 differs when geo.analytic is empty', name, mode);
end

% all modules solid: no inner set
config = sti_config('box');
config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
[~, stage2] = sti_stage2(config, 0, [NaN; 600; 500]);
design = struct('mode', 'modular_precast', 'edges', config.strip_edges, 'vs', 0, 't', NaN(3, 1), ...
    'z_ballast', -2.5, 'solid_modules', 1:3);
r = sti_realised(config, design, struct('uhpc', 2500, 'air', 1.2), stage2, struct('t_min', 0.08));
check(isempty(r.body.inner_t) && all([r.modules.V_air] == 0) && all(isnan([r.modules.t])), 'all-solid design');
fprintf('box precast, every module solid: M %.1f kg, Z_CG %.4f m, GM %.4f m | %s\n', ...
    r.props.mass_total, r.props.CG_total(3), r.props.GM_L, r.status);
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

function check(cond, varargin)
if ~cond
    error('test_sk_realised:fail', varargin{:});
end
end
