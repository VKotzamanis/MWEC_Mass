function test_evaluate_realised()
%TEST_EVALUATE_REALISED  Contract F9 on stand-in bodies (SK F1-F7 until J1).
%   Field identities with the S6/S7 inputs are asserted bitwise (same operands, same operation).
%   Oracle: an all-solid UHPC body equals sti_stage2 with every module at rho_uhpc (closed-form
%   prisms, independent assembly); differences are rounding of sums taken in another order.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
schema = mwecmass.output.export_schema();
want = setdiff({schema.final_props_fields.name}, {'stage3_status', 'stage3_check', 'realised_strips'});
for name = {'cylinder', 'box'}
    config = sti_config(name{1});
    geo = mwecmass.solid.outer_nurbs(config.ms2_model);
    config.hull_solid = geo;
    e = config.strip_edges;
    N = numel(e) - 1;
    vs = 0.5;
    t_min = 0.0762;
    t = 0.1;
    if strcmp(name{1}, 'box')
        inner = sti_inner_box(geo, t, t_min);
    else
        inner = mwecmass.solid.offset_surface(config.ms2_model, config.boundary_cache, geo, t, ...
            [e(1) e(end)], struct('t_min', t_min));
    end
    design = struct('mode', 'modular_precast', 'edges', e, 'vs', vs, 't', [NaN; t * ones(N - 1, 1)], ...
        'z_ballast', (e(2) + e(3)) / 2, 'solid_modules', 1);
    rho = struct('uhpc', 2500, 'air', 1.2);
    body = mwecmass.solid.build_body(geo, design, inner);
    bp = mwecmass.solid.body_properties(body, rho, struct());
    hs = mwecmass.solid.hydrostatics_at_draft(geo, vs, struct());
    p = mwecmass.realise.evaluate_realised(bp, hs, design, config);

    missing = setdiff(want, fieldnames(p));
    check(isempty(missing), sprintf('F9 lacks %s', strjoin(missing, ', ')));
    check(p.vertical_shift == vs && p.draft == hs.draft && p.V_sub == hs.V_sub && isequal(p.CB, hs.CB) ...
        && p.Aw == hs.Aw && p.KM == hs.KM && p.A_sub == hs.S_wet, 'S7 fields');
    check(p.mass_total == bp.total.mass && isequal(p.Inertia_Tensor, bp.total.I_cg) ...
        && p.Iyy == bp.total.I_cg(2, 2), 'S6 fields');
    check(isequal(p.CG_total, [0 0 bp.total.CG_body(3) + vs]), 'I6: CG_total(3) = CG_body(3) + vs');
    check(p.GM_L == hs.KM - p.CG_total(3), 'GM_L = KM - CG_total(3)');
    check(p.mass_buoyant_force == config.RHO_WATER * hs.V_sub && ...
        p.mass_discrepancy == p.mass_total - p.mass_buoyant_force, 'buoyancy');
    K33 = config.RHO_WATER * config.G * hs.Aw;
    check(isequal(p.K_hydro, diag([0, K33, p.mass_total * config.G * p.GM_L])) && ...
        isequal(p.K_total, p.K_hydro) && isequal(p.K_pto, zeros(3)), 'stiffness');
    [A11, A33, A55, A_full] = mwecmass.bem.interpolate_at_draft(vs, config, p.CG_total(3));
    check(isequal([p.A11 p.A33 p.A55], [A11 A33 A55]) && isequal(p.A_full, A_full), 'added mass at realised CG');
    T = mwecmass.hydrostatics.coupled_periods_by_share(diag([p.mass_total, p.mass_total, p.Iyy]) + A_full, p.K_total);
    check(p.periods.heave == T.heave && p.periods.pitch == T.pitch && isequal(p.coupled_periods, ...
        [T.surge; T.heave; T.pitch]) && p.periods.surge == Inf, 'coupled periods');
    check(p.periods.heave_uncoupled == 2 * pi * sqrt((p.mass_total + A33) / K33), 'uncoupled heave');
    rho_eff = [bp.modules.rho_eff]';
    check(isequal(p.densities_at_nodes, rho_eff) && isequal(p.realised_strip_density, rho_eff) && ...
        isequal(p.realised_strip_edges, e) && isequal([p.components.density]', rho_eff) && ...
        isequal([p.components.z_level]', (e(1:end - 1) + e(2:end)) / 2 + vs), 'module densities');
    check(isequal(p.cross_section, config.profile) && strcmp(p.fill_method, 'uhpc_fill') && ...
        strcmp(p.density_profile_source, 'realised_partition'), 'labels');
    fprintf(['%-8s M %.3f kg, Z_CG %.6f m, GM %.6f m, T_heave %.6f s (uncoupled %.6f), T_pitch %.6f s ' ...
        '(uncoupled %.6f), rho_eff %s\n'], name{1}, p.mass_total, p.CG_total(3), p.GM_L, p.periods.heave, ...
        p.periods.heave_uncoupled, p.periods.pitch, p.periods.pitch_uncoupled, mat2str(rho_eff', 6));

    % oracle: all-solid UHPC body against sti_stage2 at rho_uhpc in every module
    solid = struct('mode', 'modular_precast', 'edges', e, 'vs', vs, 't', NaN(N, 1), ...
        'z_ballast', e(1), 'solid_modules', 1:N);
    bs = mwecmass.solid.body_properties(mwecmass.solid.build_body(geo, solid, []), rho, struct());
    ps = mwecmass.realise.evaluate_realised(bs, hs, solid, config);
    f3 = sti_stage2(config, vs, rho.uhpc * ones(N, 1));
    d = [ps.mass_total / f3.mass_total - 1, ps.CG_total(3) / f3.CG_total(3) - 1, ...
        ps.Iyy / f3.Iyy - 1, ps.GM_L / f3.GM_L - 1, ps.periods.heave / f3.periods.heave - 1, ...
        ps.periods.pitch / f3.periods.pitch - 1];
    fprintf('%-8s all solid vs sti_stage2, relative differences M, Z_CG, Iyy, GM, T_heave, T_pitch: %s\n', ...
        name{1}, mat2str(d, 3));
    % N products of the same factors summed in another order: at most N rounding steps apart
    check(abs(d(1)) <= N * eps, 'all-solid mass against sti_stage2');
end
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
    error('test_evaluate_realised:fail', '%s', msg);
end
end
