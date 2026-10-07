function props = evaluate_realised(bp, hs, design, config)
%EVALUATE_REALISED  Stand-in of contract F9: final_props of a stand-in body.
%
%   props = mwecmass.realise.evaluate_realised(bp, hs, design, config)
%
%   bp: S6, hs: S7, design: S3, config.hull_solid: the S1b of a fixture deck, identified by
%   sti_closed_form('fixture', ...) through its analytic field or its hull_name, so the real F1 geo
%   of a fixture is accepted; any other geo errors mwecmass:standin:NotAnalytic. Fields and formulas
%   of contract F9 (mass, CG_total = [0 0 CG_body(3) + vs], I_cg, GM_L = KM - CG_total(3),
%   K33 = rho_w g Aw, K55 = M g GM when GM > 0, uncoupled periods 2 pi sqrt((M + A33)/K33) and
%   2 pi sqrt((I_yy + A55)/K55), Inf when K <= 1e-6, as properties_3d.m section 11), with the
%   added mass of mwecmass.bem.interpolate_at_draft at the realised CG and the coupled periods of
%   mwecmass.hydrostatics.coupled_periods_by_share, as build_realised_properties.m. props.analytic
%   marks stand-in output for the check_against_stage2 stand-in; the real F9 does not set it.

if ~isstruct(config) || ~isfield(config, 'hull_solid') || ~isstruct(config.hull_solid)
    error('mwecmass:standin:NotAnalytic', 'evaluate_realised stand-in: config.hull_solid is not a stand-in fixture geo');
end
if exist('sti_closed_form', 'file') ~= 2
    addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'fixtures'), '-end');
end
fx = sti_closed_form('fixture', config.hull_solid);
vs = design.vs;
M = bp.total.mass;
I = bp.total.I_cg;
props = struct();
props.vertical_shift = vs;
props.draft = hs.draft;
props.Aw = hs.Aw;
props.I_wp_yy = hs.I_wp_yy;
props.I_wp_xx = hs.I_wp_xx;
props.V_sub = hs.V_sub;
props.CB = hs.CB;
props.A_sub = hs.S_wet;
props.mass_buoyant_force = config.RHO_WATER * hs.V_sub;
props.KM = hs.KM;
props.mass_total = M;
props.CG_total = [0, 0, bp.total.CG_body(3) + vs];
props.Ixx = I(1, 1);
props.Iyy = I(2, 2);
props.Izz = I(3, 3);
props.Inertia_Tensor = I;
props.GM_L = hs.KM - props.CG_total(3);
props.mass_discrepancy = M - props.mass_buoyant_force;
K33 = config.RHO_WATER * config.G * hs.Aw;
K55 = 0;
if props.GM_L > 0
    K55 = M * config.G * props.GM_L;
end
props.K_hydro = diag([0, K33, K55]);
[props.A11, props.A33, props.A55, props.A_full, props.B_full] = ...
    mwecmass.bem.interpolate_at_draft(vs, config, props.CG_total(3));
props.K_pto = zeros(3);
props.K_total = props.K_hydro;
props.periods = struct('surge', Inf, 'heave', Inf, 'pitch', Inf, ...
    'heave_uncoupled', uncoupled(M + props.A33, K33), 'pitch_uncoupled', uncoupled(props.Iyy + props.A55, K55));
[T, share, Phi] = mwecmass.hydrostatics.coupled_periods_by_share( ...
    diag([M, M, props.Iyy]) + props.A_full, props.K_total);
props.periods.heave = T.heave;
props.periods.pitch = T.pitch;
props.coupled_periods = [T.surge; T.heave; T.pitch];
props.coupled_modes = Phi;
props.participation_factors = share;
props.surge_per_pitch = Phi(1, 3) / Phi(3, 3);
props.MassMatrix_CG = mwecmass.hydrostatics.build_mass_matrix(M, [0, 0, 0], I);
props.MassMatrix_Origin = mwecmass.hydrostatics.build_mass_matrix(M, props.CG_total, I);
if strcmp(design.mode, 'modular_precast')
    props.fill_method = 'uhpc_fill';
else
    props.fill_method = 'steel_fill';
end
props.density_profile_source = 'realised_partition';
props.cross_section = config.profile;
rho_eff = [bp.modules.rho_eff]';
props.densities_at_nodes = rho_eff;
props.realised_strip_density = rho_eff;
props.realised_strip_edges = design.edges(:);
e = design.edges(:);
props.components = struct('density', num2cell(rho_eff), 'z_level', num2cell((e(1:end - 1) + e(2:end)) / 2 + vs));
props.analytic = fx;
end

function T = uncoupled(m, k)
if k > 1e-6
    T = 2 * pi * sqrt(m / k);
else
    T = Inf;
end
end
