function props = evaluate_realised(bp, hs, design, config)
%EVALUATE_REALISED  final_props of a realised Stage-3 body (contract F9).
%
%   props = mwecmass.realise.evaluate_realised(bp, hs, design, config)
%
%   bp: S6 body properties (mwecmass.solid.body_properties) of the realised body; hs: S7
%   hydrostatics (mwecmass.solid.hydrostatics_at_draft) at design.vs; design: S3. config: RHO_WATER,
%   G, profile and the hydrodynamic table read by mwecmass.bem.interpolate_at_draft.
%
%   Every field of the final_props schema (mwecmass.output.export_schema) except stage3_status and
%   stage3_check, which the realisation adds. Mass, CG and inertia come from the exact body,
%   hydrostatics from hs; CG_total = [0 0 CG_body(3) + vs] (world frame; CG_x = CG_y = 0 is the
%   declared symmetric-body model, the exact x and y stay in bp), Inertia_Tensor = I_cg,
%   GM_L = KM - CG_total(3). Stiffness and periods as properties_3d.m sections 8 to 12: K33 =
%   rho_w g Aw, K55 = M g GM when GM > 0 else 0, no surge stiffness; added mass and damping at the
%   realised draft and CG (interpolate_at_draft); uncoupled periods 2 pi sqrt(m/k), Inf when
%   k <= 1e-6; periods.heave and .pitch are the coupled periods (coupled_periods_by_share).
%   Densities per module are the realised rho_eff = mass / V of each module.

vs = design.vs;
M = bp.total.mass;
I = bp.total.I_cg;
e = design.edges(:);

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
    'heave_uncoupled', period(M + props.A33, K33), ...
    'pitch_uncoupled', period(props.Iyy + props.A55, K55));
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
props.realised_strip_edges = e;
props.components = struct('density', num2cell(rho_eff), ...
    'z_level', num2cell((e(1:end - 1) + e(2:end)) / 2 + vs));
end

function T = period(m, k)
% properties_3d.m section 11: no restoring stiffness above 1e-6 gives an infinite period.
if k > 1e-6
    T = 2 * pi * sqrt(m / k);
else
    T = Inf;
end
end
