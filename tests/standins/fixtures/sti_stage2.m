function [final3d, stage2] = sti_stage2(config, vs, rho)
%STI_STAGE2  Stage-2-like result (Final3D fields) of a strip-density design on a stand-in fixture.
%
%   [final3d, stage2] = sti_stage2(config, vs, rho)
%
%   config from sti_config; vs the vertical shift [m]; rho [N x 1] bulk density of every module
%   (config.strip_edges) [kg/m^3]. One NaN entry of rho is solved for flotation, mass = rho_w
%   V_sub. Mass properties of the full modules and the hydrostatics come from sti_closed_form
%   (prisms; source in its help), the added mass from mwecmass.bem.interpolate_at_draft at the CG
%   and the periods as in properties_3d.m sections 8, 11, 12 (K33 = rho_w g Aw, K55 = M g GM,
%   coupled periods by mwecmass.hydrostatics.coupled_periods_by_share). stage2 is the S8.stage2
%   struct (vs, rho, mass, Z_CG = CG_total(3), GM, T_heave, T_pitch; world frame).

fx = sti_closed_form('fixture', config.ms2_model);
e = config.strip_edges(:);
N = numel(e) - 1;
rho = rho(:);
design = struct('mode', 'modular_precast', 'edges', e, 'vs', vs, 't', NaN(N, 1), ...
    'z_ballast', fx.z(1), 'solid_modules', 1:N);
reg = sti_closed_form('regions', fx, design, NaN(N, 1));
V = reg.uhpc.V;
hs = sti_closed_form('hydrostatics', fx, vs);
k = find(isnan(rho));
if numel(k) == 1
    others = setdiff(1:N, k);
    rho(k) = (config.RHO_WATER * hs.V_sub - rho(others)' * V(others)) / V(k);
end
M = rho' * V;
S = rho' * reg.uhpc.S;
I0 = zeros(3);
for i = 1:N
    J = reg.uhpc.J(:, :, i);
    I0 = I0 + rho(i) * (trace(J) * eye(3) - J);
end
c = S / M;
I = I0 - M * (dot(c, c) * eye(3) - c' * c);
f = struct();
f.vertical_shift = vs;
f.draft = hs.draft;
f.mass_total = M;
f.CG_total = [0, 0, c(3) + vs];
f.Inertia_Tensor = I;
f.Ixx = I(1, 1);
f.Iyy = I(2, 2);
f.Izz = I(3, 3);
f.V_sub = hs.V_sub;
f.CB = hs.CB;
f.Aw = hs.Aw;
f.I_wp_xx = hs.I_wp_xx;
f.I_wp_yy = hs.I_wp_yy;
f.A_sub = hs.S_wet;
f.KM = hs.KM;
f.GM_L = hs.KM - f.CG_total(3);
f.mass_buoyant_force = config.RHO_WATER * hs.V_sub;
f.mass_discrepancy = M - f.mass_buoyant_force;
K33 = config.RHO_WATER * config.G * hs.Aw;
K55 = 0;
if f.GM_L > 0
    K55 = M * config.G * f.GM_L;
end
f.K_hydro = diag([0, K33, K55]);
[f.A11, f.A33, f.A55, f.A_full, f.B_full] = mwecmass.bem.interpolate_at_draft(vs, config, f.CG_total(3));
f.K_pto = zeros(3);
f.K_total = f.K_hydro;
f.periods = struct('surge', Inf, 'heave', Inf, 'pitch', Inf, ...
    'heave_uncoupled', uncoupled(M + f.A33, K33), 'pitch_uncoupled', uncoupled(f.Iyy + f.A55, K55));
[T, share, Phi] = mwecmass.hydrostatics.coupled_periods_by_share(diag([M, M, f.Iyy]) + f.A_full, f.K_total);
f.periods.heave = T.heave;
f.periods.pitch = T.pitch;
f.coupled_periods = [T.surge; T.heave; T.pitch];
f.coupled_modes = Phi;
f.participation_factors = share;
f.surge_per_pitch = Phi(1, 3) / Phi(3, 3);
f.densities_at_nodes = rho;
f.components = struct('density', num2cell(rho), 'z_level', num2cell((e(1:end - 1) + e(2:end)) / 2 + vs));
f.MassMatrix_CG = mwecmass.hydrostatics.build_mass_matrix(M, [0, 0, 0], I);
f.MassMatrix_Origin = mwecmass.hydrostatics.build_mass_matrix(M, f.CG_total, I);
f.cross_section = config.profile;
final3d = f;
stage2 = struct('vs', vs, 'rho', rho, 'mass', M, 'Z_CG', f.CG_total(3), 'GM', f.GM_L, ...
    'T_heave', f.periods.heave, 'T_pitch', f.periods.pitch);
end

function T = uncoupled(m, k)
if k > 1e-6
    T = 2 * pi * sqrt(m / k);
else
    T = Inf;
end
end
