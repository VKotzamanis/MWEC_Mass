function f = stage2_objective(x, config)
%STAGE2_OBJECTIVE Sum range-normalized GM, heave, and pitch penalties.
% x=[vertical_shift, rho_1..rho_N]. Each residual is (actual-target)/max(half_range,eps), making a
% boundary penalty equal across quantities; range_penalty supplies the C1 outer curvature.
% Invalid GM or periods return config.penalty_guard. Mass and GM feasibility remain constraints.
% See docs/METHODS_ENGINE.md#optim-stage2-formulation
props = mwecmass.hydrostatics.properties_3d(x, config);

% Guard: return penalty for diverging solutions
if isnan(props.GM_L)
    f = config.penalty_guard;
    return;
end
if isnan(props.periods.heave) || isinf(props.periods.heave) || ...
   isnan(props.periods.pitch) || isinf(props.periods.pitch)
    f = config.penalty_guard;
    return;
end

% Half-ranges (computed from driver-specified ranges)
gm_half    = 0.5 * (config.gm_range(2)      - config.gm_range(1));
heave_half = 0.5 * (config.T_heave_range(2)  - config.T_heave_range(1));
pitch_half = 0.5 * (config.T_pitch_range(2)  - config.T_pitch_range(1));

% Normalised errors: r = 0 at target, r = ±1 at range boundary
r_gm    = (props.GM_L           - config.gm_target)     / gm_half;
r_heave = (props.periods.heave  - config.T_heave_goal)  / heave_half;
r_pitch = (props.periods.pitch  - config.T_pitch_goal)  / pitch_half;

% Piecewise quadratic with amplification outside |r| > 1
k_amp = config.zone_k_amp;
f = mwecmass.optim.range_penalty(r_gm, k_amp) ...
  + mwecmass.optim.range_penalty(r_heave, k_amp) ...
  + mwecmass.optim.range_penalty(r_pitch, k_amp);
end
