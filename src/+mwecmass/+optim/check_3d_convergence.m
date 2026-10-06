function [converged, quality_metrics] = check_3d_convergence(x_opt_3d, config, exitflag, output)
%CHECK_3D_CONVERGENCE Check Stage-2 density, flotation, GM, periods, and solver quality.
% x_opt_3d is [vertical_shift, rho_1..rho_N]; config supplies limits and hydrostatics.
% output must contain fmincon constrviolation, firstorderopt, and iterations.
% The wall-strip boundary is excluded from monotonicity when configured; exitflags 1, 2,
% and a bounded-residual 0 are accepted. Outputs are a logical and diagnostic struct.
props_3d  = mwecmass.hydrostatics.properties_3d(x_opt_3d, config);
densities = x_opt_3d(2:end);

% Match the constraint function by excluding the pinned wall/platform boundary.
if isfield(config, 'wall_strip_index') && ~isempty(config.wall_strip_index)
    w       = config.wall_strip_index;
    N_dens  = length(densities);
    % Build the adjacent-pair list while skipping pairs touching the wall strip.
    pairs_to_check = [];
    for ii = 1:(N_dens - 1)
        if ii == w || (ii + 1) == w
            continue;
        end
        pairs_to_check(end+1) = ii; %#ok<AGROW>
    end
    if isempty(pairs_to_check)
        monotonic_ok = true;
    else
        d = diff(densities);
        monotonic_ok = all(d(pairs_to_check) <= 0);
    end
else
    monotonic_ok = all(diff(densities) <= 0);
end

mass_balance_error = abs(props_3d.mass_total - props_3d.mass_buoyant_force);
mass_balance_ok    = mass_balance_error < 10;
GM_constraint_ok   = props_3d.GM_L > config.gm_min;
GM_margin          = props_3d.GM_L - config.gm_min;

period_heave_error = abs(props_3d.periods.heave - config.T_heave_goal) / config.T_heave_goal;
period_pitch_error = abs(props_3d.periods.pitch - config.T_pitch_goal) / config.T_pitch_goal;
periods_ok = (period_heave_error < 0.10) && (period_pitch_error < 0.10);

% Accept first-order (1), step-tolerance (2), or bounded-residual iteration-limit (0) exits.
fmincon_optimal    = (exitflag == 1);
fmincon_acceptable = (ismember(exitflag, [0, 2]) && ...
    output.constrviolation < 1e-6  && ...
    output.firstorderopt   < 1e-2);
fmincon_ok = fmincon_optimal || fmincon_acceptable;

converged = monotonic_ok && mass_balance_ok && GM_constraint_ok && fmincon_ok;

quality_metrics = struct( ...
    'monotonic',            monotonic_ok, ...
    'mass_balance',         mass_balance_ok, ...
    'mass_balance_error_kg', mass_balance_error, ...
    'GM_satisfied',         GM_constraint_ok, ...
    'GM_margin',            GM_margin, ...
    'periods_acceptable',   periods_ok, ...
    'heave_error_pct',      100*period_heave_error, ...
    'pitch_error_pct',      100*period_pitch_error, ...
    'fmincon_optimal',      fmincon_optimal, ...
    'fmincon_acceptable',   fmincon_acceptable, ...
    'exitflag',             exitflag, ...
    'firstorderopt',        output.firstorderopt, ...
    'constrviolation',      output.constrviolation, ...
    'iterations',           output.iterations);
end
