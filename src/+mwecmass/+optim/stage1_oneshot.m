function [x_opt, props_2d, conv_data, hist, iteration] = ...
        stage1_oneshot(x0, config, hist)
%STAGE1_ONESHOT Run one 2-D surrogate solve and compare it with 3-D properties.
% x0 is [draft; density nodes]; config supplies gains, bounds, and hydrostatics.
% Returns the selected design, properties, convergence data, history, and iteration=1.
    fprintf(' Running Stage 1 once (k_vol = %.1f, k_gm = %.1f)\n\n', ...
            config.k_vol, config.k_gm);

    % Mass PID: gains and limits from config (forwarded from driver)
    g = config.pid_mass_gains;
    mass_pid = mwecmass.optim.PID_Controller(g(1), g(2), g(3), ...
        'OutputLimits', config.pid_mass_limits);

    [x_opt, props_2d, exitflag, conv_data] = ...
        mwecmass.optim.solve_2d_surrogate(config, x0, mass_pid);

    if exitflag <= 0
        warning('WEC:Stage1Failed', ...
                'Stage 1 failed (exitflag=%d). Using x0.', exitflag);
        x_opt    = x0;
        props_2d = mwecmass.hydrostatics.properties_2d(x0, config);
    end

    iteration = 1;
    p3 = mwecmass.hydrostatics.properties_3d(x_opt, config);

    % Store single-iteration history
    hist.mass_errors_history          = props_2d.mass_total - p3.mass_total;
    hist.gm_errors_history            = props_2d.GM_uncorrected - p3.GM_L;
    hist.mass_correction_history      = config.k_vol;
    hist.gm_correction_history        = config.k_gm;
    hist.mass_2d_history              = props_2d.mass_total;
    hist.mass_3d_history              = p3.mass_total;
    hist.mass_2d_corrected_history    = props_2d.mass_total;
    hist.gm_2d_history                = props_2d.GM_uncorrected;
    hist.gm_3d_history                = p3.GM_L;
    hist.cg_z_2d_history              = props_2d.CG_total(3);
    hist.cg_z_3d_history              = p3.CG_total(3);

    % One-shot: single data point → R²/MAPE are undefined
    hist.R2_mass_history   = NaN;  hist.R2_GM_history   = NaN;  hist.R2_cg_history   = NaN;
    hist.MAPE_mass_history = NaN;  hist.MAPE_GM_history = NaN;  hist.MAPE_cg_history = NaN;
    hist.ME_mass_history   = NaN;  hist.ME_GM_history   = NaN;

    % Convergence metrics. Floor mass_total_denom at 1 kg to prevent Inf/NaN when mass_total is near-zero.
    mass_total_denom = max(p3.mass_total, 1);
    hist.constraints_satisfied_history = true;
    hist.solution_stable_history       = true;
    hist.mass_acceptable_history       = abs(100 * hist.mass_errors_history / mass_total_denom) < config.mass_acceptable_pct;
    % converged_history true iff fmincon succeeded (exitflag > 0); was unconditionally true, misreporting when fallback was used.
    hist.converged_history             = (exitflag > 0);
    hist.fitness_scores_history        = 1 - abs(hist.mass_errors_history / mass_total_denom);
    hist.GM_margins_history            = max(0, (p3.GM_L - config.gm_min) / config.gm_min);

    mass_gap_pct = 100 * hist.mass_errors_history / mass_total_denom;
    fprintf('\n  One-Shot Results:\n');
    fprintf('    Draft: %.3f m\n', x_opt(1));
    fprintf('    2D→3D Mass: %.1f → %.1f kg (%+.1f%%)\n', ...
            props_2d.mass_total, p3.mass_total, mass_gap_pct);
    fprintf('    2D→3D GM:   %.3f → %.3f m\n', ...
            props_2d.GM_uncorrected, p3.GM_L);
    fprintf('    2D→3D T_h:  %.2f → %.2f s\n', ...
            props_2d.periods.heave, p3.periods.heave);
    fprintf('    2D→3D T_p:  %.2f → %.2f s\n', ...
            props_2d.periods.pitch, p3.periods.pitch);
end
