function [x_opt, props_2d, conv_data, hist, ...
          converged, iteration, pid_saturated] = ...
        stage1_trained(x0, config, hist)
%STAGE1_TRAINED Iteratively correct the 2-D surrogate with volume and CG PID channels.
% k_vol responds to the 3-D/2-D submerged-volume gap; k_gm responds to the CG ratio. Both errors are
% dimensionless and updates are bounded/damped. Returns the design, properties, histories, convergence,
% iteration count, and saturation flag; optional HAMS enrichment is controlled by configuration.
% See docs/METHODS_ENGINE.md#optim-pid-correction
    fprintf(' Running Stage 1 with PID surrogate training.\n');

    % --- Unpack PID configuration from config (set by driver) ---
    bounds_kvol = config.bounds_kvol;
    bounds_kgm  = config.bounds_kgm;

    gv = config.pid_vol_gains;
    gg = config.pid_gm_gains;
    gm = config.pid_mass_gains;

    vol_pid  = mwecmass.optim.PID_Controller(gv(1), gv(2), gv(3), ...
               'OutputLimits', config.pid_vol_limits);
    gm_pid   = mwecmass.optim.PID_Controller(gg(1), gg(2), gg(3), ...
               'OutputLimits', config.pid_gm_limits);
    mass_pid = mwecmass.optim.PID_Controller(gm(1), gm(2), gm(3), ...
               'OutputLimits', config.pid_mass_limits);

    max_iters    = config.max_outer_iterations;
    stable_count = 0;

    converged     = false;
    iteration     = 0;
    pid_saturated = false;
    double_sat_count = 0;   % consecutive iterations with both PIDs saturated

    fprintf(' Bounds: k_vol [%.2f, %.2f], k_gm [%.2f, %.2f]\n', ...
            bounds_kvol, bounds_kgm);
    fprintf(' Max iterations: %d\n', max_iters);

    % Track cache size so unchanged hydro data is not rewritten at the end.
    hydro_cache_entries_before = numel(config.hydro_cache.drafts);

    while ~converged && iteration < max_iters
        iteration = iteration + 1;
        fprintf('\n  [Stage1 %2d/%d] k_vol=%.3f k_gm=%.3f ', ...
                iteration, max_iters, config.k_vol, config.k_gm);

        % --- 2D optimisation with current correction factors ---
        [x_opt, props_2d, exitflag, conv_data] = ...
            mwecmass.optim.solve_2d_surrogate(config, x0, mass_pid);

        if exitflag <= 0
            fprintf('  2D exit=%d. ', exitflag);
        end

        % Enrich live HAMS data when configured; otherwise validation uses cache interpolation.
        % Rebuild hydrodynamic fields after enrichment so the current draft is available immediately.
        vs_converged = x_opt(1);
        if strcmp(mwecmass.optim.hams_enrichment_action(config), 'run')
            [~, config.hydro_cache] = mwecmass.bem.get_or_run_hydro( ...
                vs_converged, config.hydro_cache, config, ...
                config.hams_dir, config.hams_exe, 0.01);
            config = mwecmass.bem.rebuild_config_hydro(config, config.hydro_cache);
        else
            fprintf(['  HAMS enrichment skipped (run_HAMS_MREL=false): 3D validation uses the ' ...
                     'cache interpolated at vs=%+.4f m.\n'], vs_converged);
        end

        % Guard: warn if A(inf) is all zeros (HAMS may have failed)
        if ~isempty(config.added_mass_diagonal) && max(abs(config.added_mass_diagonal(:))) < 1e-6
            warning('WEC:ZeroAddedMass', ...
                    'A(inf) = 0 at all drafts. HAMS failed — periods will be wrong.');
        end

        % --- 3D validation pass 1: get actual CG from density distribution ---
        p3 = mwecmass.hydrostatics.properties_3d(x_opt, config);

        % Re-transform hydrodynamics at the non-uniform-density CG before validation.
        if ~isempty(config.hams_dir)
            actual_cg_z = p3.CG_total(3);
            config = mwecmass.bem.retransform_at_cg( ...
                config, vs_converged, actual_cg_z);
            p3 = mwecmass.hydrostatics.properties_3d(x_opt, config);
        end

        vsub_2d = props_2d.V_sub;        vsub_3d = p3.V_sub;
        mass_2d = props_2d.mass_total;    mass_3d = p3.mass_total;
        gm_2d   = props_2d.GM;           gm_3d   = p3.GM_L;
        cg_z_2d = props_2d.CG_total(3);  cg_z_3d = p3.CG_total(3);

        vol_error_pct  = 100 * (vsub_2d - vsub_3d) / vsub_3d;
        mass_error_pct = 100 * (mass_2d - mass_3d) / mass_3d;
        gm_error       = gm_2d - gm_3d;

        % --- Append history ---
        hist.vsub_2d_history(end+1)       = vsub_2d;
        hist.vsub_3d_history(end+1)       = vsub_3d;
        hist.vol_errors_history(end+1)    = vol_error_pct;
        hist.mass_2d_history(end+1)       = mass_2d;
        hist.mass_3d_history(end+1)       = mass_3d;
        hist.mass_2d_corrected_history(end+1) = mass_2d;
        hist.gm_2d_history(end+1)         = gm_2d;
        hist.gm_3d_history(end+1)         = gm_3d;
        hist.gm_errors_history(end+1)     = gm_error;
        hist.cg_z_2d_history(end+1)       = cg_z_2d;
        hist.cg_z_3d_history(end+1)       = cg_z_3d;
        hist.kvol_history(end+1)          = config.k_vol;
        hist.kgm_history(end+1)           = config.k_gm;
        hist.mass_errors_history(end+1)   = mass_2d - mass_3d;
        hist.mass_correction_history(end+1) = config.k_vol;
        hist.gm_correction_history(end+1)   = config.k_gm;

        % --- Console report (compact) ---
        fprintf('\n    V: %+.1f%%  M: %+.1f%%  GM: %.3f→%.3f  T_h: %.1f→%.1f  T_p: %.1f→%.1f\n', ...
                vol_error_pct, mass_error_pct, gm_2d, gm_3d, ...
                props_2d.periods.heave, p3.periods.heave, ...
                props_2d.periods.pitch, p3.periods.pitch);

        % Update the volume scaling factor multiplicatively to preserve positivity.
        kvol_before = config.k_vol;
        kgm_before  = config.k_gm;

        if vsub_2d > 1e-6
            vol_err_norm = (vsub_3d - vsub_2d) / vsub_3d;
            vol_update   = vol_pid.update(vol_err_norm, 1.0);
        else
            vol_update = 0;
        end

        % Under-relaxation: cautious early, more aggressive once trend is clear
        if iteration <= config.damping_transition_iter
            damp_v = config.damping_vol_early;
        else
            damp_v = config.damping_vol_late;
        end
        kvol_new     = config.k_vol * (1 + damp_v * vol_update);
        config.k_vol = max(bounds_kvol(1), min(bounds_kvol(2), kvol_new));

        % Skip the CG-ratio update near the waterline, where the ratio is unstable.
        if abs(cg_z_2d) > config.cg_guard_floor && ...
           abs(cg_z_3d) > config.cg_guard_floor
            target_kgm      = cg_z_3d / cg_z_2d;
            gm_err_norm     = target_kgm - 1.0;
            gm_update       = gm_pid.update(gm_err_norm, 1.0);
        else
            gm_update = 0;
        end

        if iteration <= config.damping_transition_iter
            damp_g = config.damping_gm_early;
        else
            damp_g = config.damping_gm_late;
        end
        kgm_new     = config.k_gm * (1 + damp_g * gm_update);
        config.k_gm = max(bounds_kgm(1), min(bounds_kgm(2), kgm_new));

        % Accept matching errors or stable correction factors after warm-up iterations.
        delta_kvol = abs(config.k_vol - kvol_before);
        delta_kgm  = abs(config.k_gm  - kgm_before);

        vol_ok = abs(vol_error_pct) < config.vol_conv_tol_pct;
        gm_ok  = abs(gm_error)      < config.gm_conv_tol;

        if delta_kvol < config.delta_kvol_stable && ...
           delta_kgm  < config.delta_kgm_stable
            stable_count = stable_count + 1;
        else
            stable_count = 0;
        end
        kvol_stable = (stable_count >= config.stable_count_needed);

        % Saturation detection: correction factor stuck at its bound
        prox = config.pid_sat_proximity;
        pid_sat_vol = abs(config.k_vol - bounds_kvol(2)) < prox || ...
                      abs(config.k_vol - bounds_kvol(1)) < prox;
        pid_sat_gm  = abs(config.k_gm  - bounds_kgm(2)) < prox || ...
                      abs(config.k_gm  - bounds_kgm(1)) < prox;
        pid_saturated = pid_sat_vol || pid_sat_gm;

        fprintf('    PID: Δk_vol=%+.4f%s  Δk_gm=%+.4f%s\n', ...
                delta_kvol, sat_tag(pid_sat_vol), ...
                delta_kgm, sat_tag(pid_sat_gm));

        % --- Early exit: both PIDs saturated for 3 consecutive iters ---
        %  When both correction factors are at their bounds AND the
        %  solution is not changing, further iterations are wasted.
        %  Break early and let Stage 2 work from the best-effort warm-start.
        if pid_sat_vol && pid_sat_gm
            double_sat_count = double_sat_count + 1;
            if double_sat_count >= 3
                fprintf('\n  EARLY EXIT: both PIDs saturated for %d consecutive iterations.\n', ...
                        double_sat_count);
                fprintf('    k_vol=%.4f (bound), k_gm=%.4f (bound)\n', ...
                        config.k_vol, config.k_gm);
                fprintf('    Passing best-effort solution to Stage 2.\n');
                break;
            end
        else
            double_sat_count = 0;
        end

        if exitflag > 0 && ((vol_ok && gm_ok) || (kvol_stable && iteration > 3))
            converged = true;
            fprintf('\n  SURROGATE CONVERGED at iteration %d\n', iteration);
            fprintf('    V_sub err: %.1f%%, GM err: %.3f m\n', vol_error_pct, gm_error);
        elseif exitflag <= 0 && (vol_ok && gm_ok)
            % 2D→3D errors are within tolerance, but fmincon itself
            % failed (infeasible / stalled).  Do NOT declare convergence;
            % the warm-start quality is suspect.
            fprintf('\n  2D errors within tolerance but fmincon failed (exit=%d) — continuing\n', exitflag);
            x0 = x_opt;
        else
            x0 = x_opt;   % warm-start next iteration
        end

        % --- Bookkeeping for results struct ---
        hist.constraints_satisfied_history(end+1) = true;
        hist.solution_stable_history(end+1)       = kvol_stable;
        hist.mass_acceptable_history(end+1)       = ...
            abs(mass_error_pct) < config.mass_acceptable_pct;
        hist.converged_history(end+1)             = converged;
        hist.fitness_scores_history(end+1)        = 1 - abs(vol_error_pct)/100;
        hist.GM_margins_history(end+1)            = ...
            max(0, (gm_3d - config.gm_min) / config.gm_min);

        % Compute cumulative surrogate-accuracy statistics; one point is insufficient for R²/MAPE.
        n_pts = length(hist.mass_2d_history);
        if n_pts >= 2
            hist.R2_mass_history(end+1)   = compute_R2(hist.mass_3d_history, hist.mass_2d_history);
            hist.R2_GM_history(end+1)     = compute_R2(hist.gm_3d_history,   hist.gm_2d_history);
            hist.R2_cg_history(end+1)     = compute_R2(hist.cg_z_3d_history, hist.cg_z_2d_history);
            hist.MAPE_mass_history(end+1) = compute_MAPE(hist.mass_3d_history, hist.mass_2d_history);
            hist.MAPE_GM_history(end+1)   = compute_MAPE(hist.gm_3d_history,   hist.gm_2d_history);
            hist.MAPE_cg_history(end+1)   = compute_MAPE(hist.cg_z_3d_history, hist.cg_z_2d_history);
            hist.ME_mass_history(end+1)   = mean(hist.mass_2d_history - hist.mass_3d_history);
            hist.ME_GM_history(end+1)     = mean(hist.gm_2d_history   - hist.gm_3d_history);
        else
            hist.R2_mass_history(end+1)   = NaN;
            hist.R2_GM_history(end+1)     = NaN;
            hist.R2_cg_history(end+1)     = NaN;
            hist.MAPE_mass_history(end+1) = NaN;
            hist.MAPE_GM_history(end+1)   = NaN;
            hist.MAPE_cg_history(end+1)   = NaN;
            hist.ME_mass_history(end+1)   = NaN;
            hist.ME_GM_history(end+1)     = NaN;
        end
    end  % while

    % Save the enriched cache only when new entries were added and rewriting is enabled.
    hydro_cache_entries_after = numel(config.hydro_cache.drafts);
    if ~isempty(config.hams_dir) && ~isempty(config.hydro_cache_file)
        cache_gained_entries = hydro_cache_entries_after > hydro_cache_entries_before;
        cache_save_enabled   = ~isfield(config, 'output') || config.output.save.hydrodynamics.cache_rewrite;
        if cache_gained_entries && cache_save_enabled
            hydro_table = config.hydro_cache;
            mwecmass.output.save_hydro_cache(config.hydro_cache_file, hydro_table);
            fprintf('  Hydro cache saved: %d entries → %s\n', ...
                    hydro_cache_entries_after, config.hydro_cache_file);
        elseif cache_gained_entries
            fprintf('  Hydro cache gained %d entries but out.save.hydrodynamics.cache_rewrite is false; not rewritten\n', ...
                    hydro_cache_entries_after);
        else
            fprintf('  Hydro cache unchanged (%d entries); not rewritten\n', ...
                    hydro_cache_entries_after);
        end
    end

    if ~converged
        if pid_saturated
            warning('WEC:PIDSaturated', ...
                'PID saturated: k_vol=%.3f, k_gm=%.3f', ...
                config.k_vol, config.k_gm);
        else
            warning('WEC:NotConverged', ...
                'Surrogate did not converge after %d iterations', iteration);
        end
    end
end


function R2 = compute_R2(y_true, y_pred)
% COMPUTE_R2  Coefficient of determination.
%   Returns NaN when total variance is negligible (avoids 0/0).

    y_true = y_true(:);  y_pred = y_pred(:);
    SS_res = sum((y_true - y_pred).^2);
    SS_tot = sum((y_true - mean(y_true)).^2);
    if SS_tot < 1e-12, R2 = nan; else, R2 = 1 - SS_res / SS_tot; end
end

function MAPE = compute_MAPE(y_true, y_pred)
% COMPUTE_MAPE  Mean absolute percentage error (%).
%   Ignores entries where y_true ≈ 0 (Inf contribution).

    y_true = y_true(:);  y_pred = y_pred(:);
    pct = abs((y_true - y_pred) ./ y_true) * 100;
    MAPE = mean(pct(isfinite(pct)));
end

function s = sat_tag(is_saturated)
% SAT_TAG  Returns the saturation status text when PID is at its bound.
    if is_saturated, s = ' SATURATED'; else, s = ''; end
end
