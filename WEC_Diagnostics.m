classdef WEC_Diagnostics
    % WEC_DIAGNOSTICS  Post-run diagnostic reporting for the WEC pipeline.
    %
    %   Produces clear, tabulated console output that answers:
    %     1. Did the 2D surrogate converge?
    %     2. Did the 3D optimiser converge?
    %     3. Are the 2D → 3D → strip-sum properties self-consistent?
    %     4. Is the constructability realization feasible?
    %     5. Are mass/buoyancy balanced, and is the hull stable?
    %
    %   USAGE (call from WEC_Driver after the pipeline completes):
    %     WEC_Diagnostics.full_report(results, final_props);
    %
    %   All methods are static (no state).

    methods (Static)

        function full_report(results, final_props)
        % FULL_REPORT  Master diagnostic — prints all panels.
        %   Each panel is wrapped in try-catch so a single failure
        %   does not prevent the remaining panels from printing.

            config = results.config;

            try WEC_Diagnostics.print_header(config);
            catch ME, fprintf('  [Header panel failed: %s]\n', ME.message); end

            try WEC_Diagnostics.print_stage1_summary(results);
            catch ME, fprintf('  [Stage 1 panel failed: %s]\n', ME.message); end

            try WEC_Diagnostics.print_stage2_summary(results, final_props);
            catch ME, fprintf('  [Stage 2 panel failed: %s]\n', ME.message); end

            try WEC_Diagnostics.print_property_comparison(results, final_props);
            catch ME, fprintf('  [Property comparison panel failed: %s]\n', ME.message); end

            try WEC_Diagnostics.print_stability_check(final_props, config);
            catch ME, fprintf('  [Stability panel failed: %s]\n', ME.message); end

            try WEC_Diagnostics.print_density_profile(results, config);
            catch ME, fprintf('  [Density profile panel failed: %s]\n', ME.message); end

            if isfield(config, 'enable_constructability') && ...
                    config.enable_constructability && ...
                    isfield(results, 'constructability') && ...
                    ~isempty(results.constructability)
                try WEC_Diagnostics.print_constructability_report(results, config);
                catch ME, fprintf('  [Constructability panel failed: %s]\n', ME.message); end
            end

            try WEC_Diagnostics.print_verdict(results, final_props, config);
            catch ME, fprintf('  [Verdict panel failed: %s]\n', ME.message); end
        end


        %% ═══════════════════════════════════════════════════════════════
        %%  PANEL 0: HEADER
        %% ═══════════════════════════════════════════════════════════════

        function print_header(config)
            fprintf('\n');
            fprintf('╔══════════════════════════════════════════════════════════════════╗\n');
            fprintf('║              WEC PIPELINE — DIAGNOSTIC REPORT                   ║\n');
            fprintf('╠══════════════════════════════════════════════════════════════════╣\n');
            if isfield(config, 'ms2_file')
                geom_str = config.ms2_file;
            else
                geom_str = '(not stored in config)';
            end
            fprintf('║  MS2 file       : %-45s ║\n', geom_str);
            fprintf('║  Mode           : %-45s ║\n', config.stage1_mode);
            fprintf('║  Strips         : %-45s ║\n', ...
                sprintf('%d  (constructability: %s)', ...
                    config.num_ballast_sections, ...
                    ternary(config.enable_constructability, 'ON', 'OFF')));
            fprintf('║  Targets        : T_h=%.1f s, T_p=%.1f s, GM=%.2f m %18s║\n', ...
                config.T_heave_goal, config.T_pitch_goal, config.gm_target, '');
            fprintf('╚══════════════════════════════════════════════════════════════════╝\n');
        end


        %% ═══════════════════════════════════════════════════════════════
        %%  PANEL 1: STAGE 1 (2D SURROGATE) CONVERGENCE
        %% ═══════════════════════════════════════════════════════════════

        function print_stage1_summary(results)
            s1 = results.stage1_2d;
            config = results.config;

            fprintf('\n┌──────────────────────────────────────────────────────────────────┐\n');
            fprintf('│  STAGE 1: 2D SURROGATE  (%s)                              \n', upper(config.stage1_mode));
            fprintf('├──────────────────────────────────────────────────────────────────┤\n');

            fprintf('│  Iterations : %-4d    Converged: %-3s    PID saturated: %-3s\n', ...
                s1.iterations, ...
                ternary(s1.converged, 'YES', 'NO'), ...
                ternary(s1.pid_saturated, 'YES', 'no'));

            % k_vol and k_gm history
            if ~isempty(s1.mass_corrections)
                k_vol_final = s1.mass_corrections(end);
                k_gm_final  = s1.gm_corrections(end);
                fprintf('│  k_vol final: %.4f    k_gm final: %.4f\n', k_vol_final, k_gm_final);
            end

            % Final 2D→3D gaps
            if ~isempty(s1.mass_errors) && ~isempty(s1.gm_errors)
                mass_gap_kg = s1.mass_errors(end);
                gm_gap_m    = s1.gm_errors(end);

                % Mass gap as percentage of 3D mass
                if ~isempty(s1.mass_3d_history) && s1.mass_3d_history(end) > 0
                    mass_gap_pct = 100 * mass_gap_kg / s1.mass_3d_history(end);
                else
                    mass_gap_pct = NaN;
                end

                fprintf('│  Mass gap   : %+.1f kg  (%+.1f%%)\n', mass_gap_kg, mass_gap_pct);
                fprintf('│  GM gap     : %+.4f m\n', gm_gap_m);
            end

            % Surrogate accuracy
            if ~isempty(s1.R2_mass) && ~isnan(s1.R2_mass(end))
                fprintf('│  R²(mass)   : %.3f    R²(GM): %.3f    R²(CG_z): %.3f\n', ...
                    s1.R2_mass(end), s1.R2_GM(end), s1.R2_cg(end));
            end

            % Volume error
            if isfield(s1, 'convergence_data') && ...
                    isstruct(s1.convergence_data) && ...
                    isfield(s1.convergence_data, 'final_exitflag')
                fprintf('│  fmincon exit: %d\n', s1.convergence_data.final_exitflag);
            end

            fprintf('└──────────────────────────────────────────────────────────────────┘\n');
        end


        %% ═══════════════════════════════════════════════════════════════
        %%  PANEL 2: STAGE 2 (3D FMINCON) CONVERGENCE
        %% ═══════════════════════════════════════════════════════════════

        function print_stage2_summary(results, final_props)
            s2 = results.stage2_3d;
            config = results.config;

            fprintf('\n┌──────────────────────────────────────────────────────────────────┐\n');
            fprintf('│  STAGE 2: 3D HIGH-FIDELITY REFINEMENT                           \n');
            fprintf('├──────────────────────────────────────────────────────────────────┤\n');
            fprintf('│  Exit flag  : %-4d    Converged: %-3s\n', ...
                s2.exitflag, ternary(s2.converged, 'YES', 'NO'));
            fprintf('│  Cost (fval): %.6f\n', s2.fval);

            if isfield(s2, 'output') && isfield(s2.output, 'iterations')
                fprintf('│  SQP iters  : %d    FunEvals: %d\n', ...
                    s2.output.iterations, s2.output.funcCount);
            end

            % Mass balance at solution
            mass_err_pct = 100 * abs(final_props.mass_total - final_props.mass_buoyant_force) ...
                            / max(final_props.mass_total, 1);
            fprintf('│  Mass error : %.2e%%  %s\n', mass_err_pct, ...
                ternary(mass_err_pct < 0.1, '(OK)', '** VIOLATED **'));

            % True draft = keel depth below waterline
            if isfield(final_props, 'draft') && final_props.draft > 0
                draft_true = final_props.draft;
            elseif isfield(config, 'hull_z_min')
                draft_true = abs(config.hull_z_min + final_props.vertical_shift);
            else
                draft_true = 0;
            end
            fprintf('│  Draft      : %.4f m  (shift: %.4f m)\n', draft_true, final_props.vertical_shift);
            fprintf('└──────────────────────────────────────────────────────────────────┘\n');
        end


        %% ═══════════════════════════════════════════════════════════════
        %%  PANEL 3: PROPERTY COMPARISON TABLE (2D vs 3D vs Strip-Sum)
        %% ═══════════════════════════════════════════════════════════════

        function print_property_comparison(results, final_props)
            config = results.config;

            fprintf('\n┌──────────────────────────────────────────────────────────────────┐\n');
            fprintf('│  PROPERTY COMPARISON: 2D Surrogate / 3D Ground Truth / 3D Strips \n');
            fprintf('├──────────────────────────────────────────────────────────────────┤\n');

            % Recompute 2D properties at the final 3D solution
            x_final = results.stage2_3d.x_optimal;
            try
                p2d = calculate_2d_properties(x_final, config);
            catch
                p2d = [];
            end

            % Recompute 3D strip-summation (independent check)
            try
                strip_sum = WEC_Diagnostics.compute_strip_summation(x_final, config);
            catch
                strip_sum = [];
            end

            p3d = final_props;

            % Table header
            fprintf('│  %-22s %12s %12s %12s %8s\n', ...
                'Quantity', '2D Surr.', '3D Truth', '3D Strips', 'Err(%)');
            fprintf('│  %s\n', repmat('─', 1, 64));

            % Mass
            if ~isempty(p2d)
                m2 = p2d.mass_total;
            else
                m2 = NaN;
            end
            m3 = p3d.mass_total;
            if ~isempty(strip_sum)
                ms = strip_sum.mass;
            else
                ms = NaN;
            end
            WEC_Diagnostics.print_row('Mass [kg]', m2, m3, ms, '%.1f');

            % Buoyancy
            if ~isempty(p2d)
                b2 = p2d.mass_buoyant_force;
            else
                b2 = NaN;
            end
            WEC_Diagnostics.print_row('Buoyancy [kg]', b2, p3d.mass_buoyant_force, NaN, '%.1f');

            % V_sub
            if ~isempty(p2d)
                v2 = p2d.V_sub;
            else
                v2 = NaN;
            end
            WEC_Diagnostics.print_row('V_sub [m^3]', v2, p3d.V_sub, NaN, '%.4f');

            % GM
            if ~isempty(p2d)
                gm2 = p2d.GM;
                gm2_raw = p2d.GM_uncorrected;
            else
                gm2 = NaN;
                gm2_raw = NaN;
            end
            WEC_Diagnostics.print_row('GM [m]', gm2, p3d.GM_L, NaN, '%.4f');
            WEC_Diagnostics.print_row('GM_raw [m]', gm2_raw, p3d.GM_L, NaN, '%.4f');

            % CG_z
            if ~isempty(p2d)
                cg2 = p2d.CG_total(3);
            else
                cg2 = NaN;
            end
            if ~isempty(strip_sum)
                cg_s = strip_sum.CG_z;
            else
                cg_s = NaN;
            end
            WEC_Diagnostics.print_row('CG_z [m]', cg2, p3d.CG_total(3), cg_s, '%+.4f');

            % CB_z
            if ~isempty(p2d)
                cb2 = p2d.CB(3);
            else
                cb2 = NaN;
            end
            WEC_Diagnostics.print_row('CB_z [m]', cb2, p3d.CB(3), NaN, '%+.4f');

            % Iyy
            if ~isempty(p2d)
                iyy2 = p2d.Iyy;
            else
                iyy2 = NaN;
            end
            if ~isempty(strip_sum)
                iyy_s = strip_sum.Iyy;
            else
                iyy_s = NaN;
            end
            WEC_Diagnostics.print_row('Iyy [kg*m^2]', iyy2, p3d.Iyy, iyy_s, '%.1f');

            % Waterplane area
            if ~isempty(p2d)
                aw2 = p2d.Aw;
            else
                aw2 = NaN;
            end
            WEC_Diagnostics.print_row('Aw [m^2]', aw2, p3d.Aw, NaN, '%.4f');

            % Periods
            fprintf('│  %s\n', repmat('─', 1, 64));
            if ~isempty(p2d)
                th2 = p2d.periods.heave;
                tp2 = p2d.periods.pitch;
            else
                th2 = NaN;
                tp2 = NaN;
            end
            WEC_Diagnostics.print_row('T_heave [s]', th2, p3d.periods.heave, NaN, '%.3f');
            WEC_Diagnostics.print_row('T_pitch [s]', tp2, p3d.periods.pitch, NaN, '%.3f');

            if ~isinf(p3d.periods.surge) && ~isnan(p3d.periods.surge)
                if ~isempty(p2d)
                    ts2 = p2d.periods.surge;
                else
                    ts2 = NaN;
                end
                WEC_Diagnostics.print_row('T_surge [s]', ts2, p3d.periods.surge, NaN, '%.3f');
            else
                fprintf('│  %-22s %12s %12s %12s %8s\n', 'T_surge [s]', 'Inf', 'Inf', '—', '—');
            end

            % Added mass
            fprintf('│  %s\n', repmat('─', 1, 64));
            WEC_Diagnostics.print_row('A33 [kg]', NaN, p3d.A33, NaN, '%.1f');
            WEC_Diagnostics.print_row('A55 [kg*m^2]', NaN, p3d.A55, NaN, '%.1f');

            fprintf('└──────────────────────────────────────────────────────────────────┘\n');
        end


        %% ═══════════════════════════════════════════════════════════════
        %%  PANEL 4: STABILITY CHECK
        %% ═══════════════════════════════════════════════════════════════

        function print_stability_check(fp, config)
            fprintf('\n┌──────────────────────────────────────────────────────────────────┐\n');
            fprintf('│  STABILITY & PERIOD CHECK                                        \n');
            fprintf('├──────────────────────────────────────────────────────────────────┤\n');

            % GM check
            gm_ok = fp.GM_L >= config.gm_range(1) && fp.GM_L <= config.gm_range(2);
            gm_floor_ok = fp.GM_L >= config.gm_min;
            fprintf('│  GM_L = %.4f m    Target: %.2f m    Range: [%.2f, %.2f]\n', ...
                fp.GM_L, config.gm_target, config.gm_range(1), config.gm_range(2));
            fprintf('│    GM floor (%.2f m) : %s    GM in range : %s\n', ...
                config.gm_min, ...
                ternary(gm_floor_ok, 'PASS', 'FAIL'), ...
                ternary(gm_ok, 'PASS', 'WARN'));

            % CG below CB check
            cg_below_cb = fp.CG_total(3) < fp.CB(3);
            fprintf('│  CG_z = %+.4f m    CB_z = %+.4f m    CG < CB : %s\n', ...
                fp.CG_total(3), fp.CB(3), ternary(cg_below_cb, 'YES', 'NO'));

            fprintf('│  %s\n', repmat('─', 1, 64));

            % Period checks
            heave_ok = fp.periods.heave >= config.T_heave_range(1) && ...
                       fp.periods.heave <= config.T_heave_range(2);
            pitch_ok = fp.periods.pitch >= config.T_pitch_range(1) && ...
                       fp.periods.pitch <= config.T_pitch_range(2);

            fprintf('│  T_heave = %.3f s    Goal: %.1f s    Range: [%.1f, %.1f]    %s\n', ...
                fp.periods.heave, config.T_heave_goal, ...
                config.T_heave_range(1), config.T_heave_range(2), ...
                ternary(heave_ok, 'PASS', 'FAIL'));
            fprintf('│  T_pitch = %.3f s    Goal: %.1f s    Range: [%.1f, %.1f]    %s\n', ...
                fp.periods.pitch, config.T_pitch_goal, ...
                config.T_pitch_range(1), config.T_pitch_range(2), ...
                ternary(pitch_ok, 'PASS', 'FAIL'));

            % Mass balance
            mass_err_pct = 100 * abs(fp.mass_total - fp.mass_buoyant_force) / max(fp.mass_total, 1);
            fprintf('│  Mass balance: %.4e%%   %s\n', mass_err_pct, ...
                ternary(mass_err_pct < 0.1, 'PASS', 'FAIL'));

            fprintf('└──────────────────────────────────────────────────────────────────┘\n');
        end


        %% ═══════════════════════════════════════════════════════════════
        %%  PANEL 5: DENSITY PROFILE
        %% ═══════════════════════════════════════════════════════════════

        function print_density_profile(results, config)
            x_opt = results.stage2_3d.x_optimal;
            rho   = x_opt(2:end);
            N     = length(rho);
            nodes_z = config.density_nodes_z;

            fprintf('\n┌──────────────────────────────────────────────────────────────────┐\n');
            fprintf('│  DENSITY PROFILE (bottom → top)                                  \n');
            fprintf('├──────────────────────────────────────────────────────────────────┤\n');
            fprintf('│  %-6s %-10s %-10s %-8s %-8s %s\n', ...
                'Strip', 'z_node[m]', 'rho[kg/m3]', 'lb', 'ub', 'Note');
            fprintf('│  %s\n', repmat('─', 1, 58));

            for i = 1:N
                lb_i = config.ballast_density_bounds(1);
                ub_i = config.ballast_density_bounds(2);

                % Override with per-strip bounds if available
                if ~isempty(config.per_strip_density_lb)
                    lb_i = config.per_strip_density_lb(i);
                end
                if config.enable_constructability && ...
                        ~isempty(config.wall_strip_index) && ...
                        i == config.wall_strip_index
                    ub_i = config.constructability_rho_hull;
                end

                note = '';
                if config.enable_constructability && ...
                        ~isempty(config.wall_strip_index) && ...
                        i == config.wall_strip_index
                    note = 'WALL (pinned)';
                elseif abs(rho(i) - lb_i) < 1
                    note = 'at lower bound';
                elseif abs(rho(i) - ub_i) < 1
                    note = 'at upper bound';
                end

                fprintf('│  %-6d %+9.3f %10.1f %8.0f %8.0f  %s\n', ...
                    i, nodes_z(i), rho(i), lb_i, ub_i, note);
            end

            % Gradient check
            bot = mean(rho(1:min(3, N)));
            top = mean(rho(max(1, N-2):N));
            fprintf('│  %s\n', repmat('─', 1, 58));
            fprintf('│  Bottom-3 avg: %.0f kg/m^3    Top-3 avg: %.0f kg/m^3\n', bot, top);
            fprintf('│  Gradient: %s\n', ternary(bot > top, 'bottom-heavy (correct)', ...
                'TOP-HEAVY (inverted!)'));

            fprintf('└──────────────────────────────────────────────────────────────────┘\n');
        end


        %% ═══════════════════════════════════════════════════════════════
        %%  PANEL 6: CONSTRUCTABILITY
        %% ═══════════════════════════════════════════════════════════════

        function print_constructability_report(results, config)
            cstr = results.constructability;

            fprintf('\n┌──────────────────────────────────────────────────────────────────┐\n');
            fprintf('│  CONSTRUCTABILITY: UHPC + VOID REALIZATION                        \n');
            fprintf('├──────────────────────────────────────────────────────────────────┤\n');
            fprintf('│  Materials: UHPC = %.0f kg/m^3, void = %.1f kg/m^3\n', ...
                cstr.rho_hull, cstr.rho_fill);
            fprintf('│  t_min = %.1f mm,  wall height = %.2f m\n', ...
                cstr.t_min * 1000, cstr.wall_height);

            % Totals
            fprintf('│  %s\n', repmat('─', 1, 58));
            fprintf('│  Total V_hull  : %.4f m^3\n', cstr.total_V_hull);
            fprintf('│  Total V_UHPC  : %.4f m^3 (%.1f%%)\n', ...
                cstr.total_V_UHPC, 100 * cstr.total_V_UHPC / max(cstr.total_V_hull, eps));
            fprintf('│  Total V_void  : %.4f m^3 (%.1f%%)\n', ...
                cstr.total_V_void, 100 * cstr.total_V_void / max(cstr.total_V_hull, eps));
            fprintf('│  Total mass    : %.1f kg\n', cstr.total_mass);

            % Feasibility
            all_feasible = all(cstr.strip_is_feasible | cstr.strip_is_wall);
            n_infeasible = sum(~cstr.strip_is_feasible & ~cstr.strip_is_wall);

            fprintf('│  %s\n', repmat('─', 1, 58));
            fprintf('│  All strips feasible (t >= %.1f mm) : %s\n', ...
                cstr.t_min * 1000, ...
                ternary(all_feasible, 'YES', sprintf('NO (%d infeasible)', n_infeasible)));

            % Per-strip summary (compact)
            fprintf('│  %s\n', repmat('─', 1, 58));
            fprintf('│  %-5s %-8s %-8s %-7s %-8s %-7s %s\n', ...
                'Strip', 'rho_eff', 'scale', 't[mm]', 'V_UHPC', 'V_void', 'Status');

            N = length(cstr.strip_rho_eff);
            for i = 1:N
                if cstr.strip_is_wall(i)
                    status = 'WALL';
                elseif ~cstr.strip_is_feasible(i)
                    status = 'FAIL';
                else
                    status = 'OK';
                end
                fprintf('│  %-5d %7.0f %7.3f %7.1f %7.4f %7.4f  %s\n', ...
                    i, cstr.strip_rho_eff(i), cstr.strip_scale_factor(i), ...
                    cstr.strip_t_min_actual(i) * 1000, ...
                    cstr.strip_V_UHPC(i), cstr.strip_V_void(i), status);
            end

            % Verification against optimiser solution
            if isfield(cstr, 'verification')
                v = cstr.verification;
                fprintf('│  %s\n', repmat('─', 1, 58));
                fprintf('│  VERIFICATION (realised vs optimiser):\n');
                if isfield(v, 'err_mass_pct')
                    fprintf('│    Mass error   : %+.2f%%\n', v.err_mass_pct);
                end
                if isfield(v, 'err_CG_z_pct')
                    fprintf('│    CG_z error   : %+.2f%%\n', v.err_CG_z_pct);
                end
                if isfield(v, 'err_Iyy_pct')
                    fprintf('│    Iyy error    : %+.2f%%\n', v.err_Iyy_pct);
                end
                if isfield(v, 'err_GM_pct')
                    fprintf('│    GM error     : %+.2f%%\n', v.err_GM_pct);
                end
                if isfield(v, 'all_passed')
                    fprintf('│    All checks   : %s\n', ternary(v.all_passed, 'PASSED', 'FAILED'));
                end
            end

            fprintf('└──────────────────────────────────────────────────────────────────┘\n');
        end


        %% ═══════════════════════════════════════════════════════════════
        %%  PANEL 7: VERDICT
        %% ═══════════════════════════════════════════════════════════════

        function print_verdict(results, fp, config)
            fprintf('\n╔══════════════════════════════════════════════════════════════════╗\n');
            fprintf('║  VERDICT                                                         ║\n');
            fprintf('╠══════════════════════════════════════════════════════════════════╣\n');

            issues = {};

            % Check 1: Stage 1 convergence
            if ~results.stage1_2d.converged
                issues{end+1} = 'Stage 1 (2D surrogate) did NOT converge';
            end
            if results.stage1_2d.pid_saturated
                issues{end+1} = 'PID correction factors hit saturation bounds';
            end

            % Check 2: Stage 2 convergence
            if results.stage2_3d.exitflag <= 0
                issues{end+1} = sprintf('Stage 2 (3D fmincon) exit flag = %d (not converged)', ...
                    results.stage2_3d.exitflag);
            end

            % Check 3: Mass balance
            mass_err_pct = 100 * abs(fp.mass_total - fp.mass_buoyant_force) / max(fp.mass_total, 1);
            if mass_err_pct > 0.1
                issues{end+1} = sprintf('Mass balance error = %.2e%% (>0.1%%)', mass_err_pct);
            end

            % Check 4: GM
            if fp.GM_L < config.gm_min
                issues{end+1} = sprintf('GM_L = %.4f m < gm_min = %.2f m', fp.GM_L, config.gm_min);
            end
            if fp.GM_L < 0
                issues{end+1} = sprintf('GM_L = %.4f m (NEGATIVE — capsizes)', fp.GM_L);
            end

            % Check 5: Period ranges
            if fp.periods.heave < config.T_heave_range(1) || fp.periods.heave > config.T_heave_range(2)
                issues{end+1} = sprintf('T_heave = %.2f s outside [%.1f, %.1f]', ...
                    fp.periods.heave, config.T_heave_range(1), config.T_heave_range(2));
            end
            if fp.periods.pitch < config.T_pitch_range(1) || fp.periods.pitch > config.T_pitch_range(2)
                issues{end+1} = sprintf('T_pitch = %.2f s outside [%.1f, %.1f]', ...
                    fp.periods.pitch, config.T_pitch_range(1), config.T_pitch_range(2));
            end

            % Check 6: Constructability
            if config.enable_constructability && ...
                    isfield(results, 'constructability') && ...
                    ~isempty(results.constructability)
                cstr = results.constructability;
                n_infeasible = sum(~cstr.strip_is_feasible & ~cstr.strip_is_wall);
                if n_infeasible > 0
                    issues{end+1} = sprintf('%d platform strip(s) violate t_min', n_infeasible);
                end
            end

            % Print verdict
            if isempty(issues)
                fprintf('║                                                                  ║\n');
                fprintf('║    ALL CHECKS PASSED                                             ║\n');
                fprintf('║                                                                  ║\n');
                fprintf('║    The solution is converged, stable, mass-balanced, and          ║\n');
                fprintf('║    within period targets.                                        ║\n');
                if config.enable_constructability
                fprintf('║    Constructability: all strips feasible.                        ║\n');
                end
            else
                fprintf('║                                                                  ║\n');
                fprintf('║    %d ISSUE(S) FOUND:                                            ║\n', length(issues));
                fprintf('║                                                                  ║\n');
                for k = 1:length(issues)
                    line = sprintf('    %d. %s', k, issues{k});
                    fprintf('║  %-64s║\n', line);
                end
            end

            fprintf('║                                                                  ║\n');
            fprintf('║  Total time: %.2f s                                              ║\n', ...
                results.optimization_time);
            fprintf('╚══════════════════════════════════════════════════════════════════╝\n\n');
        end


        %% ═══════════════════════════════════════════════════════════════
        %%  HELPERS
        %% ═══════════════════════════════════════════════════════════════

        function print_row(label, val_2d, val_3d, val_strip, fmt)
        % PRINT_ROW  One row of the comparison table with percentage error.
            fmt_val = fmt;
            err_str = '—';

            if isfinite(val_2d) && isfinite(val_3d) && abs(val_3d) > 1e-12
                err = 100 * (val_2d - val_3d) / val_3d;
                err_str = sprintf('%+.1f', err);
            end

            s2d = WEC_Diagnostics.fmt_or_dash(val_2d, fmt_val);
            s3d = WEC_Diagnostics.fmt_or_dash(val_3d, fmt_val);
            sst = WEC_Diagnostics.fmt_or_dash(val_strip, fmt_val);

            fprintf('│  %-22s %12s %12s %12s %8s\n', label, s2d, s3d, sst, err_str);
        end


        function s = fmt_or_dash(val, fmt)
        % FMT_OR_DASH  Format a value or return '—' if NaN.
            if isnan(val) || isinf(val)
                s = '—';
            else
                s = sprintf(fmt, val);
            end
        end


        function strip_sum = compute_strip_summation(x, config)
        % COMPUTE_STRIP_SUMMATION  Independent strip-by-strip recomputation
        %   of mass, CG, and Iyy from the 3D mesh face data.
        %
        %   This is a THIRD independent computation path (distinct from both
        %   calculate_2d_properties and calculate_3d_properties) that
        %   partitions the divergence-theorem face contributions by z-bin
        %   and sums them.  If it disagrees with calculate_3d_properties,
        %   there is an inconsistency in the density assignment logic.

            densities = x(2:end);
            densities = densities(:);
            face_z   = config.topology.face_centroids(:, 3);
            face_vol = config.topology.face_volume_contribs;

            N = config.num_ballast_sections;

            % Assign density per face
            if config.enable_constructability && ~isempty(config.strip_edges)
                bins = discretize(face_z, config.strip_edges);
                bins(isnan(bins) & face_z <= config.strip_edges(1))   = 1;
                bins(isnan(bins) & face_z >= config.strip_edges(end)) = N;
                bins(isnan(bins)) = 1;
                face_rho = densities(bins);
            else
                face_rho = interp1(config.density_nodes_z, densities, ...
                    face_z, 'linear', 'extrap');
            end
            face_rho = max(config.ballast_density_bounds(1), ...
                       min(config.ballast_density_bounds(2), face_rho));

            face_mass = face_rho .* face_vol;
            strip_sum.mass = sum(face_mass);

            if abs(strip_sum.mass) > 1e-6
                cg_z_num = sum(face_z .* face_mass);
                strip_sum.CG_z = cg_z_num / strip_sum.mass + x(1);  % add draft shift
            else
                strip_sum.CG_z = 0;
            end

            % Iyy about CG
            cg_orig = [0, 0, strip_sum.CG_z - x(1)];  % in original coords
            face_cx = config.topology.face_centroids(:, 1);
            d_x = face_cx - cg_orig(1);
            d_z = face_z  - cg_orig(3);
            strip_sum.Iyy = sum(face_mass .* (d_x.^2 + d_z.^2));
        end

    end  % methods (Static)
end  % classdef


%% ═══════════════════════════════════════════════════════════════════
%%  MODULE-LEVEL HELPER
%% ═══════════════════════════════════════════════════════════════════

function r = ternary(cond, t, f)
    if cond, r = t; else, r = f; end
end