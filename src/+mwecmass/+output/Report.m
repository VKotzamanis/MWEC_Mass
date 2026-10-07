classdef Report
%REPORT Stateless console panels for one mass-distribution run.
%   Static methods read result/config structs and write formatted text to fids (default 1).
    methods (Static)
        function header(config, fids)
            if nargin < 2, fids = 1; end
            mwecmass.output.emit(fids, '\n');
            mwecmass.output.emit(fids, '╔══════════════════════════════════════════════════════════════════╗\n');
            mwecmass.output.emit(fids, '║              WEC PIPELINE — DIAGNOSTIC REPORT                   ║\n');
            mwecmass.output.emit(fids, '╠══════════════════════════════════════════════════════════════════╣\n');
            if isfield(config, 'ms2_file')
                geom_str = config.ms2_file;
            else
                geom_str = '(not stored in config)';
            end
            mwecmass.output.emit(fids, '║  MS2 file       : %-45s ║\n', geom_str);
            mwecmass.output.emit(fids, '║  Mode           : %-45s ║\n', config.stage1_mode);
            mwecmass.output.emit(fids, '║  Strips         : %-45s ║\n', ...
                sprintf('%d  (constructability: %s)', ...
                    config.num_ballast_sections, ...
                    mwecmass.internal.ternary(config.enable_constructability, 'ON', 'OFF')));
            mwecmass.output.emit(fids, '║  Targets        : T_h=%.1f s, T_p=%.1f s, GM=%.2f m %18s║\n', ...
                config.T_heave_goal, config.T_pitch_goal, config.gm_target, '');
            mwecmass.output.emit(fids, '╚══════════════════════════════════════════════════════════════════╝\n');
        end

        function stage1(results, fids)
            if nargin < 2, fids = 1; end
            s1 = results.stage1_2d;
            config = results.config;

            mwecmass.output.emit(fids, '\n┌──────────────────────────────────────────────────────────────────┐\n');
            mwecmass.output.emit(fids, '│  STAGE 1: 2D SURROGATE  (%s)                              \n', upper(config.stage1_mode));
            mwecmass.output.emit(fids, '├──────────────────────────────────────────────────────────────────┤\n');

            mwecmass.output.emit(fids, '│  Iterations : %-4d    Converged: %-3s    PID saturated: %-3s\n', ...
                s1.iterations, ...
                mwecmass.internal.ternary(s1.converged, 'YES', 'NO'), ...
                mwecmass.internal.ternary(s1.pid_saturated, 'YES', 'no'));

            if ~isempty(s1.mass_corrections)
                k_vol_final = s1.mass_corrections(end);
                k_gm_final  = s1.gm_corrections(end);
                mwecmass.output.emit(fids, '│  k_vol final: %.4f    k_gm final: %.4f\n', k_vol_final, k_gm_final);
            end

            if ~isempty(s1.mass_errors) && ~isempty(s1.gm_errors)
                mass_gap_kg = s1.mass_errors(end);
                gm_gap_m    = s1.gm_errors(end);

                if ~isempty(s1.mass_3d_history) && s1.mass_3d_history(end) > 0
                    mass_gap_pct = 100 * mass_gap_kg / s1.mass_3d_history(end);
                else
                    mass_gap_pct = NaN;
                end

                mwecmass.output.emit(fids, '│  Mass gap   : %+.1f kg  (%+.1f%%)\n', mass_gap_kg, mass_gap_pct);
                mwecmass.output.emit(fids, '│  GM gap     : %+.4f m\n', gm_gap_m);
            end

            if ~isempty(s1.R2_mass) && ~isnan(s1.R2_mass(end))
                mwecmass.output.emit(fids, '│  R²(mass)   : %.3f    R²(GM): %.3f    R²(CG_z): %.3f\n', ...
                    s1.R2_mass(end), s1.R2_GM(end), s1.R2_cg(end));
            end

            if isfield(s1, 'convergence_data') && ...
                    isstruct(s1.convergence_data) && ...
                    isfield(s1.convergence_data, 'final_exitflag')
                mwecmass.output.emit(fids, '│  fmincon exit: %d\n', s1.convergence_data.final_exitflag);
            end

            mwecmass.output.emit(fids, '└──────────────────────────────────────────────────────────────────┘\n');
        end

        function stage2(results, final_props, fids)
            if nargin < 3, fids = 1; end
            s2 = results.stage2_3d;
            config = results.config;

            mwecmass.output.emit(fids, '\n┌──────────────────────────────────────────────────────────────────┐\n');
            mwecmass.output.emit(fids, '│  STAGE 2: 3D HIGH-FIDELITY REFINEMENT                           \n');
            mwecmass.output.emit(fids, '├──────────────────────────────────────────────────────────────────┤\n');
            mwecmass.output.emit(fids, '│  Exit flag  : %-4d    Converged: %-3s\n', ...
                s2.exitflag, mwecmass.internal.ternary(s2.converged, 'YES', 'NO'));
            mwecmass.output.emit(fids, '│  Cost (fval): %.6f\n', s2.fval);

            % funcCount may not exist in older/partial fmincon output; degrade gracefully instead of throwing.
            if isfield(s2, 'output') && isfield(s2.output, 'iterations')
                if isfield(s2.output, 'funcCount')
                    mwecmass.output.emit(fids, '│  SQP iters  : %d    FunEvals: %d\n', ...
                        s2.output.iterations, s2.output.funcCount);
                else
                    mwecmass.output.emit(fids, '│  SQP iters  : %d\n', s2.output.iterations);
                end
            end

            mass_err_pct = 100 * abs(final_props.mass_total - final_props.mass_buoyant_force) ...
                            / max(final_props.mass_total, 1);
            mwecmass.output.emit(fids, '│  Mass error : %.2e%%  %s\n', mass_err_pct, ...
                mwecmass.internal.ternary(mass_err_pct < 0.1, '(OK)', '** VIOLATED **'));

            if isfield(final_props, 'draft') && final_props.draft > 0
                draft_true = final_props.draft;
            elseif isfield(config, 'hull_z_min')
                draft_true = abs(config.hull_z_min + final_props.vertical_shift);
            else
                draft_true = 0;
            end
            mwecmass.output.emit(fids, '│  Draft      : %.4f m  (shift: %.4f m)\n', draft_true, final_props.vertical_shift);
            mwecmass.output.emit(fids, '└──────────────────────────────────────────────────────────────────┘\n');
        end

        function stability(fp, config, fids)
            if nargin < 3, fids = 1; end
            mwecmass.output.emit(fids, '\n┌──────────────────────────────────────────────────────────────────┐\n');
            mwecmass.output.emit(fids, '│  STABILITY & PERIOD CHECK                                        \n');
            mwecmass.output.emit(fids, '├──────────────────────────────────────────────────────────────────┤\n');

            gm_ok = fp.GM_L >= config.gm_range(1) && fp.GM_L <= config.gm_range(2);
            gm_floor_ok = fp.GM_L >= config.gm_min;
            mwecmass.output.emit(fids, '│  GM_L = %.4f m    Target: %.2f m    Range: [%.2f, %.2f]\n', ...
                fp.GM_L, config.gm_target, config.gm_range(1), config.gm_range(2));
            mwecmass.output.emit(fids, '│    GM floor (%.2f m) : %s    GM in range : %s\n', ...
                config.gm_min, ...
                mwecmass.internal.ternary(gm_floor_ok, 'PASS', 'FAIL'), ...
                mwecmass.internal.ternary(gm_ok, 'PASS', 'WARN'));

            cg_below_cb = fp.CG_total(3) < fp.CB(3);
            mwecmass.output.emit(fids, '│  CG_z = %+.4f m    CB_z = %+.4f m    CG < CB : %s\n', ...
                fp.CG_total(3), fp.CB(3), mwecmass.internal.ternary(cg_below_cb, 'YES', 'NO'));

            mwecmass.output.emit(fids, '│  %s\n', repmat('─', 1, 64));

            heave_ok = fp.periods.heave >= config.T_heave_range(1) && ...
                       fp.periods.heave <= config.T_heave_range(2);
            pitch_ok = fp.periods.pitch >= config.T_pitch_range(1) && ...
                       fp.periods.pitch <= config.T_pitch_range(2);

            mwecmass.output.emit(fids, '│  T_heave = %.3f s    Goal: %.1f s    Range: [%.1f, %.1f]    %s\n', ...
                fp.periods.heave, config.T_heave_goal, ...
                config.T_heave_range(1), config.T_heave_range(2), ...
                mwecmass.internal.ternary(heave_ok, 'PASS', 'FAIL'));
            mwecmass.output.emit(fids, '│  T_pitch = %.3f s    Goal: %.1f s    Range: [%.1f, %.1f]    %s\n', ...
                fp.periods.pitch, config.T_pitch_goal, ...
                config.T_pitch_range(1), config.T_pitch_range(2), ...
                mwecmass.internal.ternary(pitch_ok, 'PASS', 'FAIL'));

            mass_err_pct = 100 * abs(fp.mass_total - fp.mass_buoyant_force) / max(fp.mass_total, 1);
            mwecmass.output.emit(fids, '│  Mass balance: %.4e%%   %s\n', mass_err_pct, ...
                mwecmass.internal.ternary(mass_err_pct < 0.1, 'PASS', 'FAIL'));

            mwecmass.output.emit(fids, '└──────────────────────────────────────────────────────────────────┘\n');
        end

        function density_profile(results, config, fids)
            if nargin < 3, fids = 1; end
            x_opt = results.stage2_3d.x_optimal;
            rho   = x_opt(2:end);
            N     = length(rho);
            nodes_z = config.density_nodes_z;

            mwecmass.output.emit(fids, '\n┌──────────────────────────────────────────────────────────────────┐\n');
            mwecmass.output.emit(fids, '│  DENSITY PROFILE (bottom → top)                                  \n');
            mwecmass.output.emit(fids, '├──────────────────────────────────────────────────────────────────┤\n');
            mwecmass.output.emit(fids, '│  %-6s %-10s %-10s %-8s %-8s %s\n', ...
                'Strip', 'z_node[m]', 'rho[kg/m3]', 'lb', 'ub', 'Note');
            mwecmass.output.emit(fids, '│  %s\n', repmat('─', 1, 58));

            for i = 1:N
                lb_i = config.ballast_density_bounds(1);
                ub_i = config.ballast_density_bounds(2);

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

                mwecmass.output.emit(fids, '│  %-6d %+9.3f %10.1f %8.0f %8.0f  %s\n', ...
                    i, nodes_z(i), rho(i), lb_i, ub_i, note);
            end

            bot = mean(rho(1:min(3, N)));
            top = mean(rho(max(1, N-2):N));
            mwecmass.output.emit(fids, '│  %s\n', repmat('─', 1, 58));
            mwecmass.output.emit(fids, '│  Bottom-3 avg: %.0f kg/m^3    Top-3 avg: %.0f kg/m^3\n', bot, top);
            mwecmass.output.emit(fids, '│  Gradient: %s\n', mwecmass.internal.ternary(bot > top, 'bottom-heavy (correct)', ...
                'TOP-HEAVY (inverted!)'));

            mwecmass.output.emit(fids, '└──────────────────────────────────────────────────────────────────┘\n');
        end

        function stage3(results, fids)
            if nargin < 2, fids = 1; end
            r = results.stage3;

            mwecmass.output.emit(fids, '\n┌──────────────────────────────────────────────────────────────────┐\n');
            mwecmass.output.emit(fids, '│  STAGE 3: %s REALISATION (%s)\n', upper(strrep(r.mode, '_', ' ')), upper(r.status));
            mwecmass.output.emit(fids, '├──────────────────────────────────────────────────────────────────┤\n');
            mwecmass.output.emit(fids, '│  Last step: %s    draft = %.4f m    z_ballast = %.4f m (body)\n', ...
                r.escalation, r.draft, r.design.z_ballast);
            densities = fieldnames(r.rho);
            for k = 1:numel(densities)
                mwecmass.output.emit(fids, '│  rho_%s = %.1f kg/m^3\n', densities{k}, r.rho.(densities{k}));
            end
            mwecmass.output.emit(fids, '│  %s\n', repmat('─', 1, 58));
            mwecmass.output.emit(fids, '│  %-4s %8s %9s %10s %9s %9s %9s\n', ...
                'Mod', 't[mm]', 'h_bal[m]', 'V[m^3]', 'rho_eff', 'rho_S2', 'rho_min');
            for i = 1:numel(r.modules)
                m = r.modules(i);
                mwecmass.output.emit(fids, '│  %-4d %8.2f %9.4f %10.5f %9.1f %9.1f %9.1f\n', ...
                    i, 1000 * m.t, m.h_ballast, m.V, m.rho_eff, m.rho_stage2, m.rho_floor);
            end
            mwecmass.output.emit(fids, '│  %s\n', repmat('─', 1, 58));
            mwecmass.output.emit(fids, '│  %-8s %12s %12s %9s %8s %s\n', 'Metric', 'Stage 2', 'Stage 3', 'dev[%]', 'lim[%]', 'Pass');
            for m = r.check.metrics
                mwecmass.output.emit(fids, '│  %-8s %12.5f %12.5f %9.3f %8.2f %s\n', m.name, m.stage2, m.value, ...
                    100 * m.rel_dev, 100 * m.limit, mwecmass.internal.ternary(m.pass, 'yes', 'NO'));
            end
            for q = r.check.equalities
                mwecmass.output.emit(fids, '│  Equality %-9s residual %10.3g (tol %.1g) %s\n', q.name, q.residual, q.tol, ...
                    mwecmass.internal.ternary(q.pass, 'holds', 'NOT MET'));
            end
            if ~isempty(r.reason)
                mwecmass.output.emit(fids, '│  Reason: %s\n', r.reason);
            end
            mwecmass.output.emit(fids, '└──────────────────────────────────────────────────────────────────┘\n');
        end

        function verdict(results, fp, config, fids)
            if nargin < 4, fids = 1; end
            mwecmass.output.emit(fids, '\n╔══════════════════════════════════════════════════════════════════╗\n');
            mwecmass.output.emit(fids, '║  VERDICT                                                         ║\n');
            mwecmass.output.emit(fids, '╠══════════════════════════════════════════════════════════════════╣\n');

            issues = {};

            if ~results.stage1_2d.converged
                issues{end+1} = 'Stage 1 (2D surrogate) did NOT converge';
            end
            if results.stage1_2d.pid_saturated
                issues{end+1} = 'PID correction factors hit saturation bounds';
            end

            if results.stage2_3d.exitflag <= 0
                issues{end+1} = sprintf('Stage 2 (3D fmincon) exit flag = %d (not converged)', ...
                    results.stage2_3d.exitflag);
            end

            mass_err_pct = 100 * abs(fp.mass_total - fp.mass_buoyant_force) / max(fp.mass_total, 1);
            if mass_err_pct > 0.1
                issues{end+1} = sprintf('Mass balance error = %.2e%% (>0.1%%)', mass_err_pct);
            end

            if fp.GM_L < config.gm_min
                issues{end+1} = sprintf('GM_L = %.4f m < gm_min = %.2f m', fp.GM_L, config.gm_min);
            end
            if fp.GM_L < 0
                issues{end+1} = sprintf('GM_L = %.4f m (NEGATIVE — capsizes)', fp.GM_L);
            end

            if fp.periods.heave < config.T_heave_range(1) || fp.periods.heave > config.T_heave_range(2)
                issues{end+1} = sprintf('T_heave = %.2f s outside [%.1f, %.1f]', ...
                    fp.periods.heave, config.T_heave_range(1), config.T_heave_range(2));
            end
            if fp.periods.pitch < config.T_pitch_range(1) || fp.periods.pitch > config.T_pitch_range(2)
                issues{end+1} = sprintf('T_pitch = %.2f s outside [%.1f, %.1f]', ...
                    fp.periods.pitch, config.T_pitch_range(1), config.T_pitch_range(2));
            end

            if isfield(results, 'stage3') && ~isempty(results.stage3) && ...
                    strcmp(results.stage3.status, 'failed')
                issues{end+1} = sprintf('Stage 3 failed: %s', strjoin(results.stage3.check.failed, ', '));
            end

            if isempty(issues)
                mwecmass.output.emit(fids, '║                                                                  ║\n');
                mwecmass.output.emit(fids, '║    ALL CHECKS PASSED                                             ║\n');
                mwecmass.output.emit(fids, '║                                                                  ║\n');
                mwecmass.output.emit(fids, '║    The solution is converged, stable, mass-balanced, and          ║\n');
                mwecmass.output.emit(fids, '║    within period targets.                                        ║\n');
                if isfield(results, 'stage3') && ~isempty(results.stage3)
                mwecmass.output.emit(fids, '║    Stage 3: realised design accepted against Stage 2.            ║\n');
                end
            else
                mwecmass.output.emit(fids, '║                                                                  ║\n');
                mwecmass.output.emit(fids, '║    %d ISSUE(S) FOUND:                                            ║\n', length(issues));
                mwecmass.output.emit(fids, '║                                                                  ║\n');
                for k = 1:length(issues)
                    line = sprintf('    %d. %s', k, issues{k});
                    mwecmass.output.emit(fids, '║  %-64s║\n', line);
                end
            end

            mwecmass.output.emit(fids, '║                                                                  ║\n');
            mwecmass.output.emit(fids, '║  Total time: %.2f s                                              ║\n', ...
                results.optimization_time);
            mwecmass.output.emit(fids, '╚══════════════════════════════════════════════════════════════════╝\n\n');
        end

        function exception(ME, ms2_file, hydro_table_file, fids)
            if nargin < 4, fids = 1; end
            mwecmass.output.emit(fids, '\n  OPTIMISATION FAILED: %s\n', ME.message);
            if ~isempty(ME.stack)
                mwecmass.output.emit(fids, '    in %s, line %d\n', ME.stack(1).name, ME.stack(1).line);
            end
            if ~exist(ms2_file, 'file')
                mwecmass.output.emit(fids, '  MS2 file not found: %s\n', ms2_file);
            end
            if ~isempty(hydro_table_file) && ~exist(hydro_table_file, 'file')
                mwecmass.output.emit(fids, '  HAMS hydro_table not found: %s\n', hydro_table_file);
            end
        end

        function matrix(M, prefix, fids)
            if nargin < 3, fids = 1; end
            for ii = 1:3
                mwecmass.output.emit(fids, '%s[', prefix);
                for jj = 1:3
                    mwecmass.output.emit(fids, '%12.3f', M(ii,jj));
                    if jj < 3
                        mwecmass.output.emit(fids, ', ');
                    end
                end
                mwecmass.output.emit(fids, ']\n');
            end
        end

        function final_results(props, config, fids)
            if nargin < 3, fids = 1; end
            try
                mwecmass.output.emit(fids, '\n=== FINAL 3D DESIGN RESULTS ===\n\n');
                mwecmass.output.emit(fids, '  Optimal Draft: %.3f m\n', abs(config.hull_z_min + props.vertical_shift));
                mwecmass.output.emit(fids, '  Vertical Shift: %.3f m\n', props.vertical_shift);
                if isfield(props, 'components') && ~isempty(props.components)
                    mwecmass.output.emit(fids, '\n  Density Profile:\n');
                    for i = 1:length(props.components)
                        mwecmass.output.emit(fids, '    Node %d @ Z = %6.2f m: rho = %7.0f kg/m^3\n', ...
                            i, props.components(i).z_level, props.components(i).density);
                    end
                end
                mwecmass.output.emit(fids, '\n  Mass Properties:\n');
                mwecmass.output.emit(fids, '    Total Mass:       %10.2f kg\n', props.mass_total);
                mwecmass.output.emit(fids, '    Displaced Mass:   %10.2f kg\n', props.mass_buoyant_force);
                mwecmass.output.emit(fids, '    Mass Discrepancy: %10.2e kg\n', props.mass_discrepancy);
                mwecmass.output.emit(fids, '    CG: [%.3f, %.3f, %.3f] m\n', props.CG_total);
                mwecmass.output.emit(fids, '\n  Hydrostatic Properties:\n');
                mwecmass.output.emit(fids, '    Submerged Volume: %.4f m^3\n', props.V_sub);
                mwecmass.output.emit(fids, '    Waterplane Area:  %.4f m^2\n', props.Aw);
                mwecmass.output.emit(fids, '    CB: [%.3f, %.3f, %.3f] m\n', props.CB);
                mwecmass.output.emit(fids, '    GM: %.3f m\n', props.GM_L);
                mwecmass.output.emit(fids, '\n  Natural Periods:\n');
                mwecmass.output.emit(fids, '    Heave: %.2f s (Target: %.2f s)\n', props.periods.heave, config.T_heave_goal);
                mwecmass.output.emit(fids, '    Pitch: %.2f s (Target: %.2f s)\n', props.periods.pitch, config.T_pitch_goal);
                if isinf(props.periods.surge)
                    mwecmass.output.emit(fids, '    Surge: Inf s\n');
                else
                    mwecmass.output.emit(fids, '    Surge: %.2f s\n', props.periods.surge);
                end
                mwecmass.output.emit(fids, '\n');
            catch ME
                warning('mwecmass:output:ReportFailed', 'Results reporting failed: %s', ME.message);
            end
        end

        function all_panels(results, final_props, fids)
            if nargin < 3, fids = 1; end
            config = results.config;

            try mwecmass.output.Report.header(config, fids);
            catch ME, mwecmass.output.emit(fids, '  [Header panel failed: %s]\n', ME.message); end

            try mwecmass.output.Report.stage1(results, fids);
            catch ME, mwecmass.output.emit(fids, '  [Stage 1 panel failed: %s]\n', ME.message); end

            try mwecmass.output.Report.stage2(results, final_props, fids);
            catch ME, mwecmass.output.emit(fids, '  [Stage 2 panel failed: %s]\n', ME.message); end

            try mwecmass.output.print_property_comparison(results, final_props, fids);
            catch ME, mwecmass.output.emit(fids, '  [Property comparison panel failed: %s]\n', ME.message); end

            try mwecmass.output.Report.stability(final_props, config, fids);
            catch ME, mwecmass.output.emit(fids, '  [Stability panel failed: %s]\n', ME.message); end

            try mwecmass.output.Report.density_profile(results, config, fids);
            catch ME, mwecmass.output.emit(fids, '  [Density profile panel failed: %s]\n', ME.message); end

            if isfield(results, 'stage3') && ~isempty(results.stage3)
                try mwecmass.output.Report.stage3(results, fids);
                catch ME, mwecmass.output.emit(fids, '  [Stage 3 panel failed: %s]\n', ME.message); end
            end

            try mwecmass.output.Report.verdict(results, final_props, config, fids);
            catch ME, mwecmass.output.emit(fids, '  [Verdict panel failed: %s]\n', ME.message); end
        end
    end
end
