function hydro_table = generate_draft_nodes(in, mesh_sizing, geo_config, ms2_file, hams_dir, hams_exe, cache_file)
%GENERATE_DRAFT_NODES  Run coarse and adaptive HAMS sweeps over vertical-shift nodes.
% in, mesh_sizing, and geo_config come from run/build_config; paths identify the deck, solver,
% and cache.  vertical_shift is signed body-frame [m] (z_wl = -vertical_shift), and T_heave is [s].
% The output hydro_table follows the 23-field cache contract and is unsorted; load_hydro_cache
% sorts it.  mesh_Nu/mesh_Nv are counts and wp_target_edge is [m].  Refinement adds midpoint
% runs in intervals with the largest heave-period gradient, bounded by in.bem.n_adaptive_refine.
% See docs/HAMS_MREL_ROUTE.md#function-map for draft-node conventions.
        vs_bounds = geo_config.vertical_shift_bounds;

        % n_coarse comes from the input file, not a hardcoded literal.
        n_coarse  = in.bem.n_sweep_drafts;
        vs_coarse = linspace(vs_bounds(1), vs_bounds(2), n_coarse);

        hams_opts = struct('wp_target_edge', mesh_sizing.wp_target_edge);

        fprintf('\n  Running HAMS at %d drafts: [', n_coarse);
        fprintf('%.2f ', vs_coarse);
        fprintf(']\n');

        hydro_cache = mwecmass.bem.empty_hydro_cache();
        hydro_cache.ms2_file = ms2_file;
        if exist(ms2_file, 'file')
            ms2_info = dir(ms2_file);
            hydro_cache.ms2_date = ms2_info.date;
        end

        for k = 1:n_coarse
            fprintf('\n  [%d/%d] vertical_shift = %+.4f m\n', k, n_coarse, vs_coarse(k));
            [~, hydro_cache] = mwecmass.bem.get_or_run_hydro( ...
                vs_coarse(k), hydro_cache, geo_config, hams_dir, hams_exe, 0.01, hams_opts);
        end

        % Stash mesh sizing so [L]oad can restore Nu/Nv (and the matched
        % WP lid edge) for the post-run diagnostic mesh without
        % re-parsing the .ms2 file.  Also stash the period grid so [L]oad
        % can warn if the user changed it (cache becomes stale).
        hydro_cache.mesh_Nu        = mesh_sizing.mesh_Nu;
        hydro_cache.mesh_Nv        = mesh_sizing.mesh_Nv;
        hydro_cache.panel_size     = in.geometry.panel_size;
        hydro_cache.wp_target_edge = mesh_sizing.wp_target_edge;
        hydro_cache.period_min     = in.bem.T_min;
        hydro_cache.period_max     = in.bem.T_max;
        hydro_cache.period_step    = in.bem.T_step;

        hydro_table = hydro_cache;
        save(cache_file, 'hydro_table', '-v7.3');
        fprintf('\n  Hydro cache saved: %s (%d entries)\n', cache_file, ...
                length(hydro_cache.drafts));

        % After the coarse sweep, refine intervals with the largest
        % |dT_heave/dvs| by adding midpoint HAMS runs.

        if in.bem.n_adaptive_refine > 0 && length(hydro_cache.drafts) >= 2

            fprintf('\n  ── Adaptive T_heave refinement (up to %d extra run(s)) ──\n', ...
                    in.bem.n_adaptive_refine);

            % Compute T_heave at every cached draft.
            n_cached  = length(hydro_cache.drafts);
            vs_all    = hydro_cache.drafts(:);        % vertical shift [m]
            T_h_all   = zeros(n_cached, 1);

            for kk = 1:n_cached
                vs_kk   = vs_all(kk);
                z_wl_kk = -vs_kk;   % waterline in body frame

                % Heave added mass at infinite frequency (DOF 3 = heave)
                A33_kk = hydro_cache.added_mass_inf{kk}(3, 3);

                % Aw and V_sub from geometry tables (shared z-axis)
                Aw_kk   = max(0, interp1(geo_config.Aw_table_z, ...
                                         geo_config.Aw_table, ...
                                         z_wl_kk, 'linear', 0));
                Vsub_kk = max(0, interp1(geo_config.Aw_table_z, ...
                                         geo_config.V_sub_table, ...
                                         z_wl_kk, 'linear', 0));
                M_kk = Vsub_kk * in.constants.rho_water;   % hydrostatic mass

                if Aw_kk > 1e-6 && (M_kk + A33_kk) > 0
                    K33_kk      = in.constants.rho_water * in.constants.g * Aw_kk;
                    T_h_all(kk) = 2*pi * sqrt((M_kk + A33_kk) / K33_kk);
                else
                    T_h_all(kk) = Inf;   % dry or degenerate draft
                end
            end

            % Sort by vertical shift before computing interval gradients
            [vs_srt, srt_idx] = sort(vs_all);
            T_h_srt = T_h_all(srt_idx);

            fprintf('  Coarse T_heave (sorted by draft):\n');
            fprintf('    vs=%+.3f m → T_h=%.2f s\n', ...
                    [vs_srt, T_h_srt]');

            % Compute the interval gradients.
            n_iv   = length(vs_srt) - 1;
            grad_T = zeros(n_iv, 1);
            for kk = 1:n_iv
                dvs = vs_srt(kk+1) - vs_srt(kk);
                dT  = T_h_srt(kk+1) - T_h_srt(kk);
                if dvs > 1e-6 && isfinite(dT)
                    grad_T(kk) = abs(dT / dvs);
                end
            end

            fprintf('  |dT/dvs| per interval: ');
            fprintf('%.1f ', grad_T);
            fprintf('s/m\n');

            % Bisect the highest-gradient intervals.
            [~, grad_order] = sort(grad_T, 'descend');
            n_done = 0;

            for ki = 1:n_iv
                if n_done >= in.bem.n_adaptive_refine; break; end

                iv     = grad_order(ki);
                g_iv   = grad_T(iv);

                % Threshold: skip intervals with negligible gradient
                % (< 1 s/m → adding a midpoint changes T_heave by < 0.1 s
                %  for a typical 0.1 m bisection step — not worth a HAMS run)
                if g_iv < 1.0; break; end

                vs_mid = 0.5 * (vs_srt(iv) + vs_srt(iv+1));

                % Skip if a point already exists within tolerance
                if any(abs(hydro_cache.drafts - vs_mid) < 0.01)
                    fprintf('  [skip] vs=%+.4f m already in cache.\n', vs_mid);
                    continue;
                end

                fprintf('\n  [refine %d/%d]  vs=%+.4f m  (|dT/dvs|=%.1f s/m,', ...
                        n_done+1, in.bem.n_adaptive_refine, vs_mid, g_iv);
                fprintf('  interval [%+.3f, %+.3f] m)\n', ...
                        vs_srt(iv), vs_srt(iv+1));

                [~, hydro_cache] = mwecmass.bem.get_or_run_hydro( ...
                    vs_mid, hydro_cache, geo_config, ...
                    hams_dir, hams_exe, 0.01, hams_opts);
                n_done = n_done + 1;
            end

            % Persist the updated cache.
            if n_done > 0
                hydro_table = hydro_cache;
                save(cache_file, 'hydro_table', '-v7.3');
                fprintf('\n  Hydro cache updated: %d entries total (%d refinement run(s)).\n', ...
                        length(hydro_cache.drafts), n_done);
            else
                fprintf('  No refinement needed — all interval gradients < 1 s/m.\n');
            end

        end  % adaptive refinement block
end
