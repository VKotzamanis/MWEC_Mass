classdef MWEC_WaveClimate_Plots
%MWEC_WAVECLIMATE_PLOTS  Station-only figures for a WIS climate grid.
%
%   Everything here is derived from <station>_climate_grid.mat alone — no
%   body, no BEM, no PTO, no tuning results.  Use these to characterise the
%   resource before any device is placed in it.  The WEC-coupled figures
%   live in MWEC_Tuning_Plots (fig_t1..fig_t8) and are not duplicated here.
%
%   Style follows MWEC_Tuning_Plots (Wong colorblind-safe palette, fixed cm
%   sizes) with Times New Roman as the body font; cividis carries every
%   colour-mapped field.  The T figures keep their own typeface.
%
%   Inventory
%     fig_w1_scatter    (Hs,Te) occurrence scatter table + marginals
%     fig_w2_spectrum   S_ew(omega) decomposed by sea state + J(omega)
%     fig_w3_energy     per-cell share of the resource + energy concentration
%     fig_w4_stations   cross-station spectra and resource summary
%     fig_w5_histograms T_e and H_s occurrence histograms, station per column
%     fig_w6_energy_surface  3-D surface of probability-weighted wave energy
%
%   Usage
%     cg = load('WIS_Output_WAM/ST63044_climate_grid.mat').climateGrid;
%     MWEC_WaveClimate_Plots.fig_w1_scatter(cg, struct('fig_dir', 'figs'));
%
%   Schema: STABILITY_HANDOFF_PLAN v2.1.0 (empirical-only climate grids).

    methods (Static)

        %% =================================================================
        %%  STYLE  (serif house style for the W family)
        %% =================================================================

        function s = style()
        %STYLE  Wong colorblind-safe palette on a Times New Roman body font.
        %   Same geometry as MWEC_Tuning_Plots.style so the W and T figures
        %   sit together in a document; only the typeface differs.  Keep the
        %   'tex' interpreter everywhere — 'latex' would silently swap in
        %   Computer Modern and break the typographic match.
            s = MWEC_Tuning_Plots.style();
            s.font   = MWEC_WaveClimate_Plots.serif_font();
            s.fs.tick    = 9;   s.fs.label = 11;  s.fs.title = 11;
            s.fs.sgtitle = 13;  s.fs.annot = 9;   s.fs.legend = 9;
            s.fs.cell    = 6.5;                 % in-cell scatter-table numerals
            s.c.mbar = [0, 0.4470, 0.7410];     % MATLAB default blue, as in the slide
        end

        function f = serif_font()
        %SERIF_FONT  Times New Roman where installed, nearest metric clone otherwise.
        %   listfonts costs ~1 s, so resolve once per session.
            persistent chosen
            if ~isempty(chosen), f = chosen; return; end
            preferred = {'Times New Roman', 'Times', 'Liberation Serif', ...
                         'Nimbus Roman', 'FreeSerif', 'DejaVu Serif'};
            available = listfonts;
            chosen = 'Times New Roman';
            for i = 1:numel(preferred)
                if any(strcmpi(available, preferred{i}))
                    chosen = preferred{i};  break
                end
            end
            if ~strcmp(chosen, 'Times New Roman')
                warning('MWEC_WaveClimate_Plots:font', ...
                        'Times New Roman not installed; using %s.', chosen);
            end
            f = chosen;
        end

        %% =================================================================
        %%  DERIVED QUANTITIES  (station data only)
        %% =================================================================

        function k = derive(cg, opts)
        %DERIVE  Spectral aggregates for one climate grid.
        %
        %   Mirrors MWEC_Tuning_Kernels.build_S_ew / compute_iec_moments but
        %   stays self-contained so these figures render without the tuning
        %   pipeline on the path.  Returns:
        %     Sew      probability-weighted climate spectrum  [m^2 s/rad]
        %     J        energy-flux density rho g c_g(w) Sew   [W/m per rad/s]
        %     P_wave   total resource trapz(J)                [W/m]
        %     E        per-cell contribution to P_wave        [W/m], nHs x nTe
        %     band     central 90% energy band (omega and period)
        %     IEC      spectral moments per IEC TS 62600-101
        %
        %   The omnidirectional energy flux of a sea state is
        %       J = rho g * integral c_g(w) S(w) dw            [W/m]
        %   with the group velocity from the full dispersion relation,
        %       c_g = 1/2 (1 + 2kh/sinh 2kh) w/k,   w^2 = g k tanh(kh).
        %
        %   opts.depth_model 'finite' (default) uses that c_g at the station's
        %   own water depth; 'deep' collapses it to c_g = g/(2w), the form
        %   MWEC_Tuning_Kernels.build_S_ew uses.  Deep water needs kh > pi, and
        %   at 28 m and 50 m these stations do not reach it at their energy
        %   period, where the deep form understates the flux by 6-8%.  The
        %   spectra and every IEC moment are independent of this choice; only
        %   J, E and P_wave move.
            if nargin < 2, opts = struct(); end
            if ~isfield(opts, 'depth_model') || isempty(opts.depth_model)
                opts.depth_model = 'finite';
            end
            rho = 1025; g_acc = 9.81;

            required = {'probability_grid','S_omega_grid','omega', ...
                        'Hs_centers','Te_centers'};
            for i = 1:numel(required)
                if ~isfield(cg, required{i})
                    error('MWEC_WaveClimate_Plots:derive:missingField', ...
                          'climateGrid missing required field ''%s''.', required{i});
                end
            end

            p     = cg.probability_grid;                 % nHs x nTe
            omega = cg.omega(:);                         % nOm x 1
            [N_H, N_T] = size(p);
            nOm = numel(omega);

            % nHs x nTe x nOm  ->  nOm x nCells, so every aggregate is one matrix op
            S = reshape(cg.S_omega_grid, N_H*N_T, nOm).';
            pv = p(:).';                                 % 1 x nCells

            % Empirical-only contract: an occupied cell must carry a spectrum.
            occupied = pv > 0;
            bad = occupied & ~all(isfinite(S), 1);
            if any(bad)
                error('MWEC_WaveClimate_Plots:derive:noSpectrum', ...
                      ['%d cell(s) have p>0 but a non-finite spectrum. ' ...
                       'Schema v2.1.0 is empirical-only — rebuild the grid.'], nnz(bad));
            end
            S(~isfinite(S)) = 0;                         % empty cells: p=0, no contribution

            [~, ~, depth] = MWEC_WaveClimate_Plots.ident(cg);
            [c_g, depth_used] = MWEC_WaveClimate_Plots.group_velocity( ...
                                    omega, depth, opts.depth_model, g_acc);

            k.omega = omega;
            k.Sew   = S * pv(:);                         % nOm x 1
            k.J     = rho * g_acc .* c_g .* k.Sew;       % W/m per rad/s

            % Per-cell energy flux, then its probability-weighted share
            J_cell = rho * g_acc .* trapz(omega, S .* c_g, 1);
            k.E    = reshape(pv .* J_cell, N_H, N_T);    % W/m, sums to P_wave
            k.P_wave = trapz(omega, k.J);
            k.c_g    = c_g;
            k.depth_model = depth_used;
            k.water_depth_m = depth;

            k.IEC = MWEC_WaveClimate_Plots.moments(omega, k.Sew);
            k.band = MWEC_WaveClimate_Plots.energy_band(omega, k.Sew);
            k.cdf  = MWEC_WaveClimate_Plots.energy_cdf(omega, k.Sew);

            % m0 budget check against the binned Hs (the same gate the tuning
            % pipeline applies) — a silent mismatch here means a broken grid.
            m0_target = sum(sum(p .* repmat(cg.Hs_centers(:).^2, 1, N_T))) / 16;
            k.m0_target  = m0_target;
            k.m0_err_pct = 100 * abs(k.IEC.m0 - m0_target) / max(m0_target, eps);

            k.n_cells_occupied = nnz(occupied);
            k.n_cells_total    = N_H * N_T;
        end

        function [c_g, model] = group_velocity(omega, h, model, g_acc)
        %GROUP_VELOCITY  Wave group velocity at angular frequency OMEGA.
        %   'deep'   c_g = g/(2w)                       -- valid for kh > pi
        %   'finite' c_g = 1/2 (1 + 2kh/sinh 2kh) w/k   -- valid at any depth
        %   Falls back to 'deep' when the grid carries no water depth.
            omega = omega(:);
            if strcmpi(model, 'deep') || ~isfinite(h) || h <= 0
                if ~strcmpi(model, 'deep')
                    warning('MWEC_WaveClimate_Plots:noDepth', ...
                            'No usable water depth in the grid; using deep water.');
                end
                model = 'deep';
                c_g = g_acc ./ (2 * max(omega, 1e-10));
                return
            end
            model = 'finite';
            kh = MWEC_WaveClimate_Plots.wave_number(omega, h, g_acc) * h;
            % sinh overflows past ~350; there c_g is deep-water to machine eps
            ratio = zeros(size(kh));
            small = 2*kh < 350;
            ratio(small) = 2*kh(small) ./ sinh(2*kh(small));
            n = 0.5 * (1 + ratio);
            c_g = n .* omega ./ (kh / h);
        end

        function k = wave_number(omega, h, g_acc)
        %WAVE_NUMBER  Solve w^2 = g k tanh(k h) for k.
        %   Guo's (2002) explicit initial guess, then Newton.  Three sweeps
        %   reach machine precision over 0.3-6 rad/s at 28-427 m; five is free.
            omega = omega(:);
            x  = omega.^2 * h / g_acc;
            kh = x ./ (1 - exp(-x.^1.25)).^0.4;
            k  = kh / h;
            for it = 1:5
                t  = tanh(k * h);
                f  = g_acc * k .* t - omega.^2;
                df = g_acc * t + g_acc * k * h .* (1 - t.^2);
                k  = max(k - f ./ df, 1e-12);
            end
        end

        function IEC = moments(omega, S)
        %MOMENTS  IEC TS 62600-101 omega-form spectral moments.
            omega = omega(:);  S = S(:);
            m0  = trapz(omega, S);
            m_1 = trapz(omega, S ./ max(omega, 1e-12));
            m1  = trapz(omega, S .* omega);
            m2  = trapz(omega, S .* omega.^2);
            [~, ip] = max(S);
            IEC = struct('m0', m0, 'm_1', m_1, 'm1', m1, 'm2', m2, ...
                'Hm0', 4*sqrt(max(m0, 0)), ...
                'Te',  2*pi * m_1 / max(m0, eps), ...
                'T01', 2*pi * m0  / max(m1, eps), ...
                'T02', 2*pi * sqrt(max(m0, 0)/max(m2, eps)), ...
                'Tp',  2*pi / max(omega(ip), eps), ...
                'omega_p', omega(ip), ...
                'eps_bw', sqrt(max(0, 1 - m1.^2/(max(m0, eps)*max(m2, eps)))));
        end

        function c = energy_cdf(omega, S)
        %ENERGY_CDF  Normalised cumulative energy along omega.
            c = cumtrapz(omega(:), S(:));
            c = c / max(c(end), eps);
        end

        function band = energy_band(omega, S, pct_low, pct_high)
        %ENERGY_BAND  Central energy band containing (pct_high - pct_low) of m0.
            if nargin < 3 || isempty(pct_low),  pct_low  = 0.05; end
            if nargin < 4 || isempty(pct_high), pct_high = 0.95; end
            omega = omega(:);
            c = MWEC_WaveClimate_Plots.energy_cdf(omega, S);
            % interp1 needs strictly increasing samples; flat tails carry no energy
            keep = [true; diff(c) > 0];
            omega_L = interp1(c(keep), omega(keep), pct_low,  'linear', omega(1));
            omega_H = interp1(c(keep), omega(keep), pct_high, 'linear', omega(end));
            band = struct('omega_L', omega_L, 'omega_H', omega_H, ...
                          'T_L', 2*pi/omega_H, 'T_H', 2*pi/omega_L, ...
                          'pct_low', pct_low, 'pct_high', pct_high);
        end

        %% =================================================================
        %%  W1 — Occurrence scatter table
        %% =================================================================

        function fig_w1_scatter(cg, cfg)
        %FIG_W1_SCATTER  (Hs,Te) occurrence table with Hs and Te marginals.
        %   Log colour scale — a wave climate spans 3-4 decades of occurrence,
        %   and a linear scale renders every rare-but-severe cell as one flat
        %   colour.  Unoccupied cells stay grid-grey so "never observed" never
        %   reads as "0.0%".
            if nargin < 2, cfg = struct(); end
            s = MWEC_WaveClimate_Plots.style();
            k = MWEC_WaveClimate_Plots.derive(cg, cfg);

            Ppct = cg.probability_grid * 100;
            Hs_c = cg.Hs_centers(:);  Te_c = cg.Te_centers(:);

            fig = MWEC_Tuning_Plots.new_figure(s, s.sz.full);
            tl = tiledlayout(fig, 4, 5, 'TileSpacing', 'compact', 'Padding', 'compact');

            % --- Te marginal (top) ---
            ax_t = nexttile(tl, 1, [1 4]);
            MWEC_Tuning_Plots.style_axes(ax_t, s);
            bar(ax_t, Te_c, sum(Ppct, 1), 0.92, 'FaceColor', s.c.blue, 'EdgeColor', 'none');
            ylabel(ax_t, sprintf('%% of\nrecords'), 'FontSize', s.fs.annot);
            set(ax_t, 'XTickLabel', []);
            title(ax_t, 'marginal distributions', 'FontSize', s.fs.annot, 'FontWeight', 'normal');

            ax_c = nexttile(tl, 5);  axis(ax_c, 'off');   % corner spacer

            % --- scatter table ---
            ax = nexttile(tl, 6, [3 4]);
            MWEC_Tuning_Plots.style_axes(ax, s);  hold(ax, 'on');
            set(ax, 'Color', s.c.lgrey);          % shows through where no records
            im = imagesc(ax, Te_c, Hs_c, Ppct);
            set(im, 'AlphaData', double(Ppct > 0));
            set(ax, 'YDir', 'normal', 'ColorScale', 'log');
            colormap(ax, MWEC_WaveClimate_Plots.cividis());
            c_lo = max(min(Ppct(Ppct > 0)), 1e-2);
            if isempty(c_lo), c_lo = 1e-2; end
            caxis(ax, [c_lo, max(max(Ppct(:)), 10*c_lo)]);

            MWEC_WaveClimate_Plots.label_cells(ax, Te_c, Hs_c, Ppct, 0.05, '%.1f', s);

            if isfield(cg, 'Hs_P98')
                yline(ax, cg.Hs_P98, '--', 'Color', s.c.red, 'LineWidth', s.lw.secondary, ...
                      'Label', sprintf('H_{s,P98} = %.2f m', cg.Hs_P98), ...
                      'LabelHorizontalAlignment', 'left', 'FontSize', s.fs.annot);
            end
            plot(ax, k.IEC.Te, k.IEC.Hm0, 'o', 'MarkerFaceColor', s.c.orange, ...
                 'MarkerEdgeColor', s.c.black, 'MarkerSize', s.mk.size_lg, 'LineWidth', s.lw.ref);
            text(ax, k.IEC.Te, k.IEC.Hm0, '  mean sea state', 'FontName', s.font, ...
                 'FontSize', s.fs.annot, 'VerticalAlignment', 'middle');

            xlabel(ax, 'Energy period  T_e  (s)', 'FontSize', s.fs.label);
            ylabel(ax, 'Significant wave height  H_s  (m)', 'FontSize', s.fs.label);

            annot = {sprintf('H_{m0} = %.2f m', k.IEC.Hm0); ...
                     sprintf('T_e = %.2f s', k.IEC.Te); ...
                     sprintf('T_p = %.2f s', k.IEC.Tp); ...
                     sprintf('P_{wave} = %.2f kW/m', k.P_wave/1000); ...
                     sprintf('cells occupied: %d/%d', k.n_cells_occupied, k.n_cells_total)};
            text(ax, 0.985, 0.97, annot, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'top', 'FontName', s.font, 'FontSize', s.fs.annot, ...
                 'BackgroundColor', [1 1 1], 'Margin', 4, 'EdgeColor', s.c.grey);

            cb = colorbar(ax);
            cb.Layout.Tile = 'east';          % outside the layout, clear of the marginals
            cb.Label.String = 'occurrence  (% of records)';
            cb.Label.FontSize = s.fs.label;  cb.FontSize = s.fs.tick;

            % --- Hs marginal (right) ---
            ax_r = nexttile(tl, 10, [3 1]);
            MWEC_Tuning_Plots.style_axes(ax_r, s);
            barh(ax_r, Hs_c, sum(Ppct, 2), 0.92, 'FaceColor', s.c.green, 'EdgeColor', 'none');
            set(ax_r, 'YTickLabel', []);
            xlabel(ax_r, '% of records', 'FontSize', s.fs.annot);

            linkaxes([ax_t, ax], 'x');
            linkaxes([ax_r, ax], 'y');

            title(tl, MWEC_WaveClimate_Plots.header(cg, 'W1 — occurrence scatter table'), ...
                  'FontName', s.font, 'FontSize', s.fs.sgtitle, 'FontWeight', 'bold');
            MWEC_WaveClimate_Plots.finish(fig, cg, 'fig_w1_scatter', cfg, 'image');
        end

        %% =================================================================
        %%  W2 — Climate spectrum and energy flux
        %% =================================================================

        function fig_w2_spectrum(cg, cfg)
        %FIG_W2_SPECTRUM  S_ew(omega) decomposed into its sea states, plus J(omega).
        %   The top panel stacks p_ij*S_ij for the five biggest energy
        %   contributors over the remainder, so the areas sum exactly to
        %   S_ew — it shows which sea states build the climate spectrum.
        %   Overlaying raw per-cell spectra instead would bury S_ew, which is
        %   smaller than any single cell by construction (it is a p-weighted mean).
            if nargin < 2, cfg = struct(); end
            s = MWEC_WaveClimate_Plots.style();
            k = MWEC_WaveClimate_Plots.derive(cg, cfg);

            om = k.omega;  p = cg.probability_grid;
            [N_H, N_T] = size(p);
            S = reshape(cg.S_omega_grid, N_H*N_T, numel(om)).';
            S(~isfinite(S)) = 0;
            contrib = S .* p(:).';                        % nOm x nCells, sums to S_ew

            n_top = min(5, k.n_cells_occupied);
            [~, order] = sort(k.E(:), 'descend');
            top = order(1:n_top);
            rest = setdiff(1:numel(p), top);

            stack = [sum(contrib(:, rest), 2), contrib(:, top)];
            top_cols = {s.c.red, s.c.orange, s.c.green, s.c.purple, s.c.sky};
            E_tot = sum(k.E(:));

            fig = MWEC_Tuning_Plots.new_figure(s, s.sz.full);
            tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

            % --- (a) climate spectrum, decomposed ---
            ax1 = nexttile(tl);  MWEC_Tuning_Plots.style_axes(ax1, s);
            hold(ax1, 'on');  grid(ax1, 'on');
            h = area(ax1, om, stack, 'LineStyle', 'none');
            h(1).FaceColor = s.c.lgrey;  h(1).FaceAlpha = 0.9;
            h(1).DisplayName = sprintf('all other cells (%d)', k.n_cells_occupied - n_top);
            for i = 1:n_top
                [ih, it] = ind2sub([N_H, N_T], top(i));
                h(i+1).FaceColor = top_cols{i};  h(i+1).FaceAlpha = 0.9;
                h(i+1).DisplayName = sprintf('H_s=%.2f m, T_e=%.1f s  (%.0f%% of E)', ...
                    cg.Hs_centers(ih), cg.Te_centers(it), 100*k.E(top(i))/E_tot);
            end
            plot(ax1, om, k.Sew, '-', 'Color', s.c.black, 'LineWidth', s.lw.primary, ...
                 'DisplayName', 'S_{ew}(\omega)');
            xline(ax1, k.band.omega_L, ':', 'Color', s.c.black, 'LineWidth', s.lw.ref, 'HandleVisibility', 'off');
            xline(ax1, k.band.omega_H, ':', 'Color', s.c.black, 'LineWidth', s.lw.ref, 'HandleVisibility', 'off');

            xlim(ax1, MWEC_WaveClimate_Plots.omega_view(k));
            xlabel(ax1, '\omega  (rad/s)', 'FontSize', s.fs.label);
            ylabel(ax1, 'S_{ew}  (m^2 s/rad)', 'FontSize', s.fs.label);
            title(ax1, sprintf('Climate spectrum by contributing sea state  —  90%% energy band T = [%.1f, %.1f] s', ...
                  k.band.T_L, k.band.T_H), 'FontSize', s.fs.title);
            lg = legend(ax1, 'Location', 'northeast', 'FontSize', s.fs.legend);
            set(lg, 'Box', 'off');

            % --- (b) energy-flux density + cumulative energy ---
            ax2 = nexttile(tl);  MWEC_Tuning_Plots.style_axes(ax2, s);
            hold(ax2, 'on');  grid(ax2, 'on');
            yyaxis(ax2, 'left');
            in_band = (om >= k.band.omega_L) & (om <= k.band.omega_H);
            fill(ax2, [om(in_band); flipud(om(in_band))], ...
                      [k.J(in_band)/1000; zeros(nnz(in_band), 1)], ...
                 s.c.red, 'EdgeColor', 'none', 'FaceAlpha', 0.18);
            plot(ax2, om, k.J/1000, '-', 'Color', s.c.red, 'LineWidth', s.lw.primary);
            ylabel(ax2, 'J(\omega)  (kW/m per rad/s)', 'FontSize', s.fs.label);
            set(ax2, 'YColor', s.c.red);

            yyaxis(ax2, 'right');
            plot(ax2, om, 100*k.cdf, '--', 'Color', s.c.grey, 'LineWidth', s.lw.secondary);
            ylabel(ax2, 'cumulative energy  (%)', 'FontSize', s.fs.label);
            ylim(ax2, [0, 100]);  set(ax2, 'YColor', s.c.grey);

            xlim(ax2, MWEC_WaveClimate_Plots.omega_view(k));
            xlabel(ax2, '\omega  (rad/s)', 'FontSize', s.fs.label);
            title(ax2, sprintf('Energy-flux density  —  total resource P_{wave} = %.2f kW/m', ...
                  k.P_wave/1000), 'FontSize', s.fs.title);

            title(tl, MWEC_WaveClimate_Plots.header(cg, 'W2 — climate spectrum & energy flux'), ...
                  'FontName', s.font, 'FontSize', s.fs.sgtitle, 'FontWeight', 'bold');
            MWEC_WaveClimate_Plots.finish(fig, cg, 'fig_w2_spectrum', cfg, 'vector');
        end

        %% =================================================================
        %%  W3 — Where the energy is
        %% =================================================================

        function fig_w3_energy(cg, cfg)
        %FIG_W3_ENERGY  Per-cell share of the resource + energy concentration.
        %   Occurrence and energy peak in different places: the most frequent
        %   sea state is rarely the most energetic one, because flux scales as
        %   Hs^2 Te.  W1 answers "how often", this answers "how much".
            if nargin < 2, cfg = struct(); end
            s = MWEC_WaveClimate_Plots.style();
            k = MWEC_WaveClimate_Plots.derive(cg, cfg);

            Hs_c = cg.Hs_centers(:);  Te_c = cg.Te_centers(:);
            Epct = 100 * k.E / max(sum(k.E(:)), eps);

            fig = MWEC_Tuning_Plots.new_figure(s, s.sz.wide);
            tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

            % --- energy matrix ---
            ax = nexttile(tl, 1, [1 2]);
            MWEC_Tuning_Plots.style_axes(ax, s);  hold(ax, 'on');
            set(ax, 'Color', s.c.lgrey);
            im = imagesc(ax, Te_c, Hs_c, Epct);
            set(im, 'AlphaData', double(Epct > 0));
            set(ax, 'YDir', 'normal');
            colormap(ax, MWEC_WaveClimate_Plots.cividis());
            caxis(ax, [0, max(max(Epct(:)), eps)]);

            MWEC_WaveClimate_Plots.label_cells(ax, Te_c, Hs_c, Epct, 0.5, '%.0f', s);

            plot(ax, k.IEC.Te, k.IEC.Hm0, 'o', 'MarkerFaceColor', 'none', ...
                 'MarkerEdgeColor', s.c.black, 'MarkerSize', s.mk.size_lg, 'LineWidth', s.lw.secondary);
            xlabel(ax, 'Energy period  T_e  (s)', 'FontSize', s.fs.label);
            ylabel(ax, 'Significant wave height  H_s  (m)', 'FontSize', s.fs.label);
            title(ax, sprintf('Energy contribution per sea state  (P_{wave} = %.2f kW/m)', ...
                  k.P_wave/1000), 'FontSize', s.fs.title);
            cb = colorbar(ax);
            cb.Layout.Tile = 'east';
            cb.Label.String = 'share of total wave energy  (%)';
            cb.Label.FontSize = s.fs.label;  cb.FontSize = s.fs.tick;

            % --- concentration curve ---
            ax2 = nexttile(tl);
            MWEC_Tuning_Plots.style_axes(ax2, s);  hold(ax2, 'on');  grid(ax2, 'on');
            e = sort(Epct(Epct > 0), 'descend');
            cum = cumsum(e);
            plot(ax2, 1:numel(cum), cum, '-o', 'Color', s.c.purple, ...
                 'MarkerFaceColor', s.c.purple, 'MarkerSize', 3, 'LineWidth', s.lw.secondary);
            for lv = [50, 80]
                n = find(cum >= lv, 1, 'first');
                if isempty(n), continue; end
                yline(ax2, lv, ':', 'Color', s.c.grey, 'LineWidth', s.lw.ref);
                xline(ax2, n,  ':', 'Color', s.c.grey, 'LineWidth', s.lw.ref);
                text(ax2, n, lv - 8, sprintf(' %d cells \\rightarrow %d%%', n, lv), ...
                     'FontName', s.font, 'FontSize', s.fs.annot);
            end
            ylim(ax2, [0, 102]);  xlim(ax2, [0, numel(cum) + 1]);
            xlabel(ax2, 'sea-state cells, ranked', 'FontSize', s.fs.label);
            ylabel(ax2, 'cumulative energy  (%)', 'FontSize', s.fs.label);
            title(ax2, sprintf('Concentration (%d occupied)', numel(cum)), 'FontSize', s.fs.title);

            title(tl, MWEC_WaveClimate_Plots.header(cg, 'W3 — where the wave energy is'), ...
                  'FontName', s.font, 'FontSize', s.fs.sgtitle, 'FontWeight', 'bold');
            MWEC_WaveClimate_Plots.finish(fig, cg, 'fig_w3_energy', cfg, 'image');
        end

        %% =================================================================
        %%  W4 — Cross-station comparison
        %% =================================================================

        function fig_w4_stations(cgs, cfg)
        %FIG_W4_STATIONS  Climate spectra and resource summary across stations.
        %   CGS is a cell array of climateGrid structs.  Hm0, Te and P_wave get
        %   one panel each — they are metres, seconds and kW/m, and sharing an
        %   axis between them would only make the tallest unit legible.
            if nargin < 2, cfg = struct(); end
            s = MWEC_WaveClimate_Plots.style();
            if ~iscell(cgs), cgs = {cgs}; end
            n = numel(cgs);
            if n == 0, return; end

            ks = cell(1, n);  names = cell(1, n);  legends = cell(1, n);
            for i = 1:n
                ks{i} = MWEC_WaveClimate_Plots.derive(cgs{i}, cfg);
                [sid, reg, depth] = MWEC_WaveClimate_Plots.ident(cgs{i});
                names{i}   = reg;
                legends{i} = sprintf('%s · %s  (h = %.0f m)', sid, reg, depth);
            end
            cols = {s.c.blue, s.c.orange, s.c.green, s.c.purple, s.c.red, s.c.sky};

            fig = MWEC_Tuning_Plots.new_figure(s, s.sz.full);
            tl = tiledlayout(fig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

            % --- spectra overlay ---
            ax = nexttile(tl, 1, [3 2]);
            MWEC_Tuning_Plots.style_axes(ax, s);  hold(ax, 'on');  grid(ax, 'on');
            lo = inf;  hi = -inf;
            for i = 1:n
                c = cols{mod(i-1, numel(cols)) + 1};
                plot(ax, ks{i}.omega, ks{i}.Sew, '-', 'Color', c, ...
                     'LineWidth', s.lw.primary, 'DisplayName', legends{i});
                v = MWEC_WaveClimate_Plots.omega_view(ks{i});
                lo = min(lo, v(1));  hi = max(hi, v(2));
            end
            xlim(ax, [lo, hi]);
            xlabel(ax, '\omega  (rad/s)', 'FontSize', s.fs.label);
            ylabel(ax, 'S_{ew}  (m^2 s/rad)', 'FontSize', s.fs.label);
            title(ax, 'Probability-weighted climate spectra', 'FontSize', s.fs.title);
            lg = legend(ax, 'Location', 'northeast', 'FontSize', s.fs.legend);
            set(lg, 'Box', 'off');

            % --- one metric per panel ---
            metrics = { 'Hm0',    'H_{m0}  (m)',        s.c.blue,  @(x) x.IEC.Hm0; ...
                        'Te',     'T_e  (s)',           s.c.green, @(x) x.IEC.Te; ...
                        'P_wave', 'P_{wave}  (kW/m)',   s.c.red,   @(x) x.P_wave/1000 };
            for m = 1:size(metrics, 1)
                axm = nexttile(tl, 3*m);
                MWEC_Tuning_Plots.style_axes(axm, s);  hold(axm, 'on');  grid(axm, 'on');
                v = cellfun(metrics{m, 4}, ks);
                bar(axm, 1:n, v, 0.62, 'FaceColor', metrics{m, 3}, 'EdgeColor', 'none');
                for i = 1:n
                    text(axm, i, v(i), sprintf('%.2f', v(i)), 'HorizontalAlignment', 'center', ...
                         'VerticalAlignment', 'bottom', 'FontName', s.font, 'FontSize', s.fs.annot);
                end
                ylim(axm, [0, max(v) * 1.25]);
                ylabel(axm, metrics{m, 2}, 'FontSize', s.fs.label);
                set(axm, 'XTick', 1:n, 'XLim', [0.4, n + 0.6]);
                if m == size(metrics, 1)
                    set(axm, 'XTickLabel', names, 'XTickLabelRotation', 20);
                else
                    set(axm, 'XTickLabel', []);
                end
            end

            title(tl, 'W4 — cross-station wave-climate comparison', ...
                  'FontName', s.font, 'FontSize', s.fs.sgtitle, 'FontWeight', 'bold');
            MWEC_Tuning_Plots.save_fig(fig, 'fig_w4_stations', cfg);
        end

        %% =================================================================
        %%  W5 — Marginal histograms, one column per station
        %% =================================================================

        function fig_w5_histograms(cgs, cfg)
        %FIG_W5_HISTOGRAMS  T_e and H_s occurrence histograms, station per column.
        %
        %   Axes are positioned by hand rather than by tiledlayout so the
        %   columns line up exactly and the rotated H_s interval labels get the
        %   room they need.
        %
        %   Resolution note: the bars are the climate grid's own bins —
        %   dT_e = 1.0 s and dH_s = 0.5 m as built.  The per-record T_e values
        %   are not retained in the grid, so a finer T_e histogram means
        %   lowering dbinTe in run_buildClimateGrid.m and rebuilding from the
        %   raw WIS files.
        %
        %   cfg.station_labels  optional struct mapping station id to the header
        %                       text, e.g. struct('ST63044', 'North Atlantic, NH').
            if nargin < 2, cfg = struct(); end
            if ~iscell(cgs), cgs = {cgs}; end
            n = numel(cgs);
            if n == 0, return; end
            s = MWEC_WaveClimate_Plots.style();

            width_cm = min(6.4 * n + 1.4, 34);
            fig = MWEC_Tuning_Plots.new_figure(s, [width_cm, 13.0]);

            L = 0.070;  R = 0.015;  gap = 0.052;
            colw = (1 - L - R - (n-1)*gap) / n;
            y_te = 0.575;  y_hs = 0.175;  rowh = 0.235;

            for i = 1:n
                cg = cgs{i};
                x0 = L + (i-1)*(colw + gap);
                [sid, region, depth] = MWEC_WaveClimate_Plots.ident(cg);
                label = MWEC_WaveClimate_Plots.station_label(sid, region, cfg);

                p    = cg.probability_grid;
                Te_c = cg.Te_centers(:);   Hs_c = cg.Hs_centers(:);
                Te_e = cg.Te_edges(:);     Hs_e = cg.Hs_edges(:);
                pct_Te = 100 * sum(p, 1).';       % nTe x 1
                pct_Hs = 100 * sum(p, 2);         % nHs x 1

                % --- column header ---
                annotation(fig, 'textbox', [x0 - 0.012, 0.930, colw + 0.05, 0.055], ...
                    'String', sprintf('%s:', label), 'FontName', s.font, ...
                    'FontSize', s.fs.sgtitle, 'FontWeight', 'bold', ...
                    'EdgeColor', 'none', 'VerticalAlignment', 'middle', ...
                    'Interpreter', 'tex');
                if isfinite(depth)
                    annotation(fig, 'textbox', [x0 + 0.008, 0.882, colw + 0.05, 0.050], ...
                        'String', sprintf('\\bullet  Water Depth: %.0f meters', depth), ...
                        'FontName', s.font, 'FontSize', s.fs.label, ...
                        'EdgeColor', 'none', 'VerticalAlignment', 'middle', ...
                        'Interpreter', 'tex');
                end

                % --- T_e histogram ---
                ax = axes(fig, 'Position', [x0, y_te, colw, rowh]); %#ok<LAXES>
                MWEC_WaveClimate_Plots.style_hist_axes(ax, s);
                bar(ax, Te_c, pct_Te, 1.0, 'FaceColor', s.c.mbar, 'EdgeColor', 'none');
                hold(ax, 'on');
                % continuous reconstruction of the same marginal, for shape
                Te_f = linspace(0, Te_e(end), 400).';
                dT = median(diff(Te_e));
                K = MWEC_WaveClimate_Plots.reflected_gaussian(Te_f, Te_c, 0.50*dT);
                plot(ax, Te_f, 100 * dT * (K * (pct_Te/100)), '-', ...
                     'Color', s.c.grey, 'LineWidth', s.lw.secondary);
                xlim(ax, [0, Te_e(end)]);
                ylim(ax, [0, 1.15 * max(pct_Te)]);
                xlabel(ax, 'T_e  (s)', 'FontSize', s.fs.label);
                if i == 1, ylabel(ax, 'Occurrence  (%)', 'FontSize', s.fs.label); end

                % --- H_s histogram ---
                ax = axes(fig, 'Position', [x0, y_hs, colw, rowh]); %#ok<LAXES>
                MWEC_WaveClimate_Plots.style_hist_axes(ax, s);
                bar(ax, 1:numel(Hs_c), pct_Hs, 0.86, 'FaceColor', s.c.mbar, 'EdgeColor', 'none');
                ticks = cell(numel(Hs_c), 1);
                for b = 1:numel(Hs_c)
                    ticks{b} = sprintf('[%.1f, %.1f]', Hs_e(b), Hs_e(b+1));
                end
                set(ax, 'XTick', 1:numel(Hs_c), 'XTickLabel', ticks, ...
                        'XTickLabelRotation', 30, 'XLim', [0.4, numel(Hs_c) + 0.6]);
                ylim(ax, [0, 1.15 * max(pct_Hs)]);
                xlabel(ax, 'H_s  (m)', 'FontSize', s.fs.label);
                if i == 1, ylabel(ax, 'Occurrence  (%)', 'FontSize', s.fs.label); end
            end

            MWEC_Tuning_Plots.save_fig(fig, 'fig_w5_histograms', cfg);
        end

        function style_hist_axes(ax, s)
        %STYLE_HIST_AXES  Boxed, inward-ticked axes — the slide's convention.
            set(ax, 'FontName', s.font, 'FontSize', s.fs.tick, 'LineWidth', s.lw.axes, ...
                    'TickDir', 'in', 'TickLength', [0.015, 0.015], 'Box', 'on', ...
                    'Layer', 'top', 'XGrid', 'off', 'YGrid', 'off');
        end

        %% =================================================================
        %%  W6 — 3-D energy-density surface
        %% =================================================================

        function fig_w6_energy_surface(cg, cfg, opts)
        %FIG_W6_ENERGY_SURFACE  Probability-weighted wave energy over (T_e, H_s).
        %
        %   Z is the probability-weighted energy flux per unit of the (H_s,T_e)
        %   plane, so the volume under the surface is the station's total
        %   resource in kW/m.  Each cell's flux comes from the spectral form
        %       J_ij = rho g^2 / 2 * integral S_ij(w)/w dw
        %   over that cell's stored empirical spectrum, never from the
        %   parametric rho g^2 Hs^2 Te / (64 pi).
        %
        %   The stored grid is only 4-8 H_s bins by 12-16 T_e bins, so surfing
        %   it raw gives a staircase.  Smoothing convolves the per-cell weights
        %   p_ij*J_ij with a Gaussian of one bin width, reflected about H_s = 0
        %   and T_e = 0 so no energy leaks to negative values.  Convolution
        %   preserves the integral exactly: the printed volume is the check.
        %
        %   opts.smooth        true (default) for the density surface, false to
        %                      surf the native cells in kW/m
        %   opts.sigma_factor  kernel width in bin widths (default 0.75)
        %   opts.n_Hs/opts.n_Te  fine-grid resolution (default 181 x 241)
            if nargin < 2, cfg  = struct(); end
            if nargin < 3, opts = struct(); end
            if ~isfield(opts, 'smooth'), opts.smooth = true; end
            s = MWEC_WaveClimate_Plots.style();

            fig = MWEC_Tuning_Plots.new_figure(s, [16.0, 12.5]);
            ax = axes(fig);
            set(ax, 'FontName', s.font, 'FontSize', s.fs.tick, 'LineWidth', s.lw.axes, ...
                    'Box', 'on', 'BoxStyle', 'full', 'TickDir', 'out', ...
                    'GridAlpha', 0.15, 'GridColor', [0.4 0.4 0.4]);
            hold(ax, 'on');  grid(ax, 'on');

            if opts.smooth
                q = MWEC_WaveClimate_Plots.energy_density(cg, opts);
                [T, H] = meshgrid(q.Te, q.Hs);
                Z = q.Z;
                z_label = 'p \cdot J   (kW/m per m\cdots)';
                sub = sprintf(['volume = %.3f kW/m  (\\Sigma p_{ij}J_{ij} = %.3f kW/m)   ' ...
                               '|   \\sigma = %.2f m \\times %.2f s   |   c_g: %s depth'], ...
                               q.volume, q.P_wave, q.sigma_Hs, q.sigma_Te, q.depth_model);
            else
                k = MWEC_WaveClimate_Plots.derive(cg, opts);
                [T, H] = meshgrid(cg.Te_centers(:), cg.Hs_centers(:));
                Z = k.E / 1000;
                q = struct('P_wave', sum(Z(:)), 'volume', sum(Z(:)));
                z_label = 'p \cdot J   (kW/m per cell)';
                sub = sprintf(['native %d \\times %d grid  |  \\Sigma p_{ij}J_{ij} = %.3f kW/m' ...
                               '  |  c_g: %s depth'], ...
                              size(Z, 1), size(Z, 2), q.P_wave, k.depth_model);
            end

            srf = surf(ax, T, H, Z);
            set(srf, 'EdgeColor', 'none', 'FaceColor', 'interp', ...
                     'FaceLighting', 'none', 'AmbientStrength', 1);
            colormap(ax, MWEC_WaveClimate_Plots.cividis());
            caxis(ax, [0, max(max(Z(:)), eps)]);
            zlim(ax, [0, 1.05 * max(max(Z(:)), eps)]);

            % The kernel is evaluated past the last bin edge so no mass is lost
            % at the boundary; the axes stop at the edges, where the data does.
            % A floor contour projection was tried and dropped: the surface's
            % own near-zero skirt covers the whole plane and occludes it.
            xlim(ax, [0, cg.Te_edges(end)]);  ylim(ax, [0, cg.Hs_edges(end)]);
            xlabel(ax, 'Energy period  T_e  (s)',            'FontSize', s.fs.label);
            ylabel(ax, 'Significant wave height  H_s  (m)',  'FontSize', s.fs.label);
            zlabel(ax, z_label,                              'FontSize', s.fs.label);
            view(ax, -40, 32);

            cb = colorbar(ax);
            cb.Label.String = z_label;
            cb.Label.FontSize = s.fs.label;  cb.FontSize = s.fs.tick;
            cb.Label.FontName = s.font;

            ttl = title(ax, MWEC_WaveClimate_Plots.header(cg, 'W6 — available wave energy', true));
            set(ttl, 'FontName', s.font, 'FontSize', s.fs.title, 'FontWeight', 'bold');
            sbt = subtitle(ax, sub);
            set(sbt, 'FontName', s.font, 'FontSize', s.fs.annot, 'Color', s.c.grey);

            MWEC_WaveClimate_Plots.finish(fig, cg, 'fig_w6_surface', cfg, 'image');
        end

        function q = energy_density(cg, opts)
        %ENERGY_DENSITY  Smooth (H_s,T_e) density of probability-weighted flux.
        %   Returns Z in kW/m per (m.s) on a fine grid, with trapz(trapz(Z))
        %   equal to sum(p_ij*J_ij) to within the fine grid's quadrature error.
            if nargin < 2, opts = struct(); end
            if ~isfield(opts, 'sigma_factor'), opts.sigma_factor = 0.75; end
            if ~isfield(opts, 'n_Hs'), opts.n_Hs = 181; end
            if ~isfield(opts, 'n_Te'), opts.n_Te = 241; end

            k = MWEC_WaveClimate_Plots.derive(cg, opts);
            W = k.E / 1000;                                  % kW/m per cell
            Hs_c = cg.Hs_centers(:);  Te_c = cg.Te_centers(:);
            dH = median(diff(cg.Hs_edges(:)));
            dT = median(diff(cg.Te_edges(:)));
            sigma_Hs = opts.sigma_factor * dH;
            sigma_Te = opts.sigma_factor * dT;

            % pad past the last edge so the kernel tails stay inside the domain
            Hs_f = linspace(0, cg.Hs_edges(end) + 1.5*dH, opts.n_Hs).';
            Te_f = linspace(0, cg.Te_edges(end) + 1.5*dT, opts.n_Te).';
            K_Hs = MWEC_WaveClimate_Plots.reflected_gaussian(Hs_f, Hs_c, sigma_Hs);
            K_Te = MWEC_WaveClimate_Plots.reflected_gaussian(Te_f, Te_c, sigma_Te);

            Z = K_Hs * W * K_Te.';                           % nHf x nTf
            q = struct('Te', Te_f, 'Hs', Hs_f, 'Z', Z, ...
                       'depth_model', k.depth_model, ...
                       'P_wave', sum(W(:)), ...
                       'volume', trapz(Hs_f, trapz(Te_f, Z, 2)), ...
                       'sigma_Hs', sigma_Hs, 'sigma_Te', sigma_Te, ...
                       'dH', dH, 'dT', dT, 'E_native_kWm', W);
        end

        function K = reflected_gaussian(x, mu, sigma)
        %REFLECTED_GAUSSIAN  Gaussian kernel matrix folded about the origin.
        %   K(a,b) is the density at X(a) of a unit mass at MU(b).  Folding the
        %   negative tail back keeps the total mass on [0, inf), so smoothing a
        %   cell near H_s = 0 does not lose energy off the edge of the domain.
            x = x(:);  mu = mu(:).';
            phi = @(u) exp(-0.5 * u.^2) ./ (sigma * sqrt(2*pi));
            K = phi((x - mu) / sigma) + phi((x + mu) / sigma);
        end

        function label = station_label(sid, region, cfg)
        %STATION_LABEL  Header text for a station column.
        %   cfg.station_labels overrides; otherwise the stored region is split
        %   at camel-case boundaries ('NorthAtlantic' -> 'North Atlantic').
            if nargin >= 3 && isfield(cfg, 'station_labels') && ...
                    isfield(cfg.station_labels, sid)
                label = cfg.station_labels.(sid);
                return
            end
            label = regexprep(region, '(?<=[a-z])(?=[A-Z])', ' ');
            if strcmp(label, 'unknown'), label = sid; end
        end

        %% =================================================================
        %%  HELPERS
        %% =================================================================

        function label_cells(ax, x, y, V, thresh, fmt, s)
        %LABEL_CELLS  Annotate matrix cells above THRESH, flipping the text to
        %   white over the dark end of the colormap so labels stay legible.
        %   Both figures now use cividis, which runs dark navy at the low end to
        %   yellow at the high end, so one rule covers the linear and the log
        %   axis: white below the crossover, black above it.
            lim = caxis(ax);
            is_log = strcmp(get(ax, 'ColorScale'), 'log');
            for i = 1:numel(y)
                for j = 1:numel(x)
                    if ~(V(i, j) >= thresh), continue; end
                    if is_log
                        f = (log10(max(V(i,j), lim(1))) - log10(lim(1))) / ...
                            max(log10(lim(2)) - log10(lim(1)), eps);
                    else
                        f = (V(i, j) - lim(1)) / max(lim(2) - lim(1), eps);
                    end
                    if f < 0.62, col = [1 1 1]; else, col = s.c.black; end
                    text(ax, x(j), y(i), sprintf(fmt, V(i, j)), ...
                         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                         'FontName', s.font, 'FontSize', s.fs.cell, 'Color', col);
                end
            end
        end

        function cm = cividis(n)
        %CIVIDIS  The cividis colormap (Nunez, Anderton & Renslow, 2018).
        %   Perceptually uniform like viridis, but built so a deuteranope or
        %   protanope reads the same ordering a trichromat does -- viridis is
        %   only approximately safe, cividis is optimised for it.  MATLAB does
        %   not ship it, so the reference 256-entry table is embedded verbatim:
        %   sparse anchors plus interpolation cost ~2% in the blue ramp, which
        %   is visible, and break colour-for-colour agreement with matplotlib.
            persistent T
            if isempty(T)
                T = [ ...
                  0  34  78;   0  35  79;   0  36  81;   0  37  83;   0  37  84;   0  38  86;   0  39  88;   0  40  89;
                  0  40  91;   0  41  93;   0  42  95;   0  42  97;   0  43  98;   0  44 100;   0  44 102;   0  45 104;
                  0  46 106;   0  46 108;   0  47 109;   0  48 111;   0  48 112;   0  49 112;   0  49 113;   1  50 113;
                  5  51 113;   8  51 112;  12  52 112;  15  53 112;  18  53 112;  20  54 112;  22  55 112;  24  55 111;
                 26  56 111;  28  57 111;  30  58 111;  32  58 111;  33  59 110;  35  60 110;  36  60 110;  38  61 110;
                 39  62 110;  41  63 110;  42  63 109;  43  64 109;  45  65 109;  46  65 109;  47  66 109;  49  67 109;
                 50  67 109;  51  68 109;  52  69 108;  53  69 108;  54  70 108;  56  71 108;  57  72 108;  58  72 108;
                 59  73 108;  60  74 108;  61  74 108;  62  75 108;  63  76 108;  64  76 108;  65  77 108;  66  78 108;
                 67  78 108;  68  79 108;  69  80 108;  70  81 108;  71  81 108;  72  82 108;  73  83 108;  74  83 108;
                 75  84 108;  76  85 108;  77  85 108;  78  86 108;  79  87 108;  80  87 108;  81  88 109;  82  89 109;
                 83  90 109;  84  90 109;  85  91 109;  85  92 109;  86  92 109;  87  93 109;  88  94 109;  89  94 110;
                 90  95 110;  91  96 110;  92  97 110;  93  97 110;  94  98 110;  94  99 111;  95  99 111;  96 100 111;
                 97 101 111;  98 101 111;  99 102 112; 100 103 112; 101 104 112; 101 104 112; 102 105 112; 103 106 113;
                104 106 113; 105 107 113; 106 108 113; 107 109 114; 108 109 114; 108 110 114; 109 111 114; 110 111 115;
                111 112 115; 112 113 115; 113 114 116; 114 114 116; 114 115 116; 115 116 117; 116 116 117; 117 117 117;
                118 118 118; 119 119 118; 119 119 119; 120 120 119; 121 121 119; 122 122 120; 123 122 120; 124 123 120;
                125 124 120; 126 124 120; 126 125 120; 127 126 120; 128 127 120; 129 127 120; 130 128 121; 131 129 121;
                132 130 121; 133 130 121; 134 131 121; 135 132 120; 136 133 120; 137 133 120; 138 134 120; 139 135 120;
                140 136 120; 141 136 120; 142 137 120; 143 138 120; 144 139 120; 145 139 120; 146 140 120; 146 141 120;
                147 142 120; 148 142 119; 149 143 119; 150 144 119; 151 145 119; 152 146 119; 153 146 119; 154 147 118;
                155 148 118; 156 149 118; 157 149 118; 158 150 118; 159 151 117; 160 152 117; 161 153 117; 162 153 117;
                163 154 116; 164 155 116; 165 156 116; 166 156 116; 167 157 115; 168 158 115; 169 159 115; 170 160 115;
                171 160 114; 172 161 114; 173 162 114; 174 163 113; 175 164 113; 176 165 113; 177 165 112; 179 166 112;
                180 167 111; 181 168 111; 182 169 111; 183 169 110; 184 170 110; 185 171 109; 186 172 109; 187 173 109;
                188 174 108; 189 174 108; 190 175 107; 191 176 107; 192 177 106; 193 178 106; 194 179 105; 195 179 105;
                196 180 104; 197 181 104; 198 182 103; 199 183 103; 200 184 102; 201 185 101; 203 185 101; 204 186 100;
                205 187  99; 206 188  99; 207 189  98; 208 190  98; 209 191  97; 210 192  96; 211 192  95; 212 193  95;
                213 194  94; 214 195  93; 215 196  92; 217 197  92; 218 198  91; 219 199  90; 220 200  89; 221 200  88;
                222 201  88; 223 202  87; 224 203  86; 225 204  85; 226 205  84; 228 206  83; 229 207  82; 230 208  81;
                231 209  80; 232 210  79; 233 211  78; 234 211  76; 235 212  75; 237 213  74; 238 214  73; 239 215  72;
                240 216  70; 241 217  69; 242 218  68; 243 219  66; 245 220  65; 246 221  63; 247 222  62; 248 223  60;
                249 224  58; 251 225  56; 252 226  54; 253 227  52; 254 228  52; 254 229  53; 254 230  54; 254 232  56 ] / 255;
            end
            if nargin < 1 || isempty(n), n = size(T, 1); end
            if n == size(T, 1)
                cm = T;
            else
                cm = interp1(linspace(0, 1, size(T, 1)), T, linspace(0, 1, n));
            end
        end

        function v = omega_view(k)
        %OMEGA_VIEW  Plot window around the energetic band, not the full grid.
            v = [max(0.2, 0.75 * k.band.omega_L), min(k.omega(end), 1.6 * k.band.omega_H)];
        end

        function [sid, reg, depth, n_bin, n_tot] = ident(cg)
        %IDENT  Station identity, tolerating meta/config/absent placement.
            sid = 'unknown';  reg = 'unknown';  depth = NaN;  n_bin = NaN;  n_tot = NaN;
            src = {};
            if isfield(cg, 'meta'),   src{end+1} = cg.meta;   end
            if isfield(cg, 'config'), src{end+1} = cg.config; end
            for i = 1:numel(src)
                S = src{i};
                if isfield(S, 'station_id')    && strcmp(sid, 'unknown'), sid = char(S.station_id); end
                if isfield(S, 'region')        && strcmp(reg, 'unknown'), reg = char(S.region);     end
                if isfield(S, 'water_depth_m') && ~isfinite(depth), depth = S.water_depth_m;        end
                if isfield(S, 'n_records_binned') && ~isfinite(n_bin), n_bin = S.n_records_binned;  end
                if isfield(S, 'n_records_total')  && ~isfinite(n_tot), n_tot = S.n_records_total;   end
            end
        end

        function str = header(cg, prefix, short)
        %HEADER  Common figure title: who the station is and how much data it has.
        %   SHORT drops the record count, for titles that have to fit a 3-D axes.
            if nargin < 3, short = false; end
            [sid, reg, depth, n_bin, n_tot] = MWEC_WaveClimate_Plots.ident(cg);
            str = sprintf('%s  ·  %s  ·  %s', prefix, sid, reg);
            if isfinite(depth)
                str = sprintf('%s  |  h = %.0f m', str, depth);
            end
            if ~short && isfinite(n_bin) && isfinite(n_tot)
                str = sprintf('%s,  N = %d/%d records', str, round(n_bin), round(n_tot));
            end
        end

        function finish(fig, cg, base, cfg, pdf_content)
        %FINISH  Save with the station id in the file name so runs never collide.
            sid = MWEC_WaveClimate_Plots.ident(cg);
            MWEC_Tuning_Plots.save_fig(fig, sprintf('%s_%s', base, sid), cfg, ...
                                       struct('pdf_content', pdf_content));
        end

    end
end
