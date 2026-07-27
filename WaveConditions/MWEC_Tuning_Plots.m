classdef MWEC_Tuning_Plots
%MWEC_TUNING_PLOTS  Figure-generation static methods for MWEC_Tuning v2.
%
%   Wong (2011) colorblind-safe palette, fixed cm sizes, Helvetica, line
%   weights per WEC_GM plot protocol.  Each figure saves .pdf + .png + .fig.
%   File names carry the closure tag and region so Evans / 40%-crit runs do
%   not overwrite each other.
%
%   Inventory (v2):
%     fig_t1_climate       S_ew(omega) + 90% band + IEC summary
%     fig_t2_bem           A_kk, B_kk per mode + surge-pitch coupling A15/B15
%     fig_t4_placement     coupled absorbed-power density p_k(omega).S_ew + ceiling
%     fig_t5_absorption    per-mode CW(T) vs Falnes ceiling (dipole combined)
%     fig_t6_scatter       (Hs,Te) climate + placed design periods
%     fig_c1_interference  coupled vs decoupled <P_abs> density (interference proof)
%     fig_t7_cross         cross-case / cross-pipeline eta_Falnes, CWR, <P_abs>

    methods (Static)

        %% =================================================================
        %%  STYLE
        %% =================================================================

        function s = style()
            s.c.blue   = [  0, 114, 178] / 255;
            s.c.orange = [230, 159,   0] / 255;
            s.c.green  = [  0, 158, 115] / 255;
            s.c.red    = [213,  94,   0] / 255;
            s.c.purple = [204, 121, 167] / 255;
            s.c.sky    = [ 86, 180, 233] / 255;
            s.c.yellow = [240, 228,  66] / 255;
            s.c.black  = [  0,   0,   0];
            s.c.grey   = [0.45, 0.45, 0.45];
            s.c.lgrey  = [0.82, 0.82, 0.82];

            s.sz.single = [8.8,  7.0];
            s.sz.wide   = [18.0, 7.0];
            s.sz.tall   = [8.8, 10.0];
            s.sz.full   = [18.0, 12.0];

            s.font       = 'Helvetica';
            s.fs.tick    = 8;  s.fs.label = 9;  s.fs.title = 9;
            s.fs.sgtitle = 10; s.fs.annot = 8;  s.fs.legend = 8;

            s.lw.primary = 1.5; s.lw.secondary = 1.0; s.lw.ref = 0.75; s.lw.axes = 0.75;
            s.mk.size = 6; s.mk.size_lg = 9;

            s.mode.surge = s.c.blue;  s.mode.heave = s.c.green;  s.mode.pitch = s.c.orange;
        end

        function style_axes(ax, s)
            set(ax, 'FontName', s.font, 'FontSize', s.fs.tick, 'LineWidth', s.lw.axes, ...
                    'TickDir', 'out', 'TickLength', [0.012, 0.012], ...
                    'GridAlpha', 0.18, 'GridColor', [0.5, 0.5, 0.5], 'Box', 'off');
        end

        function fig = new_figure(s, sz)
            fig = figure('Units', 'centimeters', 'Position', [2, 2, sz(1), sz(2)], ...
                         'Color', 'white', 'PaperUnits', 'centimeters', 'PaperSize', sz, ...
                         'PaperPosition', [0, 0, sz(1), sz(2)]);
        end

        function save_fig(fig, fig_name, cfg, opts)
            if nargin < 4 || isempty(opts), opts = struct(); end
            if ~isfield(opts, 'pdf_content'), opts.pdf_content = 'vector'; end
            if isfield(cfg, 'fig_dir') && ~isempty(cfg.fig_dir), out_dir = cfg.fig_dir;
            elseif isfield(cfg, 'out_dir') && ~isempty(cfg.out_dir), out_dir = cfg.out_dir;
            else, out_dir = pwd; end
            if ~isfolder(out_dir), mkdir(out_dir); end
            base = fullfile(out_dir, fig_name);
            formats = {'pdf', 'png', 'fig'};
            if isfield(cfg, 'plot') && isfield(cfg.plot, 'formats'), formats = cfg.plot.formats; end
            if any(strcmp(formats, 'pdf'))
                exportgraphics(fig, [base '.pdf'], 'ContentType', opts.pdf_content, 'BackgroundColor', 'white');
            end
            if any(strcmp(formats, 'png'))
                exportgraphics(fig, [base '.png'], 'Resolution', 300, 'BackgroundColor', 'white');
            end
            if any(strcmp(formats, 'fig')), savefig(fig, [base '.fig']); end
            if isgraphics(fig, 'figure'), close(fig); end
        end

        %% =================================================================
        %%  T1 — Climate spectrum
        %% =================================================================

        function fig_t1_climate(region, ctag, results, cfg)
            s = MWEC_Tuning_Plots.style();
            fig = MWEC_Tuning_Plots.new_figure(s, s.sz.wide);
            ax = axes(fig); MWEC_Tuning_Plots.style_axes(ax, s); hold(ax, 'on'); grid(ax, 'on');

            om = results.climate.omega;  S_ew = results.climate.S_ew;
            band = results.climate.band; part = results.climate.partition; IEC = results.climate.IEC;

            mask = (om >= band.omega_L) & (om <= band.omega_H);
            fill(ax, [om(mask); flipud(om(mask))], [S_ew(mask); zeros(sum(mask),1)], ...
                 s.c.lgrey, 'EdgeColor', 'none', 'FaceAlpha', 0.40, ...
                 'DisplayName', sprintf('90%% band [%.2f, %.2f] s', band.T_L, band.T_H));
            plot(ax, om, S_ew, '-', 'Color', s.c.blue, 'LineWidth', s.lw.primary, 'DisplayName', 'S_{ew}(\omega)');
            xline(ax, part.partition_omega, '--', 'Color', s.c.grey, 'LineWidth', s.lw.secondary, ...
                  'Label', sprintf('partition T=%.1fs', 2*pi/part.partition_omega), 'HandleVisibility', 'off');

            xlabel(ax, '\omega (rad/s)', 'FontSize', s.fs.label);
            ylabel(ax, 'S_{ew}(\omega)  (m^2 s/rad)', 'FontSize', s.fs.label);
            title(ax, sprintf('T1 — Climate spectrum  %s  [%s]', region, ctag), 'FontSize', s.fs.title);
            legend(ax, 'Location', 'northeast', 'FontSize', s.fs.legend);

            annot = {sprintf('H_{m0}=%.2f m', IEC.Hm0); sprintf('T_e=%.2f s', IEC.Te); ...
                     sprintf('T_p=%.2f s', IEC.Tp); sprintf('\\epsilon=%.3f', IEC.eps_bw)};
            text(ax, 0.98, 0.78, annot, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                 'VerticalAlignment', 'top', 'FontName', s.font, 'FontSize', s.fs.annot, ...
                 'BackgroundColor', [1 1 1 0.85], 'Margin', 4, 'EdgeColor', s.c.grey);

            MWEC_Tuning_Plots.save_fig(fig, sprintf('fig_t1_climate_%s_%s', ctag, region), cfg);
        end

        %% =================================================================
        %%  T2 — BEM verification + surge-pitch coupling
        %% =================================================================

        function fig_t2_bem(region, ctag, results, cfg)
            s = MWEC_Tuning_Plots.style();
            bem = results.bem;  om = bem.omega_BEM;
            fig = MWEC_Tuning_Plots.new_figure(s, s.sz.full);
            labels = {'Surge', 'Heave', 'Pitch'};
            cols   = {s.mode.surge, s.mode.heave, s.mode.pitch};

            for k = 1:3
                ax = subplot(3, 3, k); MWEC_Tuning_Plots.style_axes(ax, s); hold(ax, 'on'); grid(ax, 'on');
                plot(ax, om, bem.A_diag(k, :), '-', 'Color', cols{k}, 'LineWidth', s.lw.primary);
                xlabel(ax, '\omega (rad/s)', 'FontSize', s.fs.label);
                ylabel(ax, sprintf('A_{%d%d}', k, k), 'FontSize', s.fs.label);
                title(ax, labels{k}, 'FontSize', s.fs.title);

                ax = subplot(3, 3, 3 + k); MWEC_Tuning_Plots.style_axes(ax, s); hold(ax, 'on'); grid(ax, 'on');
                Bk = bem.B_diag(k, :); Bp = Bk; Bp(Bp <= 0) = NaN;
                semilogy(ax, om, Bp, '-', 'Color', cols{k}, 'LineWidth', s.lw.primary);
                set(ax, 'YScale', 'log');
                xlabel(ax, '\omega (rad/s)', 'FontSize', s.fs.label);
                ylabel(ax, sprintf('B_{%d%d}', k, k), 'FontSize', s.fs.label);
            end

            % Row 3: surge-pitch coupling (normalised cross terms)
            ax = subplot(3, 3, 7); MWEC_Tuning_Plots.style_axes(ax, s); hold(ax, 'on'); grid(ax, 'on');
            plot(ax, om, bem.coupling.A15_norm, '-', 'Color', s.c.purple, 'LineWidth', s.lw.primary);
            yline(ax, 0, ':', 'Color', s.c.grey);
            xlabel(ax, '\omega (rad/s)', 'FontSize', s.fs.label);
            ylabel(ax, 'A_{15}/\surd(A_{11}A_{33})', 'FontSize', s.fs.label);
            title(ax, 'surge-pitch added-mass coupling', 'FontSize', s.fs.title);

            ax = subplot(3, 3, 8); MWEC_Tuning_Plots.style_axes(ax, s); hold(ax, 'on'); grid(ax, 'on');
            plot(ax, om, bem.coupling.B15_norm, '-', 'Color', s.c.red, 'LineWidth', s.lw.primary);
            yline(ax, 0, ':', 'Color', s.c.grey);
            xlabel(ax, '\omega (rad/s)', 'FontSize', s.fs.label);
            ylabel(ax, 'B_{15}/\surd(B_{11}B_{33})', 'FontSize', s.fs.label);
            title(ax, 'surge-pitch radiation coupling', 'FontSize', s.fs.title);

            ax = subplot(3, 3, 9); MWEC_Tuning_Plots.style_axes(ax, s); axis(ax, 'off');
            txt = {sprintf('max |A_{15}|/\\surd = %.3f', bem.coupling.A15_norm_max); ...
                   sprintf('max |B_{15}|/\\surd = %.3f', bem.coupling.B15_norm_max); ...
                   ''; '|cross|/\surd \approx 1  \Rightarrow'; 'shared dipole channel'; ...
                   '(surge & pitch compete)'};
            text(ax, 0.05, 0.9, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
                 'FontName', s.font, 'FontSize', s.fs.annot);

            sgt = sgtitle(sprintf('T2 — BEM at CG + surge-pitch coupling  %s  [%s]', region, ctag));
            set(sgt, 'FontName', s.font, 'FontSize', s.fs.sgtitle, 'FontWeight', 'bold');
            MWEC_Tuning_Plots.save_fig(fig, sprintf('fig_t2_bem_%s_%s', ctag, region), cfg);
        end

        %% =================================================================
        %%  T4 — Coupled absorbed-power density (replaces Lorentzian)
        %% =================================================================

        function fig_t4_placement(region, ctag, results, cfg)
            s = MWEC_Tuning_Plots.style();
            bem  = results.bem;
            om   = results.climate.omega;  S_ew = results.climate.S_ew;
            band = results.climate.band;
            omega_n = results.omega_n_target;  T_n = results.T_n_target_s;
            B_PTO   = results.closure.B_PTO_per_mode;
            rho = results.meta.rho;  g = results.meta.g_acc;

            R = MWEC_Tuning_Kernels.coupled_rao(om, omega_n, bem, B_PTO, 'freq');
            pd = zeros(numel(om), 3);                       % absorbed-power density [W/(rad/s)]
            for k = 1:3
                pd(:, k) = B_PTO(k) .* om.^2 .* abs(R.X(k, :)).'.^2 .* S_ew;
            end
            F_ew      = (rho * g^2 / 2) .* S_ew ./ max(om, 1e-10);
            ceil_dens = (3 * g ./ om.^2) .* F_ew;           % Falnes-limited power density

            labels = {'Surge', 'Heave', 'Pitch'};
            cols   = {s.mode.surge, s.mode.heave, s.mode.pitch};

            fig = MWEC_Tuning_Plots.new_figure(s, s.sz.full);

            ax1 = subplot(2, 1, 1); MWEC_Tuning_Plots.style_axes(ax1, s); hold(ax1, 'on'); grid(ax1, 'on');
            plot(ax1, om, ceil_dens, '-', 'Color', s.c.grey, 'LineWidth', s.lw.secondary, ...
                 'DisplayName', 'Falnes ceiling 3g/\omega^2 \cdot F_{ew}');
            h = area(ax1, om, pd, 'LineStyle', 'none');
            for k = 1:3, h(k).FaceColor = cols{k}; h(k).FaceAlpha = 0.6; h(k).DisplayName = labels{k}; end
            plot(ax1, om, sum(pd, 2), '-', 'Color', s.c.black, 'LineWidth', s.lw.secondary, 'DisplayName', 'total');
            for k = 1:3
                xline(ax1, omega_n(k), ':', 'Color', cols{k}, 'LineWidth', s.lw.ref, 'HandleVisibility', 'off');
            end
            xlabel(ax1, '\omega (rad/s)', 'FontSize', s.fs.label);
            ylabel(ax1, 'absorbed-power density  (W per rad/s)', 'FontSize', s.fs.label);
            title(ax1, 'Coupled absorbed-power density vs Falnes ceiling', 'FontSize', s.fs.title);
            legend(ax1, 'Location', 'northeast', 'FontSize', s.fs.legend);
            xlim(ax1, [band.omega_L*0.7, band.omega_H*1.2]);

            ax2 = subplot(2, 1, 2); MWEC_Tuning_Plots.style_axes(ax2, s); hold(ax2, 'on'); grid(ax2, 'on');
            yyaxis(ax2, 'left');
            plot(ax2, om, S_ew, '-', 'Color', s.c.blue, 'LineWidth', s.lw.primary);
            ylabel(ax2, 'S_{ew}(\omega)  (m^2 s/rad)', 'FontSize', s.fs.label);
            set(ax2, 'YColor', s.c.blue);
            yyaxis(ax2, 'right');
            for k = 1:3
                plot(ax2, om, abs(R.X(k, :)).', '--', 'Color', cols{k}, 'LineWidth', s.lw.secondary, ...
                     'DisplayName', sprintf('%s |X| (T_n=%.2fs)', labels{k}, T_n(k)));
            end
            ylabel(ax2, '|X_k(\omega)|  (coupled RAO)', 'FontSize', s.fs.label);
            set(ax2, 'YColor', s.c.black);
            for k = 1:3
                xline(ax2, omega_n(k), ':', 'Color', cols{k}, 'LineWidth', s.lw.ref, 'HandleVisibility', 'off');
            end
            xlabel(ax2, '\omega (rad/s)', 'FontSize', s.fs.label);
            title(ax2, 'Coupled RAO over the climate spectrum', 'FontSize', s.fs.title);
            xlim(ax2, [band.omega_L*0.7, band.omega_H*1.2]);

            sgt = sgtitle(sprintf('T4 — Placement %s [%s]  (\\eta_F=%.4f, <P_{abs}>=%.0f W, CWR=%.1f%%)', ...
                                  region, ctag, results.eta_Falnes, results.absorbed_power_W_total, results.CWR));
            set(sgt, 'FontName', s.font, 'FontSize', s.fs.sgtitle, 'FontWeight', 'bold');
            MWEC_Tuning_Plots.save_fig(fig, sprintf('fig_t4_placement_%s_%s', ctag, region), cfg);
        end

        %% =================================================================
        %%  T5 — Capture width vs Falnes ceiling (dipole combined)
        %% =================================================================

        function fig_t5_absorption(region, ctag, results, cfg)
            s = MWEC_Tuning_Plots.style();
            A = results.absorption_vs_T;  T = A.T;  Tn = results.T_n_target_s;
            fig = MWEC_Tuning_Plots.new_figure(s, s.sz.full);

            % Panel 1: heave (monopole) vs g/omega^2
            ax = subplot(2, 1, 1); MWEC_Tuning_Plots.style_axes(ax, s); hold(ax, 'on'); grid(ax, 'on');
            plot(ax, T, A.CW_heave, '-', 'Color', s.mode.heave, 'LineWidth', s.lw.primary, 'DisplayName', 'CW_{heave}');
            plot(ax, T, A.ceil_heave, '--', 'Color', s.c.grey, 'LineWidth', s.lw.ref, 'DisplayName', 'ceiling g/\omega^2 (monopole)');
            xline(ax, Tn(2), ':', 'Color', s.mode.heave, 'LineWidth', s.lw.secondary, ...
                  'Label', sprintf('T_n=%.2fs', Tn(2)), 'HandleVisibility', 'off');
            xlabel(ax, 'Period T (s)', 'FontSize', s.fs.label); ylabel(ax, 'CW (m)', 'FontSize', s.fs.label);
            title(ax, 'Heave — monopole channel', 'FontSize', s.fs.title);
            legend(ax, 'Location', 'best', 'FontSize', s.fs.legend); xlim(ax, [T(1), T(end)]);

            % Panel 2: surge + pitch + sum vs 2g/omega^2 (combined dipole)
            ax = subplot(2, 1, 2); MWEC_Tuning_Plots.style_axes(ax, s); hold(ax, 'on'); grid(ax, 'on');
            plot(ax, T, A.CW_surge, '-', 'Color', s.mode.surge, 'LineWidth', s.lw.primary, 'DisplayName', 'CW_{surge}');
            plot(ax, T, A.CW_pitch, '-', 'Color', s.mode.pitch, 'LineWidth', s.lw.primary, 'DisplayName', 'CW_{pitch}');
            plot(ax, T, A.CW_dipole_sum, '-', 'Color', s.c.black, 'LineWidth', s.lw.secondary, 'DisplayName', 'CW_{surge}+CW_{pitch}');
            plot(ax, T, A.ceil_dipole, '--', 'Color', s.c.grey, 'LineWidth', s.lw.ref, 'DisplayName', 'ceiling 2g/\omega^2 (dipole, combined)');
            xline(ax, Tn(1), ':', 'Color', s.mode.surge, 'LineWidth', s.lw.secondary, 'HandleVisibility', 'off');
            xline(ax, Tn(3), ':', 'Color', s.mode.pitch, 'LineWidth', s.lw.secondary, 'HandleVisibility', 'off');
            xlabel(ax, 'Period T (s)', 'FontSize', s.fs.label); ylabel(ax, 'CW (m)', 'FontSize', s.fs.label);
            title(ax, 'Surge + Pitch — shared dipole channel (ceiling counted ONCE)', 'FontSize', s.fs.title);
            legend(ax, 'Location', 'best', 'FontSize', s.fs.legend); xlim(ax, [T(1), T(end)]);

            sgt = sgtitle(sprintf('T5 — Capture width vs Falnes ceiling  %s  [%s]', region, ctag));
            set(sgt, 'FontName', s.font, 'FontSize', s.fs.sgtitle, 'FontWeight', 'bold');
            MWEC_Tuning_Plots.save_fig(fig, sprintf('fig_t5_absorption_%s_%s', ctag, region), cfg);
        end

        %% =================================================================
        %%  T6 — Climate scatter + placed periods
        %% =================================================================

        function fig_t6_scatter(region, ctag, results, climateGrid, cfg)
            s = MWEC_Tuning_Plots.style();
            fig = MWEC_Tuning_Plots.new_figure(s, s.sz.wide);
            ax = axes(fig); MWEC_Tuning_Plots.style_axes(ax, s); hold(ax, 'on');

            P = climateGrid.probability_grid; Hs_c = climateGrid.Hs_centers(:); Te_c = climateGrid.Te_centers(:);
            Ppct = P * 100; Ppct(Ppct == 0) = NaN;
            imagesc(ax, Te_c, Hs_c, Ppct); set(ax, 'YDir', 'normal'); colormap(ax, parula);
            cb = colorbar(ax); cb.Label.String = 'Probability  (% of records)'; cb.Label.FontSize = s.fs.label;

            labels = {'Surge', 'Heave', 'Pitch'}; cols = {s.mode.surge, s.mode.heave, s.mode.pitch};
            for k = 1:3
                xline(ax, results.T_n_target_s(k), '--', 'Color', cols{k}, 'LineWidth', s.lw.secondary, ...
                      'Label', sprintf('%s T_n=%.2fs', labels{k}, results.T_n_target_s(k)), ...
                      'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'top', 'HandleVisibility', 'off');
            end
            xline(ax, results.climate.band.T_L, ':', 'Color', s.c.black, 'LineWidth', s.lw.ref, 'Label', 'T_L', 'HandleVisibility', 'off');
            xline(ax, results.climate.band.T_H, ':', 'Color', s.c.black, 'LineWidth', s.lw.ref, 'Label', 'T_H', 'HandleVisibility', 'off');

            xlabel(ax, 'Energy period T_e  (s)', 'FontSize', s.fs.label);
            ylabel(ax, 'Significant wave height H_s  (m)', 'FontSize', s.fs.label);
            title(ax, sprintf('T6 — %s [%s] climate + design periods', region, ctag), 'FontSize', s.fs.title);
            MWEC_Tuning_Plots.save_fig(fig, sprintf('fig_t6_scatter_%s_%s', ctag, region), cfg, struct('pdf_content', 'image'));
        end

        %% =================================================================
        %%  C1 — Coupled vs decoupled (interference evidence)
        %% =================================================================

        function fig_c1_interference(region, ctag, results, cfg)
        %FIG_C1_INTERFERENCE  Absorbed-power density with the full 3x3 dynamics
        %   vs with the surge-pitch off-diagonals zeroed.  The difference is the
        %   surge-pitch interference — destructive where coupled < decoupled.
            s = MWEC_Tuning_Plots.style();
            bem  = results.bem;
            om   = results.climate.omega;  S_ew = results.climate.S_ew;  band = results.climate.band;
            omega_n = results.omega_n_target;  B_PTO = results.closure.B_PTO_per_mode;

            % Coupled
            Rc = MWEC_Tuning_Kernels.coupled_rao(om, omega_n, bem, B_PTO, 'freq');
            % Decoupled: zero off-diagonals of A,B
            bem_d = bem;
            for q = 1:size(bem_d.A_3DOF, 3)
                Ad = diag(diag(bem_d.A_3DOF(:, :, q))); bem_d.A_3DOF(:, :, q) = Ad;
                Bd = diag(diag(bem_d.B_3DOF(:, :, q))); bem_d.B_3DOF(:, :, q) = Bd;
            end
            Rd = MWEC_Tuning_Kernels.coupled_rao(om, omega_n, bem_d, B_PTO, 'freq');

            pc = zeros(numel(om), 1); pd = zeros(numel(om), 1);
            for k = 1:3
                pc = pc + B_PTO(k) .* om.^2 .* abs(Rc.X(k, :)).'.^2 .* S_ew;
                pd = pd + B_PTO(k) .* om.^2 .* abs(Rd.X(k, :)).'.^2 .* S_ew;
            end
            Pc = trapz(om, pc);  Pd = trapz(om, pd);

            fig = MWEC_Tuning_Plots.new_figure(s, s.sz.wide);
            ax = axes(fig); MWEC_Tuning_Plots.style_axes(ax, s); hold(ax, 'on'); grid(ax, 'on');
            plot(ax, om, pd, '--', 'Color', s.c.grey, 'LineWidth', s.lw.primary, ...
                 'DisplayName', sprintf('decoupled (A_{15}=B_{15}=0):  <P>=%.0f W', Pd));
            plot(ax, om, pc, '-', 'Color', s.c.red, 'LineWidth', s.lw.primary, ...
                 'DisplayName', sprintf('coupled (full 3x3):  <P>=%.0f W', Pc));
            xline(ax, omega_n(1), ':', 'Color', s.mode.surge, 'LineWidth', s.lw.ref, 'Label', 'surge', 'HandleVisibility', 'off');
            xline(ax, omega_n(3), ':', 'Color', s.mode.pitch, 'LineWidth', s.lw.ref, 'Label', 'pitch', 'HandleVisibility', 'off');

            xlabel(ax, '\omega (rad/s)', 'FontSize', s.fs.label);
            ylabel(ax, 'absorbed-power density  (W per rad/s)', 'FontSize', s.fs.label);
            dpct = 100 * (Pc / max(Pd, eps) - 1);
            title(ax, sprintf('C1 — surge-pitch interference  %s [%s]  (coupling effect %+.1f%%)', region, ctag, dpct), ...
                  'FontSize', s.fs.title);
            legend(ax, 'Location', 'northeast', 'FontSize', s.fs.legend);
            xlim(ax, [band.omega_L*0.7, band.omega_H*1.2]);
            MWEC_Tuning_Plots.save_fig(fig, sprintf('fig_c1_interference_%s_%s', ctag, region), cfg);
        end

        %% =================================================================
        %%  T7 — Cross-case / cross-pipeline
        %% =================================================================

        function fig_t7_cross(master_summary, cfg)
        %FIG_T7_CROSS  eta_Falnes and CWR across cases, both closures overlaid.
            s = MWEC_Tuning_Plots.style();
            tags = fieldnames(master_summary);
            tag_cols = {s.c.blue, s.c.red, s.c.green, s.c.purple};

            fig = MWEC_Tuning_Plots.new_figure(s, s.sz.wide);
            ax1 = subplot(1, 2, 1); MWEC_Tuning_Plots.style_axes(ax1, s); hold(ax1, 'on'); grid(ax1, 'on');
            ax2 = subplot(1, 2, 2); MWEC_Tuning_Plots.style_axes(ax2, s); hold(ax2, 'on'); grid(ax2, 'on');

            for it = 1:numel(tags)
                S = master_summary.(tags{it});
                eta = []; pabs = []; cwr = []; names = {};
                for j = 1:numel(S)
                    if ~isfield(S(j), 'region') || isempty(S(j).region), continue; end
                    eta(end+1)  = S(j).eta_Falnes;  %#ok<AGROW>
                    pabs(end+1) = S(j).P_abs_total_W; %#ok<AGROW>
                    cwr(end+1)  = S(j).CWR; %#ok<AGROW>
                    names{end+1} = S(j).region; %#ok<AGROW>
                end
                col = tag_cols{mod(it-1, numel(tag_cols)) + 1};
                scatter(ax1, eta, pabs, 70, 'MarkerFaceColor', col, 'MarkerEdgeColor', s.c.black, ...
                        'LineWidth', s.lw.ref, 'DisplayName', tags{it});
                scatter(ax2, eta, cwr, 70, 'MarkerFaceColor', col, 'MarkerEdgeColor', s.c.black, ...
                        'LineWidth', s.lw.ref, 'DisplayName', tags{it});
                for j = 1:numel(names)
                    text(ax1, eta(j), pabs(j), ['  ' names{j}], 'FontName', s.font, 'FontSize', s.fs.annot);
                    text(ax2, eta(j), cwr(j),  ['  ' names{j}], 'FontName', s.font, 'FontSize', s.fs.annot);
                end
            end

            xlabel(ax1, '\eta_{Falnes}', 'FontSize', s.fs.label);
            ylabel(ax1, '<P_{abs}>  (W)', 'FontSize', s.fs.label);
            title(ax1, '\eta_{Falnes} vs <P_{abs}>', 'FontSize', s.fs.title);
            legend(ax1, 'Location', 'best', 'FontSize', s.fs.legend);

            xlabel(ax2, '\eta_{Falnes}', 'FontSize', s.fs.label);
            ylabel(ax2, 'CWR  (%)', 'FontSize', s.fs.label);
            title(ax2, '\eta_{Falnes} vs capture width ratio', 'FontSize', s.fs.title);
            legend(ax2, 'Location', 'best', 'FontSize', s.fs.legend);

            sgt = sgtitle('T7 — Cross-case / cross-pipeline performance');
            set(sgt, 'FontName', s.font, 'FontSize', s.fs.sgtitle, 'FontWeight', 'bold');
            MWEC_Tuning_Plots.save_fig(fig, 'fig_t7_cross', cfg);
        end

        %% =================================================================
        %%  T8 — Evans vs Critical-damping (CWR + design-period shift)
        %% =================================================================

        function fig_t8_closure_compare(master_summary, cfg)
        %FIG_T8_CLOSURE_COMPARE  Direct Evans-vs-Critical-damping comparison.
        %   Left: capture width ratio by closure (the metric).  Right: the shift
        %   in the design periods (the deliverable) between the two closures.
            s = MWEC_Tuning_Plots.style();
            if ~(isfield(master_summary, 'evans') && isfield(master_summary, 'critical_damping'))
                return;   % both pipelines must have run
            end
            SE = master_summary.evans;  SC = master_summary.critical_damping;

            regions = {}; cwrE = []; cwrC = []; dT = [];
            for k = 1:numel(cfg.cases)
                reg = cfg.cases(k).region;
                if ~any(strcmp(cfg.cases_to_run, reg)), continue; end
                e = MWEC_Tuning_Plots.match_region(SE, reg);
                c = MWEC_Tuning_Plots.match_region(SC, reg);
                if isempty(e) || isempty(c), continue; end
                regions{end+1} = reg; %#ok<AGROW>
                cwrE(end+1) = e.CWR;  cwrC(end+1) = c.CWR; %#ok<AGROW>
                dT(end+1, :) = c.T_n_target_s(:).' - e.T_n_target_s(:).'; %#ok<AGROW>
            end
            if isempty(regions), return; end
            nR = numel(regions);

            fig = MWEC_Tuning_Plots.new_figure(s, s.sz.wide);

            % Panel 1: grouped CWR bars (Evans vs Critical-damping)
            ax1 = subplot(1, 2, 1); MWEC_Tuning_Plots.style_axes(ax1, s); hold(ax1, 'on'); grid(ax1, 'on');
            b = bar(ax1, 1:nR, [cwrE(:), cwrC(:)], 'grouped');
            b(1).FaceColor = s.c.blue;  b(1).DisplayName = 'Evans';
            b(2).FaceColor = s.c.red;   b(2).DisplayName = 'Critical-damping';
            set(ax1, 'XTick', 1:nR, 'XTickLabel', regions);
            ylabel(ax1, 'CWR  (%)', 'FontSize', s.fs.label);
            title(ax1, 'Capture width ratio by closure', 'FontSize', s.fs.title);
            legend(ax1, 'Location', 'best', 'FontSize', s.fs.legend);

            % Panel 2: design-period shift Crit - Evans, per mode (the deliverable)
            ax2 = subplot(1, 2, 2); MWEC_Tuning_Plots.style_axes(ax2, s); hold(ax2, 'on'); grid(ax2, 'on');
            b2 = bar(ax2, 1:nR, dT, 'grouped');
            mcols = {s.mode.surge, s.mode.heave, s.mode.pitch};
            mlbl  = {'\DeltaT_n surge', '\DeltaT_n heave', '\DeltaT_n pitch'};
            for m = 1:3, b2(m).FaceColor = mcols{m}; b2(m).DisplayName = mlbl{m}; end
            yline(ax2, 0, '-', 'Color', s.c.grey, 'LineWidth', s.lw.ref);
            set(ax2, 'XTick', 1:nR, 'XTickLabel', regions);
            ylabel(ax2, 'T_n^{crit} - T_n^{evans}  (s)', 'FontSize', s.fs.label);
            title(ax2, 'Design-period shift between closures', 'FontSize', s.fs.title);
            legend(ax2, 'Location', 'best', 'FontSize', s.fs.legend);

            sgt = sgtitle('T8 — Evans vs Critical-damping  (CWR + design-period shift)');
            set(sgt, 'FontName', s.font, 'FontSize', s.fs.sgtitle, 'FontWeight', 'bold');
            MWEC_Tuning_Plots.save_fig(fig, 'fig_t8_closure_compare', cfg);
        end

        function e = match_region(summary, region)
        %MATCH_REGION  Return the summary entry for a region, or [] if absent.
            e = [];
            for i = 1:numel(summary)
                if isfield(summary(i), 'region') && ~isempty(summary(i).region) && strcmp(summary(i).region, region)
                    e = summary(i);  return;
                end
            end
        end

    end
end
