function plot_complete_convergence(results, config)
%PLOT_COMPLETE_CONVERGENCE Draw the multi-panel Stage-1 convergence figure (surrogate R^2/MAPE and related panels) from results.stage1_2d and save it as a PNG.
    style = mwecmass.output.figures.presentation_style(config);
    try
        s1 = results.stage1_2d;
        n_iter = s1.iterations;

        fig = mwecmass.output.figures.new_figure(style, 'tall_double_column');
        fig.Name = 'WEC Optimization Convergence';
        t = tiledlayout(fig, 2, 2);

        %% PANEL A1: Surrogate Quality (R^2 & MAPE)
        ax1 = nexttile(t);
        hold(ax1, 'on');
        if n_iter >= 2 && ~isempty(s1.R2_mass)
            iters = 1:n_iter;
            valid = ~isnan(s1.R2_mass);

            if any(valid)
                yyaxis left
                h_r2m = plot(iters(valid), s1.R2_mass(valid), '-o', 'Color', style.color.series_a);
                mwecmass.output.figures.style_line(h_r2m, style, 'curve');
                if isfield(s1, 'R2_GM') && ~isempty(s1.R2_GM)
                    h_r2g = plot(iters(valid), s1.R2_GM(valid), '-s', 'Color', style.color.series_b);
                    mwecmass.output.figures.style_line(h_r2g, style, 'curve');
                end
                h_yl1 = ylabel('$R^2$ (Goodness of Fit)');
                mwecmass.output.figures.style_text(h_yl1, style, 'label');
                ylim([0, 1.05]);
                xl = xlim;
                h_thresh = plot(xl, [0.9, 0.9], '--', 'Color', style.fill_palette.boundary, ...
                                 'HandleVisibility', 'off');
                mwecmass.output.figures.style_line(h_thresh, style, 'reference');
                set(gca, 'YColor', style.fill_palette.boundary);

                yyaxis right
                if isfield(s1, 'MAPE_mass') && ~isempty(s1.MAPE_mass)
                    h_mm = plot(iters(valid), s1.MAPE_mass(valid), '--^', 'Color', style.color.series_a);
                    mwecmass.output.figures.style_line(h_mm, style, 'curve');
                end
                if isfield(s1, 'MAPE_GM') && ~isempty(s1.MAPE_GM)
                    h_mg = plot(iters(valid), s1.MAPE_GM(valid), '--v', 'Color', style.color.series_b);
                    mwecmass.output.figures.style_line(h_mg, style, 'curve');
                end
                h_yl2 = ylabel('MAPE (\%)');
                mwecmass.output.figures.style_text(h_yl2, style, 'label');
                set(gca, 'YColor', [0.5 0.5 0.5]);

                lg1 = legend('$R^2$(Mass)', '$R^2$(GM)', 'MAPE(Mass)', 'MAPE(GM)', ...
                       'Location', 'best');
                mwecmass.output.figures.style_legend(lg1, style);
            else
                h_ins = text(0.5, 0.5, 'R$^2$ / MAPE: insufficient data', ...
                     'Units', 'normalized', 'HorizontalAlignment', 'center');
                mwecmass.output.figures.style_text(h_ins, style, 'annotation');
            end
        else
            h_few = text(0.5, 0.5, sprintf('Stage 1: %d iteration(s) (need $\\geq 2$)', n_iter), ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center');
            mwecmass.output.figures.style_text(h_few, style, 'annotation');
        end
        xlabel('Stage 1 Iteration');
        title('(A1) Surrogate Quality');
        mwecmass.output.figures.apply_axes_style(ax1, style);

        %% PANEL A2: Raw 2D-3D Errors
        ax2 = nexttile(t);
        hold(ax2, 'on');
        if n_iter >= 1 && ~isempty(s1.mass_errors)
            iters = 1:n_iter;

            yyaxis left
            h_me = plot(iters, s1.mass_errors, '-o', 'Color', style.color.series_a);
            mwecmass.output.figures.style_line(h_me, style, 'curve');
            h_yl3 = ylabel('Mass Error (kg)');
            mwecmass.output.figures.style_text(h_yl3, style, 'label');
            set(gca, 'YColor', style.color.series_a);

            yyaxis right
            h_ge = plot(iters, s1.gm_errors, '-s', 'Color', style.color.series_b);
            mwecmass.output.figures.style_line(h_ge, style, 'curve');
            h_yl4 = ylabel('GM Error (m)');
            mwecmass.output.figures.style_text(h_yl4, style, 'label');
            set(gca, 'YColor', style.color.series_b);

            lg2 = legend('Mass Error', 'GM Error', 'Location', 'best');
            mwecmass.output.figures.style_legend(lg2, style);
        else
            h_noerr = text(0.5, 0.5, 'No Stage 1 error data', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center');
            mwecmass.output.figures.style_text(h_noerr, style, 'annotation');
        end
        xlabel('Stage 1 Iteration');
        title('(A2) 2D-3D Prediction Errors');
        mwecmass.output.figures.apply_axes_style(ax2, style);

        %% PANEL B1: Stage 2 Mass & GM
        ax3 = nexttile(t);
        plot_stage2_panel( ...
            results, {'mass', 'gm'}, {'Mass Balance', 'GM'}, ...
            '(B1) Stage 2: Mass \& GM Convergence', ax3, style);

        %% PANEL B2: Stage 2 Periods
        ax4 = nexttile(t);
        plot_stage2_panel( ...
            results, {'heave', 'pitch'}, {'Heave Period', 'Pitch Period'}, ...
            '(B2) Stage 2: Period Convergence', ax4, style);

        title(t, 'WEC Optimization Convergence');
        mwecmass.output.figures.apply_layout_style(t, style);
        saved = mwecmass.output.figures.export_figure( ...
            fig, 'WEC_Complete_Convergence', config);
        fprintf('  Convergence figure saved: %s\n', saved{1});

    catch ME
        warning('mwecmass:figures:CompleteConvergenceFailed', ...
                'Complete convergence plot failed: %s', ME.message);
    end
end

function plot_stage2_panel(results, fields, labels, panel_title, ax, style)
%PLOT_STAGE2_PANEL Plot a Stage-2 iteration-history panel from results.stage2_3d.iteration_errors.
    has_data = isfield(results, 'stage2_3d') && isfield(results.stage2_3d, 'iteration_errors');
    hold(ax, 'on');
    if has_data
        err = results.stage2_3d.iteration_errors;
        markers = {'-o', '-s'};
        colors = {style.color.series_a, style.color.series_b};
        for k = 1:length(fields)
            if isfield(err, fields{k}) && ~isempty(err.(fields{k}))
                vec = err.(fields{k});
                iters_3d = 0:length(vec)-1;
                h = plot(ax, iters_3d, vec, markers{k}, 'Color', colors{k});
                mwecmass.output.figures.style_line(h, style, 'curve');
            end
        end
        ylabel(ax, 'Normalised Error');
        lg = legend(ax, labels, 'Location', 'best');
        mwecmass.output.figures.style_legend(lg, style);
    else
        h_nodata = text(ax, 0.5, 0.5, 'Stage 2 data not available', ...
             'Units', 'normalized', 'HorizontalAlignment', 'center');
        mwecmass.output.figures.style_text(h_nodata, style, 'annotation');
    end
    xlabel(ax, '3D Optimizer Iteration');
    title(ax, panel_title);
    mwecmass.output.figures.apply_axes_style(ax, style);
    hold(ax, 'off');
end
