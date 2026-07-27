classdef WEC_Visualization
    % WEC_VISUALIZATION  Plotting and reporting suite for WEC optimisation.
    %
    %   All methods are static — no instance state. Every public method
    %   wraps its body in try-catch so a plotting failure never kills the
    %   optimiser. Refactored to enforce strict publication-quality formatting
    %   (Times New Roman, LaTeX interpreters, Cividis colormaps).
    %
    %   DENSITY VISUALISATION ARCHITECTURE
    %     Density nodes are centroids of N horizontal layers.  To render
    %     true uniform-colour bands (not a smooth gradient), each plot
    %     function must assign ONE colour per layer, not per vertex.
    %
    %     2D profile plots  : hull profile is clipped into N horizontal
    %                         slab polygons; each slab is filled with its
    %                         node density via patch() + 'flat'.
    %     3D mesh plot      : density is sampled at face CENTROIDS (not
    %                         vertices); FaceVertexCData is [nFaces×1] so
    %                         'flat' assigns one colour per triangular face.
    
    properties (Constant, Access = private)
        % Publication Style Constants
        FONT_NAME       = 'Times New Roman';
        FONT_SIZE_AXIS  = 12;
        FONT_SIZE_LABEL = 12;
        FONT_SIZE_TITLE = 12;
        LW_MAIN         = 2.0;
        LW_BOUNDARY     = 1.4;
        LW_AXES         = 0.6;
        WATERLINE_COLOR = [0.15, 0.55, 0.95];
        BOUNDARY_COLOR  = [0.10, 0.10, 0.10];
    end

    methods (Static)

        %% =============================================================
        %%  CONVERGENCE PLOTS
        %% =============================================================

        function plot_complete_convergence(results, config)
            try
                s1 = results.stage1_2d;
                n_iter = s1.iterations;

                figure('Name', 'WEC Optimization Convergence', ...
                       'Color', 'white', 'Position', [100, 100, 1400, 900]);

                %% PANEL A1: Surrogate Quality (R^2 & MAPE)
                ax1 = subplot(2, 2, 1);
                hold(ax1, 'on');
                if n_iter >= 2 && ~isempty(s1.R2_mass)
                    iters = 1:n_iter;
                    valid = ~isnan(s1.R2_mass);

                    if any(valid)
                        yyaxis left
                        plot(iters(valid), s1.R2_mass(valid), 'b-o', ...
                             'LineWidth', WEC_Visualization.LW_MAIN, 'MarkerSize', 6);
                        if isfield(s1, 'R2_GM') && ~isempty(s1.R2_GM)
                            plot(iters(valid), s1.R2_GM(valid), 'r-s', ...
                                 'LineWidth', WEC_Visualization.LW_MAIN, 'MarkerSize', 6);
                        end
                        ylabel('$R^2$ (Goodness of Fit)', 'Interpreter', 'latex', ...
                               'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                        ylim([0, 1.05]);
                        xl = xlim;
                        plot(xl, [0.9, 0.9], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
                        set(gca, 'YColor', 'k');

                        yyaxis right
                        if isfield(s1, 'MAPE_mass') && ~isempty(s1.MAPE_mass)
                            plot(iters(valid), s1.MAPE_mass(valid), 'b--^', ...
                                 'LineWidth', 1.5, 'MarkerSize', 5);
                        end
                        if isfield(s1, 'MAPE_GM') && ~isempty(s1.MAPE_GM)
                            plot(iters(valid), s1.MAPE_GM(valid), 'r--v', ...
                                 'LineWidth', 1.5, 'MarkerSize', 5);
                        end
                        ylabel('MAPE (\%)', 'Interpreter', 'latex', ...
                               'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                        set(gca, 'YColor', [0.5 0.5 0.5]);
                        
                        legend('$R^2$(Mass)', '$R^2$(GM)', 'MAPE(Mass)', 'MAPE(GM)', ...
                               'Location', 'best', 'Interpreter', 'latex', 'FontName', WEC_Visualization.FONT_NAME);
                    else
                        text(0.5, 0.5, 'R$^2$ / MAPE: insufficient data', ...
                             'Interpreter', 'latex', 'Units', 'normalized', 'HorizontalAlignment', 'center');
                    end
                else
                    text(0.5, 0.5, sprintf('Stage 1: %d iteration(s) (need $\\geq 2$)', n_iter), ...
                         'Interpreter', 'latex', 'Units', 'normalized', 'HorizontalAlignment', 'center');
                end
                xlabel('Stage 1 Iteration', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                title('(A1) Surrogate Quality', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_TITLE);
                WEC_Visualization.format_axis_publication(ax1);

                %% PANEL A2: Raw 2D-3D Errors
                ax2 = subplot(2, 2, 2);
                hold(ax2, 'on');
                if n_iter >= 1 && ~isempty(s1.mass_errors)
                    iters = 1:n_iter;

                    yyaxis left
                    plot(iters, s1.mass_errors, 'b-o', 'LineWidth', WEC_Visualization.LW_MAIN, 'MarkerSize', 6);
                    ylabel('Mass Error (kg)', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                    set(gca, 'YColor', 'b');

                    yyaxis right
                    plot(iters, s1.gm_errors, 'r-s', 'LineWidth', WEC_Visualization.LW_MAIN, 'MarkerSize', 6);
                    ylabel('GM Error (m)', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                    set(gca, 'YColor', 'r');

                    legend('Mass Error', 'GM Error', 'Location', 'best', 'Interpreter', 'latex', 'FontName', WEC_Visualization.FONT_NAME);
                else
                    text(0.5, 0.5, 'No Stage 1 error data', 'Interpreter', 'latex', ...
                         'Units', 'normalized', 'HorizontalAlignment', 'center');
                end
                xlabel('Stage 1 Iteration', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                title('(A2) 2D-3D Prediction Errors', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_TITLE);
                WEC_Visualization.format_axis_publication(ax2);

                %% PANEL B1: Stage 2 Mass & GM
                ax3 = subplot(2, 2, 3);
                WEC_Visualization.plot_stage2_panel( ...
                    results, {'mass', 'gm'}, {'Mass Balance', 'GM'}, ...
                    '(B1) Stage 2: Mass \& GM Convergence', ax3);

                %% PANEL B2: Stage 2 Periods
                ax4 = subplot(2, 2, 4);
                WEC_Visualization.plot_stage2_panel( ...
                    results, {'heave', 'pitch'}, {'Heave Period', 'Pitch Period'}, ...
                    '(B2) Stage 2: Period Convergence', ax4);

                sgtitle('WEC Optimization Convergence', ...
                        'Interpreter', 'latex', 'FontName', WEC_Visualization.FONT_NAME, 'FontSize', 14);
                timestamp = datestr(now, 'yyyymmdd_HHMMSS');
                fname = sprintf('WEC_Complete_Convergence_%s.png', timestamp);
                saveas(gcf, fname);
                fprintf('  Convergence figure saved: %s\n', fname);

            catch ME
                warning('WECVisualization:CompleteConvergenceFailed', ...
                        'Complete convergence plot failed: %s', ME.message);
            end
        end

        function plot_pid_learning_convergence(stage1_data, config)
            try
                n_iter = stage1_data.iterations;
                if n_iter < 2
                    return;
                end

                figure('Name', 'Stage 1: Surrogate Learning Convergence', ...
                       'Color', 'white', 'Position', [100, 100, 1600, 900]);
                iters = 1:n_iter;

                %% PANEL A: Surrogate Quality (R^2 & MAPE)
                axA = subplot(2, 2, 1);
                hold(axA, 'on'); 

                has_r2 = isfield(stage1_data, 'R2_mass') && ~isempty(stage1_data.R2_mass);
                if has_r2
                    valid = ~isnan(stage1_data.R2_mass);
                    if any(valid)
                        yyaxis left
                        plot(iters(valid), stage1_data.R2_mass(valid), 'b-o', ...
                             'LineWidth', WEC_Visualization.LW_MAIN, 'MarkerSize', 6, 'DisplayName', '$R^2$(Mass)');
                        if isfield(stage1_data, 'R2_GM') && ~isempty(stage1_data.R2_GM)
                            plot(iters(valid), stage1_data.R2_GM(valid), 'r-s', ...
                                 'LineWidth', WEC_Visualization.LW_MAIN, 'MarkerSize', 6, 'DisplayName', '$R^2$(GM)');
                        end
                        ylabel('$R^2$ (Goodness of Fit)', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                        ylim([0, 1.05]);
                        xl = xlim;
                        plot(xl, [0.9 0.9], 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
                        set(gca, 'YColor', 'k');

                        yyaxis right
                        if isfield(stage1_data, 'MAPE_mass') && ~isempty(stage1_data.MAPE_mass)
                            plot(iters(valid), stage1_data.MAPE_mass(valid), 'b--^', ...
                                 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'MAPE(Mass)');
                        end
                        if isfield(stage1_data, 'MAPE_GM') && ~isempty(stage1_data.MAPE_GM)
                            plot(iters(valid), stage1_data.MAPE_GM(valid), 'r--v', ...
                                 'LineWidth', 2, 'MarkerSize', 5, 'DisplayName', 'MAPE(GM)');
                        end
                        ylabel('MAPE (\%)', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                        set(gca, 'YColor', [0.5 0.5 0.5]);
                    end
                end
                xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                title('(A) Surrogate Quality: $R^2$ \& MAPE', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_TITLE);
                legend('Location', 'best', 'Interpreter', 'latex', 'FontName', WEC_Visualization.FONT_NAME);
                xlim([0.5, n_iter + 0.5]);
                WEC_Visualization.format_axis_publication(axA);

                %% PANEL B: Mass & GM Errors
                axB = subplot(2, 2, 2);
                hold(axB, 'on'); 

                yyaxis left
                plot(iters, stage1_data.mass_errors, 'b-o', ...
                     'LineWidth', WEC_Visualization.LW_MAIN, 'MarkerSize', 6, 'DisplayName', 'Mass Error (Raw)');

                if isfield(stage1_data, 'mass_2d_corrected_history') && ...
                   ~isempty(stage1_data.mass_2d_corrected_history) && ~isempty(stage1_data.mass_3d_history)
                    n_plot = min(length(stage1_data.mass_2d_corrected_history), length(stage1_data.mass_3d_history));
                    if n_plot > 0
                        corr_err = stage1_data.mass_2d_corrected_history(1:n_plot) - stage1_data.mass_3d_history(1:n_plot);
                        plot(1:n_plot, corr_err, 'g-^', 'LineWidth', WEC_Visualization.LW_MAIN, 'MarkerSize', 6, 'DisplayName', 'Mass Error (Corrected)');
                    end
                end
                ylabel('Mass Error (kg)', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                set(gca, 'YColor', 'b');

                yyaxis right
                plot(iters, stage1_data.gm_errors, 'r-s', ...
                     'LineWidth', WEC_Visualization.LW_MAIN, 'MarkerSize', 6, 'DisplayName', 'GM Error');
                ylabel('GM Error (m)', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                set(gca, 'YColor', 'r');

                xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                title('(B) 2D-3D Prediction Errors', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_TITLE);
                legend('Location', 'best', 'Interpreter', 'latex', 'FontName', WEC_Visualization.FONT_NAME);
                xlim([0.5, n_iter + 0.5]);
                WEC_Visualization.format_axis_publication(axB);

                %% PANEL C: Correction Factors
                axC = subplot(2, 2, 3);
                hold(axC, 'on'); 

                plot(iters, stage1_data.mass_corrections, 'b-o', 'LineWidth', WEC_Visualization.LW_MAIN, 'MarkerSize', 5, 'DisplayName', '$k_{vol}$');
                plot(iters, stage1_data.gm_corrections, 'r-s', 'LineWidth', WEC_Visualization.LW_MAIN, 'MarkerSize', 5, 'DisplayName', '$k_{GM}$');

                xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                ylabel('Correction Factor', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                title('(C) Bias Correction Convergence', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_TITLE);
                legend('Location', 'best', 'Interpreter', 'latex', 'FontName', WEC_Visualization.FONT_NAME);
                xlim([0.5, n_iter + 0.5]);
                WEC_Visualization.format_axis_publication(axC);

                %% PANEL D: Convergence Criteria
                axD = subplot(2, 2, 4);
                hold(axD, 'on'); 
                metrics = stage1_data.convergence_metrics;

                if isfield(metrics, 'constraints_satisfied') && length(metrics.constraints_satisfied) == n_iter
                    stairs(iters, metrics.constraints_satisfied, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Constraints OK');
                    stairs(iters, metrics.solution_stable, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Stable');
                    stairs(iters, metrics.mass_acceptable, 'm-', 'LineWidth', 2.5, 'DisplayName', 'Mass $<10\%$');
                    stairs(iters, metrics.converged, 'r-', 'LineWidth', 3.5, 'DisplayName', 'CONVERGED');

                    set(gca, 'YTick', [0, 1], 'YTickLabel', {'False', 'True'});
                    ylim([-0.15, 1.25]);
                end

                xlabel('Iteration', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                ylabel('Criterion Status', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                title('(D) Convergence Criteria', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_TITLE);
                legend('Location', 'best', 'Interpreter', 'latex', 'FontName', WEC_Visualization.FONT_NAME);
                xlim([0.5, n_iter + 0.5]);
                WEC_Visualization.format_axis_publication(axD);

                sgtitle('Stage 1: Surrogate Learning Convergence', 'Interpreter', 'latex', 'FontName', WEC_Visualization.FONT_NAME, 'FontSize', 14);
                timestamp = datestr(now, 'yyyymmdd_HHMMSS');
                fname = sprintf('WEC_Stage1_Convergence_%s.png', timestamp);
                saveas(gcf, fname);

            catch ME
                warning('WECVisualization:ConvergencePlotFailed', ...
                        'PID convergence plot failed: %s', ME.message);
            end
        end

        %% =============================================================
        %%  GEOMETRY VISUALISATION
        %% =============================================================

        function visualize_3D_design(props, config)
            % FIX: density computed at face CENTROIDS (not vertices).
            % FaceVertexCData is [nFaces×1] → 'flat' assigns one colour
            % per triangular face, producing sharp horizontal band edges.
            try
                fig = figure('Name', 'Final 3D WEC Design', 'Color', 'white', 'Position', [100, 100, 1000, 800]);
                ax = axes(fig);
                hold(ax, 'on');

                % FIX (V4): high-resolution visualization mesh, independent of BEM mesh.
                % config.mesh_Nu/Nv are BEM params — VIZ_N=60 gives smooth rendering.
                VIZ_N = 60;
                if isfield(config, 'viz_Nu'), VIZ_N = config.viz_Nu; end
                opts_viz = struct('trim_wl', false, 'close_gaps', true, ...
                                  'verbose', false, 'quarter_body', false);
                mesh_viz = WEC_Panelizer.generate( ...
                    config.ms2_model, props.vertical_shift, VIZ_N, VIZ_N, opts_viz);

                densities_at_nodes = [props.components.density];

                % Convert quad panels to triangles for patch()
                panels_q = mesh_viz.panels;
                tri_faces = [panels_q(:, [1 2 3]); panels_q(:, [1 3 4])];
                verts = mesh_viz.vertices;

                face_centroids_z = ( verts(tri_faces(:,1), 3) + ...
                                     verts(tri_faces(:,2), 3) + ...
                                     verts(tri_faces(:,3), 3) ) / 3;

                face_densities = interp1(config.density_nodes_z + props.vertical_shift, ...
                    densities_at_nodes, face_centroids_z, 'nearest', 'extrap');
                face_densities = max(config.ballast_density_bounds(1), ...
                                 min(config.ballast_density_bounds(2), face_densities));

                patch(ax, 'Faces', tri_faces, ...
                      'Vertices', verts, ...
                      'FaceVertexCData', face_densities, ...
                      'FaceColor', 'flat', ...
                      'EdgeColor', 'none');
                %
                % Draw horizontal strip boundary rings (mesh cross-sections)
                node_z_shifted = config.density_nodes_z + props.vertical_shift;
                dz_half     = (node_z_shifted(2) - node_z_shifted(1)) / 2;
                z_bounds_3d = [ node_z_shifted(1) - dz_half; ...
                               (node_z_shifted(1:end-1) + node_z_shifted(2:end)) / 2; ...
                                node_z_shifted(end) + dz_half ];
                z_bounds_3d = max(z_bounds_3d, min(verts(:,3)) + 1e-4);
                z_bounds_3d = min(z_bounds_3d, max(verts(:,3)) - 1e-4);
                for k = 1:length(z_bounds_3d)
                    ring = WEC_Visualization.mesh_slice_ring(tri_faces, verts, z_bounds_3d(k));
                    if ~isempty(ring)
                        plot3(ax, ring(:,1), ring(:,2), ring(:,3), '-', ...
                              'Color', [0.1, 0.1, 0.1], 'LineWidth', 1.2);
                    end
                end
                %
                colormap(ax, WEC_Visualization.cividis_map(256));
                cbar = colorbar(ax);
                cbar.Label.String = 'Density, \rho (kg/m^3)';
                cbar.Label.Interpreter = 'tex';
                cbar.TickLabelInterpreter = 'tex';

                min_d = min(face_densities);
                max_d = max(face_densities);
                if max_d > min_d
                    caxis(ax, [min_d, max_d]);
                else
                    caxis(ax, [min_d - 1, max_d + 1]);
                end

                xlims = [min(verts(:,1)), max(verts(:,1))];
                ylims = [min(verts(:,2)), max(verts(:,2))];
                [xg, yg] = meshgrid(linspace(xlims(1), xlims(2), 10), linspace(ylims(1), ylims(2), 10));
                
                h_water = surf(ax, xg, yg, zeros(size(xg)), ...
                    'FaceColor', WEC_Visualization.WATERLINE_COLOR, 'FaceAlpha', 0.4, 'EdgeColor', 'none');

                h_cg = plot3(ax, props.CG_total(1), props.CG_total(2), props.CG_total(3), ...
                    'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
                h_cb = plot3(ax, props.CB(1), props.CB(2), props.CB(3), ...
                    'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'b');

                axis(ax, 'equal'); view(ax, 30, 25);
                xlabel(ax, '$x$ [m]', 'Interpreter', 'latex'); 
                ylabel(ax, '$y$ [m]', 'Interpreter', 'latex'); 
                zlabel(ax, '$z$ [m]', 'Interpreter', 'latex');
                
                title(ax, sprintf('Final 3D Design (GM=%.3f m, $T_h$=%.2f s, $T_p$=%.2f s)', ...
                    props.GM_L, props.periods.heave, props.periods.pitch), 'Interpreter', 'latex');
                
                legend(ax, [h_water, h_cg, h_cb], {'Waterplane', 'CG', 'CB'}, ...
                       'Location', 'bestoutside', 'Interpreter', 'none', 'FontName', WEC_Visualization.FONT_NAME);
                
                WEC_Visualization.format_axis_publication(ax);
                hold(ax, 'off');

            catch ME
                warning('WEC_Visualization:3DPlotFailed', '3D visualization failed: %s', ME.message);
            end
        end

        function visualize_2d_design(props, config, plot_title)
            % FIX: hull profile is sliced into N horizontal slab polygons.
            % Each slab is filled with its node density via a separate
            % patch() call — this gives true uniform-colour bands.
            try
                fig = figure('Name', plot_title, 'Color', 'white', 'Position', [100, 100, 1000, 600]);

                ax = subplot(1, 2, 1);
                hold(ax, 'on'); axis(ax, 'equal'); 

                % FIX (V1): smooth parametric silhouette replaces config.profile.
                % extractProfileMS2 angular-sort fails for non-convex hulls (e.g. C0).
                shifted_profile = WEC_Visualization.build_smooth_viz_profile(config) ...
                                  + [0, props.vertical_shift];

                if isfield(props, 'densities_at_nodes')
                    cmap  = WEC_Visualization.cividis_map(256);
                    d_min = min(props.densities_at_nodes);
                    d_max = max(props.densities_at_nodes);
                    if d_max <= d_min, d_max = d_min + 1; end
                    node_z_wl = config.density_nodes_z + props.vertical_shift;

                    if ~isempty(config.shell)
                        % COMPOSITE: hatched shell annulus + density-coloured core
                        inner_prof = WEC_Visualization.compute_inner_profile_2d( ...
                            shifted_profile, config.shell_thickness);
                        WEC_Visualization.draw_composite_strips( ...
                            ax, shifted_profile, inner_prof, node_z_wl, ...
                            props.densities_at_nodes, [d_min, d_max]);
                        % Inner offset boundary as dashed spline
                        WEC_Visualization.plot_inner_spline(ax, inner_prof);
                    elseif isfield(config,'enable_constructability') && ...
                            config.enable_constructability && ...
                            isfield(config,'constructability_t_min')
                        % Density strips + dashed inner void boundary at t_min offset
                        if isfield(config, 'strip_edges') && ~isempty(config.strip_edges)
                            sb = config.strip_edges(:) + props.vertical_shift;
                        else
                            sb = [];
                        end
                        WEC_Visualization.draw_density_strips( ...
                            ax, shifted_profile, node_z_wl, ...
                            props.densities_at_nodes, cmap, [d_min, d_max], sb);
                        inner_prof = WEC_Visualization.compute_inner_profile_2d( ...
                            shifted_profile, config.constructability_t_min);
                        WEC_Visualization.plot_inner_spline(ax, inner_prof);
                        text(ax, max(shifted_profile(:,1)) + 0.05, ...
                             mean([min(shifted_profile(:,2)), max(shifted_profile(:,2))]), ...
                             sprintf('%.0f mm (3 in) min wall', config.constructability_t_min*1000), ...
                             'FontSize', 9, 'Color', [0.25 0.25 0.25], ...
                             'Rotation', 90, 'VerticalAlignment', 'bottom', 'Interpreter', 'none');
                    else
                        % Compute shifted strip bounds if available
                        if isfield(config, 'strip_edges') && ~isempty(config.strip_edges)
                            sb = config.strip_edges(:) + props.vertical_shift;
                        else
                            sb = [];
                        end
                        WEC_Visualization.draw_density_strips( ...
                            ax, shifted_profile, node_z_wl, ...
                            props.densities_at_nodes, cmap, [d_min, d_max], sb);
                    end

                    colormap(ax, cmap);
                    cbar = colorbar(ax);
                    cbar.Label.String = 'Density, \rho_{core} (kg/m^3)';
                    cbar.Label.Interpreter = 'tex';
                    cbar.TickLabelInterpreter = 'tex';
                    caxis(ax, [d_min, d_max]);
                else
                    plot(ax, shifted_profile(:,1), shifted_profile(:,2), 'k-', 'LineWidth', WEC_Visualization.LW_MAIN);
                end

                % Hull outline on top
                plot(ax, shifted_profile(:,1), shifted_profile(:,2), '-', ...
                     'Color', WEC_Visualization.BOUNDARY_COLOR, 'LineWidth', WEC_Visualization.LW_BOUNDARY);

                x_hull_range = [min(shifted_profile(:,1)) - 0.2, max(shifted_profile(:,1)) + 0.2];
                h_water = plot(ax, x_hull_range, [0 0], '--', 'Color', WEC_Visualization.WATERLINE_COLOR, 'LineWidth', WEC_Visualization.LW_MAIN);

                % Wall boundary (constructability mode)
                h_wall = [];
                if isfield(config, 'enable_constructability') && config.enable_constructability && ...
                        isfield(config, 'constructability_wall_height')
                    wall_z_wl = config.hull_z_max - config.constructability_wall_height + props.vertical_shift;
                    h_wall = plot(ax, x_hull_range, [wall_z_wl wall_z_wl], '-', ...
                        'Color', [0.8, 0.2, 0.1], 'LineWidth', WEC_Visualization.LW_MAIN);
                end

                h_cg = plot(ax, props.CG_total(1), props.CG_total(3), 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
                h_cb = plot(ax, props.CB(1), props.CB(3), 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'b');

                xlabel(ax, '$x$ [m]', 'Interpreter', 'latex'); 
                ylabel(ax, '$z$ [m]', 'Interpreter', 'latex');
                title(ax, 'WEC Design Profile', 'Interpreter', 'latex');
                handles = [h_water, h_cg, h_cb];
                labels  = {'Waterline', 'CG', 'CB'};
                if ~isempty(h_wall)
                    handles = [h_water, h_wall, h_cg, h_cb];
                    labels  = {'Waterline', 'Wall boundary', 'CG', 'CB'};
                end
                legend(ax, handles, labels, ...
                       'Location', 'bestoutside', 'Interpreter', 'none', 'FontName', WEC_Visualization.FONT_NAME);
                ylim(ax, [-5, 2]);
                WEC_Visualization.format_axis_publication(ax);
                hold(ax, 'off');

                % Properties card
                ax2 = subplot(1, 2, 2);
                axis(ax2, 'off');
                txt = { ...
                    'DESIGN PROPERTIES', ...
                    '==================', ...
                    sprintf('Draft: %.3f m', abs(config.hull_z_min + props.vertical_shift)), ...
                    sprintf('Total Mass: %.0f kg', props.mass_total), ...
                    sprintf('GM: %.3f m', props.GM_L), ...
                    '', ...
                    'NATURAL PERIODS:', ...
                    sprintf('  Heave: %.2f s', props.periods.heave), ...
                    sprintf('  Pitch: %.2f s', props.periods.pitch), ...
                    sprintf('  Surge: %.2f s', props.periods.surge), ...
                    '', ...
                    'BALLAST DENSITIES:'};

                if isfield(props, 'densities_at_nodes')
                    for i = 1:length(props.densities_at_nodes)
                        z_lev = config.density_nodes_z(i) + props.vertical_shift;
                        txt{end+1} = sprintf('  Node @ Z=%.2f m: %.0f kg/m3', ...
                            z_lev, props.densities_at_nodes(i)); %#ok<AGROW>
                    end
                end

                text(ax2, 0.05, 0.95, txt, 'Units', 'normalized', ...
                     'VerticalAlignment', 'top', 'FontName', 'Courier New', ...
                     'FontSize', WEC_Visualization.FONT_SIZE_AXIS - 1, 'Interpreter', 'none');
                 
                sgtitle(plot_title, 'Interpreter', 'latex', 'FontName', WEC_Visualization.FONT_NAME, 'FontSize', 14);

            catch ME
                warning('WEC_Visualization:2DPlotFailed', '2D visualization failed: %s', ME.message);
            end
        end

        function visualize_3d_cross_section(props_3d, config, x_opt)
            % FIX: 2D profile is sliced into N horizontal slab polygons.
            % Each slab is filled with its node density — true cake layers.
            try
                fig = figure('Color', 'white', 'Position', [100, 100, 1200, 600]);

                % Use realized vertical_shift if available (post-steel/UHPC solve),
                % falling back to the optimizer value so the function works standalone.
                if isfield(props_3d, 'vertical_shift') && ~isnan(props_3d.vertical_shift)
                    draft_final = props_3d.vertical_shift;
                else
                    draft_final = x_opt(1);
                end
                % Pick the density source: prefer the AS-BUILT per-strip
                % equivalent density set by the realisation solver, falling
                % back to the optimiser's continuous densities so the
                % function still works pre-realisation.
                if isfield(props_3d, 'realised_strip_density') && ...
                        ~isempty(props_3d.realised_strip_density)
                    densities_final = props_3d.realised_strip_density(:);
                    rho_source_is_realised = true;
                else
                    densities_final = x_opt(2:end);
                    rho_source_is_realised = false;
                end

                % FIX (V3): smooth parametric silhouette replaces config.profile.
                smooth_p3 = WEC_Visualization.build_smooth_viz_profile(config);
                px = smooth_p3(:, 1);
                pz = smooth_p3(:, 2) + draft_final;
                shifted_profile = [px, pz];

                node_z_shifted = config.density_nodes_z + draft_final;

                ax1 = subplot(1, 2, 1);
                hold(ax1, 'on');

                cmap  = WEC_Visualization.cividis_map(256);
                d_min = min(densities_final);
                d_max = max(densities_final);
                if d_max <= d_min, d_max = d_min + 1; end

                if ~isempty(config.shell)
                    % COMPOSITE: hatched shell annulus + density-coloured core
                    inner_prof = WEC_Visualization.compute_inner_profile_2d( ...
                        shifted_profile, config.shell_thickness);
                    WEC_Visualization.draw_composite_strips( ...
                        ax1, shifted_profile, inner_prof, node_z_shifted, ...
                        densities_final, [d_min, d_max]);
                    WEC_Visualization.plot_inner_spline(ax1, inner_prof);
                elseif isfield(config,'enable_constructability') && ...
                        config.enable_constructability && ...
                        isfield(config,'constructability_t_min')
                    % Density strips + dashed inner void boundary at t_min offset
                    if isfield(config, 'strip_edges') && ~isempty(config.strip_edges)
                        sb2 = config.strip_edges(:) + draft_final;
                    else
                        sb2 = [];
                    end
                    WEC_Visualization.draw_density_strips( ...
                        ax1, shifted_profile, node_z_shifted, densities_final, ...
                        cmap, [d_min, d_max], sb2);
                    inner_prof = WEC_Visualization.compute_inner_profile_2d( ...
                        shifted_profile, config.constructability_t_min);
                    WEC_Visualization.plot_inner_spline(ax1, inner_prof);
                else
                    % Compute shifted strip bounds if available
                    if isfield(config, 'strip_edges') && ~isempty(config.strip_edges)
                        sb2 = config.strip_edges(:) + draft_final;
                    else
                        sb2 = [];
                    end
                    WEC_Visualization.draw_density_strips( ...
                        ax1, shifted_profile, node_z_shifted, densities_final, ...
                        cmap, [d_min, d_max], sb2);
                end

                % Hull outline on top
                plot(ax1, px, pz, '-', ...
                     'Color', WEC_Visualization.BOUNDARY_COLOR, 'LineWidth', WEC_Visualization.LW_BOUNDARY);

                colormap(ax1, cmap);
                cb = colorbar(ax1);
                cb.Label.String = 'Density, \rho_{core} (kg/m^3)';
                cb.Label.Interpreter = 'tex';
                cb.TickLabelInterpreter = 'tex';
                caxis(ax1, [d_min, d_max]);

                % Waterline — extend to full hull profile width
                x_hull_range = [min(px) - 0.2, max(px) + 0.2];
                h_wl = plot(ax1, x_hull_range, [0 0], '--', 'Color', WEC_Visualization.WATERLINE_COLOR, 'LineWidth', WEC_Visualization.LW_MAIN);

                % Wall boundary (constructability mode)
                h_wall = [];
                if isfield(config, 'enable_constructability') && config.enable_constructability && ...
                        isfield(config, 'constructability_wall_height')
                    wall_z_wl = config.hull_z_max - config.constructability_wall_height + draft_final;
                    h_wall = plot(ax1, x_hull_range, [wall_z_wl wall_z_wl], '-', ...
                        'Color', [0.8, 0.2, 0.1], 'LineWidth', WEC_Visualization.LW_MAIN);
                end

                h_cg = plot(ax1, props_3d.CG_total(1), props_3d.CG_total(3), 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 2);
                h_cb = plot(ax1, props_3d.CB(1), props_3d.CB(3), 'bs', 'MarkerSize', 12, 'MarkerFaceColor', 'b', 'LineWidth', 2);
                
                xlabel(ax1, '$x$ [m]', 'Interpreter', 'latex'); 
                ylabel(ax1, '$z$ [m]', 'Interpreter', 'latex');
                title(ax1, 'WEC Design Profile', 'Interpreter', 'latex');
                if ~isempty(h_wall)
                    legend(ax1, [h_wl, h_wall, h_cg, h_cb], {'Waterline', 'Wall boundary', 'CG', 'CB'}, ...
                        'Location', 'best', 'Interpreter', 'none', 'FontName', WEC_Visualization.FONT_NAME);
                else
                    legend(ax1, [h_wl, h_cg, h_cb], {'Waterline', 'CG', 'CB'}, ...
                        'Location', 'best', 'Interpreter', 'none', 'FontName', WEC_Visualization.FONT_NAME);
                end
                axis(ax1, 'equal'); 
                WEC_Visualization.format_axis_publication(ax1);
                hold(ax1, 'off');

                % Properties card
                ax2 = subplot(1, 2, 2);
                axis(ax2, 'off');
                ty = 0.95;
                lh = 0.04;
                FN = WEC_Visualization.FONT_NAME;
                FS = WEC_Visualization.FONT_SIZE_AXIS;

                text(ax2, 0.05, ty, 'DESIGN PROPERTIES', 'FontSize', FS+1, 'FontName', FN, 'FontWeight', 'bold', 'Interpreter', 'none');
                ty = ty - 1.5*lh;
                text(ax2, 0.05, ty, sprintf('Draft: %.3f m', abs(config.hull_z_min + draft_final)), 'FontSize', FS, 'FontName', FN, 'Interpreter', 'none');
                ty = ty - lh;
                text(ax2, 0.05, ty, sprintf('Total Mass: %.0f kg', props_3d.mass_total), 'FontSize', FS, 'FontName', FN, 'Interpreter', 'none');
                ty = ty - lh;
                text(ax2, 0.05, ty, sprintf('GM: %.3f m', props_3d.GM_L), 'FontSize', FS, 'FontName', FN, 'Interpreter', 'none');
                ty = ty - 2*lh;
                text(ax2, 0.05, ty, 'NATURAL PERIODS:', 'FontSize', FS, 'FontName', FN, 'FontWeight', 'bold', 'Interpreter', 'none');
                ty = ty - lh;
                text(ax2, 0.05, ty, sprintf('  Heave: %.2f s', props_3d.periods.heave), 'FontSize', FS, 'FontName', FN, 'Interpreter', 'none');
                ty = ty - lh;
                text(ax2, 0.05, ty, sprintf('  Pitch: %.2f s', props_3d.periods.pitch), 'FontSize', FS, 'FontName', FN, 'Interpreter', 'none');
                ty = ty - lh;
                if isinf(props_3d.periods.surge)
                    text(ax2, 0.05, ty, '  Surge: Inf s', 'FontSize', FS, 'FontName', FN, 'Interpreter', 'none');
                else
                    text(ax2, 0.05, ty, sprintf('  Surge: %.2f s', props_3d.periods.surge), 'FontSize', FS, 'FontName', FN, 'Interpreter', 'none');
                end
                ty = ty - 2*lh;
                text(ax2, 0.05, ty, 'BALLAST DENSITIES:', 'FontSize', FS, 'FontName', FN, 'FontWeight', 'bold', 'Interpreter', 'none');
                ty = ty - lh;
                for i = 1:length(densities_final)
                    z_node = config.density_nodes_z(i) + draft_final;
                    text(ax2, 0.05, ty, sprintf('  Node %d @ Z=%.2f m: rho=%.0f kg/m3', ...
                        i, z_node, densities_final(i)), 'FontSize', FS-1, 'FontName', 'Courier New', 'Interpreter', 'none');
                    ty = ty - lh;
                end

                if rho_source_is_realised
                    if isfield(props_3d, 'realisation_mode') && ...
                            strcmp(props_3d.realisation_mode, 'uhpc_fill')
                        sg_str = 'Realised UHPC + Void Hull (as-built strip densities)';
                        png_name = 'WEC_Final_3D_CrossSection_UHPC.png';
                    else
                        sg_str = 'Realised Steel-Fill Hull (as-built strip densities)';
                        png_name = 'WEC_Final_3D_CrossSection_Steel.png';
                    end
                else
                    sg_str = 'Stage 2 Final: Optimised 3D Design (continuous ballast densities)';
                    png_name = 'WEC_Final_3D_CrossSection.png';
                end
                sgtitle(sg_str, 'Interpreter', 'latex', ...
                        'FontName', WEC_Visualization.FONT_NAME, 'FontSize', 14);

                saveas(gcf, png_name);
            catch ME
                warning('WEC_Visualization:CrossSectionFailed', '3D cross-section failed: %s', ME.message);
            end
        end

        %% =============================================================
        %%  EQUIVALENT DENSITY PLOTS
        %% =============================================================

        function visualize_2d_equivalent(props, config, x_opt)
        % VISUALIZE_2D_EQUIVALENT  2D cross-section coloured by strip-average
        %   equivalent uniform density.
            try
                % Use realized vertical_shift if available (post-steel/UHPC solve),
                % falling back to the optimizer value so the function works standalone.
                if isfield(props, 'vertical_shift') && ~isnan(props.vertical_shift)
                    draft_final = props.vertical_shift;
                else
                    draft_final = x_opt(1);
                end
                densities_core  = x_opt(2:end);

                fig_name = 'WEC Equivalent Density (2D)';
                if isfield(props, 'realised_strip_density') && ...
                        ~isempty(props.realised_strip_density)
                    rho_source_is_realised = true;
                    if isfield(props, 'realisation_mode') && ...
                            strcmp(props.realisation_mode, 'uhpc_fill')
                        fig_name = 'WEC Equivalent Density (2D) — Realised UHPC';
                    else
                        fig_name = 'WEC Equivalent Density (2D) — Realised Steel';
                    end
                else
                    rho_source_is_realised = false;
                end

                fig = figure('Name', fig_name, ...
                             'Color', 'white', 'Position', [150, 150, 1000, 600]); %#ok<NASGU>

                ax = subplot(1, 2, 1);
                hold(ax, 'on'); axis(ax, 'equal');

                % FIX (V2): smooth parametric silhouette replaces config.profile.
                smooth_p2 = WEC_Visualization.build_smooth_viz_profile(config);
                px = smooth_p2(:, 1);
                pz = smooth_p2(:, 2) + draft_final;
                shifted_profile = [px, pz];
                node_z_wl = config.density_nodes_z + draft_final;

                % --- Compute equivalent densities ---
                %  Pre-realisation : derive from optimiser's per-node rho_core
                %  Post-realisation: use the as-built strip rho_eff directly
                if rho_source_is_realised
                    rho_eq = props.realised_strip_density(:);
                    if length(rho_eq) ~= length(densities_core)
                        % Shape mismatch — fall back to optimiser-derived
                        rho_eq = WEC_Visualization.compute_strip_equivalent_density( ...
                                     config, densities_core);
                        rho_source_is_realised = false;
                    end
                else
                    rho_eq = WEC_Visualization.compute_strip_equivalent_density( ...
                                 config, densities_core);
                end

                cmap  = WEC_Visualization.cividis_map(256);
                d_min = min(rho_eq);
                d_max = max(rho_eq);
                if d_max <= d_min, d_max = d_min + 1; end

                % Compute shifted strip bounds if available
                if isfield(config, 'strip_edges') && ~isempty(config.strip_edges)
                    sb3 = config.strip_edges(:) + draft_final;
                else
                    sb3 = [];
                end

                % Fill outer hull with equivalent density
                WEC_Visualization.draw_density_strips( ...
                    ax, shifted_profile, node_z_wl, rho_eq, cmap, [d_min, d_max], sb3);

                % Overlay inner offset spline to show shell boundary
                if ~isempty(config.shell)
                    inner_prof = WEC_Visualization.compute_inner_profile_2d( ...
                        shifted_profile, config.shell_thickness);
                    WEC_Visualization.plot_inner_spline(ax, inner_prof);
                end

                % Outer hull boundary
                plot(ax, px, pz, '-', ...
                     'Color', WEC_Visualization.BOUNDARY_COLOR, ...
                     'LineWidth', WEC_Visualization.LW_BOUNDARY);

                x_hull_range = [min(px) - 0.2, max(px) + 0.2];
                h_wl  = plot(ax, x_hull_range, [0 0], '--', ...
                             'Color', WEC_Visualization.WATERLINE_COLOR, ...
                             'LineWidth', WEC_Visualization.LW_MAIN);

                % Wall boundary (constructability mode)
                h_wall_eq = [];
                if isfield(config, 'enable_constructability') && config.enable_constructability && ...
                        isfield(config, 'constructability_wall_height')
                    wall_z_wl = config.hull_z_max - config.constructability_wall_height + draft_final;
                    h_wall_eq = plot(ax, x_hull_range, [wall_z_wl wall_z_wl], '-', ...
                        'Color', [0.8, 0.2, 0.1], 'LineWidth', WEC_Visualization.LW_MAIN);
                end

                h_cg  = plot(ax, props.CG_total(1), props.CG_total(3), 'ro', ...
                             'MarkerSize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 2);
                h_cb  = plot(ax, props.CB(1), props.CB(3), 'bs', ...
                             'MarkerSize', 12, 'MarkerFaceColor', 'b', 'LineWidth', 2);

                colormap(ax, cmap);
                cbar = colorbar(ax);
                cbar.Label.String = 'Density, \rho_{eq} (kg/m^3)';
                cbar.Label.Interpreter = 'tex';
                cbar.TickLabelInterpreter = 'tex';
                caxis(ax, [d_min, d_max]);

                xlabel(ax, '$x$ [m]', 'Interpreter', 'latex');
                ylabel(ax, '$z$ [m]', 'Interpreter', 'latex');
                title(ax, 'Equivalent Uniform Density per Strip', 'Interpreter', 'latex');
                handles_eq = [h_wl, h_cg, h_cb];
                labels_eq  = {'Waterline', 'CG', 'CB'};
                if ~isempty(h_wall_eq)
                    handles_eq = [h_wl, h_wall_eq, h_cg, h_cb];
                    labels_eq  = {'Waterline', 'Wall boundary', 'CG', 'CB'};
                end
                legend(ax, handles_eq, labels_eq, ...
                       'Location', 'bestoutside', 'Interpreter', 'none', ...
                       'FontName', WEC_Visualization.FONT_NAME);
                WEC_Visualization.format_axis_publication(ax);
                hold(ax, 'off');

                % --- Properties card ---
                ax2 = subplot(1, 2, 2);
                axis(ax2, 'off');
                ty = 0.95; lh = 0.045;
                FN = WEC_Visualization.FONT_NAME;
                FS = WEC_Visualization.FONT_SIZE_AXIS;
                text(ax2, 0.05, ty, 'EQUIVALENT DENSITY TABLE', ...
                     'FontSize', FS+1, 'FontName', FN, ...
                     'FontWeight', 'bold', 'Interpreter', 'none');
                ty = ty - 1.5*lh;
                text(ax2, 0.05, ty, sprintf('Draft: %.3f m', abs(config.hull_z_min + draft_final)), ...
                     'FontSize', FS, 'FontName', FN, 'Interpreter', 'none');
                ty = ty - lh;
                text(ax2, 0.05, ty, sprintf('T_heave: %.2f s  |  T_pitch: %.2f s', ...
                     props.periods.heave, props.periods.pitch), ...
                     'FontSize', FS, 'FontName', FN, 'Interpreter', 'none');
                ty = ty - 2*lh;
                text(ax2, 0.05, ty, ...
                     sprintf('%-5s  %-8s  %-10s  %-10s', ...
                             'Strip', 'Z [m]', 'rho_core', 'rho_eq'), ...
                     'FontSize', FS-1, 'FontName', 'Courier New', 'FontWeight', 'bold', 'Interpreter', 'none');
                ty = ty - lh;
                for i = 1:length(densities_core)
                    z_node = config.density_nodes_z(i) + draft_final;
                    text(ax2, 0.05, ty, ...
                         sprintf('  %-3d  %+6.3f m  %6.0f kg/m3  %6.0f kg/m3', ...
                                 i, z_node, densities_core(i), rho_eq(i)), ...
                         'FontSize', FS-1, 'FontName', 'Courier New', 'Interpreter', 'none');
                    ty = ty - lh;
                    if ty < 0.05, break; end
                end

                if rho_source_is_realised
                    if isfield(props, 'realisation_mode') && ...
                            strcmp(props.realisation_mode, 'uhpc_fill')
                        sg2_str = 'Realised UHPC + Void: as-built strip $\rho_{eq}$';
                    else
                        sg2_str = 'Realised Steel-Fill: as-built strip $\rho_{eq}$';
                    end
                else
                    sg2_str = 'Steiner-Equivalent Bulk Density (optimiser)';
                end
                sgtitle(sg2_str, 'Interpreter', 'latex', ...
                        'FontName', WEC_Visualization.FONT_NAME, 'FontSize', 14);

            catch ME
                warning('WEC_Visualization:2DEquivFailed', ...
                        '2D equivalent visualization failed: %s', ME.message);
            end
        end


        function visualize_3D_equivalent(props, config)
        % VISUALIZE_3D_EQUIVALENT  3D hull coloured by strip-average equivalent
        %   uniform density.
        %
        %   If props.realised_strip_density is present (post steel/UHPC
        %   realisation), the hull is coloured by the AS-BUILT effective
        %   density per strip; otherwise it falls back to the optimiser's
        %   continuous-density per-strip equivalent.
            try
                % Density source: realised vs optimiser
                if isfield(props, 'realised_strip_density') && ...
                        ~isempty(props.realised_strip_density)
                    rho_eq = props.realised_strip_density(:);
                    rho_source_is_realised = true;
                else
                    densities_at_nodes = [props.components.density];
                    rho_eq = WEC_Visualization.compute_strip_equivalent_density( ...
                                 config, densities_at_nodes);
                    rho_source_is_realised = false;
                end

                if rho_source_is_realised
                    if isfield(props, 'realisation_mode') && ...
                            strcmp(props.realisation_mode, 'uhpc_fill')
                        fig_name = 'WEC Equivalent Density (3D) — Realised UHPC';
                    else
                        fig_name = 'WEC Equivalent Density (3D) — Realised Steel';
                    end
                else
                    fig_name = 'WEC Equivalent Density (3D)';
                end
                fig = figure('Name', fig_name, ...
                             'Color', 'white', 'Position', [150, 150, 1000, 800]); %#ok<NASGU>
                ax = axes(fig);
                hold(ax, 'on');

                % FIX (V5): high-resolution visualization mesh, independent of BEM mesh.
                VIZ_N_eq = 60;
                if isfield(config, 'viz_Nu'), VIZ_N_eq = config.viz_Nu; end
                opts_viz = struct('trim_wl', false, 'close_gaps', true, ...
                                  'verbose', false, 'quarter_body', false);

                mesh_viz = WEC_Panelizer.generate( ...
                    config.ms2_model, props.vertical_shift, VIZ_N_eq, VIZ_N_eq, opts_viz);

                verts = mesh_viz.vertices;
                panels = mesh_viz.panels;

                % Convert quad panels to triangles for patch()
                tri_faces = [panels(:, [1 2 3]); panels(:, [1 3 4])];

                % Face centroid z-coordinates
                face_z = (verts(tri_faces(:,1), 3) + ...
                          verts(tri_faces(:,2), 3) + ...
                          verts(tri_faces(:,3), 3)) / 3;

                % Strip boundaries in world frame
                node_z_wl = config.density_nodes_z + props.vertical_shift;
                if isfield(config, 'strip_edges') && ~isempty(config.strip_edges)
                    strip_edges_wl = config.strip_edges + props.vertical_shift;
                else
                    dz_half = (node_z_wl(2) - node_z_wl(1)) / 2;
                    strip_edges_wl = [node_z_wl(1) - dz_half; ...
                                     (node_z_wl(1:end-1) + node_z_wl(2:end)) / 2; ...
                                      node_z_wl(end) + dz_half];
                end

                % Assign rho_eq to each face
                N = length(rho_eq);
                face_rho = zeros(size(tri_faces, 1), 1);
                for k = 1:N
                    mask = (face_z >= strip_edges_wl(k)) & ...
                           (face_z <  strip_edges_wl(min(k+1, length(strip_edges_wl))));
                    face_rho(mask) = rho_eq(k);
                end
                % Faces outside all strips → nearest strip value
                unassigned = (face_rho == 0);
                if any(unassigned)
                    face_rho(unassigned) = interp1(node_z_wl, rho_eq, ...
                        face_z(unassigned), 'nearest', 'extrap');
                end

                patch(ax, 'Faces', tri_faces, 'Vertices', verts, ...
                      'FaceVertexCData', face_rho, ...
                      'FaceColor', 'flat', 'EdgeColor', 'none');

                colormap(ax, WEC_Visualization.cividis_map(256));
                cbar = colorbar(ax);
                cbar.Label.String = 'Density, \rho_{eq} (kg/m^3)';
                cbar.Label.Interpreter = 'tex';
                cbar.TickLabelInterpreter = 'tex';
                d_min = min(rho_eq); d_max = max(rho_eq);
                if d_max > d_min, caxis(ax, [d_min, d_max]);
                else,             caxis(ax, [d_min-1, d_max+1]); end

                % Waterplane
                xlims = [min(verts(:,1)), max(verts(:,1))];
                ylims = [min(verts(:,2)), max(verts(:,2))];
                [xg, yg] = meshgrid(linspace(xlims(1), xlims(2), 10), ...
                                    linspace(ylims(1), ylims(2), 10));
                h_water = surf(ax, xg, yg, zeros(size(xg)), ...
                    'FaceColor', WEC_Visualization.WATERLINE_COLOR, ...
                    'FaceAlpha', 0.4, 'EdgeColor', 'none');
                h_cg = plot3(ax, props.CG_total(1), props.CG_total(2), props.CG_total(3), ...
                    'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
                h_cb = plot3(ax, props.CB(1), props.CB(2), props.CB(3), ...
                    'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'b');

                axis(ax, 'equal'); view(ax, 30, 25);
                xlabel(ax, '$x$ [m]', 'Interpreter', 'latex');
                ylabel(ax, '$y$ [m]', 'Interpreter', 'latex');
                zlabel(ax, '$z$ [m]', 'Interpreter', 'latex');
                if rho_source_is_realised
                    if isfield(props, 'realisation_mode') && ...
                            strcmp(props.realisation_mode, 'uhpc_fill')
                        title_str = sprintf( ...
                            'Realised UHPC + Void: as-built $\\rho_{eq}$ (GM=%.3f m, $T_h$=%.2f s)', ...
                            props.GM_L, props.periods.heave);
                    else
                        title_str = sprintf( ...
                            'Realised Steel-Fill: as-built $\\rho_{eq}$ (GM=%.3f m, $T_h$=%.2f s)', ...
                            props.GM_L, props.periods.heave);
                    end
                else
                    title_str = sprintf( ...
                        'Equivalent Bulk Density — Optimiser (GM=%.3f m, $T_h$=%.2f s)', ...
                        props.GM_L, props.periods.heave);
                end
                title(ax, title_str, 'Interpreter', 'latex');
                legend(ax, [h_water, h_cg, h_cb], {'Waterplane', 'CG', 'CB'}, ...
                       'Location', 'bestoutside', 'Interpreter', 'none', ...
                       'FontName', WEC_Visualization.FONT_NAME);
                WEC_Visualization.format_axis_publication(ax);
                hold(ax, 'off');

            catch ME
                warning('WEC_Visualization:3DEquivFailed', ...
                        '3D equivalent visualization failed: %s', ME.message);
            end
        end

        %% =============================================================
        %%  CONSOLE REPORTING & VALIDATION
        %% =============================================================

        function reportFinalResults(props, config)
            try
                fprintf('\n=== FINAL 3D DESIGN RESULTS ===\n\n');
                fprintf('  Optimal Draft: %.3f m\n', abs(config.hull_z_min + props.vertical_shift));
                fprintf('  Vertical Shift: %.3f m\n', props.vertical_shift);
                if isfield(props, 'components') && ~isempty(props.components)
                    fprintf('\n  Density Profile:\n');
                    for i = 1:length(props.components)
                        fprintf('    Node %d @ Z = %6.2f m: rho = %7.0f kg/m^3\n', ...
                            i, props.components(i).z_level, props.components(i).density);
                    end
                end
                fprintf('\n  Mass Properties:\n');
                fprintf('    Total Mass:       %10.2f kg\n', props.mass_total);
                fprintf('    Buoyant Force:    %10.2f kg\n', props.mass_buoyant_force);
                fprintf('    Mass Discrepancy: %10.2e kg\n', props.mass_discrepancy);
                fprintf('    CG: [%.3f, %.3f, %.3f] m\n', props.CG_total);
                fprintf('\n  Hydrostatic Properties:\n');
                fprintf('    Submerged Volume: %.4f m^3\n', props.V_sub);
                fprintf('    Waterplane Area:  %.4f m^2\n', props.Aw);
                fprintf('    CB: [%.3f, %.3f, %.3f] m\n', props.CB);
                fprintf('    GM: %.3f m\n', props.GM_L);
                fprintf('\n  Natural Periods:\n');
                fprintf('    Heave: %.2f s (Target: %.2f s)\n', props.periods.heave, config.T_heave_goal);
                fprintf('    Pitch: %.2f s (Target: %.2f s)\n', props.periods.pitch, config.T_pitch_goal);
                if isinf(props.periods.surge)
                    fprintf('    Surge: Inf s\n');
                else
                    fprintf('    Surge: %.2f s\n', props.periods.surge);
                end
                fprintf('\n');
            catch ME
                warning('WEC_Visualization:ReportFailed', 'Results reporting failed: %s', ME.message);
            end
        end

        function report_stage2_summary(stage2_data)
            fprintf('\n=== STAGE 2 SUMMARY ===\n\n');
            qm = stage2_data.quality_metrics;
            fprintf('  Convergence: %d\n', stage2_data.converged);
            fprintf('\n  Constraints:\n');
            fprintf('    Monotonic density: %d\n', qm.monotonic);
            fprintf('    Mass balance: %d (error: %.2e kg)\n', qm.mass_balance, qm.mass_balance_error_kg);
            fprintf('    GM constraint: %d (margin: %.3f m)\n', qm.GM_satisfied, qm.GM_margin);
            fprintf('\n  Solver:\n');
            fprintf('    Exit flag: %d\n', qm.exitflag);
            fprintf('    Optimal: %d\n', qm.fmincon_optimal);
            fprintf('    Constraint violation: %.2e\n', qm.constrviolation);
            fprintf('    First-order optimality: %.2e\n\n', qm.firstorderopt);
        end

        function [converged, quality_metrics] = check_3d_convergence(x_opt_3d, config, exitflag, output)
            props_3d  = calculate_3d_properties(x_opt_3d, config);
            densities = x_opt_3d(2:end);

            % OPT-MONO FIX: the optimizer's constraint function correctly
            % skips the wall-platform boundary (config.wall_strip_index).
            % Without the same exclusion here, the pinned wall strip (always
            % 2500 kg/m³) always creates a positive diff against the platform
            % strip below it (e.g. 732 → 2500), making monotonic_ok = false
            % even when the design is entirely valid.
            if isfield(config, 'wall_strip_index') && ~isempty(config.wall_strip_index)
                w       = config.wall_strip_index;
                N_dens  = length(densities);
                % Build the same pair list the constraint function uses:
                % skip any pair whose lower or upper index is the wall strip.
                pairs_to_check = [];
                for ii = 1:(N_dens - 1)
                    if ii == w || (ii + 1) == w
                        continue;   % skip wall-platform boundary
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

            % OPT-EXITFLAG FIX: fmincon exit codes for acceptable convergence:
            %   1 → First-order optimality satisfied (true optimal)
            %   2 → Step size below StepTolerance, constraints satisfied
            %       (local minimum; common outcome for SQP on nonsmooth problems)
            %   0 → Maximum iterations reached (borderline — accept only if
            %       constraint violation and first-order residual are small)
            % The original code only accepted exitflag == 0 for "acceptable",
            % missing exitflag == 2 entirely, and used a 1e-4 first-order
            % threshold that is too tight for the step-tolerance exit path
            % (typical value ~3-10e-3).
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


        %% =============================================================
        %%  HYDRODYNAMIC COEFFICIENT & RAO PLOTS
        %% =============================================================

        function plot_hydrodynamics(hydro_cache, vertical_shift, final_props, config)

            try

            FN  = WEC_Visualization.FONT_NAME;
            % FIX (D3): Removed FSA and LWA — unused variables in plot_hydrodynamics.
            FSL = WEC_Visualization.FONT_SIZE_LABEL;
            FST = WEC_Visualization.FONT_SIZE_TITLE;
            LW  = WEC_Visualization.LW_MAIN;

            c_A   = [0.12 0.47 0.71; 0.20 0.63 0.17; 0.89 0.10 0.11];
            c_B   = [0.40 0.65 0.85; 0.55 0.78 0.52; 0.95 0.50 0.50];
            c_o_A = [0.58 0.40 0.74; 1.00 0.50 0.05; 0.55 0.34 0.29];
            c_o_B = [0.75 0.60 0.85; 1.00 0.73 0.47; 0.73 0.55 0.50];

            dof6     = [1 3 5];
            dof_num  = {'1','3','5'};
            dof_name = {'Surge','Heave','Pitch'};

            % ── Find closest cached draft ─────────────────────────────
            [~, idx] = min(abs(hydro_cache.drafts - vertical_shift));
            fprintf('  Plotting hydrodynamics at vs = %+.4f m (cache idx %d)\n', ...
                    hydro_cache.drafts(idx), idx);

            A_full = hydro_cache.A{idx};
            B_full = hydro_cache.B{idx};
            omega  = hydro_cache.omega;

            if isempty(A_full) || isempty(omega)
                warning('WECVisualization:NoHydroData', ...
                        'No frequency-dependent data at vs=%+.4f.', vertical_shift);
                return;
            end

            T  = 2*pi ./ omega;
            Nf = length(omega);
            A_3x3 = A_full(dof6, dof6, :);
            B_3x3 = B_full(dof6, dof6, :);
            A_inf  = hydro_cache.A_inf{idx}(dof6, dof6);
            % FIX (L6): guard against T_band missing in older cache files.
            if isfield(hydro_cache, 'T_band') && ~isempty(hydro_cache.T_band)
                T_band = hydro_cache.T_band;
            else
                T_band = [4, 25];   % safe default: 4–25 s operational band
                warning('WEC_Visualization:NoTBand', ...
                    'hydro_cache.T_band missing — using default [4, 25] s.');
            end

            % ══════════════════════════════════════════════════════════
            %  FIGURE 1: A(w) and B(w) vs Period
            % ══════════════════════════════════════════════════════════
            fig1 = figure('Name', 'Hydrodynamic Coefficients', ...
                          'Color', 'white', 'Position', [50 50 1600 900]);

            A_units = {'[kg]','[kg]','[kg$\cdot$m$^2$]'};
            B_units = {'[kg/s]','[kg/s]','[kg$\cdot$m$^2$/s]'};

            for d = 1:3
                ax = subplot(2, 3, d);  hold(ax, 'on');
                Av = squeeze(A_3x3(d,d,:));
                Bv = squeeze(B_3x3(d,d,:));
                Ai = A_inf(d,d);

                yyaxis(ax, 'left');
                plot(ax, T, Av, '-', 'Color', c_A(d,:), 'LineWidth', LW);
                yline(ax, Ai, '--', 'Color', c_A(d,:)*0.5, 'LineWidth', 1.2, ...
                      'HandleVisibility', 'off');
                ylabel(ax, sprintf('$A_{%s%s}$ %s', dof_num{d}, dof_num{d}, A_units{d}), ...
                       'Interpreter', 'latex', 'FontSize', FSL);
                set(ax, 'YColor', c_A(d,:));

                yyaxis(ax, 'right');
                plot(ax, T, Bv, '-', 'Color', c_B(d,:), 'LineWidth', LW);
                ylabel(ax, sprintf('$B_{%s%s}$ %s', dof_num{d}, dof_num{d}, B_units{d}), ...
                       'Interpreter', 'latex', 'FontSize', FSL);
                set(ax, 'YColor', c_B(d,:));

                WEC_Visualization.shade_T_band(ax, T_band);
                xlabel(ax, '$T$ [s]', 'Interpreter', 'latex', 'FontSize', FSL);
                title(ax, sprintf('%s -- $A_\\infty$ = %.0f', dof_name{d}, Ai), ...
                      'Interpreter', 'latex', 'FontSize', FST);
                legend(ax, {'$A(\omega)$','$B(\omega)$'}, ...
                       'Interpreter', 'latex', 'FontName', FN, 'Location', 'best');
                WEC_Visualization.format_axis_publication(ax);
            end

            pairs = [1 2; 1 3; 2 3];
            for k = 1:3
                ax = subplot(2, 3, 3+k);  hold(ax, 'on');
                ii = pairs(k,1);  jj = pairs(k,2);
                Av = squeeze(A_3x3(ii,jj,:));
                Bv = squeeze(B_3x3(ii,jj,:));
                Ai = A_inf(ii,jj);

                yyaxis(ax, 'left');
                plot(ax, T, Av, '-', 'Color', c_o_A(k,:), 'LineWidth', LW);
                yline(ax, Ai, '--', 'Color', c_o_A(k,:)*0.5, 'LineWidth', 1.2, ...
                      'HandleVisibility', 'off');
                ylabel(ax, sprintf('$A_{%s%s}$', dof_num{ii}, dof_num{jj}), ...
                       'Interpreter', 'latex', 'FontSize', FSL);
                set(ax, 'YColor', c_o_A(k,:));

                yyaxis(ax, 'right');
                plot(ax, T, Bv, '-', 'Color', c_o_B(k,:), 'LineWidth', LW);
                ylabel(ax, sprintf('$B_{%s%s}$', dof_num{ii}, dof_num{jj}), ...
                       'Interpreter', 'latex', 'FontSize', FSL);
                set(ax, 'YColor', c_o_B(k,:));

                WEC_Visualization.shade_T_band(ax, T_band);
                xlabel(ax, '$T$ [s]', 'Interpreter', 'latex', 'FontSize', FSL);
                title(ax, sprintf('%s--%s -- $A_\\infty$ = %.1f', ...
                      dof_name{ii}, dof_name{jj}, Ai), ...
                      'Interpreter', 'latex', 'FontSize', FST);
                legend(ax, {'$A(\omega)$','$B(\omega)$'}, ...
                       'Interpreter', 'latex', 'FontName', FN, 'Location', 'best');
                WEC_Visualization.format_axis_publication(ax);
            end

            sgtitle(fig1, sprintf( ...
                'Hydrodynamic Coefficients at Origin ($v_s$ = %+.3f m)', ...
                vertical_shift), 'Interpreter', 'latex', 'FontName', FN, 'FontSize', 14);

            fname1 = sprintf('WEC_HydroCoeffs_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
            exportgraphics(fig1, fname1, 'Resolution', 300);
            fprintf('  Figure saved: %s\n', fname1);

            % ══════════════════════════════════════════════════════════
            %  FIGURE 2: RAO
            % ══════════════════════════════════════════════════════════
            has_Fe = isfield(hydro_cache, 'Fe') && length(hydro_cache.Fe) >= idx ...
                     && ~isempty(hydro_cache.Fe{idx});

            fig2 = figure('Name', 'Response Amplitude Operators', ...
                          'Color', 'white', 'Position', [100 50 1400 800]);

            mass   = final_props.mass_total;
            Iyy_cg = final_props.Iyy;
            z_G    = final_props.CG_total(3);
            rho_w  = config.RHO_WATER;
            g_acc  = config.G;

            Iyy_O = Iyy_cg + mass * z_G^2;

            M_O = [mass,       0,    mass*z_G;
                   0,          mass, 0;
                   mass*z_G,   0,    Iyy_O   ];

            assert(Iyy_O >= Iyy_cg - 1e-6, ...
                'Steiner produced Iyy_O < Iyy_CG — sign bug in z_G');

            C33 = rho_w * g_acc * final_props.Aw;
            z_B = final_props.CB(3);
            C55 = rho_w*g_acc*final_props.I_wp_yy ...
                + rho_w*g_acc*final_props.V_sub*z_B ...
                - mass*g_acc*z_G;
            C55 = max(C55, 0);

            C55_check = mass * g_acc * final_props.GM_L;
            if abs(C55 - C55_check) / max(abs(C55), 1) > 0.01
                fprintf('    WARNING: C55 formula (%.1f) vs mg·GM (%.1f) differ by %.1f%%\n', ...
                        C55, C55_check, abs(C55-C55_check)/max(abs(C55),1)*100);
            end

            C_O = diag([0, C33, C55]);

            if isfield(config, 'B_visc_diag') && ~isempty(config.B_visc_diag)
                B_visc = diag(config.B_visc_diag);
            else
                B_visc = zeros(3);
            end
            has_visc = any(diag(B_visc) > 0);

            if has_Fe
                Fe_3xM = hydro_cache.Fe{idx}(dof6, :);
            end

            RAO_complex = zeros(3, Nf);
            for m_idx = 1:Nf
                w = omega(m_idx);
                A_o = 0.5*(A_3x3(:,:,m_idx) + A_3x3(:,:,m_idx)');
                B_o = 0.5*(B_3x3(:,:,m_idx) + B_3x3(:,:,m_idx)');
                Z = -w^2*(M_O + A_o) + 1i*w*(B_o + B_visc) + C_O;

                if has_Fe
                    if abs(det(Z)) > 1e-30
                        X_O = Z \ Fe_3xM(:, m_idx);
                        RAO_complex(:, m_idx) = [X_O(1) + z_G * X_O(3);
                                                 X_O(2);
                                                 X_O(3)];
                    end
                else
                    if abs(det(Z)) > 1e-30
                        H = Z \ eye(3);
                        RAO_complex(1, m_idx) = H(1,1) + z_G * H(3,1);
                        RAO_complex(2, m_idx) = H(2,2);
                        RAO_complex(3, m_idx) = H(3,3);
                    end
                end
            end

            RAO   = abs(RAO_complex);
            PHASE = angle(RAO_complex) * (180/pi);

            if has_Fe
                RAO(3,:) = RAO(3,:) * (180/pi);
            end

            if has_Fe
                rao_titles = {'Surge RAO at CG $|X_1^{CG}/A|$', ...
                              'Heave RAO $|X_3/A|$', ...
                              'Pitch RAO $|X_5/A|$'};
                rao_ylabels = {'$|X_1^{CG}/A|$ [m/m]', ...
                               '$|X_3/A|$ [m/m]', ...
                               '$|X_5/A|$ [deg/m]'};
                phase_titles = {'Surge phase $\angle X_1^{CG}$', ...
                                'Heave phase $\angle X_3$', ...
                                'Pitch phase $\angle X_5$'};
            else
                rao_titles = {'Surge $|H_{11}^{CG}|$ (no $F_e$)', ...
                              'Heave $|H_{33}|$ (no $F_e$)', ...
                              'Pitch $|H_{55}|$ (no $F_e$)'};
                rao_ylabels = {'$|H_{11}^{CG}|$ [m/N]', ...
                               '$|H_{33}|$ [m/N]', ...
                               '$|H_{55}|$ [rad/(N$\cdot$m)]'};
                phase_titles = {'Surge phase $\angle H_{11}^{CG}$', ...
                                'Heave phase $\angle H_{33}$', ...
                                'Pitch phase $\angle H_{55}$'};
            end

            T_natural = [final_props.periods.surge, final_props.periods.heave, final_props.periods.pitch];
            T_nat_lbl = {'$T_{surge}$', '$T_{heave}$', '$T_{pitch}$'};

            for d = 1:3
                ax = subplot(2, 3, d);  hold(ax, 'on');
                plot(ax, T, RAO(d,:), '-', 'Color', c_A(d,:), 'LineWidth', LW);
                if isfinite(T_natural(d)) && T_natural(d) > 0 && T_natural(d) < max(T)*1.5
                    xline(ax, T_natural(d), ':', 'Color', [0.3 0.3 0.3], ...
                          'LineWidth', 1.5, ...
                          'Label', sprintf('%s = %.1f s', T_nat_lbl{d}, T_natural(d)), ...
                          'Interpreter', 'latex', 'FontSize', 10, ...
                          'LabelVerticalAlignment', 'top');
                end
                xlabel(ax, '$T$ [s]', 'Interpreter', 'latex', 'FontSize', FSL);
                ylabel(ax, rao_ylabels{d}, 'Interpreter', 'latex', 'FontSize', FSL);
                title(ax, rao_titles{d}, 'Interpreter', 'latex', 'FontSize', FST);
                WEC_Visualization.format_axis_publication(ax);
            end

            for d = 1:3
                ax = subplot(2, 3, d + 3);  hold(ax, 'on');
                plot(ax, T, PHASE(d,:), '-', 'Color', c_A(d,:), 'LineWidth', LW);
                if isfinite(T_natural(d)) && T_natural(d) > 0 && T_natural(d) < max(T)*1.5
                    xline(ax, T_natural(d), ':', 'Color', [0.3 0.3 0.3], ...
                          'LineWidth', 1.5);
                end
                yline(ax, 0, '-', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
                yline(ax, -90, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
                yline(ax, 90, '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.5);
                ylim(ax, [-180 180]);
                set(ax, 'YTick', [-180 -90 0 90 180]);
                xlabel(ax, '$T$ [s]', 'Interpreter', 'latex', 'FontSize', FSL);
                ylabel(ax, 'Phase [deg]', 'Interpreter', 'latex', 'FontSize', FSL);
                title(ax, phase_titles{d}, 'Interpreter', 'latex', 'FontSize', FST);
                WEC_Visualization.format_axis_publication(ax);
            end

            if has_visc
                visc_str = sprintf('$B_{visc}$ = [%.0f, %.0f, %.0f]', ...
                    config.B_visc_diag(1), config.B_visc_diag(2), config.B_visc_diag(3));
            else
                visc_str = 'Potential flow only ($B_{visc} = 0$)';
            end

            if has_Fe, sgt = 'RAO (from HAMS excitation force)';
            else,      sgt = 'Frequency Response (no $F_e$ data)';
            end
            sgtitle(fig2, sprintf('%s -- $v_s$ = %+.3f m, $m$ = %.0f kg -- %s', ...
                    sgt, vertical_shift, mass, visc_str), ...
                    'Interpreter', 'latex', 'FontName', FN, 'FontSize', 13);

            fname2 = sprintf('WEC_RAO_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
            exportgraphics(fig2, fname2, 'Resolution', 300);
            fprintf('  Figure saved: %s\n', fname2);

            % ══════════════════════════════════════════════════════════
            %  FIGURE 3: A_inf diagonals vs vertical_shift
            % ══════════════════════════════════════════════════════════
            N_cache = length(hydro_cache.drafts);
            if N_cache < 2
                fprintf('  Skipping A_inf vs draft plot (need >= 2 cache entries)\n');
                return;
            end

            drafts_all = hydro_cache.drafts(:);
            A_inf_diag = zeros(N_cache, 3);
            valid = true(N_cache, 1);
            for kk = 1:N_cache
                Ak = hydro_cache.A_inf{kk};
                if isempty(Ak) || max(abs(Ak(:))) < 1e-6
                    valid(kk) = false; continue;
                end
                A_inf_diag(kk,:) = [Ak(1,1), Ak(3,3), Ak(5,5)];
            end

            drafts_v = drafts_all(valid);
            A_vals_v = A_inf_diag(valid,:);
            if size(A_vals_v,1) < 2; return; end

            fig3 = figure('Name', 'A_inf vs Draft', ...
                          'Color', 'white', 'Position', [150 150 1400 450]);
            [drafts_s, si] = sort(drafts_v);
            A_s = A_vals_v(si,:);

            A_labels = {'$A_{11}^{\infty}$ [kg]', '$A_{33}^{\infty}$ [kg]', ...
                        '$A_{55}^{\infty}$ [kg$\cdot$m$^2$]'};
            for d = 1:3
                ax = subplot(1, 3, d);  hold(ax, 'on');
                plot(ax, drafts_s, A_s(:,d), 'o', 'Color', c_A(d,:), ...
                     'MarkerFaceColor', c_A(d,:), 'MarkerSize', 8, 'LineWidth', 1.5);

                p = polyfit(drafts_s, A_s(:,d), 1);
                xf = linspace(min(drafts_s), max(drafts_s), 100);
                plot(ax, xf, polyval(p, xf), '--', 'Color', c_A(d,:)*0.6, 'LineWidth', 1.5);

                A_pred = polyval(p, drafts_s);
                SS_res = sum((A_s(:,d) - A_pred).^2);
                SS_tot = sum((A_s(:,d) - mean(A_s(:,d))).^2);
                if SS_tot > 1e-12; R2 = 1 - SS_res/SS_tot; else; R2 = 1; end

                xline(ax, vertical_shift, ':', 'Color', [0.3 0.3 0.3], ...
                      'LineWidth', 1.5, 'Label', '$v_s^*$', ...
                      'Interpreter', 'latex', 'FontSize', 10);

                xlabel(ax, 'Vertical Shift [m]', 'Interpreter', 'latex', 'FontSize', FSL);
                ylabel(ax, A_labels{d}, 'Interpreter', 'latex', 'FontSize', FSL);
                title(ax, sprintf('%s -- $R^2$ = %.4f (slope = %.1f/m)', ...
                      dof_name{d}, R2, p(1)), 'Interpreter', 'latex', 'FontSize', FST);
                legend(ax, {'HAMS', 'Linear fit'}, ...
                       'Interpreter', 'latex', 'FontName', FN, 'Location', 'best');
                WEC_Visualization.format_axis_publication(ax);
            end

            sgtitle(fig3, sprintf( ...
                '$A^{\\infty}$ vs Vertical Shift (%d points -- linearity check)', ...
                length(drafts_v)), 'Interpreter', 'latex', 'FontName', FN, 'FontSize', 14);

            fname3 = sprintf('WEC_Ainf_vs_Draft_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
            exportgraphics(fig3, fname3, 'Resolution', 300);
            fprintf('  Figure saved: %s\n', fname3);

            catch ME
                warning('WECVisualization:HydroPlotFailed', ...
                        'Hydrodynamic plot failed: %s', ME.message);
            end
        end


        function plot_draft_landscape(sweep, config)
            vs   = sweep.vs;
            N    = length(vs);
            feas = sweep.feasible;

            T_h = nan(1, N);  T_p = nan(1, N);
            GM  = nan(1, N);  mass_err = nan(1, N);

            for k = 1:N
                if ~isempty(sweep.props{k}) && isfield(sweep.props{k}, 'GM_L')
                    p = sweep.props{k};
                    T_h(k)  = p.periods.heave;
                    T_p(k)  = p.periods.pitch;
                    GM(k)   = p.GM_L;
                    mass_err(k) = abs(p.mass_total - p.mass_buoyant_force) ...
                                  / max(p.mass_total, 1) * 100;
                end
            end

            fig = figure('Name', 'Draft Landscape', ...
                         'Color', 'w', 'Position', [100 100 1100 750]);

            ax1 = subplot(2, 2, 1);  hold(ax1, 'on');  grid(ax1, 'on');
            if any(feas)
                plot(ax1, vs(feas), sweep.fval(feas), 'go', ...
                     'MarkerSize', 8, 'LineWidth', 1.5);
            end
            if any(~feas)
                plot(ax1, vs(~feas), sweep.fval(~feas), 'rx', ...
                     'MarkerSize', 8, 'LineWidth', 1.5);
            end
            fi = find(feas);
            if ~isempty(fi)
                [~, best_in_feas] = min(sweep.fval(fi));
                best = fi(best_in_feas);
                plot(ax1, vs(best), sweep.fval(best), 'p', ...
                     'MarkerSize', 15, 'MarkerFaceColor', 'g', ...
                     'MarkerEdgeColor', 'k');
            end
            xlabel(ax1, 'Vertical shift [m]');
            ylabel(ax1, 'Objective $f^*$', 'Interpreter', 'latex');
            title(ax1, 'Objective landscape');
            legend(ax1, 'Feasible', 'Infeasible', 'Best', 'Location', 'best');

            ax2 = subplot(2, 2, 2);  hold(ax2, 'on');  grid(ax2, 'on');
            plot(ax2, vs, T_h, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
            plot(ax2, vs, T_p, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 4);
            yline(ax2, config.T_heave_goal, 'b--', 'LineWidth', 1);
            yline(ax2, config.T_pitch_goal, 'r--', 'LineWidth', 1);
            xl = [vs(1), vs(end)];
            patch(ax2, [xl(1) xl(2) xl(2) xl(1)], ...
                  [config.T_heave_range(1) config.T_heave_range(1) ...
                   config.T_heave_range(2) config.T_heave_range(2)], ...
                  'b', 'FaceAlpha', 0.08, 'EdgeColor', 'none');
            patch(ax2, [xl(1) xl(2) xl(2) xl(1)], ...
                  [config.T_pitch_range(1) config.T_pitch_range(1) ...
                   config.T_pitch_range(2) config.T_pitch_range(2)], ...
                  'r', 'FaceAlpha', 0.08, 'EdgeColor', 'none');
            xlabel(ax2, 'Vertical shift [m]');
            ylabel(ax2, 'Period [s]');
            title(ax2, 'Natural periods vs draft');
            legend(ax2, '$T_{heave}$', '$T_{pitch}$', ...
                   'Interpreter', 'latex', 'Location', 'best');

            ax3 = subplot(2, 2, 3);  hold(ax3, 'on');  grid(ax3, 'on');
            plot(ax3, vs, GM, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 4);
            yline(ax3, config.gm_min, 'r--', 'LineWidth', 1.5, ...
                  'Label', 'GM_{min}');
            patch(ax3, [xl(1) xl(2) xl(2) xl(1)], ...
                  [config.gm_range(1) config.gm_range(1) ...
                   config.gm_range(2) config.gm_range(2)], ...
                  'g', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
            xlabel(ax3, 'Vertical shift [m]');
            ylabel(ax3, 'GM [m]');
            title(ax3, 'Metacentric height vs draft');

            ax4 = subplot(2, 2, 4);  hold(ax4, 'on');  grid(ax4, 'on');
            plot(ax4, vs, mass_err, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 4);
            yline(ax4, 1, 'r--', 'LineWidth', 1, 'Label', '1%');
            xlabel(ax4, 'Vertical shift [m]');
            ylabel(ax4, 'Mass error [%]');
            title(ax4, 'Mass balance error vs draft');

            sgtitle(fig, sprintf('Draft Landscape Sweep (%d drafts)', N), ...
                    'FontSize', 14);

            fname = sprintf('WEC_DraftLandscape_%s.png', ...
                            datestr(now, 'yyyymmdd_HHMMSS'));
            exportgraphics(fig, fname, 'Resolution', 300);
            fprintf('  Draft landscape saved: %s\n', fname);
        end



        %% =============================================================
        %%  PANEL NORMALS DIAGNOSTIC (new — user request)
        %% =============================================================

        function plot_panel_normals(mesh, options)
        % PLOT_PANEL_NORMALS  3D hull mesh with outward panel normals as quivers.
        %
        %   WEC_Visualization.plot_panel_normals(mesh)
        %   WEC_Visualization.plot_panel_normals(mesh, options)
        %
        %   Renders the BEM panel mesh with an outward-pointing quiver arrow at
        %   each panel centroid.  Normal convention follows HAMS-MREL (cross
        %   product of panel diagonals, outward):
        %     Quad:     n̂ = normalise( cross(V3−V1, V4−V2) )
        %     Triangle: n̂ = normalise( cross(V2−V1, V3−V2) )  (V4 = V3)
        %
        %   If mesh.normals is already populated (by WEC_Panelizer), those
        %   pre-computed values are used directly.
        %
        %   TWO-PANEL FIGURE
        %     Panel 1 (left)  — 3D perspective with quiver arrows
        %     Panel 2 (right) — Histogram of n_z components
        %
        %   WHY THE n_z HISTOGRAM?
        %     For a correctly oriented BEM hull:
        %       · Hull panels span n_z ∈ (−1, +1) continuously
        %       · A spike at n_z ≈ +1 with no negatives → normals are
        %         all pointing upward (wrong sign for submerged panels)
        %       · n_z = 0 exactly → degenerate (zero-area) panel
        %     The histogram gives an immediate global view without having
        %     to rotate the 3D figure.
        %
        %   INPUTS
        %     mesh    : struct from WEC_Panelizer.generate()
        %               .panels   [nP × 4]  vertex indices (col3==col4 for tris)
        %               .vertices [nV × 3]  (x,y,z) coordinates
        %               .normals  [nP × 3]  optional precomputed unit normals
        %
        %   OPTIONS (all optional)
        %     .scale        [m]    quiver arrow length; default = 0.12 × hull diameter
        %     .subsample    [-]    plot every N-th normal; default = 1 (all)
        %     .face_alpha   [-]    panel transparency;    default = 0.30
        %     .normal_color [1×3]  arrow RGB colour;      default = [0.90 0.20 0.10]
        %     .face_color   [1×3]  panel fill RGB;        default = [0.72 0.72 0.72]
        %     .view_az      [deg]  view azimuth;          default = −35
        %     .view_el      [deg]  view elevation;        default =  25

            try
                if nargin < 2, options = struct(); end
                if ~isfield(options, 'face_alpha'),   options.face_alpha   = 0.30;              end
                if ~isfield(options, 'face_color'),   options.face_color   = [0.72 0.72 0.72];  end
                if ~isfield(options, 'normal_color'), options.normal_color = [0.90 0.20 0.10];  end
                if ~isfield(options, 'subsample'),    options.subsample    = 1;                 end
                if ~isfield(options, 'view_az'),      options.view_az      = -35;               end
                if ~isfield(options, 'view_el'),      options.view_el      =  25;               end

                verts  = mesh.vertices;
                panels = mesh.panels;
                nP     = size(panels, 1);

                % ── Auto-scale: quiver length = 12% of hull diameter ──────
                hull_dia = max(verts(:,1)) - min(verts(:,1));
                if ~isfield(options, 'scale')
                    options.scale = max(0.05, 0.12 * hull_dia);
                end

                % ── Centroids and normals ──────────────────────────────────
                centroids = zeros(nP, 3);
                normals   = zeros(nP, 3);

                has_precomputed = isfield(mesh, 'normals') && ...
                                  ~isempty(mesh.normals)   && ...
                                  size(mesh.normals, 1) == nP;

                for p = 1:nP
                    v      = panels(p, :);
                    is_tri = (v(3) == v(4));
                    if is_tri
                        centroids(p, :) = mean(verts(v(1:3), :), 1);
                    else
                        centroids(p, :) = mean(verts(v, :), 1);
                    end

                    if ~has_precomputed
                        V1 = verts(v(1), :);  V2 = verts(v(2), :);
                        V3 = verts(v(3), :);  V4 = verts(v(4), :);
                        if is_tri
                            nrm = cross(V2 - V1, V3 - V2);
                        else
                            nrm = cross(V3 - V1, V4 - V2);   % HAMS diagonal convention
                        end
                        len = norm(nrm);
                        if len > 1e-12
                            normals(p, :) = nrm / len;
                        end
                        % Zero normal left as [0 0 0] → degenerate, visible in histogram
                    end
                end

                if has_precomputed
                    normals = mesh.normals;
                    lens    = sqrt(sum(normals.^2, 2));
                    valid_n = lens > 1e-12;
                    normals(valid_n, :) = normals(valid_n, :) ./ lens(valid_n);
                end

                % ── Subsample indices ──────────────────────────────────────
                idx_plot = 1 : options.subsample : nP;
                n_show   = length(idx_plot);

                % ── Figure layout ──────────────────────────────────────────
                figure('Name', 'Panel Normals Diagnostic', ...
                       'Color', 'white', 'Position', [80, 80, 1200, 600]);

                % ── Panel 1 — 3D hull + quiver arrows ─────────────────────
                ax1 = subplot(1, 2, 1);
                hold(ax1, 'on');

                for p = 1:nP
                    v = panels(p, :);
                    if v(3) == v(4), nv = 3; else, nv = 4; end
                    vi = v(1:nv);
                    fill3(ax1, verts(vi,1), verts(vi,2), verts(vi,3), ...
                          options.face_color, ...
                          'EdgeColor', [0.45 0.45 0.45], ...
                          'FaceAlpha', options.face_alpha, ...
                          'EdgeAlpha', 0.55);
                end

                quiver3(ax1, ...
                    centroids(idx_plot,1), centroids(idx_plot,2), centroids(idx_plot,3), ...
                    normals(idx_plot,1) * options.scale, ...
                    normals(idx_plot,2) * options.scale, ...
                    normals(idx_plot,3) * options.scale, ...
                    0, 'Color', options.normal_color, ...
                    'LineWidth', 1.1, 'MaxHeadSize', 0.55);

                view(ax1, options.view_az, options.view_el);
                axis(ax1, 'equal');
                xlabel(ax1, '$x$ [m]', 'Interpreter', 'latex', ...
                       'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                ylabel(ax1, '$y$ [m]', 'Interpreter', 'latex', ...
                       'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                zlabel(ax1, '$z$ [m]', 'Interpreter', 'latex', ...
                       'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                title(ax1, sprintf('Hull mesh + outward normals (%d panels, %d shown)', ...
                      nP, n_show), 'Interpreter', 'latex', ...
                      'FontSize', WEC_Visualization.FONT_SIZE_TITLE);
                WEC_Visualization.format_axis_publication(ax1);

                % ── Panel 2 — n_z histogram ────────────────────────────────
                ax2 = subplot(1, 2, 2);
                hold(ax2, 'on');

                nz_all = normals(:, 3);
                histogram(ax2, nz_all, 50, ...
                          'FaceColor', options.normal_color, ...
                          'EdgeColor', 'none', ...
                          'FaceAlpha', 0.75);
                xline(ax2, 0, 'k--', 'LineWidth', 1.5, ...
                      'Label', '$n_z=0$', 'Interpreter', 'latex');

                n_pos  = sum(nz_all >  0.01);
                n_neg  = sum(nz_all < -0.01);
                n_zero = sum(abs(nz_all) <= 0.01);

                if n_zero > 0
                    text(ax2, 0.02, 0.96, ...
                         sprintf('\\color{red}%d degenerate panels ($|n_z| \\leq 0.01$)', n_zero), ...
                         'Units', 'normalized', 'Interpreter', 'latex', ...
                         'FontSize', WEC_Visualization.FONT_SIZE_AXIS - 1, ...
                         'VerticalAlignment', 'top');
                end

                xlabel(ax2, '$n_z$ (normal $z$-component)', ...
                       'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                ylabel(ax2, 'Panel count', ...
                       'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                title(ax2, sprintf( ...
                    '$n_z>0$: %d,  $n_z<0$: %d,  $|n_z|\\leq0.01$: %d', ...
                    n_pos, n_neg, n_zero), ...
                    'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_TITLE);
                WEC_Visualization.format_axis_publication(ax2);

                sgtitle(sprintf('BEM Panel Normals  (%d panels, quiver scale = %.3f m)', ...
                                nP, options.scale), ...
                        'Interpreter', 'latex', ...
                        'FontName', WEC_Visualization.FONT_NAME, 'FontSize', 14);

                fname = sprintf('WEC_PanelNormals_%s.png', datestr(now, 'yyyymmdd_HHMMSS'));
                exportgraphics(gcf, fname, 'Resolution', 300);
                fprintf('  Panel normals figure saved: %s\n', fname);

            catch ME
                warning('WEC_Visualization:PanelNormalsFailed', ...
                        'Panel normals plot failed: %s', ME.message);
            end
        end

        %% =============================================================
        %%  BEM MESH DIAGNOSTIC (public — called from WEC_Driver)
        %% =============================================================
        % VIZ-ACCESS FIX: plot_mesh_diagnostic was in methods(Static,Access=private)
        % but WEC_Driver calls it as WEC_Visualization.plot_mesh_diagnostic(mesh,config)
        % from outside the class.  Moved here into the public methods(Static) block.

        function plot_mesh_diagnostic(mesh, config)
            try
                figure('Name', 'BEM Mesh Diagnostic', ...
                       'Color', 'white', 'Position', [80, 80, 1500, 550]);

                verts  = mesh.vertices;
                panels = mesh.panels;
                n_p    = size(panels, 1);

                ax1 = subplot(1, 3, 1);
                hold(ax1, 'on');
                for p = 1:n_p
                    v = panels(p, :);
                    if v(3) == v(4); nv = 3; else; nv = 4; end
                    vi = v(1:nv);
                    fill3(ax1, verts(vi, 1), verts(vi, 2), verts(vi, 3), ...
                          verts(vi, 3), 'EdgeColor', [0.3 0.3 0.3], ...
                          'FaceAlpha', 0.6, 'EdgeAlpha', 0.4);
                end
                z_tol_vis = 0.01;
                wl_idx = find(abs(verts(:, 3)) < z_tol_vis);
                if ~isempty(wl_idx)
                    plot3(ax1, verts(wl_idx, 1), verts(wl_idx, 2), ...
                          zeros(length(wl_idx), 1), ...
                          'Color', WEC_Visualization.WATERLINE_COLOR, ...
                          'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 8);
                end
                xl = [min(verts(:,1))-0.5, max(verts(:,1))+0.5];
                yl = [min(verts(:,2))-0.5, max(verts(:,2))+0.5];
                fill3(ax1, [xl(1) xl(2) xl(2) xl(1)], ...
                           [yl(1) yl(1) yl(2) yl(2)], ...
                           [0 0 0 0], ...
                           WEC_Visualization.WATERLINE_COLOR, ...
                           'FaceAlpha', 0.08, 'EdgeColor', 'none');
                xlabel('$x$ [m]', 'Interpreter', 'latex');
                ylabel('$y$ [m]', 'Interpreter', 'latex');
                zlabel('$z$ [m]', 'Interpreter', 'latex');
                title(sprintf('Hull mesh (%d panels)', n_p), ...
                      'Interpreter', 'latex', ...
                      'FontSize', WEC_Visualization.FONT_SIZE_TITLE);
                view(ax1, [-35, 25]);
                axis(ax1, 'equal');
                colormap(ax1, WEC_Visualization.cividis_map(64));
                cb = colorbar(ax1);
                ylabel(cb, '$z$ [m]', 'Interpreter', 'latex');
                WEC_Visualization.format_axis_publication(ax1);

                % ── Panel 2: Waterplane lid (top view) ─────────────
                % FIX (L8): pre-initialise WP vars so Panel 3 is safe if Panel 2 throws.
                n_wp = 0;  wp_n = zeros(0,3);  wp_p = zeros(0,4);  wp_nv = zeros(0,1);

                ax2 = subplot(1, 3, 2);
                hold(ax2, 'on');

                wp_target = 0.4;
                if isfield(config, 'wp_target_edge')
                    wp_target = config.wp_target_edge;
                end

                % VIZ-MATCH FIX: replace extract_waterline_boundary +
                % generate_waterplane_mesh_structured (TFI, boundary resampled
                % at target_edge spacing — decoupled from hull nodes) with the
                % same pipeline used in run_single_hams:
                %   hull_waterline_polygon  — open-edge walk on mesh.panels,
                %                             returns exact hull mesh nodes
                %   mesh_wp_blossomquad(..., lock_boundary=true)
                %                           — CDT with hull nodes locked on
                %                             the outer boundary ring
                % This ensures the diagnostic shows exactly what HAMS received.
                boundary_xy = HAMS_Pipeline.hull_waterline_polygon(mesh);

                if ~isempty(boundary_xy) && size(boundary_xy,1) >= 3
                    [wp_n, wp_p, wp_nv] = HAMS_Pipeline.mesh_wp_blossomquad( ...
                        boundary_xy, wp_target, mesh.x_sym, mesh.y_sym, 0, true);
                    if ~isempty(wp_n), wp_n(:,3) = 0; end

                    n_wp = size(wp_p, 1);
                    for p = 1:n_wp
                        nv = wp_nv(p);
                        vi = wp_p(p, 1:nv);
                        if nv == 4
                            fc = [0.65, 0.82, 1.0];
                            ec = [0.1, 0.3, 0.7];
                        else
                            fc = [0.65, 0.95, 0.75];
                            ec = [0.05, 0.45, 0.2];
                        end
                        fill(ax2, wp_n(vi, 1), wp_n(vi, 2), fc, ...
                             'EdgeColor', ec, 'FaceAlpha', 0.6);
                    end

                    plot(ax2, boundary_xy(:, 1), boundary_xy(:, 2), ...
                         'k-', 'LineWidth', WEC_Visualization.LW_BOUNDARY);
                    plot(ax2, boundary_xy(1, 1), boundary_xy(1, 2), ...
                         'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');

                    max_edge_wp = 0;
                    for p = 1:n_wp
                        nv = wp_nv(p); vi = wp_p(p, 1:nv);
                        pts = wp_n(vi, :);
                        for e = 1:nv
                            en = mod(e, nv) + 1;
                            max_edge_wp = max(max_edge_wp, ...
                                norm(pts(en,:) - pts(e,:)));
                        end
                    end
                    title(sprintf('WP lid (%d quad + %d tri, max edge %.2f m)', ...
                          sum(wp_nv == 4), sum(wp_nv == 3), max_edge_wp), ...
                          'Interpreter', 'latex', ...
                          'FontSize', WEC_Visualization.FONT_SIZE_TITLE);
                else
                    title('WP lid (no waterline found)', ...
                          'Interpreter', 'latex');
                end

                xlabel('$x$ [m]', 'Interpreter', 'latex');
                ylabel('$y$ [m]', 'Interpreter', 'latex');
                axis(ax2, 'equal');
                WEC_Visualization.format_axis_publication(ax2);

                ax3 = subplot(1, 3, 3);
                hold(ax3, 'on');
                for p = 1:n_p
                    v = panels(p, :);
                    if v(3) == v(4); nv = 3; else; nv = 4; end
                    vi = v(1:nv);
                    fill3(ax3, verts(vi, 1), verts(vi, 2), verts(vi, 3), ...
                          [0.7 0.7 0.7], 'EdgeColor', [0.5 0.5 0.5], ...
                          'FaceAlpha', 0.25, 'EdgeAlpha', 0.3);
                end
                if ~isempty(boundary_xy)
                    for p = 1:n_wp
                        nv = wp_nv(p); vi = wp_p(p, 1:nv);
                        fill3(ax3, wp_n(vi, 1), wp_n(vi, 2), wp_n(vi, 3), ...
                              [0.3 0.5 0.9], 'EdgeColor', [0.1 0.2 0.5], ...
                              'FaceAlpha', 0.7, 'EdgeAlpha', 0.5);
                    end
                end
                if ~isempty(wl_idx)
                    plot3(ax3, verts(wl_idx, 1), verts(wl_idx, 2), ...
                          zeros(length(wl_idx), 1), '.', ...
                          'Color', WEC_Visualization.WATERLINE_COLOR, ...
                          'MarkerSize', 6);
                end
                xlabel('$x$ [m]', 'Interpreter', 'latex');
                ylabel('$y$ [m]', 'Interpreter', 'latex');
                zlabel('$z$ [m]', 'Interpreter', 'latex');
                title('Hull + WP lid (combined)', ...
                      'Interpreter', 'latex', ...
                      'FontSize', WEC_Visualization.FONT_SIZE_TITLE);
                view(ax3, [-40, 20]);
                axis(ax3, 'equal');
                WEC_Visualization.format_axis_publication(ax3);

                sgtitle('BEM Mesh Diagnostic', ...
                        'FontSize', 14, 'FontWeight', 'bold', ...
                        'FontName', WEC_Visualization.FONT_NAME);

            catch ME
                warning('WEC_Visualization:MeshDiag', ...
                        'Mesh diagnostic plot failed: %s', ME.message);
            end
        end

    end % methods (Static)

    methods (Static, Access = private)
        %
        function ring = mesh_slice_ring(faces, verts, z_cut)
            tol = 1e-9;
            n_faces = size(faces, 1);
            seg_pts = zeros(n_faces*2, 3);
            n_segs  = 0;
            edges   = [1 2; 2 3; 3 1];
            for f = 1:n_faces
                v  = verts(faces(f,:), :);
                zv = v(:, 3);
                if ~any(zv > z_cut+tol) || ~any(zv < z_cut-tol), continue; end
                cross_pts = zeros(0,3);
                for e = 1:3
                    z1 = zv(edges(e,1));  z2 = zv(edges(e,2));
                    if (z1-z_cut)*(z2-z_cut) < 0
                        t = (z_cut-z1)/(z2-z1);
                        cross_pts(end+1,:) = v(edges(e,1),:) + t*(v(edges(e,2),:)-v(edges(e,1),:)); %#ok<AGROW>
                    end
                end
                if size(cross_pts,1) == 2
                    n_segs = n_segs+1;
                    seg_pts(2*n_segs-1,:) = cross_pts(1,:);
                    seg_pts(2*n_segs,  :) = cross_pts(2,:);
                end
            end
            if n_segs == 0, ring = []; return; end
            ring = WEC_Visualization.chain_segs(seg_pts(1:2*n_segs,:));
        end

        function chain = chain_segs(seg_pts)
            n_segs = size(seg_pts,1)/2;
            used   = false(n_segs,1);
            chain  = [seg_pts(1,:); seg_pts(2,:)];
            used(1) = true;
            cur    = seg_pts(2,:);
            for iter = 1:n_segs-1
                best_d = inf; bs = 0; flip_seg = false;
                for s = 1:n_segs
                    if used(s), continue; end
                    d1 = norm(cur - seg_pts(2*s-1,:));
                    d2 = norm(cur - seg_pts(2*s,  :));
                    if d1 < best_d, best_d = d1; bs = s; flip_seg = false; end
                    if d2 < best_d, best_d = d2; bs = s; flip_seg = true;  end
                end
                if best_d > 0.5 || bs == 0, break; end
                used(bs) = true;
                if flip_seg, npt = seg_pts(2*bs-1,:); else, npt = seg_pts(2*bs,:); end
                chain(end+1,:) = npt; %#ok<AGROW>
                cur = npt;
            end
            if norm(chain(1,:)-chain(end,:)) < 0.5, chain(end+1,:) = chain(1,:); end
        end
        %
        function draw_density_strips(ax, profile, node_z, densities, cmap, clim_range, strip_bounds)
            N      = length(densities);
            node_z = node_z(:);

            % Slab boundaries: use explicit bounds if provided, else midpoints.
            % FIX (L7): guard N=1 — node_z(2) would be out-of-bounds.
            if nargin >= 7 && ~isempty(strip_bounds) && length(strip_bounds) == N + 1
                bounds = strip_bounds(:);
            elseif N > 1
                dz_half  = (node_z(2) - node_z(1)) / 2;
                bounds   = [ node_z(1) - dz_half; ...
                            (node_z(1:end-1) + node_z(2:end)) / 2; ...
                             node_z(end) + dz_half ];
            else
                % N == 1: single strip spanning the full profile z-range
                z_rng  = max(abs(node_z), 0.5);   % at least ±0.5 m half-height
                bounds = [ node_z(1) - z_rng; node_z(1) + z_rng ];
            end

            for k = 1:N
                z_lo      = bounds(k);
                z_hi      = bounds(k+1);
                rho_strip = densities(k);

                slab = WEC_Visualization.clip_profile_to_z_range(profile, z_lo, z_hi);
                if isempty(slab) || size(slab,1) < 3
                    continue;
                end

                patch(ax, slab(:,1), slab(:,2), rho_strip, ...
                      'EdgeColor', 'none', 'FaceAlpha', 0.85);
            end
            %
            for k = 1:length(bounds)
                x_cross = [];
                for j = 1:size(profile,1)
                    j2 = mod(j, size(profile,1)) + 1;
                    z1 = profile(j,2);  z2 = profile(j2,2);
                    if (z1 - bounds(k)) * (z2 - bounds(k)) <= 0 && abs(z2-z1) > 1e-12
                        t = (bounds(k) - z1) / (z2 - z1);
                        x_cross(end+1) = profile(j,1) + t*(profile(j2,1) - profile(j,1)); %#ok<AGROW>
                    end
                end
                if length(x_cross) >= 2
                    plot(ax, [min(x_cross), max(x_cross)], [bounds(k), bounds(k)], ...
                         '-', 'Color', [0.1, 0.1, 0.1], 'LineWidth', 0.8);
                end
            end
        end

        function clipped = clip_profile_to_z_range(profile, z_lo, z_hi)
            poly = profile;
            poly = WEC_Visualization.sh_clip_halfplane(poly, z_lo, +1);
            poly = WEC_Visualization.sh_clip_halfplane(poly, z_hi, -1);
            if size(poly, 1) < 3
                clipped = [];
            else
                clipped = poly;
            end
        end

        function out = sh_clip_halfplane(poly, z_bound, sign_dir)
            if isempty(poly) || size(poly, 1) < 2
                out = [];
                return;
            end
            n   = size(poly, 1);
            out = zeros(2*n, 2);
            cnt = 0;
            for i = 1:n
                j  = mod(i, n) + 1;
                zi = poly(i, 2);
                zj = poly(j, 2);
                inside_i = sign_dir * (zi - z_bound) >= -1e-12;
                inside_j = sign_dir * (zj - z_bound) >= -1e-12;
                if inside_i && inside_j
                    cnt = cnt + 1;
                    out(cnt, :) = poly(j, :);
                elseif inside_i && ~inside_j
                    t = (z_bound - zi) / (zj - zi);
                    cnt = cnt + 1;
                    out(cnt, :) = poly(i, :) + t * (poly(j, :) - poly(i, :));
                elseif ~inside_i && inside_j
                    t = (z_bound - zi) / (zj - zi);
                    cnt = cnt + 1;
                    out(cnt, :) = poly(i, :) + t * (poly(j, :) - poly(i, :));
                    cnt = cnt + 1;
                    out(cnt, :) = poly(j, :);
                end
            end
            out = out(1:cnt, :);
        end

        function inner_prof = compute_inner_profile_2d(outer_profile, t_shell)
            x = outer_profile(:,1);
            z = outer_profile(:,2);
            Nv = length(x);
            if abs(x(end)-x(1)) < 1e-12 && abs(z(end)-z(1)) < 1e-12
                x = x(1:end-1); z = z(1:end-1); Nv = Nv - 1;
            end
            if Nv < 3, inner_prof = outer_profile; return; end
            signed_area = 0.5 * sum(x .* circshift(z,-1) - circshift(x,-1) .* z);
            if signed_area < 0
                x = flipud(x); z = flipud(z);
            end
            x_off = zeros(Nv, 1);
            z_off = zeros(Nv, 1);
            for j = 1:Nv
                jm = mod(j-2, Nv) + 1;
                jp = mod(j,   Nv) + 1;
                e_prev = [x(j)-x(jm), z(j)-z(jm)];
                e_next = [x(jp)-x(j), z(jp)-z(j)];
                lp = norm(e_prev); ln = norm(e_next);
                if lp < 1e-12 || ln < 1e-12
                    x_off(j) = x(j); z_off(j) = z(j); continue;
                end
                n_prev = [-e_prev(2),  e_prev(1)] / lp;
                n_next = [-e_next(2),  e_next(1)] / ln;
                n_avg    = n_prev + n_next;
                len_avg  = norm(n_avg);
                if len_avg < 1e-12, n_avg = n_prev; len_avg = 1; end
                n_avg    = n_avg / len_avg;
                cos_half = dot(n_avg, n_prev);
                if cos_half > 0.33
                    miter = t_shell / cos_half;
                else
                    miter = t_shell * 3.0;
                end
                x_off(j) = x(j) + n_avg(1) * miter;
                z_off(j) = z(j) + n_avg(2) * miter;
            end
            inner_prof = [x_off, z_off];
        end

        function plot_inner_spline(ax, inner_prof)
            if isempty(inner_prof) || size(inner_prof,1) < 3, return; end
            x = inner_prof(:,1);
            z = inner_prof(:,2);
            x_cl = [x; x(1)];
            z_cl = [z; z(1)];
            segs = sqrt(diff(x_cl).^2 + diff(z_cl).^2);
            segs = max(segs, 1e-12);
            t_param = [0; cumsum(segs)];
            t_fine  = linspace(0, t_param(end), 400);
            x_sp = interp1(t_param, x_cl, t_fine, 'spline');
            z_sp = interp1(t_param, z_cl, t_fine, 'spline');
            plot(ax, x_sp, z_sp, '--', ...
                 'Color', [0.15, 0.15, 0.15], ...
                 'LineWidth', 1.2, ...
                 'DisplayName', 'Shell boundary');
        end

        function draw_composite_strips(ax, outer_profile, inner_profile, ...
                                       node_z, densities, clim_range)
            SHELL_GRAY    = [0.82, 0.82, 0.82];
            BOUND_COLOR   = [0.10, 0.10, 0.10];
            N      = length(densities);
            node_z = node_z(:);
            dz_half = (node_z(2) - node_z(1)) / 2;
            bounds  = [ node_z(1) - dz_half; ...
                       (node_z(1:end-1) + node_z(2:end)) / 2; ...
                        node_z(end) + dz_half ];
            outer_slabs = cell(N, 1);
            inner_slabs = cell(N, 1);
            for k = 1:N
                outer_slabs{k} = WEC_Visualization.clip_profile_to_z_range( ...
                    outer_profile, bounds(k), bounds(k+1));
                inner_slabs{k} = WEC_Visualization.clip_profile_to_z_range( ...
                    inner_profile, bounds(k), bounds(k+1));
            end
            for k = 1:N
                outer_slab = outer_slabs{k};
                if isempty(outer_slab) || size(outer_slab,1) < 3, continue; end
                patch(ax, outer_slab(:,1), outer_slab(:,2), SHELL_GRAY, ...
                      'EdgeColor', 'none', 'FaceAlpha', 1.0);
            end
            for k = 1:N
                inner_slab = inner_slabs{k};
                if isempty(inner_slab) || size(inner_slab,1) < 3, continue; end
                patch(ax, inner_slab(:,1), inner_slab(:,2), densities(k), ...
                      'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 1.0);
            end
            for k = 1:N
                z_lo = bounds(k);
                z_hi = bounds(k+1);
                outer_slab = outer_slabs{k};
                if isempty(outer_slab) || size(outer_slab,1) < 3, continue; end
                for bz = [z_lo, z_hi]
                    x_cross = [];
                    n_out = size(outer_profile,1);
                    for j = 1:n_out
                        j2 = mod(j, n_out) + 1;
                        z1 = outer_profile(j,2);  z2 = outer_profile(j2,2);
                        if (z1 - bz) * (z2 - bz) <= 0 && abs(z2-z1) > 1e-12
                            t_edge = (bz - z1) / (z2 - z1);
                            x_cross(end+1) = outer_profile(j,1) + ...
                                t_edge*(outer_profile(j2,1) - outer_profile(j,1)); %#ok<AGROW>
                        end
                    end
                    if length(x_cross) >= 2
                        line(ax, [min(x_cross), max(x_cross)], [bz, bz], ...
                             'Color', BOUND_COLOR, 'LineWidth', 0.8);
                    end
                end
            end
        end

        function rho_eq = compute_strip_equivalent_density(config, densities_core)
            N       = length(densities_core);
            rho_eq  = densities_core(:);
            if isempty(config.shell), return; end
            V_shell = config.shell.V_shell;
            V_core  = config.shell.V_core;
            V_total = V_shell + V_core;
            for i = 1:N
                if V_total(i) > 1e-12
                    rho_eq(i) = (config.shell_density * V_shell(i) + ...
                                 densities_core(i)    * V_core(i)) / V_total(i);
                end
            end
        end

        function format_axis_publication(ax)
            grid(ax, 'on'); box(ax, 'on');
            set(ax, 'GridLineStyle', ':', 'GridAlpha', 0.4, ...
                    'LineWidth', WEC_Visualization.LW_AXES, ...
                    'FontName', WEC_Visualization.FONT_NAME, ...
                    'FontSize', WEC_Visualization.FONT_SIZE_AXIS, ...
                    'TickLabelInterpreter', 'latex');
        end

        function shade_T_band(ax, T_band)
            yl = ylim(ax);
            fill(ax, [T_band(1) T_band(2) T_band(2) T_band(1)], ...
                 [yl(1) yl(1) yl(2) yl(2)], ...
                 [0.85 0.92 1.0], 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
                 'HandleVisibility', 'off');
        end

        function cmap = cividis_map(n)
            ctrl = [0.000, 0.135, 0.305; 0.000, 0.205, 0.380; 0.123, 0.263, 0.406;
                    0.253, 0.318, 0.420; 0.365, 0.373, 0.432; 0.479, 0.435, 0.440;
                    0.608, 0.518, 0.430; 0.770, 0.640, 0.380; 0.995, 0.906, 0.144];
            cmap = flipud(interp1(linspace(0, 1, 9), ctrl, linspace(0, 1, n)));
        end

        function plot_stage2_panel(results, fields, labels, panel_title, ax)
            has_data = isfield(results, 'stage2_3d') && isfield(results.stage2_3d, 'iteration_errors');
            hold(ax, 'on');
            if has_data
                err = results.stage2_3d.iteration_errors;
                markers = {'b-o', 'r-s'};
                for k = 1:length(fields)
                    if isfield(err, fields{k}) && ~isempty(err.(fields{k}))
                        vec = err.(fields{k});
                        iters_3d = 0:length(vec)-1;
                        plot(ax, iters_3d, vec, markers{k}, 'LineWidth', WEC_Visualization.LW_MAIN, 'MarkerSize', 6);
                    end
                end
                ylabel(ax, 'Normalised Error', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
                legend(ax, labels, 'Location', 'best', 'Interpreter', 'latex', 'FontName', WEC_Visualization.FONT_NAME);
            else
                text(ax, 0.5, 0.5, 'Stage 2 data not available', ...
                     'Interpreter', 'latex', 'Units', 'normalized', 'HorizontalAlignment', 'center');
            end
            xlabel(ax, '3D Optimizer Iteration', 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_LABEL);
            title(ax, panel_title, 'Interpreter', 'latex', 'FontSize', WEC_Visualization.FONT_SIZE_TITLE);
            WEC_Visualization.format_axis_publication(ax);
            hold(ax, 'off');
        end


        function boundary_xy = extract_waterline_boundary(mesh)
            boundary_xy = [];
            z_tol = 0.01;
            ec = containers.Map('KeyType','char','ValueType','int32');
            em = containers.Map('KeyType','char','ValueType','any');
            for p = 1:size(mesh.panels, 1)
                v = mesh.panels(p,:);
                if v(3)==v(4); ee=[v(1) v(2);v(2) v(3);v(3) v(1)];
                else; ee=[v(1) v(2);v(2) v(3);v(3) v(4);v(4) v(1)]; end
                for e = 1:size(ee,1)
                    key = sprintf('%d_%d', min(ee(e,:)), max(ee(e,:)));
                    if ec.isKey(key); ec(key)=ec(key)+1;
                    else; ec(key)=1; em(key)=ee(e,:); end
                end
            end
            wle = zeros(0,2); ka = ec.keys();
            for i = 1:length(ka)
                if ec(ka{i})==1
                    ev=em(ka{i});
                    if abs(mesh.vertices(ev(1),3))<z_tol && abs(mesh.vertices(ev(2),3))<z_tol
                        wle(end+1,:) = ev; %#ok<AGROW>
                    end
                end
            end
            if isempty(wle); return; end
            wvi = unique(wle(:));
            wxy = mesh.vertices(wvi,:); nw = size(wxy,1);
            g2l = containers.Map('KeyType','int32','ValueType','int32');
            for i=1:nw; g2l(wvi(i))=i; end
            le = zeros(size(wle));
            for e=1:size(wle,1); le(e,:)=[g2l(wle(e,1)),g2l(wle(e,2))]; end
            cn=(1:nw)';
            for ii=2:nw; for jj=1:ii-1
                if cn(jj)~=jj; continue; end
                if norm(wxy(ii,:)-wxy(jj,:))<1e-6; cn(ii)=jj; break; end
            end; end
            for e=1:size(le,1); le(e,:)=[cn(le(e,1)),cn(le(e,2))]; end
            le(le(:,1)==le(:,2),:)=[];
            [~,ui]=unique(sort(le,2),'rows'); le=le(ui,:);
            uid=unique(cn); nu=length(uid);
            ni=zeros(nw,1); ni(uid)=(1:nu)';
            ce=zeros(size(le));
            for e=1:size(le,1); ce(e,:)=[ni(cn(le(e,1))),ni(cn(le(e,2)))]; end
            uxy=wxy(uid,1:2);
            adj=cell(nu,1);
            for e=1:size(ce,1); adj{ce(e,1)}(end+1)=ce(e,2); adj{ce(e,2)}(end+1)=ce(e,1); end
            vis=false(nu,1); ord=zeros(nu,1); ord(1)=1; vis(1)=true;
            for step=2:nu; c=ord(step-1); nx=0;
                for ni_=1:length(adj{c}); if ~vis(adj{c}(ni_)); nx=adj{c}(ni_); break; end; end
                if nx==0; break; end; ord(step)=nx; vis(nx)=true;
            end
            nc=find(ord>0,1,'last'); ord=ord(1:nc);
            boundary_xy = uxy(ord,:);
        end

    end  % methods (Static, Access = private)

    methods (Static)

        function profile = build_smooth_viz_profile(config, n_levels)
        % BUILD_SMOOTH_VIZ_PROFILE  Parametric hull silhouette for 2D visualization.
        %
        %   profile = WEC_Visualization.build_smooth_viz_profile(config)
        %   profile = WEC_Visualization.build_smooth_viz_profile(config, n_levels)
        %
        %   Returns a smooth [2*n_levels × 2] closed polygon [x, z] representing
        %   the OUTER SILHOUETTE of the hull in the midplane (y=0 view).
        %
        %   CONSTRUCTION
        %     1. Sample n_levels z-values spanning [hull_z_min, hull_z_max].
        %     2. At each z, call extract_isocurve_at_z → waterplane contour points.
        %     3. max(pts(:,1)) gives the right half of the silhouette at that z.
        %     4. Build closed CCW polygon: right side (bottom→top) + mirrored
        %        left side (top→bottom).
        %
        %   WHY NOT config.profile?
        %     extractProfileMS2 collects boundary edge points from all surfaces and
        %     sorts them by ANGLE FROM CENTROID.  This is correct for convex polygons
        %     only.  C0 is NON-CONVEX (wide platform → narrow column): the angular
        %     sort produces a self-intersecting 'star-shaped' polygon at the
        %     platform/column junction → the spiky visual artefact seen in the plot.
        %
        %   COST: n_levels calls to extract_isocurve_at_z.  With boundary_cache
        %     (pre-built in WEC_Configuration_Builder), each call is <1 ms
        %     (pure arithmetic on cached arrays) → ~200 ms total for 300 levels.
        %
        %   FALLBACK: returns config.profile if ms2_model or boundary_cache missing.

            if nargin < 2, n_levels = 300; end

            has_model = isfield(config, 'ms2_model')      && ~isempty(config.ms2_model);
            has_cache = isfield(config, 'boundary_cache') && ~isempty(config.boundary_cache);
            has_zlims = isfield(config, 'hull_z_min')     && isfield(config, 'hull_z_max');

            if ~has_model || ~has_cache || ~has_zlims
                % Fallback: config.profile may be coarse/incorrect for non-convex hulls.
                profile = config.profile;
                return;
            end

            z_lo   = config.hull_z_min;
            z_hi   = config.hull_z_max;
            z_vals = linspace(z_lo + 1e-4, z_hi - 1e-4, n_levels)';

            n_u_cache = length(config.boundary_cache.u_samples);
            x_right   = zeros(n_levels, 1);

            for k = 1:n_levels
                pts = WEC_HydroProperties.extract_isocurve_at_z( ...
                    config.ms2_model, z_vals(k), n_u_cache, config.boundary_cache);
                if ~isempty(pts) && size(pts, 1) >= 2
                    x_right(k) = max(pts(:, 1));
                end
            end

            % Remove z-levels above/below hull (zero contour)
            valid = x_right > 1e-6;
            x_r   = x_right(valid);
            z_r   = z_vals(valid);

            if length(x_r) < 3
                profile = config.profile;   % fallback
                return;
            end

            % Build closed CCW polygon (right side bottom→top, left side top→bottom).
            % The polygon naturally closes at the keel (x_r ≈ 0 at both ends).
            profile = [x_r, z_r; -flipud(x_r), flipud(z_r)];
        end

    end
end