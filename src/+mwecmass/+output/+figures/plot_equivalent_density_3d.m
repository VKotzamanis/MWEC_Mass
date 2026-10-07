function plot_equivalent_density_3d(props, config)
%PLOT_EQUIVALENT_DENSITY_3D Draw 3D hull colour-filled by strip-average equivalent density.
% Inputs: props (realised_strip_density if available, fill_method); config (ms2_model, viz_Nu, density_nodes_z, strip_edges); vertical shift SI [m].
% Fall back to continuous optimiser densities if props.realised_strip_density is unavailable.

    opts  = mwecmass.output.figures.output_options(config);
    style = opts.style;
    try
        % Density source: realised vs optimiser
        if isfield(props, 'realised_strip_density') && ...
                ~isempty(props.realised_strip_density)
            rho_eq = props.realised_strip_density(:);
            rho_source_is_realised = true;
        else
            rho_eq = [props.components.density]';
            rho_source_is_realised = false;
        end

        if rho_source_is_realised
            if isfield(props, 'fill_method') && ...
                    strcmp(props.fill_method, 'uhpc_fill')
                fig_name = 'WEC Equivalent Density (3D) — Realised UHPC';
            else
                fig_name = 'WEC Equivalent Density (3D) — Realised Steel';
            end
        else
            fig_name = 'WEC Equivalent Density (3D)';
        end
        fig = mwecmass.output.figures.new_figure(style, 'tall_single');
        fig.Name = fig_name;
        ax = axes(fig);
        hold(ax, 'on');

        VIZ_N_eq = 60;
        if isfield(config, 'viz_Nu'), VIZ_N_eq = config.viz_Nu; end
        opts_viz = struct('trim_wl', false, 'close_gaps', true, ...
                          'verbose', false, 'quarter_body', false);

        mesh_viz = mwecmass.mesh.generate( ...
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

        colormap(ax, mwecmass.output.figures.figure_colormap(style, 256));
        cbar = colorbar(ax);
        % $...$-wrapped: style_colorbar now styles cbar.Label through style_text, which forces
        % the LaTeX interpreter (previously unset, defaulting to MATLAB's 'tex', which does not
        % require the math delimiters LaTeX does for \rho/_/^).
        cbar.Label.String = 'Density, $\rho_{eq}$ (kg/m$^3$)';
        mwecmass.output.figures.style_colorbar(cbar, style);
        d_min = min(rho_eq); d_max = max(rho_eq);
        if d_max > d_min, clim(ax, [d_min, d_max]);
        else,             clim(ax, [d_min-1, d_max+1]); end

        % Waterplane
        xlims = [min(verts(:,1)), max(verts(:,1))];
        ylims = [min(verts(:,2)), max(verts(:,2))];
        [xg, yg] = meshgrid(linspace(xlims(1), xlims(2), 10), ...
                            linspace(ylims(1), ylims(2), 10));
        h_water = surf(ax, xg, yg, zeros(size(xg)), ...
            'FaceColor', style.fill_palette.waterline, ...
            'FaceAlpha', 0.4, 'EdgeColor', 'none');
        h_cg = plot3(ax, props.CG_total(1), props.CG_total(2), props.CG_total(3), ...
            'o', 'Color', style.fill_palette.boundary, 'MarkerFaceColor', style.color.cg);
        h_cb = plot3(ax, props.CB(1), props.CB(2), props.CB(3), ...
            's', 'Color', style.fill_palette.boundary, 'MarkerFaceColor', style.color.cb);

        axis(ax, 'equal'); view(ax, 30, 25);
        xlabel(ax, '$x$ [m]');
        ylabel(ax, '$y$ [m]');
        zlabel(ax, '$z$ [m]');
        if rho_source_is_realised
            if isfield(props, 'fill_method') && ...
                    strcmp(props.fill_method, 'uhpc_fill')
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
        title(ax, title_str);
        lg = legend(ax, [h_water, h_cg, h_cb], {'Waterplane', 'CG', 'CB'}, ...
               'Location', 'bestoutside');
        mwecmass.output.figures.style_legend(lg, style);
        mwecmass.output.figures.apply_axes_style(ax, style);
        hold(ax, 'off');

        % The as-built and the optimiser-design figures are the same drawing of different
        % densities, so each keeps its own file name and its own flag.
        if rho_source_is_realised
            name_stem = 'WEC_Realised_Density_3D';
            write_figure = opts.save.stage3.realised_density_3d;
        else
            name_stem = 'WEC_Equivalent_Density_3D';
            write_figure = opts.save.stage2.density_3d;
        end
        if write_figure
            saved = mwecmass.output.figures.export_figure(fig, name_stem, config);
            fprintf('  Equivalent-density (3D) figure saved: %s\n', saved{1});
        end

    catch ME
        warning('mwecmass:figures:EquivalentDensity3DFailed', ...
                '3D equivalent visualization failed: %s', ME.message);
    end
end
