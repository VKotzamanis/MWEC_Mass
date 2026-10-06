function plot_mesh_diagnostic(mesh, config)
%PLOT_MESH_DIAGNOSTIC Plot BEM mesh, waterplane lid (CDT), and combined hull + WP mesh for verification.
% Inputs: mesh (vertices, panels, x_sym, y_sym); config (wp_target_edge default 0.4 m). Panel faces colour by z-coordinate (waterline at z=0).

    opts  = mwecmass.output.figures.output_options(config);
    style = opts.style;
    try
        fig = mwecmass.output.figures.new_figure(style, 'double_column');
        set(fig, 'Name', 'BEM Mesh Diagnostic');
        layout = tiledlayout(fig, 1, 3);

        verts  = mesh.vertices;
        panels = mesh.panels;
        n_p    = size(panels, 1);

        ax1 = nexttile(layout, 1);
        hold(ax1, 'on');
        for p = 1:n_p
            v = panels(p, :);
            if v(3) == v(4); nv = 3; else; nv = 4; end
            vi = v(1:nv);
            fill3(ax1, verts(vi, 1), verts(vi, 2), verts(vi, 3), ...
                  verts(vi, 3), 'EdgeColor', style.color.mesh_edge, ...
                  'FaceAlpha', 0.6, 'EdgeAlpha', 0.4);
        end
        z_tol_vis = 0.01;
        wl_idx = find(abs(verts(:, 3)) < z_tol_vis);
        if ~isempty(wl_idx)
            plot3(ax1, verts(wl_idx, 1), verts(wl_idx, 2), ...
                  zeros(length(wl_idx), 1), ...
                  'Color', style.fill_palette.waterline, ...
                  'LineStyle', 'none', 'Marker', '.');
        end
        xl = [min(verts(:,1))-0.5, max(verts(:,1))+0.5];
        yl = [min(verts(:,2))-0.5, max(verts(:,2))+0.5];
        fill3(ax1, [xl(1) xl(2) xl(2) xl(1)], ...
                   [yl(1) yl(1) yl(2) yl(2)], ...
                   [0 0 0 0], ...
                   style.fill_palette.waterline, ...
                   'FaceAlpha', 0.08, 'EdgeColor', 'none');
        xlabel(ax1, '$x$ [m]');
        ylabel(ax1, '$y$ [m]');
        zlabel(ax1, '$z$ [m]');
        title(ax1, sprintf('Hull mesh (%d panels)', n_p));
        view(ax1, [-35, 25]);
        axis(ax1, 'equal');
        colormap(ax1, mwecmass.output.figures.figure_colormap(style, 64));
        cb = colorbar(ax1);
        % style_colorbar now styles cb.Label too, so the explicit LaTeX
        % math-interpreter argument here is redundant with it (both resolve to
        % style.tick_label_interpreter) and is dropped.
        ylabel(cb, '$z$ [m]');
        mwecmass.output.figures.style_colorbar(cb, style);
        mwecmass.output.figures.apply_axes_style(ax1, style);

        % ── Panel 2: Waterplane lid (top view) ─────────────
        % Pre-initialize WP vars so Panel 3 is safe if Panel 2 throws.
        n_wp = 0;  wp_n = zeros(0,3);  wp_p = zeros(0,4);  wp_nv = zeros(0,1);

        ax2 = nexttile(layout, 2);
        hold(ax2, 'on');

        wp_target = 0.4;
        if isfield(config, 'wp_target_edge')
            wp_target = config.wp_target_edge;
        end

        % VIZ-MATCH FIX: replace extract_waterline_boundary +
        % generate_waterplane_mesh_structured (TFI, boundary resampled
        % at target_edge spacing — decoupled from hull nodes) with the
        % same pipeline used in run_at_draft:
        %   hull_waterline_polygon  — open-edge walk on mesh.panels,
        %                             returns exact hull mesh nodes
        %   mesh_wp_blossomquad(..., lock_boundary=true)
        %                           — CDT with hull nodes locked on
        %                             the outer boundary ring
        % This ensures the diagnostic shows exactly what HAMS received.
        boundary_xy = mwecmass.mesh.hull_waterline_polygon(mesh);

        if ~isempty(boundary_xy) && size(boundary_xy,1) >= 3
            [wp_n, wp_p, wp_nv] = mwecmass.mesh.waterplane_mesh_unstructured( ...
                boundary_xy, wp_target, mesh.x_sym, mesh.y_sym, 0, true);
            if ~isempty(wp_n), wp_n(:,3) = 0; end

            n_wp = size(wp_p, 1);
            for p = 1:n_wp
                nv = wp_nv(p);
                vi = wp_p(p, 1:nv);
                if nv == 4
                    fc = style.color.mesh_quad.fill;
                    ec = style.color.mesh_quad.edge;
                else
                    fc = style.color.mesh_tri.fill;
                    ec = style.color.mesh_tri.edge;
                end
                fill(ax2, wp_n(vi, 1), wp_n(vi, 2), fc, ...
                     'EdgeColor', ec, 'FaceAlpha', 0.6);
            end

            hb = plot(ax2, boundary_xy(:, 1), boundary_xy(:, 2), ...
                 '-', 'Color', style.fill_palette.boundary);
            mwecmass.output.figures.style_line(hb, style, 'boundary');
            plot(ax2, boundary_xy(1, 1), boundary_xy(1, 2), ...
                 'o', 'Color', style.fill_palette.boundary, ...
                 'MarkerFaceColor', style.fill_palette.boundary);

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
            title(ax2, sprintf('WP lid (%d quad + %d tri, max edge %.2f m)', ...
                  sum(wp_nv == 4), sum(wp_nv == 3), max_edge_wp));
        else
            title(ax2, 'WP lid (no waterline found)');
        end

        xlabel(ax2, '$x$ [m]');
        ylabel(ax2, '$y$ [m]');
        axis(ax2, 'equal');
        mwecmass.output.figures.apply_axes_style(ax2, style);

        ax3 = nexttile(layout, 3);
        hold(ax3, 'on');
        for p = 1:n_p
            v = panels(p, :);
            if v(3) == v(4); nv = 3; else; nv = 4; end
            vi = v(1:nv);
            fill3(ax3, verts(vi, 1), verts(vi, 2), verts(vi, 3), ...
                  style.fill_palette.shell, 'EdgeColor', style.color.mesh_edge, ...
                  'FaceAlpha', 0.25, 'EdgeAlpha', 0.3);
        end
        if ~isempty(boundary_xy)
            for p = 1:n_wp
                nv = wp_nv(p); vi = wp_p(p, 1:nv);
                fill3(ax3, wp_n(vi, 1), wp_n(vi, 2), wp_n(vi, 3), ...
                      style.color.mesh_quad.fill, 'EdgeColor', style.color.mesh_quad.edge, ...
                      'FaceAlpha', 0.7, 'EdgeAlpha', 0.5);
            end
        end
        if ~isempty(wl_idx)
            plot3(ax3, verts(wl_idx, 1), verts(wl_idx, 2), ...
                  zeros(length(wl_idx), 1), '.', ...
                  'Color', style.fill_palette.waterline);
        end
        xlabel(ax3, '$x$ [m]');
        ylabel(ax3, '$y$ [m]');
        zlabel(ax3, '$z$ [m]');
        title(ax3, 'Hull + WP lid (combined)');
        view(ax3, [-40, 20]);
        axis(ax3, 'equal');
        mwecmass.output.figures.apply_axes_style(ax3, style);

        % FontWeight is not one of the grep-checked literal categories and is kept for emphasis.
        title(layout, 'BEM Mesh Diagnostic', 'FontWeight', 'bold');
        mwecmass.output.figures.apply_layout_style(layout, style);

        if opts.save.hydrodynamics.mesh_diagnostic
            saved = mwecmass.output.figures.export_figure( ...
                fig, 'WEC_Mesh_Diagnostic', config);
            fprintf('  Mesh diagnostic figure saved: %s\n', saved{1});
        end

    catch ME
        warning('mwecmass:figures:MeshDiag', ...
                'Mesh diagnostic plot failed: %s', ME.message);
    end
end
