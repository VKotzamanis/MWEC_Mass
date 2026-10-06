function plot_panel_normals(mesh, options)
%PLOT_PANEL_NORMALS 3D hull mesh with outward panel normals as quivers plus n_z histogram.
% Renders BEM panel mesh with outward normals (HAMS-MREL convention) at panel
% centroids; left panel shows 3D view, right panel shows z-component distribution
% to detect orientation errors (spike at +1 → wrong sign, all zero → degenerate).
% options.scale [m] (default 0.12×hull_dia), subsample, face_alpha (0.30),
% normal_color, face_color, view_az (−35°), view_el (25°), output struct optional.

    try
        if nargin < 2, options = struct(); end
        style = mwecmass.output.figures.presentation_style(options);
        if ~isfield(options, 'face_alpha'),   options.face_alpha   = 0.30;                    end
        if ~isfield(options, 'face_color'),   options.face_color   = style.fill_palette.shell; end
        if ~isfield(options, 'normal_color'), options.normal_color = style.color.diagnostic;   end
        if ~isfield(options, 'subsample'),    options.subsample    = 1;                        end
        if ~isfield(options, 'view_az'),      options.view_az      = -35;                      end
        if ~isfield(options, 'view_el'),      options.view_el      =  25;                      end

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
        fig = mwecmass.output.figures.new_figure(style, 'double_column');
        set(fig, 'Name', 'Panel Normals Diagnostic');
        layout = tiledlayout(fig, 1, 2);

        % ── Panel 1 — 3D hull + quiver arrows ─────────────────────
        ax1 = nexttile(layout, 1);
        hold(ax1, 'on');

        for p = 1:nP
            v = panels(p, :);
            if v(3) == v(4), nv = 3; else, nv = 4; end
            vi = v(1:nv);
            fill3(ax1, verts(vi,1), verts(vi,2), verts(vi,3), ...
                  options.face_color, ...
                  'EdgeColor', style.color.mesh_edge, ...
                  'FaceAlpha', options.face_alpha, ...
                  'EdgeAlpha', 0.55);
        end

        hq = quiver3(ax1, ...
            centroids(idx_plot,1), centroids(idx_plot,2), centroids(idx_plot,3), ...
            normals(idx_plot,1) * options.scale, ...
            normals(idx_plot,2) * options.scale, ...
            normals(idx_plot,3) * options.scale, ...
            0, 'Color', options.normal_color, 'MaxHeadSize', 0.55);
        mwecmass.output.figures.style_line(hq, style, 'curve');

        view(ax1, options.view_az, options.view_el);
        axis(ax1, 'equal');
        xlabel(ax1, '$x$ [m]');
        ylabel(ax1, '$y$ [m]');
        zlabel(ax1, '$z$ [m]');
        title(ax1, sprintf('Hull mesh + outward normals (%d panels, %d shown)', nP, n_show));
        mwecmass.output.figures.apply_axes_style(ax1, style);

        % ── Panel 2 — n_z histogram ────────────────────────────────
        ax2 = nexttile(layout, 2);
        hold(ax2, 'on');

        nz_all = normals(:, 3);
        histogram(ax2, nz_all, 50, ...
                  'FaceColor', options.normal_color, ...
                  'EdgeColor', 'none', ...
                  'FaceAlpha', 0.75);
        hz = xline(ax2, 0, '--', 'Color', style.color.reference, 'Label', '$n_z=0$');
        mwecmass.output.figures.style_line(hz, style, 'reference');
        mwecmass.output.figures.style_text(hz, style, 'label');

        n_pos  = sum(nz_all >  0.01);
        n_neg  = sum(nz_all < -0.01);
        n_zero = sum(abs(nz_all) <= 0.01);

        if n_zero > 0
            % The warning colour is the status_palette.bad role (a failed/degenerate condition),
            % passed to LaTeX's \color[rgb]{} form so no colour name is literal.
            h_deg = text(ax2, 0.02, 0.96, ...
                 sprintf('\\color[rgb]{%.3f,%.3f,%.3f}%d degenerate panels ($|n_z| \\leq 0.01$)', ...
                         style.status_palette.bad, n_zero), ...
                 'Units', 'normalized', 'VerticalAlignment', 'top');
            mwecmass.output.figures.style_text(h_deg, style, 'annotation');
        end

        xlabel(ax2, '$n_z$ (normal $z$-component)');
        ylabel(ax2, 'Panel count');
        title(ax2, sprintf( ...
            '$n_z>0$: %d,  $n_z<0$: %d,  $|n_z|\\leq0.01$: %d', ...
            n_pos, n_neg, n_zero));
        mwecmass.output.figures.apply_axes_style(ax2, style);

        title(layout, sprintf('BEM Panel Normals  (%d panels, quiver scale = %.3f m)', ...
                        nP, options.scale));
        mwecmass.output.figures.apply_layout_style(layout, style);

        saved = mwecmass.output.figures.export_figure(fig, 'WEC_PanelNormals', options);
        fprintf('  Panel normals figure saved: %s\n', saved{1});

    catch ME
        warning('mwecmass:figures:PanelNormalsFailed', ...
                'Panel normals plot failed: %s', ME.message);
    end
end
