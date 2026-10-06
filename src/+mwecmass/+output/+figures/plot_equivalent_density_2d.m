function plot_equivalent_density_2d(props, config, x_opt)
%PLOT_EQUIVALENT_DENSITY_2D Draw XZ cross-section colour-filled by strip-average equivalent density.
% Inputs: props (hull properties); config (density_nodes_z, shell, strip_edges, hull_z_min, shell_thickness); x_opt (vertical shift, core densities). Waterline at z=0; SI units.
% Fall back to continuous optimiser densities if props.realised_strip_density is unavailable.

    opts  = mwecmass.output.figures.output_options(config);
    style = opts.style;
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
            if isfield(props, 'fill_method') && ...
                    strcmp(props.fill_method, 'uhpc_fill')
                fig_name = 'WEC Equivalent Density (2D) — Realised UHPC';
            else
                fig_name = 'WEC Equivalent Density (2D) — Realised Steel';
            end
        else
            rho_source_is_realised = false;
        end

        fig = mwecmass.output.figures.new_figure(style, 'tall_double_column');
        fig.Name = fig_name;
        t = tiledlayout(fig, 1, 2);

        ax = nexttile(t);
        hold(ax, 'on'); axis(ax, 'equal');

        smooth_p2 = mwecmass.output.figures.build_silhouette_profile(config);
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
                rho_eq = mwecmass.output.figures.compute_strip_equivalent_density( ...
                             config, densities_core);
                rho_source_is_realised = false;
            end
        else
            rho_eq = mwecmass.output.figures.compute_strip_equivalent_density( ...
                         config, densities_core);
        end

        cmap  = mwecmass.output.figures.figure_colormap(style, 256);
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
        mwecmass.output.figures.draw_density_strips( ...
            ax, shifted_profile, node_z_wl, rho_eq, cmap, [d_min, d_max], sb3, style);

        % Overlay inner offset spline to show shell boundary
        if ~isempty(config.shell)
            inner_prof = mwecmass.output.figures.compute_inner_profile( ...
                shifted_profile, config.shell_thickness);
            mwecmass.output.figures.plot_inner_spline(ax, inner_prof, style);
        end

        % Outer hull boundary
        h_hull = plot(ax, px, pz, '-', 'Color', style.fill_palette.boundary);
        mwecmass.output.figures.style_line(h_hull, style, 'boundary');

        x_hull_range = [min(px) - 0.2, max(px) + 0.2];
        h_wl  = plot(ax, x_hull_range, [0 0], '--', 'Color', style.fill_palette.waterline);
        mwecmass.output.figures.style_line(h_wl, style, 'reference');

        % Wall boundary (constructability realisation type: modular_precast)
        h_wall_eq = [];
        if isfield(config, 'enable_constructability') && config.enable_constructability && ...
                isfield(config, 'constructability_wall_height')
            wall_z_wl = config.hull_z_max - config.constructability_wall_height + draft_final;
            h_wall_eq = plot(ax, x_hull_range, [wall_z_wl wall_z_wl], '-', ...
                'Color', style.fill_palette.wall_boundary);
            mwecmass.output.figures.style_line(h_wall_eq, style, 'boundary');
        end

        h_cg  = plot(ax, props.CG_total(1), props.CG_total(3), 'o', ...
                     'Color', style.fill_palette.boundary, 'MarkerFaceColor', style.color.cg);
        h_cb  = plot(ax, props.CB(1), props.CB(3), 's', ...
                     'Color', style.fill_palette.boundary, 'MarkerFaceColor', style.color.cb);
        mwecmass.output.figures.style_line(h_cg, style, 'reference');
        mwecmass.output.figures.style_line(h_cb, style, 'reference');

        colormap(ax, cmap);
        cbar = colorbar(ax);
        % $...$-wrapped: style_colorbar now styles cbar.Label through style_text, which forces
        % the LaTeX interpreter (previously unset, defaulting to MATLAB's 'tex', which does not
        % require the math delimiters LaTeX does for \rho/_/^).
        cbar.Label.String = 'Density, $\rho_{eq}$ (kg/m$^3$)';
        mwecmass.output.figures.style_colorbar(cbar, style);
        clim(ax, [d_min, d_max]);

        xlabel(ax, '$x$ [m]');
        ylabel(ax, '$z$ [m]');
        title(ax, 'Equivalent Uniform Density per Strip');
        handles_eq = [h_wl, h_cg, h_cb];
        labels_eq  = {'Waterline', 'CG', 'CB'};
        if ~isempty(h_wall_eq)
            handles_eq = [h_wl, h_wall_eq, h_cg, h_cb];
            labels_eq  = {'Waterline', 'Wall boundary', 'CG', 'CB'};
        end
        lg = legend(ax, handles_eq, labels_eq, 'Location', 'bestoutside');
        mwecmass.output.figures.style_legend(lg, style);
        mwecmass.output.figures.apply_axes_style(ax, style);
        hold(ax, 'off');

        % --- Properties card ---
        % Free-standing text() routed through style_text. The 'T_heave',
        % 'rho_core', 'rho_eq' strings below carry real underscores; forcing the LaTeX
        % interpreter (style_text always does, in place of this card's former 'none') needs them
        % escaped to '\_' so they still print literally rather than
        % render as a subscript trigger. The mono-font numeric rows keep style.mono_font_name via
        % style_text's optional font_name override, so the table stays column-aligned; escaping
        % 'rho_core'/'rho_eq' to 'rho\_core'/'rho\_eq' costs one extra character in that one
        % header cell, a minor known misalignment against the %-10s field width below it.
        ax2 = nexttile(t);
        axis(ax2, 'off');
        ty = 0.95; lh = 0.045;
        h_t1 = text(ax2, 0.05, ty, 'EQUIVALENT DENSITY TABLE', 'FontWeight', 'bold');
        mwecmass.output.figures.style_text(h_t1, style, 'annotation');
        ty = ty - 1.5*lh;
        h_t2 = text(ax2, 0.05, ty, sprintf('Draft: %.3f m', abs(config.hull_z_min + draft_final)));
        mwecmass.output.figures.style_text(h_t2, style, 'annotation');
        ty = ty - lh;
        h_t3 = text(ax2, 0.05, ty, escape_latex_local(sprintf('T_heave: %.2f s  |  T_pitch: %.2f s', ...
             props.periods.heave, props.periods.pitch)));
        mwecmass.output.figures.style_text(h_t3, style, 'annotation');
        ty = ty - 2*lh;
        h_t4 = text(ax2, 0.05, ty, ...
             escape_latex_local(sprintf('%-5s  %-8s  %-10s  %-10s', ...
                     'Strip', 'Z [m]', 'rho_core', 'rho_eq')), 'FontWeight', 'bold');
        mwecmass.output.figures.style_text(h_t4, style, 'annotation', style.mono_font_name);
        ty = ty - lh;
        for i = 1:length(densities_core)
            z_node = config.density_nodes_z(i) + draft_final;
            h_row = text(ax2, 0.05, ty, ...
                 sprintf('  %-3d  %+6.3f m  %6.0f kg/m3  %6.0f kg/m3', ...
                         i, z_node, densities_core(i), rho_eq(i)));
            mwecmass.output.figures.style_text(h_row, style, 'annotation', style.mono_font_name);
            ty = ty - lh;
            if ty < 0.05, break; end
        end

        if rho_source_is_realised
            if isfield(props, 'fill_method') && ...
                    strcmp(props.fill_method, 'uhpc_fill')
                sg2_str = 'Realised UHPC + Void: as-built strip $\rho_{eq}$';
            else
                sg2_str = 'Realised Steel-Fill: as-built strip $\rho_{eq}$';
            end
        else
            sg2_str = 'Steiner-Equivalent Bulk Density (optimiser)';
        end
        title(t, sg2_str);
        mwecmass.output.figures.apply_layout_style(t, style);

        if opts.save.stage1.density_2d
            saved = mwecmass.output.figures.export_figure( ...
                fig, 'WEC_Equivalent_Density_2D', config);
            fprintf('  Equivalent-density (2D) figure saved: %s\n', saved{1});
        end

    catch ME
        warning('mwecmass:figures:EquivalentDensity2DFailed', ...
                '2D equivalent visualization failed: %s', ME.message);
    end
end

function txt = escape_latex_local(txt)
% ESCAPE_LATEX_LOCAL  Escape the LaTeX-reserved characters this file's properties-card strings
%   can contain (%, _) so they print literally once style_text forces the LaTeX interpreter,
%   applied uniformly to every properties-card string this file feeds to text().
    txt = strrep(txt, '%', '\%');
    txt = strrep(txt, '_', '\_');
end
