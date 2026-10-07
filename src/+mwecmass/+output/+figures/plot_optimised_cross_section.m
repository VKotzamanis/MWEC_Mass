function plot_optimised_cross_section(props_3d, config, x_opt, realised)
%PLOT_OPTIMISED_CROSS_SECTION Plot the hull y = 0 section colour-filled by the density of each module.
% Inputs: props_3d (vertical_shift, realised_strip_density when realised, CG_total, CB); config
% (density_nodes_z, strip_edges, hull_solid, output options); x_opt (draft, Stage-2 densities);
% realised (optional, results.stage3). The hull outline is the exact y = 0 section of the outer
% patches (config.hull_solid). With realised, the outline, the voids (dashed), the ballast level and
% the status come from the sections of the realised solid (realised_section_data) and the fill is
% the realised density of each module; without it the fill is the Stage-2 density. Waterline at
% z = 0; SI units.

    if nargin < 4
        realised = [];
    end
    style = mwecmass.output.figures.presentation_style(config);
    try
        fig = mwecmass.output.figures.new_figure(style, 'tall_single');
        t = tiledlayout(fig, 1, 1);

        if isfield(props_3d, 'vertical_shift') && ~isnan(props_3d.vertical_shift)
            draft_final = props_3d.vertical_shift;
        else
            draft_final = x_opt(1);
        end
        if isfield(props_3d, 'realised_strip_density') && ~isempty(props_3d.realised_strip_density)
            densities_final = props_3d.realised_strip_density(:);
        else
            densities_final = x_opt(2:end);
        end

        voids = [];
        if ~isempty(realised)
            data = mwecmass.output.figures.realised_section_data(realised);
            if ~isempty(data.omitted)
                warning('mwecmass:figures:SectionsOmitted', ...
                    '%d section heights had no usable section and are left out of the figure (first: z = %g, %s).', ...
                    numel(data.omitted), data.omitted(1).z, data.omitted(1).reason);
            end
            profile_body = data.outline.profile;
            voids = data.polygons(strcmp({data.polygons.role}, 'void'));
        else
            if ~isfield(config, 'hull_solid') || isempty(config.hull_solid)
                error('mwecmass:figures:NoHullSolid', ...
                    'config.hull_solid (the exact hull geometry) is required to draw the hull outline.');
            end
            profile_body = mwecmass.output.figures.hull_outline_data(config.hull_solid).profile;
        end
        px = profile_body(:, 1);
        pz = profile_body(:, 2) + draft_final;
        shifted_profile = [px, pz];
        node_z_shifted = config.density_nodes_z + draft_final;

        ax1 = nexttile(t);
        hold(ax1, 'on');

        cmap  = mwecmass.output.figures.figure_colormap(style, 256);
        d_min = min(densities_final);
        d_max = max(densities_final);
        if d_max <= d_min, d_max = d_min + 1; end

        if isfield(config, 'strip_edges') && ~isempty(config.strip_edges)
            sb2 = config.strip_edges(:) + draft_final;
        else
            sb2 = [];
        end
        mwecmass.output.figures.draw_density_strips( ...
            ax1, shifted_profile, node_z_shifted, densities_final, cmap, [d_min, d_max], sb2, style);

        h_void = [];
        for k = 1:numel(voids)
            h_v = plot(ax1, voids(k).xz([1:end, 1], 1), voids(k).xz([1:end, 1], 2) + draft_final, '--', ...
                'Color', style.fill_palette.inner_boundary, 'LineWidth', style.line_width.boundary);
            if isempty(h_void)
                h_void = h_v;
            else
                set(h_v, 'HandleVisibility', 'off');
            end
        end

        h_hull = plot(ax1, px([1:end, 1]), pz([1:end, 1]), '-', 'Color', style.fill_palette.boundary);
        mwecmass.output.figures.style_line(h_hull, style, 'boundary');

        colormap(ax1, cmap);
        cb = colorbar(ax1);
        cb.Label.String = 'Uniform Segment Density, $\rho_{segment}$ [kg/m$^{3}$]';
        mwecmass.output.figures.style_colorbar(cb, style);
        clim(ax1, [d_min, d_max]);

        x_hull_range = [min(px) - 0.2, max(px) + 0.2];
        h_wl = plot(ax1, x_hull_range, [0 0], '--', 'Color', style.fill_palette.waterline);
        mwecmass.output.figures.style_line(h_wl, style, 'reference');
        h_wl.LineWidth = 2 * style.line_width.reference;

        h_zb = [];
        if ~isempty(realised) && realised.design.z_ballast > config.hull_z_min && ...
                realised.design.z_ballast < config.hull_z_max
            h_zb = plot(ax1, x_hull_range, realised.design.z_ballast * [1 1] + draft_final, '-', ...
                'Color', style.fill_palette.ballast_level);
            mwecmass.output.figures.style_line(h_zb, style, 'reference');
        end

        % Stage-2 wall boundary of the modular-precast mode; a realised design shows its solid
        % modules as built.
        if isempty(realised) && isfield(config, 'enable_constructability') && config.enable_constructability && ...
                isfield(config, 'constructability_wall_height')
            wall_z_wl = config.hull_z_max - config.constructability_wall_height + draft_final;
            h_wall = plot(ax1, x_hull_range, [wall_z_wl wall_z_wl], '-', 'Color', style.fill_palette.wall_boundary);
            mwecmass.output.figures.style_line(h_wall, style, 'boundary');
        end

        h_cg = plot(ax1, props_3d.CG_total(1), props_3d.CG_total(3), 'o', ...
            'Color', style.fill_palette.boundary, 'MarkerFaceColor', style.color.cg);
        h_cb = plot(ax1, props_3d.CB(1), props_3d.CB(3), 's', ...
            'Color', style.fill_palette.boundary, 'MarkerFaceColor', style.color.cb);
        mwecmass.output.figures.style_line(h_cg, style, 'reference');
        mwecmass.output.figures.style_line(h_cb, style, 'reference');

        handles = [h_wl, h_cg, h_cb];
        labels = {'Waterline, $Z = 0$ m', 'Center of Gravity, $Z_{CG}$', 'Center of Buoyancy, $Z_{CB}$'};
        if ~isempty(h_zb)
            handles(end + 1) = h_zb;
            labels{end + 1} = 'Ballast level, $z_{\mathrm{ballast}}$';
        end
        if ~isempty(h_void)
            handles(end + 1) = h_void;
            labels{end + 1} = 'Void boundary of the realised solid';
        end
        lg1 = legend(ax1, handles, labels, 'Location', 'north', 'NumColumns', 2);
        mwecmass.output.figures.style_legend(lg1, style);
        axis(ax1, 'equal');
        mwecmass.output.figures.apply_axes_style(ax1, style);
        set(ax1, 'Box', 'on', 'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off');

        if ~isempty(realised)
            if strcmp(realised.mode, 'modular_precast')
                sg_str = 'Stage 3: Modular Precast Construction';
                name_stem = 'WEC_Final_3D_CrossSection_UHPC';
            else
                sg_str = 'Stage 3: Thin Shell Construction';
                name_stem = 'WEC_Final_3D_CrossSection_Steel';
            end
            if strcmp(data.status, 'accepted')
                colour = style.status_palette.ok;
            else
                colour = style.status_palette.bad;
            end
            h_status = text(ax1, 0.02, 0.02, strrep(strrep(data.status_lines, '_', '\_'), '%', '\%'), ...
                'Units', 'normalized', 'VerticalAlignment', 'bottom', 'BackgroundColor', 'w', 'Margin', 1, ...
                'Color', colour, 'FontWeight', 'bold');
            mwecmass.output.figures.style_text(h_status, style, 'annotation');
        else
            sg_str = 'Stage 2: Material Unaware Solution';
            name_stem = 'WEC_Final_3D_CrossSection';
        end
        hold(ax1, 'off');
        xlabel(ax1, {'X [m]', sg_str});
        ylabel(ax1, 'Z [m]');
        fig.Name = sg_str;
        mwecmass.output.figures.apply_layout_style(t, style);

        saved = mwecmass.output.figures.export_figure(fig, name_stem, config);
        fprintf('  Cross-section figure saved: %s\n', saved{1});
    catch ME
        warning('mwecmass:figures:CrossSectionFailed', '3D cross-section failed: %s', ME.message);
    end
end
