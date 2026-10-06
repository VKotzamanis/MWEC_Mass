function plot_optimised_cross_section(props_3d, config, x_opt)
%PLOT_OPTIMISED_CROSS_SECTION Plot optimised hull XZ cross-section colour-filled by layer density.
% Inputs: props_3d (vertical_shift, realised_strip_density if post-realisation, fill_method); config (density_nodes_z, shell, constructability settings); x_opt (draft, core densities). Waterline at z=0; SI units.

    style = mwecmass.output.figures.presentation_style(config);
    try
        % Single panel: the second tile used to hold a numeric "properties card" duplicating the
        % console report; that panel is removed -- the numbers it showed are
        % written to Output/<type>/final_results.log instead. 'tall_single' replaces
        % 'tall_double_column', which was sized for the two-tile layout.
        fig = mwecmass.output.figures.new_figure(style, 'tall_single');
        t = tiledlayout(fig, 1, 1);

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

        % Stored profiles contain body-frame XZ points. Their historical
        % centroid-angle ordering crosses the concave neck of this hull;
        % order a local copy along each side before any plotting consumer.
        if isfield(props_3d, 'cross_section') && ~isempty(props_3d.cross_section)
            profile_body = props_3d.cross_section;
        elseif isfield(config, 'profile') && ~isempty(config.profile)
            profile_body = config.profile;
        else
            error('mwecmass:figures:CrossSectionProfileMissing', ...
                'A stored body-frame XZ cross_section or config.profile is required.');
        end
        profile_body = ordered_plot_profile(profile_body);
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

        if ~isempty(config.shell)
            % COMPOSITE: hatched shell annulus + density-coloured core
            inner_prof = mwecmass.output.figures.compute_inner_profile( ...
                shifted_profile, config.shell_thickness);
            mwecmass.output.figures.draw_composite_strips( ...
                ax1, shifted_profile, inner_prof, node_z_shifted, ...
                densities_final, [d_min, d_max], style);
            mwecmass.output.figures.plot_inner_spline(ax1, inner_prof, style);
        elseif isfield(config,'enable_constructability') && ...
                config.enable_constructability && ...
                isfield(config,'constructability_t_min')
            % Density strips + dashed inner void boundary at t_min offset
            if isfield(config, 'strip_edges') && ~isempty(config.strip_edges)
                sb2 = config.strip_edges(:) + draft_final;
            else
                sb2 = [];
            end
            mwecmass.output.figures.draw_density_strips( ...
                ax1, shifted_profile, node_z_shifted, densities_final, ...
                cmap, [d_min, d_max], sb2, style);
            inner_prof = mwecmass.output.figures.compute_inner_profile( ...
                shifted_profile, config.constructability_t_min);
            mwecmass.output.figures.plot_inner_spline(ax1, inner_prof, style);
        else
            % Compute shifted strip bounds if available
            if isfield(config, 'strip_edges') && ~isempty(config.strip_edges)
                sb2 = config.strip_edges(:) + draft_final;
            else
                sb2 = [];
            end
            mwecmass.output.figures.draw_density_strips( ...
                ax1, shifted_profile, node_z_shifted, densities_final, ...
                cmap, [d_min, d_max], sb2, style);
        end

        % Hull outline on top
        h_hull = plot(ax1, px, pz, '-', 'Color', style.fill_palette.boundary);
        mwecmass.output.figures.style_line(h_hull, style, 'boundary');

        colormap(ax1, cmap);
        cb = colorbar(ax1);
        % $...$-wrapped: style_colorbar now styles cb.Label through style_text, which forces the
        % LaTeX interpreter (previously unset, defaulting to MATLAB's 'tex', which does not
        % require the math delimiters LaTeX does for \rho/_/^).
        cb.Label.String = 'Uniform Segment Density, $\rho_{segment}$ [kg/m$^{3}$]';
        mwecmass.output.figures.style_colorbar(cb, style);
        clim(ax1, [d_min, d_max]);

        % Waterline — extend to full hull profile width
        x_hull_range = [min(px) - 0.2, max(px) + 0.2];
        h_wl = plot(ax1, x_hull_range, [0 0], '--', 'Color', style.fill_palette.waterline);
        mwecmass.output.figures.style_line(h_wl, style, 'reference');
        h_wl.LineWidth = 2 * style.line_width.reference;

        % Wall boundary (constructability realisation type: modular_precast)
        h_wall = [];
        if ~(rho_source_is_realised && isfield(props_3d, 'fill_method') && ...
                strcmp(props_3d.fill_method, 'uhpc_fill')) && ...
                isfield(config, 'enable_constructability') && config.enable_constructability && ...
                isfield(config, 'constructability_wall_height')
            wall_z_wl = config.hull_z_max - config.constructability_wall_height + draft_final;
            h_wall = plot(ax1, x_hull_range, [wall_z_wl wall_z_wl], '-', ...
                'Color', style.fill_palette.wall_boundary);
            mwecmass.output.figures.style_line(h_wall, style, 'boundary');
        end

        h_cg = plot(ax1, props_3d.CG_total(1), props_3d.CG_total(3), 'o', ...
            'Color', style.fill_palette.boundary, 'MarkerFaceColor', style.color.cg);
        h_cb = plot(ax1, props_3d.CB(1), props_3d.CB(3), 's', ...
            'Color', style.fill_palette.boundary, 'MarkerFaceColor', style.color.cb);
        mwecmass.output.figures.style_line(h_cg, style, 'reference');
        mwecmass.output.figures.style_line(h_cb, style, 'reference');

        xlabel(ax1, 'X [m]');
        ylabel(ax1, 'Z [m]');
        if rho_source_is_realised && isfield(props_3d, 'fill_method') && ...
                strcmp(props_3d.fill_method, 'uhpc_fill')
            h_construct_proxy = plot(ax1, nan, nan, '--', ...
                'Color', style.fill_palette.inner_boundary, ...
                'LineWidth', style.line_width.boundary);
            t_min_mm = config.constructability_t_min * 1000;
            lg1 = legend(ax1, [h_wl, h_cg, h_cb, h_construct_proxy], ...
                {'Waterline, $Z = 0$ m', 'Center of Gravity, $Z_{CG}$', ...
                 'Center of Buoyancy, $Z_{CB}$', ...
                 sprintf('Constructability, $t_{min}$ = %.3g [mm]', t_min_mm)}, ...
                'Location', 'north', 'NumColumns', 2);
        else
            lg1 = legend(ax1, [h_wl, h_cg, h_cb], ...
                {'Waterline, $Z = 0$ m', 'Center of Gravity, $Z_{CG}$', ...
                 'Center of Buoyancy, $Z_{CB}$'}, 'Location', 'north', ...
                'NumColumns', 2);
        end
        mwecmass.output.figures.style_legend(lg1, style);
        axis(ax1, 'equal');
        mwecmass.output.figures.apply_axes_style(ax1, style);
        set(ax1, 'Box', 'on', 'XGrid', 'off', 'YGrid', 'off', ...
            'XMinorGrid', 'off', 'YMinorGrid', 'off');
        hold(ax1, 'off');

        if rho_source_is_realised
            if isfield(props_3d, 'fill_method') && ...
                    strcmp(props_3d.fill_method, 'uhpc_fill')
                sg_str = 'Stage 2: Modular Precast Construction';
                name_stem = 'WEC_Final_3D_CrossSection_UHPC';
            else
                sg_str = 'Stage 2: Thin Shell Construction';
                name_stem = 'WEC_Final_3D_CrossSection_Steel';
            end
        else
            sg_str = 'Stage 2: Material Unaware Solution';
            name_stem = 'WEC_Final_3D_CrossSection';
        end
        xlabel(ax1, {'X [m]', sg_str});
        fig.Name = sg_str;
        mwecmass.output.figures.apply_layout_style(t, style);

        saved = mwecmass.output.figures.export_figure(fig, name_stem, config);
        fprintf('  Cross-section figure saved: %s\n', saved{1});
    catch ME
        warning('mwecmass:figures:CrossSectionFailed', '3D cross-section failed: %s', ME.message);
    end
end

function profile = ordered_plot_profile(points)
%ORDERED_PLOT_PROFILE Traverse the stored symmetric, single-valued hull sides.
% Reordering and exact deduplication only: coordinates remain unchanged.
% Unsupported profiles are rejected instead of inventing a new boundary.
    profile_error = 'mwecmass:figures:CrossSectionProfileInvalid';
    if ~isnumeric(points) || ~isreal(points) || ~ismatrix(points) || ...
            size(points, 2) ~= 2 || any(~isfinite(points(:)))
        error(profile_error, 'The stored hull profile must be finite real numeric Nx2 [X,Z] data.');
    end
    points = unique(points, 'rows', 'stable');
    if size(points, 1) < 6
        error(profile_error, 'The stored hull profile has too few distinct boundary points.');
    end

    right = sortrows(points(points(:, 1) > 0, :), 2);
    left = sortrows(points(points(:, 1) < 0, :), 2);
    centre = sortrows(points(points(:, 1) == 0, :), 2);
    if size(right, 1) < 2 || size(left, 1) < 2 || size(centre, 1) ~= 2
        error(profile_error, ...
            'Expected two hull sides and exactly two centreline endpoints (keel and top).');
    end
    if any(diff(right(:, 2)) <= 0) || any(diff(left(:, 2)) <= 0)
        error(profile_error, ...
            'Each stored hull side must have a single boundary point at each Z level.');
    end
    if centre(1, 2) >= min([right(:, 2); left(:, 2)]) || ...
            centre(2, 2) <= max([right(:, 2); left(:, 2)])
        error(profile_error, 'Centreline endpoints must bound both hull sides in Z.');
    end
    tolerance = 1e-10 * max(1, max(abs(points(:))));
    if size(right, 1) ~= size(left, 1) || ...
            any(abs(right(:, 1) + left(:, 1)) > tolerance) || ...
            any(abs(right(:, 2) - left(:, 2)) > tolerance)
        error(profile_error, ...
            'Stored hull sides are not symmetric paired XZ samples; ordering is ambiguous.');
    end

    % Right side ascends from keel to top; left side descends back to keel.
    % Strictly monotone sides of opposite X sign cannot cross one another.
    % Endpoint checks above make the top and keel joins local continuations.
    profile = [centre(1, :); right; centre(2, :); flipud(left); centre(1, :)];
    edge_lengths = hypot(diff(profile(:, 1)), diff(profile(:, 2)));
    signed_area = 0.5 * sum(profile(1:end-1, 1) .* profile(2:end, 2) - ...
        profile(2:end, 1) .* profile(1:end-1, 2));
    if any(~isfinite(edge_lengths)) || any(edge_lengths <= 0) || ...
            ~isfinite(signed_area) || signed_area <= 0
        error(profile_error, 'The ordered hull boundary has invalid adjacency or closure.');
    end
end
