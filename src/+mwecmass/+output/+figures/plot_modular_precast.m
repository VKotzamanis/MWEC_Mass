function plot_modular_precast(cstr, config)
%PLOT_MODULAR_PRECAST Plot UHPC + void hull constructability: XZ elevation and per-strip plan-view cross-sections.
% Inputs: cstr (strip geometry, offsets, void volumes, fill method); config (material properties, visualization options). Void is perpendicular offset (t_strip), then volume-contracted to match V_void/V_total.

    style = mwecmass.output.figures.presentation_style(config);

    N = length(cstr.strip_z_lo);

    % Per-strip thickness (Phase 1b)
    if isfield(cstr, 't_offset_strip') && ~isempty(cstr.t_offset_strip)
        t_strip = cstr.t_offset_strip;
    else
        t_strip = cstr.t_UHPC * ones(N, 1);
    end
    if isfield(cstr, 'strip_is_solid') && ~isempty(cstr.strip_is_solid)
        is_solid = cstr.strip_is_solid(:);
    else
        is_solid = false(N, 1);
    end
    if isfield(cstr, 'strip_is_wall') && ~isempty(cstr.strip_is_wall)
        is_wall = cstr.strip_is_wall(:);
    else
        is_wall = false(N, 1);
    end

    vs       = cstr.vertical_shift;
    t_min_mm = cstr.t_min * 1000;

    % Smooth world-frame hull silhouette.  The shared builder returns an
    % already ordered profile when its geometry source is available; only a
    % raw config.profile fallback needs the plotting-only side ordering below.
    profile_is_raw = false;
    try
        prof_body = mwecmass.output.figures.build_silhouette_profile(config);
        if isfield(config, 'profile') && ~isempty(config.profile) && ...
                isequal(prof_body, config.profile)
            profile_is_raw = true;
        end
    catch
        if isfield(config, 'profile') && ~isempty(config.profile)
            prof_body = config.profile;
            profile_is_raw = true;
        else
            prof_body = [];
        end
    end
    if isempty(prof_body)
        warning('mwecmass:figures:plot_modular_precast:NoProfile', ...
                'No smooth profile available — falling back to outer contours.');
        prof_body = build_profile_from_outer_contours(cstr);
    end
    if profile_is_raw
        profile_error = 'mwecmass:figures:plot_modular_precast:ProfileInvalid';
        if ~isnumeric(prof_body) || ~isreal(prof_body) || ...
                ~ismatrix(prof_body) || size(prof_body, 2) ~= 2 || ...
                any(~isfinite(prof_body(:)))
            error(profile_error, ...
                'The fallback hull profile must be finite real numeric Nx2 [X,Z] data.');
        end
        points = unique(prof_body, 'rows', 'stable');
        if size(points, 1) < 6
            error(profile_error, 'The fallback hull profile has too few distinct boundary points.');
        end
        right  = sortrows(points(points(:, 1) > 0, :), 2);
        left   = sortrows(points(points(:, 1) < 0, :), 2);
        centre = sortrows(points(points(:, 1) == 0, :), 2);
        if size(right, 1) < 2 || size(left, 1) < 2 || size(centre, 1) ~= 2
            error(profile_error, ...
                'Expected two hull sides and exactly two centreline endpoints (keel and top).');
        end
        if any(diff(right(:, 2)) <= 0) || any(diff(left(:, 2)) <= 0)
            error(profile_error, ...
                'Each stored hull side must have a single boundary point at each Z level.');
        end
        tolerance = 1e-10 * max(1, max(abs(points(:))));
        if centre(1, 2) >= min([right(:, 2); left(:, 2)]) || ...
                centre(2, 2) <= max([right(:, 2); left(:, 2)]) || ...
                size(right, 1) ~= size(left, 1) || ...
                any(abs(right(:, 1) + left(:, 1)) > tolerance) || ...
                any(abs(right(:, 2) - left(:, 2)) > tolerance)
            error(profile_error, ...
                'Fallback hull profile sides are not valid symmetric paired XZ samples.');
        end
        % Exact approved local plotting order: right side bottom-to-top,
        % then left side top-to-bottom, closing at the keel.
        prof_body = [centre(1, :); right; centre(2, :); flipud(left); centre(1, :)];
        edge_lengths = hypot(diff(prof_body(:, 1)), diff(prof_body(:, 2)));
        signed_area = 0.5 * sum(prof_body(1:end-1, 1) .* prof_body(2:end, 2) - ...
            prof_body(2:end, 1) .* prof_body(1:end-1, 2));
        if any(~isfinite(edge_lengths)) || any(edge_lengths <= 0) || ...
                ~isfinite(signed_area) || signed_area <= 0
            error(profile_error, 'The ordered fallback hull boundary is invalid.');
        end
    end
    px_outer = prof_body(:, 1);
    pz_outer = prof_body(:, 2) + vs;

    %% ═════════════════════════════════════════════════════════
    %%  FIGURE 1 — XZ midplane elevation
    %% ═════════════════════════════════════════════════════════
    % A taller double-column canvas gives the composition legend, elevation,
    % caption, and east-side feature legend their own layout-managed space.
    fig1 = mwecmass.output.figures.new_figure(style, 'tall_double_column');
    set(fig1, 'Name', 'Constructability: XZ midplane');
    tl1 = tiledlayout(fig1, 4, 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
    ax_comp = nexttile(tl1, 1);
    axis(ax_comp, 'off');
    hold(ax_comp, 'on');
    ax1 = nexttile(tl1, 2, [3 1]);
    hold(ax1, 'on');

    x_min_hull = min(px_outer);
    x_max_hull = max(px_outer);

    % ── Inner offset polygons computed on the FULL profile ────
    % CRITICAL: offset the entire hull profile (NOT a strip-clipped
    % piece).  Offsetting a clipped piece adds the strip's horizontal
    % top/bottom edges to the polygon, and the offset of those fake
    % edges produces the "shelves" we keep seeing at strip
    % boundaries.  Offsetting the full profile gives a uniform
    % perpendicular wall thickness on the actual hull surface.
    %
    % Cache one inner polygon per unique thickness.  Use
    % polyshape/polybuffer when available (handles convex AND
    % concave regions correctly via Minkowski erosion); fall back
    % to plot_steel_solve's offset_polygon otherwise.
    inner_keys   = {};
    inner_polys  = {};
    for i = 1:N
        if is_wall(i) || is_solid(i), continue; end
        ti = t_strip(i);
        if ~isfinite(ti) || ti <= 0, continue; end
        key = sprintf('%.6f', ti);
        if any(strcmp(inner_keys, key)), continue; end
        inner_keys{end+1} = key;                                  %#ok<AGROW>
        inner_polys{end+1} = compute_inner_offset_local( ...
                                  px_outer, pz_outer, ti);        %#ok<AGROW>
    end
    inner_lookup = @(t) inner_polys{find(strcmp(inner_keys, ...
                                               sprintf('%.6f', t)), 1)};

    % ── Per-strip render ─────────────────────────────────────
    for i = 1:N
        z_lo_w = cstr.strip_z_lo(i) + vs;
        z_hi_w = cstr.strip_z_hi(i) + vs;

        % Clip the FULL outer profile to this strip
        [x_top, z_top]   = mwecmass.internal.clip_z(px_outer, pz_outer, z_lo_w, 'above');
        [x_clip, z_clip] = mwecmass.internal.clip_z(x_top, z_top,        z_hi_w, 'below');
        if length(x_clip) < 3, continue; end

        if is_wall(i)
            fc = style.fill_palette.solid_material;
        elseif is_solid(i)
            fc = style.fill_palette.fill_material;
        else
            fc = style.fill_palette.jacket_material;
        end
        patch(ax1, x_clip, z_clip, fc, ...
              'EdgeColor', 'none', 'FaceAlpha', 1.0, ...
              'HandleVisibility', 'off');

        % Annular strips: clip the precomputed inner-offset polygon
        % (full hull, perpendicular t) to this strip's z-range,
        % then HORIZONTALLY contract the resulting void polygon so
        % its 2D area equals the volume target
        %     A_void_target = A_strip_XZ * V_void(i)/V_total(i).
        % This accounts for the 3" shell already taking some UHPC
        % volume; the void shrinks symmetrically inward in x to
        % make room for the remaining UHPC required by the strip.
        if ~is_wall(i) && ~is_solid(i) && isfinite(t_strip(i)) && t_strip(i) > 0
            P_in = inner_lookup(t_strip(i));
            if ~isempty(P_in) && size(P_in, 1) >= 3
                [xi_top, zi_top]   = mwecmass.internal.clip_z( ...
                                        P_in(:,1), P_in(:,2), z_lo_w, 'above');
                [xi_clip, zi_clip] = mwecmass.internal.clip_z( ...
                                        xi_top, zi_top, z_hi_w, 'below');
                if length(xi_clip) >= 3
                    % Volume-based target void area in this strip
                    V_u = 0; V_v = 0;
                    if ~isempty(cstr.strip_V_UHPC), V_u = cstr.strip_V_UHPC(i); end
                    if ~isempty(cstr.strip_V_void), V_v = cstr.strip_V_void(i); end
                    V_t = V_u + V_v;
                    A_strip_XZ = polyarea(x_clip, z_clip);
                    A_void_offset = polyarea(xi_clip, zi_clip);
                    if V_t > 1e-12 && A_strip_XZ > 1e-9 && A_void_offset > 1e-9
                        A_void_target = A_strip_XZ * (V_v / V_t);
                        if A_void_offset > A_void_target
                            % Horizontal scale toward x-centroid
                            s_h = A_void_target / A_void_offset;
                            cx  = mean(xi_clip);
                            xi_clip = (xi_clip - cx) * s_h + cx;
                        end
                    end
                    h_void_strip = patch(ax1, xi_clip, zi_clip, style.fill_palette.void, ...
                          'EdgeColor', style.fill_palette.inner_boundary, ...
                          'LineStyle', '--', ...
                          'HandleVisibility', 'off');
                    mwecmass.output.figures.style_line(h_void_strip, style, 'boundary');
                    mwecmass.output.figures.draw_hatch_strips(ax1, xi_clip, zi_clip, ...
                                     style.hatch_spacing, style.fill_palette.hatch);
                end
            end
        end
    end

    % Outer hull outline (over the patches)
    h_outer = plot(ax1, [px_outer; px_outer(1)], [pz_outer; pz_outer(1)], ...
         '-', 'Color', style.fill_palette.boundary, ...
         'HandleVisibility', 'off');
    mwecmass.output.figures.style_line(h_outer, style, 'boundary');

    % Reference lines and strip boundaries. Subtle grey dividers between strips: sourced from
    % style.line_width.grid (0.6 pt), the closest existing role by value and by intent (a thin
    % non-data guide line); style_line's role enum (curve/boundary/reference) has no divider
    % role and 'reference' would double this width to 1.2 pt, so these stay a direct set() call
    % rather than be misclassified into a role that changes their weight.
    x_range = [x_min_hull - 0.30, x_max_hull + 0.30];
    for i = 1:N
        z_b = cstr.strip_z_lo(i) + vs;
        line(ax1, x_range, [z_b z_b], 'Color', [0.55 0.55 0.55], ...
             'LineStyle', ':', 'LineWidth', style.line_width.grid, ...
             'HandleVisibility', 'off');
    end
    line(ax1, x_range, ...
         [cstr.strip_z_hi(N) + vs, cstr.strip_z_hi(N) + vs], ...
         'Color', [0.55 0.55 0.55], 'LineStyle', ':', 'LineWidth', style.line_width.grid, ...
         'HandleVisibility', 'off');

    h_wl = plot(ax1, x_range, [0 0], '--', ...
                'Color', style.fill_palette.waterline);
    mwecmass.output.figures.style_line(h_wl, style, 'reference');

    % Legend proxies
    h_uhpc = patch(ax1, NaN, NaN, style.fill_palette.jacket_material,       'EdgeColor', 'none');
    h_fs   = patch(ax1, NaN, NaN, style.fill_palette.fill_material, 'EdgeColor', 'none');
    h_wall = patch(ax1, NaN, NaN, style.fill_palette.solid_material,       'EdgeColor', 'none');
    h_void = patch(ax1, NaN, NaN, style.fill_palette.void, ...
                   'EdgeColor', style.fill_palette.inner_boundary, 'LineStyle', '--');
    mwecmass.output.figures.style_line(h_void, style, 'boundary');

    % Per-module annotation (right of the hull), preserving the existing
    % bottom-to-top strip numbering and avoiding an extra label for the
    % upper exterior boundary line.
    x_annot = x_range(2) + 0.10;
    for i = 1:N
        z_mid = 0.5 * (cstr.strip_z_lo(i) + cstr.strip_z_hi(i)) + vs;
        h_lbl = text(ax1, x_annot, z_mid, sprintf('Module %d', i), ...
             'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
             'Color', [0.15 0.15 0.15], 'Clipping', 'on');
        mwecmass.output.figures.style_text(h_lbl, style, 'annotation');
        h_lbl.FontSize = 8;
    end

    rho_material_txt = num2str(cstr.rho_hull, '%.0f');
    shell_thickness_txt = num2str(t_min_mm, '%.1f');
    air_density_txt = num2str(cstr.rho_fill, '%.1f');
    xlabel(ax1, {'X [m]', ...
        ['Stage 3: Modular Precast Construction, $\rho_{\mathrm{material}} = ' ...
         rho_material_txt '$ [kg/m$^3$]']});
    ylabel(ax1, 'Z [m]');

    lg1 = legend(ax1, [h_uhpc, h_fs, h_wall, h_void, h_wl], ...
           {['Precast Shell, $t_{\min} = ' shell_thickness_txt '$ [mm]'], ...
            'Solid Module (Ballast)', ...
            'Solid Module (Wall)', ...
            ['Hollow Volume, $\rho_{\mathrm{Air}} = ' air_density_txt ...
             '$ [kg/m$^3$]'], ...
            'Waterline, $Z = 0$ m'}, ...
           'Location', 'eastoutside', 'AutoUpdate', 'off');
    mwecmass.output.figures.style_legend(lg1, style);
    lg1.FontSize = min(style.font_size.legend, 8.5);
    axis(ax1, 'equal');
    mwecmass.output.figures.apply_axes_style(ax1, style);
    set(ax1, 'TickDir', 'out', 'Box', 'on', ...
        'XGrid', 'off', 'YGrid', 'off', ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off');
    right_margin = max(1.25, 0.40 * (x_max_hull - x_min_hull));
    xlim(ax1, [x_range(1), x_annot + right_margin]);
    z_span = max(pz_outer) - min(pz_outer);
    z_pad = max(0.05 * z_span, eps(max(abs(pz_outer))));
    ylim(ax1, [min(pz_outer) - z_pad, max(pz_outer) + z_pad]);

    % The composition legend lives on the invisible first-tile axes.  Its
    % NaN-only proxies cannot affect the elevation's data limits.
    h_comp = gobjects(N, 1);
    comp_labels = cell(N, 1);
    for i = 1:N
        h_comp(i) = plot(ax_comp, NaN, NaN, 'LineStyle', 'none', ...
                         'Marker', 'none', 'HandleVisibility', 'on');
        V_u = 0;
        V_v = 0;
        if ~isempty(cstr.strip_V_UHPC), V_u = cstr.strip_V_UHPC(i); end
        if ~isempty(cstr.strip_V_void), V_v = cstr.strip_V_void(i); end
        if V_v > 1e-12
            ratio_pct = 100 * V_u / V_v;
            comp_labels{i} = ['Module ' num2str(i) ...
                              ': $V_{\mathrm{material}}/V_{\mathrm{air}} = ' ...
                              num2str(ratio_pct, '%.1f') '\%$'];
        elseif V_u > 1e-12
            comp_labels{i} = ['Module ' num2str(i) ...
                              ': $V_{\mathrm{material}}/V_{\mathrm{air}} = \infty\%$'];
        else
            comp_labels{i} = ['Module ' num2str(i) ...
                              ': $V_{\mathrm{material}}/V_{\mathrm{air}} = ' ...
                              '\mathrm{undefined}$'];
        end
    end
    lg_comp = legend(ax_comp, h_comp, comp_labels, ...
                     'Location', 'north', ...
                     'NumColumns', ceil(N / 2), ...
                     'Orientation', 'vertical', 'AutoUpdate', 'off');
    mwecmass.output.figures.style_legend(lg_comp, style);
    lg_comp.FontSize = 7.5;
    lg_comp.Box = 'off';
    lg_comp.ItemTokenSize = [6 5];
    hold(ax_comp, 'off');

    try
        saved_xz = mwecmass.output.figures.export_figure(fig1, ...
            'WEC_Constructability_XZ', config, ...
            mwecmass.output.output_dir('modular_precast'));
        fprintf('      Figure saved: %s\n', saved_xz{1});
    catch
    end

    %% ═════════════════════════════════════════════════════════
    %%  FIGURE 2 — Per-strip plan-view cross-sections
    %%  Plot the retained bottom (z = strip_z_lo) and top
    %%  (z = strip_z_hi) module views.  Degenerate endpoints are
    %%  skipped automatically; Strip 1 Top and structural-wall
    %%  views are excluded explicitly below.
    %%  Void is drawn as solid white (no hatching).
    %% ═════════════════════════════════════════════════════════

    % Build the (strip_idx, end_label, contour, z) panel list.  The first
    % strip's top view and both views of the structural-wall strip are
    % intentionally omitted; the existing area test still removes any
    % other degenerate endpoint contour.
    panels = struct('i', {}, 'module_no', {}, 'end_label', {}, ...
                    'contour', {}, 'z', {});
    min_area_keep = 1e-4;   % m^2 — drop degenerate end-cross-sections
    for i = 1:N
        conts = cstr.contours_outer{i};
        if isempty(conts), continue; end
        bot_c = conts{1};
        if ~is_wall(i) && ~isempty(bot_c) && size(bot_c, 1) >= 3 && ...
                polyarea(bot_c(:,1), bot_c(:,2)) > min_area_keep
            panels(end+1) = struct( ...
                'i', i, 'module_no', 0, 'end_label', 'Bottom', ...
                'contour', bot_c, 'z', cstr.strip_z_lo(i)); %#ok<AGROW>
        end
        top_c = conts{end};
        if i ~= 1 && ~is_wall(i) && ...
                ~isempty(top_c) && size(top_c, 1) >= 3 && ...
                polyarea(top_c(:,1), top_c(:,2)) > min_area_keep
            panels(end+1) = struct( ...
                'i', i, 'module_no', 0, 'end_label', 'Top', ...
                'contour', top_c, 'z', cstr.strip_z_hi(i)); %#ok<AGROW>
        end
    end

    % Number only retained modules, in the same bottom-to-top strip order
    % used by cstr.  Keep the original strip index in panels(p).i so all
    % geometry and volume lookups continue to use the authoritative data.
    if isempty(panels)
        retained_strip_ids = [];
    else
        retained_strip_ids = unique([panels.i], 'stable');
        for p = 1:length(panels)
            panels(p).module_no = find(retained_strip_ids == panels(p).i, 1);
        end
    end

    n_panels = length(panels);
    if n_panels == 0
        fprintf('      Figure 2 skipped: no usable strip cross-sections.\n');
    else
        % Modules 1--3 occupy columns 1--3, with top views above their
        % matching bottom views.
        n_cols = 3;
        tile_indices = zeros(1, n_panels);
        for p = 1:n_panels
            if strcmp(panels(p).end_label, 'Top')
                tile_indices(p) = panels(p).module_no;
            else
                tile_indices(p) = n_cols + panels(p).module_no;
            end
        end
        fig2 = mwecmass.output.figures.new_figure(style, 'tall_double_column');
        set(fig2, 'Name', 'Constructability: Strip Plan-View');
        tl2 = tiledlayout(fig2, 2, 3, ...
            'TileSpacing', 'compact', 'Padding', 'compact');
        title(tl2, 'Precast Construction: Plan-View of the Hollow Modules');
        mwecmass.output.figures.apply_layout_style(tl2, style);

        [~, panel_order] = sort(tile_indices);
        for q = 1:n_panels
            p = panel_order(q);
            ax = nexttile(tl2, tile_indices(p));
            hold(ax, 'on');
            info  = panels(p);
            i     = info.i;
            outer = info.contour;

            xo = outer(:, 1);  yo = outer(:, 2);
            if abs(xo(end)-xo(1)) > 1e-10 || abs(yo(end)-yo(1)) > 1e-10
                xo = [xo; xo(1)];  yo = [yo; yo(1)]; %#ok<AGROW> -- xo/yo are re-sliced fresh from outer(:,1:2) each panel iteration (two lines above); this closes the polygon with at most one extra point, not an accumulating loop.
            end

            if is_wall(i)
                fc = style.fill_palette.solid_material;
            elseif is_solid(i)
                fc = style.fill_palette.fill_material;
            else
                fc = style.fill_palette.jacket_material;
            end
            patch(ax, xo, yo, fc, ...
                  'EdgeColor', 'none', 'FaceAlpha', 1.0);

            % Void overlay — perpendicular offset by t_strip(i),
            % then 2D-contracted about its centroid until its
            % area matches A_outer_plan * V_void(i)/V_total(i).
            % Filled SOLID WHITE (no hatching).
            if ~is_wall(i) && ~is_solid(i) && ...
                    isfinite(t_strip(i)) && t_strip(i) > 0
                P_in = compute_inner_offset_local(xo, yo, t_strip(i));
                if ~isempty(P_in) && size(P_in, 1) >= 3
                    xi = P_in(:,1);  yi = P_in(:,2);
                    A_outer_plan = polyarea(xo, yo);
                    V_u = 0; V_v = 0;
                    if ~isempty(cstr.strip_V_UHPC), V_u = cstr.strip_V_UHPC(i); end
                    if ~isempty(cstr.strip_V_void), V_v = cstr.strip_V_void(i); end
                    V_t = V_u + V_v;
                    A_in_now = polyarea(xi, yi);
                    if V_t > 1e-12 && A_outer_plan > 1e-9 && A_in_now > 1e-9
                        A_void_target = A_outer_plan * (V_v / V_t);
                        if A_in_now > A_void_target
                            s_p = sqrt(A_void_target / A_in_now);
                            cx = mean(xi);  cy = mean(yi);
                            xi = (xi - cx) * s_p + cx;
                            yi = (yi - cy) * s_p + cy;
                        end
                    end
                    h_void_p = patch(ax, xi, yi, style.fill_palette.void, ...
                          'EdgeColor', style.fill_palette.inner_boundary, 'LineStyle', '--');
                    mwecmass.output.figures.style_line(h_void_p, style, 'boundary');
                end
            end

            h_outer_p = plot(ax, xo, yo, '-', 'Color', style.fill_palette.boundary);
            mwecmass.output.figures.style_line(h_outer_p, style, 'boundary');

            % The two endpoint views share the authoritative volume totals
            % integrated for their original strip.
            V_u = 0; V_v = 0;
            if ~isempty(cstr.strip_V_UHPC), V_u = cstr.strip_V_UHPC(i); end
            if ~isempty(cstr.strip_V_void), V_v = cstr.strip_V_void(i); end

            h_volume_material = plot(ax, NaN, NaN, '-', ...
                'Color', [0 0 0], 'LineWidth', style.line_width.boundary);
            h_volume_air = plot(ax, NaN, NaN, '--', ...
                'Color', [0.55 0.55 0.55], ...
                'LineWidth', style.line_width.boundary);
            lg_panel = legend(ax, [h_volume_material, h_volume_air], ...
                {sprintf('$V_{\\mathrm{Material}} = %.3g$ [m$^{3}$]', V_u), ...
                 sprintf('$V_{\\mathrm{Air}} = %.3g$ [m$^{3}$]', V_v)}, ...
                'Location', 'northoutside', 'NumColumns', 2, ...
                'Orientation', 'horizontal', 'AutoUpdate', 'off');
            mwecmass.output.figures.style_legend(lg_panel, style);
            set(lg_panel, 'FontSize', 6.5, 'ItemTokenSize', [6 5], ...
                'Box', 'on');

            x_pad = max(0.05 * (max(xo) - min(xo)), eps(max(abs(xo))));
            y_pad = max(0.05 * (max(yo) - min(yo)), eps(max(abs(yo))));
            xlim(ax, [min(xo) - x_pad, max(xo) + x_pad]);
            ylim(ax, [min(yo) - y_pad, max(yo) + y_pad]);
            axis(ax, 'equal');
            xlabel(ax, {'X [m]', sprintf('Module %d, %s View', ...
                                         info.module_no, info.end_label)});
            ylabel(ax, 'Y [m]');
            mwecmass.output.figures.apply_axes_style(ax, style, 'small_multiple');
            set(ax, 'Box', 'on', 'TickDir', 'out', ...
                'XGrid', 'off', 'YGrid', 'off', 'ZGrid', 'off', ...
                'XMinorGrid', 'off', 'YMinorGrid', 'off', ...
                'ZMinorGrid', 'off');
            hold(ax, 'off');
        end
    end

    if n_panels > 0
        try
            saved_strips = mwecmass.output.figures.export_figure(fig2, ...
                'WEC_Constructability_Strips', config, ...
                mwecmass.output.output_dir('modular_precast'));
            fprintf('      Figure saved: %s\n', saved_strips{1});
        catch
        end
    end
end

function txt = escape_latex_local(txt)
% ESCAPE_LATEX_LOCAL  Escape the LaTeX-reserved characters this file's generated strings can
%   contain (%, _) so they print literally once apply_axes_style/style_text force the LaTeX
%   interpreter, applied uniformly to every generated string in this file rather than only the
%   one title that first needed it.
    txt = strrep(txt, '%', '\%');
    txt = strrep(txt, '_', '\_');
end

function P_in = compute_inner_offset_local(x_outer, y_outer, t)
% COMPUTE_INNER_OFFSET_LOCAL  Inner offset of polygon (x_outer, y_outer)
%   by perpendicular distance t.
%
%   Method order (best→worst, falls through on failure):
%     1. polyshape + polybuffer(-t)   — Minkowski erosion, the
%        literature-standard for offsetting polygons with mixed convex/
%        concave regions (handles self-intersection correctly).  Available
%        in MATLAB R2018b+.  This is the GEOMETRICALLY CORRECT choice.
%     2. the thin-shell realisation.offset_vertices_raw(x, y, t) — vertex-normal
%        offset with miter limit.  Same routine plot_steel_solve.m uses.
%        Robust on convex polygons; can produce small artefacts on sharp
%        concave corners but rarely catastrophic for hull-like shapes.
%
%   Returns Nx2 polygon (closed-implicitly, last vertex != first), or []
%   if the offset polygon collapsed.
    P_in = [];
    if length(x_outer) < 3 || ~isfinite(t) || t <= 0, return; end
    x_outer = x_outer(:);  y_outer = y_outer(:);

    % --- Method 1: polybuffer ---------------------------------------------
    try
        % polyshape silently warns on duplicate / collinear points; suppress.
        ws = warning('off', 'MATLAB:polyshape:repairedBySimplify');
        ps = polyshape(x_outer, y_outer);
        ps_in = polybuffer(ps, -t, 'JointType', 'miter', 'MiterLimit', 10);
        warning(ws);
        if ~isempty(ps_in.Vertices) && area(ps_in) > 1e-10
            % polybuffer may return multiple regions; keep the largest.
            R = regions(ps_in);
            if numel(R) > 1
                [~, kbig] = max(arrayfun(@area, R));
                ps_in = R(kbig);
            end
            V = ps_in.Vertices;
            % polyshape vertices may include NaN row separators; drop them.
            V = V(all(~isnan(V), 2), :);
            if size(V, 1) >= 3
                P_in = V;
                return;
            end
        end
    catch
        % polyshape/polybuffer not available or failed — fall through.
    end

    % --- Method 2: vertex-normal offset (steel's method) ------------------
    try
        [xi, yi] = mwecmass.internal.offset_polygon(x_outer, y_outer, t);
        if length(xi) >= 3
            A_outer = polyarea(x_outer, y_outer);
            A_inner = polyarea(xi, yi);
            if A_inner > 1e-10 && A_inner < A_outer
                P_in = [xi(:), yi(:)];
            end
        end
    catch
    end
end

function prof = build_profile_from_outer_contours(cstr)
% BUILD_PROFILE_FROM_OUTER_CONTOURS  Last-resort fallback when
%   mwecmass.output.figures.build_silhouette_profile is unavailable and
%   config.profile is empty.  Reconstructs a coarse XZ silhouette
%   from cstr.contours_outer by extracting (x_min, x_max) at each
%   z-sample.  Body frame.
    N = length(cstr.strip_z_lo);
    pts_R = [];  pts_L = [];
    for i = 1:N
        ci = cstr.contours_outer{i};
        if isempty(ci), continue; end
        n_z = length(ci);
        z_b = linspace(cstr.strip_z_lo(i), cstr.strip_z_hi(i), n_z)';
        for k = 1:n_z
            pk = ci{k};
            if ~isempty(pk) && size(pk, 1) >= 3
                pts_R(end+1, :) = [max(pk(:,1)), z_b(k)];   %#ok<AGROW>
                pts_L(end+1, :) = [min(pk(:,1)), z_b(k)];   %#ok<AGROW>
            end
        end
    end
    if isempty(pts_R)
        prof = [];
        return;
    end
    [~, ord_R] = sort(pts_R(:,2));
    [~, ord_L] = sort(pts_L(:,2), 'descend');
    prof = [pts_R(ord_R, :); pts_L(ord_L, :)];
end
