function fig = plot_steel_solve(config, steel_data)
%PLOT_STEEL_SOLVE 2D cross-section of realised steel-fill geometry with zones, outline, fill level.
% Renders solid-steel zone (z≤z_fill), jacket annulus, air cavity with 45° hatching, and
% outlines (hull, inner offset for thickness, waterline at draft, z_fill cut). Zone
% annotations (ρ values). NaN-resilient: when solver returns
% infeasible, uses optimiser frame and draws red banner. Uses polygon offset from
% mwecmass.internal.offset_polygon and silhouette from build_silhouette_profile.

    style = mwecmass.output.figures.presentation_style(config);

    %% ── Validation ────────────────────────────────────────────────────
    if isempty(steel_data) || ~isstruct(steel_data)
        error('plot_steel_solve:NoSteelData', 'steel_data is empty or not a struct.');
    end
    required = {'t_steel','z_fill','draft','vertical_shift','feasible', ...
                'M_steel','M_air','M_total','V_steel','V_air', ...
                'GM_realised','T_heave_realised','T_pitch_realised', ...
                'CG_z_world','CB_z_world','rho_steel','rho_air', ...
                'targets','residuals','mass_balance_error_pct', ...
                'vs_optimiser','draft_optimiser','t_min','t_min_active'};
    for k = 1:length(required)
        if ~isfield(steel_data, required{k})
            error('plot_steel_solve:MissingField', ...
                  'steel_data.%s missing.', required{k});
        end
    end
    if ~isfield(config, 'profile') || isempty(config.profile)
        error('plot_steel_solve:NoProfile', ...
              'config.profile missing — needed for hull silhouette.');
    end

    %% ── Plotting frame (NaN-safe) ─────────────────────────────────────
    if isfinite(steel_data.vertical_shift)
        vs_plot      = steel_data.vertical_shift;
        draft_plot   = steel_data.draft;
        is_real      = true;
    else
        vs_plot      = steel_data.vs_optimiser;
        draft_plot   = steel_data.draft_optimiser;
        is_real      = false;
    end

    %% ── Smooth silhouette in WORLD frame (waterline = 0) ─────────────
    % Preserve the builder's ordered profile; only the legacy raw fallback
    % needs plotting-only side ordering before offsetting and clipping.
    profile_is_raw = false;
    try
        prof_body = mwecmass.output.figures.build_silhouette_profile(config);
        if isequal(prof_body, config.profile)
            profile_is_raw = true;
        end
    catch
        prof_body = config.profile;
        profile_is_raw = true;
    end
    if profile_is_raw
        profile_error = 'mwecmass:figures:plot_steel_solve:ProfileInvalid';
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
        % Traverse right side bottom-to-top, then left side top-to-bottom.
        % Reordering and exact deduplication only; coordinate values are kept.
        prof_body = [centre(1, :); right; centre(2, :); flipud(left); centre(1, :)];
        edge_lengths = hypot(diff(prof_body(:, 1)), diff(prof_body(:, 2)));
        signed_area = 0.5 * sum(prof_body(1:end-1, 1) .* prof_body(2:end, 2) - ...
            prof_body(2:end, 1) .* prof_body(1:end-1, 2));
        if any(~isfinite(edge_lengths)) || any(edge_lengths <= 0) || ...
                ~isfinite(signed_area) || signed_area <= 0
            error(profile_error, 'The ordered fallback hull boundary is invalid.');
        end
    end
    profile_world = prof_body + [0, vs_plot];
    px_outer = profile_world(:, 1);
    pz_outer = profile_world(:, 2);
    z_fill_world = steel_data.z_fill + vs_plot;

    %% ── Inner jacket offset (the THICKNESS visualisation) ────────────
    if isfinite(steel_data.t_steel) && steel_data.t_steel > 1e-6
        [px_inner, pz_inner] = mwecmass.internal.offset_polygon( ...
                                   px_outer, pz_outer, steel_data.t_steel);
    else
        px_inner = [];  pz_inner = [];
    end
    has_inner = length(px_inner) >= 3;

    %% ── Sutherland–Hodgman clips at z_fill (no polyshape) ────────────
    [x_below_o, z_below_o] = clip_polygon_below_z(px_outer, pz_outer, z_fill_world);
    [x_above_o, z_above_o] = clip_polygon_above_z(px_outer, pz_outer, z_fill_world);
    if has_inner
        [x_above_i, z_above_i] = clip_polygon_above_z(px_inner, pz_inner, z_fill_world);
    else
        x_above_i = [];  z_above_i = [];
    end

    %% ── Figure ───────────────────────────────────────────────────────
    fig = mwecmass.output.figures.new_figure(style, 'tall_double_column');
    % This single-panel diagnostic needs additional vertical plot-box room so
    % axis equal can display the full hull at a useful physical width beneath
    % the 2x2 legend. Keep the shared preset unchanged for other figures.
    set(fig, 'Units', 'centimeters');
    fig.Position(3:4) = [17, 15];
    set(fig, 'Name', 'Steel-Fill Solver');

    %% ============= 2D cross-section =====================
    ax1 = axes(fig);
    hold(ax1, 'on'); axis(ax1, 'equal');

    % ─ PASS 1: filled regions ──────────────────────────────────────
    if length(x_below_o) >= 3
        patch(ax1, x_below_o, z_below_o, style.fill_palette.solid_material, ...
              'EdgeColor', 'none', 'FaceAlpha', 1.0, ...
              'HandleVisibility', 'off');
    end
    if length(x_above_o) >= 3
        patch(ax1, x_above_o, z_above_o, style.fill_palette.jacket_material, ...
              'EdgeColor', 'none', 'FaceAlpha', 1.0, ...
              'HandleVisibility', 'off');
    end
    if length(x_above_i) >= 3
        patch(ax1, x_above_i, z_above_i, style.fill_palette.void, ...
              'EdgeColor', 'none', 'FaceAlpha', 1.0, ...
              'HandleVisibility', 'off');

        % ─ PASS 2: 45° hatch over the air cavity ────────────────
        try
            draw_hatch(ax1, x_above_i, z_above_i, style.fill_palette.hatch, style.hatch_spacing, ...
                       style.line_width.hatch);
        catch ME
            warning('plot_steel_solve:HatchFailed', ...
                    'Hatch rendering failed (non-fatal): %s', ME.message);
        end
    end

    % ─ PASS 3: outlines ────────────────────────────────────────────
    h_hull = plot(ax1, [px_outer; px_outer(1)], [pz_outer; pz_outer(1)], ...
                  '-', 'Color', style.fill_palette.boundary, ...
                  'HandleVisibility', 'off');
    mwecmass.output.figures.style_line(h_hull, style, 'boundary');

    h_inner = [];
    if has_inner
        h_inner = plot(ax1, [px_inner; px_inner(1)], [pz_inner; pz_inner(1)], ...
                       '--', 'Color', style.fill_palette.inner_boundary, ...
                       'HandleVisibility', 'off');
        mwecmass.output.figures.style_line(h_inner, style, 'boundary');
    end

    % Keep the waterline and x-limits close to the realised hull so the
    % equal-aspect axes use the available width beneath the legend.
    hull_x_limits = [min(px_outer), max(px_outer)];
    hull_x_margin = 0.05 * max(diff(hull_x_limits), eps);
    x_range = hull_x_limits + [-hull_x_margin, hull_x_margin];
    h_wl    = plot(ax1, x_range, [0, 0], '--', ...
                   'Color', style.fill_palette.waterline);
    mwecmass.output.figures.style_line(h_wl, style, 'reference');

    latex_slash = char(92);

    % Infeasibility banner
    if ~is_real
        yl = ylim(ax1);
        h_infeas = text(ax1, mean(x_range), yl(2) - 0.05*(yl(2)-yl(1)), ...
             'INFEASIBLE — geometry rendered with optimiser draft', ...
             'HorizontalAlignment', 'center', 'Color', style.status_palette.bad, ...
             'FontWeight', 'bold');
        mwecmass.output.figures.style_text(h_infeas, style, 'annotation');
    end

    % ─ PASS 4: legend (build conditionally) ─────────────────────────
    h_solid_proxy  = patch(ax1, NaN, NaN, style.fill_palette.solid_material,  'EdgeColor', 'none');
    h_jacket_proxy = patch(ax1, NaN, NaN, style.fill_palette.jacket_material, 'EdgeColor', 'none');
    h_air_proxy    = patch(ax1, NaN, NaN, style.fill_palette.void, ...
                           'EdgeColor', style.fill_palette.hatch, 'LineStyle', '--');

    rho_material_txt = num2str(steel_data.rho_fill, '%.0f');
    shell_thickness_txt = num2str(steel_data.t_steel, '%.4f');
    air_density_txt = num2str(steel_data.rho_air, '%.1f');
    leg_h = [h_solid_proxy, h_jacket_proxy, h_air_proxy, h_wl];
    leg_l = {['Solid Ballast, $' latex_slash 'rho_{' latex_slash ...
              'mathrm{material}} = ' rho_material_txt '$ [kg/m$^3$]'], ...
             ['Thin Shell, $t_{plate} = ' shell_thickness_txt '$ [m]'], ...
             ['Hollow Volume, $' latex_slash 'rho_{' latex_slash ...
              'mathrm{Air}} = ' air_density_txt '$ [kg/m$^3$]'], ...
             'Waterline, $Z = 0$ m'};

    lg = legend(ax1, leg_h, leg_l, 'Location', 'northoutside', ...
        'NumColumns', 2, 'AutoUpdate', 'off');
    mwecmass.output.figures.style_legend(lg, style);
    lg.FontSize = min(style.font_size.legend, 8.5);

    xlabel(ax1, {'X [m]', ...
        ['Stage 3: Thin Shell Construction, $' latex_slash 'rho_{' latex_slash ...
         'mathrm{material}} = ' ...
         rho_material_txt '$ [kg/m$^3$]']});
    ylabel(ax1, 'Z [m]');
    mwecmass.output.figures.apply_axes_style(ax1, style);
    set(ax1, 'Box', 'on', 'XGrid', 'off', 'YGrid', 'off', ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off');
    hull_z_limits = [min([pz_outer; 0]), max([pz_outer; 0])];
    hull_z_margin = 0.05 * max(diff(hull_z_limits), eps);
    xlim(ax1, hull_x_limits + [-1.1, 1.1] * hull_x_margin);
    zlim(ax1, hull_z_limits + [-1.1, 1.1] * hull_z_margin);
    hold(ax1, 'off');

    %% ── Save ────────────────────────────────────────────────────────
    save_figure(fig, config);
end


% =======================================================================
%  POLYGON CLIPPING (Sutherland–Hodgman against horizontal half-plane)
% =======================================================================

function [xc, zc] = clip_polygon_below_z(x, z, z_cut)
% Returns the polygon clipped to z <= z_cut.
    [xc, zc] = mwecmass.internal.clip_z(x, z, z_cut, 'below');
end

function [xc, zc] = clip_polygon_above_z(x, z, z_cut)
% Returns the polygon clipped to z >= z_cut.
    [xc, zc] = mwecmass.internal.clip_z(x, z, z_cut, 'above');
end


% =======================================================================
%  HATCH (45° diagonal lines clipped to a polygon)
% =======================================================================

function draw_hatch(ax, x_poly, z_poly, color, spacing, line_width)
% Sweep lines x − z = c and clip to polygon via inpolygon.
% line_width [pt] from style.line_width.hatch (role-sourced, not hardcoded).
    if length(x_poly) < 3, return; end

    x_lo = min(x_poly);  x_hi = max(x_poly);
    z_lo = min(z_poly);  z_hi = max(z_poly);
    if (x_hi - x_lo) < 1e-9 || (z_hi - z_lo) < 1e-9, return; end

    c_min = x_lo - z_hi;
    c_max = x_hi - z_lo;
    c_vals = c_min : spacing : c_max;

    for ci = 1:length(c_vals)
        c = c_vals(ci);
        z_line = linspace(z_lo, z_hi, 200)';
        x_line = z_line + c;

        in = inpolygon(x_line, z_line, x_poly, z_poly);
        if ~any(in), continue; end

        d = diff(in);
        starts = find(d == 1) + 1;
        stops  = find(d == -1);
        if in(1),   starts = [1; starts];          end %#ok<AGROW> -- starts is rebuilt from find(d==1) fresh each ci iteration; this prepends at most one element, not an accumulating loop.
        if in(end), stops  = [stops; length(in)];  end %#ok<AGROW> -- stops is rebuilt from find(d==-1) fresh each ci iteration; this appends at most one element, not an accumulating loop.
        n_seg = min(length(starts), length(stops));
        for s = 1:n_seg
            plot(ax, x_line(starts(s):stops(s)), ...
                     z_line(starts(s):stops(s)), '-', ...
                 'Color', color, 'LineWidth', line_width, ...
                 'HandleVisibility', 'off');
        end
    end
end


% =======================================================================
%  SAVE FIGURE
% =======================================================================

function save_figure(fig, config)
    out_dir = mwecmass.output.output_dir('thin_shell');
    try
        [png_files, ~] = mwecmass.output.figures.export_figure( ...
            fig, 'Steel_Solve', config, out_dir);
        fprintf('      Steel-fill plot saved: %s\n', png_files{1});
    catch ME
        warning('plot_steel_solve:SaveFailed', ...
                'Could not save figure: %s', ME.message);
    end
end
