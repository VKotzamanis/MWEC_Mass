function fig = plot_steel_solve(config, steel_data)
% PLOT_STEEL_SOLVE  2D diagnostic for the steel-fill solver.
%
%   fig = plot_steel_solve(config, steel_data)
%
%   One figure, two panels (matches WEC_Visualization.visualize_2d_design):
%
%     LEFT  — 2D midplane cross-section.  Solid-steel zone (z<=z_fill)
%             filled in dark grey.  Steel jacket annulus above z_fill
%             filled in lighter grey.  Air cavity inside the jacket
%             rendered white with 45° diagonal hatching.  Inner jacket
%             offset is overlaid as a dashed dark-grey line — this is
%             where the thickness is visible.  Hull outline (solid black),
%             dashed-blue waterline at the re-solved draft, dash-dot
%             orange z_fill cut, CG/CB markers, and a "SOLID STEEL" /
%             "AIR" label centred in each zone.
%
%     RIGHT — Numeric card: solution (t_steel, z_fill, draft), mass
%             breakdown (V_steel/V_air/M_steel/M_air/M_total + balance
%             error), and target-vs-realised table for GM and natural
%             periods, color-coded (green ≤ 5%, amber ≤ 15%, red > 15%).
%
%   NaN-resilient.  When the solver returned infeasible (vertical_shift
%   = NaN), the geometry is rendered using steel_data.vs_optimiser as
%   the plotting frame and a red banner is drawn on the figure.
%
%   See also: WEC_Shell_Offset.solve, WEC_Visualization.visualize_2d_design

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
    % Prefer WEC_Visualization.build_smooth_viz_profile (handles non-convex
    % hulls cleanly).  Fall back to raw config.profile if access-restricted.
    try
        prof_body = WEC_Visualization.build_smooth_viz_profile(config);
    catch
        prof_body = config.profile;
    end
    profile_world = prof_body + [0, vs_plot];
    px_outer = profile_world(:, 1);
    pz_outer = profile_world(:, 2);
    z_fill_world = steel_data.z_fill + vs_plot;

    %% ── Inner jacket offset (the THICKNESS visualisation) ────────────
    if isfinite(steel_data.t_steel) && steel_data.t_steel > 1e-6
        [px_inner, pz_inner] = WEC_Shell_Offset.offset_vertices_raw( ...
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

    %% ── Style constants (match WEC_Constructable_Hull / Visualization) ─
    STEEL_DARK    = [0.45, 0.46, 0.50];   % solid steel below z_fill
    STEEL_LIGHT   = [0.74, 0.76, 0.80];   % jacket annulus above z_fill
    AIR_WHITE     = [1.00, 1.00, 1.00];
    HATCH_COLOR   = [0.50, 0.52, 0.58];
    HATCH_SPACING = 0.06;
    BOUND_COLOR   = [0.08, 0.08, 0.08];
    INNER_COLOR   = [0.30, 0.30, 0.32];
    WL_COLOR      = [0.15, 0.55, 0.95];
    ZFILL_COLOR   = [0.95, 0.55, 0.10];
    OK_GREEN      = [0.10, 0.55, 0.10];
    WARN_AMBER    = [0.80, 0.55, 0.10]; %#ok<NASGU>
    BAD_RED       = [0.80, 0.15, 0.15];
    FN            = 'Times New Roman';
    FN_MONO       = 'Consolas';

    %% ── Figure ───────────────────────────────────────────────────────
    fig = figure('Name', 'Steel-Fill Solver', 'Color', 'w', ...
                 'Position', [80, 80, 1500, 760]);

    %% ============= LEFT PANEL: 2D cross-section =====================
    ax1 = subplot(1, 2, 1);
    hold(ax1, 'on'); axis(ax1, 'equal');

    % ─ PASS 1: filled regions ──────────────────────────────────────
    if length(x_below_o) >= 3
        patch(ax1, x_below_o, z_below_o, STEEL_DARK, ...
              'EdgeColor', 'none', 'FaceAlpha', 1.0, ...
              'HandleVisibility', 'off');
    end
    if length(x_above_o) >= 3
        patch(ax1, x_above_o, z_above_o, STEEL_LIGHT, ...
              'EdgeColor', 'none', 'FaceAlpha', 1.0, ...
              'HandleVisibility', 'off');
    end
    if length(x_above_i) >= 3
        patch(ax1, x_above_i, z_above_i, AIR_WHITE, ...
              'EdgeColor', 'none', 'FaceAlpha', 1.0, ...
              'HandleVisibility', 'off');

        % ─ PASS 2: 45° hatch over the air cavity ────────────────
        try
            draw_hatch(ax1, x_above_i, z_above_i, HATCH_COLOR, HATCH_SPACING);
        catch ME
            warning('plot_steel_solve:HatchFailed', ...
                    'Hatch rendering failed (non-fatal): %s', ME.message);
        end
    end

    % ─ PASS 3: outlines ────────────────────────────────────────────
    h_hull = plot(ax1, [px_outer; px_outer(1)], [pz_outer; pz_outer(1)], ...
                  '-', 'Color', BOUND_COLOR, 'LineWidth', 1.8);

    h_inner = [];
    if has_inner
        h_inner = plot(ax1, [px_inner; px_inner(1)], [pz_inner; pz_inner(1)], ...
                       '--', 'Color', INNER_COLOR, 'LineWidth', 1.3);
    end

    x_range = [min(px_outer) - 0.30, max(px_outer) + 0.30];
    h_zfill = plot(ax1, x_range, [z_fill_world, z_fill_world], '-.', ...
                   'Color', ZFILL_COLOR, 'LineWidth', 2.0);
    h_wl    = plot(ax1, x_range, [0, 0], '--', ...
                   'Color', WL_COLOR, 'LineWidth', 2.0);

    h_cg = [];  h_cb = [];
    if isfinite(steel_data.CG_z_world)
        h_cg = plot(ax1, 0, steel_data.CG_z_world, 'ro', ...
                    'MarkerSize', 12, 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
    end
    if isfinite(steel_data.CB_z_world)
        h_cb = plot(ax1, 0, steel_data.CB_z_world, 'bs', ...
                    'MarkerSize', 12, 'MarkerFaceColor', 'b', 'LineWidth', 1.5);
    end

    % ─ PASS 4: zone annotations (centered text in each region) ──────
    if length(x_below_o) >= 3
        z_anno = 0.5 * (max(z_below_o) + min(z_below_o));
        text(ax1, 0, z_anno, ...
             sprintf('SOLID STEEL\n\\rho = %.0f kg/m^3', steel_data.rho_steel), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'Color', 'w', 'FontWeight', 'bold', 'FontSize', 11, ...
             'FontName', FN, 'Interpreter', 'tex');
    end
    if length(x_above_i) >= 3
        z_anno_air = 0.5 * (max(z_above_i) + min(z_above_i));
        text(ax1, 0, z_anno_air, ...
             sprintf('AIR\n\\rho = %.0f kg/m^3', steel_data.rho_air), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
             'Color', [0.30 0.30 0.30], 'FontWeight', 'bold', 'FontSize', 10, ...
             'FontName', FN, 'Interpreter', 'tex', ...
             'BackgroundColor', AIR_WHITE, 'Margin', 3);
    end

    % t_steel callout at the top jacket edge — tiny arrow + "t = X mm"
    if has_inner && length(x_above_o) >= 3
        z_top  = max(z_above_o) - 0.15 * (max(z_above_o) - min(z_above_o));
        x_jacket_at_top_outer = max(px_outer);
        x_jacket_at_top_inner = max(px_inner);
        if x_jacket_at_top_outer > x_jacket_at_top_inner
            x_mid_jacket = 0.5 * (x_jacket_at_top_outer + x_jacket_at_top_inner);
            text(ax1, x_mid_jacket + 0.08, z_top, ...
                 sprintf('t_{steel} = %.1f mm', steel_data.t_steel * 1000), ...
                 'FontSize', 10, 'FontName', FN, 'Color', INNER_COLOR, ...
                 'Rotation', 90, 'VerticalAlignment', 'bottom', ...
                 'Interpreter', 'tex');
        end
    end

    % Infeasibility banner
    if ~is_real
        yl = ylim(ax1);
        text(ax1, mean(x_range), yl(2) - 0.05*(yl(2)-yl(1)), ...
             'INFEASIBLE — geometry rendered with optimiser draft', ...
             'HorizontalAlignment', 'center', 'Color', BAD_RED, ...
             'FontWeight', 'bold', 'FontSize', 12, 'FontName', FN);
    end

    % ─ PASS 5: legend (build conditionally) ─────────────────────────
    h_solid_proxy  = patch(ax1, NaN, NaN, STEEL_DARK,  'EdgeColor', 'none');
    h_jacket_proxy = patch(ax1, NaN, NaN, STEEL_LIGHT, 'EdgeColor', 'none');
    h_air_proxy    = patch(ax1, NaN, NaN, AIR_WHITE, ...
                           'EdgeColor', HATCH_COLOR, 'LineStyle', '--');

    leg_h = [h_solid_proxy, h_jacket_proxy, h_air_proxy, h_hull];
    if steel_data.t_min_active
        jacket_label = sprintf('Steel jacket  (t = %.4f m  *** at min %.4f m ***)', ...
                               steel_data.t_steel, steel_data.t_min);
    else
        jacket_label = sprintf('Steel jacket  (t = %.4f m,  min = %.4f m)', ...
                               steel_data.t_steel, steel_data.t_min);
    end
    leg_l = {'Solid steel  (z < z_{fill})', ...
             jacket_label, ...
             'Air cavity', ...
             'Hull outer surface'};
    if ~isempty(h_inner)
        leg_h(end+1) = h_inner;
        leg_l{end+1} = 'Jacket inner offset (dashed)';
    end
    leg_h(end+1) = h_zfill;
    leg_l{end+1} = sprintf('z_{fill} = %.3f m', steel_data.z_fill);
    leg_h(end+1) = h_wl;
    leg_l{end+1} = sprintf('Waterline (draft = %.3f m)', draft_plot);
    if ~isempty(h_cg), leg_h(end+1) = h_cg; leg_l{end+1} = 'CG'; end
    if ~isempty(h_cb), leg_h(end+1) = h_cb; leg_l{end+1} = 'CB'; end

    legend(ax1, leg_h, leg_l, ...
           'Location', 'eastoutside', 'FontName', FN, ...
           'FontSize', 10, 'Interpreter', 'tex');

    xlabel(ax1, '$x$ [m]', 'Interpreter', 'latex', 'FontSize', 13, 'FontName', FN);
    ylabel(ax1, '$z$ [m]   (world frame, waterline = 0)', ...
           'Interpreter', 'latex', 'FontSize', 13, 'FontName', FN);
    title(ax1, 'Steel-Fill Realised Geometry', ...
          'FontSize', 14, 'FontName', FN);
    grid(ax1, 'on');
    set(ax1, 'GridLineStyle', ':', 'GridAlpha', 0.3, 'Box', 'on', ...
             'FontName', FN, 'FontSize', 11, ...
             'TickLabelInterpreter', 'latex');
    hold(ax1, 'off');

    %% ============= RIGHT PANEL: properties card =====================
    ax2 = subplot(1, 2, 2);
    axis(ax2, 'off');

    if steel_data.feasible
        status_str = 'FEASIBLE';
        status_col = OK_GREEN;
    else
        status_str = 'INFEASIBLE — closest mass-balanced point';
        status_col = BAD_RED;
    end

    % Build the text block.  Use a monospace font so columns align.
    txt = {};
    txt{end+1} = 'STEEL-FILL SOLUTION';
    txt{end+1} = '====================';
    if steel_data.t_min_active
        txt{end+1} = sprintf('  t_steel    = %.5f m  (%.2f in)  *** AT t_min BOUND ***', ...
                              steel_data.t_steel, steel_data.t_steel/0.0254);
    else
        txt{end+1} = sprintf('  t_steel    = %.5f m  (%.2f in)', ...
                              steel_data.t_steel, steel_data.t_steel/0.0254);
    end
    txt{end+1} = sprintf('  t_min      = %.5f m  (%.2f in, fabrication floor)', ...
                          steel_data.t_min, steel_data.t_min/0.0254);
    txt{end+1} = sprintf('  z_fill     = %.4f m  (body frame)', steel_data.z_fill);
    txt{end+1} = sprintf('  draft      = %s m', fmt_or_dash(steel_data.draft, '%.4f'));
    if isfield(steel_data, 't_max')
        txt{end+1} = sprintf('  t_max      = %.4f m  (search bound)', steel_data.t_max);
    end
    txt{end+1} = '';
    txt{end+1} = 'MASS BREAKDOWN';
    txt{end+1} = '====================';
    txt{end+1} = sprintf('  V_steel    = %.4f m^3', steel_data.V_steel);
    txt{end+1} = sprintf('  V_air      = %.4f m^3', steel_data.V_air);
    txt{end+1} = sprintf('  M_steel    = %.1f kg  (rho %.0f)', ...
                          steel_data.M_steel, steel_data.rho_steel);
    txt{end+1} = sprintf('  M_air      = %.1f kg  (rho %.0f)', ...
                          steel_data.M_air, steel_data.rho_air);
    txt{end+1} = sprintf('  M_total    = %.1f kg', steel_data.M_total);
    txt{end+1} = sprintf('  M_target   = %.1f kg  (err %s)', ...
                          steel_data.targets.mass, ...
                          fmt_pct(steel_data.residuals.dmass_pct));
    txt{end+1} = sprintf('  mass-bal   = %.4f%%  (constraint, hard)', ...
                          steel_data.mass_balance_error_pct);
    txt{end+1} = '';
    txt{end+1} = 'TARGETS vs REALISED';
    txt{end+1} = '====================';
    txt{end+1} = sprintf('  %-10s %10s %10s %10s', ...
                          'Quantity', 'Target', 'Realised', 'Residual');
    txt{end+1} = sprintf('  %-10s %10.4f %10s %10s', ...
                          'GM (m)', steel_data.targets.GM, ...
                          fmt_or_dash(steel_data.GM_realised, '%.4f'), ...
                          fmt_pct(steel_data.residuals.dGM_pct));
    txt{end+1} = sprintf('  %-10s %10.3f %10s %10s', ...
                          'T_heave(s)', steel_data.targets.T_heave, ...
                          fmt_or_dash(steel_data.T_heave_realised, '%.3f'), ...
                          fmt_pct(steel_data.residuals.dT_heave_pct));
    txt{end+1} = sprintf('  %-10s %10.3f %10s %10s', ...
                          'T_pitch(s)', steel_data.targets.T_pitch, ...
                          fmt_or_dash(steel_data.T_pitch_realised, '%.3f'), ...
                          fmt_pct(steel_data.residuals.dT_pitch_pct));

    text(ax2, 0.02, 0.97, txt, 'Units', 'normalized', ...
         'VerticalAlignment', 'top', 'FontName', FN_MONO, ...
         'FontSize', 11, 'Interpreter', 'none');

    % Status line — bold, colour-coded
    text(ax2, 0.02, 0.13, sprintf('STATUS: %s', status_str), ...
         'Units', 'normalized', 'FontName', FN, 'FontSize', 13, ...
         'FontWeight', 'bold', 'Color', status_col, 'Interpreter', 'none');

    % Bound-hit alert line in red, just below the STATUS line
    if steel_data.t_min_active
        text(ax2, 0.02, 0.07, ...
             sprintf('t_steel AT FABRICATION FLOOR (t_min = %.4f m / %.2f in)', ...
                     steel_data.t_min, steel_data.t_min/0.0254), ...
             'Units', 'normalized', 'FontName', FN, 'FontSize', 12, ...
             'FontWeight', 'bold', 'Color', BAD_RED, 'Interpreter', 'none');
    end

    if ~is_real
        text(ax2, 0.02, 0.03, ...
             '(geometry rendered with optimiser draft)', ...
             'Units', 'normalized', 'FontName', FN, 'FontSize', 10, ...
             'FontAngle', 'italic', 'Color', [0.45 0.45 0.45], ...
             'Interpreter', 'none');
    end

    %% ── Super-title ─────────────────────────────────────────────────
    sgtitle(sprintf(['Steel-Fill Solver Result   |   ' ...
                     't_{steel} = %.4f m   z_{fill} = %.3f m   ' ...
                     'draft = %.3f m   M_{total} = %.0f kg'], ...
                    steel_data.t_steel, steel_data.z_fill, draft_plot, ...
                    steel_data.M_total), ...
            'FontName', FN, 'FontSize', 13, 'FontWeight', 'bold');

    %% ── Save ────────────────────────────────────────────────────────
    save_figure(fig);
end


% =======================================================================
%  POLYGON CLIPPING (Sutherland–Hodgman against horizontal half-plane)
% =======================================================================

function [xc, zc] = clip_polygon_below_z(x, z, z_cut)
% Returns the polygon clipped to z <= z_cut.
    [xc, zc] = clip_halfspace(x, z, z_cut, true);
end

function [xc, zc] = clip_polygon_above_z(x, z, z_cut)
% Returns the polygon clipped to z >= z_cut.
    [xc, zc] = clip_halfspace(x, z, z_cut, false);
end

function [xc, zc] = clip_halfspace(x, z, z_cut, keep_below)
    n = length(x);
    if n < 3
        xc = [];  zc = [];
        return;
    end
    xc = zeros(2*n, 1);
    zc = zeros(2*n, 1);
    cnt = 0;

    for i = 1:n
        j = mod(i, n) + 1;
        zi = z(i);  zj = z(j);
        if keep_below
            in_i = (zi <= z_cut);
            in_j = (zj <= z_cut);
        else
            in_i = (zi >= z_cut);
            in_j = (zj >= z_cut);
        end

        if in_i
            cnt = cnt + 1;
            xc(cnt) = x(i);
            zc(cnt) = zi;
        end
        if in_i ~= in_j && abs(zj - zi) > 1e-14
            t  = (z_cut - zi) / (zj - zi);
            cnt = cnt + 1;
            xc(cnt) = x(i) + t * (x(j) - x(i));
            zc(cnt) = z_cut;
        end
    end

    xc = xc(1:cnt);
    zc = zc(1:cnt);
end


% =======================================================================
%  HATCH (45° diagonal lines clipped to a polygon)
% =======================================================================

function draw_hatch(ax, x_poly, z_poly, color, spacing)
% Reuses Constructable_Hull's algorithm: sweep lines x − z = c and clip
% to the polygon via inpolygon.
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
        if in(1),   starts = [1; starts];          end
        if in(end), stops  = [stops; length(in)];  end
        n_seg = min(length(starts), length(stops));
        for s = 1:n_seg
            plot(ax, x_line(starts(s):stops(s)), ...
                     z_line(starts(s):stops(s)), '-', ...
                 'Color', color, 'LineWidth', 0.4, ...
                 'HandleVisibility', 'off');
        end
    end
end


% =======================================================================
%  STRING FORMATTING
% =======================================================================

function s = fmt_or_dash(v, fmt)
    if isfinite(v)
        s = sprintf(fmt, v);
    else
        s = '   —   ';
    end
end

function s = fmt_pct(v)
    if isfinite(v)
        s = sprintf('%+.2f%%', v);
    else
        s = '   —   ';
    end
end


% =======================================================================
%  SAVE FIGURE
% =======================================================================

function save_figure(fig)
    out_dir = fullfile(fileparts(mfilename('fullpath')), 'Plots');
    if ~exist(out_dir, 'dir')
        [ok, msg] = mkdir(out_dir);
        if ~ok
            warning('plot_steel_solve:MkdirFailed', ...
                    'Could not create Plots directory: %s', msg);
            out_dir = fileparts(mfilename('fullpath'));
        end
    end
    png_path = fullfile(out_dir, 'Steel_Solve.png');
    fig_path = fullfile(out_dir, 'Steel_Solve.fig');
    try
        saveas(fig, png_path);
        savefig(fig, fig_path);
        fprintf('      Steel-fill plot saved: %s\n', png_path);
    catch ME
        warning('plot_steel_solve:SaveFailed', ...
                'Could not save figure: %s', ME.message);
    end
end
