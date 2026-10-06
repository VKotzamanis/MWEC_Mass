function plot_draft_landscape(sweep, config)
%PLOT_DRAFT_LANDSCAPE Stage-1 objective and hydrostatic convergence report.
% The sweep variable is vertical_shift; physical draft is read from each
% evaluated property record (vertical shift is not used as the plot variable).
% Tier-1 properties_2d records do not provide draft, while refined
% properties_3d records may provide p.draft.  The fallback below matches
% optim/run.m and properties_3d.m: draft = abs(hull_z_min + vertical_shift).

    style = mwecmass.output.figures.presentation_style(config);

    vs = sweep.vs(:)';
    feasible = logical(sweep.feasible(:)');
    fval = sweep.fval(:)';
    N = numel(vs);
    if numel(feasible) ~= N || numel(fval) ~= N || numel(sweep.props) ~= N
        error('WEC:DraftLandscapeSizeMismatch', ...
            'sweep.vs, sweep.feasible, sweep.fval, and sweep.props must have equal lengths.');
    end

    draft = nan(1, N);
    T_h = nan(1, N);
    T_p = nan(1, N);
    GM = nan(1, N);
    mass_err = nan(1, N);
    for k = 1:N
        if ~isempty(sweep.props{k})
            p = sweep.props{k};
            if isfield(p, 'draft')
                draft(k) = p.draft;
            end
            if isfield(p, 'periods')
                T_h(k) = p.periods.heave;
                T_p(k) = p.periods.pitch;
            end
            if isfield(p, 'GM_L')
                GM(k) = p.GM_L;
            end
            if isfield(p, 'mass_total') && isfield(p, 'mass_buoyant_force')
                mass_err(k) = abs(p.mass_total - p.mass_buoyant_force) ...
                    / max(p.mass_total, 1) * 100;
            end
        end
    end

    missing_draft = ~isfinite(draft);
    if any(missing_draft)
        if ~isfield(config, 'hull_z_min') || ~isscalar(config.hull_z_min) || ...
                ~isfinite(config.hull_z_min)
            error('WEC:DraftLandscapeMissingDraft', ...
                'Physical draft is missing and config.hull_z_min is unavailable.');
        end
        draft(missing_draft) = abs(config.hull_z_min + vs(missing_draft));
    end

    fig = mwecmass.output.figures.new_figure(style, 'tall_double_column');
    fig.Name = 'Stage 1 Convergence Report';
    tiled = tiledlayout(fig, 2, 2);
    draft_limits = local_limits(draft);

    % Objective: logarithmic X axis; invalid and nonpositive values are omitted.
    ax1 = nexttile(tiled);
    hold(ax1, 'on');
    valid_obj = isfinite(fval) & fval > 0 & isfinite(draft);
    h_feas = gobjects(0);
    h_infeas = gobjects(0);
    h_best = gobjects(0);
    if any(valid_obj & feasible)
        h_feas = plot(ax1, fval(valid_obj & feasible), draft(valid_obj & feasible), ...
            'o', 'Color', style.status_palette.ok);
    end
    if any(valid_obj & ~feasible)
        h_infeas = plot(ax1, fval(valid_obj & ~feasible), draft(valid_obj & ~feasible), ...
            'x', 'Color', style.status_palette.bad);
    end
    feasible_indices = find(valid_obj & feasible);
    if ~isempty(feasible_indices)
        [~, best_local] = min(fval(feasible_indices));
        best = feasible_indices(best_local);
        h_best = plot(ax1, fval(best), draft(best), 'p', ...
            'MarkerFaceColor', style.status_palette.ok, ...
            'MarkerEdgeColor', style.fill_palette.boundary);
    end
    mwecmass.output.figures.style_line([h_feas h_infeas h_best], style, 'curve');
    set(ax1, 'XScale', 'log', 'YLim', draft_limits);
    xlabel(ax1, {'Stage 1 Objective, $J(x^{2d})$ [-]', 'Convergence Feasibility'});
    ylabel(ax1, 'WEC Draft, $d$ [m]');
    local_legend(ax1, [h_feas h_infeas h_best], ...
        {'Feasible', 'Infeasible', 'Best'}, style);
    local_finalize_axes(ax1, style);

    % Natural periods: period is X and physical draft is Y.
    ax2 = nexttile(tiled);
    hold(ax2, 'on');
    period_limits = local_limits([T_h T_p config.T_heave_goal config.T_pitch_goal ...
        config.T_heave_range(:)' config.T_pitch_range(:)']);
    patch(ax2, [config.T_heave_range(1) config.T_heave_range(2) ...
        config.T_heave_range(2) config.T_heave_range(1)], ...
        [draft_limits(1) draft_limits(1) draft_limits(2) draft_limits(2)], ...
        style.color.series_a, 'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    patch(ax2, [config.T_pitch_range(1) config.T_pitch_range(2) ...
        config.T_pitch_range(2) config.T_pitch_range(1)], ...
        [draft_limits(1) draft_limits(1) draft_limits(2) draft_limits(2)], ...
        style.color.series_b, 'FaceAlpha', 0.08, 'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
    h_th = plot(ax2, T_h, draft, '-o', 'Color', style.color.series_a);
    h_tp = plot(ax2, T_p, draft, '-s', 'Color', style.color.series_b);
    h_th_goal = xline(ax2, config.T_heave_goal, '--', 'Color', style.color.series_a);
    h_tp_goal = xline(ax2, config.T_pitch_goal, '--', 'Color', style.color.series_b);
    mwecmass.output.figures.style_line([h_th h_tp], style, 'curve');
    mwecmass.output.figures.style_line([h_th_goal h_tp_goal], style, 'reference');
    set(ax2, 'XLim', period_limits, 'YLim', draft_limits);
    xlabel(ax2, {'Natural Periods, $T_{n}$ [s]', ...
        'Natural Periods, $T_{n}$, VS Design Periods, $T_{d}$'});
    ylabel(ax2, 'WEC Draft, $d$ [m]');
    local_legend(ax2, [h_th h_tp h_th_goal h_tp_goal], ...
        {'$T^{\mathrm{Heave}}_{n}$', '$T^{\mathrm{Pitch}}_{n}$', ...
         '$T^{\mathrm{Heave}}_{d}$', '$T^{\mathrm{Pitch}}_{d}$'}, style, 2);
    local_finalize_axes(ax2, style);

    % Pitch stability: retain the minimum-GM threshold and remove the old fill.
    ax3 = nexttile(tiled);
    hold(ax3, 'on');
    h_gm = plot(ax3, GM, draft, '-o', 'Color', style.fill_palette.boundary);
    h_gm_min = xline(ax3, config.gm_min, '--', 'Color', style.status_palette.bad);
    mwecmass.output.figures.style_line(h_gm, style, 'curve');
    mwecmass.output.figures.style_line(h_gm_min, style, 'reference');
    set(ax3, 'YLim', draft_limits);
    xlabel(ax3, {'Metacentric Height, $GM$ [m]', 'Hydrostatic Stability, Pitch'});
    ylabel(ax3, 'WEC Draft, $d$ [m]');
    local_legend(ax3, [h_gm h_gm_min], {'$GM$', '$GM_{\min}=0.2\,\mathrm{m}$'}, style);
    local_finalize_axes(ax3, style);

    % Heave stability: logarithmic mass-balance error, without a threshold callout.
    ax4 = nexttile(tiled);
    hold(ax4, 'on');
    mass_plot = mass_err;
    mass_plot(~isfinite(mass_plot) | mass_plot <= 0) = NaN;
    h_mass = plot(ax4, mass_plot, draft, '-o', 'Color', style.fill_palette.boundary);
    mwecmass.output.figures.style_line(h_mass, style, 'curve');
    set(ax4, 'XScale', 'log', 'YLim', draft_limits);
    xlabel(ax4, {'Mass Balance Error [\%]', 'Hydrostatic Stability, Heave'});
    ylabel(ax4, 'WEC Draft, $d$ [m]');
    local_finalize_axes(ax4, style);

    title(tiled, 'Stage 1 Convergence Report');
    mwecmass.output.figures.apply_layout_style(tiled, style);
    saved = mwecmass.output.figures.export_figure(fig, 'WEC_DraftLandscape', config);
    fprintf('  Draft landscape saved: %s\n', saved{1});
end

function limits = local_limits(values)
    values = values(isfinite(values));
    if isempty(values)
        limits = [0 1];
        return;
    end
    lo = min(values);
    hi = max(values);
    if lo == hi
        padding = max(abs(lo) * 0.05, 0.05);
    else
        padding = 0.05 * (hi - lo);
    end
    limits = [lo - padding, hi + padding];
end

function local_legend(ax, handles, labels, style, num_columns)
    if nargin < 5
        num_columns = 1;
    end
    keep = isgraphics(handles);
    if any(keep)
        lg = legend(ax, handles(keep), labels(keep), 'Location', 'best');
        lg.NumColumns = num_columns;
        mwecmass.output.figures.style_legend(lg, style);
    end
end

function local_finalize_axes(ax, style)
    mwecmass.output.figures.apply_axes_style(ax, style);
    set(ax, 'Box', 'on', 'XGrid', 'off', 'YGrid', 'off', ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off');
end
