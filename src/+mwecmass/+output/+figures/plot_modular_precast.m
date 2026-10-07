function plot_modular_precast(realised, config)
%PLOT_MODULAR_PRECAST Stage-3 UHPC realisation: elevation and per-module plan views of the realised solid.
% Inputs: realised (results.stage3: the realised or closest-fail design with its body, status and
% per-module volumes); config (output options). Both figures draw only the exact sections of the
% realised body (mwecmass.output.figures.realised_section_data): the y = 0 elevation with the walls,
% the voids, the ballast level and the solid modules as they are built, and the plan section at the
% bottom and top of every module that is not a solid module. The status of the realisation is
% written on the elevation.

    style = mwecmass.output.figures.presentation_style(config);
    if ~isstruct(realised) || ~isfield(realised, 'mode') || ~strcmp(realised.mode, 'modular_precast')
        error('mwecmass:figures:BadRealised', 'plot_modular_precast needs the modular_precast realisation (results.stage3).');
    end
    data = mwecmass.output.figures.realised_section_data(realised);
    warn_omitted(data);
    vs = data.vs;
    N = numel(realised.modules);

    %% Figure 1: elevation in the y = 0 plane
    fig1 = mwecmass.output.figures.new_figure(style, 'tall_double_column');
    set(fig1, 'Name', 'Stage 3: Modular Precast Elevation');
    tl1 = tiledlayout(fig1, 4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    ax_comp = nexttile(tl1, 1);
    axis(ax_comp, 'off');
    hold(ax_comp, 'on');
    ax1 = nexttile(tl1, 2, [3 1]);
    hold(ax1, 'on');

    draw_polygons(ax1, data.polygons, vs, style);
    profile = data.outline.profile;
    h_outline = plot(ax1, profile([1:end, 1], 1), profile([1:end, 1], 2) + vs, '-', ...
        'Color', style.fill_palette.boundary, 'HandleVisibility', 'off');
    mwecmass.output.figures.style_line(h_outline, style, 'boundary');

    x_range = [min(profile(:, 1)) - 0.30, max(profile(:, 1)) + 0.30];
    for z = data.edges(:)'
        line(ax1, x_range, [z z] + vs, 'Color', [0.55 0.55 0.55], 'LineStyle', ':', ...
            'LineWidth', style.line_width.grid, 'HandleVisibility', 'off');
    end
    h_wl = plot(ax1, x_range, [0 0], '--', 'Color', style.fill_palette.waterline);
    mwecmass.output.figures.style_line(h_wl, style, 'reference');
    handles = h_wl;
    labels = {'Waterline, $Z = 0$ m'};
    if data.z_ballast > data.z_range(1) && data.z_ballast < data.z_range(2)
        h_zb = plot(ax1, x_range, data.z_ballast * [1 1] + vs, '-', 'Color', style.fill_palette.ballast_level);
        mwecmass.output.figures.style_line(h_zb, style, 'reference');
        handles(end + 1) = h_zb;
        labels{end + 1} = 'Ballast level, $z_{\mathrm{ballast}}$';
    end

    x_annot = x_range(2) + 0.10;
    for i = 1:N
        z_mid = 0.5 * (data.edges(i) + data.edges(i + 1)) + vs;
        h_lbl = text(ax1, x_annot, z_mid, module_label(realised.modules(i), i), ...
            'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
            'Color', [0.15 0.15 0.15], 'Clipping', 'on');
        mwecmass.output.figures.style_text(h_lbl, style, 'annotation');
        h_lbl.FontSize = 8;
    end

    [h_mat, l_mat] = material_proxies(ax1, data.polygons, style, realised);
    lg1 = legend(ax1, [h_mat, handles], [l_mat, labels], 'Location', 'eastoutside', 'AutoUpdate', 'off');
    mwecmass.output.figures.style_legend(lg1, style);
    lg1.FontSize = min(style.font_size.legend, 8.5);
    xlabel(ax1, {'X [m]', 'Stage 3: Modular Precast Construction'});
    ylabel(ax1, 'Z [m]');
    axis(ax1, 'equal');
    mwecmass.output.figures.apply_axes_style(ax1, style);
    set(ax1, 'TickDir', 'out', 'Box', 'on', 'XGrid', 'off', 'YGrid', 'off', ...
        'XMinorGrid', 'off', 'YMinorGrid', 'off');
    z_world = profile(:, 2) + vs;
    z_pad = 0.05 * (max(z_world) - min(z_world));
    xlim(ax1, [x_range(1), x_annot + max(1.25, 0.40 * (max(profile(:, 1)) - min(profile(:, 1))))]);
    ylim(ax1, [min(z_world) - z_pad, max(z_world) + z_pad]);
    draw_status(ax1, data, style);

    h_comp = gobjects(N, 1);
    comp_labels = cell(N, 1);
    for i = 1:N
        h_comp(i) = plot(ax_comp, NaN, NaN, 'LineStyle', 'none', 'Marker', 'none', 'HandleVisibility', 'on');
        comp_labels{i} = composition_label(realised.modules(i), i);
    end
    lg_comp = legend(ax_comp, h_comp, comp_labels, 'Location', 'north', 'NumColumns', ceil(N / 2), ...
        'Orientation', 'vertical', 'AutoUpdate', 'off');
    mwecmass.output.figures.style_legend(lg_comp, style);
    lg_comp.FontSize = 7.5;
    lg_comp.Box = 'off';
    lg_comp.ItemTokenSize = [6 5];
    hold(ax_comp, 'off');
    save_figure(fig1, 'WEC_Constructability_XZ', config);

    %% Figure 2: plan sections at the bottom and top of the modules that are not solid modules
    panels = data.strip_panels;
    if isempty(panels)
        fprintf('      Figure 2 skipped: no plan section of a hollow or ballast module.\n');
        return
    end
    n_cols = max([panels.col]);
    fig2 = mwecmass.output.figures.new_figure(style, 'tall_double_column');
    set(fig2, 'Name', 'Stage 3: Modular Precast Plan Views');
    tl2 = tiledlayout(fig2, 2, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl2, 'Precast Construction: Plan Sections of the Modules');
    mwecmass.output.figures.apply_layout_style(tl2, style);
    for q = 1:numel(panels)
        pn = panels(q);
        item = data.plan(pn.plan);
        ax = nexttile(tl2, (pn.row - 1) * n_cols + pn.col);
        hold(ax, 'on');
        patch(ax, item.outer(:, 1), item.outer(:, 2), role_color(item.role, style), 'EdgeColor', 'none');
        if ~item.solid
            h_void = patch(ax, item.inner(:, 1), item.inner(:, 2), style.fill_palette.void, ...
                'EdgeColor', style.fill_palette.inner_boundary, 'LineStyle', '--');
            mwecmass.output.figures.style_line(h_void, style, 'boundary');
        end
        h_out = plot(ax, item.outer([1:end, 1], 1), item.outer([1:end, 1], 2), '-', ...
            'Color', style.fill_palette.boundary);
        mwecmass.output.figures.style_line(h_out, style, 'boundary');
        m = realised.modules(pn.module);
        h_u = plot(ax, NaN, NaN, '-', 'Color', [0 0 0], 'LineWidth', style.line_width.boundary);
        h_a = plot(ax, NaN, NaN, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', style.line_width.boundary);
        lg = legend(ax, [h_u, h_a], ...
            {sprintf('$V_{\\mathrm{UHPC}} = %.3g$ [m$^{3}$]', m.V_uhpc), ...
             sprintf('$V_{\\mathrm{air}} = %.3g$ [m$^{3}$]', m.V_air)}, ...
            'Location', 'northoutside', 'NumColumns', 2, 'Orientation', 'horizontal', 'AutoUpdate', 'off');
        mwecmass.output.figures.style_legend(lg, style);
        set(lg, 'FontSize', 6.5, 'ItemTokenSize', [6 5], 'Box', 'on');
        pad_x = 0.05 * (max(item.outer(:, 1)) - min(item.outer(:, 1)));
        pad_y = 0.05 * (max(item.outer(:, 2)) - min(item.outer(:, 2)));
        xlim(ax, [min(item.outer(:, 1)) - pad_x, max(item.outer(:, 1)) + pad_x]);
        ylim(ax, [min(item.outer(:, 2)) - pad_y, max(item.outer(:, 2)) + pad_y]);
        axis(ax, 'equal');
        xlabel(ax, {'X [m]', sprintf('Module %d, %s view', pn.module, pn.end)});
        ylabel(ax, 'Y [m]');
        mwecmass.output.figures.apply_axes_style(ax, style, 'small_multiple');
        set(ax, 'Box', 'on', 'TickDir', 'out', 'XGrid', 'off', 'YGrid', 'off', ...
            'XMinorGrid', 'off', 'YMinorGrid', 'off');
        hold(ax, 'off');
    end
    save_figure(fig2, 'WEC_Constructability_Strips', config);
end

function draw_polygons(ax, polygons, vs, style)
% Material first, voids on top with their dashed boundary and hatch.
    for k = 1:numel(polygons)
        p = polygons(k);
        if strcmp(p.role, 'void')
            continue
        end
        patch(ax, p.xz(:, 1), p.xz(:, 2) + vs, role_color(p.role, style), 'EdgeColor', 'none', ...
            'FaceAlpha', 1.0, 'HandleVisibility', 'off');
    end
    for k = 1:numel(polygons)
        p = polygons(k);
        if ~strcmp(p.role, 'void')
            continue
        end
        h = patch(ax, p.xz(:, 1), p.xz(:, 2) + vs, style.fill_palette.void, ...
            'EdgeColor', style.fill_palette.inner_boundary, 'LineStyle', '--', 'HandleVisibility', 'off');
        mwecmass.output.figures.style_line(h, style, 'boundary');
        mwecmass.output.figures.draw_hatch_strips(ax, p.xz(:, 1), p.xz(:, 2) + vs, ...
            style.hatch_spacing, style.fill_palette.hatch);
    end
end

function c = role_color(role, style)
    switch role
        case 'solid_module'
            c = style.fill_palette.solid_material;
        case 'ballast'
            c = style.fill_palette.ballast_material;
        case 'wall'
            c = style.fill_palette.jacket_material;
        case 'void'
            c = style.fill_palette.void;
        otherwise
            error('mwecmass:figures:UnknownRole', 'Unknown section role ''%s''.', role);
    end
end

function [h, labels] = material_proxies(ax, polygons, style, realised)
% Legend entries for the roles that occur in the figure, in a fixed order.
    present = unique({polygons.role});
    order = {'solid_module', 'ballast', 'wall', 'void'};
    names = {'Solid Module (Wall)', 'Ballast (Solid UHPC)', 'UHPC', ...
             sprintf('Void, $\\rho_{\\mathrm{air}} = %.1f$ [kg/m$^3$]', realised.rho.air)};
    h = gobjects(1, 0);
    labels = {};
    for q = 1:numel(order)
        if ~any(strcmp(present, order{q}))
            continue
        end
        if strcmp(order{q}, 'void')
            hq = patch(ax, NaN, NaN, style.fill_palette.void, 'EdgeColor', style.fill_palette.inner_boundary, ...
                'LineStyle', '--');
            mwecmass.output.figures.style_line(hq, style, 'boundary');
        else
            hq = patch(ax, NaN, NaN, role_color(order{q}, style), 'EdgeColor', 'none');
        end
        h(end + 1) = hq; %#ok<AGROW>
        labels{end + 1} = names{q}; %#ok<AGROW>
    end
end

function txt = module_label(m, i)
    if isfinite(m.t)
        txt = sprintf('Module %d, $t$ = %.1f mm', i, 1000 * m.t);
    else
        txt = sprintf('Module %d', i);
    end
end

function txt = composition_label(m, i)
    if m.V_air > 0
        txt = sprintf('Module %d: $V_{\\mathrm{UHPC}}/V_{\\mathrm{air}} = %.1f\\%%$', i, 100 * m.V_uhpc / m.V_air);
    else
        txt = sprintf('Module %d: solid, $V_{\\mathrm{UHPC}} = %.3g$ m$^3$', i, m.V_uhpc);
    end
end

function draw_status(ax, data, style)
    if strcmp(data.status, 'accepted')
        colour = style.status_palette.ok;
    else
        colour = style.status_palette.bad;
    end
    h = text(ax, 0.02, 0.02, strrep(strrep(data.status_lines, '_', '\_'), '%', '\%'), 'Units', 'normalized', ...
        'VerticalAlignment', 'bottom', 'BackgroundColor', 'w', 'Margin', 1, 'Color', colour, 'FontWeight', 'bold');
    mwecmass.output.figures.style_text(h, style, 'annotation');
end

function warn_omitted(data)
    if ~isempty(data.omitted)
        warning('mwecmass:figures:SectionsOmitted', ...
            '%d section heights had no usable section and are left out of the figure (first: z = %g, %s).', ...
            numel(data.omitted), data.omitted(1).z, data.omitted(1).reason);
    end
end

function save_figure(fig, stem, config)
    try
        saved = mwecmass.output.figures.export_figure(fig, stem, config, mwecmass.output.output_dir('modular_precast'));
        fprintf('      Figure saved: %s\n', saved{1});
    catch err
        warning('mwecmass:figures:SaveFailed', 'Could not save %s: %s', stem, err.message);
    end
end
