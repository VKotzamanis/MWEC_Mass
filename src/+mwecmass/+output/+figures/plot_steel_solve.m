function fig = plot_steel_solve(realised, config)
%PLOT_STEEL_SOLVE Stage-3 thin-shell realisation: elevation of the realised solid with its ballast level.
% Inputs: realised (results.stage3: the realised or closest-fail design with its body and status);
% config (output options). Draws the exact y = 0 section of the realised body
% (mwecmass.output.figures.realised_section_data): solid ballast below z_ballast, the shell of
% uniform thickness above it with its air void, the waterline, the ballast level and the status of
% the realisation. Nothing is offset or clipped here. Returns the figure handle.

    style = mwecmass.output.figures.presentation_style(config);
    if ~isstruct(realised) || ~isfield(realised, 'mode') || ~strcmp(realised.mode, 'thin_shell')
        error('mwecmass:figures:BadRealised', 'plot_steel_solve needs the thin_shell realisation (results.stage3).');
    end
    data = mwecmass.output.figures.realised_section_data(realised);
    if ~isempty(data.omitted)
        warning('mwecmass:figures:SectionsOmitted', ...
            '%d section heights had no usable section and are left out of the figure (first: z = %g, %s).', ...
            numel(data.omitted), data.omitted(1).z, data.omitted(1).reason);
    end
    vs = data.vs;
    profile = data.outline.profile;

    fig = mwecmass.output.figures.new_figure(style, 'tall_double_column');
    % This single-panel figure needs vertical room so axis equal can show the full hull beneath
    % the two-column legend above the axes; the shared size preset stays unchanged for other figures.
    set(fig, 'Units', 'centimeters');
    fig.Position(3:4) = [17, 15];
    set(fig, 'Name', 'Stage 3: Thin Shell');
    ax = axes(fig);
    hold(ax, 'on');
    axis(ax, 'equal');

    for k = 1:numel(data.polygons)
        p = data.polygons(k);
        if ~strcmp(p.role, 'void')
            patch(ax, p.xz(:, 1), p.xz(:, 2) + vs, role_color(p.role, style), 'EdgeColor', 'none', ...
                'FaceAlpha', 1.0, 'HandleVisibility', 'off');
        end
    end
    for k = 1:numel(data.polygons)
        p = data.polygons(k);
        if strcmp(p.role, 'void')
            patch(ax, p.xz(:, 1), p.xz(:, 2) + vs, style.fill_palette.void, 'EdgeColor', 'none', ...
                'FaceAlpha', 1.0, 'HandleVisibility', 'off');
            h_in = plot(ax, p.xz([1:end, 1], 1), p.xz([1:end, 1], 2) + vs, '--', ...
                'Color', style.fill_palette.inner_boundary, 'HandleVisibility', 'off');
            mwecmass.output.figures.style_line(h_in, style, 'boundary');
            mwecmass.output.figures.draw_hatch_strips(ax, p.xz(:, 1), p.xz(:, 2) + vs, ...
                style.hatch_spacing, style.fill_palette.hatch);
        end
    end
    h_hull = plot(ax, profile([1:end, 1], 1), profile([1:end, 1], 2) + vs, '-', ...
        'Color', style.fill_palette.boundary, 'HandleVisibility', 'off');
    mwecmass.output.figures.style_line(h_hull, style, 'boundary');

    x_limits = [min(profile(:, 1)), max(profile(:, 1))];
    x_margin = 0.05 * max(diff(x_limits), eps);
    x_range = x_limits + [-x_margin, x_margin];
    h_wl = plot(ax, x_range, [0, 0], '--', 'Color', style.fill_palette.waterline);
    mwecmass.output.figures.style_line(h_wl, style, 'reference');
    handles = [];
    labels = {};
    if data.z_ballast > data.z_range(1) && data.z_ballast < data.z_range(2)
        h_zb = plot(ax, x_range, data.z_ballast * [1 1] + vs, '-', 'Color', style.fill_palette.ballast_level);
        mwecmass.output.figures.style_line(h_zb, style, 'reference');
        handles = h_zb;
        labels = {'Ballast level, $z_{\mathrm{ballast}}$'};
    end

    present = unique({data.polygons.role});
    t = realised.design.t(isfinite(realised.design.t));
    if any(strcmp(present, 'ballast'))
        h_b = patch(ax, NaN, NaN, role_color('ballast', style), 'EdgeColor', 'none');
        handles = [h_b, handles];
        labels = [{sprintf('Solid Ballast, $\\rho_{\\mathrm{ballast}} = %.0f$ [kg/m$^3$]', realised.rho.ballast)}, labels];
    end
    if any(strcmp(present, 'wall'))
        h_s = patch(ax, NaN, NaN, role_color('wall', style), 'EdgeColor', 'none');
        handles = [handles, h_s];
        labels = [labels, {sprintf('Thin Shell, $t$ = %.1f [mm]', 1000 * t(1))}];
    end
    if any(strcmp(present, 'void'))
        h_a = patch(ax, NaN, NaN, style.fill_palette.void, 'EdgeColor', style.fill_palette.hatch, 'LineStyle', '--');
        handles = [handles, h_a];
        labels = [labels, {sprintf('Hollow Volume, $\\rho_{\\mathrm{air}} = %.1f$ [kg/m$^3$]', realised.rho.air)}];
    end
    handles = [handles, h_wl];
    labels = [labels, {'Waterline, $Z = 0$ m'}];
    lg = legend(ax, handles, labels, 'Location', 'northoutside', 'NumColumns', 2, 'AutoUpdate', 'off');
    mwecmass.output.figures.style_legend(lg, style);
    lg.FontSize = min(style.font_size.legend, 8.5);

    if strcmp(data.status, 'accepted')
        colour = style.status_palette.ok;
    else
        colour = style.status_palette.bad;
    end
    z_world = profile(:, 2) + vs;
    z_limits = [min([z_world; 0]), max([z_world; 0])];
    z_margin = 0.05 * max(diff(z_limits), eps);
    h_status = text(ax, 0.02, 0.02, strrep(strrep(data.status_lines, '_', '\_'), '%', '\%'), 'Units', 'normalized', ...
        'VerticalAlignment', 'bottom', 'BackgroundColor', 'w', 'Margin', 1, 'Color', colour, 'FontWeight', 'bold');
    mwecmass.output.figures.style_text(h_status, style, 'annotation');

    xlabel(ax, {'X [m]', 'Stage 3: Thin Shell Construction'});
    ylabel(ax, 'Z [m]');
    mwecmass.output.figures.apply_axes_style(ax, style);
    set(ax, 'Box', 'on', 'XGrid', 'off', 'YGrid', 'off', 'XMinorGrid', 'off', 'YMinorGrid', 'off');
    xlim(ax, x_limits + [-1.1, 1.1] * x_margin);
    ylim(ax, z_limits + [-1.1, 1.1] * z_margin);
    hold(ax, 'off');

    out_dir = mwecmass.output.output_dir('thin_shell');
    try
        files = mwecmass.output.figures.export_figure(fig, 'Steel_Solve', config, out_dir);
        fprintf('      Steel-fill plot saved: %s\n', files{1});
    catch err
        warning('mwecmass:figures:SaveFailed', 'Could not save Steel_Solve: %s', err.message);
    end
end

function c = role_color(role, style)
    switch role
        case 'ballast'
            c = style.fill_palette.solid_material;
        case 'wall'
            c = style.fill_palette.jacket_material;
        otherwise
            error('mwecmass:figures:UnknownRole', 'Unknown section role ''%s'' for a thin-shell figure.', role);
    end
end
