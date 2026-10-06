function hull_at_draft(result_file, output_folder, N_grid, view_az, view_el)
%HULL_AT_DRAFT Visualise the hull at a given vertical shift, for comparison against the WAMIT pre-processor.
% Inputs: result_file (path to .mat with results/final_props); output_folder (default
% Output/diagnostics/); N_grid (default 40, panel density); view_az (default 35 deg); view_el
% (default 22 deg). Produces a 3D plot in body-fixed frame with translucent waterline plane at
% z_body = -vertical_shift, matching WAMIT GDF conventions. Writes PNG to output_folder.

    if nargin < 3 || isempty(N_grid),  N_grid  = 40; end
    if nargin < 4 || isempty(view_az), view_az = 35; end
    if nargin < 5 || isempty(view_el), view_el = 22; end
    if nargin < 2 || isempty(output_folder)
        output_folder = mwecmass.output.output_dir('diagnostics');
    end
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end

    [results, final_props] = mwecmass.output.load_results(result_file);
    vs = final_props.vertical_shift;      % [m] realised design draft (world-frame vertical shift)
    ms2_path = results.config.ms2_file;   % absolute path, as written by build_config.m
    if ~isfile(ms2_path)
        error('hull_at_draft:MissingDeck', 'Hull-deck file not found: %s', ms2_path);
    end

    %% ═══════════════════════════════════════════════════════════════════
    %%  BUILD MESH (full hull, NOT trimmed)
    %% ═══════════════════════════════════════════════════════════════════

    parser = mwecmass.geometry.MS2Parser.parse(ms2_path);

    pan_opts = struct( ...
        'trim_wl',         false, ...  % keep full hull (wet + dry)
        'close_gaps',      false, ...
        'cosine_spacing',  false, ...
        'wl_conform',      false, ...
        'verbose',         false, ...
        'quarter_body',    false, ...
        'half_body',       false);

    mesh = mwecmass.mesh.generate(parser, vs, N_grid, N_grid, pan_opts);

    % Panelizer returns world coords.  mwecmass.mesh.generate.m:612 does
    %     S(:, :, 3) = S(:, :, 3) + draft;
    % i.e. world_z = body_z + vs.  Inverse: body_z = world_z − vs.
    verts_body      = mesh.vertices;
    verts_body(:,3) = verts_body(:,3) - vs;         % world → body
    z_wl_body       = -vs;                          % waterline in body frame

    %% ═══════════════════════════════════════════════════════════════════
    %%  PANEL CENTROID & WET/DRY CLASSIFICATION
    %% ═══════════════════════════════════════════════════════════════════

    n_panels = size(mesh.panels, 1);
    cent_z   = zeros(n_panels, 1);
    for k = 1:n_panels
        vidx = mesh.panels(k, :);
        cent_z(k) = mean(verts_body(vidx, 3));
    end

    wet_mask = cent_z <= z_wl_body;        % panel below WL → wetted
    n_wet    = sum(wet_mask);
    n_dry    = n_panels - n_wet;

    z_min = min(verts_body(:,3));  z_max = max(verts_body(:,3));

    fprintf('\n  Hull z-extent (body frame) : [%+.3f, %+.3f] m\n', z_min, z_max);
    fprintf('  vertical_shift             : %+0.4f m\n', vs);
    fprintf('  Waterline (body frame)     : z_body = %+0.4f m\n', z_wl_body);
    fprintf('  Panels total / wet / dry   : %d / %d / %d\n', n_panels, n_wet, n_dry);
    fprintf('  Wetted z-span              : [%+.3f, %+.3f] m  (= %.3f m of hull below WL)\n', ...
            z_min, z_wl_body, z_wl_body - z_min);

    %% ═══════════════════════════════════════════════════════════════════
    %%  PLOT
    %% ═══════════════════════════════════════════════════════════════════

    c_wet = [0.20 0.45 0.78];          % submerged: blue
    c_dry = [0.92 0.84 0.65];          % above WL : wheat

    fig = figure('Color', 'white', 'Position', [60 60 1300 800]);

    ax_main = subplot(1, 2, 1);  hold(ax_main, 'on');
    % Wet patches (one patch call per group → fast)
    patch(ax_main, ...
          'Vertices', verts_body, ...
          'Faces',    mesh.panels(wet_mask, :), ...
          'FaceColor', c_wet, 'EdgeColor', [0.15 0.20 0.30], ...
          'EdgeAlpha', 0.35, 'FaceAlpha', 0.95, 'LineWidth', 0.4);
    patch(ax_main, ...
          'Vertices', verts_body, ...
          'Faces',    mesh.panels(~wet_mask, :), ...
          'FaceColor', c_dry, 'EdgeColor', [0.45 0.40 0.30], ...
          'EdgeAlpha', 0.35, 'FaceAlpha', 0.95, 'LineWidth', 0.4);

    % Waterline plane
    xr = [min(verts_body(:,1))-0.5, max(verts_body(:,1))+0.5];
    yr = [min(verts_body(:,2))-0.5, max(verts_body(:,2))+0.5];
    patch(ax_main, ...
          'Vertices', [xr(1) yr(1) z_wl_body; xr(2) yr(1) z_wl_body; ...
                       xr(2) yr(2) z_wl_body; xr(1) yr(2) z_wl_body], ...
          'Faces', 1:4, ...
          'FaceColor', [0.45 0.65 0.92], 'FaceAlpha', 0.20, ...
          'EdgeColor', [0.20 0.35 0.70], 'LineWidth', 1.2);

    % Annotate
    text(ax_main, xr(2), yr(2), z_wl_body, ...
         sprintf('  WL: z_{body}=%+.4f m', z_wl_body), ...
         'FontSize', 10, 'Color', [0.20 0.35 0.70], 'VerticalAlignment', 'middle');

    xlabel(ax_main, 'x_{body} [m]', 'FontSize', 12);
    ylabel(ax_main, 'y_{body} [m]', 'FontSize', 12);
    zlabel(ax_main, 'z_{body} [m]', 'FontSize', 12);
    title(ax_main, sprintf('3D view  (v_s = %+0.4f m)', vs), 'FontSize', 13);
    axis(ax_main, 'equal');  grid(ax_main, 'on');  box(ax_main, 'on');
    view(ax_main, view_az, view_el);
    camlight(ax_main, 'headlight');  lighting(ax_main, 'gouraud');
    set(ax_main, 'FontName', 'Times New Roman', 'FontSize', 11);

    % Side view (xz) to make the waterline cut crystal clear
    ax_side = subplot(1, 2, 2);  hold(ax_side, 'on');
    patch(ax_side, ...
          'Vertices', verts_body, 'Faces', mesh.panels(wet_mask, :), ...
          'FaceColor', c_wet, 'EdgeColor', [0.15 0.20 0.30], ...
          'EdgeAlpha', 0.25, 'FaceAlpha', 0.95, 'LineWidth', 0.4);
    patch(ax_side, ...
          'Vertices', verts_body, 'Faces', mesh.panels(~wet_mask, :), ...
          'FaceColor', c_dry, 'EdgeColor', [0.45 0.40 0.30], ...
          'EdgeAlpha', 0.25, 'FaceAlpha', 0.95, 'LineWidth', 0.4);

    % Waterline horizontal line in xz-view
    yline(ax_side, z_wl_body, '-', 'Color', [0.20 0.35 0.70], 'LineWidth', 1.8, ...
          'Label', sprintf('z_{body} = %+.4f m', z_wl_body), ...
          'LabelHorizontalAlignment', 'right', 'LabelVerticalAlignment', 'top');
    % z = 0 marker (the .ms2 body origin) for reference
    yline(ax_side, 0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0, ...
          'Label', 'z_{body} = 0 (origin)', ...
          'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom');

    xlabel(ax_side, 'x_{body} [m]', 'FontSize', 12);
    ylabel(ax_side, 'z_{body} [m]', 'FontSize', 12);
    title(ax_side, 'Side view (xz, looking down y)', 'FontSize', 13);
    axis(ax_side, 'equal');  grid(ax_side, 'on');  box(ax_side, 'on');
    view(ax_side, 0, 0);
    set(ax_side, 'FontName', 'Times New Roman', 'FontSize', 11);

    sgtitle(fig, sprintf( ...
        'Hull at v_s = %+0.4f m  —  body-fixed frame  —  blue = wet, beige = dry', vs), ...
        'FontName', 'Times New Roman', 'FontSize', 14);

    %% ═══════════════════════════════════════════════════════════════════
    %%  SAVE
    %% ═══════════════════════════════════════════════════════════════════

    fname = sprintf('WAMIT_GeomVerify_vs%+0.4f_%s.png', vs, char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
    exportgraphics(fig, fullfile(output_folder, fname), 'Resolution', 300);
    close(fig);   % 16 GB machine: never leave figures open
    fprintf('\n  Figure saved: %s\n', fullfile(output_folder, fname));

    [~, deck_name, deck_ext] = fileparts(ms2_path);
    fprintf('\n  HOW TO VERIFY AGAINST WAMIT\n');
    fprintf('    1. Open your WAMIT pre-processor / viewer with %s%s''s .gdf loaded.\n', deck_name, deck_ext);
    fprintf('    2. Set XTRIMWL = %+0.4f in the .cfg and reload.\n', -vs);
    fprintf('    3. The wetted region in WAMIT should match the BLUE patches here.\n');
    fprintf('       The dry region (BEIGE) should be above WAMIT''s trim plane.\n');
end
