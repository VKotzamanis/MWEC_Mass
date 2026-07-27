%% PLOT_GEOMETRY_AT_DRAFT  Visualise the C0 hull at a given vertical shift
%
%  Produces a 3D plot of the FULL hull (wet + dry) in the .ms2 NATURAL
%  body-fixed frame (i.e., the frame your WAMIT GDF is drawn in).  Adds
%  a translucent waterline plane at z_body = -vertical_shift so you can
%  visually compare what the WAMIT setup is "seeing" against the body
%  pictured in your WAMIT pre-processor.
%
%  CONVENTIONS (must match WAMIT setup)
%    XBODY      = (0, 0, 0)        — body origin = global origin
%    XTRIMWL    = -vertical_shift  — waterline trim plane in WAMIT
%    Hence the body is drawn UNTRANSLATED in the .ms2 natural frame;
%    only the WL plane moves.
%
%  USAGE
%    Edit `vs` in §1 below, then run.
%
%  OUTPUT
%    WAMIT_GeomVerify_vs<+/-X.XXXX>_<stamp>.png  (300 dpi)

clear; clc; close all;

%% ═══════════════════════════════════════════════════════════════════
%%  §1  USER CONFIG
%% ═══════════════════════════════════════════════════════════════════

vs           = 1.4286;            % vertical_shift [m]  ← edit this
ms2_file     = 'C0.ms2';
N_grid       = 40;                % Nu = Nv = N_grid (panel density)
view_az      = 35;                % azimuth [deg]
view_el      = 22;                % elevation [deg]

%% ═══════════════════════════════════════════════════════════════════
%%  §2  BUILD MESH (full hull, NOT trimmed)
%% ═══════════════════════════════════════════════════════════════════

suite_dir = fileparts(mfilename('fullpath'));
if ~isempty(suite_dir); cd(suite_dir); end

parser = WEC_MS2_Parser.parse(ms2_file);

pan_opts = struct( ...
    'trim_wl',         false, ...  % keep full hull (wet + dry)
    'close_gaps',      false, ...
    'cosine_spacing',  false, ...
    'wl_conform',      false, ...
    'verbose',         false, ...
    'quarter_body',    false, ...
    'half_body',       false);

mesh = WEC_Panelizer.generate(parser, vs, N_grid, N_grid, pan_opts);

% Panelizer returns world coords.  WEC_Panelizer.m:612 does
%     S(:, :, 3) = S(:, :, 3) + draft;
% i.e. world_z = body_z + vs.  Inverse: body_z = world_z − vs.
verts_body      = mesh.vertices;
verts_body(:,3) = verts_body(:,3) - vs;         % world → body
z_wl_body       = -vs;                          % waterline in body frame

%% ═══════════════════════════════════════════════════════════════════
%%  §3  PANEL CENTROID & WET/DRY CLASSIFICATION
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
%%  §4  PLOT
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
set(ax_main, 'FontName', 'Times New Roman', 'FontSize', 11, ...
             'GridLineStyle', ':', 'GridAlpha', 0.4, ...
             'TickDir', 'out', 'Layer', 'top');

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
set(ax_side, 'FontName', 'Times New Roman', 'FontSize', 11, ...
             'GridLineStyle', ':', 'GridAlpha', 0.4, ...
             'TickDir', 'out', 'Layer', 'top');

sgtitle(fig, sprintf( ...
    'C0 hull at v_s = %+0.4f m  —  body-fixed frame  —  blue = wet, beige = dry', vs), ...
    'FontName', 'Times New Roman', 'FontSize', 14);

%% ═══════════════════════════════════════════════════════════════════
%%  §5  SAVE
%% ═══════════════════════════════════════════════════════════════════

% Written to the shared Plots/ directory alongside every other suite figure.
fname = sprintf('WAMIT_GeomVerify_vs%+0.4f', vs);
WEC_Visualization.save_figure(fig, fname);

fprintf('\n  HOW TO VERIFY AGAINST WAMIT\n');
fprintf('    1. Open your WAMIT pre-processor / viewer with C0.gdf loaded.\n');
fprintf('    2. Set XTRIMWL = %+0.4f in the .cfg and reload.\n', -vs);
fprintf('    3. The wetted region in WAMIT should match the BLUE patches here.\n');
fprintf('       The dry region (BEIGE) should be above WAMIT''s trim plane.\n');
