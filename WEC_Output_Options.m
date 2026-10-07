function out = WEC_Output_Options()

%WEC_OUTPUT_OPTIONS Define output selection, export formats, and shared figure styling.
% Syntax: out = WEC_Output_Options(); no inputs.
% Output groups save/export/style control reports and figures; sizes are [width height] cm,
% line/font widths are points, RGB values lie in [0,1], and hatch spacing is metres.

%% Outputs selected by pipeline stage
out.save.stage1.density_2d = false;            % bool: save the Stage-1 equivalent-density 2-D figure.
out.save.stage1.draft_landscape = true;        % bool: save the Stage-1 skipped-mode draft landscape.

out.save.stage2.density_3d = false;            % bool: save the Stage-2 optimiser-density 3-D figure.
out.save.stage2.cross_section = true;          % bool: save the Stage-2 final cross-section figure.
out.save.stage2.summary_log = true;            % bool: print the Stage-2 summary log.

out.save.convergence = true;                   % bool: optimiser convergence history over stages 1 and 2.

out.save.stage3.steel_solve = true;            % bool: save the thin-shell steel-solve diagnostic figure.
out.save.stage3.steel_solve_log = true;        % bool: print the thin-shell steel-solve diagnostic log.
out.save.stage3.precast_midplane = true;       % bool: save the modular-precast midplane figure.
out.save.stage3.precast_strips = true;         % bool: save the modular-precast strip-plan figure.
out.save.stage3.realised_density_3d = false;   % bool: save the realised-density 3-D figure.
out.save.stage3.final_results_log = true;      % bool: print the final-results log.
out.save.stage3.diagnostic_panels_log = true;  % bool: print the diagnostic-panels log.

out.save.hydrodynamics.coefficients = true;          % bool: save added-mass and radiation-damping curves.
out.save.hydrodynamics.raos = true;                   % bool: save response-amplitude-operator curves.
out.save.hydrodynamics.added_mass_vs_draft = true;    % bool: save infinite-frequency added mass by draft.
out.save.hydrodynamics.mesh_diagnostic = false;       % bool: save the boundary-element mesh diagnostic.
out.save.hydrodynamics.panel_normals = true;           % bool: save panel-normal directions.
out.save.hydrodynamics.cache_rewrite = true;           % bool: rewrite an enriched hydrodynamic cache.

out.save.results_mat = true;                    % bool: write the results and final-properties MAT-file.

out.save.diagnostics.hull_at_draft = false;       % bool: save the hull-at-draft geometry verification figure.
out.save.diagnostics.uhpc_mass_balance = false;    % bool: save the UHPC mass-balance four-panel figure (modular_precast only).
out.save.diagnostics.stage_animations = false;     % bool: save the per-stage optimisation animation GIFs (slow).

%% File delivery
out.export.formats = {'png', 'fig', 'pdf'};      % cellstr: export each figure in these formats.
out.export.dpi = 450;                            % dpi: raster PNG resolution.
out.export.pdf_content = 'vector';               % char: PDF content type passed to exportgraphics.
out.export.timestamp_filenames = false;          % bool: append yyyyMMdd_HHmmss to exported figure names.
out.console_echo = true;                         % bool: allow enabled console reports to print.
out.output_dir = 'Output';                       % char: repository-relative folder for run output.

%% Shared presentation by role
out.style.font_name = 'Times New Roman';         % char: font used for figure text.
out.style.mono_font_name = 'Courier New';        % char: monospace font used for numeric cards.
out.style.font_size.title = 14;                  % pt: figure and axes-title text.
out.style.font_size.axes = 12;                   % pt: axis-label text.
out.style.font_size.tick_label = 12;             % pt: axis tick-label text, equal to axis labels.
out.style.font_size.colorbar = 10.8;             % pt: colorbar tick-label text.
out.style.font_size.legend = 9.8;                % pt: legend-entry text.
out.style.font_size.annotation = 11;             % pt: in-figure annotation text.
out.style.font_size.small_multiple = 8;          % pt: title, tick, and axis-label text within one panel.
out.style.colormap = 'cividis';                  % char: density colormap ('cividis', 'parula', or 'jet').

out.style.line_width.main = 2.0;                 % pt: main-curve width used by existing figures.
out.style.line_width.curve = 2.0;                % pt: primary plotted-curve width.
out.style.line_width.boundary = 1.4;             % pt: hull and material-boundary width.
out.style.line_width.reference = 1.2;            % pt: reference-curve width.
out.style.line_width.axes = 0.6;                 % pt: axes-box and tick width.
out.style.line_width.grid = 0.6;                 % pt: grid-line width.
out.style.line_width.hatch = 0.4;                % pt: void-region hatch-stroke width.
out.style.marker.size = 6;                       % pt: plotted marker size.

out.style.figure_size.single_column = [8.5, 6.65]; % cm: [width height]; read by plot_steel_solve.m.
out.style.figure_size.double_column = [17, 6.65];  % cm: [width height] used by plot_mesh_diagnostic.m.
out.style.figure_size.tall_single = [12.5, 10];    % cm: [width height] for a tall single panel.
out.style.figure_size.tall_double_column = [17, 12.5]; % cm: [width height] for a tall multi-panel layout.
out.style.layout.tile_spacing = 'compact';       % char: tiledlayout TileSpacing setting.
out.style.layout.padding = 'compact';            % char: tiledlayout Padding setting.

out.style.grid.alpha = 0.4;                      % [-]: grid-line opacity from 0 to 1.
out.style.grid.line_style = ':';                 % char: grid-line style.
out.style.tick_label_interpreter = 'latex';      % char: interpreter for tick labels and axes text.
out.style.hatch_spacing = 0.06;                  % m: void-region hatch-line separation.

out.style.status_palette.ok = [0.10, 0.55, 0.10]; % RGB: acceptable realisation-card status colour.
out.style.status_palette.bad = [0.80, 0.15, 0.15]; % RGB: failed realisation-card status colour.

out.style.color.cg = [1.00, 0.00, 0.00];         % RGB: centre-of-gravity marker fill.
out.style.color.cb = [0.00, 0.00, 1.00];         % RGB: centre-of-buoyancy marker fill; same three readers as color.cg.
out.style.color.series_a = [0.00, 0.00, 1.00];   % RGB: first curve in a two-series comparison panel.
out.style.color.series_b = [1.00, 0.00, 0.00];   % RGB: second curve in a two-series comparison panel.

out.style.fill_palette.solid_material = [0.45, 0.46, 0.50]; % RGB: solid structural material.
out.style.fill_palette.jacket_material = [0.74, 0.76, 0.80]; % RGB: jacket structural material.
out.style.fill_palette.ballast_material = [0.55, 0.55, 0.60]; % RGB: filled strip material.
out.style.fill_palette.shell = [0.82, 0.82, 0.82]; % RGB: shell annulus material.
out.style.fill_palette.void = [1.00, 1.00, 1.00]; % RGB: void or air region.
out.style.fill_palette.hatch = [0.50, 0.52, 0.58]; % RGB: void-region hatch lines.
out.style.fill_palette.boundary = [0.10, 0.10, 0.10]; % RGB: outer material boundary.
out.style.fill_palette.inner_boundary = [0.30, 0.30, 0.32]; % RGB: inner shell boundary.
out.style.fill_palette.waterline = [0.15, 0.55, 0.95]; % RGB: still-water line.
out.style.fill_palette.ballast_level = [0.95, 0.55, 0.10]; % RGB: fill-level line.
out.style.fill_palette.wall_boundary = [0.80, 0.15, 0.10]; % RGB: precast wall boundary.

out.style.color.dof = [0.12 0.47 0.71; 0.20 0.63 0.17; 0.89 0.10 0.11];           % RGB rows for surge, heave, and pitch.
out.style.color.dof_secondary = [0.40 0.65 0.85; 0.55 0.78 0.52; 0.95 0.50 0.50]; % RGB secondary rows for surge, heave, and pitch.
out.style.color.dof_pair = [0.58 0.40 0.74; 1.00 0.50 0.05; 0.55 0.34 0.29];      % RGB rows for DOF pairs (1,3), (1,5), and (3,5).
out.style.color.dof_pair_secondary = [0.75 0.60 0.85; 1.00 0.73 0.47; 0.73 0.55 0.50]; % RGB secondary rows for DOF pairs (1,3), (1,5), and (3,5).
out.style.color.reference = [0.3, 0.3, 0.3];       % RGB neutral colour for thresholds, asymptotes, and natural periods.
out.style.color.band = [0.85, 0.92, 1.0];          % RGB fill colour for operational-period bands.
out.style.color.diagnostic = [0.90, 0.20, 0.10];   % RGB highlight colour for mesh and geometry diagnostics.
out.style.color.mesh_edge = [0.35, 0.35, 0.35];    % RGB BEM-panel edge outline colour.
out.style.color.mesh_quad.fill = [0.65, 0.82, 1.0]; % RGB: waterplane-lid quadrilateral-panel fill; read by plot_mesh_diagnostic.m.
out.style.color.mesh_quad.edge = [0.10, 0.30, 0.70]; % RGB: waterplane-lid quadrilateral-panel edge; read by plot_mesh_diagnostic.m.
out.style.color.mesh_tri.fill = [0.65, 0.95, 0.75]; % RGB: waterplane-lid triangular-panel fill; read by plot_mesh_diagnostic.m.
out.style.color.mesh_tri.edge = [0.05, 0.45, 0.20]; % RGB: waterplane-lid triangular-panel edge; read by plot_mesh_diagnostic.m.
end
