function config = WEC_Configuration_Builder(params)
% WEC_CONFIGURATION_BUILDER  Build the unified config struct from driver params.
%
%   config = WEC_CONFIGURATION_BUILDER(params)
%
%   Reads every tuneable parameter from the `params` struct defined in
%   WEC_Driver.m.  No physics constants or optimisation settings are
%   hardcoded here — change them in the driver.
%
%   This function:
%     1. Copies physical constants and targets straight from params.
%     2. Parses the MS2 parametric model and extracts 2D midplane profile.
%     3. Builds a Y-span lookup table from spline cross-sections.
%     4. Reads and transforms WAMIT added-mass / damping data.
%     5. Computes an initial density guess (smooth tanh profile).
%     6. Validates internal consistency of the assembled config.
%
%   INPUT
%     params : struct defined in WEC_Driver.m (single source of truth)
%
%   OUTPUT
%     config : struct consumed by every kernel in the pipeline
%
%   DEPENDENCIES
%     WEC_MS2_Parser      — parametric geometry parser
%     WEC_Core_Functions   — spline cross-sections, polygon operations
%     WEC_File_IO         — WAMIT reader + coordinate transformation
%
%   See also: WEC_Driver, WEC_Main_Optimizer

try
    fprintf('\n========== WEC CONFIGURATION BUILDER ==========\n');

    %% §1  PHYSICAL CONSTANTS  ─────────────────────────────────────────
    %  Copied directly from driver.  Units: SI throughout.

    config.RHO_WATER    = params.RHO_WATER;       % [kg/m^3]
    config.G            = params.G;                % [m/s^2]
    config.seabed_depth = params.seabed_depth;     % [m]
    config.wamit_L      = params.wamit_L;          % [m]
    config.ms2_file     = params.ms2_file;         % for diagnostic reporting

    fprintf('  Constants: rho = %.0f kg/m^3, g = %.2f m/s^2\n', ...
            config.RHO_WATER, config.G);

    %% §2  GEOMETRY — MS2 PARAMETRIC MODEL  ──────────────────────────────
    %  Parse the .ms2 file into evaluatable B-spline surfaces and curves.
    %  Extract the midplane profile for the 2D surrogate.
    %  Compute total hull volume via strip integration on spline cross-sections.
    %
    %  VARIABLE DICTIONARY (outputs stored in config)
    %    ms2_model           : WEC_MS2_Parser object (evaluatable geometry)
    %    profile             : [M×2]  midplane polygon (x, z) for 2D model
    %    topology.strip_data : struct with pre-computed strip geometry
    %    total_wec_volume    : [m³]   enclosed hull volume (sanity check)
    %
    %  WHY splines instead of mesh?
    %    The parametric surfaces give exact cross-section polygons at any
    %    z-level — no mesh-density dependence, no sawtooth clipping error.
    %    Volume and mass properties are computed by strip integration of
    %    A(z) via cubic spline fitting + analytic integration.  This is
    %    O(h^4) accurate vs O(h^2) for trapz on a mesh.

    ms2_file = params.ms2_file;
    if ~exist(ms2_file, 'file')
        error('WEC:FileNotFound', 'MS2 file not found: %s', ms2_file);
    end

    fprintf('  Parsing geometry: %s\n', ms2_file);
    config.ms2_model = WEC_MS2_Parser.parse(ms2_file);
    config.profile   = WEC_Core_Functions.extractProfileMS2(config.ms2_model);

    fprintf('  Visible surfaces: %d\n', length(config.ms2_model.visible_surfs));
    fprintf('  Profile vertices: %d\n', size(config.profile, 1));

    % Store panelizer grid density (used only for .pnl export, not in optimiser)
    config.mesh_Nu = params.mesh_Nu;
    config.mesh_Nv = params.mesh_Nv;

    % Hull z-extents from ACTUAL SURFACE EVALUATION (v5.0)
    %
    %  model.extents comes from the .ms2 file header — the bounding box
    %  of ALL entities including B-spline control points.  B-splines
    %  APPROXIMATE control points (don't pass through interior ones),
    %  so control points can lie well outside the actual hull surface.
    %  For C0, a control point at z = -3.0 gives hull_z_min = -3.0,
    %  but the surface only reaches z ≈ -2.5.  This wastes a density
    %  strip on empty space below the keel.
    %
    %  FIX: Sample all visible surfaces on a coarse grid and extract
    %  the actual min/max z.  Cost: ~0.05 s (negligible at config time).
    [hull_z_lo, hull_z_hi] = WEC_Core_Functions.compute_surface_z_range(config.ms2_model);
    hull_extents_z = [hull_z_lo, hull_z_hi];
    hull_z_min = hull_extents_z(1);
    hull_z_max = hull_extents_z(2);
    fprintf('  Surface z-range: [%.4f, %.4f] m (file header: [%.4f, %.4f])\n', ...
            hull_z_lo, hull_z_hi, ...
            min(config.ms2_model.extents(3), config.ms2_model.extents(6)), ...
            max(config.ms2_model.extents(3), config.ms2_model.extents(6)));
    hp_full = WEC_HydroProperties.compute(config.ms2_model, ...
                  struct('n_quad', 20, 'verbose', false));
    config.total_wec_volume = hp_full.volume;
    config.hull_centroid    = hp_full.centroid;

    % Second volume moments (v8.1) — for inertia tensor in HAMS pipeline
    config.hull_int_x2 = hp_full.int_x2;
    config.hull_int_y2 = hp_full.int_y2;
    config.hull_int_z2 = hp_full.int_z2;

    fprintf('  Hull volume: %.4f m^3 (parametric divergence theorem)\n', config.total_wec_volume);

    %% §2b  PRECOMPUTE BOUNDARY CACHE + WATERPLANE AREA TABLE (v8.0)
    %
    %  BOUNDARY CACHE (new in v8.0):
    %    Evaluates all source surface boundary curves ONCE at n_u
    %    u-samples.  Stores the results in a cache struct that
    %    extract_isocurve_at_z reuses for every z-level query.
    %    Cost: ~0.15 s.  Eliminates all redundant eval_curve/eval_snake
    %    calls from the Aw table construction, r_min computation,
    %    and any standalone compute_submerged/compute_strip calls.
    %
    %  AW TABLE (updated from v7.0):
    %    Same adaptive sampling strategy (coarse + refinement at high
    %    |dAw/dz|), but contour extraction now uses type-dispatched
    %    algebraic solvers instead of brute-force bisection sweeps.
    %    Validated: 0.00% Aw error vs reference, 1249× speedup on C0.
    %
    %  USAGE: calculate_3d_properties interpolates Aw from this table
    %    and passes it as Aw_override to compute_submerged.  No contour
    %    tracing in the optimizer loop.

    fprintf('  Precomputing boundary cache...\n');
    n_u_cache = 100;
    config.boundary_cache = WEC_HydroProperties.precompute_boundary_cache( ...
                                config.ms2_model, n_u_cache);
    fprintf('    Sources: %d, mirrors: %d, u-samples: %d\n', ...
            length(config.boundary_cache.sources), ...
            length(config.boundary_cache.mirrors), n_u_cache);

    fprintf('  Building waterplane area table (adaptive, fast iso-z)...\n');

    % Parametric resolution (v9.0): n_coarse from target dz spacing
    if isfield(params, 'aw_table_dz') && params.aw_table_dz > 0
        n_coarse = max(20, ceil((hull_z_max - hull_z_min) / params.aw_table_dz));
    else
        n_coarse = 40;  % backward compatibility
    end
    z_coarse = linspace(hull_z_min + 1e-4, hull_z_max - 1e-4, n_coarse)';
    Aw_coarse  = zeros(n_coarse, 1);
    Ixx_coarse = zeros(n_coarse, 1);
    Iyy_coarse = zeros(n_coarse, 1);
    P_coarse   = zeros(n_coarse, 1);

    for k = 1:n_coarse
        wl_pts = WEC_HydroProperties.extract_isocurve_at_z( ...
                     config.ms2_model, z_coarse(k), n_u_cache, ...
                     config.boundary_cache);
        if ~isempty(wl_pts) && size(wl_pts, 1) >= 3
            [Aw_coarse(k), Ixx_coarse(k), Iyy_coarse(k)] = ...
                WEC_HydroProperties.waterplane_properties(wl_pts);
            x_cl = [wl_pts(:,1); wl_pts(1,1)];
            y_cl = [wl_pts(:,2); wl_pts(1,2)];
            P_coarse(k) = sum(sqrt(diff(x_cl).^2 + diff(y_cl).^2));
        end
    end

    % Pass 2: adaptive refinement where Aw changes rapidly
    dAw = abs(diff(Aw_coarse));
    dz_c = diff(z_coarse);
    grad_Aw = dAw ./ dz_c;
    grad_thresh = 0.3 * max(grad_Aw);  % refine top 30% gradient

    z_refine   = [];
    Aw_refine  = [];
    Ixx_refine = [];
    Iyy_refine = [];
    P_refine   = [];

    for k = 1:length(grad_Aw)
        if grad_Aw(k) > grad_thresh
            z_mid = linspace(z_coarse(k), z_coarse(k+1), 5)';
            z_mid = z_mid(2:end-1);
            for j = 1:length(z_mid)
                wl_pts = WEC_HydroProperties.extract_isocurve_at_z( ...
                             config.ms2_model, z_mid(j), n_u_cache, ...
                             config.boundary_cache);
                Aw_j = 0; Ixx_j = 0; Iyy_j = 0; P_j = 0;
                if ~isempty(wl_pts) && size(wl_pts, 1) >= 3
                    [Aw_j, Ixx_j, Iyy_j] = ...
                        WEC_HydroProperties.waterplane_properties(wl_pts);
                    x_cl = [wl_pts(:,1); wl_pts(1,1)];
                    y_cl = [wl_pts(:,2); wl_pts(1,2)];
                    P_j  = sum(sqrt(diff(x_cl).^2 + diff(y_cl).^2));
                end
                z_refine(end+1)   = z_mid(j); %#ok<AGROW>
                Aw_refine(end+1)  = Aw_j;     %#ok<AGROW>
                Ixx_refine(end+1) = Ixx_j;    %#ok<AGROW>
                Iyy_refine(end+1) = Iyy_j;    %#ok<AGROW>
                P_refine(end+1)   = P_j;       %#ok<AGROW>
            end
        end
    end

    % Merge coarse + refined, sort by z
    z_all_aw  = [z_coarse(:); z_refine(:)];
    Aw_all    = [Aw_coarse(:);  Aw_refine(:)];
    Ixx_all   = [Ixx_coarse(:); Ixx_refine(:)];
    Iyy_all   = [Iyy_coarse(:); Iyy_refine(:)];
    P_all     = [P_coarse(:);   P_refine(:)];

    [z_all_aw, sort_idx] = sort(z_all_aw);
    Aw_all  = Aw_all(sort_idx);
    Ixx_all = Ixx_all(sort_idx);
    Iyy_all = Iyy_all(sort_idx);
    P_all   = P_all(sort_idx);

    % Add boundary values: Aw=0 (and P=0) at exact hull limits
    z_all_aw  = [hull_z_min; z_all_aw; hull_z_max];
    Aw_all    = [0; Aw_all;  0];
    Ixx_all   = [0; Ixx_all; 0];
    Iyy_all   = [0; Iyy_all; 0];
    P_all     = [0; P_all;   0];

    config.Aw_table_z     = z_all_aw;
    config.Aw_table       = Aw_all;
    config.I_wp_xx_table  = Ixx_all;
    config.I_wp_yy_table  = Iyy_all;
    config.P_table        = P_all;

    fprintf('    Aw table: %d points (%d coarse + %d refined)\n', ...
            length(z_all_aw), n_coarse, length(z_refine));
    fprintf('    Aw range: [%.4f, %.4f] m²\n', min(Aw_all), max(Aw_all));

    %% §2c  PRECOMPUTE V_sub AND CB_z TABLES  ─────────────────────────
    %  V_sub and CB_z from cumulative Aw integration (v9.0).
    %
    %  PREVIOUS METHOD (v8.x):
    %    Called compute_submerged (divergence theorem on partial parametric
    %    surfaces) at each z-level.  This produces 13–33% errors for
    %    offset-axis hull families because GL quadrature cannot resolve
    %    the step discontinuity at the waterline crossing.
    %
    %  NEW METHOD (v9.0):
    %    V_sub(z_wl) = ∫_{z_min}^{z_wl} A(z) dz         (trapz)
    %    CB_z(z_wl)  = ∫_{z_min}^{z_wl} z·A(z) dz / V_sub
    %
    %    Integrates the smooth, validated Aw(z) table.  No surface
    %    integrals, no orientation checks, no cancellation issues.
    %    Accuracy limited only by Aw table resolution (controlled by
    %    params.aw_table_dz).

    fprintf('  Precomputing V_sub / CB_z tables (%d z-levels, Aw-trapz)...\n', ...
            length(z_all_aw));

    n_vtab      = length(z_all_aw);
    V_sub_table = zeros(n_vtab, 1);
    CB_z_table  = zeros(n_vtab, 1);

    for k_tab = 1:n_vtab
        z_wl_k = z_all_aw(k_tab);

        % Cumulative integral of A(z) from hull bottom to z_wl_k
        mask = z_all_aw <= z_wl_k;
        z_int = z_all_aw(mask);
        A_int = Aw_all(mask);

        if length(z_int) < 2
            V_sub_table(k_tab) = 0;
            CB_z_table(k_tab)  = hull_z_min;
            continue;
        end

        V_sub_table(k_tab) = trapz(z_int, A_int);

        if V_sub_table(k_tab) > 1e-10
            CB_z_table(k_tab) = trapz(z_int, z_int .* A_int) / V_sub_table(k_tab);
        else
            V_sub_table(k_tab) = 0;
            CB_z_table(k_tab)  = hull_z_min;
        end
    end

    config.V_sub_table = V_sub_table;
    config.CB_z_table  = CB_z_table;

    fprintf('    V_sub range: [%.4f, %.4f] m^3\n', min(V_sub_table), max(V_sub_table));
    fprintf('    CB_z  range: [%.4f, %.4f] m\n',   min(CB_z_table),  max(CB_z_table));

    % §2c-override: Replace the divergence-theorem hull volume with the
    % Aw-trapz maximum.  The DT value (from hp_full.volume) overcounts
    % for hulls whose visible surfaces do NOT form a fully closed
    % boundary.  E1.ms2 Ellipsoid is an open bowl — the DT includes
    % a virtual cap volume (~4 m³) that does not physically exist.
    %
    % V_sub_max = integral Aw(z) dz is always correct: it counts only
    % the cross-sectional area present at each z-level.  For a fully
    % closed hull both values agree; for open hulls the Aw-trapz value
    % is the physically meaningful displaced volume.
    %
    % The centroid (config.hull_centroid) and inertia integrals from
    % the DT are retained — they remain accurate for open hulls.
    V_sub_max_Aw = max(V_sub_table);
    if abs(V_sub_max_Aw - config.total_wec_volume) / ...
            max(config.total_wec_volume, 1e-6) > 0.05
        fprintf('    NOTE: DT volume (%.4f m3) differs from Aw-trapz (%.4f m3) by %.1f%%.\n', ...
                config.total_wec_volume, V_sub_max_Aw, ...
                abs(V_sub_max_Aw - config.total_wec_volume) / ...
                config.total_wec_volume * 100);
        fprintf('    Using Aw-trapz value (open-hull / non-closed surface).\n');
    end
    config.total_wec_volume = V_sub_max_Aw;

    % ── §2d  PRECOMPUTE WETTED SURFACE AREA TABLE ─────────────────────
    %  S_wet(z_wl) = ∬_{z ≤ z_wl} ||Su×Sv|| du dv  (hull sides only)
    %
    %  Built on the same z-grid as Aw / V_sub so that
    %  calculate_3d_properties can interpolate with interp1 on Aw_table_z.
    %
    %  IMPORTANT: the augmentation block below (strip edges → 81 pts) also
    %  extends S_wet_table, keeping it the same length as Aw_table_z.
    %
    %  Definition: wetted hull area EXCLUDING the waterplane.
    %  Standard use: frictional drag, viscous BEM corrections, Re-scaling.
    %
    %  Guard: only if config.ms2_model is present (WAMIT-only PATH B/C has
    %  no geometry). Absence sets config.S_wet_table = [].

    if isfield(config, 'ms2_model') && ~isempty(config.ms2_model)
        fprintf('  Precomputing S_wet table (%d z-levels, GL n=20)...\n', n_vtab);
        S_wet_table = zeros(n_vtab, 1);
        for k_tab = 1:n_vtab
            z_wl_k = z_all_aw(k_tab);
            if z_wl_k <= hull_z_min
                S_wet_table(k_tab) = 0;
            else
                S_wet_table(k_tab) = WEC_HydroProperties.compute_wetted_surface_area( ...
                    config.ms2_model, z_wl_k);
            end
        end
        config.S_wet_table = S_wet_table;
        fprintf('    S_wet range: [%.4f, %.4f] m^2\n', ...
                min(S_wet_table), max(S_wet_table));
    else
        config.S_wet_table = [];   % WAMIT-only or geometry-only: no MS2 geometry
        fprintf('  S_wet table: skipped (no ms2_model)\n');
    end

    % ── Validation: monotonicity and bounds check ──
    dV = diff(V_sub_table);
    n_violations = sum(dV < -1e-6);
    if n_violations > 0
        warning('WEC:VsubTableNonMono', ...
                'V_sub table has %d non-monotonic intervals (max dV = %.4e m^3).', ...
                n_violations, min(dV));
    end
    if min(CB_z_table) < hull_z_min - 0.1 || max(CB_z_table) > hull_z_max + 0.1
        warning('WEC:CBzOutOfBounds', ...
                'CB_z table [%.3f, %.3f] exceeds hull bounds [%.3f, %.3f].', ...
                min(CB_z_table), max(CB_z_table), hull_z_min, hull_z_max);
    end

    % Topology struct — repurposed for strip geometry data
    config.topology = struct();

    %% §3  BALLAST CONFIGURATION  ─────────────────────────────────────
    %  Density nodes: one per horizontal strip.
    %
    %  DEFAULT (constructability disabled):
    %    N nodes placed uniformly from hull bottom to hull top.
    %
    %  CONSTRUCTABILITY ENABLED (wall_position = 'top', default):
    %    Strip N = wall (Z_max − wall_height to Z_max), pinned at rho_UHPC.
    %    Strips 1..N-1 = platform (Z_min to wall_z_boundary), uniform spacing.
    %
    %  CONSTRUCTABILITY ENABLED (wall_position = 'bottom', capsized _180):
    %    Strip 1 = wall (Z_min to Z_min + wall_height), pinned at rho_UHPC.
    %    Strips 2..N = platform (wall_z_boundary to Z_max), uniform spacing.
    %
    %  Strip boundaries are defined FIRST, then node positions are derived
    %  as strip centroids.  This guarantees the wall boundary lands exactly
    %  at wall_z_boundary — the midpoint formula would not.
    %
    %  WHY store strip_edges in config?
    %    When constructability is enabled, downstream code
    %    (WEC_Constructable_Hull.realize) needs the exact strip boundaries
    %    to define the wall/platform split.  Recomputing from node midpoints
    %    would shift the wall boundary by half a strip.

    config.num_ballast_sections   = params.num_ballast_sections;
    config.ballast_density_bounds = params.ballast_density_bounds;
    config.max_density_ratio      = params.max_density_ratio;

    hull_z_min = hull_extents_z(1);
    hull_z_max = hull_extents_z(2);
    config.hull_z_min = hull_z_min;
    config.hull_z_max = hull_z_max;

    if params.enable_constructability
        % ── Non-uniform strip layout: N-1 platform + 1 wall ─────
        %
        %  Wall position:
        %    'top'    — wall is the LAST strip:  [Z_max − wall_h, Z_max]
        %               Platform strips fill [Z_min, Z_max − wall_h].
        %               Default for upright hulls (e.g. C0.ms2).
        %
        %    'bottom' — wall is the FIRST strip: [Z_min, Z_min + wall_h]
        %               Platform strips fill [Z_min + wall_h, Z_max].
        %               For capsized hulls (e.g. C0_180.ms2).
        %
        %  In both cases: N+1 edges → N strips, wall pinned at rho_UHPC.

        N = config.num_ballast_sections;

        % Default to 'top' for backward compatibility
        if isfield(params, 'wall_position')
            wall_pos = params.wall_position;
        else
            wall_pos = 'top';
        end
        config.wall_position = wall_pos;

        wall_h = params.constructability_wall_height;

        switch wall_pos
            case 'top'
                % ── Wall at top (upright hull) ────────────────────
                wall_z_boundary = hull_z_max - wall_h;

                assert(wall_z_boundary > hull_z_min + 0.01, ...
                       'WEC:WallTooTall', ...
                       'wall_height (%.2f m) leaves no platform region. Hull spans [%.2f, %.2f] m.', ...
                       wall_h, hull_z_min, hull_z_max);

                % N-1 platform strips below, 1 wall strip on top
                platform_edges = linspace(hull_z_min, wall_z_boundary, N)';
                config.strip_edges = [platform_edges; hull_z_max];

                config.wall_strip_index = N;   % last strip = wall

                fprintf('  Ballast: %d strips (%d platform + 1 wall at TOP)\n', N, N - 1);
                fprintf('    Platform: [%.2f, %.2f] m  (%d strips, uniform)\n', ...
                        hull_z_min, wall_z_boundary, N - 1);
                fprintf('    Wall:     [%.2f, %.2f] m  (1 strip, solid UHPC)\n', ...
                        wall_z_boundary, hull_z_max);

            case 'bottom'
                % ── Wall at bottom (capsized hull) ────────────────
                wall_z_boundary = hull_z_min + wall_h;

                assert(wall_z_boundary < hull_z_max - 0.01, ...
                       'WEC:WallTooTall', ...
                       'wall_height (%.2f m) leaves no platform region. Hull spans [%.2f, %.2f] m.', ...
                       wall_h, hull_z_min, hull_z_max);

                % 1 wall strip at bottom, N-1 platform strips above
                platform_edges = linspace(wall_z_boundary, hull_z_max, N)';
                config.strip_edges = [hull_z_min; platform_edges];

                config.wall_strip_index = 1;   % first strip = wall

                fprintf('  Ballast: %d strips (1 wall at BOTTOM + %d platform)\n', N, N - 1);
                fprintf('    Wall:     [%.2f, %.2f] m  (1 strip, solid UHPC)\n', ...
                        hull_z_min, wall_z_boundary);
                fprintf('    Platform: [%.2f, %.2f] m  (%d strips, uniform)\n', ...
                        wall_z_boundary, hull_z_max, N - 1);

            otherwise
                error('WEC:InvalidWallPosition', ...
                      'wall_position must be ''top'' or ''bottom'', got ''%s''.', wall_pos);
        end

        % Node positions = strip centroids (same for both cases)
        config.density_nodes_z = zeros(N, 1);
        for i = 1:N
            config.density_nodes_z(i) = 0.5 * (config.strip_edges(i) + config.strip_edges(i+1));
        end

        fprintf('    Density bounds: [%.0f, %.0f] kg/m^3 (platform only)\n', ...
                config.ballast_density_bounds(1), config.ballast_density_bounds(2));
    else
        % ── Default: uniform layout ──────────────────────────────
        config.density_nodes_z = linspace(hull_z_min, hull_z_max, ...
                                          config.num_ballast_sections)';
        config.wall_strip_index = [];
        config.strip_edges      = [];   % recompute from midpoints downstream

        fprintf('  Ballast: %d nodes over [%.2f, %.2f] m, rho in [%.0f, %.0f] kg/m^3\n', ...
                config.num_ballast_sections, hull_z_min, hull_z_max, ...
                config.ballast_density_bounds(1), config.ballast_density_bounds(2));
    end

    %% §3b  STEEL-FILL SOLVER (post-optimisation parameters)  ───────
    %  Forward parameters consumed by WEC_Shell_Offset.solve, which runs
    %  ONCE after Stage-2 fmincon converges (see WEC_Main_Optimizer §4c).
    %  Nothing is computed here — the inverse solve needs final_props.
    %
    %  config.shell is set to [] for downstream backward compatibility:
    %  any code that still gates on `~isempty(config.shell)` will take
    %  its else branch (the homogeneous-strip path), which is correct
    %  because the optimiser sees a homogeneous hull (no shell mass) and
    %  the steel jacket / fill is realised purely as a post-step.

    config.enable_steel_solve      = params.enable_steel_solve;
    config.rho_steel               = params.rho_steel;
    config.rho_air                 = params.rho_air;
    config.steel_t_init            = params.steel_t_init;
    config.steel_t_min             = params.steel_t_min;
    config.steel_max_slope_factor  = params.steel_max_slope_factor;
    config.steel_n_z_grid          = params.steel_n_z_grid;
    config.shell                   = [];   % retired forward shell — see comment above

    % HAMS period-grid forwards (consumed by HAMS_Pipeline.default_hams_params
    % via the config-override path).  Optional — defaults apply if absent.
    if isfield(params, 'hams_T_min')
        config.hams_T_min  = params.hams_T_min;
    end
    if isfield(params, 'hams_T_max')
        config.hams_T_max  = params.hams_T_max;
    end
    if isfield(params, 'hams_T_step')
        config.hams_T_step = params.hams_T_step;
    end

    if config.enable_steel_solve
        fprintf('  Steel-fill solver: ENABLED (post-optimisation, see Stage 4)\n');
        fprintf('    rho_steel = %.0f kg/m^3   rho_air = %.0f kg/m^3\n', ...
                config.rho_steel, config.rho_air);
        fprintf('    t_steel initial guess: %.4f m (%.2f in)\n', ...
                config.steel_t_init, config.steel_t_init / 0.0254);
        fprintf('    t_steel min (fabrication floor): %.5f m (%.2f in)\n', ...
                config.steel_t_min, config.steel_t_min / 0.0254);
    else
        fprintf('  Steel-fill solver: DISABLED\n');
    end

    %% §3c  CONSTRUCTABILITY POST-PROCESSING  ─────────────────────────
    %  Forward all constructability parameters to config.
    %  When disabled, all fields exist but are ignored downstream.

    config.enable_constructability       = params.enable_constructability;
    config.constructability_rho_hull     = params.constructability_rho_hull;
    config.constructability_rho_fill     = params.constructability_rho_fill;
    config.constructability_t_min        = params.constructability_t_min;
    config.constructability_wall_height  = params.constructability_wall_height;
    config.constructability_n_sub        = params.constructability_n_sub;

    % UHPC global-solve tuning (forwarded only when present in params)
    if isfield(params, 'uhpc_t_init')
        config.uhpc_t_init = params.uhpc_t_init;
    end
    if isfield(params, 'uhpc_max_slope_factor')
        config.uhpc_max_slope_factor = params.uhpc_max_slope_factor;
    end
    if isfield(params, 'uhpc_n_z_grid')
        config.uhpc_n_z_grid = params.uhpc_n_z_grid;
    end

    if config.enable_constructability
        fprintf('  Constructability post-processing: ENABLED\n');
        fprintf('    Hull material:  %.0f kg/m³ (UHPC)\n', config.constructability_rho_hull);
        fprintf('    Void fill:      %.1f kg/m³ (air)\n', config.constructability_rho_fill);
        fprintf('    Min thickness:  %.4f m (%.1f in)\n', ...
                config.constructability_t_min, config.constructability_t_min / 0.0254);
        fprintf('    Wall height:    %.2f m\n', config.constructability_wall_height);
    else
        fprintf('  Constructability post-processing: DISABLED\n');
    end

    %% §3d  CONSTRUCTABILITY — PER-STRIP DENSITY BOUNDS  ──────────────
    %
    %  MODEL (Phase 2 / Increment 2 — uniform-thickness perpendicular
    %  offset shell, replaces the prior homothetic-scaling annular void)
    %
    %    Each platform strip is realised as a UHPC shell of uniform
    %    perpendicular thickness t (≥ t_min for structural feasibility),
    %    optionally with a UHPC fill in the cavity.  The minimum strip
    %    mass — and hence the minimum strip effective density — occurs at
    %    t = t_min with no fill:
    %
    %      V_shell(strip, t_min) = π · ∫_{z_lo}^{z_hi} (r² − r_inner²) dz
    %      r_inner(z)            = max(0, r(z) − t_min · √(1 + r'(z)²))
    %
    %    ρ_min_strip = (V_shell · ρ_hull + (V_strip − V_shell) · ρ_fill) / V_strip
    %
    %    The cap r_inner = 0 wherever t_min · √(1+r'²) ≥ r captures the
    %    fact that very steep or very narrow z-sections cannot support a
    %    uniform-thickness inner offset and degenerate to solid UHPC.
    %    See WEC_HydroProperties.compute_perpendicular_shell_volume.
    %
    %  WHY THIS MODEL (vs the prior homothetic ρ_min)
    %    The homothetic model scales the entire inner contour by a single
    %    s ∈ [0, 1].  When the strip spans a region where r changes
    %    rapidly (the bulb-to-neck shoulder), s gets capped by the worst
    %    z-sample, forcing a thick UHPC puck across the WHOLE strip even
    %    where the geometry has plenty of room for a thin shell.  The
    %    offset-shell model has uniform PERPENDICULAR thickness and
    %    handles the shoulder pointwise — solid where r' is huge, thin
    %    elsewhere.  ρ_min(strip) drops by 5×–10× on shoulder strips.
    %
    %  CONTRACT
    %    This bound MUST match what WEC_Constructable_Hull.realize will
    %    physically produce in Increment 3 (where the realiser is
    %    rewritten to use the same offset-shell model).  Until then, the
    %    OLD homothetic realiser will cap densities below the perpendicular
    %    s_max bound (still computed below for that legacy path) — so
    %    densities the optimiser picks below the perpendicular ρ_min will
    %    be silently pulled up by the realiser, producing measurable mass
    %    drift in the verification table.  That drift IS the test signal
    %    we want, and it disappears when the realiser is replaced.
    %
    %  PARAMETER CHAIN
    %    wall_height → strip_edges → {r(z), r'(z)} → V_shell(t_min)
    %         ↑           ↑                ↑              ↑
    %      user set    config           geometry     t_min (user set)
    %                                                    ↓
    %                                       ρ_min(strip) → lb(i)
    %
    %  IMPLEMENTATION
    %    For each platform strip, sample z-levels at the realiser's
    %    cadence (config.constructability_n_sub points), call
    %    compute_rmin_at_z to build r(z), then call
    %    compute_perpendicular_shell_volume which integrates V_shell at
    %    t_min via central differences for r' and the dome-tip cap.
    %    The perpendicular s_max is also retained for backward compat
    %    with the OLD homothetic realiser (will be removed in Increment 3).

    if config.enable_constructability
        fprintf('  Computing per-strip density bounds from hull geometry...\n');

        N            = config.num_ballast_sections;
        rho_hull_c   = config.constructability_rho_hull;
        rho_fill_c   = config.constructability_rho_fill;
        t_min_c      = config.constructability_t_min;
        w_idx        = config.wall_strip_index;
        % CONTRACT: mirror the realiser's z-sampling density.  The realiser
        % uses config.constructability_n_sub (default 100) per strip — if
        % we sample more sparsely the realiser will find a tighter s_max
        % that the optimiser bound never saw, re-opening the relaxed-vs-
        % true feasibility-set gap this whole subsystem closes.
        n_rmin_sub   = config.constructability_n_sub;

        per_strip_lb         = ones(N, 1) * config.ballast_density_bounds(1);
        per_strip_rmin       = zeros(N, 1);
        per_strip_smax       = zeros(N, 1);
        per_strip_m_min      = zeros(N, 1);
        per_strip_V_shell    = zeros(N, 1);  % at t_min, perpendicular offset
        per_strip_A_outer    = zeros(N, 1);  % lateral surface area

        for i = 1:N
            if i == w_idx
                % Wall strip: pinned to rho_hull, no geometry scan needed
                per_strip_lb(i)         = rho_hull_c;
                per_strip_rmin(i)       = Inf;
                per_strip_smax(i)       = 0;
                per_strip_V_shell(i)    = NaN;  % wall is solid — no shell concept
                per_strip_A_outer(i)    = NaN;
                continue;
            end

            z_lo = config.strip_edges(i);
            z_hi = config.strip_edges(i + 1);

            % Mirror the realiser's z-sampling EXACTLY
            % (WEC_Constructable_Hull.realize lines 335,349-350): include
            % both endpoints, with the top sample nudged 1mm below z_hi to
            % avoid the wall-strip surface that starts exactly at z_hi.
            % Excluding endpoints (the prior approach) skipped tight
            % shoulders at the strip boundary and made the central-
            % difference r' use different neighbours than the realiser.
            n_rmin_sub_i = max(n_rmin_sub, ceil((z_hi - z_lo) / 0.01));
            z_samples = linspace(z_lo, z_hi, n_rmin_sub_i);
            z_samples(end) = z_hi - 1e-3;

            % ── Pass 1: collect r_min_k and A_k at each z-sample ────
            %  We keep the global r_min purely as a diagnostic output
            %  (per_strip_rmin).  s_max is now computed via the
            %  PERPENDICULAR formula below — strictly tighter than the
            %  prior radial 1−t_min/r which ignored the slope r'.
            n_k     = length(z_samples);
            r_min_k_arr  = zeros(n_k, 1);
            A_k_arr      = zeros(n_k, 1);
            r_min_i      = Inf;

            for k = 1:n_k
                [r_k, ~] = WEC_HydroProperties.compute_rmin_at_z( ...
                        config.ms2_model, z_samples(k), 100, ...
                        config.boundary_cache);

                if r_k > 1e-10 && r_k < r_min_i
                    r_min_i = r_k;
                end

                r_min_k_arr(k) = r_k;
                A_k_arr(k) = max(0, interp1(config.Aw_table_z, ...
                                 config.Aw_table, z_samples(k), ...
                                 'linear', 0));
            end

            per_strip_rmin(i) = r_min_i;

            % ── Pass 2: profile slope rp(z) via central differences ──
            %  Mirrors WEC_Constructable_Hull.realize lines 409-421.
            %  Required by the perpendicular wall-thickness formula below.
            rp_k_arr = zeros(n_k, 1);
            if n_k >= 2
                for k_rp = 1:n_k
                    if k_rp == 1
                        dz_rp = z_samples(2) - z_samples(1);
                        rp_k_arr(k_rp) = (r_min_k_arr(2) - r_min_k_arr(1)) / dz_rp;
                    elseif k_rp == n_k
                        dz_rp = z_samples(n_k) - z_samples(n_k-1);
                        rp_k_arr(k_rp) = (r_min_k_arr(n_k) - r_min_k_arr(n_k-1)) / dz_rp;
                    else
                        dz_rp = z_samples(k_rp+1) - z_samples(k_rp-1);
                        rp_k_arr(k_rp) = (r_min_k_arr(k_rp+1) - r_min_k_arr(k_rp-1)) / dz_rp;
                    end
                end
            end
            % Dome-tip correction: zero rp at samples where r is below t_min
            % (the surface has effectively reached its closing point and the
            % one-sided finite difference produces a spurious large |r'|).
            % Mirrors WEC_Constructable_Hull.realize line 436.
            rp_k_arr(r_min_k_arr < t_min_c) = 0;

            % ── Pass 3: per-z perpendicular s_max → strip s_max → ρ_min ──
            %  PERPENDICULAR wall thickness on a revolution profile:
            %      t_perp(z) = (1 − s)·r·sqrt(1 + r'²) / (1 + s·r'²)
            %  Setting t_perp = t_min and solving for s:
            %      s_max(z) = (r·L − t_min) / (r·L + t_min·r'²),  L = √(1+r'²)
            %  Matches WEC_Constructable_Hull.realize lines 455-481 EXACTLY.
            %  At r' = 0 this reduces to the old radial 1 − t_min/r.
            %
            %  The realiser uses a UNIFORM s_i across the whole strip, capped
            %  at the GLOBAL MIN of s_max(z) over z-samples.  So the strip's
            %  minimum achievable ρ_eff is:
            %      ρ_min_strip = ρ_hull − s_max_strip² · (ρ_hull − ρ_fill)
            %  (NOT the volume-weighted local mean — that was wrong.  The
            %  realiser's actual mass formula is uniform: m_strip = V·ρ_eff,
            %  see WEC_Constructable_Hull.realize lines 506-510.)
            % Degenerate-strip guard: if no valid r samples exist (the
            % B-spline evaluator found no surface or the strip closes
            % completely), force solid (s_max=0, rho=rho_hull) to mirror
            % WEC_Constructable_Hull.realize:451-453.  Without this guard
            % the loop below leaves s_max_strip = 1 and the strip's lb
            % collapses to rho_fill, the OPPOSITE of the realiser's
            % behaviour.
            valid_mask = (r_min_k_arr > 1e-10) & isfinite(r_min_k_arr);
            if ~any(valid_mask)
                s_max_strip = 0;
            else
                s_max_strip = 1.0;
                for k = 1:n_k
                    r_k  = r_min_k_arr(k);
                    rp_k = rp_k_arr(k);

                    if r_k < 1e-10 || ~isfinite(r_k)
                        continue;   % degenerate sample — does not constrain s
                    end

                    L_k = r_k * sqrt(1 + rp_k^2);
                    if L_k <= t_min_c
                        s_max_k = 0;
                    else
                        s_max_k = (L_k - t_min_c) / (L_k + t_min_c * rp_k^2);
                    end
                    s_max_k = max(0, min(1, s_max_k));

                    s_max_strip = min(s_max_strip, s_max_k);
                end
            end

            per_strip_smax(i) = s_max_strip;

            % ── Increment-2 SWITCH: offset-shell ρ_min ──────────────
            % Replace the homothetic perpendicular-s_max formula with a
            % uniform-thickness offset-shell formula (mirrors the
            % WEC_Shell_Offset realisation model).  Strip mass at ρ_min:
            %   V_shell·ρ_hull  +  V_int·ρ_fill
            % where V_shell is the perpendicular-offset shell volume at
            % t = t_min, V_int = V_strip − V_shell, and r_inner is capped
            % at zero where t·sqrt(1+r'²) ≥ r (locally solid).
            %
            % This is strictly looser than the homothetic ρ_min for any
            % strip whose worst z forces s_max(strip) → 0 — the dominant
            % bottleneck on hulls with sharp shoulders.  See
            % WEC_HydroProperties.compute_perpendicular_shell_volume for
            % the derivation and accuracy claim.
            V_strip_i_est = trapz(z_samples(:), A_k_arr(:));
            try
                [V_shell_at_tmin, A_outer_i, ~] = ...
                    WEC_HydroProperties.compute_perpendicular_shell_volume( ...
                        config.ms2_model, z_lo, z_hi, t_min_c, ...
                        n_rmin_sub_i, config.boundary_cache);
            catch ME
                warning('WEC_Configuration_Builder:OffsetShellFailed', ...
                    'Strip %d: compute_perpendicular_shell_volume failed (%s). Falling back to perpendicular-formula ρ_min.', ...
                    i, ME.message);
                V_shell_at_tmin = V_strip_i_est * (rho_hull_c - rho_fill_c) / rho_hull_c;  % ρ_min ≈ perpendicular fallback
                A_outer_i = NaN;
            end

            per_strip_V_shell(i) = V_shell_at_tmin;
            per_strip_A_outer(i) = A_outer_i;

            if V_shell_at_tmin >= V_strip_i_est
                rho_min_offset = rho_hull_c;     % strip too narrow → solid
            else
                V_int_i        = V_strip_i_est - V_shell_at_tmin;
                rho_min_offset = (V_shell_at_tmin * rho_hull_c + V_int_i * rho_fill_c) / V_strip_i_est;
            end
            rho_min_offset = max(rho_fill_c, min(rho_hull_c, rho_min_offset));

            per_strip_lb(i) = max(config.ballast_density_bounds(1), rho_min_offset);

            % Diagnostic: minimum strip mass at this lower bound
            per_strip_m_min(i) = per_strip_lb(i) * V_strip_i_est;
        end

        config.per_strip_density_lb       = per_strip_lb(:)';
        config.per_strip_rmin             = per_strip_rmin(:)';
        config.per_strip_smax             = per_strip_smax(:)';
        config.per_strip_m_min            = per_strip_m_min(:)';
        config.per_strip_V_shell_at_tmin  = per_strip_V_shell(:)';
        config.per_strip_A_outer          = per_strip_A_outer(:)';

        fprintf('    Strip   z_lo     z_hi    A_outer[m²]  V_shell[m³]  rho_min[kg/m3]   (offset-shell @ t_min)\n');
        fprintf('    %s\n', repmat('-', 1, 80));
        for i = 1:N
            if i == w_idx
                fprintf('    %-5d  %+6.3f   %+6.3f   WALL (pinned to %.0f kg/m³)\n', ...
                        i, config.strip_edges(i), config.strip_edges(i+1), rho_hull_c);
            else
                fprintf('    %-5d  %+6.3f   %+6.3f   %10.4f   %10.4f   %10.1f\n', ...
                        i, config.strip_edges(i), config.strip_edges(i+1), ...
                        per_strip_A_outer(i), per_strip_V_shell(i), per_strip_lb(i));
            end
        end
    else
        config.per_strip_density_lb       = [];
        config.per_strip_rmin             = [];
        config.per_strip_smax             = [];
        config.per_strip_m_min            = [];
        config.per_strip_V_shell_at_tmin  = [];
        config.per_strip_A_outer          = [];
    end

    %% §3e  MINIMUM ACHIEVABLE MASS (constructability mode)  ──────────
    %
    %  PURPOSE
    %    Compute the minimum total hull mass when every strip is at its
    %    constructability-enforced lower density bound.  This is used:
    %      (a) as a feasibility check (can the hull float at all?), and
    %      (b) as an inequality constraint in both Stage 1 and Stage 2
    %          optimizers to prevent exploration of infeasible density
    %          regions.
    %
    %  COMPUTATION
    %    m_min = Σ_i  V_strip_i × per_strip_density_lb(i)
    %
    %    Each V_strip_i is computed from the 3D parametric geometry via
    %    WEC_HydroProperties.compute_strip (divergence theorem), not from
    %    the 2D surrogate.  This ensures the constraint is accurate
    %    regardless of k_vol.
    %
    %  FEASIBILITY CHECK
    %    If m_min > ρ_water × V_total, the body cannot float even when
    %    fully submerged.  The optimization problem is infeasible.
    %    This is detected here and reported as an error BEFORE entering
    %    the optimizer loop.

    if config.enable_constructability && ~isempty(config.strip_edges) ...
            && ~isempty(config.per_strip_density_lb)

        fprintf('  Computing minimum achievable mass (V_strip × per_strip_density_lb, offset-shell)...\n');
        N_strips = config.num_ballast_sections;
        m_min = 0;
        m_max = 0;

        for i = 1:N_strips
            z_lo_i = config.strip_edges(i);
            z_hi_i = config.strip_edges(i + 1);

            if z_hi_i <= z_lo_i + 1e-10
                continue;
            end

            % Strip volume via divergence theorem (canonical, matches
            % §3f and the realiser's compute_strip call).
            strip_i = WEC_HydroProperties.compute_strip( ...
                config.ms2_model, z_lo_i, z_hi_i, ...
                struct('n_quad', 16, ...
                       'Aw_table_z', config.Aw_table_z, ...
                       'Aw_table', config.Aw_table));
            V_strip_i = strip_i.V;

            % Offset-shell ρ_min represents uniform-s mass density across
            % the strip, so m_strip_min = V_strip · ρ_min(strip).  Drop the
            % §3d trapz-based per_strip_m_min fallback — V_strip from
            % divergence theorem is the canonical strip volume.
            m_min = m_min + V_strip_i * config.per_strip_density_lb(i);
            m_max = m_max + V_strip_i * config.ballast_density_bounds(2);
        end

        config.m_min_constructability = m_min;
        config.m_max_constructability = m_max;

        % Buoyancy limits: maximum buoyancy = fully submerged hull
        max_buoyancy = config.RHO_WATER * config.total_wec_volume;

        fprintf('    Minimum achievable mass: %.1f kg\n', m_min);
        fprintf('    Maximum achievable mass: %.1f kg\n', m_max);
        fprintf('    Maximum buoyancy (fully submerged): %.1f kg\n', max_buoyancy);

        % Feasibility check: can the hull float at all?
        if m_min > max_buoyancy
            error('WEC:Infeasible', ...
                ['Constructability minimum mass (%.1f kg) exceeds maximum ' ...
                 'buoyancy (%.1f kg = rho_water × V_total).\n' ...
                 '  The hull cannot float with these constructability ' ...
                 'constraints.\n' ...
                 '  Reduce wall_height (currently %.2f m), reduce t_min ' ...
                 '(currently %.4f m),\n' ...
                 '  or increase hull volume (currently %.4f m^3).'], ...
                m_min, max_buoyancy, ...
                config.constructability_wall_height, ...
                config.constructability_t_min, ...
                config.total_wec_volume);
        end

        % Tight feasibility warning
        buoyancy_margin = (max_buoyancy - m_min) / max_buoyancy * 100;
        if buoyancy_margin < 20
            warning('WEC:TightFeasibility', ...
                ['Buoyancy margin is only %.1f%%. ' ...
                 'Convergence may be difficult. ' ...
                 'Consider relaxing constructability constraints.'], ...
                buoyancy_margin);
        end

        fprintf('    Buoyancy margin: %.1f%%\n', buoyancy_margin);
    else
        config.m_min_constructability = 0;
        config.m_max_constructability = inf;
    end

    %% §3f  PRECOMPUTE STRIP GEOMETRY (v4.3)  ─────────────────────────
    %
    %  Strip volumes, centroids, and inertias depend only on the hull
    %  geometry and strip_edges — NOT on density or draft.  Computing
    %  them once here eliminates all property evaluation calls from
    %  the optimisation loop.
    %
    %  Uses WEC_HydroProperties.compute_strip (parametric divergence
    %  theorem with source/mirror dedup) — the same method as
    %  compute_submerged, applied to horizontal strips.
    %
    %  calculate_3d_properties uses these as:
    %    mass_i       = config.strip_V(i) × ρ_i
    %    cg_z_num    += mass_i × config.strip_CB_z(i)
    %    Iyy_total   += ρ_i × config.strip_Iyy(i)

    N = config.num_ballast_sections;

    % Determine strip edges (for default mode, constructability already set)
    if ~config.enable_constructability || isempty(config.strip_edges)
        if N > 1
            dz_se = config.density_nodes_z(2) - config.density_nodes_z(1);
        else
            dz_se = hull_z_max - hull_z_min;
        end
        se = zeros(N + 1, 1);
        se(1) = config.density_nodes_z(1) - dz_se/2;
        for i = 2:N
            se(i) = 0.5 * (config.density_nodes_z(i-1) + config.density_nodes_z(i));
        end
        se(N+1) = config.density_nodes_z(N) + dz_se/2;
        se(1)   = max(se(1), hull_z_min);
        se(end) = min(se(end), hull_z_max);
        config.strip_edges = se;
    end

    %% §3f-pre  AUGMENT Aw TABLE WITH EXACT STRIP BOUNDARY EVALUATIONS
    %
    %  The Aw_table was built on an adaptive grid of 45 z-levels.  Strip
    %  boundary z-values (e.g. z = −0.9 m, the platform-to-column shoulder)
    %  may fall BETWEEN table grid points.  The compute_strip fast path uses
    %  interp1 to get A(z_lo) and A(z_hi), but interpolation is inaccurate
    %  near rapid A(z) transitions (Aw jumps from 0.6 to 10.4 m² over 0.34 m
    %  at the shoulder).
    %
    %  FIX: evaluate A(z) exactly at every strip edge using the parametric
    %  boundary definition (extract_isocurve_at_z), then insert those points
    %  into the Aw_table grid so interp1 finds an exact match.
    %
    %  V_sub_table is on the same z-grid (config.Aw_table_z); augmenting both
    %  together keeps them consistent for calculate_3d_properties, which uses
    %    interp1(config.Aw_table_z, config.V_sub_table, z_wl)
    %  V_sub at the new strip-edge z-values is computed with compute_submerged.

    n_edge_added  = 0;

    for ie = 1:length(config.strip_edges)
        z_edge = config.strip_edges(ie);

        % Skip if a table point already exists within numerical tolerance
        if any(abs(config.Aw_table_z - z_edge) < 1e-8)
            continue;
        end

        % Exact isocurve at the strip edge
        wl_edge = WEC_HydroProperties.extract_isocurve_at_z( ...
                      config.ms2_model, z_edge, n_u_cache, config.boundary_cache);
        if ~isempty(wl_edge) && size(wl_edge, 1) >= 3
            [Aw_edge, Ixx_edge, Iyy_edge] = ...
                WEC_HydroProperties.waterplane_properties(wl_edge);
        else
            Aw_edge  = 0;
            Ixx_edge = 0;
            Iyy_edge = 0;
        end

        % V_sub and CB_z at strip edge from cumulative Aw integration (v9.0)
        %   Same trapz method as §2c — no compute_submerged.
        mask_edge = config.Aw_table_z <= z_edge;
        z_int_e   = config.Aw_table_z(mask_edge);
        A_int_e   = config.Aw_table(mask_edge);

        % CRITICAL: re-sort the masked subset.  Earlier iterations of this
        % loop append new edge points to the BOTTOM of Aw_table_z (line 926
        % below) without re-sorting until line 951.  When ie>=2, the masked
        % subset (z_int_e, A_int_e) inherits that out-of-order tail, and
        % trapz on unsorted x produces wrong-sign Δx contributions —
        % corrupting Vsub_edge / CBz_edge for every edge after the first
        % and yielding a non-monotonic V_sub_table.
        [z_int_e, sort_e] = sort(z_int_e);
        A_int_e           = A_int_e(sort_e);

        if length(z_int_e) >= 2
            Vsub_edge = trapz(z_int_e, A_int_e);
            if Vsub_edge > 1e-10
                CBz_edge = trapz(z_int_e, z_int_e .* A_int_e) / Vsub_edge;
            else
                Vsub_edge = 0;
                CBz_edge  = hull_z_min;
            end
        else
            Vsub_edge = 0;
            CBz_edge  = hull_z_min;
        end

        % Insert into tables
        config.Aw_table_z    = [config.Aw_table_z;    z_edge  ];
        config.Aw_table      = [config.Aw_table;      Aw_edge ];
        config.I_wp_xx_table = [config.I_wp_xx_table; Ixx_edge];
        config.I_wp_yy_table = [config.I_wp_yy_table; Iyy_edge];
        config.V_sub_table   = [config.V_sub_table;   Vsub_edge];
        config.CB_z_table    = [config.CB_z_table;    CBz_edge ];

        % Augment S_wet_table at this strip edge to keep it the same
        % length as Aw_table_z.  Without this, interp1 in
        % calculate_3d_properties crashes with "X and V must be of the
        % same length" because Aw_table_z grows but S_wet_table does not.
        if ~isempty(config.S_wet_table)
            if z_edge <= hull_z_min
                Swet_edge = 0;
            else
                Swet_edge = WEC_HydroProperties.compute_wetted_surface_area( ...
                                config.ms2_model, z_edge);
            end
            config.S_wet_table = [config.S_wet_table; Swet_edge];
        end

        n_edge_added = n_edge_added + 1;
    end

    % Re-sort all tables by z (insertion above was unsorted)
    [config.Aw_table_z, sort_ie] = sort(config.Aw_table_z);
    config.Aw_table      = config.Aw_table(sort_ie);
    config.I_wp_xx_table = config.I_wp_xx_table(sort_ie);
    config.I_wp_yy_table = config.I_wp_yy_table(sort_ie);
    config.V_sub_table   = config.V_sub_table(sort_ie);
    config.CB_z_table    = config.CB_z_table(sort_ie);
    if ~isempty(config.S_wet_table)
        config.S_wet_table = config.S_wet_table(sort_ie);
    end

    % RE-COMPUTE V_sub_table and CB_z_table from scratch on the FINAL sorted
    % grid.  Required because the §3f-pre per-edge trapz (lines 906–923) used
    % a DIFFERENT (possibly partial / earlier-state) grid for each entry,
    % producing values that are no longer mutually consistent after sorting.
    % Using identical cumulative-trapz semantics as §2c guarantees a strictly
    % non-decreasing V_sub_table after augmentation.
    n_vtab2 = length(config.Aw_table_z);
    for k_tab = 1:n_vtab2
        z_wl_k = config.Aw_table_z(k_tab);
        mask   = config.Aw_table_z <= z_wl_k;
        z_int  = config.Aw_table_z(mask);
        A_int  = config.Aw_table(mask);

        if length(z_int) < 2
            config.V_sub_table(k_tab) = 0;
            config.CB_z_table(k_tab)  = hull_z_min;
            continue;
        end

        Vk = trapz(z_int, A_int);
        if Vk > 1e-10
            config.V_sub_table(k_tab) = Vk;
            config.CB_z_table(k_tab)  = trapz(z_int, z_int .* A_int) / Vk;
        else
            config.V_sub_table(k_tab) = 0;
            config.CB_z_table(k_tab)  = hull_z_min;
        end
    end

    % Hard guard: assert monotonicity after the rebuild.  If this trips,
    % the Aw_table itself has issues (negative entries or duplicate z's).
    dVsub = diff(config.V_sub_table);
    if any(dVsub < -1e-9)
        bad = find(dVsub < -1e-9, 1);
        error('WEC_Configuration_Builder:VsubNonMonotonic', ...
            'V_sub_table decreased by %.3em³ between z=%.4f and z=%.4f after rebuild — Aw_table likely has negative or duplicate-z entries.', ...
            -dVsub(bad), config.Aw_table_z(bad), config.Aw_table_z(bad+1));
    end

    fprintf('  Aw/V_sub tables augmented: %d strip edge(s) added → %d total points\n', ...
            n_edge_added, length(config.Aw_table_z));
    fprintf('  V_sub_table rebuilt on final grid; monotonic check passed.\n');

    fprintf('  Precomputing strip geometry (%d strips, table-based trapz)...\n', N);
    config.strip_V       = zeros(N, 1);
    config.strip_CB_z    = zeros(N, 1);
    config.strip_Iyy     = zeros(N, 1);
    config.strip_Ixx     = zeros(N, 1);
    config.strip_Izz     = zeros(N, 1);

    strip_quad_opts = struct('n_quad', 16, ...
                             'Aw_table_z', config.Aw_table_z, ...
                             'Aw_table', config.Aw_table, ...
                             'I_wp_xx_table', config.I_wp_xx_table, ...
                             'I_wp_yy_table', config.I_wp_yy_table);
    for i = 1:N
        z_lo_i = config.strip_edges(i);
        z_hi_i = config.strip_edges(i + 1);

        if z_hi_i <= z_lo_i + 1e-10
            continue;
        end

        strip_i = WEC_HydroProperties.compute_strip( ...
            config.ms2_model, z_lo_i, z_hi_i, strip_quad_opts);

        config.strip_V(i)    = strip_i.V;
        config.strip_CB_z(i) = strip_i.CB_z;
        config.strip_Iyy(i)  = strip_i.Iyy;
        config.strip_Ixx(i)  = strip_i.Ixx;
        config.strip_Izz(i)  = strip_i.Izz;

        fprintf('    Strip %2d [%+6.3f, %+6.3f]: V=%.4f m³, z̄=%.3f m\n', ...
                i, z_lo_i, z_hi_i, strip_i.V, strip_i.CB_z);
    end

    V_sum = sum(config.strip_V);
    fprintf('  Strip volume sum: %.4f m³ (hull total: %.4f m³, diff: %.2f%%)\n', ...
            V_sum, config.total_wec_volume, ...
            abs(V_sum - config.total_wec_volume) / config.total_wec_volume * 100);

    %% §4  2D SURROGATE MODEL  ────────────────────────────────────────
    %  CORRECTED EFFECTIVE WIDTH TABLE (v4.2)
    %
    %  The 2D surrogate computes strip volume as:
    %
    %    V_strip_2D(z) = profile_chord(z) × y_span(z) × dz × k_vol
    %
    %  where profile_chord(z) = max(x) − min(x) on the midplane profile.
    %
    %  PROBLEM (v4.0–v4.1):
    %    y_span was the raw transverse width: max(y) − min(y).  For hulls
    %    with non-rectangular cross-sections (e.g., C0's wide platform),
    %    chord × raw_width ≠ actual polygon area.  The per-strip volume
    %    distribution was wrong even when k_vol corrected the TOTAL volume.
    %    This caused mass balance infeasibility (exitflag = -2) because
    %    different strips got incorrect fractions of the total volume.
    %
    %  FIX (v4.2):
    %    Compute a corrected effective width so that the 2D strip area
    %    exactly matches the 3D cross-section area at every z-level:
    %
    %      y_span_corrected(z) = A_3D(z) / profile_chord(z)
    %
    %    Then:  chord(z) × y_span_corrected(z) = A_3D(z)  exactly.
    %
    %    This eliminates the need for k_vol to correct per-strip volume
    %    distribution — only a residual total-volume correction remains
    %    (from strip height discretization), which is typically < 1%.
    %
    %  CONSTRUCTION
    %    1. Evaluate 3D cross-section polygon at each z-level → A_3D(z).
    %    2. Find the 2D profile chord at each z-level → chord(z).
    %    3. y_span(z) = A_3D(z) / chord(z) when chord > 0.
    %    4. Fall back to raw Y-span when chord = 0 (degenerate).
    %    5. Fill zero-width levels with nearest non-zero neighbour.

    config.n_density_strips = params.n_density_strips;

    % Effective width from parser extents (full hull Y-span) — fallback
    y_extents = [config.ms2_model.extents(2), config.ms2_model.extents(5)];
    config.effective_width = max(y_extents) - min(y_extents);

    n_z_levels = params.n_z_levels;
    config.y_span_z_levels = linspace(hull_z_min, hull_z_max, n_z_levels)';
    config.y_span_table    = zeros(n_z_levels, 1);

    fprintf('  Building area-corrected Y-span table (%d levels, divergence theorem)...\n', n_z_levels);

    % Compute A(z) via thin-strip divergence theorem: A ≈ V_strip / dz.
    % This uses the trusted WEC_HydroProperties pipeline with correct
    % source/mirror dedup — no cross-section polygon assembly.
    h_thin = 0.5 * (hull_z_max - hull_z_min) / n_z_levels;  % half-strip thickness
    yspan_quad_opts = struct('n_quad', 12);

    for kz = 1:n_z_levels
        z_k = config.y_span_z_levels(kz);

        % --- 3D cross-section area from thin-strip divergence theorem ---
        z_lo_k = max(hull_z_min, z_k - h_thin);
        z_hi_k = min(hull_z_max, z_k + h_thin);
        dz_k   = z_hi_k - z_lo_k;

        A_3D_k = 0;
        if dz_k > 1e-10
            strip_k = WEC_HydroProperties.compute_strip( ...
                          config.ms2_model, z_lo_k, z_hi_k, yspan_quad_opts);
            if strip_k.V > 1e-10
                A_3D_k = strip_k.V / dz_k;
            end
        end

        % --- 2D profile chord at this z-level ---
        chord_k = 0;
        try
            isects = WEC_Core_Functions.find_waterline_intersections( ...
                         config.profile, z_k);
            if ~isempty(isects) && size(isects, 1) >= 2
                chord_k = max(isects(:, 1)) - min(isects(:, 1));
            end
        catch
            % Degenerate z-level — no intersection
        end

        % --- Corrected effective width ---
        if chord_k > 1e-6 && A_3D_k > 1e-10
            % Area-corrected: makes 2D strip area = 3D polygon area
            config.y_span_table(kz) = A_3D_k / chord_k;
        elseif A_3D_k > 1e-10
            % Fallback: use A as-is (degenerate profile chord)
            config.y_span_table(kz) = A_3D_k;
        end
    end

    % Fill gaps with nearest non-zero neighbour
    nonzero_idx = find(config.y_span_table > 0);
    if ~isempty(nonzero_idx)
        for kz = 1:n_z_levels
            if config.y_span_table(kz) == 0
                [~, nearest] = min(abs(nonzero_idx - kz));
                config.y_span_table(kz) = config.y_span_table(nonzero_idx(nearest));
            end
        end
    else
        config.y_span_table(:) = config.effective_width;
    end

    config.eff_w_floor = params.eff_w_floor;

    fprintf('  Y-span table: %d levels, corrected width range [%.4f, %.4f] m\n', ...
            n_z_levels, min(config.y_span_table), max(config.y_span_table));
    fprintf('  Strip count: %d\n', config.n_density_strips);

    %% §5  HYDRODYNAMIC DATA  ─────────────────────────────────────────
    %  Two paths:
    %    (A) HAMS hydro_table — precomputed by HAMS_DraftSweep.  Contains
    %        A(∞) and B_avg at N drafts, 6×6 at origin.  Extracted to 3×3
    %        [surge, heave, pitch] and transformed from origin to CG here.
    %    (B) WAMIT .1 files — legacy path.  Reads 3 pre-existing cases.
    %
    %  Both paths populate the SAME config fields:
    %    config.wamit_drafts, config.wamit_A, config.wamit_A_full, config.wamit_B_full
    %  so that ALL downstream code (interpolate_wamit_added_mass,
    %  calculate_3d_properties, coupled eigenvalue analysis) is UNCHANGED.

    use_hams_table = isfield(params, 'hydro_table_file') && ...
                     ~isempty(params.hydro_table_file) && ...
                     exist(params.hydro_table_file, 'file');

    if isfield(params, 'hydro_table_file') && ~isempty(params.hydro_table_file) ...
            && ~exist(params.hydro_table_file, 'file')
        warning('WEC:HydroTableNotFound', ...
                'hydro_table_file specified but not found: %s\nFalling back to WAMIT .1 files.', ...
                params.hydro_table_file);
    end

    if use_hams_table
        % ── PATH A: Load precomputed HAMS hydro_table ─────────────
        fprintf('  Loading HAMS hydro_table: %s\n', params.hydro_table_file);
        ht_loaded = load(params.hydro_table_file, 'hydro_table');
        config.hydro_table = ht_loaded.hydro_table;

        N_ht = length(config.hydro_table.drafts);
        config.wamit_drafts = config.hydro_table.drafts(:);
        config.wamit_z_cg   = config.hydro_table.z_cg(:);
        config.wamit_A      = zeros(N_ht, 3);
        config.wamit_B      = zeros(N_ht, 3);
        config.wamit_A_full = cell(N_ht, 1);
        config.wamit_B_full = cell(N_ht, 1);

        idx_3dof = [1, 3, 5];  % surge, heave, pitch from 6-DOF

        for i = 1:N_ht
            A_inf_6x6 = config.hydro_table.A_inf{i};   % 6×6 at origin
            B_avg_6x6 = config.hydro_table.B_avg{i};   % 6×6 at origin
            z_cg_i    = config.hydro_table.z_cg(i);     % CG_z in global frame

            % Extract 3×3 [surge, heave, pitch] submatrix
            A_3x3_origin = A_inf_6x6(idx_3dof, idx_3dof);
            B_3x3_origin = B_avg_6x6(idx_3dof, idx_3dof);

            % Transform from origin to CG (same function as legacy path)
            [A_cg, B_cg] = WEC_File_IO.transform_hydrodynamic_matrices_3x3( ...
                A_3x3_origin, B_3x3_origin, z_cg_i);

            config.wamit_A(i, :)   = [A_cg(1,1), A_cg(2,2), A_cg(3,3)];
            config.wamit_B(i, :)   = [B_cg(1,1), B_cg(2,2), B_cg(3,3)];
            config.wamit_A_full{i} = A_cg;
            config.wamit_B_full{i} = B_cg;
        end

        % Sort by draft (ascending)
        [config.wamit_drafts, si] = sort(config.wamit_drafts);
        config.wamit_z_cg   = config.wamit_z_cg(si);
        config.wamit_A      = config.wamit_A(si, :);
        config.wamit_B      = config.wamit_B(si, :);
        config.wamit_A_full = config.wamit_A_full(si);
        config.wamit_B_full = config.wamit_B_full(si);

        num_cases = N_ht;
        config.hydro_ready = true;
        fprintf('  HAMS hydro_table: %d drafts, A(inf) + B_avg at CG\n', N_ht);
        for i = 1:N_ht
            fprintf('    draft=%+.2f m, A_diag=[%.1f, %.1f, %.1f]\n', ...
                    config.wamit_drafts(i), config.wamit_A(i,:));
        end

    else
        % ── PATH B: Load WAMIT .1 files (legacy) ─────────────────
        if isfield(params, 'wamit_cases') && ~isempty(params.wamit_cases)
            wamit_cases = params.wamit_cases;
            num_cases   = size(wamit_cases, 1);

            if ~iscell(wamit_cases) || size(wamit_cases, 2) < 3
                error('WEC:BadInput', ...
                      'wamit_cases must be a cell array with columns {prefix, draft, z_cg}');
            end

            config.wamit_drafts = zeros(num_cases, 1);
            config.wamit_z_cg   = zeros(num_cases, 1);
            config.wamit_A      = zeros(num_cases, 3);
            config.wamit_B      = zeros(num_cases, 3);
            config.wamit_A_full = cell(num_cases, 1);
            config.wamit_B_full = cell(num_cases, 1);

            fprintf('  Loading %d WAMIT cases:\n', num_cases);

            for i = 1:num_cases
                prefix = wamit_cases{i, 1};
                draft  = wamit_cases{i, 2};
                z_cg   = wamit_cases{i, 3};

                filename = [prefix, '.1'];
                if ~exist(filename, 'file')
                    error('WEC:FileNotFound', 'WAMIT file not found: %s', filename);
                end

                wamit_data = WEC_File_IO.read_wamit_data( ...
                    prefix, config.RHO_WATER, config.G, config.wamit_L, z_cg);

                config.wamit_drafts(i) = draft;
                config.wamit_z_cg(i)   = z_cg;

                config.wamit_A(i, :)   = [wamit_data.A(1,1), ...
                                           wamit_data.A(2,2), ...
                                           wamit_data.A(3,3)];
                config.wamit_B(i, :)   = [wamit_data.B_rad_avg(1,1), ...
                                           wamit_data.B_rad_avg(2,2), ...
                                           wamit_data.B_rad_avg(3,3)];
                config.wamit_A_full{i} = wamit_data.A;
                config.wamit_B_full{i} = wamit_data.B_rad_avg;

                fprintf('    %s: draft = %+.2f m, z_cg = %.4f m\n', prefix, draft, z_cg);
                fprintf('      A_diag = [%.1f, %.1f, %.1f],  A13 = %.1f\n', ...
                        config.wamit_A(i,:), wamit_data.A(1,3));
            end

            [config.wamit_drafts, si] = sort(config.wamit_drafts);
            config.wamit_z_cg   = config.wamit_z_cg(si);
            config.wamit_A      = config.wamit_A(si, :);
            config.wamit_B      = config.wamit_B(si, :);
            config.wamit_A_full = config.wamit_A_full(si);
            config.wamit_B_full = config.wamit_B_full(si);

            config.hydro_ready = true;
            num_cases = length(config.wamit_drafts);
            fprintf('  WAMIT data loaded (full 3x3 at CG)\n');

        else
            % ── PATH C: No hydro data (geometry-only config) ──────
            %  Used when building config for the initial HAMS sweep.
            %  Downstream code must check config.hydro_ready before
            %  using wamit_* fields.
            fprintf('  Hydrodynamic data: SKIPPED (geometry-only config)\n');
            config.wamit_drafts = [];
            config.wamit_z_cg   = [];
            config.wamit_A      = [];
            config.wamit_B      = [];
            config.wamit_A_full = {};
            config.wamit_B_full = {};
            config.hydro_ready  = false;
            num_cases = 0;
        end
    end

    %% §6  PTO / MOORING  ─────────────────────────────────────────────

    config.enable_pto_effects = params.enable_pto_effects;
    config.K33_pto       = params.K33_pto;
    config.K55_pto       = params.K55_pto;
    config.K11_pto       = params.K11_pto;
    config.pto_angle_deg = params.pto_angle_deg;
    config.ht            = params.ht;
    config.bt            = params.bt;

    fprintf('  PTO mode: %s\n', ...
            ternary(config.enable_pto_effects, 'augmented', 'free-floating'));

    %% §7  OPTIMISATION TARGETS  ──────────────────────────────────────

    config.T_heave_goal  = params.T_heave_goal;
    config.T_pitch_goal  = params.T_pitch_goal;
    config.T_surge_goal  = params.T_surge_goal;
    config.T_heave_range = params.T_heave_range;
    config.T_pitch_range = params.T_pitch_range;
    config.T_surge_range = params.T_surge_range;

    config.gm_min    = params.gm_min;
    config.gm_range  = params.gm_range;
    config.gm_target = params.gm_target;

    fprintf('  Targets: T_h = %.1f s, T_p = %.1f s, GM = %.2f m\n', ...
            config.T_heave_goal, config.T_pitch_goal, config.gm_target);

    %% §8  OPTIMISATION BOUNDS  ────────────────────────────────────────
    %  Auto-compute from hull z-extents (v6.0 — full-body sweep).
    %
    %  INTENT
    %    The HAMS sweep must sample the waterline at evenly-spaced levels
    %    throughout the full hull height, from just below the top to just
    %    above the keel.  The design vector relationship is:
    %
    %      z_wl = -vs   (waterline in body frame)
    %
    %    so the bounds are:
    %
    %      vs_min = -hull_z_max + δ    (waterline δ below the top)
    %      vs_max = -hull_z_min - δ    (waterline δ above the keel)
    %
    %  WHY NOT the v5.0 centroid cap?
    %    v5.0 replaced vs_max = -hull_z_min - δ  with  vs_max = -centroid_z - δ
    %    to prevent the waterline from reaching the keel tip (V_sub → 0).
    %    That concern is real but is already handled at the physics level:
    %    screen_draft returns fval=Inf/feasible=false when V_sub < 1e-6, and
    %    calculate_3d_properties clamps the interpolation to the table range.
    %    There is no need to exclude those drafts from the bounds.
    %
    %    The centroid cap is NOT geometry-neutral.  For an inverted hull
    %    (C0_180), the volume centroid sits near the TOP of the body frame,
    %    so the cap clips vs_max to cover only ~23% of the hull height.
    %    Every HAMS draft then lands in the narrow apex-cone region, which
    %    produces degenerate waterplane meshes and NaN A(∞) for all runs.
    %    For C0 (upright) the centroid is nearer the keel so the error is
    %    smaller (67% coverage) and was not previously noticed.
    %
    %  NEAR-KEEL DRAFTS
    %    At large vs (waterline near keel), V_sub ≈ 0 and the HAMS mesh
    %    has only a few panels.  HAMS may fail or return A(∞)≈0.  Both
    %    outcomes are safe:
    %      • run_single_hams returns status='failed', A_inf=zeros — stored in cache
    %      • screen_draft flags those entries as infeasible (fval=Inf)
    %      • fmincon lb/ub still include the full range; SQP avoids the
    %        infeasible keel region through the mass-balance constraint
    %    The adaptive T_heave refinement skips intervals where T_h=Inf
    %    (isfinite(dT) guard), so no wasted HAMS runs are triggered there.

    if isfield(params, 'vertical_shift_bounds') && ~isempty(params.vertical_shift_bounds)
        config.vertical_shift_bounds = params.vertical_shift_bounds;
        fprintf('  Vertical shift bounds: [%.2f, %.2f] m (manual override)\n', ...
                config.vertical_shift_bounds);
    else
        delta_vs = 0.1;
        vs_min = -hull_z_max + delta_vs;   % waterline δ below hull top
        vs_max = -hull_z_min - delta_vs;   % waterline δ above hull keel

        config.vertical_shift_bounds = [vs_min, vs_max];
        fprintf('  Vertical shift bounds: [%.2f, %.2f] m (auto, wl ∈ [%.2f, %.2f] m)\n', ...
                vs_min, vs_max, -vs_max, -vs_min);
        fprintf('    Hull height: %.2f m  |  sweep coverage: %.0f%%\n', ...
                hull_z_max - hull_z_min, ...
                (vs_max - vs_min) / (hull_z_max - hull_z_min) * 100);
    end

    % CG target: mean of WAMIT z_cg values.
    % WHY mean?  Each WAMIT case was run at a different draft, so each
    % has a different CG assumption.  The mean gives a representative
    % reference that the optimiser can use for warm-starting or
    % diagnostics, without biasing toward any single draft case.
    if ~isempty(config.wamit_z_cg)
        config.z_cg_target = mean(config.wamit_z_cg);
    else
        config.z_cg_target = NaN;
    end

    %% §9  PID CORRECTION FACTORS (initial values)  ────────────────────
    %  Both start at unity = no correction.  The PID loop in
    %  WEC_Main_Optimizer drives them away from 1.0 during 'trained'
    %  mode.  In 'oneshot' and 'skip' modes they stay at 1.0.
    %
    %  mass_correction_factor = k_vol  (naming kept for backward compat)
    %  gm_correction_factor   = k_gm

    config.mass_correction_factor = params.k_vol_init;
    config.gm_correction_factor   = params.k_gm_init;

    % Runtime fields expected by calculate_2d_properties.
    % Stage 1 PID loop updates these; when stage1_mode = 'skip' they
    % stay at initial values.
    config.k_vol = params.k_vol_init;
    config.k_gm  = params.k_gm_init;

    %% §10  STAGE 1 PID TUNING  ───────────────────────────────────────
    %  Forward all PID gains, limits, damping, and convergence criteria
    %  so that WEC_Main_Optimizer reads them from config, not literals.

    config.stage1_mode             = params.stage1_mode;
    config.k_vol_init              = params.k_vol_init;
    config.k_gm_init               = params.k_gm_init;

    config.pid_mass_gains          = params.pid_mass_gains;
    config.pid_vol_gains           = params.pid_vol_gains;
    config.pid_gm_gains            = params.pid_gm_gains;
    config.pid_mass_limits         = params.pid_mass_limits;
    config.pid_vol_limits          = params.pid_vol_limits;
    config.pid_gm_limits           = params.pid_gm_limits;

    config.bounds_kvol             = params.bounds_kvol;
    config.bounds_kgm              = params.bounds_kgm;

    config.damping_vol_early       = params.damping_vol_early;
    config.damping_vol_late        = params.damping_vol_late;
    config.damping_gm_early        = params.damping_gm_early;
    config.damping_gm_late         = params.damping_gm_late;
    config.damping_transition_iter = params.damping_transition_iter;

    config.vol_conv_tol_pct        = params.vol_conv_tol_pct;
    config.gm_conv_tol             = params.gm_conv_tol;
    config.delta_kvol_stable       = params.delta_kvol_stable;
    config.delta_kgm_stable        = params.delta_kgm_stable;
    config.stable_count_needed     = params.stable_count_needed;
    config.mass_acceptable_pct     = params.mass_acceptable_pct;
    config.cg_guard_floor          = params.cg_guard_floor;
    config.pid_sat_proximity       = params.pid_sat_proximity;

    %% §11  OBJECTIVE-FUNCTION SHAPING  ────────────────────────────────

    config.period_delta  = params.period_delta;
    config.zone_k_amp    = params.zone_k_amp;
    config.penalty_guard = params.penalty_guard;

    %% §12  SOLVER SETTINGS  ──────────────────────────────────────────

    config.max_outer_iterations = params.max_outer_iterations;
    config.convergence_tol      = params.convergence_tol;
    config.n_sweep_refine       = params.n_sweep_refine;

    %% §13  AUTOCAD INERTIA VALIDATION  ────────────────────────────────

    config.autocad_Ixx = params.autocad_Ixx;
    config.autocad_Iyy = params.autocad_Iyy;
    config.autocad_Izz = params.autocad_Izz;
    config.autocad_discrepancy_pct = params.autocad_discrepancy_pct;

    %% §14  INITIAL CONDITIONS  ────────────────────────────────────────
    %  Smooth tanh density profile: heavy at bottom, light at top.
    %  Average density targets neutral buoyancy (rho_water).
    %
    %  CONSTRUCTION
    %    1. Normalise node z-positions to [0, 1] (z_norm).
    %    2. Apply tanh: heavy at z_norm=0 (keel), light at z_norm=1 (deck).
    %    3. Scale so that mean(density) = rho_core_target (see below).
    %    4. Check adjacent-ratio constraint; reduce steepness if violated.
    %    5. Clamp to physical bounds.
    %
    %  VARIABLE DICTIONARY
    %    z_norm          : [N×1]  normalised node heights, 0 = keel, 1 = deck
    %    density_profile : [N×1]  unitless shape (mean = 1 after normalisation)
    %    adj_ratios      : [N-1×1] rho(i) / rho(i+1) for constraint check
    %
    %  SHELL CORRECTION
    %    When the shell is enabled the design variables are CORE densities,
    %    not bulk densities.  Mass equilibrium requires:
    %
    %      rho_shell * V_shell + rho_core_avg * V_core = rho_water * V_hull
    %
    %    Targeting rho_water for the core ignores the shell mass, which is
    %    already 2400 * V_shell kg.  The initial guess is then far too heavy,
    %    and fmincon drives draft negative (more submersion) to recover
    %    buoyancy -- pushing the waterplane onto the narrow column and sending
    %    T_heave to ~19 s.
    %
    %    Correct target:
    %      rho_core_target = (rho_water * V_hull - rho_shell * V_shell) / V_core

    % Initial vertical shift: put the waterline at the hull volume centroid.
    %
    %   vs = -centroid_z  →  z_wl = centroid_z  (waterline at centroid)
    %
    %  This guarantees ≈ 50% of hull volume is submerged for any hull
    %  orientation — a physically realistic warm-start for both upright
    %  (C0) and inverted (C0_180) geometries.
    %
    %  Using mean(vs_bounds) was broken by the v6.0 full-body sweep fix:
    %  the midpoint of the new wider bounds sits near the keel for
    %  inverted hulls (V_sub ≈ 0 there), producing an infeasible x0.
    %  The centroid criterion is geometry-neutral.
    %
    %  Clamp to vs_bounds so the initial point is always inside the
    %  optimisation domain.
    vs_at_centroid = -config.hull_centroid(3);
    config.initial_vertical_shift = max(config.vertical_shift_bounds(1), ...
                                    min(config.vertical_shift_bounds(2), vs_at_centroid));

    rho_avg_target = config.RHO_WATER;
    z_norm = (config.density_nodes_z - config.density_nodes_z(1)) / ...
             (config.density_nodes_z(end) - config.density_nodes_z(1) + eps);

    shape = params.tanh_shape_param;
    density_profile = 0.5 + 0.5 * tanh(shape * (0.5 - z_norm));
    density_profile = density_profile / mean(density_profile);
    config.initial_densities = rho_avg_target * density_profile;

    % Check adjacent-ratio constraint; reduce steepness if violated
    adj_ratios = config.initial_densities(1:end-1) ./ ...
                (config.initial_densities(2:end) + eps);
    max_adj = max(adj_ratios);

    if max_adj > config.max_density_ratio
        fprintf('  Initial profile ratio %.2f > limit %.1f — reducing steepness\n', ...
                max_adj, config.max_density_ratio);
        shape = params.tanh_shape_fallback;
        density_profile = 0.5 + 0.5 * tanh(shape * (0.5 - z_norm));
        density_profile = density_profile / mean(density_profile);
        config.initial_densities = rho_avg_target * density_profile;
        adj_ratios = config.initial_densities(1:end-1) ./ ...
                    (config.initial_densities(2:end) + eps);
        max_adj = max(adj_ratios);
    end

    % Clamp to physical bounds and force row vector
    config.initial_densities = max(config.ballast_density_bounds(1), ...
        min(config.ballast_density_bounds(2), config.initial_densities));
    config.initial_densities = config.initial_densities(:)';

    % Clamp to per-strip density lower bounds (constructability mode)
    if config.enable_constructability && ~isempty(config.per_strip_density_lb)
        for i = 1:length(config.initial_densities)
            config.initial_densities(i) = max(config.initial_densities(i), ...
                                              config.per_strip_density_lb(i));
        end
    end

    % Pin wall strip to rho_UHPC (constructability mode)
    if config.enable_constructability && ~isempty(config.wall_strip_index)
        w_idx = config.wall_strip_index;
        config.initial_densities(w_idx) = params.constructability_rho_hull;
        fprintf('  Wall strip %d pinned to %.0f kg/m^3 (UHPC)\n', ...
                w_idx, params.constructability_rho_hull);
    end

    fprintf('  Initial vertical shift: %.2f m\n', config.initial_vertical_shift);
    fprintf('  Initial density: [%.0f ... %.0f] kg/m^3, max adj ratio: %.2f\n', ...
            config.initial_densities(1), config.initial_densities(end), max_adj);

    %% §15  VALIDATION  ────────────────────────────────────────────────

    % Hydro data validation (skip when geometry-only config)
    if config.hydro_ready
        assert(length(config.wamit_drafts) == size(config.wamit_A, 1), ...
               'WAMIT draft count does not match added-mass rows');

        for i = 1:num_cases
            assert(isequal(size(config.wamit_A_full{i}), [3, 3]), ...
                   'wamit_A_full{%d} is not 3x3', i);
            assert(isequal(size(config.wamit_B_full{i}), [3, 3]), ...
                   'wamit_B_full{%d} is not 3x3', i);
        end
    end

    assert(config.num_ballast_sections >= 1, ...
           'Need at least 1 ballast section');
    assert(config.max_density_ratio >= 1, ...
           'max_density_ratio must be >= 1');
    assert(config.gm_range(1) >= 0, ...
           'Minimum GM must be non-negative');
    assert(config.gm_range(2) > config.gm_range(1), ...
           'Max GM must exceed min GM');

    if config.gm_target < config.gm_range(1) || ...
       config.gm_target > config.gm_range(2)
        warning('WEC:GMOutOfRange', ...
                'gm_target (%.2f) is outside gm_range [%.2f, %.2f]', ...
                config.gm_target, config.gm_range(1), config.gm_range(2));
    end

    %  Mutual exclusion: steel-fill and constructability are both
    %  post-optimisation realisers — running both would overwrite each
    %  other's final_props with no well-defined ordering.
    assert(~(config.enable_steel_solve && config.enable_constructability), ...
           'Mutually exclusive — set exactly one to true in WEC_Driver.m');

    if config.enable_constructability
        assert(config.num_ballast_sections >= 3, ...
               'WEC:TooFewStrips', ...
               'Constructability mode requires num_ballast_sections >= 3 (1 wall + 2 platform minimum).');
        assert(~isempty(config.wall_strip_index), ...
               'WEC:NoWallIndex', ...
               'Constructability enabled but wall_strip_index was not set.');
    end

    fprintf('  Validation passed\n');

    %% §16  HAMS RUNTIME PATHS  ───────────────────────────────────────
    %  Forward HAMS directory, executable, and cache file from params.
    %  The outer loop in WEC_Main_Optimizer uses these to run HAMS at
    %  converged vertical_shifts and enrich the hydro cache on-the-fly.

    if isfield(params, 'hams_dir') && ~isempty(params.hams_dir)
        config.hams_dir         = params.hams_dir;
        config.hams_exe         = params.hams_exe;
        config.hydro_cache_file = params.hydro_cache_file;

        % Load hydro cache into config for outer-loop access
        if isfield(params, 'hydro_cache_file') && ...
                exist(params.hydro_cache_file, 'file')
            loaded = load(params.hydro_cache_file, 'hydro_table');
            config.hydro_cache = loaded.hydro_table;
        else
            config.hydro_cache = HAMS_Pipeline.empty_hydro_cache();
        end
    else
        config.hams_dir         = '';
        config.hams_exe         = '';
        config.hydro_cache_file = '';
        config.hydro_cache      = HAMS_Pipeline.empty_hydro_cache();
    end

    fprintf('=================================================\n\n');

catch ME
    fprintf('\n  CONFIGURATION FAILED: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('    in %s, line %d\n', ME.stack(1).name, ME.stack(1).line);
    end
    rethrow(ME);
end

end

%% ─────────────────────────────────────────────────────────────────────
function r = ternary(cond, t, f)
    if cond, r = t; else, r = f; end
end