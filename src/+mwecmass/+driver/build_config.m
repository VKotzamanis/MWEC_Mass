function config = build_config(in, hydro_table, mesh_sizing)
%BUILD_CONFIG Build the unified hydrostatic/geometry config struct consumed by the whole pipeline.
% Assemble geometry, ballast, hydrostatic, hydrodynamic, optimisation, and
% realisation settings into one config struct.  All author-set values come
% from in; hydro_table is an optional loaded BEM table and mesh_sizing carries
% optional panel counts.  Empty hydro_table produces a geometry-only config.
% The model is free-floating and has no PTO or mooring fields.
%

try
    fprintf('\n========== WEC CONFIGURATION BUILDER ==========\n');

    %% Paths and mode
    %  The only place an author-set value is reshaped rather than copied. Everything else in
    %  this function is a straight copy out of `in`, so a value has exactly one source.
    %
    %  Four fileparts() calls: this file sits three directories below the repository root
    %  (src/+mwecmass/+driver/build_config.m), so the chain walks file -> +driver -> +mwecmass
    %  -> src -> root. The file names in `in` are bare names resolved under Input/ (the deck and
    %  the BEM cache) or repository-relative paths (the HAMS-MREL workspace and binary).
    repo_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));

    ms2_file        = fullfile(repo_root, 'Input', in.files.ms2_file);      % hull geometry deck
    bem_cache_file  = fullfile(repo_root, 'Input', in.bem.bem_cache_file);  % BEM cache, read only
    hams_dir        = fullfile(repo_root, in.bem.hams_dir);                 % HAMS-MREL workspace
    hams_exe        = fullfile(repo_root, in.bem.hams_exe);                 % HAMS-MREL binary
    % Write target of the HAMS-MREL route's own BEM cache, under Output/hams_mrel/ and never
    % under Input/: the draft-node sweep saves into it and the optimiser rewrites it when a run
    % has added draft nodes. Path value only; nothing reaches it while in.bem.run_HAMS_MREL is
    % false.
    hams_cache_file = fullfile(repo_root, in.bem.hams_dir, in.bem.hams_cache_file);

    % Geometry cache folder: repository-relative unless absolute. Empty or absent switches the
    % cache off.
    if isfield(in.files, 'geometry_cache_dir') && ~isempty(in.files.geometry_cache_dir)
        geometry_cache_dir = in.files.geometry_cache_dir;
        if isempty(regexp(geometry_cache_dir, '^([A-Za-z]:|[\\/])', 'once'))
            geometry_cache_dir = fullfile(repo_root, geometry_cache_dir);
        end
    else
        geometry_cache_dir = '';
    end

    % The realisation type selects the post-optimisation model; the stages read it as the two
    % mutually-exclusive booleans below ("mode" is reserved for the Stage-1 draft search).
    switch in.materials.realisation_type
        case 'thin_shell'
            enable_steel_solve      = true;
            enable_constructability = false;
        case 'modular_precast'
            enable_steel_solve      = false;
            enable_constructability = true;
        case 'preliminary'
            enable_steel_solve      = false;
            enable_constructability = false;
    end

    % Wall position: read off the deck's file name, not an author-set field.
    if contains(ms2_file, '_180')
        wall_position = 'bottom';
    else
        wall_position = 'top';
    end

    if nargin < 3 || isempty(mesh_sizing), mesh_sizing = struct(); end

    %% Physical constants
    %  Copied directly from the input struct.  Units: SI throughout.

    config.RHO_WATER    = in.constants.rho_water;  % [kg/m^3]
    config.G            = in.constants.g;          % [m/s^2]
    % config.seabed_depth is not set: in.context.seabed_depth no longer exists and the suite
    % carries one site water depth only, in.bem.water_depth.
    config.wamit_L      = in.constants.wamit_L;    % [m]
    config.ms2_file     = ms2_file;                % for diagnostic reporting

    fprintf('  Constants: rho = %.0f kg/m^3, g = %.2f m/s^2\n', ...
            config.RHO_WATER, config.G);

    %% Geometry products
    %  The expensive geometry steps run in compute_geometry_products at the end of this file: the
    %  deck's surface extents and volume integrals, the hydrostatic tables with their boundary
    %  cache, the strip layout and per-strip density bounds, the strip geometry and the Y-span
    %  table. That function receives only the values geometry_inputs collects, so the deck and
    %  those values are its whole input. Its result holds plain arrays and structs; it is saved
    %  under in.files.geometry_cache_dir and reloaded when the key of
    %  mwecmass.internal.geometry_cache (deck bytes, those values, runtime version and the source
    %  of every file the steps can call) is unchanged.
    %
    %  VARIABLE DICTIONARY (outputs stored in config)
    %    ms2_model           : mwecmass.geometry.MS2Parser object (evaluatable geometry)
    %    profile             : [Mx2]  midplane polygon (x, z) for 2D model
    %    total_wec_volume    : [m^3]  displaced volume of record

    % Store panelizer grid density (used only for .pnl export, not in optimiser). Guarded:
    % the counts come from mwecmass.mesh.panel_grid_counts on the HAMS-MREL path and from the
    % loaded cache's own metadata on the WAMIT path; when the caller supplies neither, [] here is
    % a documented default, not a fabricated value.
    if isfield(mesh_sizing, 'mesh_Nu'), config.mesh_Nu = mesh_sizing.mesh_Nu; else, config.mesh_Nu = []; end
    if isfield(mesh_sizing, 'mesh_Nv'), config.mesh_Nv = mesh_sizing.mesh_Nv; else, config.mesh_Nv = []; end

    g = geometry_inputs(in, enable_constructability, wall_position);
    key = '';
    hit = false;
    if ~isempty(geometry_cache_dir) && exist(ms2_file, 'file')   % a missing deck errors in the parse
        [~, deck_stem] = fileparts(ms2_file);
        cache_file = fullfile(geometry_cache_dir, ...
                              [deck_stem '_' in.materials.realisation_type '_geometry.mat']);
        key = mwecmass.internal.geometry_cache('key', ms2_file, g, fullfile(repo_root, 'src'));
        [hit, products, reason] = mwecmass.internal.geometry_cache('load', cache_file, key);
    end
    if hit
        fprintf('  Geometry products: reloaded from %s\n', cache_file);
        ms2_model = mwecmass.geometry.MS2Parser.parse(ms2_file);
    else
        if isempty(key)
            fprintf('  Geometry products: cache off, computing\n');
        else
            fprintf('  Geometry products: %s, computing\n', reason);
        end
        [products, ms2_model] = compute_geometry_products(ms2_file, g);
        if ~isempty(key) && mwecmass.internal.geometry_cache('save', cache_file, key, products)
            fprintf('  Geometry products: saved to %s\n', cache_file);
        end
    end

    config.ms2_model = ms2_model;
    product_names = fieldnames(products.config);
    for k = 1:numel(product_names)
        config.(product_names{k}) = products.config.(product_names{k});
    end
    for k = 1:size(products.warnings, 1)
        warning(products.warnings{k, 1}, '%s', products.warnings{k, 2});
    end
    hull_z_min = config.hull_z_min;
    hull_z_max = config.hull_z_max;

    % Topology struct — repurposed for strip geometry data
    config.topology = struct();

    %% Ballast configuration
    %  The strip layout (density node heights, strip edges, wall strip index) comes with the
    %  geometry products; the three values below are plain copies.

    config.num_ballast_sections   = in.geometry.num_ballast_sections;
    config.ballast_density_bounds = mode_density_bounds(in);
    config.max_density_ratio      = in.bounds.max_density_ratio;

    %% Thin-shell solver parameters
    %  Forward parameters consumed by mwecmass.realise.thin_shell.solve, which runs
    %  The inverse solve runs once after Stage-2 optimisation converges.
    %  Nothing is computed here — the inverse solve needs final_props.
    %
    %  config.shell remains empty because shell and ballast are applied after optimisation.

    config.enable_steel_solve      = enable_steel_solve;
    config.rho_steel               = in.materials.thin_shell.rho_shell;
    config.rho_air                 = in.materials.thin_shell.rho_air;
    % config.rho_shell is the shell-region density of the two-density thin-shell split; the
    % input file carries one shell density, so both fields read in.materials.thin_shell.rho_shell.
    config.rho_shell                = in.materials.thin_shell.rho_shell;
    config.rho_ballast             = thin_shell_rho_ballast(in);
    config.steel_t_init            = in.materials.thin_shell.t_init;
    config.steel_t_min             = in.materials.thin_shell.t_min;
    config.steel_max_slope_factor  = in.materials.thin_shell.max_slope_factor;
    config.steel_n_z_grid          = in.materials.thin_shell.n_z_grid;
    config.shell                   = [];   % shell and ballast are applied after optimisation

    % HAMS period-grid forwards (consumed by mwecmass.bem.hams_mrel.default_hams_params
    % via the config-override path).  Optional — defaults apply if absent.
    if isfield(in.bem, 'T_min')
        config.period_min  = in.bem.T_min;
    end
    if isfield(in.bem, 'T_max')
        config.period_max  = in.bem.T_max;
    end
    if isfield(in.bem, 'T_step')
        config.period_step = in.bem.T_step;
    end
    % Site water depth: guarded the same way as the T_min/T_max/T_step forwards above.
    % mwecmass.bem.hams_mrel.default_hams_params(config) reads config.water_depth; when it is
    % absent (e.g. an empty_hydro_cache() call, which passes no config) that
    % function's own nargin<1 branch supplies config = struct() and has no fallback of its own --
    % The normal driver supplies water depth; an empty cache remains usable without it.
    if isfield(in.bem, 'water_depth')
        config.water_depth = in.bem.water_depth;
    end

    if config.enable_steel_solve
        fprintf('  Steel-fill solver: ENABLED (post-optimisation, see Stage 3: material realisation)\n');
        % Print the three densities consumed by the two-density thin-shell solve.
        fprintf('    rho_shell = %.0f kg/m^3   rho_ballast = %.0f kg/m^3   rho_air = %.0f kg/m^3\n', ...
                config.rho_shell, config.rho_ballast, config.rho_air);
        fprintf('    t_steel initial guess: %.4f m (%.2f in)\n', ...
                config.steel_t_init, config.steel_t_init / 0.0254);
        fprintf('    t_steel min (fabrication floor): %.5f m (%.2f in)\n', ...
                config.steel_t_min, config.steel_t_min / 0.0254);
    else
        fprintf('  Steel-fill solver: DISABLED\n');
    end

    %% Constructability parameters
    %  Forward all constructability parameters to config.
    %  When disabled, all fields exist but are ignored downstream.

    config.enable_constructability       = enable_constructability;
    config.constructability_rho_hull     = in.materials.modular_precast.rho_hull;
    config.constructability_rho_air      = in.materials.modular_precast.rho_air;
    config.constructability_t_min        = in.materials.modular_precast.t_min;
    config.constructability_wall_height  = in.materials.modular_precast.wall_height;
    config.constructability_n_sub        = in.materials.modular_precast.n_sub;

    % UHPC global-solve tuning (forwarded only when the input struct carries the field)
    if isfield(in.materials.modular_precast, 't_init')
        config.uhpc_t_init = in.materials.modular_precast.t_init;
    end
    if isfield(in.materials.modular_precast, 'max_slope_factor')
        config.uhpc_max_slope_factor = in.materials.modular_precast.max_slope_factor;
    end
    if isfield(in.materials.modular_precast, 'n_z_grid')
        config.uhpc_n_z_grid = in.materials.modular_precast.n_z_grid;
    end

    if config.enable_constructability
        fprintf('  Constructability post-processing: ENABLED\n');
        fprintf('    Hull material:  %.0f kg/m³ (UHPC)\n', config.constructability_rho_hull);
        fprintf('    Void air:       %.1f kg/m³\n', config.constructability_rho_air);
        fprintf('    Min thickness:  %.4f m (%.1f in)\n', ...
                config.constructability_t_min, config.constructability_t_min / 0.0254);
        fprintf('    Wall height:    %.2f m\n', config.constructability_wall_height);
    else
        fprintf('  Constructability post-processing: DISABLED\n');
    end

    %% Minimum constructable mass
    %
    %  The strip volumes times their lower density bounds were summed with the geometry products.
    %  Compare this minimum mass with full-submergence buoyancy before solving.

    if config.enable_constructability && ~isempty(config.strip_edges) ...
            && ~isempty(config.per_strip_density_lb)
        m_min = products.report.m_min_constructability;
        m_max = products.report.m_max_constructability;
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
    end

    config.n_density_strips = in.geometry.num_ballast_sections;   % alias used by the 2-D model
    config.eff_w_floor = in.geometry.eff_w_floor;

    fprintf('  Y-span table: %d levels, corrected width range [%.4f, %.4f] m\n', ...
            in.geometry.n_z_levels, min(config.y_span_table), max(config.y_span_table));
    fprintf('  Strip count: %d\n', config.n_density_strips);
    %% Hydrodynamic data
    %  Two paths:
    %    (A) HAMS hydro_table — precomputed by HAMS_DraftSweep.  Contains
    %        A(∞) and B_avg at N drafts, 6×6 at origin.  Extracted to 3×3
    %        [surge, heave, pitch] and transformed from origin to CG here.
    %    (B) WAMIT input files, when supplied by the loaded table.
    %
    %  Both paths populate the SAME config fields:
    %    config.hydro_drafts, config.hydro_z_cg, config.added_mass_diagonal,
    %    config.added_mass_full, config.radiation_damping_diagonal, config.radiation_damping_full
    %  so that ALL downstream code (mwecmass.bem.interpolate_at_draft,
    %  properties_3d, coupled eigenvalue analysis) is UNCHANGED.

    %% Attach the supplied hydrodynamic table
    %  hydro_table is supplied by the caller (mwecmass.bem.load_hydro_cache or
    %  mwecmass.bem.hams_mrel.run, selected by mwecmass.driver.run's own Stage-1 branch), so this
    %  function never reads a WAMIT case list itself. The second argument is optional: when
    %  omitted or empty, the geometry-only branch below runs, for callers that build a config
    %  before any hydro_table exists (the HAMS draft-node sweep and the WAMIT input writer).
    if nargin < 2 || isempty(hydro_table)
        % Geometry-only branch: no hydrodynamic data is attached.
        fprintf('  Hydrodynamic data: SKIPPED (geometry-only config)\n');
        config.hydro_drafts               = [];
        config.hydro_z_cg                 = [];
        config.added_mass_diagonal        = [];
        config.radiation_damping_diagonal = [];
        config.added_mass_full            = {};
        config.radiation_damping_full     = {};
        config.hydro_ready  = false;
        num_cases = 0;
    else
    config.hydro_table = hydro_table;

        N_ht = length(config.hydro_table.drafts);
        config.hydro_drafts               = config.hydro_table.drafts(:);   % [m] draft (vertical-shift) nodes, signed z up
        config.hydro_z_cg                 = config.hydro_table.z_cg(:);     % [m] vertical CG each node's matrices are referred to
        config.added_mass_diagonal        = zeros(N_ht, 3);                 % [A11 kg, A33 kg, A55 kg m^2] per draft node
        config.radiation_damping_diagonal = zeros(N_ht, 3);                 % [B11 N s/m, B33 N s/m, B55 N m s/rad], band-averaged
        config.added_mass_full            = cell(N_ht, 1);                  % 3x3 [surge, heave, pitch] added mass per node; surge-pitch terms kg m
        config.radiation_damping_full     = cell(N_ht, 1);                  % 3x3 band-averaged radiation damping per node

        idx_3dof = [1, 3, 5];  % surge, heave, pitch from 6-DOF

        for i = 1:N_ht
            A_inf_6x6 = config.hydro_table.added_mass_inf{i};   % 6×6 at origin
            B_avg_6x6 = config.hydro_table.radiation_damping_band_avg{i};   % 6×6 at origin
            z_cg_i    = config.hydro_table.z_cg(i);     % CG_z in global frame

            % Extract 3×3 [surge, heave, pitch] submatrix
            A_3x3_origin = A_inf_6x6(idx_3dof, idx_3dof);
            B_3x3_origin = B_avg_6x6(idx_3dof, idx_3dof);

            % Transform from origin to CG.
            [A_cg, B_cg] = mwecmass.bem.transform_to_cg( ...
                A_3x3_origin, B_3x3_origin, z_cg_i);

            config.added_mass_diagonal(i, :)        = [A_cg(1,1), A_cg(2,2), A_cg(3,3)];
            config.radiation_damping_diagonal(i, :) = [B_cg(1,1), B_cg(2,2), B_cg(3,3)];
            config.added_mass_full{i}               = A_cg;
            config.radiation_damping_full{i}        = B_cg;
        end

        % Sort by draft (ascending)
        [config.hydro_drafts, si]         = sort(config.hydro_drafts);
        config.hydro_z_cg                 = config.hydro_z_cg(si);
        config.added_mass_diagonal        = config.added_mass_diagonal(si, :);
        config.radiation_damping_diagonal = config.radiation_damping_diagonal(si, :);
        config.added_mass_full            = config.added_mass_full(si);
        config.radiation_damping_full     = config.radiation_damping_full(si);

        num_cases = N_ht;
        config.hydro_ready = true;
        fprintf('  HAMS hydro_table: %d drafts, A(inf) + B_avg at CG\n', N_ht);
        for i = 1:N_ht
            fprintf('    draft=%+.2f m, A_diag=[%.1f, %.1f, %.1f]\n', ...
                    config.hydro_drafts(i), config.added_mass_diagonal(i,:));
        end
    end

    %% Free-floating model
    %  The suite models a free-floating body: no PTO or mooring fields exist on config.


    %% Optimisation targets

    config.T_heave_goal  = in.targets.T_heave_goal;
    config.T_pitch_goal  = in.targets.T_pitch_goal;
    config.T_heave_range = in.targets.T_heave_range;
    config.T_pitch_range = in.targets.T_pitch_range;

    config.gm_min    = in.targets.gm_min;
    config.gm_range  = in.targets.gm_range;
    config.gm_target = in.targets.gm_target;

    fprintf('  Targets: T_h = %.1f s, T_p = %.1f s, GM = %.2f m\n', ...
            config.T_heave_goal, config.T_pitch_goal, config.gm_target);

    %% Optimisation bounds
    %  Auto-compute from the full hull z-extents.
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
    %  Use the full-body range rather than a centroid cap.
    %    A centroid cap can exclude valid waterline positions.
    %    to prevent the waterline from reaching the keel tip (V_sub → 0).
    %    That concern is real but is already handled at the physics level:
    %    stage1_screen_draft returns fval=Inf/feasible=false when V_sub < 1e-6, and
    %    properties_3d clamps the interpolation to the table range.
    %    There is no need to exclude those drafts from the bounds.
    %
    %    The centroid cap is NOT geometry-neutral.  For an inverted hull
    %    (C0_180), the volume centroid sits near the TOP of the body frame,
    %    so the cap clips vs_max to cover only ~23% of the hull height.
    %    Every HAMS draft then lands in the narrow apex-cone region, which
    %    produces degenerate waterplane meshes and NaN A(∞) for all runs.
    %    For C0 (upright) the centroid is nearer the keel so the error is
    %    Upright and inverted hulls then receive different sweep coverage.
    %
    %  NEAR-KEEL DRAFTS
    %    At large vs (waterline near keel), V_sub ≈ 0 and the HAMS mesh
    %    has only a few panels.  HAMS may fail or return A(∞)≈0.  Both
    %    outcomes are safe:
    %      • run_at_draft returns status='failed', A_inf=zeros — stored in cache
    %      • stage1_screen_draft flags those entries as infeasible (fval=Inf)
    %      • fmincon lb/ub still include the full range; SQP avoids the
    %        infeasible keel region through the mass-balance constraint
    %    The adaptive T_heave refinement skips intervals where T_h=Inf
    %    (isfinite(dT) guard), so no wasted HAMS runs are triggered there.

    if isfield(in.bounds, 'vertical_shift_bounds') && ~isempty(in.bounds.vertical_shift_bounds)
        config.vertical_shift_bounds = in.bounds.vertical_shift_bounds;
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
    % Each hydrodynamic case can use a different draft CG, so the mean
    % has a different CG assumption.  The mean gives a representative
    % reference that the optimiser can use for warm-starting or
    % diagnostics, without biasing toward any single draft case.
    if ~isempty(config.hydro_z_cg)
        config.z_cg_target = mean(config.hydro_z_cg);
    else
        config.z_cg_target = NaN;
    end

    %% PID correction factors
    %  Both start at unity = no correction.  The PID loop in
    %  mwecmass.optim.run drives them away from 1.0 during 'trained'
    %  mode.  In 'oneshot' and 'skip' modes they stay at 1.0.
    %
    %  mass_correction_factor = k_vol  (naming kept for backward compat)
    %  gm_correction_factor   = k_gm

    config.mass_correction_factor = in.pid.k_vol_init;
    config.gm_correction_factor   = in.pid.k_gm_init;

    % Runtime fields expected by properties_2d.
    % Stage 1 PID loop updates these; when stage1_mode = 'skip' they
    % stay at initial values.
    config.k_vol = in.pid.k_vol_init;
    config.k_gm  = in.pid.k_gm_init;

    %% Stage-1 PID tuning
    %  Forward all PID gains, limits, damping, and convergence criteria
    %  so that mwecmass.optim.run reads them from config, not literals.

    config.stage1_mode             = in.pid.stage1_mode;
    config.k_vol_init              = in.pid.k_vol_init;
    config.k_gm_init               = in.pid.k_gm_init;

    config.pid_mass_gains          = in.pid.mass_gains;
    config.pid_vol_gains           = in.pid.vol_gains;
    config.pid_gm_gains            = in.pid.gm_gains;
    config.pid_mass_limits         = in.pid.mass_limits;
    config.pid_vol_limits          = in.pid.vol_limits;
    config.pid_gm_limits           = in.pid.gm_limits;

    config.bounds_kvol             = in.pid.bounds_kvol;
    config.bounds_kgm              = in.pid.bounds_kgm;

    config.damping_vol_early       = in.pid.damping_vol_early;
    config.damping_vol_late        = in.pid.damping_vol_late;
    config.damping_gm_early        = in.pid.damping_gm_early;
    config.damping_gm_late         = in.pid.damping_gm_late;
    config.damping_transition_iter = in.pid.damping_transition_iter;

    config.vol_conv_tol_pct        = in.pid.vol_conv_tol_pct;
    config.gm_conv_tol             = in.pid.gm_conv_tol;
    config.delta_kvol_stable       = in.pid.delta_kvol_stable;
    config.delta_kgm_stable        = in.pid.delta_kgm_stable;
    config.stable_count_needed     = in.pid.stable_count_needed;
    config.mass_acceptable_pct     = in.pid.mass_acceptable_pct;
    config.cg_guard_floor          = in.pid.cg_guard_floor;
    config.pid_sat_proximity       = in.pid.sat_proximity;

    %% Objective shaping

    config.zone_k_amp    = in.objective.zone_k_amp;
    config.penalty_guard = in.objective.penalty_guard;

    %% Solver settings
    %  Stage-2 fmincon algorithm (default 'sqp' — production setting).
    if isfield(in.solver, 'stage2_algorithm') && ~isempty(in.solver.stage2_algorithm)
        config.stage2_algorithm = in.solver.stage2_algorithm;
    else
        config.stage2_algorithm = 'sqp';
    end

    config.max_outer_iterations = in.solver.max_outer_iterations;
    config.n_sweep_refine       = in.solver.n_sweep_refine;

    %% CAD inertia validation

    config.autocad_Ixx = in.validation.autocad_Ixx;
    config.autocad_Iyy = in.validation.autocad_Iyy;
    config.autocad_Izz = in.validation.autocad_Izz;
    config.autocad_discrepancy_pct = in.validation.autocad_discrepancy_pct;

    %% Initial conditions
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
    %  The bounds midpoint can place an inverted hull near the keel:
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

    shape = in.pid.tanh_shape_param;
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
        shape = in.pid.tanh_shape_fallback;
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

    % Clamp to per-strip density lower bounds (constructability realisation type)
    if config.enable_constructability && ~isempty(config.per_strip_density_lb)
        for i = 1:length(config.initial_densities)
            config.initial_densities(i) = max(config.initial_densities(i), ...
                                              config.per_strip_density_lb(i));
        end
    end

    % Pin wall strip to rho_UHPC (constructability realisation type)
    if config.enable_constructability && ~isempty(config.wall_strip_index)
        w_idx = config.wall_strip_index;
        config.initial_densities(w_idx) = in.materials.modular_precast.rho_hull;
        fprintf('  Wall strip %d pinned to %.0f kg/m^3 (UHPC)\n', ...
                w_idx, in.materials.modular_precast.rho_hull);
    end

    fprintf('  Initial vertical shift: %.2f m\n', config.initial_vertical_shift);
    fprintf('  Initial density: [%.0f ... %.0f] kg/m^3, max adj ratio: %.2f\n', ...
            config.initial_densities(1), config.initial_densities(end), max_adj);

    %% Validation

    % Hydro data validation (skip when geometry-only config)
    if config.hydro_ready
        assert(length(config.hydro_drafts) == size(config.added_mass_diagonal, 1), ...
               'Cached draft-node count does not match added-mass rows');

        for i = 1:num_cases
            assert(isequal(size(config.added_mass_full{i}), [3, 3]), ...
                   'added_mass_full{%d} is not 3x3', i);
            assert(isequal(size(config.radiation_damping_full{i}), [3, 3]), ...
                   'radiation_damping_full{%d} is not 3x3', i);
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
           'Mutually exclusive — in.materials.realisation_type selects exactly one (WEC_User_Input.m)');

    if config.enable_constructability
        assert(config.num_ballast_sections >= 3, ...
               'WEC:TooFewStrips', ...
               'Constructability realisation type requires num_ballast_sections >= 3 (1 wall + 2 platform minimum).');
        assert(~isempty(config.wall_strip_index), ...
               'WEC:NoWallIndex', ...
               'Constructability enabled but wall_strip_index was not set.');
    end

    fprintf('  Validation passed\n');

    %% HAMS runtime paths
    %  Forward the HAMS-MREL directory, executable and cache file resolved above.
    %  The outer loop in mwecmass.optim.run uses these to run HAMS at
    %  converged vertical_shifts and enrich the hydro cache on-the-fly.
    %  config.hydro_cache_file is that enrichment's WRITE target (under
    %  Output/hams_mrel/); config.hydro_cache below is READ from the WAMIT
    %  cache under Input/ (in.bem.bem_cache_file), which this suite never
    %  writes.

    % config.run_HAMS_MREL: the flag mwecmass.optim.hams_enrichment_action reads to decide
    % whether mwecmass.optim.stage1_trained runs its per-iteration HAMS-enrichment call or
    % skips it in favour of cache interpolation. Optional -- defaults to false (WAMIT-cache path)
    % for an input struct that does not carry the switch; same isfield-guarded-optional pattern as
    % period_min/period_max/period_step above.
    if isfield(in.bem, 'run_HAMS_MREL')
        config.run_HAMS_MREL = in.bem.run_HAMS_MREL;
    else
        config.run_HAMS_MREL = false;
    end

    if isfield(in.bem, 'hams_dir') && ~isempty(in.bem.hams_dir)
        config.hams_dir         = hams_dir;
        config.hams_exe         = hams_exe;
        config.hydro_cache_file = hams_cache_file;

        % Load the BEM cache into config for outer-loop access. The source is bem_cache_file, the
        % cache this suite reads; config.hydro_cache_file above is the file the optimiser WRITES
        % when a run has enriched the cache, and with the HAMS-MREL route supplying its own cache
        % under Output/hams_mrel/ the two are no longer the same file. Reading the write target
        % here would leave config.hydro_cache empty on the WAMIT-cache path, so the read keeps its
        % own field.
        if exist(bem_cache_file, 'file')
            loaded = load(bem_cache_file, 'hydro_table');
            config.hydro_cache = loaded.hydro_table;
            if config.hydro_cache.water_depth ~= config.water_depth, error('mwecmass:driver:WaterDepthMismatch', 'Cache water_depth (%.15g m) differs from config.water_depth (%.15g m).', config.hydro_cache.water_depth, config.water_depth); end
        else
            config.hydro_cache = mwecmass.bem.empty_hydro_cache();
        end
    else
        config.hams_dir         = '';
        config.hams_exe         = '';
        config.hydro_cache_file = '';
        config.hydro_cache      = mwecmass.bem.empty_hydro_cache();
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

function g = geometry_inputs(in, enable_constructability, wall_position)
%GEOMETRY_INPUTS The author inputs compute_geometry_products reads, and nothing else.
% The wall inputs exist only for the modular-precast realisation; the wall position (read off the
% deck's file name) only matters with them.
    g.num_ballast_sections    = in.geometry.num_ballast_sections;
    g.n_z_levels              = in.geometry.n_z_levels;
    if isfield(in.geometry, 'aw_table_dz')
        g.aw_table_dz         = in.geometry.aw_table_dz;   % [m] target z-spacing of the coarse grid
    else
        g.aw_table_dz         = [];                        % absent -> the builder's own fallback grid
    end
    g.ballast_density_bounds  = mode_density_bounds(in);
    g.enable_constructability = enable_constructability;
    if enable_constructability
        g.wall_position = wall_position;
        g.wall_height   = in.materials.modular_precast.wall_height;
        g.rho_hull      = in.materials.modular_precast.rho_hull;
        g.rho_air       = in.materials.modular_precast.rho_air;
        g.t_min         = in.materials.modular_precast.t_min;
        g.n_sub         = in.materials.modular_precast.n_sub;
    else
        g.wall_position = '';
        g.wall_height   = [];
        g.rho_hull      = [];
        g.rho_air       = [];
        g.t_min         = [];
        g.n_sub         = [];
    end
end

function bounds = mode_density_bounds(in)
%MODE_DENSITY_BOUNDS Per-strip density bounds [lo, hi] of the realisation type.
% The lower bound is the author input; the upper bound is the solid density of the type's
% material (steel ballast for thin shell, UHPC for modular precast). Only 'preliminary', which
% has no material, takes the upper bound from the input.
    bounds = in.bounds.ballast_density_bounds;
    switch in.materials.realisation_type
        case 'thin_shell'
            bounds(2) = thin_shell_rho_ballast(in);
        case 'modular_precast'
            bounds(2) = in.materials.modular_precast.rho_hull;
    end
end

function rho_ballast = thin_shell_rho_ballast(in)
%THIN_SHELL_RHO_BALLAST Solid ballast density; the shell density when the input leaves it unset.
    if isfield(in.materials.thin_shell, 'rho_ballast') && ~isempty(in.materials.thin_shell.rho_ballast)
        rho_ballast = in.materials.thin_shell.rho_ballast;
    else
        rho_ballast = in.materials.thin_shell.rho_shell;
    end
end

function [products, ms2_model] = compute_geometry_products(ms2_file, g)
%COMPUTE_GEOMETRY_PRODUCTS Run the expensive geometry steps of build_config.
% Reads the deck and the values in g only. products.config holds the config fields these steps
% set, products.report the minimum and maximum constructable mass, and products.warnings {id, message} rows
% raised by this function itself, re-issued when the products are reloaded. Warnings raised
% inside the library functions it calls (the hydrostatic tables, the z-crossing search, the
% midplane profile, the cap contribution) are shown on a fresh build only: lastwarn keeps one
% warning and evalc would hide the progress output of a long build. Everything in products is a
% plain array or struct; the parsed model is returned separately.

    geo = struct();
    notes = cell(0, 2);

    hull = mwecmass.driver.parse_hull_deck(ms2_file);
    geo.ms2_model = hull.ms2_model;
    geo.profile   = hull.profile;

    hull_extents_z = hull.extents_z;   % [m] [z_min, z_max] of the actual hull surfaces
    hull_z_min = hull_extents_z(1);
    hull_z_max = hull_extents_z(2);

    geo.total_wec_volume = hull.volume;     % [m^3] divergence theorem; hydrostatic tables may replace it
    geo.hull_centroid    = hull.centroid;   % [1x3] [m] body frame, z up

    % Store second volume moments for the HAMS inertia tensor.
    geo.hull_int_x2 = hull.int_x2;          % [m^5]
    geo.hull_int_y2 = hull.int_y2;          % [m^5]
    geo.hull_int_z2 = hull.int_z2;          % [m^5]

    %% Hydrostatic tables
    %  mwecmass.driver.build_hydrostatic_tables builds the waterplane-area, waterplane-inertia,
    %  perimeter, submerged-volume, centre-of-buoyancy and wetted-area tables on one adaptive
    %  z-grid, together with the boundary cache they are evaluated from. The displaced volume of
    %  record comes back with them: for a hull whose visible surfaces do not form a closed
    %  boundary the Aw-trapz maximum replaces the divergence-theorem value (that function's own
    %  override note gives the reason).

    tables = mwecmass.driver.build_hydrostatic_tables(geo.ms2_model, hull_z_min, hull_z_max, ...
                 g.aw_table_dz, geo.total_wec_volume);

    geo.boundary_cache   = tables.boundary_cache;
    geo.Aw_table_z       = tables.Aw_table_z;      % [m]
    geo.Aw_table         = tables.Aw_table;        % [m^2]
    geo.I_wp_xx_table    = tables.I_wp_xx_table;   % [m^4]
    geo.I_wp_yy_table    = tables.I_wp_yy_table;   % [m^4]
    geo.P_table          = tables.P_table;         % [m]
    geo.V_sub_table      = tables.V_sub_table;     % [m^3]
    geo.CB_z_table       = tables.CB_z_table;      % [m]
    geo.S_wet_table      = tables.S_wet_table;     % [m^2], [] without a parsed deck
    geo.total_wec_volume = tables.total_wec_volume; % [m^3] displaced volume of record

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
    %  Keep strip_edges so downstream realisation uses the same wall boundary.
    %    When constructability is enabled, downstream code
    %    (mwecmass.realise.modular_precast.solve_and_extract) needs the exact strip boundaries
    %    to define the wall/platform split.  Recomputing from node midpoints
    %    would shift the wall boundary by half a strip.

    geo.hull_z_min = hull_z_min;
    geo.hull_z_max = hull_z_max;

    if g.enable_constructability
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

        N = g.num_ballast_sections;

        wall_pos = g.wall_position;   % 'bottom' for a capsized _180 deck, 'top' otherwise
        geo.wall_position = wall_pos;

        wall_h = g.wall_height;

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
                geo.strip_edges = [platform_edges; hull_z_max];

                geo.wall_strip_index = N;   % last strip = wall

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
                geo.strip_edges = [hull_z_min; platform_edges];

                geo.wall_strip_index = 1;   % first strip = wall

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
        geo.density_nodes_z = zeros(N, 1);
        for i = 1:N
            geo.density_nodes_z(i) = 0.5 * (geo.strip_edges(i) + geo.strip_edges(i+1));
        end

        fprintf('    Density bounds: [%.0f, %.0f] kg/m^3 (platform only)\n', ...
                g.ballast_density_bounds(1), g.ballast_density_bounds(2));
    else
        % ── Default: uniform layout ──────────────────────────────
        geo.density_nodes_z = linspace(hull_z_min, hull_z_max, ...
                                          g.num_ballast_sections)';
        geo.wall_strip_index = [];
        geo.strip_edges      = [];   % recompute from midpoints downstream

        fprintf('  Ballast: %d nodes over [%.2f, %.2f] m, rho in [%.0f, %.0f] kg/m^3\n', ...
                g.num_ballast_sections, hull_z_min, hull_z_max, ...
                g.ballast_density_bounds(1), g.ballast_density_bounds(2));
    end

    %% Per-strip constructability bounds
    %
    %  Platform strips use a perpendicular offset shell of thickness t_min.
    %  The shell volume and effective-density lower bound are evaluated pointwise;
    %  a zero inner radius represents a locally solid section.
    %  The perpendicular s_max is retained for diagnostic reporting.

    if g.enable_constructability
        fprintf('  Computing per-strip density bounds from hull geometry...\n');

        N            = g.num_ballast_sections;
        rho_hull_c   = g.rho_hull;
        rho_air_c    = g.rho_air;
        t_min_c      = g.t_min;
        w_idx        = geo.wall_strip_index;
        % CONTRACT: mirror the realiser's z-sampling density.  The realiser
        % uses config.constructability_n_sub (g.n_sub here, default 100) per strip — if
        % we sample more sparsely the realiser will find a tighter s_max
        % that the optimiser bound never saw, re-opening the relaxed-vs-
        % true feasibility-set gap this whole subsystem closes.
        n_rmin_sub   = g.n_sub;

        per_strip_lb         = ones(N, 1) * g.ballast_density_bounds(1);
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

            z_lo = geo.strip_edges(i);
            z_hi = geo.strip_edges(i + 1);

            % Mirror the realiser's z-sampling EXACTLY
            % Include
            % both endpoints, with the top sample nudged 1mm below z_hi to
            % avoid the wall-strip surface that starts exactly at z_hi.
            % Including endpoints captures tight
            % shoulders at the strip boundary and made the central-
            % difference r' use different neighbours than the realiser.
            n_rmin_sub_i = max(n_rmin_sub, ceil((z_hi - z_lo) / 0.01));
            z_samples = linspace(z_lo, z_hi, n_rmin_sub_i);
            z_samples(end) = z_hi - 1e-3;

            % ── Pass 1: collect r_min_k and A_k at each z-sample ────
            %  We keep the global r_min purely as a diagnostic output
            %  (per_strip_rmin).  s_max is now computed via the
            %  PERPENDICULAR formula below — strictly tighter than the
            %  radial 1−t_min/r while accounting for the profile slope r'.
            n_k     = length(z_samples);
            r_min_k_arr  = zeros(n_k, 1);
            A_k_arr      = zeros(n_k, 1);
            r_min_i      = Inf;

            for k = 1:n_k
                [r_k, ~] = mwecmass.geometry.compute_rmin_at_z( ...
                        geo.ms2_model, z_samples(k), 100, ...
                        geo.boundary_cache);

                if r_k > 1e-10 && r_k < r_min_i
                    r_min_i = r_k;
                end

                r_min_k_arr(k) = r_k;
                A_k_arr(k) = max(0, interp1(geo.Aw_table_z, ...
                                 geo.Aw_table, z_samples(k), ...
                                 'linear', 0));
            end

            per_strip_rmin(i) = r_min_i;

            % ── Pass 2: profile slope rp(z) via central differences ──
            % Use endpoint and central finite differences for the profile slope.
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
            % Suppress the finite-difference spike at the dome tip.
            rp_k_arr(r_min_k_arr < t_min_c) = 0;

            % The exported lower bound uses the offset-shell rule computed below.
            % Pass 3's own output, s_max_strip / per_strip_smax(i), feeds only the diagnostic
            % console table below and is not read by any bound; the ρ_min_strip formula shown
            % below is therefore not the value geo.per_strip_density_lb takes.
            % ── Pass 3: per-z perpendicular s_max → strip s_max → ρ_min (diagnostic only) ──
            %  PERPENDICULAR wall thickness on a revolution profile:
            %      t_perp(z) = (1 − s)·r·sqrt(1 + r'²) / (1 + s·r'²)
            %  Setting t_perp = t_min and solving for s:
            %      s_max(z) = (r·L − t_min) / (r·L + t_min·r'²),  L = √(1+r'²)
            % Apply the perpendicular-thickness constraint.
            %  At r' = 0 this reduces to the radial 1 − t_min/r.
            %
            %  The realiser uses a UNIFORM s_i across the whole strip, capped
            %  at the GLOBAL MIN of s_max(z) over z-samples.  So this pass's diagnostic
            %  minimum achievable ρ_eff is:
            %      ρ_min_strip = ρ_hull − s_max_strip² · (ρ_hull − ρ_air)
            %  (NOT the volume-weighted local mean — that was wrong.  The
            %  realiser's actual mass formula is uniform: m_strip = V·ρ_eff,
            %  The realiser uses uniform effective strip density.  This diagnostic value has no
            %  consumer: geo.per_strip_density_lb is set below from the offset-shell rule,
            %  not from this formula.
            % Degenerate-strip guard: if no valid r samples exist (the
            % B-spline evaluator found no surface or the strip closes
            % completely), force solid (s_max=0, rho=rho_hull) to mirror
            % mwecmass.realise.modular_precast.solve_and_extract:451-453.  Without this guard
            % the loop below leaves s_max_strip = 1 and the strip's lb
            % collapses to rho_air, the OPPOSITE of the realiser's
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

            % ── Offset-shell lower density ──────────────
            % Replace the homothetic perpendicular-s_max formula with a
            % uniform-thickness offset-shell formula (mirrors the
            % the thin-shell realisation realisation model).  Strip mass at ρ_min:
            %   V_shell·ρ_hull  +  V_int·ρ_air
            % where V_shell is the perpendicular-offset shell volume at
            % t = t_min, V_int = V_strip − V_shell, and r_inner is capped
            % at zero where t·sqrt(1+r'²) ≥ r (locally solid).
            %
            % This is strictly looser than the homothetic ρ_min for any
            % strip whose worst z forces s_max(strip) → 0 — the dominant
            % bottleneck on hulls with sharp shoulders.  See
            % mwecmass.hydrostatics.compute_perpendicular_shell_volume for
            % the derivation and accuracy claim.
            V_strip_i_est = trapz(z_samples(:), A_k_arr(:));
            try
                [V_shell_at_tmin, A_outer_i, ~] = ...
                    mwecmass.hydrostatics.compute_perpendicular_shell_volume( ...
                        geo.ms2_model, z_lo, z_hi, t_min_c, ...
                        n_rmin_sub_i, geo.boundary_cache);
            catch ME
                notes(end+1, :) = {'mwecmass:driver:OffsetShellFailed', sprintf( ...
                    'Strip %d: compute_perpendicular_shell_volume failed (%s). Falling back to perpendicular-formula ρ_min.', ...
                    i, ME.message)};
                V_shell_at_tmin = V_strip_i_est * (rho_hull_c - rho_air_c) / rho_hull_c;  % ρ_min ≈ perpendicular fallback
                A_outer_i = NaN;
            end

            per_strip_V_shell(i) = V_shell_at_tmin;
            per_strip_A_outer(i) = A_outer_i;

            if V_shell_at_tmin >= V_strip_i_est
                rho_min_offset = rho_hull_c;     % strip too narrow → solid
            else
                V_int_i        = V_strip_i_est - V_shell_at_tmin;
                rho_min_offset = (V_shell_at_tmin * rho_hull_c + V_int_i * rho_air_c) / V_strip_i_est;
            end
            rho_min_offset = max(rho_air_c, min(rho_hull_c, rho_min_offset));

            per_strip_lb(i) = max(g.ballast_density_bounds(1), rho_min_offset);

            % Diagnostic: minimum strip mass at this lower bound
            per_strip_m_min(i) = per_strip_lb(i) * V_strip_i_est;
        end

        geo.per_strip_density_lb       = per_strip_lb(:)';
        % Diagnostic arrays remain local; only the density lower bound is exported.
        fprintf('    Strip   z_lo     z_hi    A_outer[m²]  V_shell[m³]  rho_min[kg/m3]   (offset-shell @ t_min)\n');
        fprintf('    %s\n', repmat('-', 1, 80));
        for i = 1:N
            if i == w_idx
                fprintf('    %-5d  %+6.3f   %+6.3f   WALL (pinned to %.0f kg/m³)\n', ...
                        i, geo.strip_edges(i), geo.strip_edges(i+1), rho_hull_c);
            else
                fprintf('    %-5d  %+6.3f   %+6.3f   %10.4f   %10.4f   %10.1f\n', ...
                        i, geo.strip_edges(i), geo.strip_edges(i+1), ...
                        per_strip_A_outer(i), per_strip_V_shell(i), per_strip_lb(i));
            end
        end
    else
        geo.per_strip_density_lb       = [];
    end

    %% Minimum constructable mass
    %
    %  Sum canonical 3-D strip volumes times their lower density bounds. build_config compares
    %  the minimum mass with full-submergence buoyancy before solving.

    m_min = 0;
    m_max = 0;
    if g.enable_constructability && ~isempty(geo.strip_edges) ...
            && ~isempty(geo.per_strip_density_lb)

        fprintf('  Computing minimum achievable mass (V_strip × per_strip_density_lb, offset-shell)...\n');
        N_strips = g.num_ballast_sections;

        for i = 1:N_strips
            z_lo_i = geo.strip_edges(i);
            z_hi_i = geo.strip_edges(i + 1);

            if z_hi_i <= z_lo_i + 1e-10
                continue;
            end

            % Strip volume via the divergence theorem, matching the realiser's compute_strip call.
            strip_i = mwecmass.hydrostatics.compute_strip( ...
                geo.ms2_model, z_lo_i, z_hi_i, ...
                struct('n_quad', 16, ...
                       'Aw_table_z', geo.Aw_table_z, ...
                       'Aw_table', geo.Aw_table));
            V_strip_i = strip_i.V;

            % Offset-shell ρ_min is uniform across the strip, so use the
            % divergence-theorem volume for m_strip_min.
            m_min = m_min + V_strip_i * geo.per_strip_density_lb(i);
            m_max = m_max + V_strip_i * g.ballast_density_bounds(2);
        end
    end

    %% Strip geometry and augmented tables
    %  mwecmass.driver.build_strip_geometry_tables evaluates the waterplane area exactly at every
    %  strip boundary, inserts those points into the z-tables, rebuilds V_sub and CB_z on the
    %  final grid and integrates each strip there. The uniform strip layout is derived inside it
    %  when the wall-pinned layout has not already fixed the edges (geo.strip_edges empty).

    strip_layout = struct();
    strip_layout.num_ballast_sections    = g.num_ballast_sections;      % [-] strip count
    strip_layout.enable_constructability = g.enable_constructability;   % wall-pinned layout
    strip_layout.strip_edges             = geo.strip_edges;               % [m], [] = derive
    strip_layout.density_nodes_z         = geo.density_nodes_z;           % [m] strip centroids

    [tables, strips] = mwecmass.driver.build_strip_geometry_tables(geo.ms2_model, tables, ...
                           strip_layout, hull_z_min, hull_z_max);

    geo.strip_edges    = strips.strip_edges;
    geo.Aw_table_z     = tables.Aw_table_z;
    geo.Aw_table       = tables.Aw_table;
    geo.I_wp_xx_table  = tables.I_wp_xx_table;
    geo.I_wp_yy_table  = tables.I_wp_yy_table;
    geo.V_sub_table    = tables.V_sub_table;
    geo.CB_z_table     = tables.CB_z_table;
    geo.S_wet_table    = tables.S_wet_table;
    geo.strip_V        = strips.strip_V;      % [m^3]
    geo.strip_CB_z     = strips.strip_CB_z;   % [m]
    geo.strip_Iyy      = strips.strip_Iyy;    % [m^5]
    geo.strip_Ixx      = strips.strip_Ixx;    % [m^5]
    geo.strip_Izz      = strips.strip_Izz;    % [m^5]

    %% Two-dimensional surrogate
    %  The effective width matches each 3-D cross-sectional area.
    %
    %  The 2D surrogate computes strip volume as:
    %
    %    V_strip_2D(z) = profile_chord(z) × y_span(z) × dz × k_vol
    %
    %  where profile_chord(z) = max(x) − min(x) on the midplane profile.
    %
    %  A raw transverse width can mismatch a non-rectangular cross-section:
    %    y_span was the raw transverse width: max(y) − min(y).  For hulls
    %    with non-rectangular cross-sections (e.g., C0's wide platform),
    %    chord × raw_width ≠ actual polygon area.  The per-strip volume
    %    distribution was wrong even when k_vol corrected the TOTAL volume.
    %    This caused mass balance infeasibility (exitflag = -2) because
    %    different strips got incorrect fractions of the total volume.
    %
    %  Use the area-to-chord ratio at each z level:
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

    % Effective width from parser extents (full hull Y-span) — fallback
    y_extents = [geo.ms2_model.extents(2), geo.ms2_model.extents(5)];
    geo.effective_width = max(y_extents) - min(y_extents);

    n_z_levels = g.n_z_levels;
    geo.y_span_z_levels = linspace(hull_z_min, hull_z_max, n_z_levels)';
    geo.y_span_table    = zeros(n_z_levels, 1);

    fprintf('  Building area-corrected Y-span table (%d levels, divergence theorem)...\n', n_z_levels);

    % Compute A(z) via thin-strip divergence theorem: A ≈ V_strip / dz.
    % This uses the trusted mwecmass.hydrostatics pipeline with correct
    % source/mirror dedup — no cross-section polygon assembly.
    h_thin = 0.5 * (hull_z_max - hull_z_min) / n_z_levels;  % half-strip thickness
    yspan_quad_opts = struct('n_quad', 12);

    for kz = 1:n_z_levels
        z_k = geo.y_span_z_levels(kz);

        % --- 3D cross-section area from thin-strip divergence theorem ---
        z_lo_k = max(hull_z_min, z_k - h_thin);
        z_hi_k = min(hull_z_max, z_k + h_thin);
        dz_k   = z_hi_k - z_lo_k;

        A_3D_k = 0;
        if dz_k > 1e-10
            strip_k = mwecmass.hydrostatics.compute_strip( ...
                          geo.ms2_model, z_lo_k, z_hi_k, yspan_quad_opts);
            if strip_k.V > 1e-10
                A_3D_k = strip_k.V / dz_k;
            end
        end

        % --- 2D profile chord at this z-level ---
        chord_k = 0;
        try
            isects = mwecmass.geometry.find_waterline_intersections( ...
                         geo.profile, z_k);
            if ~isempty(isects) && size(isects, 1) >= 2
                chord_k = max(isects(:, 1)) - min(isects(:, 1));
            end
        catch
            % Degenerate z-level — no intersection
        end

        % --- Corrected effective width ---
        if chord_k > 1e-6 && A_3D_k > 1e-10
            % Area-corrected: makes 2D strip area = 3D polygon area
            geo.y_span_table(kz) = A_3D_k / chord_k;
        elseif A_3D_k > 1e-10
            % Fallback: use A as-is (degenerate profile chord)
            geo.y_span_table(kz) = A_3D_k;
        end
    end

    % Fill gaps with nearest non-zero neighbour
    nonzero_idx = find(geo.y_span_table > 0);
    if ~isempty(nonzero_idx)
        for kz = 1:n_z_levels
            if geo.y_span_table(kz) == 0
                [~, nearest] = min(abs(nonzero_idx - kz));
                geo.y_span_table(kz) = geo.y_span_table(nonzero_idx(nearest));
            end
        end
    else
        geo.y_span_table(:) = geo.effective_width;
    end

    products.config = rmfield(geo, 'ms2_model');
    products.report.m_min_constructability = m_min;
    products.report.m_max_constructability = m_max;
    products.warnings = notes;
    ms2_model = geo.ms2_model;
end
