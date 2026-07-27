%% WEC_DRIVER  Central parameter definition and execution script.
%
%  This script is the single source of truth for every tunable parameter
%  in the WEC mass-distribution optimisation pipeline.  No numerical
%  constant that affects results should be hardcoded anywhere else in the
%  suite; every kernel file reads its values from the `params` struct
%  assembled here and forwarded through `config`.
%
%  ARCHITECTURE
%  ────────────────────────────────────────────────────────────────────
%  WEC_Driver  (this file)
%    ├─ §1  Physical constants & environment
%    ├─ §2  Input files (MS2 geometry, WAMIT cases)
%    ├─ §3  Ballast configuration
%    ├─ §4  2D surrogate model parameters
%    ├─ §5  PTO / mooring configuration
%    ├─ §6  Optimisation targets (periods, GM)
%    ├─ §7  Optimisation variable bounds
%    ├─ §8  Stage 1 PID surrogate-training settings
%    ├─ §9  Objective-function shaping
%    ├─ §10 Solver settings
%    ├─ §11 AutoCAD inertia validation (optional)
%    ├─ §12 Initial-condition profile shape
%    └─ §13 Execute pipeline
%  ────────────────────────────────────────────────────────────────────
%
%  DEPENDENCIES (kernel files, called by WEC_Main_Optimizer)
%    WEC_Configuration_Builder  – merges `params` into `config`
%    WEC_Main_Optimizer          – Stage 1 + Stage 2 orchestration
%    run_2d_optimizer            – Stage 1 inner fmincon loop
%    calculate_2d_properties     – 2D strip-theory surrogate
%    calculate_3d_properties     – 3D mesh ground-truth model
%    zone_penalty                – C1-continuous period penalty
%    WEC_PID_Controller          – discrete PID with anti-windup
%    WEC_Core_Functions          – geometry & hydrostatic utilities
%    WEC_File_IO                 – WAMIT reader + coord transform
%    WEC_Visualization           – plotting & reporting
%
%  USAGE
%    >> WEC_Driver     % edit parameters below, then run
%
%  Author:  [Your name / WEC Optimisation Team]
%  Date:    2025-xx-xx
%  Version: 7.1 — portable paths, auto-detect wall position, params compliance

clear; clc; close all;

%% ═══════════════════════════════════════════════════════════════════
%%  §0  SUITE DIRECTORY (portable path detection)
%% ═══════════════════════════════════════════════════════════════════

suite_dir = fileparts(mfilename('fullpath'));
if ~isempty(suite_dir)
    cd(suite_dir);
end
fprintf('Suite directory: %s\n', pwd);

%% ═══════════════════════════════════════════════════════════════════
%%  §1  PHYSICAL CONSTANTS & ENVIRONMENT
%% ═══════════════════════════════════════════════════════════════════
%  SI units throughout.  Change only for fresh-water or scaled-model tests.
%
%  NOTE ON g:
%    HAMS-MREL (WavDynMods.f90) uses g = 9.80665 m/s² (ISO standard gravity).
%    Using 9.81 in the hydrostatic model and 9.80665 in HAMS introduces a
%    systematic 0.01% offset in C55 stiffness and mass-balance comparisons.
%    Set to 9.80665 here for consistency.  Document explicitly if overriding.
%
%  NOTE ON rho_water:
%    HAMS internal density is hardcoded at 1000 kg/m³ (fresh water).
%    The pipeline post-scales HAMS outputs to 1025 kg/m³ via RHO_WATER.
%    Do NOT change RHO_WATER without also verifying the HAMS scaling step.

params.RHO_WATER    = 1025;       % [kg/m^3]  seawater density
params.G            = 9.80665;    % [m/s^2]   ISO standard gravity (matches HAMS)
params.seabed_depth = -28;        % [m]       seabed elevation (for tether anchoring)
params.wamit_L      = 1.0;        % [m]       WAMIT reference length

%% ═══════════════════════════════════════════════════════════════════
%%  §2  INPUT FILES
%% ═══════════════════════════════════════════════════════════════════

params.ms2_file    = 'C1.ms2';
params.panel_size  = 0.1;    % [m]  target BEM panel edge length (mean).
                              %      Drives Nu, Nv via WEC_Mesh_Sizing
                              %      when the user requests [R]egenerate
                              %      of the hydro cache.  Smaller →
                              %      finer mesh, slower HAMS runs.
params.quarter_body = false;  % [-]  true = mesh only Q1 with [1,1] symmetry (4× fewer panels)

% IRFR NOTE: wp_target_edge must be COMPARABLE to the hull panel edge at the waterline,
% not dramatically finer. Hull panel edge ≈ 2π×R_wp/(4×Nu).
% For C1 (R_wp≈4.5m, Nu=20): edge ≈ 0.35m → target_edge = 0.3m gives ratio 1.2 (good).
% Setting target_edge too fine (e.g. 0.1m vs 0.59m hull edge) creates a pressure
% discontinuity at the waterplane that DEGRADES irregular frequency suppression.
params.wp_target_edge = 0.1;  % [m]  WP lid panel edge — match hull edge at waterline

% HAMS WAVE-PERIOD GRID (matches the 2D / 3D / climate-analysis bands)
%   Input to HAMS is period in seconds (Input_frequency_type=4), uniformly
%   stepped from T_min to T_max in steps of T_step.  HAMS output is in
%   omega (rad/s) — the downstream parser convention.
%
%   Resolution is period-uniform: very fine at long T (low ω, design band),
%   coarser at short T (high ω, IRFR tail).  Refine panel_size to push
%   irregular frequencies above the short-T end, not by adding more samples.
%
%   IMPORTANT: changing any of these values requires regenerating the
%   hydro cache via [R] at the [L]/[R] prompt in §13a.
params.hams_T_min  = 3.0;     % [s]  shortest analysis period
params.hams_T_max  = 20.0;    % [s]  longest analysis period
params.hams_T_step = 0.5;     % [s]  period step  (37 freqs at default range)

% VISUALISATION MESH (separate from BEM mesh)
% 3D density renders use VIZ_N×VIZ_N panels on the full body.
% 60×60 = 3600 panels (smooth rendering). Does NOT affect HAMS run time.
% To override: set params.viz_mesh_N and add the forwarding line below.
% WEC_Configuration_Builder §16 should forward: config.viz_Nu = params.viz_mesh_N;
params.viz_mesh_N   = 60;     % [-]  3D visualisation mesh resolution per direction

% --- OPTION A: HAMS hydro_table (preferred) ---
params.hydro_table_file = '';   % e.g. 'C1_hydro_table.mat'

% --- OPTION B: WAMIT .1 files (legacy) ---
% params.wamit_cases = {
%     'SR0_C1',  1.0, -0.4626
%     'SR1_C1',  0.0, -1.4626
%     'SR2_C1', -1.0, -2.4626
% };

%% ═══════════════════════════════════════════════════════════════════
%%  §3  BALLAST CONFIGURATION
%% ═══════════════════════════════════════════════════════════════════

params.num_ballast_sections    = 5;
params.ballast_density_bounds  = [20, 2500];  % [kg/m^3]
params.max_density_ratio       = 100.0;

%% ═══════════════════════════════════════════════════════════════════
%%  §3b  STEEL-FILL SOLVER (post-optimisation)
%% ═══════════════════════════════════════════════════════════════════

params.enable_steel_solve      = true;
params.rho_steel               = 7500;       % [kg/m³]  structural steel
params.rho_air                 = 100;        % [kg/m³]  air + trapped gas (small positive to avoid singularity)
params.steel_t_init            = 0.02;       % [m]  initial thickness guess
params.steel_t_min             = 0.01905;    % [m]  3/4" splash-zone heavy-duty plate (fabrication floor)
params.steel_max_slope_factor  = 5.0;
params.steel_n_z_grid          = 300;
% Note: The realisation solvers (steel and UHPC) use the same range-normalised
% objective form as Stage 2: phi(r_heave) + phi(r_pitch), where each residual
% is divided by the half-range of the corresponding target band (set in §6).
% Mass is enforced as an equality constraint (buoyancy balance) and GM as an
% inequality (GM ≥ gm_min) inside WEC_Shell_Offset.solve / solve_constructable.
% No per-target weights are required — the range normalization makes every
% term peak at 1.0 at its respective band boundary.

%% ═══════════════════════════════════════════════════════════════════
%%  §3c  CONSTRUCTABILITY MODE (optional)
%% ═══════════════════════════════════════════════════════════════════
%  MUTUAL EXCLUSION: enable_steel_solve and enable_constructability
%  cannot both be true.
%
%  constructability_n_sub: z-samples per strip for area integration.
%    100 samples gives <0.1% error on C1 taper regions; 20 is too coarse
%    near the keel where the cross-section changes rapidly over 0.01 m.

params.enable_constructability       = false;
params.constructability_rho_hull     = 2500;      % [kg/m³] UHPC
params.constructability_rho_fill     = 1.2;       % [kg/m³] air at STP
params.constructability_t_min        = 0.0762;    % [m]     3 inches
params.constructability_wall_height  = 1.8;       % [m]     wall region height
params.constructability_n_sub        = 100;        % [-]     z-samples per strip (100 recommended)

%  UHPC global-solve tuning (used by WEC_Constructable_Hull.realize).
%  These override the steel-solve defaults when constructability is enabled.
%  t_init is auto-computed as steel_t_init × (rho_steel / rho_UHPC) if not set;
%  provide explicitly here to override that heuristic.
params.uhpc_t_init           = params.steel_t_init * ...
                               (params.rho_steel / params.constructability_rho_hull);
params.uhpc_max_slope_factor = params.steel_max_slope_factor;
params.uhpc_n_z_grid         = params.steel_n_z_grid;

% Aw-table resolution: controls accuracy of V_sub / CB_z tables.
%   dz = 0.05 m → ~72 points for a 3.6 m hull (< 1% V_sub error).
%   Finer dz = higher accuracy in V_sub / CB_z interpolation.
params.aw_table_dz = 0.05;  % [m] target z-spacing for Aw table

if contains(params.ms2_file, '_180')
    params.wall_position = 'bottom';
else
    params.wall_position = 'top';
end

%% ═══════════════════════════════════════════════════════════════════
%%  §4  2D SURROGATE MODEL PARAMETERS
%% ═══════════════════════════════════════════════════════════════════
%  eff_w_floor: minimum effective width to prevent division-by-zero
%  in the 2D strip-theory surrogate.  Consumed by WEC_Configuration_Builder.

params.n_density_strips  = params.num_ballast_sections;
params.n_z_levels        = 200;
params.eff_w_floor       = 0.05;  % [m]  consumed by Configuration_Builder §4

%% ═══════════════════════════════════════════════════════════════════
%%  §5  PTO / MOORING CONFIGURATION
%% ═══════════════════════════════════════════════════════════════════

params.enable_pto_effects = 0;
params.K33_pto            = 0;
params.K55_pto            = 0;
params.K11_pto            = 0;
params.pto_angle_deg      = 45;
params.ht                 = 1.0;
params.bt                 = 0.5;

%% ═══════════════════════════════════════════════════════════════════
%%  §6  OPTIMISATION TARGETS
%% ═══════════════════════════════════════════════════════════════════

params.T_heave_goal  = 7.77;
params.T_pitch_goal  = 3.89;
params.T_surge_goal  = 8.95;

params.T_heave_range = [7.0, 10.0];
params.T_pitch_range = [3.0, 5.0];
params.T_surge_range = [7.0, 9.0];

params.gm_min    = 0.2;
params.gm_range  = [0.2, 0.7];
params.gm_target = 0.5;
%% ═══════════════════════════════════════════════════════════════════
%%  §7  OPTIMISATION VARIABLE BOUNDS
%% ═══════════════════════════════════════════════════════════════════

params.vertical_shift_bounds = [];   % auto-computed in Config Builder §8

%% ═══════════════════════════════════════════════════════════════════
%%  §8  STAGE 1 — PID SURROGATE-TRAINING SETTINGS
%% ═══════════════════════════════════════════════════════════════════

params.stage1_mode = 'sweep';

params.k_vol_init = 1.0;
params.k_gm_init  = 1.0;

params.pid_mass_gains   = [0.8,  0.002, 0.05];
params.pid_vol_gains    = [0.6,  0.01,  0.03];
params.pid_gm_gains     = [0.3,  0.02,  0.02];

params.pid_mass_limits  = [0.5,  5.0];
params.pid_vol_limits   = [-1.0, 1.0];
params.pid_gm_limits    = [-1.0, 1.0];

params.bounds_kvol = [0.30, 3.00];
params.bounds_kgm  = [0.50, 3.00];

params.damping_vol_early = 0.5;
params.damping_vol_late  = 0.7;
params.damping_gm_early  = 0.4;
params.damping_gm_late   = 0.5;
params.damping_transition_iter = 2;

params.vol_conv_tol_pct    = 3.0;
params.gm_conv_tol         = 0.05;
params.delta_kvol_stable   = 0.005;
params.delta_kgm_stable    = 0.01;
params.stable_count_needed = 2;
params.mass_acceptable_pct = 10;
params.cg_guard_floor      = 0.01;
params.pid_sat_proximity   = 0.01;

%% ═══════════════════════════════════════════════════════════════════
%%  §9  OBJECTIVE-FUNCTION SHAPING
%% ═══════════════════════════════════════════════════════════════════

params.period_delta  = 3.0;
params.zone_k_amp    = 5.0;
params.penalty_guard = 1e4;

%% ═══════════════════════════════════════════════════════════════════
%%  §10  SOLVER SETTINGS
%% ═══════════════════════════════════════════════════════════════════

params.max_outer_iterations = 50;
params.convergence_tol      = 1e-2;
params.n_sweep_refine       = 4;

%% ═══════════════════════════════════════════════════════════════════
%%  §10b  HAMS SWEEP DENSITY  (FIX L1 — moved from hardcoded §13a)
%% ═══════════════════════════════════════════════════════════════════
%  n_hams_sweep_drafts: number of draft levels for the initial HAMS coarse
%  sweep (§13a).  Each point requires one complete HAMS run (~5–30 s per
%  run depending on mesh density).
%
%  GUIDANCE
%    8  drafts: adequate for smooth A(∞) vs draft curves; ~2 min total
%   12  drafts: recommended for complex geometries with nonlinear A(∞)
%   20  drafts: overkill for most geometries; retain only if A(∞) vs draft
%               shows strong nonlinearity (R² < 0.99 on linear fit)
%
%  WHY this belongs in params (not hardcoded in §13a):
%    The sweep density directly affects the accuracy of the hydro_table
%    interpolation used throughout the optimiser loop.  Hardcoding it
%    in the execution block violates the single-source-of-truth principle.

params.n_hams_sweep_drafts    = 8;   % [-]  recommended starting value

%  n_hams_adaptive_refine: number of EXTRA HAMS runs inserted in the
%  highest |dT_heave/dvs| interval(s) after the coarse sweep.
%
%  WHY ADAPTIVE RATHER THAN A LARGER n_hams_sweep_drafts?
%    Uniform refinement spreads runs evenly across vs_bounds.  For C1-class
%    hulls, the T_heave transition (12 s → 3 s) is compressed into ~0.35 m
%    of vertical shift at the platform-column shoulder.  Uniform doubling
%    from 8→16 adds at most 2 extra points there; adaptive refinement
%    concentrates all N extra runs in exactly that gap.
%
%  GUIDANCE
%    0  — disabled (coarse sweep only; acceptable if geometry is convex)
%    3  — one bisection per top-3 gradient intervals; typical default
%    5  — recommended for C1-class with stepped hull cross-section
%   10  — use only if T_heave is still non-monotone after 5-run refinement
%
%  RUNTIME: each refinement point = 1 full HAMS run (~5–30 s on mesh
%           of 400–800 panels).  5 extra runs ≈ 2–3 min overhead.
%
%  NOTE: refinement runs only during fresh cache generation (choose 'R'
%        at the cache prompt).  To add refinement to an existing cache,
%        either delete <ms2_name>_hydro_cache.mat and re-run, or
%        temporarily set this to 0 to skip and load the old cache.

params.n_hams_adaptive_refine = 5;   % [-]  extra runs in high-gradient zone

%% ═══════════════════════════════════════════════════════════════════
%%  §11  AUTOCAD INERTIA VALIDATION (optional)
%% ═══════════════════════════════════════════════════════════════════

params.autocad_Ixx = [];
params.autocad_Iyy = [];
params.autocad_Izz = [];
params.autocad_discrepancy_pct = 5;

%% ═══════════════════════════════════════════════════════════════════
%%  §12  INITIAL-CONDITION PROFILE SHAPE
%% ═══════════════════════════════════════════════════════════════════

params.tanh_shape_param     = 3.0;
params.tanh_shape_fallback  = 1.5;

%% ═══════════════════════════════════════════════════════════════════
%%  §13  EXECUTE PIPELINE
%% ═══════════════════════════════════════════════════════════════════

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════╗\n');
fprintf('║     WEC OPTIMISATION PIPELINE — DRIVER v7.1      ║\n');
fprintf('╚══════════════════════════════════════════════════╝\n');
fprintf('  MS2 file    : %s\n',   params.ms2_file);
fprintf('  Wall pos.   : %s\n',   params.wall_position);

%% §13a  HYDRODYNAMIC DATA — CHECK / GENERATE  ─────────────────────
%  FIX (L1): n_coarse now reads from params.n_hams_sweep_drafts.
%  FIX (L5): HAMS directory now resolved via OS detection.

[~, ms2_name] = fileparts(params.ms2_file);
cache_file       = [ms2_name '_hydro_cache.mat'];   % [R]egenerate target (HAMS sweep writes here)
wamit_cache_file = [ms2_name '_wamit_cache.mat'];   % [L]oad source — WAMIT-derived hydro_table
                                                    % (same schema as the HAMS cache: top-level
                                                    %  struct `hydro_table` with drafts, A_inf,
                                                    %  B_avg, A, B, Fe, omega, z_cg, V_sub, mass,
                                                    %  T_band, hams_params, mesh_Nu/Nv,
                                                    %  panel_size, wp_target_edge, hams_T_*).
                                                    %  After a successful [L]oad cache_file is
                                                    %  reassigned to this path so all downstream
                                                    %  references (params.hydro_table_file, the
                                                    %  post-run replot, etc.) follow the file
                                                    %  that was actually consumed.

% OS-aware HAMS directory resolution (FIX L5)
%
%   hams_dir  — per-OS I/O workspace (Input/, Output/).  Mirrored layout
%               under HAMS_MREL/<OS>/ so Windows and Linux runs don't
%               clobber each other's files.
%   hams_exe  — absolute path to the HAMS-MREL solver binary.  On Windows
%               this lives inside HAMS_MREL/Windows alongside the
%               workspace; on Linux it's the Fedora-compiled binary built
%               in HAMS-MREL_Fedora/src/ (separate repo clone).
if ispc
    hams_sub = 'Windows';
    hams_dir = fullfile(pwd, 'HAMS_MREL', hams_sub);
    hams_exe = fullfile(hams_dir, 'HAMS_MREL.exe');
elseif ismac
    hams_sub = 'macOS';
    hams_dir = fullfile(pwd, 'HAMS_MREL', hams_sub);
    hams_exe = fullfile(hams_dir, 'HAMS_MREL');
else
    hams_sub = 'Linux';
    hams_dir = fullfile(pwd, 'HAMS_MREL', hams_sub);
    hams_exe = fullfile(pwd, 'HAMS-MREL_Fedora', 'src', 'HAMS_MREL');
end

if ~exist(hams_exe, 'file')
    warning('WEC:HAMSNotFound', ...
            'HAMS executable not found: %s\nFalling back to WAMIT if available.', hams_exe);
    hams_available = false;
else
    hams_available = true;
end

if hams_available
    % [L]oad path now consumes the WAMIT-derived cache; [R]egenerate still
    % writes the HAMS cache (cache_file).  Probe the WAMIT cache first so
    % the prompt reflects what L will actually do.
    if exist(wamit_cache_file, 'file')
        info = dir(wamit_cache_file);
        loaded_check = load(wamit_cache_file, 'hydro_table');
        n_cached = length(loaded_check.hydro_table.drafts);
        fprintf('  Hydro cache : %s (%.1f KB, %d entries, %s)  [WAMIT]\n', ...
                wamit_cache_file, info.bytes/1024, n_cached, info.date);

        if isfield(loaded_check.hydro_table, 'ms2_file') && ...
                ~isempty(loaded_check.hydro_table.ms2_file)
            if ~strcmp(loaded_check.hydro_table.ms2_file, params.ms2_file)
                fprintf('  WARNING: Cache was built for %s, current is %s\n', ...
                        loaded_check.hydro_table.ms2_file, params.ms2_file);
            end
        end

        choice = input('  [L]oad existing (WAMIT cache) / [R]egenerate (HAMS sweep) ? ', 's');
        if isempty(choice); choice = 'L'; end
    elseif exist(cache_file, 'file')
        % WAMIT cache absent but a HAMS cache is on disk.  L is configured
        % to load the WAMIT cache, so we cannot honour it — force R.
        fprintf(['  WAMIT cache not found: %s\n' ...
                 '  (HAMS cache %s is present but [L]oad now targets ' ...
                 'the WAMIT cache.)\n  Forcing [R]egenerate.\n'], ...
                wamit_cache_file, cache_file);
        choice = 'R';
    else
        fprintf('  Hydro cache : not found. Will generate.\n');
        choice = 'R';
    end

    if upper(choice(1)) == 'R'
        % ── Mesh sizing preprocess (only on [R]) ──────────────
        %   Loads the .ms2 geometry once, measures per-surface arc
        %   lengths, and picks global (Nu, Nv) so the worst-case
        %   surface has mean panel edge ≈ params.panel_size and mean
        %   aspect ratio ≈ 1.  See WEC_Mesh_Sizing.m.  The chosen
        %   values are written back to params and into the hydro
        %   cache so a later [L]oad restores them for the diagnostic
        %   mesh without re-running this step.
        parser_pre = WEC_MS2_Parser.parse(params.ms2_file);
        [Nu, Nv, sizing_report] = WEC_Mesh_Sizing.compute( ...
            parser_pre, params.panel_size);
        if Nu > 40 || Nv > 40
            fprintf(['  WARNING: Computed (Nu, Nv) = (%d, %d) ' ...
                     'exceeds the 40-per-direction soft cap.\n'], Nu, Nv);
            fprintf('           Total panels per source surface: %d.\n', (Nu-1)*(Nv-1));
            fprintf(['           HAMS runtime grows ~ O(N_panels^2).  ' ...
                    'Consider a larger params.panel_size.\n']);
            resp = input('  Continue with this mesh density? [y/N] ', 's');
            if isempty(resp) || ~strcmpi(resp(1), 'y')
                error('WEC:MeshSizingAborted', ...
                    'Aborted at user request — increase params.panel_size and re-run.');
            end
        end
        params.mesh_Nu        = Nu;
        params.mesh_Nv        = Nv;
        params.wp_target_edge = sizing_report.wp_target_edge;
        fprintf('  WP target edge   : %.3f m (matched to panel_size)\n', ...
                params.wp_target_edge);

        fprintf('\n  Building geometry config for HAMS sweep...\n');

        params_geo = params;
        params_geo.hydro_table_file = '';
        geo_config = WEC_Configuration_Builder(params_geo);

        vs_bounds = geo_config.vertical_shift_bounds;

        % FIX (L1): n_coarse comes from params, not a hardcoded literal.
        n_coarse  = params.n_hams_sweep_drafts;
        vs_coarse = linspace(vs_bounds(1), vs_bounds(2), n_coarse);

        hams_opts = struct('wp_target_edge', params.wp_target_edge, ...
                           'quarter_body', params.quarter_body);

        fprintf('\n  Running HAMS at %d drafts: [', n_coarse);
        fprintf('%.2f ', vs_coarse);
        fprintf(']\n');

        hydro_cache = HAMS_Pipeline.empty_hydro_cache();
        hydro_cache.ms2_file = params.ms2_file;
        if exist(params.ms2_file, 'file')
            ms2_info = dir(params.ms2_file);
            hydro_cache.ms2_date = ms2_info.date;
        end

        for k = 1:n_coarse
            fprintf('\n  [%d/%d] vertical_shift = %+.4f m\n', k, n_coarse, vs_coarse(k));
            [~, hydro_cache] = HAMS_Pipeline.get_or_run_hams( ...
                vs_coarse(k), hydro_cache, geo_config, hams_dir, hams_exe, 0.01, hams_opts);
        end

        % Stash mesh sizing so [L]oad can restore Nu/Nv (and the matched
        % WP lid edge) for the post-run diagnostic mesh without
        % re-parsing the .ms2 file.  Also stash the period grid so [L]oad
        % can warn if the user changed it (cache becomes stale).
        hydro_cache.mesh_Nu        = params.mesh_Nu;
        hydro_cache.mesh_Nv        = params.mesh_Nv;
        hydro_cache.panel_size     = params.panel_size;
        hydro_cache.wp_target_edge = params.wp_target_edge;
        hydro_cache.hams_T_min     = params.hams_T_min;
        hydro_cache.hams_T_max     = params.hams_T_max;
        hydro_cache.hams_T_step    = params.hams_T_step;

        hydro_table = hydro_cache; %#ok<NASGU>
        save(cache_file, 'hydro_table', '-v7.3');
        fprintf('\n  Hydro cache saved: %s (%d entries)\n', cache_file, ...
                length(hydro_cache.drafts));

        %% §13a-refine  ADAPTIVE T_HEAVE REFINEMENT  ─────────────────
        %  After the coarse sweep, compute |dT_heave/dvs| between every
        %  adjacent pair of draft points, then insert up to
        %  params.n_hams_adaptive_refine extra HAMS runs at the midpoints
        %  of the highest-gradient intervals.
        %
        %  RATIONALE — why T_heave is the refinement signal:
        %    T_heave = 2π√((M+A33)/(ρgAw)).  Near a platform-column
        %    shoulder, Aw and A33 can change 5–15× within 0.3–0.5 m of
        %    vertical shift.  The resulting steep T_heave gradient means a
        %    single uniform coarse point may span T_heave = 3–13 s in one
        %    interval; the optimiser then interpolates over a cliff and
        %    misses the operating point entirely.  Targeting that interval
        %    directly (bisection) resolves A33 and Aw exactly where
        %    needed, at the cost of N extra HAMS runs only.
        %
        %  ANALOGY: identical adaptive logic to §3f-pre (Aw_table strip
        %  boundary injection) in WEC_Configuration_Builder.

        if params.n_hams_adaptive_refine > 0 && length(hydro_cache.drafts) >= 2

            fprintf('\n  ── Adaptive T_heave refinement (up to %d extra run(s)) ──\n', ...
                    params.n_hams_adaptive_refine);

            % ── Step 1: compute T_heave at every cached draft ─────────
            n_cached  = length(hydro_cache.drafts);
            vs_all    = hydro_cache.drafts(:);        % vertical shift [m]
            T_h_all   = zeros(n_cached, 1);

            for kk = 1:n_cached
                vs_kk   = vs_all(kk);
                z_wl_kk = -vs_kk;   % waterline in body frame

                % Heave added mass at infinite frequency (DOF 3 = heave)
                A33_kk = hydro_cache.A_inf{kk}(3, 3);

                % Aw and V_sub from geometry tables (shared z-axis)
                Aw_kk   = max(0, interp1(geo_config.Aw_table_z, ...
                                         geo_config.Aw_table, ...
                                         z_wl_kk, 'linear', 0));
                Vsub_kk = max(0, interp1(geo_config.Aw_table_z, ...
                                         geo_config.V_sub_table, ...
                                         z_wl_kk, 'linear', 0));
                M_kk = Vsub_kk * params.RHO_WATER;   % hydrostatic mass

                if Aw_kk > 1e-6 && (M_kk + A33_kk) > 0
                    K33_kk      = params.RHO_WATER * params.G * Aw_kk;
                    T_h_all(kk) = 2*pi * sqrt((M_kk + A33_kk) / K33_kk);
                else
                    T_h_all(kk) = Inf;   % dry or degenerate draft
                end
            end

            % Sort by vertical shift before computing interval gradients
            [vs_srt, srt_idx] = sort(vs_all);
            T_h_srt = T_h_all(srt_idx);

            fprintf('  Coarse T_heave (sorted by draft):\n');
            fprintf('    vs=%+.3f m → T_h=%.2f s\n', ...
                    [vs_srt, T_h_srt]');

            % ── Step 2: |dT_heave/dvs| per interval ───────────────────
            n_iv   = length(vs_srt) - 1;
            grad_T = zeros(n_iv, 1);
            for kk = 1:n_iv
                dvs = vs_srt(kk+1) - vs_srt(kk);
                dT  = T_h_srt(kk+1) - T_h_srt(kk);
                if dvs > 1e-6 && isfinite(dT)
                    grad_T(kk) = abs(dT / dvs);
                end
            end

            fprintf('  |dT/dvs| per interval: ');
            fprintf('%.1f ', grad_T);
            fprintf('s/m\n');

            % ── Step 3: bisect top-gradient intervals ──────────────────
            [~, grad_order] = sort(grad_T, 'descend');
            n_done = 0;

            for ki = 1:n_iv
                if n_done >= params.n_hams_adaptive_refine; break; end

                iv     = grad_order(ki);
                g_iv   = grad_T(iv);

                % Threshold: skip intervals with negligible gradient
                % (< 1 s/m → adding a midpoint changes T_heave by < 0.1 s
                %  for a typical 0.1 m bisection step — not worth a HAMS run)
                if g_iv < 1.0; break; end

                vs_mid = 0.5 * (vs_srt(iv) + vs_srt(iv+1));

                % Skip if a point already exists within tolerance
                if any(abs(hydro_cache.drafts - vs_mid) < 0.01)
                    fprintf('  [skip] vs=%+.4f m already in cache.\n', vs_mid);
                    continue;
                end

                fprintf('\n  [refine %d/%d]  vs=%+.4f m  (|dT/dvs|=%.1f s/m,', ...
                        n_done+1, params.n_hams_adaptive_refine, vs_mid, g_iv);
                fprintf('  interval [%+.3f, %+.3f] m)\n', ...
                        vs_srt(iv), vs_srt(iv+1));

                [~, hydro_cache] = HAMS_Pipeline.get_or_run_hams( ...
                    vs_mid, hydro_cache, geo_config, ...
                    hams_dir, hams_exe, 0.01, hams_opts);
                n_done = n_done + 1;
            end

            % ── Step 4: persist updated cache ─────────────────────────
            if n_done > 0
                hydro_table = hydro_cache;
                save(cache_file, 'hydro_table', '-v7.3');
                fprintf('\n  Hydro cache updated: %d entries total (%d refinement run(s)).\n', ...
                        length(hydro_cache.drafts), n_done);
            else
                fprintf('  No refinement needed — all interval gradients < 1 s/m.\n');
            end

        end  % adaptive refinement block
        %% ────────────────────────────────────────────────────────────

        any_valid = false;
        for kk = 1:length(hydro_cache.A_inf)
            if max(abs(hydro_cache.A_inf{kk}(:))) > 1e-6
                any_valid = true; break;
            end
        end
        if ~any_valid
            error('WEC:AllHAMSFailed', ...
                  ['All %d HAMS runs failed — A(inf) is zero everywhere.\n', ...
                   'Most likely cause: Hydrostatic.in format mismatch.\n', ...
                   'Delete %s and re-run after fixing.'], ...
                  length(hydro_cache.drafts), cache_file);
        end
    else
        fprintf('  Loading existing WAMIT hydro cache.\n');
        % Redirect every downstream consumer (params.hydro_table_file,
        % the post-run replot at §13b, restore_mesh_sizing_from_cache, etc.)
        % to the WAMIT cache by reassigning cache_file.  Schema is identical
        % to the HAMS cache (same `hydro_table` struct and field set), so
        % WEC_Configuration_Builder, the 6×6→3×3 congruence transform in
        % WEC_File_IO, and WEC_Visualization.plot_hydrodynamics consume it
        % unchanged.
        cache_file = wamit_cache_file;
        [params.mesh_Nu, params.mesh_Nv, wp_edge_from_cache] = ...
            restore_mesh_sizing_from_cache(loaded_check.hydro_table, params);
        if ~isempty(wp_edge_from_cache)
            params.wp_target_edge = wp_edge_from_cache;
        end
    end

    params.hydro_table_file  = cache_file;
    params.hams_dir          = hams_dir;
    params.hams_exe          = hams_exe;
    params.hydro_cache_file  = cache_file;

    fprintf('  Hydro data  : %s\n', cache_file);
elseif exist(wamit_cache_file, 'file') || exist(cache_file, 'file')
    % HAMS solver not installed but a cache exists — load and proceed
    % without prompting.  Prefer the WAMIT cache (matching the [L]oad
    % policy enforced above); fall back to the HAMS cache only if WAMIT
    % is absent.
    if exist(wamit_cache_file, 'file')
        cache_file = wamit_cache_file;
    end
    info         = dir(cache_file);
    loaded_check = load(cache_file, 'hydro_table');
    n_cached     = length(loaded_check.hydro_table.drafts);
    fprintf('  HAMS executable not installed — loading existing hydro cache.\n');
    fprintf('  Hydro cache : %s (%.1f KB, %d entries, %s)\n', ...
            cache_file, info.bytes/1024, n_cached, info.date);
    [params.mesh_Nu, params.mesh_Nv, wp_edge_from_cache] = ...
        restore_mesh_sizing_from_cache(loaded_check.hydro_table, params);
    if ~isempty(wp_edge_from_cache)
        params.wp_target_edge = wp_edge_from_cache;
    end
    params.hydro_table_file  = cache_file;
    params.hams_dir          = hams_dir;
    params.hams_exe          = hams_exe;
    params.hydro_cache_file  = cache_file;
    fprintf('  Hydro data  : cache (HAMS exe not installed) (%s)\n', cache_file);
elseif isfield(params, 'wamit_cases') && ~isempty(params.wamit_cases)
    fprintf('  Hydro data  : WAMIT (%d cases)\n', size(params.wamit_cases, 1));
else
    error('WEC:NoHydroData', ...
          ['No hydrodynamic data available.\n', ...
           'Either install HAMS at %s\n', ...
           'or uncomment params.wamit_cases in §2.'], hams_exe);
end

fprintf('  Stage 1 mode: %s\n',   params.stage1_mode);
fprintf('  PTO mode    : %s\n\n', ternary(params.enable_pto_effects, ...
                                          'PTO-augmented', 'Free-floating'));

try
    [results, final_props] = WEC_Main_Optimizer(params);

    WEC_Diagnostics.full_report(results, final_props);

    % Use optimiser draft for HAMS mesh diagnostic (mesh was panelised there)
    if isfield(results, 'Final3D')
        fp_for_mesh = results.Final3D;
    else
        fp_for_mesh = final_props;
    end

    if hams_available && exist(cache_file, 'file')
        loaded_cache = load(cache_file, 'hydro_table');
        config_final = results.config;
        WEC_Visualization.plot_hydrodynamics(loaded_cache.hydro_table, ...
                           fp_for_mesh.vertical_shift, ...
                           final_props, config_final);
    end

    if hams_available && isfield(params, 'ms2_file') && exist(params.ms2_file, 'file')
        try
            parser_viz = WEC_MS2_Parser.parse(params.ms2_file);
            draft_viz  = fp_for_mesh.vertical_shift;
            % VIZ-MATCH FIX: run_single_hams hardcodes quarter_body=false (full
            % body, ISX=0 ISY=0) for every HAMS run regardless of params.quarter_body.
            % Force the same here so the diagnostic plot shows exactly the mesh
            % that was handed to HAMS — not a quarter-body approximation.
            pan_opts_viz = struct('trim_wl', true, 'close_gaps', false, ...
                                  'cosine_spacing', false, 'verbose', false, ...
                                  'quarter_body', false);
            mesh_viz = WEC_Panelizer.generate(parser_viz, draft_viz, ...
                           params.mesh_Nu, params.mesh_Nv, pan_opts_viz);
            config_viz = struct('wp_target_edge', params.wp_target_edge);

            % BEM mesh diagnostic: panels, waterplane lid, waterline
            WEC_Visualization.plot_mesh_diagnostic(mesh_viz, config_viz);

            % Panel normals diagnostic: quiver plot + n_z histogram
            WEC_Visualization.plot_panel_normals(mesh_viz);

        catch ME
            warning('WEC:MeshPlotFailed', 'Mesh/normal plot failed: %s', ME.message);
        end
    end

    if strcmp(params.stage1_mode, 'skip') && ...
            isfield(results.stage1_2d, 'convergence_data') && ...
            isfield(results.stage1_2d.convergence_data, 'sweep')
        WEC_Visualization.plot_draft_landscape( ...
            results.stage1_2d.convergence_data.sweep, results.config);
    end

catch ME
    print_error(ME, params);
    rethrow(ME);
end

fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf('║   SCRIPT COMPLETED                               ║\n');
fprintf('╚══════════════════════════════════════════════════╝\n\n');

%% ═══════════════════════════════════════════════════════════════════
%%  LOCAL FUNCTIONS
%% ═══════════════════════════════════════════════════════════════════
%  DEAD CODE NOTES (D1, D2):
%    print_summary and print_physics_diagnostic are never called in the
%    execution flow above.  WEC_Diagnostics.full_report is the active
%    reporting path.  These functions are retained for standalone debugging
%    but are NOT wired into the pipeline.

function print_summary(results, fp) %#ok<DEFNU>  % DEAD IN PRODUCTION (D1)
% PRINT_SUMMARY  Console report — NOT called from pipeline; use
% WEC_Diagnostics.full_report instead.

    fprintf('\n╔══════════════════════════════════════════════════╗\n');
    fprintf('║   OPTIMISATION COMPLETED SUCCESSFULLY            ║\n');
    fprintf('╚══════════════════════════════════════════════════╝\n\n');

    fprintf('  Optimal draft:       %8.3f m\n',  fp.vertical_shift);
    fprintf('  Total mass:          %8.1f kg\n', fp.mass_total);
    fprintf('  Buoyant force:       %8.1f kg\n', fp.mass_buoyant_force);
    fprintf('  Mass error:          %8.2f%%\n',  ...
            abs(fp.mass_discrepancy / fp.mass_total) * 100);
    fprintf('  GM_L:                %8.3f m\n',  fp.GM_L);
    fprintf('  T_heave:             %8.2f s\n',  fp.periods.heave);
    fprintf('  T_pitch:             %8.2f s\n',  fp.periods.pitch);

    if isinf(fp.periods.surge)
        fprintf('  T_surge:                  Inf s  (free-floating)\n');
    else
        fprintf('  T_surge:             %8.2f s\n', fp.periods.surge);
    end

    fprintf('\n  Stage 1 iterations:  %d  (converged: %s)\n', ...
            results.stage1_2d.iterations, ...
            ternary(results.stage1_2d.converged, 'yes', 'no'));
    fprintf('  Stage 2 exit flag:   %d\n',  results.stage2_3d.exitflag);
    fprintf('  Stage 2 cost:        %.4f\n', results.stage2_3d.fval);
    fprintf('  Total time:          %.2f s\n', results.optimization_time);
end

function print_physics_diagnostic(results, fp) %#ok<DEFNU>  % DEAD IN PRODUCTION (D2)
% PRINT_PHYSICS_DIAGNOSTIC  CG-vs-CB check — NOT called from pipeline.

    config = results.config;

    fprintf('\n╔══════════════════════════════════════════════════╗\n');
    fprintf('║   PHYSICS DIAGNOSTIC: BALLAST DISTRIBUTION       ║\n');
    fprintf('╚══════════════════════════════════════════════════╝\n\n');

    fprintf('  CG_z:  %+.4f m\n', fp.CG_total(3));
    fprintf('  CB_z:  %+.4f m\n', fp.CB(3));
    fprintf('  GM_L:   %.4f m\n', fp.GM_L);

    if fp.CG_total(3) < fp.CB(3)
        fprintf('  CG is %.4f m below CB  (correct)\n', fp.CB(3) - fp.CG_total(3));
    else
        fprintf('  WARNING: CG is %.4f m ABOVE CB  (top-heavy)\n', ...
                fp.CG_total(3) - fp.CB(3));
    end

    if isfield(results.stage2_3d, 'x_optimal')
        rho = results.stage2_3d.x_optimal(2:end);
        n   = length(rho);
        bot = mean(rho(1:min(3,n)));
        top = mean(rho(max(1,n-2):n));
        fprintf('\n  Bottom-3 avg density: %.0f kg/m^3\n', bot);
        fprintf('  Top-3    avg density: %.0f kg/m^3\n', top);
        if bot > top
            fprintf('  Gradient: bottom-heavy (correct)\n');
        else
            fprintf('  WARNING: top-heavy gradient (inverted!)\n');
        end
    end

    mass_err = abs(fp.mass_total - fp.mass_buoyant_force) / fp.mass_total * 100;
    fprintf('\n  Mass balance error: %.2e%%  %s\n', mass_err, ...
            ternary(mass_err < 0.1, '(satisfied)', '(violated)'));
    in_range = fp.GM_L >= config.gm_range(1) && fp.GM_L <= config.gm_range(2);
    fprintf('  GM in [%.2f, %.2f] m?  %s\n', ...
            config.gm_range(1), config.gm_range(2), ...
            ternary(in_range, 'yes', 'no'));
end

function print_error(ME, params)
    fprintf('\n  OPTIMISATION FAILED: %s\n', ME.message);
    if ~isempty(ME.stack)
        fprintf('    in %s, line %d\n', ME.stack(1).name, ME.stack(1).line);
    end
    if ~exist(params.ms2_file, 'file')
        fprintf('  MS2 file not found: %s\n', params.ms2_file);
    end
    if isfield(params, 'hydro_table_file') && ~isempty(params.hydro_table_file)
        if ~exist(params.hydro_table_file, 'file')
            fprintf('  HAMS hydro_table not found: %s\n', params.hydro_table_file);
        end
    elseif isfield(params, 'wamit_cases') && ~isempty(params.wamit_cases)
        for i = 1:size(params.wamit_cases, 1)
            f = [params.wamit_cases{i,1}, '.1'];
            if ~exist(f, 'file')
                fprintf('  WAMIT file not found: %s\n', f);
            end
        end
    end
end

function r = ternary(cond, t, f)
    if cond, r = t; else, r = f; end
end

function [Nu, Nv, wp_edge] = restore_mesh_sizing_from_cache(hydro_table, params)
    % Restore (Nu, Nv, wp_target_edge) from a loaded hydro cache for
    % downstream diagnostic mesh generation.  Falls back to (12, 12) for
    % caches built before WEC_Mesh_Sizing was introduced.
    %
    % If `params` is supplied, compares cached period-grid fields against
    % the current params.hams_T_min/max/step and warns on mismatch — the
    % cache is stale and downstream periods won't match the user's intent.
    if isfield(hydro_table, 'mesh_Nu') && isfield(hydro_table, 'mesh_Nv')
        Nu = hydro_table.mesh_Nu;
        Nv = hydro_table.mesh_Nv;
        fprintf('  Mesh sizing : Nu=%d, Nv=%d (restored from cache)\n', Nu, Nv);
    else
        Nu = 12;
        Nv = 12;
        warning('WEC:OldCacheNoMeshSize', ...
            ['Hydro cache pre-dates the panel_size sizer. ' ...
             'Using fallback (Nu, Nv) = (12, 12) for the diagnostic mesh.\n' ...
             '         Re-running with [R]egenerate will adopt params.panel_size.']);
    end
    if isfield(hydro_table, 'wp_target_edge') && ~isempty(hydro_table.wp_target_edge)
        wp_edge = hydro_table.wp_target_edge;
        fprintf('  WP target edge : %.3f m (restored from cache)\n', wp_edge);
    else
        wp_edge = [];   % signals caller to leave params.wp_target_edge unchanged
    end

    % Period-grid sanity check (best-effort — old caches lack these fields)
    if nargin >= 2 && isstruct(params)
        cache_keys   = {'hams_T_min', 'hams_T_max', 'hams_T_step'};
        cache_labels = {'T_min', 'T_max', 'T_step'};
        mismatches = {};
        for k = 1:length(cache_keys)
            ck = cache_keys{k};
            if isfield(hydro_table, ck) && ~isempty(hydro_table.(ck)) && ...
                    isfield(params, ck) && ~isempty(params.(ck))
                if abs(hydro_table.(ck) - params.(ck)) > 1e-6
                    mismatches{end+1} = sprintf( ...
                        '%s: cache=%.3f s, params=%.3f s', ...
                        cache_labels{k}, hydro_table.(ck), params.(ck)); %#ok<AGROW>
                end
            end
        end
        if ~isempty(mismatches)
            warning('WEC:CachePeriodGridMismatch', ...
                ['Loaded hydro cache uses a DIFFERENT period grid than ' ...
                 'params.hams_T_min/max/step:\n         %s\n' ...
                 '         Cached hydrodynamic coefficients are sampled on the ' ...
                 'OLD grid; downstream consumers will interpolate but cannot ' ...
                 'recover data outside the cached omega range.\n' ...
                 '         Choose [R] at the next run to regenerate on the ' ...
                 'current grid.'], strjoin(mismatches, '; '));
        end
    end
end