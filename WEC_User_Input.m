function in = WEC_User_Input()

%WEC_USER_INPUT Define the single author-editable input struct and optionally run the suite.
% Syntax: in = WEC_User_Input(); with no output, execute the configured pipeline.
% Output: in struct containing physical, geometry, material, target, solver, and site settings.
% Side effects in no-output mode: changes to the repository folder/path and calls driver.run.

in = struct();

%% ------------------------------------------------------------ physical constants ---------
in.constants.rho_water = 1025;      % kg/m^3, seawater density.
in.constants.g         = 9.80665;   % m/s^2, standard gravity used by the solver.
in.constants.wamit_L   = 1.0;       % m, WAMIT reference length.

%% ------------------------------------------------------------ input files -----------------
in.files.ms2_file = 'C1.ms2';       % filename, hull geometry deck, resolved under Input/.

%% ------------------------------------------------------------ geometry discretisation -----
in.geometry.panel_size  = 0.1;      % m, target mean BEM panel edge length.
in.geometry.wp_target_edge = 0.1;   % m, target waterplane-lid panel edge length.
in.geometry.aw_table_dz = 0.05;     % m, target z-spacing for the waterplane-area table.
in.geometry.num_ballast_sections = 5;  % count, number of ballast strips; must be positive.
in.geometry.n_z_levels  = 200;      % count, z-levels for the 2D surrogate's y-span table.
in.geometry.eff_w_floor = 0.05;     % m, effective-width floor for the 2D surrogate.

%% ------------------------------------------------------------ realisation type and materials
in.materials.realisation_type = 'thin_shell';  % one of 'preliminary' | 'thin_shell' | 'modular_precast'.

in.materials.thin_shell.rho_shell = 7500;   % kg/m^3, structural steel wall density.
in.materials.thin_shell.rho_air   = 1.2;    % kg/m^3, air density used for void volume.
in.materials.thin_shell.rho_ballast = in.materials.thin_shell.rho_shell;  % kg/m^3, density of the solid ballast below z_ballast.
assert(in.materials.thin_shell.rho_ballast >= in.materials.thin_shell.rho_shell, ...
    'WEC_User_Input:ThinShellSeedNotMonotone', ...
    ['thin_shell.rho_ballast (%g) must be >= ', ...
     'thin_shell.rho_shell (%g). The analytic z_ballast warm-start seed in ', ...
     '+thin_shell/solve.m (the weighted-cumulative inversion of ', ...
     'the two-density ballast model) integrates ', ...
     '(rho_ballast-rho_shell)*(A_o-A_i) + (rho_ballast-rho_air)*A_i, which is only guaranteed ', ...
     'non-negative (monotone cumulative integral, required for the seed''s interp1 to be ', ...
     'well-posed) when rho_ballast >= rho_shell. A polymer shell with a ', ...
     'steel ballast satisfies this; the thin-shell mode requires this ordering.'], ...
    in.materials.thin_shell.rho_ballast, in.materials.thin_shell.rho_shell);
in.materials.thin_shell.t_init          = 0.02;   % m, initial shell-thickness guess.
in.materials.thin_shell.t_min           = 0.025;  % m, minimum shell thickness for splash-zone plate.
in.materials.thin_shell.max_slope_factor = 5.0;   % dimensionless factor, thickness-taper limit.
in.materials.thin_shell.n_z_grid        = 300;    % count, z-grid resolution for the shell solve.

in.materials.modular_precast.rho_hull    = 2500;   % kg/m^3, UHPC.
in.materials.modular_precast.rho_air     = 1.2;    % kg/m^3, air density in precast voids.
in.materials.modular_precast.t_min       = 0.0762; % m, 3 inches.
in.materials.modular_precast.wall_height = 1.8;    % m, wall-region height.
in.materials.modular_precast.n_sub       = 100;    % count, z-samples per strip.
in.materials.modular_precast.t_init = in.materials.thin_shell.t_init * ...
    (in.materials.thin_shell.rho_shell / in.materials.modular_precast.rho_hull);
in.materials.modular_precast.max_slope_factor = in.materials.thin_shell.max_slope_factor;
in.materials.modular_precast.n_z_grid = in.materials.thin_shell.n_z_grid;

%% ------------------------------------------------------------ bounds ----------------------
in.bounds.ballast_density_bounds = [20, 2500];  % kg/m^3, [lo, hi] per-strip density bound.
in.bounds.max_density_ratio      = 100.0;       % dimensionless ratio; asserted >=1 downstream.
in.bounds.vertical_shift_bounds  = [];          % m, [lo, hi] draft-shift bound; [] selects auto.

%% ------------------------------------------------------------ optimisation targets and ranges
in.targets.T_heave_goal   = 7.77;         % s, heave natural-period target.
in.targets.T_pitch_goal   = 3.89;         % s, pitch natural-period target.
in.targets.T_heave_range  = [7.0, 10.0];  % s, [lo, hi] heave-period band for range penalties.
in.targets.T_pitch_range  = [3.0, 5.0];   % s, [lo, hi] pitch-period band.
in.targets.gm_min    = 0.2;   % m, GM_L constraint floor (fmincon inequality GM_L >= gm_min).
in.targets.gm_range  = [0.2, 0.7];  % m, [lo, hi] GM_L reporting/consistency band.
in.targets.gm_target = 0.5;   % m, GM_L design target (documentation/objective reference).


%% ------------------------------------------------------------ PID calibration -------------
in.pid.stage1_mode = 'sweep';   % enum, Stage-1 draft-search mode.
in.pid.k_vol_init = 1.0;    % dimensionless correction factor, initial volume-PID gain scale.
in.pid.k_gm_init  = 1.0;    % dimensionless, initial GM-PID gain scale.
in.pid.mass_gains = [0.8, 0.002, 0.05];   % [P, I, D], mass-correction PID gains.
in.pid.vol_gains  = [0.6, 0.01, 0.03];    % [P, I, D], volume-correction PID gains.
in.pid.gm_gains   = [0.3, 0.02, 0.02];    % [P, I, D], GM-correction PID gains.
in.pid.mass_limits = [0.5, 5.0];    % [lo, hi], PID output saturation limits (mass correction).
in.pid.vol_limits  = [-1.0, 1.0];   % [lo, hi], PID output saturation limits (volume).
in.pid.gm_limits   = [-1.0, 1.0];   % [lo, hi], PID output saturation limits (GM).
in.pid.bounds_kvol = [0.30, 3.00];  % [lo, hi], k_vol search bounds.
in.pid.bounds_kgm  = [0.50, 3.00];  % [lo, hi], k_gm search bounds.
in.pid.damping_vol_early = 0.5;   % dimensionless, early-iteration volume-PID damping.
in.pid.damping_vol_late  = 0.7;   % dimensionless, late-iteration volume-PID damping.
in.pid.damping_gm_early  = 0.4;   % dimensionless, early-iteration GM-PID damping.
in.pid.damping_gm_late   = 0.5;   % dimensionless, late-iteration GM-PID damping.
in.pid.damping_transition_iter = 2;   % count, iteration at which early->late damping switches.
in.pid.vol_conv_tol_pct = 3.0;    % %, volume-PID convergence tolerance.
in.pid.gm_conv_tol      = 0.05;   % m, GM-PID convergence tolerance.
in.pid.delta_kvol_stable = 0.005; % dimensionless, k_vol stability-detection delta.
in.pid.delta_kgm_stable  = 0.01;  % dimensionless, k_gm stability-detection delta.
in.pid.stable_count_needed = 2;   % count, consecutive stable iterations required.
in.pid.mass_acceptable_pct = 10;  % %, acceptable mass-balance error band.
in.pid.cg_guard_floor      = 0.01;  % m, CG-guard floor.
in.pid.sat_proximity       = 0.01;  % dimensionless, PID saturation-proximity threshold.
in.pid.tanh_shape_param    = 3.0;   % dimensionless, tanh-profile shape for the initial guess.
in.pid.tanh_shape_fallback = 1.5;   % dimensionless, fallback tanh-profile shape parameter.

%% ------------------------------------------------------------ objective shaping -----------
in.objective.zone_k_amp     = 5.0;   % dimensionless, curvature-amplification factor outside transition zones.
in.objective.penalty_guard  = 1e4;   % dimensionless, fallback penalty for failed evaluations.

%% ------------------------------------------------------------ solver settings -------------
in.solver.max_outer_iterations = 50;      % count, outer PID/optimiser iteration cap.
in.solver.n_sweep_refine       = 4;       % count, Stage-1 sweep-refinement passes.
in.solver.stage2_algorithm     = 'sqp';   % enum, fmincon Stage-2 algorithm name.

in.validation.autocad_Ixx = [];  % kg*m^2, external-CAD Ixx cross-check value; [] = unset.
in.validation.autocad_Iyy = [];  % kg*m^2.
in.validation.autocad_Izz = [];  % kg*m^2.
in.validation.autocad_discrepancy_pct = 5;  % %, allowed CAD-vs-computed inertia discrepancy.

%% ------------------------------------------------------------ BEM source ------------------
in.bem.run_HAMS_MREL = false;   % bool, BEM source switch: false = WAMIT cache, true = HAMS-MREL.

in.bem.bem_cache_file = 'C1_wamit_cache.mat';  % filename, resolved under Input/.

in.bem.hams_cache_file = 'C1_hams_cache.mat';  % filename, resolved under in.bem.hams_dir.

in.bem.water_depth = 74;   % m, BEM water depth positive downward; -1 denotes deep water.

in.bem.hams_dir = fullfile('Output', 'hams_mrel');   % path, HAMS-MREL runtime I/O workspace.
in.bem.hams_exe = fullfile('Input', 'HAMS_MREL', 'HAMS_MREL');  % path, HAMS-MREL solver binary.
in.bem.T_min = 3.0;    % s, shortest BEM analysis period.
in.bem.T_max = 20.0;   % s, longest BEM analysis period.
in.bem.T_step = 0.5;   % s, BEM period step (37 frequencies at this default range).
in.bem.n_sweep_drafts = 8;    % count, coarse draft-node sweep density.
in.bem.n_adaptive_refine = 5;  % count, extra HAMS runs in the high-gradient draft zone.


%% ------------------------------------------------------------ plotting and export ---------
in.plots.viz_mesh_N = 60;   % count, 3D visualisation mesh resolution per direction.

%% ------------------------------------------------------------ site and hydrodynamic context
in.context.T_heave_goal  = in.targets.T_heave_goal;    % s. Mirror; see in.targets above.
in.context.T_pitch_goal  = in.targets.T_pitch_goal;    % s. Mirror; see in.targets above.
in.context.T_heave_range = in.targets.T_heave_range;   % s [lo,hi]. Mirror; see in.targets above.
in.context.T_pitch_range = in.targets.T_pitch_range;   % s [lo,hi]. Mirror; see in.targets above.
in.context.gm_target     = in.targets.gm_target;       % m. Mirror; see in.targets above.
in.context.T_surge_goal  = 8.95;         % s, surge natural-period target (metadata only).
in.context.T_surge_range = [7.0, 9.0];   % s [lo, hi], documentation-only surge-period band.
in.context.WIS_station = '';       % char, optional WIS station ID such as 'ST73135'.
in.context.data_year    = '';      % string, optional WIS data year or range.

%% ------------------------------------------------------------ run the suite ----------------
if nargout == 0
    repo_root = fileparts(mfilename('fullpath'));
    cd(repo_root);
    addpath(fullfile(repo_root, 'src'));
    out = WEC_Output_Options();
    mwecmass.driver.run(in, out);
end

end
