function results = run_MWEC_design_stage2_proto(climate_grid_path, final_props_path, opts)
%RUN_MWEC_DESIGN_STAGE2_PROTO  PROTOTYPE Stage-2 hardware-design tool.
%
%   *** SCOPE WARNING ***
%   This is NOT a Stage-1 placement tool.  It performs a joint
%   passive-feasible optimisation over (omega_n, alpha, B_PTO) and derives
%   (K_PTO, ell, K_moor) from the placement.  Output includes per-tether
%   geometry and damping/stiffness scalars -- those are Stage-2 design
%   choices.  Saved here for later inspection; not on the Stage-1
%   deliverable path.  The proper Stage-1 placement-only script is
%   run_MWEC_tuning.m.
%
%   Original docstring follows:
%
%RUN_MWEC_TUNING  Stage-1 natural-period placement target generator.
%
%   Given (i) a WIS WAM climate grid (probability-weighted, empirical 2-D
%   spectra) and (ii) a Stage-1 hull/mass-optimisation 'final_props' struct
%   plus the BEM hydro cache it was built against, this script returns
%   the three target natural periods (surge, heave, pitch) that maximise
%   broadband absorbed-power coverage of the site climate, subject to
%   user-supplied Q values and the Budal--Falnes weighting.
%
%   This script does NOT produce AEP, capture width, or absorbed power
%   in W -- by design.  Its sole deliverable is the (T_n^surge,
%   T_n^heave, T_n^pitch) placement target plus a per-period theoretical
%   absorption diagnostic.
%
%   STAGE-1 CLOSURE (Evans-optimal, no PTO design variables)
%     Evans (1976, J. Fluid Mech. 77:1, eq. 11): maximum absorbed power at
%     resonance requires B_PTO = B_rad(omega_n) -- resistance matching.
%     Applied here as a uniform scalar derived from heave:
%       B_PTO_scalar = B_22(omega_n,heave_hydro)        [N.s/m]
%       B_PTO_rot    = B_PTO_scalar * k_gyr^2           [N.m.s/rad]  for pitch
%                      k_gyr = sqrt(Iyy/M) [m]
%     The three Q scalars are then:
%       Q_heave = omega_n,h * (M+A_22(omega_n,h)) / (B_22(omega_n,h) + B_PTO_scalar)
%                 = Q_rad,heave / 2  (Evans optimal)
%       Q_pitch = omega_n,p * (Iyy+A_33(omega_n,p)) / (B_33(omega_n,p) + B_PTO_rot)
%       Q_surge = omega_c  * (M+A_11(omega_c))  / (B_11(omega_c)  + B_PTO_scalar)
%                 evaluated at energy-band centroid omega_c (K_hyd,11=0, no hydrostatic T_n)
%
%   USAGE
%     results = run_MWEC_tuning()                                  % interactive
%     results = run_MWEC_tuning(climate_grid_path, final_props_path)
%     results = run_MWEC_tuning(climate_grid_path, final_props_path, opts)
%
%   INPUTS
%     climate_grid_path  path to <station>_climate_grid.mat
%                          (schema 2.x with 'wam2d' source cells).
%     final_props_path   path to WEC_Mass_*.mat carrying 'final_props'.
%     opts (struct)      see DEFAULTS below.
%
%   DEFAULTS (override via opts):
%     .AB_model             'freq_dependent' (default) | 'constant'
%                             'constant' replaces A(omega) by A_inf and B(omega)
%                             by mean of positive-only B_kk(omega) -- useful
%                             for hulls with BEM irregularities (negative B at
%                             irregular frequencies); see computeQfromBEM.
%     .objective            'falnes' -- weighting for the placement objective.
%     .falnes_weights       [2; 1; 2]   axisymmetric deep-water Falnes [nu_s;nu_h;nu_p]
%                                        nu_k = CW_max,k / (g/w^2): monopole=1, dipole=2
%     .pct_low              0.05     lower CDF percentile for placement band
%     .pct_high             0.95     upper CDF percentile for placement band
%     .T_min                0        physical period floor (s); 0 disables
%     .plot                 true     produce diagnostic figures
%     .absorption_T_grid    2:0.1:10 (s) periods for absorption diagnostic
%     .L_ref                []       reference length for B-normalisation (m);
%                                    [] -> y-extent of final_props.cross_section.
%     .verbose              true
%
%   OUTPUT (struct)
%     .T_n_target_s         [T_surge; T_heave; T_pitch]  PRIMARY OUTPUT (s)
%     .omega_n_target       [3 x 1]  rad/s
%     .Q_used               [3 x 1]  three derived Q scalars
%     .Q_rad_at_target      [3 x 1]  BEM radiation-only Q at the target T_n's
%     .closure              struct   B_PTO_scalar, B_crit_heave, zeta_h, etc.
%     .eta_F                scalar   broadband-capture metric in [0,1]
%     .eta_per_mode         [3 x 1]  per-mode contribution to eta_F
%     .climate              struct   S_ew, omega, band, IEC summary
%     .absorption_vs_T      struct   per-T % absorption (three normalisations)
%     .bem                  struct   computeQfromBEM output
%     .opt_info             struct   placement-optimiser diagnostics
%     .meta                 struct   provenance
%
%   See also: computeQfromBEM, run_lorentzian_placement, build_S_ew_wam,
%             compute_90pct_bandwidth, lorentzian_objective.

    %% ====================================================================
    %  PATH SETUP
    %  ====================================================================
    script_dir = fileparts(mfilename('fullpath'));
    sub_dir    = fullfile(script_dir, 'MATLAB Script for Plotting');
    if exist(sub_dir, 'dir'), addpath(sub_dir); end
    if exist(script_dir, 'dir'), addpath(script_dir); end

    %% ====================================================================
    %  ARGUMENTS / FILE PICKERS
    %  ====================================================================
    if nargin < 1 || isempty(climate_grid_path)
        [fn, fp] = uigetfile( ...
            {'*_climate_grid.mat','Climate grid (*_climate_grid.mat)'; ...
             '*.mat','All .mat files'}, ...
            'Select climate grid (WAM-derived)');
        if isequal(fn, 0); fprintf('No climate grid selected.\n'); results = []; return; end
        climate_grid_path = fullfile(fp, fn);
    end

    if nargin < 2 || isempty(final_props_path)
        [fn, fp] = uigetfile( ...
            {'WEC_Mass_*.mat','Stage-1 final props (WEC_Mass_*.mat)'; ...
             '*.mat','All .mat files'}, ...
            'Select Stage-1 final_props .mat');
        if isequal(fn, 0); fprintf('No final_props selected.\n'); results = []; return; end
        final_props_path = fullfile(fp, fn);
    end

    if nargin < 3 || isempty(opts), opts = struct(); end
    opts = apply_defaults(opts);

    fprintf('\n========================================\n');
    fprintf('  MWEC TUNING -- Stage 1 placement target\n');
    fprintf('  Climate    : %s\n', short(climate_grid_path));
    fprintf('  Final props: %s   (BEM hydro_table read inline)\n', short(final_props_path));
    fprintf('  Closure    : Evans-optimal (B_PTO = B_22 at heave resonance)\n');
    fprintf('  Objective  : %s   (nu = [%.1f, %.1f, %.1f])\n', ...
            opts.objective, opts.falnes_weights(1), opts.falnes_weights(2), ...
            opts.falnes_weights(3));
    fprintf('  AB model   : %s\n', opts.AB_model);
    fprintf('========================================\n');

    %% ====================================================================
    %  STEP 1 -- LOAD CLIMATE GRID AND BUILD S_ew(omega)
    %  ====================================================================
    fprintf('\nSTEP 1: Loading climate grid and building S_ew(omega)...\n');
    cg_struct = load(climate_grid_path);
    if ~isfield(cg_struct, 'climateGrid')
        error('run_MWEC_tuning:cg', ...
              'Expected variable ''climateGrid'' in %s', climate_grid_path);
    end
    climateGrid = cg_struct.climateGrid;
    required_cg = {'probability_grid','source_grid','S_omega_grid', ...
                   'omega','Hs_centers','Te_centers'};
    for k = 1:numel(required_cg)
        if ~isfield(climateGrid, required_cg{k})
            error('run_MWEC_tuning:cgFields', ...
                  'climateGrid missing required field ''%s''.', required_cg{k});
        end
    end
    if isfield(climateGrid, 'meta') && isfield(climateGrid.meta, 'schema_version')
        fprintf('   Schema version: %s\n', climateGrid.meta.schema_version);
    end
    if isfield(climateGrid, 'meta') && isfield(climateGrid.meta, 'station_id')
        station_id = climateGrid.meta.station_id;
    else
        [~, station_id] = fileparts(climate_grid_path);
    end

    omega = climateGrid.omega(:);
    [~, S_ew, sew_meta] = build_S_ew_wam(climateGrid, omega);
    fprintf('   omega grid: [%.2f, %.2f] rad/s (N=%d)\n', ...
            omega(1), omega(end), numel(omega));
    fprintf('   S_ew m_0   = %.5f m^2  (target %.5f, err %.2f%%)\n', ...
            sew_meta.m0_Sew, sew_meta.m0_target, sew_meta.m0_err_pct);

    % --- IEC TS 62600-101 omega-form spectral moments ---
    % Per MWEC_Stage1_Theory_v2.tex eqs. (Hm0)-(eps); all in omega-form.
    IEC = compute_IEC_descriptors(omega, S_ew);
    fprintf('   IEC: Hm0=%.3f m  Te=%.2f s  T01=%.2f s  T02=%.2f s  Tp=%.2f s  eps=%.3f\n', ...
            IEC.Hm0, IEC.Te, IEC.T01, IEC.T02, IEC.Tp, IEC.eps_bw);

    %% ====================================================================
    %  STEP 2 -- 90% ENERGY BAND OF S_ew (placement bounds)
    %  ====================================================================
    fprintf('\nSTEP 2: 90%% energy band of S_ew...\n');
    [omega_L, omega_H, B90, ~, T_H, T_L] = compute_90pct_bandwidth( ...
        S_ew, omega, opts.pct_low, opts.pct_high, opts.T_min);

    %% ====================================================================
    %  STEP 2.5 -- SPECTRAL PARTITION (per-mode placement bounds)
    %  ====================================================================
    fprintf('\nSTEP 2.5: Spectral partition for per-mode placement bounds...\n');
    [omega_c, partition_omega, spectral_info] = compute_band_partition( ...
        S_ew, omega, omega_L, omega_H);
    omega_bounds = [omega_L,          partition_omega; ...  % surge: low zone
                    omega_L,          omega_H;          ...  % heave: full band
                    partition_omega,  omega_H          ];    % pitch: high zone
    if spectral_info.bimodal
        fprintf('   [BIMODAL] Valley/peak ratio=%.2f at T=%.2f s  =>  partition at T=%.2f s.\n', ...
                spectral_info.valley_ratio, 2*pi/spectral_info.valley_omega, ...
                2*pi/partition_omega);
    else
        fprintf('   Unimodal: partition at spectral centroid T=%.2f s.\n', ...
                2*pi/omega_c);
    end
    fprintf('   Surge zone : T in [%.2f, %.2f] s\n', 2*pi/partition_omega, 2*pi/omega_L);
    fprintf('   Heave zone : T in [%.2f, %.2f] s  (full band)\n', 2*pi/omega_H, 2*pi/omega_L);
    fprintf('   Pitch zone : T in [%.2f, %.2f] s\n', 2*pi/omega_H, 2*pi/partition_omega);

    %% ====================================================================
    %  STEP 3 -- BEM-DERIVED RADIATION Q PER MODE
    %  ====================================================================
    fprintf('\nSTEP 3: Computing radiation Q from BEM...\n');
    % When AB_model='constant', restrict the positive-B mean to the
    % climate's 90% energy band (avoids picking up zero-radiation tails
    % at very low or very high omega that would under-damp the average).
    if strcmpi(opts.AB_model, 'constant')
        bem_avg_band = [omega_L, omega_H];
    else
        bem_avg_band = [];
    end
    bem = computeQfromBEM(final_props_path, ...
                          'verbose',  opts.verbose, ...
                          'AB_model', opts.AB_model, ...
                          'avg_band', bem_avg_band);

    %% ====================================================================
    %  STEP 4 -- PLACEMENT + DAMPING + GEOMETRY  (objective-dependent)
    %  ====================================================================
    if strcmpi(opts.objective, 'eta_power_passive')
        fprintf('\nSTEP 4: Joint passive-feasible eta_power optimisation...\n');
        fprintf('   Variables: (omega_n_surge, omega_n_heave, omega_n_pitch, alpha, B_PTO)\n');
        fprintf('   Derived from placement: K_PTO, ell, K_moor  (matrix-structure tied)\n');
        fprintf('   Constraints: passive K_PTO>=0, mooring alpha in [%g,%g] deg, hull ell <= ell_max\n', ...
                opts.alpha_bounds(1), opts.alpha_bounds(2));
        fprintf('   Objective : eta_power = sum nu_k * 4 beta/(1+beta)^2 * Lorentzian / Falnes denom\n');

        placement_data = struct('S_ew', S_ew, 'omega', omega, 'omega_L', omega_L, ...
            'omega_H', omega_H, 'omega_bounds', omega_bounds, 'omega_c', omega_c);
        [closure, omega_n_opt, eta_opt, opt_info] = joint_passive_placement( ...
            bem, final_props_path, placement_data, opts);
        T_n_opt = 2*pi ./ omega_n_opt;
        Q_scalars = closure.Q_k;
        Q_rad_at_target = diag(bem.Q_rad_func(omega_n_opt));
        feasibility = struct('pass', true(3,1), 'margin', NaN(3,1));   % passive by construction

        % Provide legacy field aliases (some downstream code consumes these)
        closure_initial  = closure;        % no separate "initial" step in joint solve
        opt_info_initial = opt_info;

        if opts.verbose
            mode_labels = {'Surge', 'Heave', 'Pitch'};
            fprintf('\n   Mode    T_n (s)   omega_n (r/s)   beta_k     eta_match    Q_k\n');
            for kk = 1:3
                fprintf('   %-6s  %8.3f   %8.4f      %8.3f   %8.4f    %7.3f\n', ...
                    mode_labels{kk}, T_n_opt(kk), omega_n_opt(kk), ...
                    closure.beta(kk), closure.eta_match(kk), closure.Q_k(kk));
            end
            fprintf('   eta_power = %.4f   (matching-corrected broadband absorption)\n', eta_opt);
            fprintf('   K_PTO = %.4g N/m  (per-tether spring; passive -- must be >=0)\n', closure.K_PTO_scalar);
            fprintf('   K_moor = %.4g N/m (surge mooring spring; passive)\n', closure.K_moor);
            fprintf('   alpha = %.2f deg   ell = %.4f m   B_PTO = %.2f N.s/m\n', ...
                rad2deg(closure.alpha), closure.ell, closure.B_PTO_scalar);
        end

    else
        % --- Legacy 'falnes' pipeline (spectral-coverage placement) ---------
        fprintf('\nSTEP 3.5: Evans-optimal closure (legacy: B_PTO = B_22 at heave resonance)...\n');
        closure = compute_stage1_closure(bem, final_props_path, omega_c, opts.verbose);
        Q_scalars = [closure.Q_surge; closure.Q_heave; closure.Q_pitch];

        fprintf('\nSTEP 4: Lorentzian placement optimisation [Q fixed, objective=%s]...\n', opts.objective);
        [omega_n_opt, eta_opt, opt_info] = run_lorentzian_placement( ...
            S_ew, omega, omega_L, omega_H, Q_scalars, ...
            'objective',      opts.objective, ...
            'falnes_weights', opts.falnes_weights, ...
            'omega_bounds',   omega_bounds, ...
            'verbose',        opts.verbose);
        T_n_opt = 2*pi ./ omega_n_opt;

        fprintf('\nSTEP 4.5: Q self-consistency correction at placed omega_n...\n');
        omega_n_initial  = omega_n_opt;
        T_n_initial      = T_n_opt;
        Q_initial        = Q_scalars;
        eta_initial      = eta_opt;
        opt_info_initial = opt_info;
        closure_initial  = closure;
        [Q_scalars, closure] = correct_Q_at_placed(omega_n_opt, bem, closure);
        delta_Q      = Q_scalars - Q_initial;
        delta_Q_pct  = 100 * delta_Q ./ max(abs(Q_initial), eps);
        mode_labels  = {'Surge', 'Heave', 'Pitch'};
        fprintf('   Mode    Q_initial   Q_corrected   delta_Q   delta_%%\n');
        for kk = 1:3
            fprintf('   %-6s  %8.3f    %8.3f    %+8.3f   %+6.1f%%\n', ...
                    mode_labels{kk}, Q_initial(kk), Q_scalars(kk), delta_Q(kk), delta_Q_pct(kk));
        end
        fprintf('   B_PTO_scalar: %.4f -> %.4f N.s/m  (%+.1f%%)\n', ...
                closure_initial.B_PTO_scalar, closure.B_PTO_scalar, ...
                100*(closure.B_PTO_scalar / max(closure_initial.B_PTO_scalar, eps) - 1));

        [omega_n_opt, eta_opt, opt_info] = run_lorentzian_placement( ...
            S_ew, omega, omega_L, omega_H, Q_scalars, ...
            'objective',      opts.objective, ...
            'falnes_weights', opts.falnes_weights, ...
            'omega_bounds',   omega_bounds, ...
            'verbose',        false);
        T_n_opt = 2*pi ./ omega_n_opt;
        delta_T     = T_n_opt - T_n_initial;
        delta_T_pct = 100 * delta_T ./ max(T_n_initial, eps);
        fprintf('   eta_F: %.4f -> %.4f  (delta %+.4f)\n', eta_initial, eta_opt, eta_opt - eta_initial);

        Q_rad_at_target = bem.Q_rad_func(omega_n_opt);
        Q_rad_at_target = diag(Q_rad_at_target);
        feasibility = check_Q_feasibility(Q_scalars, Q_rad_at_target);
    end

    %% ====================================================================
    %  STEP 5 -- ABSORPTION-VS-T DIAGNOSTIC
    %  ====================================================================
    fprintf('\nSTEP 5: Per-period absorption diagnostic (T in [%.1f, %.1f] s)...\n', ...
            opts.absorption_T_grid(1), opts.absorption_T_grid(end));
    [L_ref, L_ref_source] = resolve_L_ref(final_props_path, opts.L_ref);
    % bem.K_hydro is already a [3x1] vector (set in computeQfromBEM.m) --
    % pass through as-is.  Sanity-check the shape.
    K_hydro_diag = bem.K_hydro(:);
    assert(numel(K_hydro_diag) == 3, ...
           'bem.K_hydro must be a 3-vector; got size %s', mat2str(size(K_hydro_diag)));
    absorption_vs_T = compute_absorption_vs_T(opts.absorption_T_grid, ...
                                              bem, closure, K_hydro_diag, ...
                                              opts.falnes_weights, L_ref);
    fprintf('   L_ref = %.3f m  (%s)\n', L_ref, L_ref_source);
    if strcmpi(opts.objective, 'eta_power_passive')
        fprintf('   Diagnostic uses FULL RAO with passive design stiffness (K_hyd + K_PTO + K_moor).\n');
        fprintf('   Placed T_n match hardware resonances (passive feasibility enforced).\n');
    else
        fprintf('   Diagnostic uses FULL RAO on the BARE hull (K_hyd only) + closure damping.\n');
        fprintf('   Hull''s actual resonances: heave at T=%.2f s, pitch at T=%.2f s\n', ...
                2*pi/closure.omega_n_h_hydro, 2*pi/closure.omega_n_p_hydro);
        fprintf('   (placed T_n are design targets, NOT used in the absorption diagnostic).\n');
    end

    %% ====================================================================
    %  STEP 6 -- HEADLINE REPORT
    %  ====================================================================
    print_headline(T_n_opt, omega_n_opt, Q_scalars, Q_rad_at_target, ...
                   eta_opt, opt_info.eta_per_mode, bem, T_L, T_H, ...
                   closure, feasibility, opts);

    %% ====================================================================
    %  STEP 7 -- DIAGNOSTIC FIGURES
    %  ====================================================================
    if opts.plot
        % Forward per-mode matching factor when running the eta_power pipeline,
        % so Fig 5 can plot the matching-corrected integrand (legacy gets ones).
        if strcmpi(opts.objective, 'eta_power_passive') && isfield(closure, 'eta_match')
            eta_match_plot = closure.eta_match(:);
        else
            eta_match_plot = ones(3, 1);
        end
        figdata = struct( ...
            'station_id',       station_id, ...
            'S_ew',             S_ew, ...
            'omega',            omega, ...
            'omega_L',          omega_L, ...
            'omega_H',          omega_H, ...
            'omega_c',          omega_c, ...
            'omega_bounds',     omega_bounds, ...
            'spectral_info',    spectral_info, ...
            'omega_n_opt',      omega_n_opt, ...
            'T_n_opt',          T_n_opt, ...
            'Q_used',           Q_scalars, ...
            'Q_rad_at_target',  Q_rad_at_target, ...
            'eta_F',            eta_opt, ...
            'eta_match',        eta_match_plot, ...
            'eta_per_mode',     opt_info.eta_per_mode, ...
            'falnes_weights',   opts.falnes_weights, ...
            'bem',              bem, ...
            'closure',          closure, ...
            'absorption_vs_T',  absorption_vs_T, ...
            'L_ref',            L_ref, ...
            'climateGrid',      climateGrid, ...
            'opts',             opts);
        try
            plot_placement_diagnostics(figdata);
        catch ME_plot
            warning('run_MWEC_tuning:plotFail', 'Plotting failed: %s', ME_plot.message);
        end
    end

    %% ====================================================================
    %  PACK OUTPUT
    %  ====================================================================
    results = struct();
    results.T_n_target_s    = T_n_opt;
    results.omega_n_target  = omega_n_opt;
    results.Q_used          = Q_scalars;
    results.Q_rad_at_target = Q_rad_at_target;
    results.closure         = closure;
    results.closure_initial = closure_initial;
    results.feasibility     = feasibility;
    % Legacy "initial vs corrected" fields -- populated only by the legacy
    % falnes pipeline.  The eta_power_passive joint solve has no separate
    % initial step (self-consistency is built in), so these become NaN.
    if exist('T_n_initial', 'var')
        results.T_n_initial      = T_n_initial;
        results.omega_n_initial  = omega_n_initial; %#ok<NODEF>
        results.Q_initial        = Q_initial;       %#ok<NODEF>
        results.eta_F_initial    = eta_initial;     %#ok<NODEF>
        results.Q_correction_pct = delta_Q_pct;     %#ok<NODEF>
        results.T_correction_pct = delta_T_pct;     %#ok<NODEF>
    else
        results.T_n_initial      = nan(3,1);
        results.omega_n_initial  = nan(3,1);
        results.Q_initial        = nan(3,1);
        results.eta_F_initial    = NaN;
        results.Q_correction_pct = nan(3,1);
        results.T_correction_pct = nan(3,1);
    end
    results.eta_F           = eta_opt;     % primary objective value (eta_power if passive pipeline)
    if strcmpi(opts.objective, 'eta_power_passive')
        results.eta_power          = eta_opt;
        results.beta_per_mode      = closure.beta;
        results.eta_match_per_mode = closure.eta_match;
        results.passive_feasible   = true;
    else
        % Keep output schema invariant: legacy 'falnes' run did not compute
        % these, so populate NaN sentinels (callers can branch on isfinite).
        results.eta_power          = NaN;
        results.beta_per_mode      = nan(3,1);
        results.eta_match_per_mode = nan(3,1);
        results.passive_feasible   = false;   % legacy did not enforce K_PTO>=0
    end
    results.eta_per_mode    = opt_info.eta_per_mode;
    results.climate         = struct( ...
        'station_id',    station_id, ...
        'S_ew',          S_ew, ...
        'omega',         omega, ...
        'omega_L',       omega_L,      'omega_H',       omega_H, ...
        'T_L',           T_L,          'T_H',           T_H, ...
        'B90',           B90, ...
        'omega_c',       omega_c, ...
        'omega_bounds',  omega_bounds, ...
        'spectral_info', spectral_info, ...
        'm0_check',      sew_meta, ...
        'IEC',           IEC);
    results.absorption_vs_T = absorption_vs_T;
    results.L_ref           = L_ref;
    results.L_ref_source    = L_ref_source;
    results.bem             = bem;
    results.opt_info        = opt_info;
    results.meta            = struct( ...
        'climate_grid_path', climate_grid_path, ...
        'final_props_path',  final_props_path, ...
        'opts',              opts);
end

% =====================================================================
function opts = apply_defaults(opts)
    if ~isfield(opts, 'AB_model'),            opts.AB_model = 'freq_dependent'; end
    % Objective: 'eta_power_passive' (default, post-2026 pipeline) joint-optimises
    % (omega_n, alpha, B_PTO) under passive K_PTO >= 0, matrix-structure B_PTO,
    % and matching-efficiency-weighted broadband absorption.  'falnes' is the
    % legacy spectral-coverage objective (peak-normalised Lorentzian) kept for
    % comparison only -- it does NOT include the per-mode 4 beta/(1+beta)^2 factor.
    if ~isfield(opts, 'objective'),           opts.objective = 'eta_power_passive'; end
    if ~isfield(opts, 'falnes_weights'),      opts.falnes_weights = [2; 1; 2]; end
    if ~isfield(opts, 'pct_low'),             opts.pct_low = 0.05; end
    if ~isfield(opts, 'pct_high'),            opts.pct_high = 0.95; end
    if ~isfield(opts, 'T_min'),               opts.T_min = 0; end
    if ~isfield(opts, 'plot'),                opts.plot = true; end
    if ~isfield(opts, 'verbose'),             opts.verbose = true; end
    % New options for eta_power_passive
    if ~isfield(opts, 'alpha_bounds'),        opts.alpha_bounds = [30, 60]; end   % deg
    if ~isfield(opts, 'B_PTO_bounds'),        opts.B_PTO_bounds = [1, 1e5]; end   % N.s/m
    if ~isfield(opts, 'nstarts'),             opts.nstarts = 16; end
    % Whether to enforce surge/pitch spectral partition zones during placement.
    % Default false for eta_power_passive (let the optimiser roam under the
    % matrix coupling) and true for legacy 'falnes' (preserves prior behaviour).
    if ~isfield(opts, 'use_partition_zones'), opts.use_partition_zones = strcmpi(opts.objective, 'falnes'); end
    % Backwards-compatibility: opts.damping_ratio_heave was removed when the
    % closure was updated to Evans-optimal (B_PTO = B_rad,heave directly).
    if isfield(opts, 'damping_ratio_heave')
        warning('run_MWEC_tuning:obsoleteOpt', ...
                ['opts.damping_ratio_heave is no longer used. ' ...
                 'Stage-1 closure now uses Evans-optimal B_PTO = B_22(omega_n,heave).']);
    end
    % Per-period absorption-vs-T diagnostic
    if ~isfield(opts, 'absorption_T_grid'),   opts.absorption_T_grid = 2:0.1:10; end
    % Reference length for normalisation (B): per-mode capture width / L_ref.
    % Default: y-extent of final_props.cross_section (col 1 span), since the
    % raft's longitudinal axis (y) faces the wave crest in this convention.
    % Set [] for auto-detect from cross_section, scalar to override.
    if ~isfield(opts, 'L_ref'),               opts.L_ref = []; end
    opts.falnes_weights = opts.falnes_weights(:);
    opts.absorption_T_grid = opts.absorption_T_grid(:);
end

% =====================================================================
function s = short(p)
    [~, base, ext] = fileparts(p);
    s = [base, ext];
end

% =====================================================================
function IEC = compute_IEC_descriptors(omega, S_ew)
%COMPUTE_IEC_DESCRIPTORS  omega-form spectral moments and IEC period
%   descriptors per MWEC_Stage1_Theory_v2.tex eqs. (Hm0)-(eps).
    omega = omega(:); S_ew = S_ew(:);
    m0  = trapz(omega, S_ew);
    m_1 = trapz(omega, S_ew ./ max(omega, 1e-12));
    m1  = trapz(omega, S_ew .* omega);
    m2  = trapz(omega, S_ew .* omega.^2);
    [~, ip] = max(S_ew);
    omega_p = omega(ip);
    IEC = struct( ...
        'm0',     m0, ...
        'm_1',    m_1, ...
        'm1',     m1, ...
        'm2',     m2, ...
        'Hm0',    4*sqrt(max(m0, 0)), ...
        'Te',     2*pi * m_1 / max(m0, eps), ...
        'T01',    2*pi * m0  / max(m1, eps), ...
        'T02',    2*pi * sqrt(max(m0, 0)/max(m2, eps)), ...
        'Tp',     2*pi / max(omega_p, eps), ...
        'omega_p',omega_p, ...
        'eps_bw', sqrt(max(0, 1 - m1.^2/(max(m0, eps)*max(m2, eps)))));
end

% =====================================================================
function closure = compute_stage1_closure(bem, final_props_path, omega_c, verbose)
%COMPUTE_STAGE1_CLOSURE  Evans-optimal B_PTO scalar + three derived Q scalars.
%
%   Evans (1976, J. Fluid Mech. 77:1, eq. 11): maximum absorbed power at
%   resonance requires B_PTO = B_rad(omega_n) -- resistance matching.
%   Applied as a single scalar derived from heave:
%     B_PTO_scalar = B_22(omega_n,heave_hydro)        [N.s/m]
%
%   PITCH UNIT CONVERSION (F1 fix)
%     B_PTO [N.s/m] cannot be added to B_33 [N.m.s/rad] directly.
%     Conversion via radius of gyration k_gyr = sqrt(Iyy/M) [m]:
%       B_PTO_rot = B_PTO_scalar * k_gyr^2            [N.m.s/rad]
%     Stage-2 replaces k_gyr with the actual PTO moment arm.
%
%   Q_SURGE (F2 fix)
%     K_hyd,11 = 0 so there is no hydrostatic natural frequency for surge.
%     Q_surge is evaluated at the energy-weighted band centroid omega_c,
%     using actual BEM A_11(omega_c) and B_11(omega_c).  B_PTO_scalar
%     dominates because B_11(omega_c) << B_22(omega_n,h) for a raft.
%
%   Q formulae:
%     Q_heave = omega_n,h * (M+A_22(omega_n,h)) / (B_22(omega_n,h) + B_PTO_scalar)
%     Q_pitch = omega_n,p * (Iyy+A_33(omega_n,p)) / (B_33(omega_n,p) + B_PTO_rot)
%     Q_surge = omega_c  * (M+A_11(omega_c))  / (B_11(omega_c)  + B_PTO_scalar)

    FP = load(final_props_path);
    fp = FP.final_props;
    M   = fp.mass_total;
    Iyy = fp.Iyy;
    K22 = fp.K_hydro(2,2);

    % --- Evans-optimal PTO scalar (F3 fix: B_22 at omega_n, not A_inf) ---
    if ~isfield(fp, 'B_full') || ~isfinite(fp.B_full(2,2)) || fp.B_full(2,2) <= 0
        error('run_MWEC_tuning:noB22', ...
              ['final_props.B_full(2,2) = %g. Evans closure requires B_22 > 0 ' ...
               'at heave hydrostatic resonance.'], fp.B_full(2,2));
    end
    B_PTO_scalar = fp.B_full(2,2);   % Evans (1976): B_PTO = B_rad,heave

    % --- Rotational equivalent for pitch (F1 fix) ---
    k_gyr     = sqrt(max(Iyy / max(M, eps), 0));   % radius of gyration [m]
    B_PTO_rot = B_PTO_scalar * k_gyr^2;            % [N.m.s/rad]

    % --- Heave at hydrostatic resonance ---
    T_n_h_hydro = fp.periods.heave;
    omega_n_h   = 2*pi / T_n_h_hydro;
    A22_h       = fp.A_full(2,2);
    B22_h       = fp.B_full(2,2);
    Q_heave     = omega_n_h * (M + A22_h) / (B22_h + B_PTO_scalar);

    % --- Pitch at hydrostatic resonance (dimensionally consistent) ---
    T_n_p_hydro = fp.periods.pitch;
    omega_n_p   = 2*pi / T_n_p_hydro;
    A33_p       = fp.A_full(3,3);
    B33_p       = fp.B_full(3,3);
    Q_pitch     = omega_n_p * (Iyy + A33_p) / (B33_p + B_PTO_rot);

    % --- Surge at band centroid omega_c (F2 fix) ---
    A11_c = max(bem.A_kk_func(1, omega_c), 0);
    B11_c = max(bem.B_kk_func(1, omega_c), 0);   % suppress negative BEM noise
    Q_surge = omega_c * (M + A11_c) / (B11_c + B_PTO_scalar);

    closure = struct( ...
        'M',              M, ...
        'Iyy',            Iyy, ...
        'K22',            K22, ...
        'k_gyr',          k_gyr, ...
        'B_PTO_scalar',   B_PTO_scalar, ...
        'B_PTO_rot',      B_PTO_rot, ...
        'omega_n_h_hydro',omega_n_h, ...
        'omega_n_p_hydro',omega_n_p, ...
        'A22_at_omh',     A22_h, ...
        'B22_at_omh',     B22_h, ...
        'A33_at_omp',     A33_p, ...
        'B33_at_omp',     B33_p, ...
        'omega_c',        omega_c, ...
        'A11_at_omc',     A11_c, ...
        'B11_at_omc',     B11_c, ...
        'Q_surge',        Q_surge, ...
        'Q_heave',        Q_heave, ...
        'Q_pitch',        Q_pitch);

    if verbose
        fprintf('   M=%.1f kg   Iyy=%.1f kg.m^2   k_gyr=%.4f m\n', M, Iyy, k_gyr);
        fprintf('   B_PTO_scalar = B_22(omega_n,h) = %.4f N.s/m  (Evans 1976)\n', B_PTO_scalar);
        fprintf('   B_PTO_rot    = B_PTO * k_gyr^2 = %.4f N.m.s/rad  (pitch, F1 fix)\n', B_PTO_rot);
        fprintf('   omega_n_heave = %.4f rad/s (T=%.3f s)  A22=%.1f  B22=%.4f\n', ...
                omega_n_h, 2*pi/omega_n_h, A22_h, B22_h);
        fprintf('   omega_n_pitch = %.4f rad/s (T=%.3f s)  A33=%.1f  B33=%.5f\n', ...
                omega_n_p, 2*pi/omega_n_p, A33_p, B33_p);
        fprintf('   omega_c (band centroid) = %.4f rad/s (T=%.3f s)\n', omega_c, 2*pi/omega_c);
        fprintf('   A_11(omega_c)=%.1f  B_11(omega_c)=%.5f  B_11/B_PTO=%.4f\n', ...
                A11_c, B11_c, B11_c/max(B_PTO_scalar, eps));
        fprintf('   Q_surge = %.3f\n', Q_surge);
        fprintf('   Q_heave = %.3f  (Q_rad,heave/2 = %.3f)\n', ...
                Q_heave, omega_n_h*(M+A22_h)/(2*B22_h));
        fprintf('   Q_pitch = %.3f\n', Q_pitch);
        if Q_heave < 2 || Q_pitch < 2
            fprintf(['   [CAVEAT] Q<2 for at least one mode: Lorentzian approximation ' ...
                     'near its validity limit (theory sec:lorentzian).\n']);
        end
    end
end

% =====================================================================
function feasibility = check_Q_feasibility(Q_used, Q_rad)
%CHECK_Q_FEASIBILITY  Assert Q_des <= Q_rad per theory eq. (Q_constraint).
%   With B_PTO_scalar >= 0 in the closure, the closed-loop Q is bounded
%   above by Q_rad mathematically; this asserts the inequality numerically
%   to catch BEM-noise pathologies (negative B at irregular frequencies).
    labels = {'surge', 'heave', 'pitch'};
    ok = true(3, 1);
    margin = nan(3, 1);
    tol = 1e-3;
    for k = 1:3
        if isfinite(Q_rad(k)) && Q_rad(k) > 0
            margin(k) = Q_rad(k) - Q_used(k);
            if Q_used(k) > Q_rad(k) * (1 + tol)
                ok(k) = false;
                warning('run_MWEC_tuning:Qfeasibility', ...
                        ['Q_used(%s)=%.3f exceeds Q_rad=%.3f at placed omega_n -- ' ...
                         'likely BEM-noise (negative B near irregular frequency). ' ...
                         'Consider opts.AB_model=''constant''.'], ...
                        labels{k}, Q_used(k), Q_rad(k));
            end
        end
    end
    feasibility = struct( ...
        'pass',   ok, ...
        'margin', margin);
end

% =====================================================================
function [L_ref, src] = resolve_L_ref(final_props_path, override)
%RESOLVE_L_REF  Reference length for capture-width normalisation.
%   Reads column 1 of final_props.cross_section [m], which stores
%   y-coordinates in the raft convention (y = longitudinal axis, parallel
%   to the incident wave crest; wave propagates in x).
%   Col-1 span = raft beam width seen by an incoming wave.
    if ~isempty(override) && isfinite(override) && override > 0
        L_ref = double(override);
        src   = 'opts.L_ref (manual)';
        return;
    end
    try
        FP = load(final_props_path);
        cs = FP.final_props.cross_section;
        L_ref = max(cs(:,1)) - min(cs(:,1));
        src   = 'final_props.cross_section col 1 (y-extent)';
    catch ME
        warning('run_MWEC_tuning:LrefFallback', ...
                'L_ref auto-detect failed (%s); using L_ref = 1.0 m.', ME.message);
        L_ref = 1.0;
        src   = 'fallback (cross_section missing)';
    end
end

% =====================================================================
function tbl = compute_absorption_vs_T(T_grid, bem, closure, K_hydro_diag, nu, L_ref)
%COMPUTE_ABSORPTION_VS_T  Per-period absorbed-power diagnostic using the
%   full BEM-based RAO of the Stage-1 hull (theory eqs. RAO, CW_formula).
%
%   TWO REGIMES, switched by closure provenance:
%     (a) eta_power_passive pipeline -- closure provides per-mode B_PTO_kk
%         (matrix-structure entries) and K_total_diag = K_hyd + K_PTO + K_moor.
%         The diagnostic then shows the FULL passive design absorbing at
%         the placed T_n (heave/pitch at hydrostatic, surge at K_moor).
%     (b) Legacy 'falnes' pipeline -- closure has only B_PTO_scalar and the
%         bare-hull K_hyd.  Heave resonates at T_n,heave,hydro; pitch at
%         T_n,pitch,hydro; surge K_hyd,11=0 so pure damped mass response.
%         The placed T_n are design targets NOT enforced here.
%
%   Per-mode capture width (theory eq. CW_formula), units: metres:
%       l_k(omega) = 2 * B_PTO * omega^3 * |Fe_k|^2 / (rho * g^2 * |D_k|^2)
%       D_k(omega) = -omega^2*(M_k+A_kk(omega))
%                     + i*omega*(B_kk(omega)+B_PTO_scalar)
%                     + K_hyd,kk
%
%   Three "% absorbed" normalisations:
%     (A) pct_A_k = l_k / sum_j(nu_j * g/omega^2) * 100
%                  -- fraction of Falnes-bounded ceiling at this frequency.
%     (B) pct_B_k = l_k / L_ref * 100   -- vs raft y-extent.
%     (C) pct_C_k = l_k / lambda * 100  -- vs wavelength.

    rho   = 1025;
    g_acc = 9.81;

    T_grid     = T_grid(:);
    omega_grid = 2*pi ./ T_grid;
    lambda_grid= g_acc .* T_grid.^2 ./ (2*pi);
    nu         = nu(:);
    K_hydro_diag = K_hydro_diag(:);

    M_diag = bem.M_diag(:);     % [M; M; Iyy]
    % Pitch unit-conversion lever arm: use body radius of gyration k_gyr unless
    % the closure provides a true k_gyr (legacy) -- the matrix structure (new
    % pipeline) uses k_gyr only as a coordinate transform, not as a PTO arm.
    if isfield(closure, 'k_gyr') && isfinite(closure.k_gyr) && closure.k_gyr > 0
        k_gyr = closure.k_gyr;
    else
        k_gyr = sqrt(max(M_diag(3) / max(M_diag(1), eps), 0));
    end
    B_PTO  = closure.B_PTO_scalar;        % legacy scalar (for backward compat)
    % Per-mode damping entries (matrix structure if available, else uniform)
    if isfield(closure, 'B_PTO_kk') && numel(closure.B_PTO_kk) == 3
        B_PTO_kk_in = closure.B_PTO_kk(:);   % [N.s/m; N.s/m; N.m.s/rad]
    else
        % Legacy: uniform scalar in surge/heave, k_gyr-converted in pitch
        B_PTO_kk_in = [B_PTO; B_PTO; B_PTO * k_gyr^2];
    end
    % Per-mode total stiffness (hydrostatic + PTO + mooring) -- passive pipeline
    % populates K_total_diag in the closure; legacy falls back to bare hull.
    if isfield(closure, 'K_total_diag') && numel(closure.K_total_diag) == 3
        K_total_diag = closure.K_total_diag(:);
    else
        K_total_diag = K_hydro_diag;        % legacy: bare hull
    end

    N = numel(T_grid);
    Akk_mat = zeros(N, 3);  Bkk_mat = zeros(N, 3);
    Fek_mat = zeros(N, 3) + 0i;
    Dk_mat  = zeros(N, 3) + 0i;
    Hk_sq   = zeros(N, 3);
    l_mode  = zeros(N, 3);

    for k = 1:3
        Ak  = bem.A_kk_func(k, omega_grid);   Ak(~isfinite(Ak)) = 0;
        Bk  = bem.B_kk_func(k, omega_grid);   Bk(~isfinite(Bk)) = 0;
        Bk(Bk < 0) = 0;
        Fek = bem.Fe_kk_func(k, omega_grid);  Fek(~isfinite(Fek)) = 0;

        % Store raw BEM values (before any pitch transformation)
        Akk_mat(:, k) = Ak;
        Bkk_mat(:, k) = Bk;
        Fek_mat(:, k) = Fek;

        if k == 3
            % Pitch: transform to equivalent translational oscillator via k_gyr.
            % Fe3 [N.m] / k_gyr -> [N]; A33/k_gyr^2 -> [kg]; B33/k_gyr^2 -> [N.s/m].
            % B_PTO,33 [N.m.s/rad] / k_gyr^2 -> [N.s/m].  Same conversion for K_33.
            % Result: l_pitch has units [m], directly comparable to l_surge/l_heave.
            k2       = k_gyr^2;
            Fek_use  = Fek / k_gyr;
            Ak_use   = Ak  / k2;
            Bk_use   = Bk  / k2;
            M_k      = M_diag(k) / k2;          % Iyy/k_gyr^2 = M_body
            K_k      = K_total_diag(k) / k2;    % rotational -> translational
            B_PTO_k  = B_PTO_kk_in(k) / k2;     % rotational -> translational
        else
            Fek_use  = Fek;
            Ak_use   = Ak;
            Bk_use   = Bk;
            M_k      = M_diag(k);
            K_k      = K_total_diag(k);
            B_PTO_k  = B_PTO_kk_in(k);
        end

        % Impedance: BARE HULL (no PTO stiffness; B_PTO is closure damping only)
        Dk = -omega_grid.^2 .* (M_k + Ak_use) ...
             + 1i .* omega_grid .* (Bk_use + B_PTO_k) ...
             + K_k;
        Dk_mat(:, k) = Dk;

        Dsq = abs(Dk).^2;
        Hk_sq(:, k)  = abs(Fek_use).^2 ./ max(Dsq, eps);
        l_mode(:, k) = 2 .* B_PTO_k .* omega_grid.^3 .* abs(Fek_use).^2 ./ ...
                       max(rho * g_acc^2 .* Dsq, eps);
    end

    falnes_total_capture = (sum(nu) .* g_acc) ./ omega_grid.^2;
    pct_A = l_mode ./ falnes_total_capture .* 100;
    pct_B = l_mode ./ L_ref .* 100;
    pct_C = l_mode ./ lambda_grid .* 100;
    pct_A_total = sum(pct_A, 2);
    pct_B_total = sum(pct_B, 2);
    pct_C_total = sum(pct_C, 2);

    J_unit = rho * g_acc^2 ./ (4 * omega_grid);

    tbl = struct( ...
        'T',           T_grid, ...
        'omega',       omega_grid, ...
        'lambda',      lambda_grid, ...
        'J_unit_W_per_m', J_unit, ...
        'A_kk',        Akk_mat, ...
        'B_kk',        Bkk_mat, ...
        'Fe_kk',       Fek_mat, ...
        'Dk',          Dk_mat, ...
        'H_sq',        Hk_sq, ...
        'l_surge_m',   l_mode(:,1), ...
        'l_heave_m',   l_mode(:,2), ...
        'l_pitch_m',   l_mode(:,3), ...
        'l_total_m',   sum(l_mode, 2), ...
        'pct_A_surge', pct_A(:,1), ...
        'pct_A_heave', pct_A(:,2), ...
        'pct_A_pitch', pct_A(:,3), ...
        'pct_A_total', pct_A_total, ...
        'pct_B_surge', pct_B(:,1), ...
        'pct_B_heave', pct_B(:,2), ...
        'pct_B_pitch', pct_B(:,3), ...
        'pct_B_total', pct_B_total, ...
        'pct_C_surge', pct_C(:,1), ...
        'pct_C_heave', pct_C(:,2), ...
        'pct_C_pitch', pct_C(:,3), ...
        'pct_C_total', pct_C_total, ...
        'B_PTO_scalar', B_PTO, ...
        'L_ref',       L_ref, ...
        'nu',          nu, ...
        'note',        ['Full-RAO capture width: bare hull + Evans B_PTO closure. ' ...
                        'Pitch transformed via k_gyr: Fe3/k_gyr [N], A33/k_gyr^2 [kg], ' ...
                        'B33/k_gyr^2 [N.s/m]; l_pitch [m] directly comparable to l_surge/l_heave.']);
end

% =====================================================================
function print_headline(T_n, omega_n, Q_used, Q_rad, eta, eta_pm, ...
                         bem, T_L, T_H, closure, feasibility, opts)
    fprintf('\n========================================\n');
    fprintf('  STAGE-1 PLACEMENT TARGETS\n');
    fprintf('  Energy band (90%%): T in [%.2f, %.2f] s\n', T_L, T_H);

    is_passive = strcmpi(opts.objective, 'eta_power_passive');
    if is_passive
        fprintf('\n  Closure (eta_power passive): joint (omega_n, alpha, B_PTO) optimisation\n');
        fprintf('   alpha    = %.2f deg\n', rad2deg(closure.alpha));
        fprintf('   ell      = %.4f m   (derived from S_3/S_2 under matrix coupling)\n', closure.ell);
        fprintf('   B_PTO    = %.3f N.s/m   (per-tether scalar)\n', closure.B_PTO_scalar);
        fprintf('   K_PTO    = %.4g N/m    (per-tether spring; PASSIVE => >=0)\n', closure.K_PTO_scalar);
        fprintf('   K_moor   = %.4g N/m    (surge mooring spring; PASSIVE => >=0)\n', closure.K_moor);
        fprintf('   B_PTO matrix diag : [%.3g, %.3g, %.3g]  [Ns/m,Ns/m,Nms/rad]\n', closure.B_PTO_kk);
        fprintf('   K_PTO matrix diag : [%.3g, %.3g, %.3g]  [N/m,N/m,Nm/rad]\n', closure.K_PTO_kk);
        fprintf('   Off-diag B_13 = %.4g  (rank-1 sub-block, unavoidable kinematic coupling)\n', closure.B13);
    else
        fprintf('\n  Closure (Evans 1976, legacy): B_PTO_scalar = B_22(T_n,heave) = %.4f N.s/m\n', ...
                closure.B_PTO_scalar);
        fprintf('           B_PTO_rot (pitch) = %.4f N.m.s/rad  (k_gyr = %.4f m)\n', ...
                closure.B_PTO_rot, closure.k_gyr);
    end

    fprintf('\n  Mode    T_n target (s)  omega_n (rad/s)   Q_used    Q_rad@target\n');
    labels = {'Surge', 'Heave', 'Pitch'};
    for k = 1:3
        if isfinite(Q_rad(k))
            qrad_str = sprintf('%9.2f', Q_rad(k));
        else
            qrad_str = '      Inf';   % surge B_11 -> 0 at long T
        end
        fprintf('  %-6s   %8.3f         %.4f         %6.3f   %s\n', ...
                labels{k}, T_n(k), omega_n(k), Q_used(k), qrad_str);
    end

    if is_passive
        fprintf('\n  Per-mode matching diagnostics  (4 beta/(1+beta)^2 = SDOF resonant capture):\n');
        fprintf('  Mode    beta_k       eta_match    contrib to eta_power\n');
        for k = 1:3
            fprintf('  %-6s  %10.3g    %8.4f      %.4f\n', labels{k}, ...
                closure.beta(k), closure.eta_match(k), eta_pm(k));
        end
        fprintf('\n  eta_power = %.4f   (matching-corrected Falnes-weighted broadband absorption)\n', eta);
    else
        fprintf('\n  Broadband capture metric eta_%s = %.4f\n', opts.objective, eta);
        fprintf('  Per-mode contribution    : surge=%.3f  heave=%.3f  pitch=%.3f\n', ...
                eta_pm(1), eta_pm(2), eta_pm(3));
    end

    fprintf('\n  Hydrostatic T_n (BEM)    : heave=%.3f s  pitch=%.3f s\n', ...
            bem.T_n_hydrostatic(2), bem.T_n_hydrostatic(3));
    fprintf('  Q_rad at hydrostatic T_n : heave=%.2f    pitch=%.2f\n', ...
            bem.Q_at_T_n(2), bem.Q_at_T_n(3));

    % Target-vs-hydrostatic Tn shift implied for Stage 2 (no PTO talk here;
    % Stage 1 only reports the gap, not the mechanism that would close it).
    for k = 2:3
        Tnt = T_n(k);  Tnh = bem.T_n_hydrostatic(k);
        if isfinite(Tnh) && Tnh > 0 && abs(Tnt - Tnh)/Tnh > 0.05
            fprintf('  [Stage-2 gap] %s: target T_n=%.2f s vs hydrostatic %.2f s (Delta = %+.2f s)\n', ...
                    labels{k}, Tnt, Tnh, Tnt - Tnh);
        end
    end

    % Implied mooring stiffness for surge target
    A11_at_surge = bem.A_kk_func(1, omega_n(1));
    if ~isfinite(A11_at_surge), A11_at_surge = 0; end
    K_moor_implied = (closure.M + max(A11_at_surge, 0)) * omega_n(1)^2;
    fprintf('\n  [Surge mooring] Target T_n=%.2f s requires K_moor = %.1f N/m\n', ...
            T_n(1), K_moor_implied);
    fprintf('                  (M + A_11(omega_n,surge) = %.1f kg)\n', ...
            closure.M + max(A11_at_surge, 0));

    % Caveats
    if Q_used(2) < 2 || Q_used(3) < 2
        fprintf(['  [CAVEAT] Q<2 for heave or pitch: Lorentzian approximation\n' ...
                 '           near its validity limit (theory sec:lorentzian).\n']);
    end

    if any(~feasibility.pass)
        fprintf('  [WARNING] Q feasibility check failed (see warnings above).\n');
    end
    fprintf('========================================\n');
end

% =====================================================================
function plot_placement_diagnostics(D)
%PLOT_PLACEMENT_DIAGNOSTICS  Five figures for Stage-1 placement output.

    labels = {'Surge', 'Heave', 'Pitch'};
    colors = {[0.85 0.30 0.30], [0.30 0.65 0.30], [0.25 0.40 0.85]};
    omega = D.omega;
    S_ew  = D.S_ew;
    mask_band = (omega >= D.omega_L) & (omega <= D.omega_H);
    T = 2*pi ./ max(omega, 1e-9);

    % ---------- Fig 1: S_ew with placed Lorentzians ----------
    f1 = figure('Name', sprintf('[%s] Fig 1: S_ew with placed modes', D.station_id), ...
                'Color', 'w', 'Position', [80 80 1100 460]);   %#ok<NASGU>

    subplot(1,2,1);
    plot(omega, S_ew, 'b-', 'LineWidth', 1.5); hold on;
    fill([omega(mask_band); flipud(omega(mask_band))], ...
         [S_ew(mask_band); zeros(sum(mask_band),1)], ...
         [0.75 0.88 1.0], 'EdgeColor','none', 'FaceAlpha',0.35, ...
         'DisplayName','90% energy band');
    for k = 1:3
        omk = D.omega_n_opt(k);
        det = omega./omk - omk./omega;
        Pk  = 1 ./ (1 + D.Q_used(k)^2 .* det.^2);
        % Scale each Lorentzian so its peak equals S_ew(omega_nk): no omega^3 factor,
        % so tails decay naturally with P_k rather than growing artificially.
        Sew_at_omk = interp1(omega, S_ew, omk, 'pchip', 0);
        plot(omega, Pk .* Sew_at_omk, '--', ...
             'Color', colors{k}, 'LineWidth', 1.2, ...
             'DisplayName', sprintf('%s P_k', labels{k}));
        xline(omk, ':', 'Color', colors{k}, 'LineWidth', 1.2, 'HandleVisibility','off');
        % T_n labels placed at staggered y-positions inside the plot to avoid
        % colliding with the top period axis added below.
        y_label_frac = [0.65, 0.55, 0.45];
        text(omk, max(S_ew)*y_label_frac(k), ...
             sprintf(' %s T_n=%.2fs', labels{k}(1), D.T_n_opt(k)), ...
             'Color', colors{k}, 'FontWeight', 'bold', ...
             'BackgroundColor', [1 1 1 0.7]);
    end
    xline(D.omega_L, 'k--', 'HandleVisibility','off');
    xline(D.omega_H, 'k--', 'HandleVisibility','off');
    if isfield(D, 'spectral_info') && D.spectral_info.bimodal
        xline(D.spectral_info.partition_omega, 'm-.', 'LineWidth', 1.0, ...
              'Label', sprintf('partition T=%.1fs', 2*pi/D.spectral_info.partition_omega), ...
              'HandleVisibility','off');
    end
    hold off;
    xlabel('\omega (rad/s)'); ylabel('S_{ew}(\omega)  (m^2 s/rad)');
    title(sprintf('S_{ew} and three-mode Lorentzian placement -- %s', D.station_id));
    legend('show', 'Location', 'northeast');
    grid on;
    add_period_axis(gca);

    % Right panel: same in period domain
    subplot(1,2,2);
    plot(T, S_ew, 'b-', 'LineWidth', 1.5); hold on;
    for k = 1:3
        xline(D.T_n_opt(k), ':', 'Color', colors{k}, 'LineWidth', 1.5, ...
              'Label', sprintf('%s %.2fs', labels{k}(1), D.T_n_opt(k)), ...
              'LabelHorizontalAlignment','center', 'HandleVisibility','off');
    end
    xline(2*pi/D.omega_L, 'k--', 'HandleVisibility','off', 'Label', 'T_H');
    xline(2*pi/D.omega_H, 'k--', 'HandleVisibility','off', 'Label', 'T_L');
    hold off;
    xlim([0, 2*(2*pi/D.omega_L)]);
    xlabel('Period T (s)'); ylabel('S_{ew}  (m^2 s/rad)');
    title('Same, period axis');
    grid on;

    % ---------- Fig 2: Q_rad(omega) per mode with target markers ----------
    f2 = figure('Name', sprintf('[%s] Fig 2: BEM Q_{rad}(\\omega) per mode', D.station_id), ...
                'Color', 'w', 'Position', [120 120 1000 420]);   %#ok<NASGU>
    om_bem = D.bem.omega_BEM;
    Q_rad  = D.bem.Q_rad;
    for k = 1:3
        semilogy(om_bem, Q_rad(k, :), '-', 'Color', colors{k}, 'LineWidth', 1.6, ...
                 'DisplayName', sprintf('%s Q_{rad}', labels{k})); hold on;
    end
    for k = 1:3
        yline(D.Q_used(k), '--', 'Color', colors{k}, 'LineWidth', 1.0, ...
              'Label', sprintf('Q_{des,%s}=%.2f', labels{k}(1), D.Q_used(k)), ...
              'HandleVisibility', 'off');
    end
    for k = 1:3
        if isfinite(D.Q_rad_at_target(k))
            plot(D.omega_n_opt(k), D.Q_rad_at_target(k), 'o', ...
                 'MarkerFaceColor', colors{k}, 'MarkerEdgeColor', 'k', ...
                 'MarkerSize', 9, 'HandleVisibility', 'off');
        end
        xline(D.omega_n_opt(k), ':', 'Color', colors{k}, 'LineWidth', 0.8, ...
              'HandleVisibility', 'off');
    end
    xline(D.omega_L, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    xline(D.omega_H, 'k--', 'LineWidth', 0.8, 'HandleVisibility','off');
    hold off;
    xlabel('\omega (rad/s)'); ylabel('Q_{rad}');
    title(sprintf('Radiation-only Q vs \\omega -- %s', D.station_id));
    legend('show', 'Location', 'best'); grid on;
    add_period_axis(gca);

    % ---------- Fig 3: BEM verification (raw A_kk and B_kk per mode) ----
    plot_bem_verification(D, labels, colors);

    % ---------- Fig 4: scatter (Hs, Te) with target T_n lines ----------
    plot_scatter_with_targets(D, labels, colors);

    % ---------- Fig 5: Falnes-weighted spectrum + stacked eta_F integrand ----
    plot_falnes_panel(D, labels, colors);

    % ---------- Fig 6: per-period absorption diagnostic (T in [2,10] s) ----
    if isfield(D, 'absorption_vs_T') && ~isempty(D.absorption_vs_T)
        plot_absorption_vs_T(D, labels, colors);
    end
end

% =====================================================================
function plot_bem_verification(D, labels, colors)
%PLOT_BEM_VERIFICATION  Fig 3: raw BEM A_kk(omega) and B_kk(omega) per mode.
%   2x3 grid: top row A, bottom row B.  Marks:
%     - A_inf horizontal line on each A panel (if present)
%     - hydrostatic T_n vertical for heave and pitch
%     - log y-axis for B (wide dynamic range; pitch B is ~100x smaller)
%
%   PHYSICS CHECKS this plot catches:
%     - B_kk -> 0 at low and high omega (correct physical limits)
%     - A_kk -> A_inf (finite constant) at high omega
%     - Pitch B_33 visibly << heave B_22 at the same frequency
%   If these shapes are wrong, every downstream calculation is corrupted.

    om_bem = D.bem.omega_BEM;
    A_diag = D.bem.A_diag;
    B_diag = D.bem.B_diag;
    A_inf  = nan(3, 1);
    if isfield(D.bem, 'AB_const_summary') && ...
       isfield(D.bem.AB_const_summary, 'A_inf_3DOF')
        A_inf = D.bem.AB_const_summary.A_inf_3DOF;
    end
    M_diag = D.bem.M_diag;
    Tn_hydro = D.bem.T_n_hydrostatic;   % [Inf; T_h; T_p]

    figure('Name', sprintf('[%s] Fig 3: BEM verification (A_{kk}, B_{kk})', D.station_id), ...
           'Color', 'w', 'Position', [120 60 1280 720]);

    A_units = {'kg', 'kg', 'kg.m^2'};
    B_units = {'N.s/m', 'N.s/m', 'N.m.s/rad'};

    for k = 1:3
        % --- Top row: A_kk(omega) ---
        subplot(2, 3, k);
        plot(om_bem, A_diag(k, :), '-', 'Color', colors{k}, ...
             'LineWidth', 1.8, 'DisplayName', sprintf('A_{%d%d}(\\omega)', k, k));
        hold on;
        if isfinite(A_inf(k))
            yline(A_inf(k), '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.2, ...
                  'Label', sprintf('A_\\infty = %.0f', A_inf(k)), ...
                  'HandleVisibility', 'off');
        end
        if k > 1 && isfinite(Tn_hydro(k)) && Tn_hydro(k) > 0
            omn = 2*pi / Tn_hydro(k);
            if omn > om_bem(1) && omn < om_bem(end)
                xline(omn, ':', 'Color', colors{k}, 'LineWidth', 1.2, ...
                      'Label', sprintf('T_{n,hydro}=%.2fs', Tn_hydro(k)), ...
                      'HandleVisibility', 'off');
            end
        end
        hold off;
        xlabel('\omega (rad/s)');
        ylabel(sprintf('A_{%d%d}  (%s)', k, k, A_units{k}));
        title(sprintf('%s: added mass A_{%d%d}(\\omega)  [M=%g %s]', ...
                      labels{k}, k, k, M_diag(k), A_units{k}));
        grid on;
        xlim([om_bem(1), om_bem(end)]);

        % --- Bottom row: B_kk(omega) ---
        subplot(2, 3, 3 + k);
        % Plot positive B with log y; show negative B as red markers
        Bk = B_diag(k, :);
        Bk_pos = Bk;  Bk_pos(Bk_pos <= 0) = NaN;
        semilogy(om_bem, Bk_pos, '-', 'Color', colors{k}, ...
                 'LineWidth', 1.8, 'DisplayName', sprintf('B_{%d%d}(\\omega)', k, k));
        hold on;
        if any(Bk < 0)
            neg_om = om_bem(Bk < 0);
            neg_B  = abs(Bk(Bk < 0));
            semilogy(neg_om, neg_B, 'rx', 'MarkerSize', 8, ...
                     'DisplayName', '|B|<0 (BEM noise)');
        end
        if k > 1 && isfinite(Tn_hydro(k)) && Tn_hydro(k) > 0
            omn = 2*pi / Tn_hydro(k);
            if omn > om_bem(1) && omn < om_bem(end)
                xline(omn, ':', 'Color', colors{k}, 'LineWidth', 1.2, ...
                      'HandleVisibility', 'off');
            end
        end
        hold off;
        xlabel('\omega (rad/s)');
        ylabel(sprintf('B_{%d%d}  (%s)', k, k, B_units{k}));
        title(sprintf('%s: radiation damping B_{%d%d}(\\omega)', labels{k}, k, k));
        grid on;
        xlim([om_bem(1), om_bem(end)]);
        if any(Bk < 0)
            legend('show', 'Location', 'best');
        end
    end

    sgtitle(sprintf(['BEM verification: A_{kk}(\\omega) and B_{kk}(\\omega) ' ...
                     'all three modes -- %s'], D.station_id));
end

% =====================================================================
function plot_scatter_with_targets(D, labels, colors)
%PLOT_SCATTER_WITH_TARGETS  Fig 4: (Hs, Te) probability heatmap with the
%   three target periods overlaid as vertical lines.
%
%   Catches:
%     - Targets falling at nonsensical Te values (period-domain mismatch,
%       moments computed in f-form used with omega-form formulas)
%     - Surge target at long periods where energetic events concentrate
%     - Heave target near the high-probability cluster
%     - Pitch target in the storm-period regime

    if ~isfield(D, 'climateGrid') || isempty(D.climateGrid)
        warning('run_MWEC_tuning:noScatter', ...
                'climateGrid not in figdata; skipping scatter diagram.');
        return;
    end
    cg = D.climateGrid;
    if ~isfield(cg, 'probability_grid') || ...
       ~isfield(cg, 'Hs_centers')       || ...
       ~isfield(cg, 'Te_centers')
        warning('run_MWEC_tuning:scatterFields', ...
                'climateGrid missing scatter fields; skipping Fig 4.');
        return;
    end
    P    = cg.probability_grid;
    Hs_c = cg.Hs_centers(:);
    Te_c = cg.Te_centers(:);

    figure('Name', sprintf('[%s] Fig 4: scatter (Hs, Te) + target periods', D.station_id), ...
           'Color', 'w', 'Position', [160 100 950 540]);

    % Display in pct so the colourbar is more readable
    Ppct = P * 100;
    Ppct(Ppct == 0) = NaN;   % white background for empty cells
    imagesc(Te_c, Hs_c, Ppct);
    set(gca, 'YDir', 'normal');
    colormap(parula);
    cb = colorbar;
    cb.Label.String = 'Probability  (% of records)';
    hold on;

    % Target Tn lines
    yl = ylim;
    for k = 1:3
        xline(D.T_n_opt(k), '--', 'Color', colors{k}, 'LineWidth', 2.0, ...
              'Label', sprintf('%s T_n=%.2fs', labels{k}, D.T_n_opt(k)), ...
              'LabelOrientation', 'horizontal', ...
              'LabelVerticalAlignment', 'top', ...
              'HandleVisibility', 'off');
    end
    % 90% energy band edges in T-domain
    T_L = 2*pi / D.omega_H;
    T_H = 2*pi / D.omega_L;
    xline(T_L, 'k:', 'LineWidth', 1.0, 'Label', 'T_L (90% band)', ...
          'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
    xline(T_H, 'k:', 'LineWidth', 1.0, 'Label', 'T_H (90% band)', ...
          'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
    if isfield(D, 'spectral_info') && D.spectral_info.bimodal
        T_partition = 2*pi / D.spectral_info.partition_omega;
        xline(T_partition, 'm-.', 'LineWidth', 1.2, ...
              'Label', sprintf('partition T=%.1fs', T_partition), ...
              'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
    end
    ylim(yl);
    hold off;

    xlabel('Energy period T_e  (s)');
    ylabel('Significant wave height H_s  (m)');
    title(sprintf('Climate scatter with placed T_n targets -- %s', D.station_id));
    grid on;

    % Annotation: sanity check on placement vs scatter
    [~, imax] = max(P(:));
    [iH, iT]  = ind2sub(size(P), imax);
    text(0.02, 0.98, sprintf(['Modal cell: H_s=%.2f m, T_e=%.2f s (p=%.2f%%)\n' ...
                              'Placed T_n: surge %.2f s | heave %.2f s | pitch %.2f s'], ...
                              Hs_c(iH), Te_c(iT), P(iH,iT)*100, ...
                              D.T_n_opt(1), D.T_n_opt(2), D.T_n_opt(3)), ...
         'Units', 'normalized', 'VerticalAlignment', 'top', ...
         'BackgroundColor', [1 1 1 0.8], 'FontSize', 9, 'Margin', 4);
end

% =====================================================================
function plot_falnes_panel(D, labels, colors)
%PLOT_FALNES_PANEL  Fig 5: two-panel figure.
%   Top:    S_ew(omega)/omega^3 with 90% band, placed omega_n markers.
%   Bottom: stacked eta_F integrand per mode:
%             g_k(omega) = nu_k * P_k(omega) * S_ew(omega)/omega^3 / denom
%             denom      = (sum nu) * int_{wL}^{wH} S_ew/omega^3 domega
%           Area under the stacked total in the placement band = eta_F.

    omega = D.omega;
    S_ew  = D.S_ew;
    mask_band = (omega >= D.omega_L) & (omega <= D.omega_H);

    nu = D.falnes_weights;
    if numel(nu) ~= 3
        warning('plot_falnes_panel:badNu', ...
                'falnes_weights malformed (numel=%d); falling back to [2;1;2].', numel(nu));
        nu = [2; 1; 2];
    end
    nu_total = sum(nu);

    SF = S_ew ./ max(omega, 1e-9).^3;
    denom = nu_total * trapz(omega(mask_band), SF(mask_band));
    if denom <= 0, denom = eps; end

    % Per-mode matching factor 4 beta/(1+beta)^2 (1 if not provided -- legacy)
    if isfield(D, 'eta_match') && numel(D.eta_match) == 3
        eta_match = D.eta_match(:);
        objective_label = '\eta_{power}';
    else
        eta_match = ones(3, 1);
        objective_label = '\eta_{F}';
    end

    % Per-mode normalised integrand (so sum_k area_in_band = eta_F or eta_power)
    integrand = zeros(numel(omega), 3);
    for k = 1:3
        omk = D.omega_n_opt(k);  Qk = D.Q_used(k);
        det = omega./omk - omk./omega;
        Pk  = 1 ./ (1 + Qk^2 .* det.^2);
        integrand(:, k) = nu(k) .* eta_match(k) .* Pk .* SF ./ denom;
    end
    integrand_total = sum(integrand, 2);

    % Verify eta_F / eta_power (sanity print -- should match D.eta_F)
    eta_band = trapz(omega(mask_band), integrand_total(mask_band));

    figure('Name', sprintf('[%s] Fig 5: Falnes weighting + %s integrand', D.station_id, objective_label), ...
           'Color', 'w', 'Position', [200 60 1100 760]);

    % --- Top panel: raw Falnes-weighted spectrum ---
    subplot(2, 1, 1);
    plot(omega, SF, 'b-', 'LineWidth', 1.6, 'DisplayName', 'S_{ew}(\omega)/\omega^3'); hold on;
    fill([omega(mask_band); flipud(omega(mask_band))], ...
         [SF(mask_band); zeros(sum(mask_band),1)], ...
         [0.85 0.92 1.0], 'EdgeColor','none', 'FaceAlpha',0.35, ...
         'DisplayName', '90% energy band');
    for k = 1:3
        xline(D.omega_n_opt(k), ':', 'Color', colors{k}, 'LineWidth', 1.5, ...
              'Label', sprintf('%s', labels{k}(1)), ...
              'HandleVisibility','off');
    end
    hold off;
    xlabel('\omega (rad/s)');
    ylabel('S_{ew}(\omega)/\omega^3  (m^2 s^4 / rad^4)');
    title('(top)  Falnes-weighted spectrum -- the function the optimiser integrates');
    legend('show', 'Location', 'northeast');
    grid on;
    add_period_axis(gca);

    % --- Bottom panel: stacked eta_F / eta_power integrand ---
    subplot(2, 1, 2);
    Y = integrand';   % [3 x N]
    h = area(omega, Y', 'LineStyle', 'none');
    h(1).FaceColor = colors{1}; h(1).FaceAlpha = 0.6;
    h(2).FaceColor = colors{2}; h(2).FaceAlpha = 0.6;
    h(3).FaceColor = colors{3}; h(3).FaceAlpha = 0.6;
    h(1).DisplayName = sprintf('%s (\\nu_1=%g)', labels{1}, nu(1));
    h(2).DisplayName = sprintf('%s (\\nu_2=%g)', labels{2}, nu(2));
    h(3).DisplayName = sprintf('%s (\\nu_3=%g)', labels{3}, nu(3));
    hold on;
    plot(omega, integrand_total, 'k-', 'LineWidth', 1.6, ...
         'DisplayName', 'Total integrand');
    % Mark the 90% band edges and integration area
    yl = ylim;
    fill([D.omega_L D.omega_L D.omega_H D.omega_H], ...
         [0 yl(2) yl(2) 0], [0.95 0.95 0.95], ...
         'EdgeColor', 'none', 'FaceAlpha', 0.10, 'HandleVisibility','off');
    for k = 1:3
        xline(D.omega_n_opt(k), ':', 'Color', colors{k}, 'LineWidth', 1.5, ...
              'HandleVisibility','off');
    end
    hold off;
    xlabel('\omega (rad/s)');
    ylabel(sprintf('g_k(\\omega) -- normalised integrand of %s', objective_label));
    title(sprintf(['(bottom)  Stacked %s integrand -- area in band ' ...
                   '= %s = %.4f (check)'], objective_label, objective_label, eta_band));
    legend('show', 'Location', 'northeast');
    grid on;
    add_period_axis(gca);

    % No sgtitle (clashes with the period-axis tick labels from
    % add_period_axis on the top subplot).  Station_id is already in the
    % figure Name and in the subplot titles.
end

% =====================================================================
function plot_absorption_vs_T(D, labels, colors)
%PLOT_ABSORPTION_VS_T  Figure 7: per-mode RAO-based absorption (%) at
%   each wave period T, using the full BEM Fe(omega) and the Stage-1
%   closure damping.  Three normalisations:
%     (A) share of Falnes-bounded total ceiling at that frequency
%     (B) capture width / L_ref (raft y-extent)
%     (C) capture width / wavelength

    A = D.absorption_vs_T;
    figure('Name', sprintf('[%s] Fig 6: full-RAO absorption vs T', D.station_id), ...
           'Color', 'w', 'Position', [320 80 1280 760]);

    panels = { ...
        struct('tag','A', 'title',sprintf('(A) Share of Falnes max ceiling (denom = sum nu * g/omega^2; sum nu = %g)', sum(A.nu)), ...
               'ylab','% of Falnes ceiling', ...
               'sk',A.pct_A_surge,'hk',A.pct_A_heave,'pk',A.pct_A_pitch,'tot',A.pct_A_total), ...
        struct('tag','B', 'title',sprintf('(B) Capture width / L_{ref}   (L_{ref} = %.3f m, raft y-extent)', D.L_ref), ...
               'ylab','% of L_{ref}-width wave power', ...
               'sk',A.pct_B_surge,'hk',A.pct_B_heave,'pk',A.pct_B_pitch,'tot',A.pct_B_total), ...
        struct('tag','C', 'title','(C) Capture width / wavelength', ...
               'ylab','% of wavelength (= l_k / \lambda)', ...
               'sk',A.pct_C_surge,'hk',A.pct_C_heave,'pk',A.pct_C_pitch,'tot',A.pct_C_total) };

    for ip = 1:numel(panels)
        P = panels{ip};
        subplot(3, 1, ip);

        % Stacked-area: surge bottom, heave middle, pitch top
        Y = [P.sk(:)'; P.hk(:)'; P.pk(:)'];
        h = area(A.T, Y', 'LineStyle', 'none');
        h(1).FaceColor = colors{1}; h(1).FaceAlpha = 0.55;
        h(2).FaceColor = colors{2}; h(2).FaceAlpha = 0.55;
        h(3).FaceColor = colors{3}; h(3).FaceAlpha = 0.55;
        h(1).DisplayName = labels{1};
        h(2).DisplayName = labels{2};
        h(3).DisplayName = labels{3};
        hold on;

        % Total line on top
        plot(A.T, P.tot, 'k-', 'LineWidth', 1.6, 'DisplayName', 'Total');

        % Mark the three placed T_n
        ymax_obs = max([P.tot(:); 100]) * 1.05;
        for k = 1:3
            xline(D.T_n_opt(k), ':', 'Color', colors{k}, 'LineWidth', 1.4, ...
                  'Label', sprintf('T_{n,%s}=%.2fs', labels{k}(1), D.T_n_opt(k)), ...
                  'LabelVerticalAlignment','bottom', ...
                  'LabelHorizontalAlignment','center', ...
                  'HandleVisibility','off');
        end

        % Theory ceilings for panel A (dashed horizontal lines)
        if strcmp(P.tag, 'A')
            nu_tot = sum(A.nu);
            yline(A.nu(1)/nu_tot*100, ':', 'Color', colors{1}, ...
                  'Label', sprintf('surge max %.1f%%', A.nu(1)/nu_tot*100), ...
                  'HandleVisibility','off');
            yline(A.nu(2)/nu_tot*100, ':', 'Color', colors{2}, ...
                  'Label', sprintf('heave max %.1f%%', A.nu(2)/nu_tot*100), ...
                  'HandleVisibility','off');
            yline(100, 'k:', 'Label','axisymmetric ideal 100%', ...
                  'HandleVisibility','off');
        end

        hold off;
        xlabel('Wave period T (s)');
        ylabel(P.ylab);
        title(P.title);
        xlim([A.T(1), A.T(end)]);
        if strcmp(P.tag, 'A')
            ylim([0, max(105, ymax_obs)]);
        else
            ylim([0, ymax_obs]);
        end
        grid on;
        if ip == 1
            legend('show', 'Location', 'northeast');
        end
    end

    sgtitle(sprintf('Per-period absorption (full-RAO, Stage-1 closure) -- %s', D.station_id));

    % --- Print a small console table at round T values ---
    T_print = [3, 4, 5, 6, 8, 10];
    fprintf('\n   Per-T RAO-based absorption (panel B, %% of L_{ref}=%.2f m):\n', D.L_ref);
    fprintf('     T(s)   surge%%   heave%%   pitch%%   total%%\n');
    for tt = T_print
        [~, ii] = min(abs(A.T - tt));
        fprintf('     %4.1f   %6.2f   %6.2f   %6.2f   %6.2f\n', ...
                A.T(ii), A.pct_B_surge(ii), A.pct_B_heave(ii), ...
                A.pct_B_pitch(ii), A.pct_B_total(ii));
    end
end

% =====================================================================
function ax_top = add_period_axis(ax_bot)
%ADD_PERIOD_AXIS  Add a top x-axis showing period T = 2*pi/omega.
%   Call AFTER all plotting on ax_bot is finished.  Returns the handle
%   of the transparent twin axes.  Standard period ticks (3--30 s) are
%   placed where they fall inside the current x-limits of ax_bot.  The
%   title of ax_bot, if any, is moved to ax_top so it sits above the
%   period-axis tick labels (avoiding the label--title collision).

    if nargin < 1 || isempty(ax_bot), ax_bot = gca; end
    xl = xlim(ax_bot);
    if xl(1) <= 0, xl(1) = max(xl(2)*1e-3, 1e-3); end

    % Shrink ax_bot height slightly so labels + title fit above it
    pos = ax_bot.Position;
    new_h = pos(4) * 0.90;
    ax_bot.Position = [pos(1), pos(2), pos(3), new_h];

    ax_top = axes('Position',       ax_bot.Position, ...
                  'XAxisLocation',  'top', ...
                  'YAxisLocation',  'right', ...
                  'Color',          'none', ...
                  'YTick',          [], ...
                  'YColor',         'none', ...
                  'Box',            'off', ...
                  'XLim',           xl, ...
                  'HitTest',        'off', ...
                  'PickableParts',  'none');

    % Period ticks restricted to operationally useful waves (3--30 s).
    % Beyond 2.1 rad/s (~T<3 s) the ticks crowd together and add no info.
    T_candidates  = [30, 20, 15, 12, 10, 8, 6, 5, 4, 3];
    om_candidates = 2*pi ./ T_candidates;
    keep = (om_candidates >= xl(1)) & (om_candidates <= xl(2));
    om_ticks = om_candidates(keep);
    T_ticks  = T_candidates(keep);
    % Limit to ~7 ticks for readability
    if numel(om_ticks) > 7
        idx = round(linspace(1, numel(om_ticks), 7));
        om_ticks = om_ticks(idx);
        T_ticks  = T_ticks(idx);
    end
    if ~isempty(om_ticks)
        set(ax_top, 'XTick', om_ticks, ...
                    'XTickLabel', arrayfun(@(t) sprintf('%g', t), T_ticks, ...
                                            'UniformOutput', false));
    end
    xlabel(ax_top, 'Period T (s)');

    % Move title from ax_bot to ax_top so it lands above the period axis
    % rather than overlapping the period tick labels.
    if ~isempty(ax_bot.Title.String)
        title_text = ax_bot.Title.String;
        title_intp = ax_bot.Title.Interpreter;
        title(ax_bot, '');
        title(ax_top, title_text, 'Interpreter', title_intp);
    end
end

% =====================================================================
function [omega_c, partition_omega, info] = compute_band_partition( ...
        S_ew, omega, omega_L, omega_H)
%COMPUTE_BAND_PARTITION  Energy-weighted centroid + optional bimodal partition.
%
%   Unimodal: partition_omega = omega_c (spectral centroid in [omega_L,omega_H]).
%   Bimodal:  partition_omega = valley frequency between the two dominant peaks.
%             Bimodality declared when valley / higher-peak < 0.60.
%
%   No Signal Processing Toolbox required: manual Gaussian smoothing +
%   local-maxima search.
%
%   OUTPUTS
%     omega_c         scalar   energy-weighted centroid of S_ew on [omega_L,omega_H]
%     partition_omega scalar   partition frequency (rad/s)
%     info            struct   .bimodal, .valley_ratio, .valley_omega,
%                              .peak_omegas, .partition_omega, .omega_c

    omega = omega(:);  S_ew = S_ew(:);
    mask  = (omega >= omega_L) & (omega <= omega_H);
    om_b  = omega(mask);
    Sb    = S_ew(mask);

    % Energy-weighted centroid
    m0_b = trapz(om_b, Sb);
    if m0_b <= 0
        omega_c = (omega_L + omega_H) / 2;
    else
        omega_c = trapz(om_b, Sb .* om_b) / m0_b;
    end

    % Manual Gaussian smoothing (sigma ~ 5% of N_b, min 2 points)
    N_b   = numel(om_b);
    sigma = max(0.05 * N_b, 2);
    hw    = min(ceil(2.5 * sigma), floor((N_b - 1) / 2));
    kk    = (-hw:hw)';
    kern  = exp(-0.5 * (kk / sigma).^2);
    kern  = kern / sum(kern);
    Sb_sm = conv(Sb, kern, 'same');

    % Local maxima: strict inequality, height > 10% of smoothed peak
    min_height = 0.10 * max(Sb_sm);
    is_max = false(N_b, 1);
    for i = 2:N_b-1
        if Sb_sm(i) > Sb_sm(i-1) && Sb_sm(i) > Sb_sm(i+1) && Sb_sm(i) >= min_height
            is_max(i) = true;
        end
    end
    pk_idx  = find(is_max);
    pk_vals = Sb_sm(pk_idx);

    bimodal      = false;
    valley_ratio = NaN;
    valley_omega = NaN;
    peak_omegas  = [];

    if numel(pk_idx) >= 2
        [~, order] = sort(pk_vals, 'descend');
        idx_a = pk_idx(order(1));
        idx_b = pk_idx(order(2));
        if idx_a > idx_b, [idx_a, idx_b] = deal(idx_b, idx_a); end

        [val_val, val_rel] = min(Sb_sm(idx_a:idx_b));
        val_abs = idx_a + val_rel - 1;
        ratio   = val_val / max(pk_vals(order(1)), pk_vals(order(2)));

        if ratio < 0.60
            bimodal      = true;
            valley_ratio = ratio;
            valley_omega = om_b(val_abs);
            peak_omegas  = [om_b(idx_a); om_b(idx_b)];
        end
    end

    if bimodal
        partition_omega = valley_omega;
    else
        partition_omega = omega_c;
    end

    info = struct( ...
        'bimodal',          bimodal, ...
        'valley_ratio',     valley_ratio, ...
        'valley_omega',     valley_omega, ...
        'partition_omega',  partition_omega, ...
        'peak_omegas',      peak_omegas, ...
        'omega_c',          omega_c);
end

% =====================================================================
function [Q_corr, closure_corr] = correct_Q_at_placed(omega_n_placed, bem, closure)
%CORRECT_Q_AT_PLACED  Re-evaluate Evans-optimal Q at the placed natural frequencies.
%
%   The initial Q scalars are computed at the hydrostatic BEM frequencies or
%   band centroid.  This function re-evaluates them at the actual placed omega_n
%   using the full frequency-dependent BEM A_kk(omega) and B_kk(omega).
%
%   Evans-optimal B_PTO is also updated to B_22(omega_n,heave_placed).

    omega_n_placed = omega_n_placed(:);
    M     = closure.M;
    Iyy   = closure.Iyy;
    k_gyr = closure.k_gyr;

    % Evans-optimal B_PTO at placed heave frequency
    A22_p = max(bem.A_kk_func(2, omega_n_placed(2)), 0);
    B22_p = max(bem.B_kk_func(2, omega_n_placed(2)), 0);
    if B22_p <= 0
        warning('correct_Q_at_placed:B22zero', ...
                'B_22 at placed heave omega_n=%.4f is %.4g; reverting to initial B_PTO.', ...
                omega_n_placed(2), B22_p);
        B22_p = closure.B_PTO_scalar;
    end
    B_PTO_s = B22_p;
    B_PTO_r = B_PTO_s * k_gyr^2;

    % Surge at placed omega_n1
    A11_p = max(bem.A_kk_func(1, omega_n_placed(1)), 0);
    B11_p = max(bem.B_kk_func(1, omega_n_placed(1)), 0);
    Q_surge = omega_n_placed(1) * (M + A11_p) / (B11_p + B_PTO_s);

    % Heave at placed omega_n2  (Evans optimal: Q = Q_rad,heave / 2)
    Q_heave = omega_n_placed(2) * (M + A22_p) / (B22_p + B_PTO_s);

    % Pitch at placed omega_n3
    A33_p = max(bem.A_kk_func(3, omega_n_placed(3)), 0);
    B33_p = max(bem.B_kk_func(3, omega_n_placed(3)), 0);
    Q_pitch = omega_n_placed(3) * (Iyy + A33_p) / (B33_p + B_PTO_r);

    Q_corr = [Q_surge; Q_heave; Q_pitch];

    closure_corr = closure;
    closure_corr.B_PTO_scalar  = B_PTO_s;
    closure_corr.B_PTO_rot     = B_PTO_r;
    closure_corr.A22_at_placed = A22_p;
    closure_corr.B22_at_placed = B22_p;
    closure_corr.A11_at_placed = A11_p;
    closure_corr.B11_at_placed = B11_p;
    closure_corr.A33_at_placed = A33_p;
    closure_corr.B33_at_placed = B33_p;
    closure_corr.Q_surge       = Q_surge;
    closure_corr.Q_heave       = Q_heave;
    closure_corr.Q_pitch       = Q_pitch;
end

% closed_loop_heave_Q (heave-SDOF + tether PTO sensitivity) deleted in
% Stage-1 refactor: PTO geometry (K_PTO, alpha) is a Stage-2 concern and
% no longer enters the Stage-1 deliverable.

% =====================================================================
function [closure, omega_n_opt, eta_opt, opt_info] = ...
        joint_passive_placement(bem, final_props_path, pdata, opts)
%JOINT_PASSIVE_PLACEMENT  5-variable passive-feasible eta_power optimisation.
%
%   Variables x = [omega_n,1; omega_n,2; omega_n,3; alpha; B_PTO].
%   Geometry (alpha, ell) and stiffness (K_PTO, K_moor) are derived from the
%   placement via the symmetric-tether matrix-structure stiffness equations.
%   Passive constraints: K_PTO >= 0 (S_2 >= 0), K_moor >= 0 (S_1 >= S_2 cot^2 a),
%   pitch S_3 >= 0, hull ell <= b_t,max sin a + h_t,max cos a.

    FP = load(final_props_path);
    fp = FP.final_props;
    M   = fp.mass_total;
    Iyy = fp.Iyy;
    K22 = fp.K_hydro(2,2);
    K33 = fp.K_hydro(3,3);

    % Hull bounds
    cs = fp.cross_section;
    bt_max = max(abs(cs(:,1)));
    ht_max = max(cs(:,2));

    % Hydrostatic Tn (passive lower bounds for heave/pitch omega_n)
    om_hydro_h = 2*pi / bem.T_n_hydrostatic(2);
    om_hydro_p = 2*pi / bem.T_n_hydrostatic(3);
    k_gyr      = sqrt(max(Iyy / max(M, eps), 0));

    % Bundle for objective / constraint
    omega   = pdata.omega(:);
    S_ew    = pdata.S_ew(:);
    omega_L = pdata.omega_L;
    omega_H = pdata.omega_H;
    mask    = (omega >= omega_L) & (omega <= omega_H);
    SF      = S_ew ./ max(omega, 1e-9).^3;
    nu      = opts.falnes_weights(:);
    denom   = sum(nu) * trapz(omega(mask), SF(mask));
    if denom <= 0, denom = eps; end

    d.bem = bem; d.M = M; d.Iyy = Iyy; d.M_diag = [M; M; Iyy];
    d.K22_hyd = K22; d.K33_hyd = K33;
    d.bt_max = bt_max; d.ht_max = ht_max;
    d.omega = omega; d.SF = SF; d.mask = mask; d.denom = denom; d.nu = nu;

    % Bounds.  Heave/pitch lb = hydrostatic (passive cap on T_n).  Surge in 90% band.
    % If opts.use_partition_zones is true, surge/pitch zones from omega_bounds apply.
    if opts.use_partition_zones && isfield(pdata, 'omega_bounds') && size(pdata.omega_bounds,1) == 3
        lb1 = pdata.omega_bounds(1,1);   ub1 = pdata.omega_bounds(1,2);
        lb2 = max(pdata.omega_bounds(2,1), om_hydro_h);  ub2 = pdata.omega_bounds(2,2);
        lb3 = max(pdata.omega_bounds(3,1), om_hydro_p);  ub3 = pdata.omega_bounds(3,2);
    else
        lb1 = omega_L;   ub1 = omega_H;
        lb2 = om_hydro_h; ub2 = omega_H;
        lb3 = om_hydro_p; ub3 = omega_H;
    end
    lb = [lb1; lb2; lb3; deg2rad(opts.alpha_bounds(1)); opts.B_PTO_bounds(1)];
    ub = [ub1; ub2; ub3; deg2rad(opts.alpha_bounds(2)); opts.B_PTO_bounds(2)];

    % Multistart seeds (8 hand-picked + remaining random)
    rng(42);
    nstarts = max(opts.nstarts, 1);
    seeds = zeros(nstarts, 5);
    seeds(1, :) = [0.5,   max(om_hydro_h*1.001, lb2),  max(om_hydro_p*1.001, lb3),  deg2rad(45),   500];
    seeds(2, :) = [0.7,   max(om_hydro_h*1.001, lb2),  max(om_hydro_p*1.001, lb3),  deg2rad(45),  2000];
    seeds(3, :) = [1.0,   max(om_hydro_h*1.001, lb2),  max(om_hydro_p*1.001, lb3),  deg2rad(45),  5000];
    seeds(4, :) = [0.5,   1.5,                          max(om_hydro_p*1.001, lb3),  deg2rad(45),  3000];
    seeds(5, :) = [0.5,   max(om_hydro_h*1.001, lb2),   1.7,                         deg2rad(45),  1000];
    seeds(6, :) = [0.5,   max(om_hydro_h*1.001, lb2),  max(om_hydro_p*1.001, lb3),  deg2rad(35),   500];
    seeds(7, :) = [0.5,   max(om_hydro_h*1.001, lb2),  max(om_hydro_p*1.001, lb3),  deg2rad(55),   500];
    seeds(8, :) = [0.7,   1.5,                          1.7,                         deg2rad(50),  2000];
    for s = 9:nstarts
        seeds(s,1) = lb(1) + (ub(1)-lb(1))*rand;
        seeds(s,2) = lb(2) + (ub(2)-lb(2))*rand;
        seeds(s,3) = lb(3) + (ub(3)-lb(3))*rand;
        seeds(s,4) = lb(4) + (ub(4)-lb(4))*rand;
        seeds(s,5) = exp(log(max(lb(5),1)) + (log(ub(5))-log(max(lb(5),1)))*rand);
    end
    % Clip seeds to bounds (defensive)
    for s = 1:nstarts
        seeds(s,:) = max(min(seeds(s,:), ub'), lb');
    end

    % Optimise
    obj_fn = @(x) -compute_eta_power_local(x, d);
    con_fn = @(x) feasibility_local(x, d);
    optsfm = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
        'TolFun', 1e-10, 'TolX', 1e-8, 'MaxIterations', 1500, ...
        'MaxFunctionEvaluations', 8000, 'ConstraintTolerance', 1e-8);

    bestE = -inf; bestX = nan(5,1);
    multistart = nan(nstarts, 7);  % [eta, x(1..5), exitflag]
    for s = 1:nstarts
        try
            [x_s, neg_e_s, ex_s] = fmincon(obj_fn, seeds(s,:)', [], [], [], [], ...
                lb, ub, con_fn, optsfm);
            e_s = -neg_e_s;
            multistart(s,:) = [e_s, x_s(1), x_s(2), x_s(3), x_s(4), x_s(5), ex_s];
            if e_s > bestE
                bestE = e_s; bestX = x_s;
            end
        catch ME
            if opts.verbose
                fprintf('   [joint_passive_placement] start %d failed: %s\n', s, ME.message);
            end
        end
    end

    if ~isfinite(bestE)
        error('joint_passive_placement:noFeasible', ...
              'All multistart restarts failed to find a passive-feasible point.');
    end

    % Diagnostics at the optimum
    [eta_opt, deriv] = compute_eta_power_local(bestX, d);
    omega_n_opt = bestX(1:3);
    alpha_opt   = bestX(4);
    B_PTO_opt   = bestX(5);

    % Per-mode contribution to eta_power
    contrib = zeros(3,1);
    for k = 1:3
        det = omega./omega_n_opt(k) - omega_n_opt(k)./omega;
        Pk  = 1 ./ (1 + deriv.Q_k(k)^2 .* det.^2);
        contrib(k) = nu(k) * deriv.eta_match(k) * trapz(omega(mask), Pk(mask) .* SF(mask)) / denom;
    end

    % Pack closure struct.  Maintain legacy fields where meaningful.
    closure = struct();
    closure.M          = M;
    closure.Iyy        = Iyy;
    closure.K22        = K22;
    closure.K33        = K33;
    closure.k_gyr      = k_gyr;
    closure.alpha      = alpha_opt;
    closure.ell        = deriv.ell;
    closure.B_PTO_scalar = B_PTO_opt;
    closure.K_PTO_scalar = deriv.K_PTO;
    closure.K_moor       = deriv.K_moor;
    closure.B_PTO_kk     = deriv.B_PTO_kk;   % [N.s/m; N.s/m; N.m.s/rad]
    closure.K_PTO_kk     = deriv.K_PTO_kk;
    closure.B_PTO_rot    = deriv.B_PTO_kk(3);
    closure.K_total_diag = [deriv.K_PTO_kk(1) + deriv.K_moor;
                            K22 + deriv.K_PTO_kk(2);
                            K33 + deriv.K_PTO_kk(3)];
    closure.B13          = deriv.B13;
    closure.K13          = deriv.K13;
    closure.beta         = deriv.beta;
    closure.eta_match    = deriv.eta_match;
    closure.Q_k          = deriv.Q_k;
    closure.A_at_placed  = deriv.A_at;
    closure.B_rad_at_placed = deriv.B_rad;
    closure.S1           = deriv.S1;
    closure.S2           = deriv.S2;
    closure.S3           = deriv.S3;
    closure.omega_n_h_hydro = om_hydro_h;
    closure.omega_n_p_hydro = om_hydro_p;
    closure.omega_c         = pdata.omega_c;
    closure.passive_feasible = true;
    closure.objective    = 'eta_power_passive';
    % Legacy aliases consumed by some plot/print helpers:
    closure.A22_at_omh = bem.A_kk_func(2, om_hydro_h);
    closure.B22_at_omh = max(bem.B_kk_func(2, om_hydro_h), 0);
    closure.A33_at_omp = bem.A_kk_func(3, om_hydro_p);
    closure.B33_at_omp = max(bem.B_kk_func(3, om_hydro_p), 0);
    closure.A11_at_omc = bem.A_kk_func(1, pdata.omega_c);
    closure.B11_at_omc = max(bem.B_kk_func(1, pdata.omega_c), 0);
    closure.Q_surge = deriv.Q_k(1);
    closure.Q_heave = deriv.Q_k(2);
    closure.Q_pitch = deriv.Q_k(3);

    % opt_info: mirror the legacy run_lorentzian_placement output for compat
    opt_info = struct();
    opt_info.eta_per_mode = contrib;
    opt_info.multistart   = multistart;
    opt_info.objective    = 'eta_power_passive';
    opt_info.note         = 'Joint (omega_n, alpha, B_PTO) optimisation under passive K_PTO matrix structure.';
end

% =====================================================================
function [eta, deriv] = compute_eta_power_local(x, d)
%COMPUTE_ETA_POWER_LOCAL  Inner objective for joint_passive_placement.
    om_n  = x(1:3); a = x(4); B_PTO = x(5);
    A_at = max([d.bem.A_kk_func(1, om_n(1));
                d.bem.A_kk_func(2, om_n(2));
                d.bem.A_kk_func(3, om_n(3))], 0);
    S1 = om_n(1)^2 * (d.M   + A_at(1));
    S2 = om_n(2)^2 * (d.M   + A_at(2)) - d.K22_hyd;
    S3 = om_n(3)^2 * (d.Iyy + A_at(3)) - d.K33_hyd;
    if S2 <= 0 || S3 < 0
        eta = -1e3;
        deriv = struct('K_PTO',NaN,'ell',NaN,'K_moor',NaN, ...
            'B_PTO_kk',nan(3,1),'K_PTO_kk',nan(3,1),'B13',NaN,'K13',NaN, ...
            'beta',nan(3,1),'eta_match',nan(3,1),'Q_k',nan(3,1), ...
            'S1',S1,'S2',S2,'S3',S3,'A_at',A_at,'B_rad',nan(3,1));
        return;
    end
    K_PTO  = S2 / (2*sin(a)^2);
    ell    = sin(a) * sqrt(max(S3,0)/S2);
    K_moor = S1 - 2*K_PTO*cos(a)^2;
    g       = [2*cos(a)^2; 2*sin(a)^2; 2*ell^2];
    B_PTO_kk = B_PTO * g;
    K_PTO_kk = K_PTO * g;
    B_rad  = max([d.bem.B_kk_func(1, om_n(1));
                  d.bem.B_kk_func(2, om_n(2));
                  d.bem.B_kk_func(3, om_n(3))], 0);
    beta      = B_PTO_kk ./ max(B_rad, eps);
    eta_match = 4*beta ./ (1+beta).^2;
    Q_k       = om_n .* (d.M_diag + A_at) ./ max(B_rad + B_PTO_kk, eps);
    eta = 0;
    for k = 1:3
        det = d.omega./om_n(k) - om_n(k)./d.omega;
        Pk  = 1 ./ (1 + Q_k(k)^2 .* det.^2);
        eta = eta + d.nu(k) * eta_match(k) * trapz(d.omega(d.mask), Pk(d.mask) .* d.SF(d.mask));
    end
    eta = eta / d.denom;
    deriv = struct('K_PTO',K_PTO,'ell',ell,'K_moor',K_moor, ...
        'B_PTO_kk',B_PTO_kk,'K_PTO_kk',K_PTO_kk, ...
        'B13',-2*B_PTO*ell*cos(a),'K13',-2*K_PTO*ell*cos(a), ...
        'beta',beta,'eta_match',eta_match,'Q_k',Q_k, ...
        'S1',S1,'S2',S2,'S3',S3,'A_at',A_at,'B_rad',B_rad);
end

% =====================================================================
function [c, ceq] = feasibility_local(x, d)
%FEASIBILITY_LOCAL  Inner passive-feasibility constraint for joint solve.
    om_n = x(1:3); a = x(4);
    A_at = max([d.bem.A_kk_func(1, om_n(1));
                d.bem.A_kk_func(2, om_n(2));
                d.bem.A_kk_func(3, om_n(3))], 0);
    S1 = om_n(1)^2 * (d.M   + A_at(1));
    S2 = om_n(2)^2 * (d.M   + A_at(2)) - d.K22_hyd;
    S3 = om_n(3)^2 * (d.Iyy + A_at(3)) - d.K33_hyd;
    c1 = -S2;                                                       % heave passive
    c2 = -S3;                                                       % pitch passive
    c3 = S2*cos(a)^2 / max(sin(a)^2, eps) - S1;                     % K_moor >= 0
    if S2 > 0 && S3 >= 0
        ell      = sin(a) * sqrt(max(S3,0)/S2);
        ell_max  = d.bt_max*sin(a) + d.ht_max*cos(a);
        c4 = ell - ell_max;
    else
        c4 = -1;   % infeasibility already flagged by c1/c2
    end
    c   = [c1; c2; c3; c4];
    ceq = [];
end
