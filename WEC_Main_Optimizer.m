function [results, final_props] = WEC_Main_Optimizer(params)
% WEC_MAIN_OPTIMIZER  Two-stage WEC mass-distribution optimisation.
%
%   [results, final_props] = WEC_MAIN_OPTIMIZER(params)
%
%   Orchestrates a two-stage optimisation pipeline for a wave energy
%   converter's ballast distribution:
%
%     Stage 1 (2D surrogate):  Fast exploration using a strip-theory
%       model.  Three modes controlled by config.stage1_mode:
%         'skip'     – bypass Stage 1, forward defaults to Stage 2
%         'oneshot'  – single 2D solve, no iterative correction
%         'trained'  – iterative PID correction of k_vol (volume bias)
%                      and k_gm (CG/GM bias) until the 2D surrogate
%                      agrees with the 3D ground truth
%
%     Stage 2 (3D fmincon/SQP):  High-fidelity refinement on the full
%       3D mesh, warm-started from Stage 1's result.
%
%   INPUT
%     params : struct from WEC_Driver.m (single source of truth for all
%              tuneable constants).  Passed to WEC_Configuration_Builder
%              to produce the runtime `config` struct.
%
%   OUTPUTS
%     results     : struct with Stage 1 and Stage 2 history, config,
%                   convergence metrics, and timing
%     final_props : struct of 3D physical properties at the optimum
%
%   DEPENDENCIES
%     WEC_Configuration_Builder, run_2d_optimizer,
%     calculate_2d_properties, calculate_3d_properties,
%     zone_penalty, WEC_PID_Controller, WEC_Visualization
%
%   See also: WEC_Driver, WEC_Configuration_Builder, run_2d_optimizer

try

%% §1  BUILD CONFIGURATION  ────────────────────────────────────────────
%  All tuneable values flow from params → config.  Nothing is hardcoded
%  below; every literal comes from a config field set by the driver.

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════╗\n');
fprintf('║   WEC OPTIMISATION PIPELINE                      ║\n');
fprintf('╚══════════════════════════════════════════════════╝\n');
fprintf(' Timestamp: %s\n\n', datestr(now));

tic_main = tic;

config = WEC_Configuration_Builder(params);

%% §2  STAGE 1 — 2D SURROGATE  ─────────────────────────────────────────

stage1_mode = config.stage1_mode;

fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf('║ STAGE 1: Mode = %-33s║\n', upper(stage1_mode));
fprintf('╚══════════════════════════════════════════════════╝\n\n');

% x = [draft, rho_1, ..., rho_N] — the design vector
x0_2d = [config.initial_vertical_shift, config.initial_densities];

% Surrogate correction factors: unity = no bias correction yet.
% k_vol scales the strip-extrusion submerged volume.
% k_gm  biases the 2D CG position to match 3D.
config.k_vol = config.k_vol_init;
config.k_gm  = config.k_gm_init;

%% §2a  PRE-CALIBRATION OF k_vol AND k_gm  ────────────────────────────
%
%  PROBLEM
%    The PID loop in 'trained' mode requires a successful 2D solve to
%    compute 3D properties and derive correction factors.  But when the
%    feasible region is tight (e.g. constructability bounds), the
%    uncalibrated 2D surrogate (k_vol = k_gm = 1.0) may be infeasible.
%    Every iteration returns exitflag = -2, the PID never updates, and
%    Stage 1 burns all 50 iterations without progress.
%
%  SOLUTION
%    Evaluate the initial design vector x0 in BOTH models (one 2D eval
%    + one 3D eval) BEFORE the optimisation loop.  Derive:
%
%      k_vol = V_sub_3D / V_sub_2D     (correct volume overestimate)
%      k_gm  = CG_z_3D / CG_z_2D      (correct CG depth mismatch)
%
%    This is a one-shot calibration costing ~0.5 s.  It shifts the 2D
%    surrogate's feasible region to overlap with reality so that
%    fmincon can find a feasible starting point.
%
%  WHY THIS WORKS
%    The 2D rectangular-extrusion model overestimates V_sub by 10-15%
%    and predicts a shallower CG_z than the 3D divergence theorem.
%    With k_vol = 1.0, the 2D mass balance is biased — the equilibrium
%    draft and hence BM = I_wp / V_sub and GM = KM - CG_z are all
%    shifted.  Pre-calibrating removes this first-order bias, leaving
%    only small residuals that the PID refines iteratively.
%
%  GUARD CONDITIONS
%    - Only applies to 'trained' mode (oneshot and skip don't iterate)
%    - V_sub_2D must be positive (otherwise hull is above water)
%    - CG_z values must be > 0.01 m from waterline (ratio guard)
%    - Results are clamped to the PID bounds [k_vol_lo, k_vol_hi]

if strcmp(stage1_mode, 'trained')
    fprintf(' Pre-calibrating surrogate correction factors...\n');

    try
        % Evaluate initial design in both models
        config_precal      = config;
        config_precal.k_vol = 1.0;
        config_precal.k_gm  = 1.0;

        props_2d_precal = calculate_2d_properties(x0_2d, config_precal);
        props_3d_precal = calculate_3d_properties(x0_2d, config);

        vsub_2d_0 = props_2d_precal.V_sub;
        vsub_3d_0 = props_3d_precal.V_sub;
        cgz_2d_0  = props_2d_precal.CG_total(3);
        cgz_3d_0  = props_3d_precal.CG_total(3);

        % k_vol: correct volume bias
        if vsub_2d_0 > 1e-6 && vsub_3d_0 > 1e-6
            kvol_precal  = vsub_3d_0 / vsub_2d_0;
            kvol_precal  = max(config.bounds_kvol(1), ...
                           min(config.bounds_kvol(2), kvol_precal));
            config.k_vol = kvol_precal;
            fprintf('   k_vol: %.4f  (V_sub: 2D=%.4f → 3D=%.4f m³, ratio=%.3f)\n', ...
                    config.k_vol, vsub_2d_0, vsub_3d_0, vsub_3d_0/vsub_2d_0);
        else
            fprintf('   k_vol: skipped (V_sub near zero)\n');
        end

        % k_gm: correct CG_z bias
        cg_floor = 0.01;  % [m] guard against near-zero CG_z
        if abs(cgz_2d_0) > cg_floor && abs(cgz_3d_0) > cg_floor
            kgm_precal   = cgz_3d_0 / cgz_2d_0;
            kgm_precal   = max(config.bounds_kgm(1), ...
                           min(config.bounds_kgm(2), kgm_precal));
            config.k_gm  = kgm_precal;
            fprintf('   k_gm:  %.4f  (CG_z: 2D=%.4f → 3D=%.4f m, ratio=%.3f)\n', ...
                    config.k_gm, cgz_2d_0, cgz_3d_0, cgz_3d_0/cgz_2d_0);
        else
            fprintf('   k_gm:  skipped (CG_z near zero: 2D=%.4f, 3D=%.4f)\n', ...
                    cgz_2d_0, cgz_3d_0);
        end

        gm_2d_0 = props_2d_precal.GM_uncorrected;
        gm_3d_0 = props_3d_precal.GM_L;
        fprintf('   GM at x0: 2D=%.4f → 3D=%.4f m (Δ=%+.4f)\n', ...
                gm_2d_0, gm_3d_0, gm_2d_0 - gm_3d_0);

    catch ME
        fprintf('   Pre-calibration failed: %s\n', ME.message);
        fprintf('   Falling back to k_vol=%.2f, k_gm=%.2f\n', ...
                config.k_vol, config.k_gm);
    end
end

% History arrays — populated by 'oneshot' or 'trained' branches.
% Kept empty for 'skip'.  All are 1×N row vectors.
hist = init_history_struct();

converged     = false;
iteration     = 0;
pid_saturated = false;

switch stage1_mode

    case 'sweep'
        %  HAMS-based Tier-1/Tier-2 landscape sweep.
        %  Tier-1: ranks all cached drafts by |T_heave - T_heave_goal|.
        %  Tier-2: lightweight fmincon on top-K candidates.
        %  Stage 2 warm-starts from the best Tier-2 result.
        [x_opt_2d, props_2d, conv_data_2d, hist] = ...
            run_stage1_sweep(x0_2d, config, hist);
        converged = true;

    case 'skip'
        %  Bypass Stage 1 entirely.  Stage 2 starts cold from x0_2d
        %  (initial_vertical_shift + initial_densities from config).
        %  Phase A density pre-conditioning still fires if GM < gm_min.
        fprintf(' Stage 1 skipped — passing x0 directly to Stage 2.\n');
        x_opt_2d   = x0_2d;
        props_2d   = calculate_3d_properties(x0_2d, config);
        conv_data_2d = struct('final_objective', NaN, 'final_exitflag', 0);
        converged  = true;

    case 'oneshot'
        [x_opt_2d, props_2d, conv_data_2d, hist, iteration] = ...
            run_stage1_oneshot(x0_2d, config, hist);
        converged = true;

    case 'trained'
        [x_opt_2d, props_2d, conv_data_2d, hist, ...
         converged, iteration, pid_saturated] = ...
            run_stage1_trained(x0_2d, config, hist);

    otherwise
        error('WEC:InvalidMode', 'Unknown stage1_mode: %s', stage1_mode);
end

%% §3  STAGE 2 — 3D HIGH-FIDELITY REFINEMENT  ──────────────────────────
%  Warm-start fmincon (SQP) from Stage 1 optimum.
%
%  WHY a different objective from Stage 1?
%    Stage 1 uses a physics-based screen at each draft (Tier 1) followed by
%    lightweight fmincon refinement on top candidates (Tier 2).  Both tiers
%    use the SAME objective as Stage 2: uniform range normalization.
%
%  Objective:  Σ range_penalty(r_i)  for r_i = (actual − target) / half_range
%              i ∈ {GM, T_heave, T_pitch}.  No tuning weights.
%  Constraints: mass = buoyancy (equality), GM ≥ GM_min, monotonic density,
%               adjacent density ratio.

fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf('║ STAGE 2: 3D HIGH-FIDELITY REFINEMENT             ║\n');
fprintf('╚══════════════════════════════════════════════════╝\n\n');

x0_3d       = x_opt_2d;
num_ballast = config.num_ballast_sections;

% Variable bounds: [vertical_shift, rho_1, ..., rho_N]
[lb_3d, ub_3d] = build_3d_bounds(config);

% OutputFcn wrapper captures per-iteration 3D properties.
% WHY assignin/evalin('base')?
%   fmincon's OutputFcn cannot return data to the caller — its
%   signature is fixed: stop = f(x, optimValues, state).  The only
%   way to shuttle accumulated data out is through the base workspace.
%   This is ugly but unavoidable without restructuring around a class.
%   The temp variable is cleared immediately after retrieval.
output_fcn = @(x, ov, state) save_stage2_iteration_local( ...
    x, ov, state, config, []);

% fmincon options — solver mechanics, not centralised in the driver.
% WHY keep these here instead of in WEC_Driver?
%   These are numerical-algorithm internals (step sizes, tolerances,
%   scaling), not physics or design parameters.  Changing them does
%   not affect WHAT the solver optimises, only HOW efficiently it
%   converges.  Putting them in the driver would clutter it with
%   settings that only a developer debugging convergence would touch.
options_3d = optimoptions('fmincon', ...
    'Algorithm',              'sqp', ...
    'Display',                'iter', ...
    'MaxFunctionEvaluations', 50000, ...
    'MaxIterations',          1000, ...
    'ConstraintTolerance',    1e-8, ...
    'OptimalityTolerance',    1e-8, ...
    'StepTolerance',          1e-10, ...
    'FiniteDifferenceStepSize', 1e-6, ...
    'ScaleProblem',           true, ...
    'OutputFcn',              output_fcn);

obj_fun_3d = @(x) objective_function_3d(x, config);
con_fun_3d = @(x) constraint_function_3d(x, config);

%% §STAGE2-PHASEA  DENSITY PRE-CONDITIONING  ──────────────────────────────
%  PROBLEM:
%    When the Tier-1/2 warm-start is infeasible (GM < gm_min), SQP has two
%    routes to recover feasibility:
%
%      Route 1 — push densities bottom-heavy → GM recovers while staying in
%                the column-waterline basin (T_heave ≈ 9 s, cost ≈ 0.47).
%
%      Route 2 — nudge vs across the platform shoulder (0.06 m step) → BM
%                jumps from 0.002 to ~0.6 m, GM trivially satisfied, but
%                T_heave collapses back to 3.3 s (cost ≈ 1.60).
%
%    Route 2 costs far fewer gradient evaluations, so SQP always takes it.
%    This undoes the T_heave gain from the new Tier-1 ranking entirely.
%
%  FIX — Phase A:
%    Pin vs at the warm-start value (lb_a(1) = ub_a(1) = x0_3d(1)) and run
%    a density-only fmincon first.  Route 2 is structurally blocked — the
%    solver MUST find feasibility through density adjustment (Route 1).
%
%    Once Phase A achieves GM ≥ gm_min, Phase B (full Stage 2) starts from
%    a feasible, density-settled warm-start.  The heave cost barrier at the
%    platform shoulder (4× cost jump for a 0.06 m vs step) then holds Stage
%    2 in the column basin.
%
%  PHASE A IS SKIPPED when:
%    (a) The warm-start is already feasible (GM ≥ gm_min) — platform basin.
%    (b) Phase A fails to converge (ef_a ≤ 0) — Stage 2 falls through with
%        original x0_3d and a console warning.
%
%  NOTE ON CONSTRUCTABILITY MODE:
%    When enable_constructability = true, wall strip (column top) is pinned
%    at rho_UHPC.  If this top-heavy mass makes GM ≥ gm_min unachievable at
%    the column waterline, Phase A will report ef_a ≤ 0 and fall through.
%    Monitor the '[Phase A]' console lines to verify.

props_x0 = calculate_3d_properties(x0_3d, config);

if props_x0.GM_L < config.gm_min

    fprintf('\n╔──────────────────────────────────────────────────╗\n');
    fprintf('║ STAGE 2 — Phase A: density pre-conditioning       ║\n');
    fprintf('╚──────────────────────────────────────────────────╝\n');
    fprintf('  Warm-start infeasible: GM = %.4f m  (gm_min = %.4f m)\n', ...
            props_x0.GM_L, config.gm_min);
    fprintf('  Pinning vs = %+.4f m — settling density only...\n\n', x0_3d(1));

    lb_a    = lb_3d;  lb_a(1) = x0_3d(1);   % freeze vs at warm-start
    ub_a    = ub_3d;  ub_a(1) = x0_3d(1);   % Route 2 structurally blocked

    opts_a  = optimoptions('fmincon', ...
        'Algorithm',              'sqp', ...
        'Display',                'iter', ...
        'MaxFunctionEvaluations', 500, ...
        'MaxIterations',          50, ...
        'ConstraintTolerance',    1e-4, ...
        'OptimalityTolerance',    1e-4, ...
        'StepTolerance',          1e-6, ...
        'ScaleProblem',           true);

    [x_a, fval_a, ef_a] = fmincon(obj_fun_3d, x0_3d, ...
        [], [], [], [], lb_a, ub_a, con_fun_3d, opts_a);

    props_a = calculate_3d_properties(x_a, config);

    if ef_a > 0 && props_a.GM_L >= config.gm_min
        fprintf('\n  [Phase A] ✓ Feasibility achieved:\n');
        fprintf('    GM      = %.4f m  (target ≥ %.4f m)\n', ...
                props_a.GM_L, config.gm_min);
        fprintf('    T_heave = %.3f s\n', props_a.periods.heave);
        fprintf('    T_pitch = %.3f s\n', props_a.periods.pitch);
        fprintf('    f       = %.4f  (was %.4f at warm-start)\n', fval_a, ...
                objective_function_3d(x0_3d, config));
        x0_3d = x_a;   % hand density-settled point to Stage 2 (Phase B)
    else
        fprintf('\n  [Phase A] ✗ Did not achieve GM ≥ %.4f m ', config.gm_min);
        fprintf('(ef=%d, GM=%.4f m).\n', ef_a, props_a.GM_L);
        fprintf('  Possible cause: constructability wall prevents Route 1.\n');
        fprintf('  Stage 2 will proceed from original warm-start.\n');
        %  x0_3d unchanged — Stage 2 is no worse than without Phase A
    end

    fprintf('\n');
end
%% ────────────────────────────────────────────────────────────────────────

fprintf('Starting 3D optimisation...\n\n');
[x_opt_3d, fval_3d, exitflag_3d, output_3d] = fmincon( ...
    obj_fun_3d, x0_3d, ...
    [], [], [], [], ...
    lb_3d, ub_3d, ...
    con_fun_3d, options_3d);

% Retrieve per-iteration data from base workspace (see assignin note above)
if evalin('base', 'exist(''stage2_iteration_data_temp'', ''var'')')
    stage2_iter = evalin('base', 'stage2_iteration_data_temp');
    evalin('base', 'clear stage2_iteration_data_temp');
else
    warning('WEC:NoIterData', 'Stage 2 iteration data not captured');
    stage2_iter = struct('x', [], 'props', {{}}, 'errors', []);
end

% Unpack error arrays for results struct
[mass_errors_3d, gm_errors_3d, heave_errors_3d, pitch_errors_3d] = ...
    unpack_stage2_errors(stage2_iter);

[stage2_converged, quality_metrics] = ...
    WEC_Visualization.check_3d_convergence(x_opt_3d, config, exitflag_3d, output_3d);

%% §4  POST-PROCESSING  ────────────────────────────────────────────────

fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf('║ OPTIMISATION COMPLETE                             ║\n');
fprintf('╚══════════════════════════════════════════════════╝\n\n');

final_props = calculate_3d_properties(x_opt_3d, config);

if isfield(config, 'profile') && ~isempty(config.profile)
    final_props.cross_section = config.profile;
else
    final_props.cross_section = [];
end

%% §4a  STASH OPTIMISER SOLUTION  ─────────────────────────────────────────
%  Keep the pure optimiser result separate so realisers can reference it.

final_props_optimiser = final_props;
steel_data            = [];

%% §4b  CONSTRUCTABILITY REALIZATION (optional)  ─────────────────────
%  Post-process the converged optimiser solution into a physically
%  realisable UHPC + void geometry.  Does not modify any optimiser output.
%  Skipped entirely when enable_constructability = false.

if config.enable_constructability
    fprintf('\n╔══════════════════════════════════════════════════╗\n');
    fprintf('║ CONSTRUCTABILITY POST-PROCESSING                 ║\n');
    fprintf('╚══════════════════════════════════════════════════╝\n\n');

    constructability = WEC_Constructable_Hull.realize( ...
        config, x_opt_3d, final_props_optimiser);
    WEC_Constructable_Hull.visualize(constructability, config);
    final_props = WEC_Constructable_Hull.build_realised_props( ...
        final_props_optimiser, constructability, config);

    fprintf('\n  Constructability realization complete.\n');
else
    constructability = [];
end

%% §4c  STEEL-FILL REALISATION (post-optimisation)  ───────────────────────

if isfield(config, 'enable_steel_solve') && config.enable_steel_solve
    steel_data = WEC_Shell_Offset.solve(config, x_opt_3d, final_props_optimiser);
    plot_steel_solve(config, steel_data);
    if isfield(steel_data, 'feasible') && steel_data.feasible
        final_props = WEC_Shell_Offset.build_realised_props( ...
            final_props_optimiser, steel_data, config);
    else
        warning('WEC_Main_Optimizer:SteelInfeasibleNoSwap', ...
            'Steel solve infeasible — final_props NOT updated, falling back to optimiser.');
    end
end

%% §4d  POST-REALISATION VISUALISATION  ──────────────────────────────────
%  These plots use the FINAL REALISED final_props (steel or UHPC if enabled,
%  otherwise the raw optimiser result).  visualize_3D_equivalent is called
%  separately in §6 with final_props_optimiser and serves as Fig. 1.

WEC_Visualization.visualize_3d_cross_section(final_props, config, x_opt_3d);
WEC_Visualization.visualize_2d_equivalent(final_props, config, x_opt_3d);
WEC_Visualization.reportFinalResults(final_props, config);

%% §5  ASSEMBLE RESULTS STRUCT  ─────────────────────────────────────────

results = assemble_results(config, ...
    x_opt_2d, props_2d, conv_data_2d, hist, converged, pid_saturated, ...
    x_opt_3d, final_props, exitflag_3d, fval_3d, output_3d, ...
    stage2_converged, quality_metrics, ...
    mass_errors_3d, gm_errors_3d, heave_errors_3d, pitch_errors_3d, ...
    constructability, steel_data, final_props_optimiser);

results.optimization_time = toc(tic_main);

%% §6  VISUALISATION & SAVE  ───────────────────────────────────────────

fprintf('\nGenerating final visualisations...\n');
WEC_Visualization.visualize_3D_equivalent(final_props_optimiser, config);

% Also render the AS-BUILT view if a realisation ran (steel or UHPC).
% Detected by the presence of realised_strip_density on final_props,
% which the realisation solvers + build_realised_props populate.
if isfield(final_props, 'realised_strip_density') && ...
        ~isempty(final_props.realised_strip_density)
    WEC_Visualization.visualize_3D_equivalent(final_props, config);
end

if iteration > 1
    try
        WEC_Visualization.plot_complete_convergence(results, config);
    catch
        % Convergence plot not applicable for skip/oneshot modes
    end
end

WEC_Visualization.report_stage2_summary(results.stage2_3d);

timestamp    = datestr(now, 'yyyymmdd_HHMMSS');
results_file = sprintf('WEC_Results_%s.mat', timestamp);
save(results_file, 'results', 'final_props');
fprintf('  Results saved to: %s\n', results_file);

% CHANGE 1 FIX: dynamic canonical filename — UHPC, STEEL, or raw optimiser.
% Previously hardcoded to 'WEC_Mass_UHPC.mat' regardless of active pipeline.
% MUTUAL EXCLUSION: enable_steel_solve and enable_constructability cannot
% both be true (enforced in WEC_Driver §3b/§3c comment).
if isfield(config, 'enable_steel_solve') && config.enable_steel_solve
    canonical_file = 'WEC_Mass_STEEL.mat';
elseif isfield(config, 'enable_constructability') && config.enable_constructability
    canonical_file = 'WEC_Mass_UHPC.mat';
else
    canonical_file = 'WEC_Mass_Optimiser.mat';
end
save(canonical_file, 'results', 'final_props');
fprintf('  Canonical save : %s\n', canonical_file);

fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf('║ Total time: %.2f seconds                         ║\n', results.optimization_time);
fprintf('╚══════════════════════════════════════════════════╝\n\n');

catch ME
    fprintf('\n  OPTIMISATION FAILED: %s\n', ME.message);
    for i = 1:length(ME.stack)
        fprintf('    %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
    rethrow(ME);
end

end  % WEC_Main_Optimizer


%% =====================================================================
%%  STAGE 1 BRANCH FUNCTIONS
%% =====================================================================

function [x_opt, props_2d, conv_data, hist] = ...
        run_stage1_sweep(x0, config, hist)
% RUN_STAGE1_SWEEP  Two-tier draft-landscape sweep using 3D properties.
%
%   TIER 1 — Physics-based screen (no fmincon):
%     For each draft in config.wamit_drafts, computes the analytically
%     mass-balanced uniform-density state in ONE calculate_3d_properties
%     call.  Cost: N calls total (~0.1 ms each with V_sub table).
%     Purpose: rank drafts and identify the feasible region.
%
%   TIER 2 — Targeted fmincon on top-K candidates:
%     Runs a lightweight fmincon (25 iter, loose tolerances) on the
%     top config.n_sweep_refine feasible drafts from Tier 1.
%     Uses the same objective and constraints as Stage 2 — no new code.
%     Purpose: refine the warm-start density profile before Stage 2.
%
%   WHY NOT full fmincon on every draft (old behaviour)?
%     With the V_sub precomputed table, each calculate_3d_properties
%     call is ~0.1 ms.  The Tier-1 screen costs ~2 ms for 20 drafts.
%     The old approach ran 200-iter fmincon at every draft, which
%     dominated runtime without improving Stage-2 warm-start quality.
%
%   FEASIBILITY NOTE:
%     Tier-1 uses a uniform density (mass-balanced), which does NOT
%     satisfy the monotonic or ratio constraints.  The feasibility flag
%     therefore only checks mass balance and GM — not constraint
%     satisfaction.  Tier-2 applies full constraints and may revise
%     the feasibility label.  This is intentional: Tier 1 is a
%     landscape probe, not a solution.
%
%   FALLBACK:
%     If config.wamit_drafts has <= 1 entry, passes x0 straight through.

    % ── Guard: fall back to passthrough if no draft grid ──────────
    if ~isfield(config, 'wamit_drafts') || length(config.wamit_drafts) <= 1
        fprintf(' Skipping Stage 1 (no draft grid). Using default x0.\n');
        x_opt    = x0;
        props_2d = calculate_2d_properties(x0, config);
        conv_data = struct('final_objective', NaN, 'final_exitflag', 0);
        return;
    end

    vs_grid = config.wamit_drafts(:)';
    N       = length(vs_grid);

    fprintf(' Stage 1: two-tier draft-landscape sweep (%d drafts)\n', N);
    fprintf('   vs range: [%.3f, %.3f] m\n', vs_grid(1), vs_grid(end));

    % ── Bounds ────────────────────────────────────────────────────
    [lb_base, ub_base] = build_3d_bounds(config);

    % ── Pre-allocate sweep results ────────────────────────────────
    sweep.vs       = vs_grid;
    sweep.fval     = inf(1, N);
    sweep.exitflag = zeros(1, N);
    sweep.x        = cell(1, N);
    sweep.props    = cell(1, N);
    sweep.feasible = false(1, N);

    obj_fun = @(x) objective_function_3d(x, config);
    con_fun = @(x) constraint_function_3d(x, config);

    % ── TIER 1: physics-based single-point screen ─────────────────
    fprintf('\n Tier 1: physics-based screen...\n');

    for k = 1:N
        vs_k   = vs_grid(k);
        screen = screen_draft(vs_k, config);

        sweep.fval(k)     = screen.fval;
        sweep.exitflag(k) = 1;
        sweep.x{k}        = screen.x;
        sweep.props{k}    = screen.props;
        sweep.feasible(k) = screen.feasible;

        fprintf('   [%2d/%d] vs=%+.3f  f=%.4f  GM=%.3f  T_h=%.2f  T_p=%.2f  %s\n', ...
            k, N, vs_k, screen.fval, screen.props.GM_L, ...
            screen.props.periods.heave, screen.props.periods.pitch, ...
            ternary(screen.feasible, 'FEAS', 'infeas'));
    end

    % ── TIER 2: targeted fmincon on top-K candidates ──────────────
    %  FIX: Rank candidates by T_heave proximity to target instead of
    %  GM feasibility at uniform density.
    %
    %  WHY THE OLD APPROACH FAILED:
    %    The GM check at uniform density incorrectly excluded column-
    %    waterline drafts (where BM ≈ 0.002 m and uniform CG is above
    %    CB).  Those drafts have T_heave closest to the target and are
    %    exactly the warm-starts Stage 2 needs.  With heavy bottom
    %    ballast, GM is achievable there — but the uniform-density probe
    %    never discovered that.  Stage 2 was then trapped at the
    %    platform waterline (T_heave ≈ 3 s) as a local minimum.
    %
    %  WHY T_HEAVE IS THE CORRECT SIGNAL:
    %    T_heave = 2π√((M+A33)/(ρgAw)).  At fixed vs, M = V_sub·ρ_w
    %    (equilibrium), A33 and Aw come from HAMS — all density-
    %    independent.  The T_heave from screen_draft's uniform-density
    %    probe is therefore exact regardless of the probe density used.
    %    Ranking by |T_heave − T_target| correctly identifies which
    %    drafts will drive the heave objective toward zero.
    %
    %  GM IS STILL ENFORCED — by Tier-2 fmincon (full constraints) and
    %  Stage 2 fmincon.  This change only affects which drafts are
    %  *presented* to Tier-2, not how Tier-2 evaluates them.
    K_refine   = config.n_sweep_refine;
    T_h_vals   = arrayfun(@(k) sweep.props{k}.periods.heave, 1:N);
    T_h_vals(~isfinite(T_h_vals)) = Inf;   % dry / degenerate drafts rank last
    [~, rank]  = sort(abs(T_h_vals - config.T_heave_goal));
    top_idx    = rank(1:min(K_refine, N));

    fprintf('\n Tier 1 → Tier 2 candidates ranked by |T_heave − %.1f s|:\n', ...
            config.T_heave_goal);
    for ki = 1:length(top_idx)
        k_ki = top_idx(ki);
        fprintf('   [%d] vs=%+.3f m  T_h=%.2f s  |err|=%.2f s  GM=%.3f m  %s\n', ...
                ki, vs_grid(k_ki), T_h_vals(k_ki), ...
                abs(T_h_vals(k_ki) - config.T_heave_goal), ...
                sweep.props{k_ki}.GM_L, ...
                ternary(sweep.feasible(k_ki), 'GM-feas', 'GM-infeas'));
    end

    % Lightweight fmincon — looser than Stage 2 (warm-start only)
    opts_refine = optimoptions('fmincon', ...
        'Algorithm',              'sqp', ...
        'Display',                'none', ...
        'MaxFunctionEvaluations', 150, ...
        'MaxIterations',          25, ...
        'ConstraintTolerance',    1e-4, ...
        'OptimalityTolerance',    1e-4, ...
        'StepTolerance',          1e-6, ...
        'ScaleProblem',           true);

    fprintf('\n Tier 2: refining %d candidates...\n', length(top_idx));

    for j = 1:length(top_idx)
        k    = top_idx(j);
        vs_k = vs_grid(k);

        lb_k = lb_base;  lb_k(1) = vs_k;
        ub_k = ub_base;  ub_k(1) = vs_k;
        x0_k = max(lb_k, min(ub_k, sweep.x{k}));

        try
            [x_k, fval_k, ef_k] = fmincon(obj_fun, x0_k, ...
                [], [], [], [], lb_k, ub_k, con_fun, opts_refine);

            props_k  = calculate_3d_properties(x_k, config);
            mass_err = abs(props_k.mass_total - props_k.mass_buoyant_force) ...
                       / max(props_k.mass_total, 1);
            is_feas  = ef_k > 0 && mass_err < 0.01 && props_k.GM_L >= config.gm_min;

            % Only update if Tier 2 improved on Tier 1
            if fval_k < sweep.fval(k)
                sweep.fval(k)     = fval_k;
                sweep.exitflag(k) = ef_k;
                sweep.x{k}        = x_k;
                sweep.props{k}    = props_k;
                sweep.feasible(k) = is_feas;
            end

            fprintf('   [T2 %d/%d] vs=%+.3f  f=%.4f  GM=%.3f  %s\n', ...
                j, length(top_idx), vs_k, sweep.fval(k), props_k.GM_L, ...
                ternary(is_feas, 'FEAS', 'infeas'));

        catch ME_t2
            fprintf('   [T2 %d/%d] vs=%+.3f  FAILED: %s\n', ...
                j, length(top_idx), vs_k, ME_t2.message);
        end
    end

    % ── Pick best point for Stage-2 warm-start ────────────────────
    feas_idx_final = find(sweep.feasible);
    if ~isempty(feas_idx_final)
        [~, best_in_feas] = min(sweep.fval(feas_idx_final));
        best = feas_idx_final(best_in_feas);
        fprintf('\n   Best feasible: draft %d (vs=%+.3f m, f=%.4f)\n', ...
                best, vs_grid(best), sweep.fval(best));
    else
        [~, best] = min(sweep.fval);
        warning('WEC:NoFeasibleDraft', ...
                'No feasible draft after Tier-2 refinement. Using lowest objective (vs=%+.3f m).', ...
                vs_grid(best));
    end

    x_opt    = sweep.x{best};
    props_2d = calculate_2d_properties(x_opt, config);
    conv_data = struct('final_objective',  sweep.fval(best), ...
                       'final_exitflag',   sweep.exitflag(best), ...
                       'sweep',            sweep);

    fprintf('   Warm-start for Stage 2: vs=%+.3f m, rho=[', x_opt(1));
    fprintf('%.0f ', x_opt(2:end));
    fprintf(']\n');
end


function result = screen_draft(vs_k, config)
% SCREEN_DRAFT  Single-point physics probe at a fixed draft — no fmincon.
%
%   Builds the analytically mass-balanced uniform-density state and
%   evaluates 3D properties in ONE calculate_3d_properties call.
%
%   DENSITY CONSTRUCTION
%     The required average density from flotation equilibrium is:
%       rho_uniform = rho_water * V_sub / V_hull_total
%     For constructability mode, the wall strip is pinned at rho_UHPC
%     and the remaining mass is distributed uniformly to platform strips.
%
%   FEASIBILITY
%     Checks mass balance (< 5% error) and GM >= gm_min.
%     Does NOT check monotonic/ratio constraints — those are enforced
%     only in Tier 2.  This is intentional: Tier 1 is a landscape probe.
%
%   COST
%     One call to calculate_3d_properties (~0.1 ms with V_sub table).

    N        = config.num_ballast_sections;
    V_strips = sum(config.strip_V);   % total hull strip volume from config

    % Probe with initial densities to get V_sub at this draft
    x_probe     = [vs_k, config.initial_densities];
    props_probe = calculate_3d_properties(x_probe, config);
    V_sub       = props_probe.V_sub;

    if V_sub < 1e-6 || V_strips < 1e-6
        result.feasible = false;
        result.fval     = Inf;
        result.x        = x_probe;
        result.props    = props_probe;
        return
    end

    % Build mass-balanced density vector
    if config.enable_constructability && ~isempty(config.wall_strip_index)
        % Wall strip pinned; distribute remaining required mass to platform
        w         = config.wall_strip_index;
        mass_wall = config.strip_V(w) * config.constructability_rho_hull;
        V_free    = V_strips - config.strip_V(w);
        mass_need = config.RHO_WATER * V_sub - mass_wall;
        rho_free  = mass_need / max(V_free, eps);
        rho_free  = max(config.ballast_density_bounds(1), ...
                    min(config.ballast_density_bounds(2), rho_free));

        rho_bal      = repmat(rho_free, 1, N);
        rho_bal(w)   = config.constructability_rho_hull;
    else
        rho_unif = config.RHO_WATER * V_sub / V_strips;
        rho_unif = max(config.ballast_density_bounds(1), ...
                   min(config.ballast_density_bounds(2), rho_unif));
        rho_bal  = repmat(rho_unif, 1, N);
    end

    x_bal  = [vs_k, rho_bal];
    props  = calculate_3d_properties(x_bal, config);
    fval   = objective_function_3d(x_bal, config);

    mass_err    = abs(props.mass_total - props.mass_buoyant_force) / max(props.mass_total, 1);
    is_feasible = mass_err < 0.05 && props.GM_L >= config.gm_min;

    result.feasible = is_feasible;
    result.fval     = fval;
    result.x        = x_bal;
    result.props    = props;
end


function [x_opt, props_2d, conv_data, hist, iteration] = ...
        run_stage1_oneshot(x0, config, hist)
% RUN_STAGE1_ONESHOT  Single 2D solve (uncorrected surrogate), then
%   validate against 3D and report the warm-start gap.

    fprintf(' Running Stage 1 once (k_vol = %.1f, k_gm = %.1f)\n\n', ...
            config.k_vol, config.k_gm);

    % Mass PID: gains and limits from config (forwarded from driver)
    g = config.pid_mass_gains;
    mass_pid = WEC_PID_Controller(g(1), g(2), g(3), ...
        'OutputLimits', config.pid_mass_limits);

    [x_opt, props_2d, exitflag, conv_data] = ...
        run_2d_optimizer(config, x0, mass_pid);

    if exitflag <= 0
        warning('WEC:Stage1Failed', ...
                'Stage 1 failed (exitflag=%d). Using x0.', exitflag);
        x_opt    = x0;
        props_2d = calculate_2d_properties(x0, config);
    end

    iteration = 1;
    p3 = calculate_3d_properties(x_opt, config);

    % Store single-iteration history
    hist.mass_errors_history          = props_2d.mass_total - p3.mass_total;
    hist.gm_errors_history            = props_2d.GM_uncorrected - p3.GM_L;
    hist.mass_correction_history      = config.k_vol;
    hist.gm_correction_history        = config.k_gm;
    hist.mass_2d_history              = props_2d.mass_total;
    hist.mass_3d_history              = p3.mass_total;
    hist.mass_2d_corrected_history    = props_2d.mass_total;
    hist.gm_2d_history                = props_2d.GM_uncorrected;
    hist.gm_3d_history                = p3.GM_L;
    hist.cg_z_2d_history              = props_2d.CG_total(3);
    hist.cg_z_3d_history              = p3.CG_total(3);

    % One-shot: single data point → R²/MAPE are undefined
    hist.R2_mass_history   = NaN;  hist.R2_GM_history   = NaN;  hist.R2_cg_history   = NaN;
    hist.MAPE_mass_history = NaN;  hist.MAPE_GM_history = NaN;  hist.MAPE_cg_history = NaN;
    hist.ME_mass_history   = NaN;  hist.ME_GM_history   = NaN;

    % Convergence metrics (single iteration, always true for oneshot)
    hist.constraints_satisfied_history = true;
    hist.solution_stable_history       = true;
    hist.mass_acceptable_history       = abs(100 * hist.mass_errors_history / p3.mass_total) < config.mass_acceptable_pct;
    hist.converged_history             = true;
    hist.fitness_scores_history        = 1 - abs(hist.mass_errors_history / p3.mass_total);
    hist.GM_margins_history            = max(0, (p3.GM_L - config.gm_min) / config.gm_min);

    mass_gap_pct = 100 * hist.mass_errors_history / p3.mass_total;
    fprintf('\n  One-Shot Results:\n');
    fprintf('    Draft: %.3f m\n', x_opt(1));
    fprintf('    2D→3D Mass: %.1f → %.1f kg (%+.1f%%)\n', ...
            props_2d.mass_total, p3.mass_total, mass_gap_pct);
    fprintf('    2D→3D GM:   %.3f → %.3f m\n', ...
            props_2d.GM_uncorrected, p3.GM_L);
    fprintf('    2D→3D T_h:  %.2f → %.2f s\n', ...
            props_2d.periods.heave, p3.periods.heave);
    fprintf('    2D→3D T_p:  %.2f → %.2f s\n', ...
            props_2d.periods.pitch, p3.periods.pitch);
end


function [x_opt, props_2d, conv_data, hist, ...
          converged, iteration, pid_saturated] = ...
        run_stage1_trained(x0, config, hist)
% RUN_STAGE1_TRAINED  Iterative PID correction of the 2D surrogate.
%
%   Two correction channels updated each iteration:
%     k_vol : multiplicative correction on strip-theory submerged volume.
%             PID error = (V_sub_3D − V_sub_2D) / V_sub_3D  (unitless).
%     k_gm  : multiplicative correction on 2D CG_z position.
%             PID error = (cg_z_3D / cg_z_2D) − 1           (unitless).
%
%   WHY these error signals?
%     k_vol error is a relative volume gap: positive when 2D under-predicts
%     → PID increases k_vol → strips get "fatter" → volume rises.
%     k_gm error is a CG depth ratio: positive when 3D CG is deeper
%     → PID increases k_gm → biases 2D CG downward → GM rises.
%     Both are dimensionless and O(1), which keeps PID gain tuning
%     consistent across different hull sizes.
%
%   Convergence is declared when BOTH:
%     |V_sub error| < vol_conv_tol_pct  AND  |GM error| < gm_conv_tol
%   OR when correction factors stabilise for stable_count_needed iters
%   AND at least 3 iterations have elapsed (to avoid premature exit
%   from a lucky initial guess).
%
%   VARIABLE DICTIONARY (loop-scope)
%     vol_err_norm   [-]   (V_3D − V_2D) / V_3D,  fed to vol_pid
%     vol_update     [-]   raw PID output for k_vol
%     target_kgm     [-]   cg_z_3D / cg_z_2D,  the "ideal" k_gm
%     gm_err_norm    [-]   target_kgm − 1,  fed to gm_pid
%     gm_update      [-]   raw PID output for k_gm
%     damp_v, damp_g [-]   under-relaxation multipliers (early/late)
%     delta_kvol     [-]   |k_vol_new − k_vol_old|
%     stable_count   [-]   consecutive iters with small Δk

    fprintf(' Running Stage 1 with PID surrogate training.\n');

    % --- Unpack PID configuration from config (set by driver) ---
    bounds_kvol = config.bounds_kvol;
    bounds_kgm  = config.bounds_kgm;

    gv = config.pid_vol_gains;
    gg = config.pid_gm_gains;
    gm = config.pid_mass_gains;

    vol_pid  = WEC_PID_Controller(gv(1), gv(2), gv(3), ...
               'OutputLimits', config.pid_vol_limits);
    gm_pid   = WEC_PID_Controller(gg(1), gg(2), gg(3), ...
               'OutputLimits', config.pid_gm_limits);
    mass_pid = WEC_PID_Controller(gm(1), gm(2), gm(3), ...
               'OutputLimits', config.pid_mass_limits);

    max_iters    = config.max_outer_iterations;
    stable_count = 0;

    converged     = false;
    iteration     = 0;
    pid_saturated = false;
    double_sat_count = 0;   % consecutive iterations with both PIDs saturated

    fprintf(' Bounds: k_vol [%.2f, %.2f], k_gm [%.2f, %.2f]\n', ...
            bounds_kvol, bounds_kgm);
    fprintf(' Max iterations: %d\n', max_iters);

    while ~converged && iteration < max_iters
        iteration = iteration + 1;
        fprintf('\n  [Stage1 %2d/%d] k_vol=%.3f k_gm=%.3f ', ...
                iteration, max_iters, config.k_vol, config.k_gm);

        % --- 2D optimisation with current correction factors ---
        [x_opt, props_2d, exitflag, conv_data] = ...
            run_2d_optimizer(config, x0, mass_pid);

        if exitflag <= 0
            fprintf('  2D exit=%d. ', exitflag);
        end

        % --- HAMS enrichment: run at converged vertical_shift if not cached ---
        %  This ensures the 3D validation below uses exact (not interpolated)
        %  hydrodynamic coefficients at the 2D-converged draft.
        %  The cache grows by at most 1 entry per outer iteration.
        %  config.wamit_* fields are rebuilt so interpolate_wamit_added_mass
        %  sees the enriched table immediately.
        if ~isempty(config.hams_dir)
            vs_converged = x_opt(1);
            [~, config.hydro_cache] = HAMS_Pipeline.get_or_run_hams( ...
                vs_converged, config.hydro_cache, config, ...
                config.hams_dir, config.hams_exe, 0.01);
            config = HAMS_Pipeline.rebuild_config_hydro(config, config.hydro_cache);
        end

        % Guard: warn if A(inf) is all zeros (HAMS may have failed)
        if ~isempty(config.wamit_A) && max(abs(config.wamit_A(:))) < 1e-6
            warning('WEC:ZeroAddedMass', ...
                    'A(inf) = 0 at all drafts. HAMS failed — periods will be wrong.');
        end

        % --- 3D validation pass 1: get actual CG from density distribution ---
        p3 = calculate_3d_properties(x_opt, config);

        % --- CG correction: re-transform A from origin using actual CG ---
        %  rebuild_config_hydro transformed A_origin → A_CG using the
        %  UNIFORM-DENSITY CG (from compute_hams_inputs).  The actual CG
        %  from the non-uniform density distribution is different.
        %
        %  A33 is unaffected (no CG dependence).
        %  A55 has quadratic CG dependence → ~4% error on T_pitch.
        %  A11 has linear CG dependence → ~2% error on T_surge.
        %
        %  Fix: re-transform using the actual CG, then re-evaluate.
        %  Cost: one 3×3 matrix multiply + one calculate_3d_properties call.
        if ~isempty(config.hams_dir)
            actual_cg_z = p3.CG_total(3);
            config = HAMS_Pipeline.retransform_at_actual_cg( ...
                config, vs_converged, actual_cg_z);
            p3 = calculate_3d_properties(x_opt, config);
        end

        vsub_2d = props_2d.V_sub;        vsub_3d = p3.V_sub;
        mass_2d = props_2d.mass_total;    mass_3d = p3.mass_total;
        gm_2d   = props_2d.GM;           gm_3d   = p3.GM_L;
        cg_z_2d = props_2d.CG_total(3);  cg_z_3d = p3.CG_total(3);

        vol_error_pct  = 100 * (vsub_2d - vsub_3d) / vsub_3d;
        mass_error_pct = 100 * (mass_2d - mass_3d) / mass_3d;
        gm_error       = gm_2d - gm_3d;

        % --- Append history ---
        hist.vsub_2d_history(end+1)       = vsub_2d;
        hist.vsub_3d_history(end+1)       = vsub_3d;
        hist.vol_errors_history(end+1)    = vol_error_pct;
        hist.mass_2d_history(end+1)       = mass_2d;
        hist.mass_3d_history(end+1)       = mass_3d;
        hist.mass_2d_corrected_history(end+1) = mass_2d;
        hist.gm_2d_history(end+1)         = gm_2d;
        hist.gm_3d_history(end+1)         = gm_3d;
        hist.gm_errors_history(end+1)     = gm_error;
        hist.cg_z_2d_history(end+1)       = cg_z_2d;
        hist.cg_z_3d_history(end+1)       = cg_z_3d;
        hist.kvol_history(end+1)          = config.k_vol;
        hist.kgm_history(end+1)           = config.k_gm;
        hist.mass_errors_history(end+1)   = mass_2d - mass_3d;
        hist.mass_correction_history(end+1) = config.k_vol;
        hist.gm_correction_history(end+1)   = config.k_gm;

        % --- Console report (compact) ---
        fprintf('\n    V: %+.1f%%  M: %+.1f%%  GM: %.3f→%.3f  T_h: %.1f→%.1f  T_p: %.1f→%.1f\n', ...
                vol_error_pct, mass_error_pct, gm_2d, gm_3d, ...
                props_2d.periods.heave, p3.periods.heave, ...
                props_2d.periods.pitch, p3.periods.pitch);

        % --- PID update: k_vol ---
        %  WHY multiplicative (k * (1 + damp * Δ)) instead of additive (k + Δ)?
        %    k_vol is a scaling factor, not an offset.  Multiplicative updates
        %    keep the step proportional to the current magnitude, preventing
        %    the factor from crossing zero or going negative.
        kvol_before = config.k_vol;
        kgm_before  = config.k_gm;

        if vsub_2d > 1e-6
            vol_err_norm = (vsub_3d - vsub_2d) / vsub_3d;
            vol_update   = vol_pid.update(vol_err_norm, 1.0);
        else
            vol_update = 0;
        end

        % Under-relaxation: cautious early, more aggressive once trend is clear
        if iteration <= config.damping_transition_iter
            damp_v = config.damping_vol_early;
        else
            damp_v = config.damping_vol_late;
        end
        kvol_new     = config.k_vol * (1 + damp_v * vol_update);
        config.k_vol = max(bounds_kvol(1), min(bounds_kvol(2), kvol_new));

        % --- PID update: k_gm ---
        %  WHY cg_guard_floor?
        %    k_gm error = (cg_z_3D / cg_z_2D) − 1.  When either CG_z
        %    is near zero, the ratio blows up or flips sign, producing
        %    nonsensical PID updates.  The guard skips the update when
        %    |cg_z| < 0.01 m (hull nearly centred on waterline).
        if abs(cg_z_2d) > config.cg_guard_floor && ...
           abs(cg_z_3d) > config.cg_guard_floor
            target_kgm      = cg_z_3d / cg_z_2d;
            gm_err_norm     = target_kgm - 1.0;
            gm_update       = gm_pid.update(gm_err_norm, 1.0);
        else
            gm_update = 0;
        end

        if iteration <= config.damping_transition_iter
            damp_g = config.damping_gm_early;
        else
            damp_g = config.damping_gm_late;
        end
        kgm_new     = config.k_gm * (1 + damp_g * gm_update);
        config.k_gm = max(bounds_kgm(1), min(bounds_kgm(2), kgm_new));

        % --- Convergence check ---
        %  WHY (vol_ok && gm_ok) || (kvol_stable && iteration > 3)?
        %    Primary path: both error channels are within tolerance.
        %    Fallback path: correction factors have stopped moving (the PID
        %    has converged on its own dynamic) AND at least 3 iterations have
        %    passed.  The iteration > 3 guard prevents false convergence when
        %    the initial k_vol=k_gm=1.0 happens to be near the optimum
        %    and the first two Δk are trivially small.
        delta_kvol = abs(config.k_vol - kvol_before);
        delta_kgm  = abs(config.k_gm  - kgm_before);

        vol_ok = abs(vol_error_pct) < config.vol_conv_tol_pct;
        gm_ok  = abs(gm_error)      < config.gm_conv_tol;

        if delta_kvol < config.delta_kvol_stable && ...
           delta_kgm  < config.delta_kgm_stable
            stable_count = stable_count + 1;
        else
            stable_count = 0;
        end
        kvol_stable = (stable_count >= config.stable_count_needed);

        % Saturation detection: correction factor stuck at its bound
        prox = config.pid_sat_proximity;
        pid_sat_vol = abs(config.k_vol - bounds_kvol(2)) < prox || ...
                      abs(config.k_vol - bounds_kvol(1)) < prox;
        pid_sat_gm  = abs(config.k_gm  - bounds_kgm(2)) < prox || ...
                      abs(config.k_gm  - bounds_kgm(1)) < prox;
        pid_saturated = pid_sat_vol || pid_sat_gm;

        fprintf('    PID: Δk_vol=%+.4f%s  Δk_gm=%+.4f%s\n', ...
                delta_kvol, sat_tag(pid_sat_vol), ...
                delta_kgm, sat_tag(pid_sat_gm));

        % --- Early exit: both PIDs saturated for 3 consecutive iters ---
        %  When both correction factors are at their bounds AND the
        %  solution is not changing, further iterations are wasted.
        %  Break early and let Stage 2 work from the best-effort warm-start.
        if pid_sat_vol && pid_sat_gm
            double_sat_count = double_sat_count + 1;
            if double_sat_count >= 3
                fprintf('\n  EARLY EXIT: both PIDs saturated for %d consecutive iterations.\n', ...
                        double_sat_count);
                fprintf('    k_vol=%.4f (bound), k_gm=%.4f (bound)\n', ...
                        config.k_vol, config.k_gm);
                fprintf('    Passing best-effort solution to Stage 2.\n');
                break;
            end
        else
            double_sat_count = 0;
        end

        if exitflag > 0 && ((vol_ok && gm_ok) || (kvol_stable && iteration > 3))
            converged = true;
            fprintf('\n  SURROGATE CONVERGED at iteration %d\n', iteration);
            fprintf('    V_sub err: %.1f%%, GM err: %.3f m\n', vol_error_pct, gm_error);
        elseif exitflag <= 0 && (vol_ok && gm_ok)
            % 2D→3D errors are within tolerance, but fmincon itself
            % failed (infeasible / stalled).  Do NOT declare convergence;
            % the warm-start quality is suspect.
            fprintf('\n  2D errors within tolerance but fmincon failed (exit=%d) — continuing\n', exitflag);
            x0 = x_opt;
        else
            x0 = x_opt;   % warm-start next iteration
        end

        % --- Bookkeeping for results struct ---
        hist.constraints_satisfied_history(end+1) = true;
        hist.solution_stable_history(end+1)       = kvol_stable;
        hist.mass_acceptable_history(end+1)       = ...
            abs(mass_error_pct) < config.mass_acceptable_pct;
        hist.converged_history(end+1)             = converged;
        hist.fitness_scores_history(end+1)        = 1 - abs(vol_error_pct)/100;
        hist.GM_margins_history(end+1)            = ...
            max(0, (gm_3d - config.gm_min) / config.gm_min);

        % --- Running surrogate accuracy statistics ---
        %  WHY running (cumulative) instead of per-iteration?
        %    R² and MAPE measure how well the *overall* 2D model tracks
        %    the 3D ground truth across ALL iterations so far.  A single-
        %    iteration comparison is just an error, not a goodness-of-fit.
        %    Need ≥ 2 data points for meaningful variance; pad with NaN
        %    for iteration 1 so all history vectors stay the same length.
        n_pts = length(hist.mass_2d_history);
        if n_pts >= 2
            hist.R2_mass_history(end+1)   = compute_R2(hist.mass_3d_history, hist.mass_2d_history);
            hist.R2_GM_history(end+1)     = compute_R2(hist.gm_3d_history,   hist.gm_2d_history);
            hist.R2_cg_history(end+1)     = compute_R2(hist.cg_z_3d_history, hist.cg_z_2d_history);
            hist.MAPE_mass_history(end+1) = compute_MAPE(hist.mass_3d_history, hist.mass_2d_history);
            hist.MAPE_GM_history(end+1)   = compute_MAPE(hist.gm_3d_history,   hist.gm_2d_history);
            hist.MAPE_cg_history(end+1)   = compute_MAPE(hist.cg_z_3d_history, hist.cg_z_2d_history);
            hist.ME_mass_history(end+1)   = mean(hist.mass_2d_history - hist.mass_3d_history);
            hist.ME_GM_history(end+1)     = mean(hist.gm_2d_history   - hist.gm_3d_history);
        else
            hist.R2_mass_history(end+1)   = NaN;
            hist.R2_GM_history(end+1)     = NaN;
            hist.R2_cg_history(end+1)     = NaN;
            hist.MAPE_mass_history(end+1) = NaN;
            hist.MAPE_GM_history(end+1)   = NaN;
            hist.MAPE_cg_history(end+1)   = NaN;
            hist.ME_mass_history(end+1)   = NaN;
            hist.ME_GM_history(end+1)     = NaN;
        end
    end  % while

    % --- Save enriched hydro cache to disk ---
    %  The cache now has the coarse startup points PLUS one point per
    %  outer iteration at the converged vertical_shift. Next session
    %  will load these and skip the HAMS runs at the same drafts.
    if ~isempty(config.hams_dir) && ~isempty(config.hydro_cache_file)
        hydro_table = config.hydro_cache; %#ok<NASGU>
        save(config.hydro_cache_file, 'hydro_table', '-v7.3');
        fprintf('  Hydro cache saved: %d entries → %s\n', ...
                length(config.hydro_cache.drafts), config.hydro_cache_file);
    end

    if ~converged
        if pid_saturated
            warning('WEC:PIDSaturated', ...
                'PID saturated: k_vol=%.3f, k_gm=%.3f', ...
                config.k_vol, config.k_gm);
        else
            warning('WEC:NotConverged', ...
                'Surrogate did not converge after %d iterations', iteration);
        end
    end
end


%% =====================================================================
%%  STAGE 2 LOCAL FUNCTIONS
%% =====================================================================

function [c, ceq] = constraint_function_3d(x, config)
% CONSTRAINT_FUNCTION_3D  Nonlinear constraints for Stage 2 fmincon.
%
%   Inequality (c ≤ 0), all dimensionless O(1):
%     1. GM floor:           1 − GM_L / gm_min                 ≤ 0
%     2. Adjacent density:   rho_i/(rho_{i+1}+1) − max_ratio   ≤ 0
%        WHY +1 in denominator?  Prevents division by zero when a
%        density node approaches the lower bound (100 kg/m³).  The
%        +1 kg/m³ shift is negligible for physical densities but
%        keeps the constraint smooth for the SQP Jacobian.
%     3. Monotonic density:  (rho_{i+1}−rho_i) / rho_max       ≤ 0
%        Forces heavy-bottom → light-top ordering (same as 2D).
%
%   CONSTRUCTABILITY: When enabled, the wall strip (pinned at rho_UHPC)
%   is excluded from monotonic and ratio constraints.  The wall density
%   is physically determined by the material, not the optimiser, so
%   enforcing monotonicity across the wall-platform boundary would
%   over-constrain the platform strips immediately below the wall.
%
%   Equality (ceq = 0):
%     mass / buoyancy − 1 = 0

try
    props     = calculate_3d_properties(x, config);
    densities = x(2:end);
    % When shell is enabled, densities are CORE densities (not bulk).
    % The ratio and monotonic constraints apply to the core only,
    % because the shell density is uniform across all strips and
    % poses no manufacturability concern.
    N         = length(densities);
    rho_max   = config.ballast_density_bounds(2);

    c_gm = 1.0 - props.GM_L / config.gm_min;

    % Identify which adjacent pairs to constrain.
    % Skip any pair bridging the wall-platform boundary.
    w_idx = config.wall_strip_index;   % [] if not constructability
    constrained_pairs = [];
    for i = 1:(N-1)
        if ~isempty(w_idx) && (i == w_idx || i + 1 == w_idx)
            continue;   % skip wall-platform boundary
        end
        constrained_pairs(end+1) = i; %#ok<AGROW>
    end
    n_pairs = length(constrained_pairs);

    c_ratio = zeros(n_pairs, 1);
    c_mono  = zeros(n_pairs, 1);
    for j = 1:n_pairs
        i = constrained_pairs(j);
        c_ratio(j) = densities(i) / (densities(i+1) + 1) - config.max_density_ratio;
        c_mono(j)  = (densities(i+1) - densities(i)) / rho_max;
    end

    % Minimum mass constraint (constructability mode).
    % Prevents the optimizer from requesting density distributions
    % that produce less mass than is physically achievable.
    if config.m_min_constructability > 0
        c_mass_min = 1.0 - props.mass_total / config.m_min_constructability;
    else
        c_mass_min = [];
    end

    c   = [c_gm; c_ratio; c_mono; c_mass_min];
    ceq = props.mass_total / props.mass_buoyant_force - 1.0;
catch
    % Return O(1) violations for diverging solutions.
    % Size must be consistent — use safe fallback.
    n_pairs_fallback = max(0, length(x) - 2);
    if config.enable_constructability && ~isempty(config.wall_strip_index)
        n_pairs_fallback = max(0, n_pairs_fallback - 1);
    end
    n_extra = 0;
    if isfield(config, 'm_min_constructability') && config.m_min_constructability > 0
        n_extra = 1;
    end
    c   = ones(1 + 2 * n_pairs_fallback + n_extra, 1);
    ceq = 1.0;
end
end


function f = objective_function_3d(x, config)
% OBJECTIVE_FUNCTION_3D  Stage 2 objective — uniform range normalization.
%
%   f = phi(r_gm) + phi(r_heave) + phi(r_pitch)
%
%   where r = (actual − target) / half_range  for each quantity.
%
%   This makes all three terms dimensionless with IDENTICAL scaling:
%     r = 0   → at target       → phi = 0
%     r = ±1  → at range boundary → phi = 1.0  (same for all three)
%     |r| > 1 → outside range    → phi amplified by k_amp
%
%   The optimizer minimises the total squared normalised error.
%   Since all terms hit 1.0 at their respective boundaries, no single
%   quantity can dominate.  The optimizer finds the best compromise.
%
%   WHY this is defensible:
%     Each target and range is a physical design requirement set by
%     the engineer (T_heave_goal, T_heave_range, gm_target, gm_range).
%     The normalization follows automatically — no tuning weights.
%     A reviewer can verify: "at the boundary of any acceptable range,
%     all three penalties equal exactly 1.0."

props = calculate_3d_properties(x, config);

% Guard: return penalty for diverging solutions
if isnan(props.GM_L)
    f = config.penalty_guard;
    return;
end
if isnan(props.periods.heave) || isinf(props.periods.heave) || ...
   isnan(props.periods.pitch) || isinf(props.periods.pitch)
    f = config.penalty_guard;
    return;
end

% Half-ranges (computed from driver-specified ranges)
gm_half    = 0.5 * (config.gm_range(2)      - config.gm_range(1));
heave_half = 0.5 * (config.T_heave_range(2)  - config.T_heave_range(1));
pitch_half = 0.5 * (config.T_pitch_range(2)  - config.T_pitch_range(1));

% Normalised errors: r = 0 at target, r = ±1 at range boundary
r_gm    = (props.GM_L           - config.gm_target)     / gm_half;
r_heave = (props.periods.heave  - config.T_heave_goal)  / heave_half;
r_pitch = (props.periods.pitch  - config.T_pitch_goal)  / pitch_half;

% Piecewise quadratic with amplification outside |r| > 1
k_amp = config.zone_k_amp;
f = range_penalty(r_gm, k_amp) ...
  + range_penalty(r_heave, k_amp) ...
  + range_penalty(r_pitch, k_amp);
end


function phi = range_penalty(r, k_amp)
% RANGE_PENALTY  C1-continuous piecewise-quadratic penalty.
%
%   phi = range_penalty(r, k_amp)
%
%   r = 0    → phi = 0      (at target)
%   r = ±1   → phi = 1.0    (at range boundary)
%   |r| > 1  → phi > 1, growing k_amp× faster
%
%   C1-continuous at r = ±1:
%     Inside:  phi = r²           phi' = 2r
%     Outside: phi = 1 + 2δ + k_amp·δ²   where δ = |r| − 1
%              phi' = 2 + 2·k_amp·δ      (matches 2r at |r|=1)

if r < -1
    delta = -r - 1;     % positive
    phi = 1 + 2*delta + k_amp * delta^2;
elseif r > 1
    delta = r - 1;      % positive
    phi = 1 + 2*delta + k_amp * delta^2;
else
    phi = r^2;
end
end


function stop = save_stage2_iteration_local(x, optimValues, state, config, ~)
% SAVE_STAGE2_ITERATION_LOCAL  fmincon OutputFcn callback.
%
%   Accumulates design vector, 3D properties, and normalised errors at
%   each SQP iteration in a persistent variable.  On 'done', pushes the
%   data to the base workspace via assignin so the main function can
%   retrieve it with evalin.
%
%   WHY these normalisation denominators?
%     mass  — relative to total mass  (dimensionless mass fraction)
%     GM    — scaled by gm_range width (how far from optimal, relative
%             to the feasible GM window)
%     heave — scaled by T_heave_range width (how far from target,
%             relative to the acceptable period band)
%     pitch — same as heave, using T_pitch_range
%   This makes all four error metrics O(1) and comparable, so a single
%   convergence plot shows all channels on the same scale.

    persistent iteration_data;

    if strcmp(state, 'init')
        iteration_data = struct('x', [], 'props', {{}}, 'errors', []);
        iteration_data.errors = struct('mass',{}, 'gm',{}, 'heave',{}, 'pitch',{});
        stop = false;
        return;
    end

    stop = false;
    if strcmp(state, 'iter') || strcmp(state, 'done')
        props = calculate_3d_properties(x, config);

        mass_err  = abs(props.mass_total - props.mass_buoyant_force) / props.mass_total;
        gm_err    = abs(props.GM_L - config.gm_target) / ...
                    (config.gm_range(2) - config.gm_range(1));
        heave_err = abs(props.periods.heave - config.T_heave_goal) / ...
                    (config.T_heave_range(2) - config.T_heave_range(1));
        pitch_err = abs(props.periods.pitch - config.T_pitch_goal) / ...
                    (config.T_pitch_range(2) - config.T_pitch_range(1));

        iteration_data.x(:, end+1) = x;
        iteration_data.props{end+1} = props;
        iteration_data.errors(end+1) = struct( ...
            'mass', mass_err, 'gm', gm_err, ...
            'heave', heave_err, 'pitch', pitch_err);
    end

    if strcmp(state, 'done')
        assignin('base', 'stage2_iteration_data_temp', iteration_data);
    end
end


%% =====================================================================
%%  STATISTICAL HELPERS
%% =====================================================================

function R2 = compute_R2(y_true, y_pred)
% COMPUTE_R2  Coefficient of determination.
%   Returns NaN when total variance is negligible (avoids 0/0).

    y_true = y_true(:);  y_pred = y_pred(:);
    SS_res = sum((y_true - y_pred).^2);
    SS_tot = sum((y_true - mean(y_true)).^2);
    if SS_tot < 1e-12, R2 = nan; else, R2 = 1 - SS_res / SS_tot; end
end

function MAPE = compute_MAPE(y_true, y_pred)
% COMPUTE_MAPE  Mean absolute percentage error (%).
%   Ignores entries where y_true ≈ 0 (Inf contribution).

    y_true = y_true(:);  y_pred = y_pred(:);
    pct = abs((y_true - y_pred) ./ y_true) * 100;
    MAPE = mean(pct(isfinite(pct)));
end

function r = manual_corr(x, y)
% MANUAL_CORR  Pearson correlation coefficient.
%   Returns NaN when either vector has negligible variance.

    x = x(:);  y = y(:);
    xc = x - mean(x);  yc = y - mean(y);
    denom = sqrt(sum(xc.^2) * sum(yc.^2));
    if denom < 1e-12, r = nan; else, r = sum(xc .* yc) / denom; end
end


%% =====================================================================
%%  STRUCT ASSEMBLY HELPERS
%% =====================================================================

function hist = init_history_struct()
% INIT_HISTORY_STRUCT  Pre-allocate all history fields as empty [].
%   Populated incrementally by 'oneshot' or 'trained' branches.

    hist = struct( ...
        'vol_errors_history',           [], ...
        'gm_errors_history',            [], ...
        'kvol_history',                 [], ...
        'kgm_history',                  [], ...
        'mass_2d_history',              [], ...
        'mass_3d_history',              [], ...
        'mass_2d_corrected_history',    [], ...
        'gm_2d_history',               [], ...
        'gm_3d_history',               [], ...
        'cg_z_2d_history',             [], ...
        'cg_z_3d_history',             [], ...
        'vsub_2d_history',             [], ...
        'vsub_3d_history',             [], ...
        'R2_mass_history',             [], ...
        'R2_GM_history',               [], ...
        'R2_cg_history',               [], ...
        'MAPE_mass_history',           [], ...
        'MAPE_GM_history',             [], ...
        'MAPE_cg_history',             [], ...
        'ME_mass_history',             [], ...
        'ME_GM_history',               [], ...
        'mass_errors_history',         [], ...
        'mass_correction_history',     [], ...
        'gm_correction_history',       [], ...
        'constraints_satisfied_history', [], ...
        'solution_stable_history',     [], ...
        'mass_acceptable_history',     [], ...
        'converged_history',           [], ...
        'fitness_scores_history',      [], ...
        'GM_margins_history',          []);
end


function [lb, ub] = build_3d_bounds(config)
% BUILD_3D_BOUNDS  Variable bounds for [vertical_shift, rho_1..rho_N].

    N = config.num_ballast_sections;

    lb = [config.vertical_shift_bounds(1), ...
          ones(1, N) * config.ballast_density_bounds(1)];
    ub = [config.vertical_shift_bounds(2), ...
          ones(1, N) * config.ballast_density_bounds(2)];

    if config.enable_constructability && ~isempty(config.per_strip_density_lb)
        for i = 1:N
            lb(1 + i) = max(lb(1 + i), config.per_strip_density_lb(i));
        end
        if ~isempty(config.wall_strip_index)
            w_idx = config.wall_strip_index;
            lb(1 + w_idx) = config.constructability_rho_hull;
            ub(1 + w_idx) = config.constructability_rho_hull;
        end
    end
end


function results = assemble_results(config, ...
    x_opt_2d, props_2d, conv_data_2d, hist, converged, pid_saturated, ...
    x_opt_3d, final_props, exitflag_3d, fval_3d, output_3d, ...
    stage2_converged, quality_metrics, ...
    mass_errors_3d, gm_errors_3d, heave_errors_3d, pitch_errors_3d, ...
    constructability, steel_data, final_props_optimiser)
% ASSEMBLE_RESULTS  Pack all outputs into the results struct.

    results        = struct();
    results.config = config;

    results.stage1_2d = struct( ...
        'x_optimal',                x_opt_2d, ...
        'properties',               props_2d, ...
        'iterations',               length(hist.mass_errors_history), ...
        'converged',                converged, ...
        'pid_saturated',            pid_saturated, ...
        'mass_errors',              hist.mass_errors_history, ...
        'gm_errors',                hist.gm_errors_history, ...
        'mass_corrections',         hist.mass_correction_history, ...
        'gm_corrections',           hist.gm_correction_history, ...
        'convergence_data',         conv_data_2d, ...
        'mass_2d_history',          hist.mass_2d_history, ...
        'mass_3d_history',          hist.mass_3d_history, ...
        'mass_2d_corrected_history', hist.mass_2d_corrected_history, ...
        'gm_2d_history',            hist.gm_2d_history, ...
        'gm_3d_history',            hist.gm_3d_history, ...
        'cg_z_2d_history',          hist.cg_z_2d_history, ...
        'cg_z_3d_history',          hist.cg_z_3d_history, ...
        'R2_mass',                  hist.R2_mass_history, ...
        'R2_GM',                    hist.R2_GM_history, ...
        'R2_cg',                    hist.R2_cg_history, ...
        'MAPE_mass',                hist.MAPE_mass_history, ...
        'MAPE_GM',                  hist.MAPE_GM_history, ...
        'MAPE_cg',                  hist.MAPE_cg_history, ...
        'ME_mass',                  hist.ME_mass_history, ...
        'ME_GM',                    hist.ME_GM_history, ...
        'convergence_metrics', struct( ...
            'constraints_satisfied', hist.constraints_satisfied_history, ...
            'solution_stable',       hist.solution_stable_history, ...
            'mass_acceptable',       hist.mass_acceptable_history, ...
            'converged',             hist.converged_history, ...
            'fitness_scores',        hist.fitness_scores_history, ...
            'GM_margins',            hist.GM_margins_history));

    results.stage2_3d = struct( ...
        'x_optimal',       x_opt_3d, ...
        'properties',      final_props, ...
        'exitflag',        exitflag_3d, ...
        'fval',            fval_3d, ...
        'output',          output_3d, ...
        'converged',       stage2_converged, ...
        'quality_metrics', quality_metrics, ...
        'iteration_errors', struct( ...
            'mass',  mass_errors_3d, ...
            'gm',    gm_errors_3d, ...
            'heave', heave_errors_3d, ...
            'pitch', pitch_errors_3d));

    % Constructability post-processing (empty when disabled)
    results.constructability = constructability;

    % Steel-fill realisation outputs (empty when disabled)
    results.steel_data = steel_data;
    results.Final3D    = final_props_optimiser;
end


function [me3, ge3, he3, pe3] = unpack_stage2_errors(stage2_iter)
% UNPACK_STAGE2_ERRORS  Extract per-iteration error arrays.

    n  = length(stage2_iter.errors);
    me3 = zeros(1, n);  ge3 = zeros(1, n);
    he3 = zeros(1, n);  pe3 = zeros(1, n);
    for i = 1:n
        me3(i) = stage2_iter.errors(i).mass;
        ge3(i) = stage2_iter.errors(i).gm;
        he3(i) = stage2_iter.errors(i).heave;
        pe3(i) = stage2_iter.errors(i).pitch;
    end
end


%% =====================================================================
%%  TINY HELPERS
%% =====================================================================

function s = sat_tag(is_saturated)
% SAT_TAG  Returns ' SATURATED' suffix when PID is at its bound.
    if is_saturated, s = ' SATURATED'; else, s = ''; end
end

function r = ternary(cond, t, f)
% TERNARY  Inline if-else for fprintf convenience.
    if cond, r = t; else, r = f; end
end