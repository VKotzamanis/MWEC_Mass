function [opt_results, x_opt, iteration] = run(config, ~)
%RUN Execute the two-stage ballast optimisation (2-D surrogate then 3-D fmincon/SQP).
% config is the runtime driver configuration; the unused second input is retained for dispatch parity.
% Outputs are opt_results, x_opt=[vertical_shift,rho_1..rho_N], and the Stage-1 iteration count.
% Realisation is a later pipeline stage; this function returns un-realised Stage-2 properties.
try

%% STAGE 1: 2-D DENSITY SEARCH

fprintf('\n');
fprintf('╔══════════════════════════════════════════════════╗\n');
fprintf('║   WEC OPTIMISATION PIPELINE                      ║\n');
fprintf('╚══════════════════════════════════════════════════╝\n');
fprintf(' Timestamp: %s\n\n', char(datetime('now', 'Format', 'dd-MMM-yyyy HH:mm:ss')));

tic_main = tic;

stage1_mode = config.stage1_mode;

fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf('║ STAGE 1: 2-D DENSITY SEARCH — Mode = %-33s║\n', upper(stage1_mode));
fprintf('╚══════════════════════════════════════════════════╝\n\n');

% x = [draft, rho_1, ..., rho_N] — the design vector
x0_2d = [config.initial_vertical_shift, config.initial_densities];

% Surrogate correction factors: unity = no bias correction yet.
% k_vol scales the strip-extrusion submerged volume.
% k_gm  biases the 2D CG position to match 3D.
config.k_vol = config.k_vol_init;
config.k_gm  = config.k_gm_init;

%% PRE-CALIBRATION OF k_vol AND k_gm
% In trained mode, evaluate the initial design in both models and clamp the
% volume and CG ratios to their configured PID bounds before iterating.

if strcmp(stage1_mode, 'trained')
    fprintf(' Pre-calibrating surrogate correction factors...\n');

    try
        % Evaluate initial design in both models
        config_precal      = config;
        config_precal.k_vol = 1.0;
        config_precal.k_gm  = 1.0;

        props_2d_precal = mwecmass.hydrostatics.properties_2d(x0_2d, config_precal);
        props_3d_precal = mwecmass.hydrostatics.properties_3d(x0_2d, config);

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

iteration     = 0;
pid_saturated = false;

switch stage1_mode

    case 'sweep'
        %  HAMS-based Tier-1/Tier-2 landscape sweep.
        %  Tier-1: ranks all cached drafts by |T_heave - T_heave_goal|.
        %  Tier-2: lightweight fmincon on top-K candidates.
        %  Stage 2 warm-starts from the best Tier-2 result.
        [x_opt_2d, props_2d, conv_data_2d, hist] = ...
            mwecmass.optim.stage1_sweep(x0_2d, config, hist);
        converged = true;

    case 'skip'
        %  Bypass Stage 1 entirely.  Stage 2 starts cold from x0_2d
        %  (initial_vertical_shift + initial_densities from config).
        %  Phase A density pre-conditioning still fires if GM < gm_min.
        fprintf(' Stage 1 skipped — passing x0 directly to Stage 2.\n');
        x_opt_2d   = x0_2d;
        props_2d   = mwecmass.hydrostatics.properties_3d(x0_2d, config);
        conv_data_2d = struct('final_objective', NaN, 'final_exitflag', 0);
        converged  = true;

    case 'oneshot'
        [x_opt_2d, props_2d, conv_data_2d, hist, iteration] = ...
            mwecmass.optim.stage1_oneshot(x0_2d, config, hist);
        converged = true;

    case 'trained'
        [x_opt_2d, props_2d, conv_data_2d, hist, ...
         converged, iteration, pid_saturated] = ...
            mwecmass.optim.stage1_trained(x0_2d, config, hist);

    otherwise
        error('WEC:InvalidMode', 'Unknown stage1_mode: %s', stage1_mode);
end

%% STAGE 2: 3-D OPTIMISATION
% fmincon (SQP) from two starts: the Stage-1 result and a bottom-filled design at the Stage-1
% draft. The start whose result has the lower objective is kept and both are logged. Objective and
% constraints are defined by stage2_objective and stage2_constraints.

fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf('║ STAGE 2: 3-D OPTIMISATION                        ║\n');
fprintf('╚══════════════════════════════════════════════════╝\n\n');

% Variable bounds: [vertical_shift, rho_1, ..., rho_N]
[lb_3d, ub_3d] = mwecmass.optim.stage2_bounds(config);

% OutputFcn wrapper captures per-iteration 3D properties.
% The OutputFcn stores iteration data in a temporary base-workspace variable
% because fmincon callbacks cannot return accumulated data to their caller.
output_fcn = @(x, ov, state) save_stage2_iteration_local( ...
    x, ov, state, config, []);

% fmincon options — solver mechanics, not centralised in the driver.
% Solver tolerances and scaling are local numerical settings.
% Algorithm from config (default 'sqp', set in WEC_User_Input.m).
if isfield(config, 'stage2_algorithm') && ~isempty(config.stage2_algorithm)
    stage2_algo = config.stage2_algorithm;
else
    stage2_algo = 'sqp';
end
fprintf('Stage 2 algorithm: %s\n', stage2_algo);

constraint_tol = 1e-8;
options_3d = optimoptions('fmincon', ...
    'Algorithm',              stage2_algo, ...
    'Display',                'iter', ...
    'MaxFunctionEvaluations', 50000, ...
    'MaxIterations',          1000, ...
    'ConstraintTolerance',    constraint_tol, ...
    'OptimalityTolerance',    1e-8, ...
    'StepTolerance',          1e-10, ...
    'FiniteDifferenceStepSize', 1e-6, ...
    'ScaleProblem',           true, ...
    'OutputFcn',              output_fcn);

obj_fun_3d = @(x) mwecmass.optim.stage2_objective(x, config);
con_fun_3d = @(x) mwecmass.optim.stage2_constraints(x, config);

start_x0 = {x_opt_2d, ...
            mwecmass.optim.stage2_bottom_filled_start(x_opt_2d(1), config, lb_3d, ub_3d)};
start_labels = {'Stage-1 result', 'bottom-filled'};
n_starts = numel(start_x0);
stage2_runs = cell(1, n_starts);
for k = 1:n_starts
    fprintf('\n  ── Stage 2 start %d of %d: %s ──\n', k, n_starts, start_labels{k});
    stage2_runs{k} = solve_stage2_start(start_x0{k}, start_labels{k}, config, ...
        lb_3d, ub_3d, options_3d, obj_fun_3d, con_fun_3d);
end

% A start that ends with a constraint violation above the solver's own tolerance does not win
% on objective; when no start ends feasible, the smallest violation wins.
run_fval      = cellfun(@(r) r.fval, stage2_runs);
run_violation = cellfun(@(r) r.output.constrviolation, stage2_runs);
run_feasible  = run_violation <= constraint_tol;
if any(run_feasible)
    candidates = find(run_feasible);
    [~, j_best] = min(run_fval(candidates));
    kept = candidates(j_best);
else
    [~, kept] = min(run_violation);
end

fprintf('\n  Stage 2 starts (constraint tolerance %.0e)\n', constraint_tol);
fprintf('  %-15s %12s %12s %9s %12s %6s %10s %8s\n', 'start', 'f(x0)', 'fval', ...
        'exitflag', 'violation', 'iter', 'mass [kg]', 'GM [m]');
for k = 1:n_starts
    r = stage2_runs{k};
    fprintf('  %-15s %12.5g %12.5g %9d %12.3g %6d %10.1f %8.4f%s\n', r.label, r.f0, r.fval, ...
            r.exitflag, r.output.constrviolation, r.output.iterations, ...
            r.props.mass_total, r.props.GM_L, mwecmass.internal.ternary(k == kept, '   <- kept', ''));
end
fprintf('\n');

x_opt_3d      = stage2_runs{kept}.x;
fval_3d       = stage2_runs{kept}.fval;
exitflag_3d   = stage2_runs{kept}.exitflag;
output_3d     = stage2_runs{kept}.output;
stage2_iter   = stage2_runs{kept}.iter;

% Unpack error arrays for results struct
[mass_errors_3d, gm_errors_3d, heave_errors_3d, pitch_errors_3d] = ...
    unpack_stage2_errors(stage2_iter);

[stage2_converged, quality_metrics] = ...
    mwecmass.optim.check_3d_convergence(x_opt_3d, config, exitflag_3d, output_3d);

%% POST-PROCESSING

fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf('║ OPTIMISATION COMPLETE                             ║\n');
fprintf('╚══════════════════════════════════════════════════╝\n\n');

final_props = mwecmass.hydrostatics.properties_3d(x_opt_3d, config);

if isfield(config, 'profile') && ~isempty(config.profile)
    final_props.cross_section = config.profile;
else
    final_props.cross_section = [];
end

%% STASH OPTIMISER SOLUTION
%  Keep the pure optimiser result separate so realisers can reference it.

final_props_optimiser = final_props;
steel_data            = [];

%% REALISATION DEFERRED TO REALISE
%  Realisation runs in realise: mwecmass.driver.run calls
%  mwecmass.realise.<type>.run(config, x_opt, opt_results) after this function returns, so
%  calling it here too would run it twice. `constructability` is [] here and `steel_data` is
%  already empty; realise overwrites both for the two realising types, and the
%  'preliminary' type leaves them at [].
constructability = [];

%% ASSEMBLE RESULTS STRUCT

opt_results = mwecmass.optim.report_assemble_results(config, ...
    x_opt_2d, props_2d, conv_data_2d, hist, converged, pid_saturated, ...
    x_opt_3d, final_props, exitflag_3d, fval_3d, output_3d, ...
    stage2_converged, quality_metrics, ...
    mass_errors_3d, gm_errors_3d, heave_errors_3d, pitch_errors_3d, ...
    constructability, steel_data, final_props_optimiser);

%% STAGE-2 DESIGN TRAJECTORY
% Persist the accepted SQP design and properties for diagnostics without rerunning.
%
%  x is [1+N x n_iter]: row 1 = vertical_shift [m], rows 2..N+1 = rho [kg/m^3].
%  Sign convention: z_wl = -vertical_shift; draft = |hull_z_min + vs|.
n_it2 = numel(stage2_iter.props);
traj  = struct('x', stage2_iter.x, 'n', n_it2, ...
    'vs',         zeros(1, n_it2), ...   % [m]
    'draft',      zeros(1, n_it2), ...   % [m]
    'mass_total', zeros(1, n_it2), ...   % [kg]
    'mass_buoy',  zeros(1, n_it2), ...   % [kg]
    'GM',         zeros(1, n_it2), ...   % [m]
    'CG_z',       zeros(1, n_it2), ...   % [m] world frame
    'T_heave',    zeros(1, n_it2), ...   % [s]
    'T_pitch',    zeros(1, n_it2), ...   % [s]
    'fval',       zeros(1, n_it2));      % [-] objective
for it_t = 1:n_it2
    p_t = stage2_iter.props{it_t};
    traj.vs(it_t)         = p_t.vertical_shift;
    traj.draft(it_t)      = p_t.draft;
    traj.mass_total(it_t) = p_t.mass_total;
    traj.mass_buoy(it_t)  = p_t.mass_buoyant_force;
    traj.GM(it_t)         = p_t.GM_L;
    traj.CG_z(it_t)       = p_t.CG_total(3);
    traj.T_heave(it_t)    = p_t.periods.heave;
    traj.T_pitch(it_t)    = p_t.periods.pitch;
    traj.fval(it_t)       = mwecmass.optim.stage2_objective(stage2_iter.x(:, it_t), config);
end
opt_results.stage2_3d.trajectory = traj;

% One record per Stage-2 start; kept marks the start whose result is reported above.
stage2_starts = struct('label', start_labels, 'x0', start_x0, ...
    'f0',             cellfun(@(r) r.f0, stage2_runs, 'UniformOutput', false), ...
    'x',              cellfun(@(r) r.x, stage2_runs, 'UniformOutput', false), ...
    'fval',           num2cell(run_fval), ...
    'exitflag',       cellfun(@(r) r.exitflag, stage2_runs, 'UniformOutput', false), ...
    'constrviolation', num2cell(run_violation), ...
    'iterations',     cellfun(@(r) r.output.iterations, stage2_runs, 'UniformOutput', false), ...
    'kept',           num2cell((1:n_starts) == kept));
opt_results.stage2_3d.starts = stage2_starts;

opt_results.optimization_time = toc(tic_main);

fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf('║ optimisation (mwecmass.optim.run) time: %.2f seconds ║\n', opt_results.optimization_time);
fprintf('╚══════════════════════════════════════════════════╝\n\n');

% x_opt_3d is the Stage-2 fmincon output assigned above; this line exposes it as this function's
% second output (x_opt), which mwecmass.driver.run's realise call, realise_fn(config, x_opt,
% opt_results), and the report below read back as x_opt_3d. A rename, not a new computation.
x_opt = x_opt_3d;

% `iteration` (the third output) is the local variable set inside the switch above: 0 at
% initialisation, overwritten by the 'oneshot' and 'trained' branches. It is returned so the
% post-run reporting section can keep its "iteration > 1" gate on the convergence figure.

catch ME
    fprintf('\n  OPTIMISATION FAILED: %s\n', ME.message);
    for i = 1:length(ME.stack)
        fprintf('    %s (line %d)\n', ME.stack(i).name, ME.stack(i).line);
    end
    rethrow(ME);
end

end  % run


%% STAGE 2 SOLVE FROM ONE START

function r = solve_stage2_start(x0, label, config, lb, ub, options, obj_fun, con_fun)
% SOLVE_STAGE2_START  Run the Stage-2 fmincon solve from x0 and collect its result.
%   If x0 violates the GM floor, densities are first settled with the vertical shift pinned
%   (Phase A) and the full solve starts from that point. r holds the start (label, x0, f0), the
%   result (x, fval, exitflag, output, props) and the per-iteration data of the solve (iter).

    r = struct('label', label, 'x0', x0, 'f0', obj_fun(x0));
    props_x0 = mwecmass.hydrostatics.properties_3d(x0, config);
    fprintf('  start x0: vs = %+.4f m, rho = [%s] kg/m^3, f = %.5g, GM = %.4f m\n', ...
            x0(1), sprintf(' %.1f', x0(2:end)), r.f0, props_x0.GM_L);

    %% DENSITY PRE-CONDITIONING
    % If the warm start violates the GM floor, first optimize densities with vertical shift held
    % constant, then use that feasible point for the full solve.
    if props_x0.GM_L < config.gm_min

        fprintf('\n╔──────────────────────────────────────────────────╗\n');
        fprintf('║ STAGE 2 — Phase A: density pre-conditioning       ║\n');
        fprintf('╚──────────────────────────────────────────────────╝\n');
        fprintf('  Warm-start infeasible: GM = %.4f m  (gm_min = %.4f m)\n', ...
                props_x0.GM_L, config.gm_min);
        fprintf('  Pinning vs = %+.4f m — settling density only...\n\n', x0(1));

        lb_a    = lb;  lb_a(1) = x0(1);   % freeze vs at warm-start
        ub_a    = ub;  ub_a(1) = x0(1);   % Route 2 structurally blocked

        opts_a  = optimoptions('fmincon', ...
            'Algorithm',              'sqp', ...
            'Display',                'iter', ...
            'MaxFunctionEvaluations', 500, ...
            'MaxIterations',          50, ...
            'ConstraintTolerance',    1e-4, ...
            'OptimalityTolerance',    1e-4, ...
            'StepTolerance',          1e-6, ...
            'ScaleProblem',           true);

        [x_a, fval_a, ef_a] = fmincon(obj_fun, x0, ...
            [], [], [], [], lb_a, ub_a, con_fun, opts_a);

        props_a = mwecmass.hydrostatics.properties_3d(x_a, config);

        if ef_a > 0 && props_a.GM_L >= config.gm_min
            fprintf('\n  [Phase A] ✓ Feasibility achieved:\n');
            fprintf('    GM      = %.4f m  (target ≥ %.4f m)\n', ...
                    props_a.GM_L, config.gm_min);
            fprintf('    T_heave = %.3f s\n', props_a.periods.heave);
            fprintf('    T_pitch = %.3f s\n', props_a.periods.pitch);
            fprintf('    f       = %.4f  (was %.4f at warm-start)\n', fval_a, r.f0);
            x0 = x_a;   % hand density-settled point to Stage 2 (Phase B)
        else
            fprintf('\n  [Phase A] ✗ Did not achieve GM ≥ %.4f m ', config.gm_min);
            fprintf('(ef=%d, GM=%.4f m).\n', ef_a, props_a.GM_L);
            fprintf('  Possible cause: constructability wall prevents Route 1.\n');
            fprintf('  Stage 2 will proceed from original warm-start.\n');
            %  x0 unchanged — Stage 2 is no worse than without Phase A
        end

        fprintf('\n');
    end
    fprintf('Starting 3D optimisation...\n\n');
    [r.x, r.fval, r.exitflag, r.output] = fmincon(obj_fun, x0, ...
        [], [], [], [], lb, ub, con_fun, options);
    r.props = mwecmass.hydrostatics.properties_3d(r.x, config);

    % Retrieve per-iteration data from base workspace (see assignin note in Stage 2)
    if evalin('base', 'exist(''stage2_iteration_data_temp'', ''var'')')
        r.iter = evalin('base', 'stage2_iteration_data_temp');
        evalin('base', 'clear stage2_iteration_data_temp');
    else
        warning('WEC:NoIterData', 'Stage 2 iteration data not captured');
        r.iter = struct('x', [], 'props', {{}}, 'errors', []);
    end
end


%% STAGE 2 fmincon OUTPUT CALLBACK

function stop = save_stage2_iteration_local(x, ~, state, config, ~)
% SAVE_STAGE2_ITERATION_LOCAL  fmincon OutputFcn callback.
%
%   Accumulates design vector, 3D properties, and normalised errors at
%   each SQP iteration in a persistent variable.  On 'done', pushes the
%   data to the base workspace via assignin so the main function can
%   retrieve it with evalin.
%
%   Errors use each configured target/range so the recorded channels are comparable.

    persistent iteration_data;

    if strcmp(state, 'init')
        iteration_data = struct('x', [], 'props', {{}}, 'errors', []);
        iteration_data.errors = struct('mass',{}, 'gm',{}, 'heave',{}, 'pitch',{});
        % Reset early to avoid reading previous run's stale data if this run errors before first 'iter'.
        assignin('base', 'stage2_iteration_data_temp', iteration_data);
        stop = false;
        return;
    end

    stop = false;
    if strcmp(state, 'iter') || strcmp(state, 'done')
        props = mwecmass.hydrostatics.properties_3d(x, config);

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

        % Flush every iterate to preserve partial results on error, not just at 'done'.
        assignin('base', 'stage2_iteration_data_temp', iteration_data);
    end
end


%% STRUCT ASSEMBLY HELPER

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


function [me3, ge3, he3, pe3] = unpack_stage2_errors(stage2_iter)
% UNPACK_STAGE2_ERRORS  Extract per-iteration error arrays.
% Stays local here (single call site in Sec.3 above, well under the 40-line local-function threshold).

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
