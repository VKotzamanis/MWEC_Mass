function [x_optimal, final_props, exitflag, convergence_data] = run_2d_optimizer(config, x0, mass_pid)
% RUN_2D_OPTIMIZER - Stage 1: find draft + density distribution
%
% Objective (J_2D) — three dimensionless, unit-weighted terms:
%
%   f = Phi(T_heave) + Phi(T_pitch) + (GM * k_gm / GM_pref - 1)^2
%
% where Phi is a C1-continuous zone penalty (see zone_penalty.m):
%   Inside  [T_target ± delta]:  Phi = (T/T_target - 1)^2
%   Outside range:               curvature amplified by factor k_amp
%
% Constraints (all dimensionless, O(1)):
%   Equality:   mass / (buoyancy * k_mass) - 1 = 0
%   Inequality: 1 - GM*k_gm/GM_min <= 0, monotonic density, density ratio
%
% Corrections k_mass, k_gm appear only here (from PID loop), not in the
% 2D property calculator.
%
% x = [draft, rho_1, ..., rho_N]
%
% INPUTS:
%   config: Configuration struct (includes k_mass, k_gm, gm_pref)
%   x0: [1+N x 1] Initial guess [draft; densities]
%   mass_pid: PID controller (kept for interface consistency)
%
% OUTPUTS:
%   x_optimal, final_props, exitflag, convergence_data

    try
        %% ===== INITIALIZE CONVERGENCE TRACKING =====
        convergence_data = struct();
        convergence_data.fmincon_iterations = [];
        convergence_data.mass_errors = [];
        convergence_data.objective_values = [];
        convergence_data.gm_values = [];
        convergence_data.constraint_violations = [];
        convergence_data.variable_changes = [];
        convergence_data.convergence_rates = [];
        convergence_data.best_feasible = struct('x', x0, 'props', [], 'feasibility', inf);
        
        %% ===== SETUP OPTIMIZATION VARIABLES =====
        % Decision vector: x = [vertical_shift; density_node_1; ...; density_node_N]
        num_densities = config.num_ballast_sections;
        
        % Bounds
        lb = [config.vertical_shift_bounds(1), ...
              ones(1, num_densities) * config.ballast_density_bounds(1)];
        ub = [config.vertical_shift_bounds(2), ...
              ones(1, num_densities) * config.ballast_density_bounds(2)];

        % Per-strip density lower bounds (constructability mode).
        % Enforces t_min by preventing the optimizer from requesting
        % densities that cannot be physically realised.
        if config.enable_constructability && ~isempty(config.per_strip_density_lb)
            for i = 1:num_densities
                lb(1 + i) = max(lb(1 + i), config.per_strip_density_lb(i));
            end
            % Pin wall strip
            if ~isempty(config.wall_strip_index)
                w_idx = config.wall_strip_index;
                lb(1 + w_idx) = config.constructability_rho_hull;
                ub(1 + w_idx) = config.constructability_rho_hull;
            end
        end
        
        % No linear inequality constraints (density ratio handled in nonlinear constraints)
        A_ineq = [];
        b_ineq = [];
        
        fprintf('(2D solve... ');
        
        %% ===== CONFIGURE OPTIMIZER =====
        % Custom output function for convergence tracking
        output_fcn = @(x, optimValues, state) capture_convergence_data(...
            x, optimValues, state, config, convergence_data);
        
        options = optimoptions('fmincon', ...
            'Algorithm', 'sqp', ...
            'Display', 'off', ...
            'MaxFunctionEvaluations', 10000, ...
            'ConstraintTolerance', 1e-3, ...
            'StepTolerance', 1e-4, ...
            'OptimalityTolerance', 1e-5, ...
            'MaxIterations', 1000, ...
            'OutputFcn', output_fcn);
        
        %% ===== RUN OPTIMIZATION =====
        [x_optimal_col, fval_final, exitflag] = fmincon(...
            @objective_func, x0', ...   % Objective
            A_ineq, b_ineq, ...          % Linear inequality
            [], [], ...                  % Linear equality
            lb', ub', ...                % Bounds
            @constraint_func, ...        % Nonlinear constraints
            options);
        
        if isempty(x_optimal_col)
            warning('run_2d_optimizer:EmptyResult', 'fmincon returned empty, using x0');
            x_optimal = x0;
        else
            x_optimal = x_optimal_col';
        end
        
        %% ===== CALCULATE FINAL PROPERTIES =====
        final_props = calculate_2d_properties(x_optimal, config);
        
        %% ===== UPDATE CONVERGENCE DATA =====
        convergence_data.final_objective = fval_final;
        convergence_data.final_exitflag = exitflag;
        
        if final_props.mass_buoyant_force > 1e-6
            convergence_data.final_mass_error = ...
                abs(final_props.mass_discrepancy / final_props.mass_buoyant_force);
        else
            convergence_data.final_mass_error = inf;
        end
        
        convergence_data.total_fmincon_iterations = length(convergence_data.fmincon_iterations);
        
        % Calculate final convergence rate
        if length(convergence_data.mass_errors) > 1
            final_rate = calculate_convergence_rate(convergence_data.mass_errors);
            convergence_data.final_convergence_rate = final_rate;
        end
        
        %% ===== USE BEST FEASIBLE IF FINAL IS INFEASIBLE =====
        if ~isempty(convergence_data.best_feasible.props)
            final_feasibility = calculate_feasibility_score(final_props, config);
            
            if convergence_data.best_feasible.feasibility < final_feasibility
                fprintf('    ⚠ Using best feasible solution from iteration %d\n', ...
                        length(convergence_data.best_feasible.x));
                x_optimal = convergence_data.best_feasible.x;
                final_props = convergence_data.best_feasible.props;
            end
        end
        
        %% ===== REPORT FINAL RESULTS =====
        fprintf('exit=%d, draft=%.3f, GM=%.3f) ', exitflag, x_optimal(1), final_props.GM);
        
    catch ME
        warning('run_2d_optimizer:OptimizationFailed', ...
                'Optimization failed: %s. Returning initial guess.', ME.message);
        
        x_optimal = x0;
        final_props = calculate_2d_properties(x0, config);
        exitflag = -99;
        convergence_data.final_objective = inf;
        convergence_data.final_exitflag = exitflag;
        convergence_data.final_mass_error = inf;
        convergence_data.total_fmincon_iterations = 0;
    end
    
    
    %% ===== NESTED OBJECTIVE FUNCTION =====
    function f = objective_func(x)
        try
            props = calculate_2d_properties(x', config);
            
            % Guard: no waterplane or invalid GM
            if isnan(props.GM) || props.Aw < 1e-6
                f = 1e4;
                return;
            end
            
            % Guard: heave/pitch must be finite
            % (Surge is Inf in free-floating mode — physically correct)
            if isnan(props.periods.heave) || isinf(props.periods.heave) || ...
                    isnan(props.periods.pitch) || isinf(props.periods.pitch)
                f = 1e4;
                return;
            end
            
            %% === J(2D): uniform range-normalised objective ===
            %
            %  f = phi(r_gm) + phi(r_heave) + phi(r_pitch)
            %
            %  r = (actual − target) / half_range for each quantity.
            %  All three terms hit exactly 1.0 at their range boundary.
            %  Matches the Stage 2 (3D) objective normalization.
            
            k_amp = get_field(config, 'zone_k_amp', 3.0);
            
            % Half-ranges
            gm_half    = 0.5 * (config.gm_range(2)      - config.gm_range(1));
            heave_half = 0.5 * (config.T_heave_range(2)  - config.T_heave_range(1));
            pitch_half = 0.5 * (config.T_pitch_range(2)  - config.T_pitch_range(1));
            
            % Normalised errors (r = 0 at target, r = ±1 at boundary)
            gm_pref = get_field(config, 'gm_target', 0.20);
            r_gm    = (props.GM              - gm_pref)             / max(gm_half, 1e-6);
            r_heave = (props.periods.heave   - config.T_heave_goal) / max(heave_half, 1e-6);
            r_pitch = (props.periods.pitch   - config.T_pitch_goal) / max(pitch_half, 1e-6);
            
            f = range_penalty_2d(r_gm, k_amp) ...
              + range_penalty_2d(r_heave, k_amp) ...
              + range_penalty_2d(r_pitch, k_amp);
            
        catch ME
            warning('objective_func:EvalFailed', 'Objective evaluation failed: %s', ME.message);
            f = 1e4;
        end
    end
    
    
    %% ===== NESTED CONSTRAINT FUNCTION =====
    function [c, ceq] = constraint_func(x)
        try
            props = calculate_2d_properties(x', config);
            
            %% INEQUALITY CONSTRAINTS (c ≤ 0), all normalized to O(1)
            densities = x(2:end);
            % When shell is enabled, densities are CORE densities (not bulk).
            % The ratio and monotonic constraints apply to the core only,
            % because the shell density is uniform and poses no
            % manufacturability concern.
            
            % NOTE: k_vol is applied inside calculate_2d_properties.
            % No separate mass correction needed in constraint.
            
            % Identify platform densities (exclude pinned wall strip).
            w_idx = config.wall_strip_index;   % [] when not constructability
            if ~isempty(w_idx)
                platform_densities = densities([1:w_idx-1, w_idx+1:end]);
            else
                platform_densities = densities;
            end
            
            % 1. Density ratio constraint (dimensionless) — platform only
            max_density = max(platform_densities);
            min_density_candidates = platform_densities(platform_densities > 50);
            if ~isempty(min_density_candidates)
                min_density = min(min_density_candidates);
            else
                min_density = 100;
            end
            c_density_ratio = (max_density / min_density) / config.max_density_ratio - 1.0;
            
            % 2. GM floor (dimensionless): 1 - GM/GM_min ≤ 0
            % NOTE: props.GM already includes the k_gm correction on CG_z
            % (applied in calculate_2d_properties §5).  Do NOT multiply
            % by k_gm again — that double-applies the correction.
            c_gm = 1.0 - props.GM / config.gm_min;
            
            % 3. Monotonic density (normalized by rho_max)
            %    Skip wall-platform boundary (wall is pinned, not optimised).
            rho_max = config.ballast_density_bounds(2);
            constrained_pairs = [];
            for i = 1:(length(densities)-1)
                if ~isempty(w_idx) && (i == w_idx || i + 1 == w_idx)
                    continue;   % skip wall-platform boundary
                end
                constrained_pairs(end+1) = i; %#ok<AGROW>
            end
            c_monotonic = zeros(length(constrained_pairs), 1);
            for j = 1:length(constrained_pairs)
                i = constrained_pairs(j);
                c_monotonic(j) = (densities(i+1) - densities(i)) / rho_max;
            end
            
            % Minimum mass constraint (constructability mode).
            if isfield(config, 'm_min_constructability') && ...
                    config.m_min_constructability > 0
                c_mass_min = 1.0 - props.mass_total / config.m_min_constructability;
            else
                c_mass_min = [];
            end

            c = [c_density_ratio; c_gm; c_monotonic; c_mass_min];
            
            %% EQUALITY CONSTRAINT (dimensionless): mass/buoyancy - 1 = 0
            % k_vol is applied inside calculate_2d_properties, so
            % mass and buoyancy are already volume-corrected.
            if props.mass_buoyant_force > 1e-6
                ceq = props.mass_total / props.mass_buoyant_force - 1.0;
            else
                ceq = 1.0;  % No buoyancy → infeasible
            end
            
        catch ME
            warning('constraint_func:EvalFailed', 'Constraint evaluation failed: %s', ME.message);
            n_pairs_fb = max(0, length(densities) - 1);
            if ~isempty(config.wall_strip_index)
                n_pairs_fb = max(0, n_pairs_fb - 1);
            end
            n_extra = 0;
            if isfield(config, 'm_min_constructability') && config.m_min_constructability > 0
                n_extra = 1;
            end
            c = ones(2 + n_pairs_fb + n_extra, 1);   % c_density_ratio + c_gm + c_monotonic + c_mass_min
            ceq = 1.0;
        end
    end
%
end


%% ===== HELPER FUNCTION: CAPTURE CONVERGENCE DATA =====

function stop = capture_convergence_data(x, optimValues, state, config, convergence_data)
% CAPTURE_CONVERGENCE_DATA - OutputFcn callback for fmincon
%
% Tracks iteration history and updates best feasible solution
    
    stop = false;
    
    try
        if strcmp(state, 'iter') || strcmp(state, 'init') || strcmp(state, 'done')
            % Calculate properties at current iterate
            props = calculate_2d_properties(x', config);
            
            % Store iteration data
            convergence_data.fmincon_iterations(end+1) = optimValues.iteration;
            convergence_data.objective_values(end+1) = optimValues.fval;
            convergence_data.gm_values(end+1) = props.GM;
            
            % Mass error (relative)
            if props.mass_buoyant_force > 1e-6
                mass_error = abs(props.mass_discrepancy / props.mass_buoyant_force);
            else
                mass_error = inf;
            end
            convergence_data.mass_errors(end+1) = mass_error;
            
            % Constraint violation (GM only, mass is equality)
            if props.GM < config.gm_min
                violation = config.gm_min - props.GM;
            else
                violation = 0;
            end
            convergence_data.constraint_violations(end+1) = violation;
            
            % Variable change and convergence rate
            if length(convergence_data.fmincon_iterations) > 1
                convergence_data.variable_changes(end+1) = optimValues.stepsize;
                
                if length(convergence_data.mass_errors) >= 3
                    conv_rate = calculate_convergence_rate(...
                        convergence_data.mass_errors(max(1, end-2):end));
                    convergence_data.convergence_rates(end+1) = conv_rate;
                else
                    convergence_data.convergence_rates(end+1) = nan;
                end
            else
                convergence_data.variable_changes(end+1) = nan;
                convergence_data.convergence_rates(end+1) = nan;
            end
            
            % Track best feasible solution
            feasibility_score = calculate_feasibility_score(props, config);
            if feasibility_score < convergence_data.best_feasible.feasibility
                convergence_data.best_feasible.x = x';
                convergence_data.best_feasible.props = props;
                convergence_data.best_feasible.feasibility = feasibility_score;
            end
            
            % Console output suppressed for clean PID output.
            % Data is still tracked in convergence_data for post-analysis.
        end
        
    catch ME
        warning('capture_convergence_data:Failed', 'Data capture failed: %s', ME.message);
    end
end


%% ===== HELPER FUNCTION: CALCULATE CONVERGENCE RATE =====

function rate = calculate_convergence_rate(error_history)
% CALCULATE_CONVERGENCE_RATE - Estimate error reduction rate
%
% Uses log-linear fit of recent error history to estimate convergence rate.
% Handles edge cases: insufficient data, zero variance, negative errors.
%
% INPUT:
%   error_history: [N×1] vector of errors over iterations
%
% OUTPUT:
%   rate: Convergence rate (positive = reducing, NaN = insufficient data)

    try
        % Need at least 2 points for rate calculation
        if length(error_history) < 2
            rate = nan;
            return;
        end
        
        % Use recent history (last 5 points or all available)
        n = length(error_history);
        recent_start = max(1, n - 4);
        
        % Take logarithm (protect against zeros)
        log_errors = log(max(error_history(recent_start:end), 1e-12));
        
        % Check if we have enough variance for fit
        if length(log_errors) >= 3 && var(log_errors) > 1e-12
            x = (1:length(log_errors))';
            p = polyfit(x, log_errors, 1);  % Linear fit in log space
            rate = abs(p(1));  % Slope magnitude
        elseif length(log_errors) == 2
            % Simple difference for 2 points
            rate = abs(log_errors(end) - log_errors(end-1));
        else
            % Insufficient variance or data
            rate = 0;
        end
        
    catch ME
        warning('calculate_convergence_rate:Failed', 'Rate calculation failed: %s', ME.message);
        rate = nan;
    end
end


%% ===== HELPER FUNCTION: CALCULATE FEASIBILITY SCORE =====

function feasibility_score = calculate_feasibility_score(props, config)
% CALCULATE_FEASIBILITY_SCORE - Quantify constraint satisfaction
%
% Uses the same dimensionless terms as J(2D):
%   phi(r_gm) + phi(r_heave) + phi(r_pitch) + mass_error^2
%
% Lower score = more feasible.

    try
        gm_pref = get_field(config, 'gm_target', 0.20);
        k_amp   = get_field(config, 'zone_k_amp', 3.0);
        
        gm_half    = 0.5 * (config.gm_range(2)     - config.gm_range(1));
        heave_half = 0.5 * (config.T_heave_range(2) - config.T_heave_range(1));
        pitch_half = 0.5 * (config.T_pitch_range(2) - config.T_pitch_range(1));
        
        score = 0;
        
        % GM term (same normalization as objective)
        if ~isnan(props.GM)
            r_gm = (props.GM - gm_pref) / max(gm_half, 1e-6);
            score = score + range_penalty_2d(r_gm, k_amp);
        else
            score = score + 1.0;
        end
        
        % Mass balance (relative error squared)
        if props.mass_buoyant_force > 1e-6
            mass_error = abs(props.mass_discrepancy / props.mass_buoyant_force);
            score = score + mass_error^2;
        else
            score = score + 1.0;
        end
        
        % Period terms (same normalization as objective)
        if isinf(props.periods.heave) || isnan(props.periods.heave) || props.periods.heave > 50
            score = score + 1.0;
        else
            r_h = (props.periods.heave - config.T_heave_goal) / max(heave_half, 1e-6);
            score = score + range_penalty_2d(r_h, k_amp);
        end
        
        if isinf(props.periods.pitch) || isnan(props.periods.pitch) || props.periods.pitch > 20
            score = score + 1.0;
        else
            r_p = (props.periods.pitch - config.T_pitch_goal) / max(pitch_half, 1e-6);
            score = score + range_penalty_2d(r_p, k_amp);
        end
        
        feasibility_score = score;
        
    catch ME
        warning('calculate_feasibility_score:Failed', 'Feasibility scoring failed: %s', ME.message);
        feasibility_score = inf;
    end
end


%% ===== HELPER: SAFE CONFIG FIELD ACCESS =====
function val = get_field(s, field_name, default_val)
% GET_FIELD - Return struct field if it exists, otherwise default.
    if isfield(s, field_name)
        val = s.(field_name);
    else
        val = default_val;
    end
end


function phi = range_penalty_2d(r, k_amp)
% RANGE_PENALTY_2D  C1-continuous piecewise-quadratic penalty.
%   Same as range_penalty in WEC_Main_Optimizer.
%   r = 0 → phi = 0;  r = ±1 → phi = 1;  |r| > 1 → amplified.
    if r < -1
        delta = -r - 1;
        phi = 1 + 2*delta + k_amp * delta^2;
    elseif r > 1
        delta = r - 1;
        phi = 1 + 2*delta + k_amp * delta^2;
    else
        phi = r^2;
    end
end