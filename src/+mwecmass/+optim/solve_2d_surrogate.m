function [x_optimal, final_props, exitflag, convergence_data] = solve_2d_surrogate(config, x0, ~)
%SOLVE_2D_SURROGATE Solve the Stage-1 2-D surrogate with fmincon/SQP.
% x0=[draft;rho_1..rho_N]; config supplies targets, ranges, correction factors, and constraints.
% The objective uses range-normalized GM/heave/pitch penalties; mass balance and GM are constrained.
% Returns x_optimal, final_props, exitflag, and convergence_data. Invalid evaluations use a guard penalty.
    try
        %% Initialize convergence tracking
        convergence_data = struct();
        convergence_data.fmincon_iterations = [];
        convergence_data.mass_errors = [];
        convergence_data.objective_values = [];
        convergence_data.gm_values = [];
        convergence_data.constraint_violations = [];
        convergence_data.variable_changes = [];
        convergence_data.convergence_rates = [];
        convergence_data.best_feasible = struct('x', x0, 'props', [], 'feasibility', inf);

        %% Set up optimization variables
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
            for i = 1:num_densities %#ok<FXUP> -- reused loop index
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

        %% Configure optimizer
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

        %% Run optimization
        [x_optimal_col, fval_final, exitflag] = fmincon(...
            @objective_func, x0', ...   % Objective
            A_ineq, b_ineq, ...          % Linear inequality
            [], [], ...                  % Linear equality
            lb', ub', ...                % Bounds
            @constraint_func, ...        % Nonlinear constraints
            options);

        if isempty(x_optimal_col)
            warning('solve_2d_surrogate:EmptyResult', 'fmincon returned empty, using x0');
            x_optimal = x0;
        else
            x_optimal = x_optimal_col';
        end

        %% Calculate final properties
        final_props = mwecmass.hydrostatics.properties_2d(x_optimal, config);

        %% Update convergence data
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

        %% Use best feasible point when needed
        if ~isempty(convergence_data.best_feasible.props)
            final_feasibility = calculate_feasibility_score(final_props, config);

            if convergence_data.best_feasible.feasibility < final_feasibility
                fprintf('    ⚠ Using best feasible solution from iteration %d\n', ...
                        length(convergence_data.best_feasible.x));
                x_optimal = convergence_data.best_feasible.x;
                final_props = convergence_data.best_feasible.props;
            end
        end

        %% Report final results
        fprintf('exit=%d, draft=%.3f, GM=%.3f) ', exitflag, x_optimal(1), final_props.GM);

    catch ME
        warning('solve_2d_surrogate:OptimizationFailed', ...
                'Optimization failed: %s. Returning initial guess.', ME.message);

        x_optimal = x0;
        final_props = mwecmass.hydrostatics.properties_2d(x0, config);
        exitflag = -99;
        convergence_data.final_objective = inf;
        convergence_data.final_exitflag = exitflag;
        convergence_data.final_mass_error = inf;
        convergence_data.total_fmincon_iterations = 0;
    end


    %% Nested objective function
    function f = objective_func(x)
        try
            props = mwecmass.hydrostatics.properties_2d(x', config);

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

            %% J(2D): uniform range-normalised objective
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


    %% Nested constraint function
    function [c, ceq] = constraint_func(x)
        try
            props = mwecmass.hydrostatics.properties_2d(x', config);

            %% INEQUALITY CONSTRAINTS (c ≤ 0), all normalized to O(1)
            densities = x(2:end);

            % NOTE: k_vol is applied inside properties_2d.
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
            % (applied in the properties_2d objective).  Do NOT multiply
            % by k_gm again — that double-applies the correction.
            % config.gm_min is also the Stage-2 floor for props.GM_L (+optim/run.m), which has
            % no k_gm bias; the same threshold is applied to a biased quantity here and an
            % unbiased one there. Documentation-only note, no value change.
            c_gm = 1.0 - props.GM / config.gm_min;

            c = [c_density_ratio; c_gm];

            % Equality constraint: mass/buoyancy - 1 = 0
            % k_vol is applied inside properties_2d, so
            % mass and buoyancy already include the volume adjustment.
            if props.mass_buoyant_force > 1e-6
                ceq = props.mass_total / props.mass_buoyant_force - 1.0;
            else
                ceq = 1.0;  % No buoyancy → infeasible
            end

        catch ME
            warning('constraint_func:EvalFailed', 'Constraint evaluation failed: %s', ME.message);
            c = ones(2, 1);   % c_density_ratio + c_gm
            ceq = 1.0;
        end
    end
%
end


%% Capture convergence data

function stop = capture_convergence_data(x, optimValues, state, config, convergence_data)
% CAPTURE_CONVERGENCE_DATA - OutputFcn callback for fmincon
%
% Tracks iteration history and updates best feasible solution

    stop = false;

    try
        if strcmp(state, 'iter') || strcmp(state, 'init') || strcmp(state, 'done')
            % Calculate properties at current iterate
            props = mwecmass.hydrostatics.properties_2d(x', config);

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


%% Calculate convergence rate

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


%% Calculate feasibility score

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


%% Read an optional config field
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
%   Same as mwecmass.optim.range_penalty.
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
