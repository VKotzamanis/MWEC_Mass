classdef WEC_PID_Controller < handle
% WEC_PID_CONTROLLER  Discrete PID controller and penalty shaping for WEC optimisation.
%
%   This handle class provides two services to the optimisation pipeline:
%
%   (A) FEEDBACK CONTROL — discrete-time PID with anti-windup
%       Used by the Stage 1 trained loop to iteratively correct k_vol
%       and k_gm until the 2D surrogate matches the 3D ground truth.
%
%   (B) PENALTY SHAPING — C1-continuous zone penalty for period targets
%       Used by the fmincon objective functions in Stage 1 and Stage 2
%       to penalise natural periods outside an acceptable band.
%
%   WHY a handle class?
%     The PID controller carries state (integral accumulator, previous
%     error) between calls.  A handle class lets WEC_Main_Optimizer hold
%     a reference and call update() repeatedly without copying the
%     object each time.
%
%   WHY unitless error input?
%     The PID operates on normalised (dimensionless) errors:
%       error = (target / actual) − 1
%     This prevents gain saturation from raw SI quantities.  Example:
%       raw mass error = 8 892 kg  →  Kp × 8892 = 7113  (SATURATES)
%       normalised     = 0.099     →  Kp × 0.099 = 0.079 (correct)
%
%   USAGE (PID)
%     pid = WEC_PID_Controller(Kp, Ki, Kd, 'OutputLimits', [lo, hi]);
%     err = (target / actual) - 1.0;          % unitless
%     u   = pid.update(err, dt);
%     pid.reset();
%
%   USAGE (zone penalty)
%     phi = WEC_PID_Controller.zone_penalty(T, T_goal, [T_lo T_hi], k);
%
%   PROPERTIES
%     Kp, Ki, Kd       [-]  PID gains (set by driver, forwarded via config)
%     error_sum        [-]  integral accumulator
%     error_prev       [-]  previous error (for derivative)
%     output_prev      [-]  previous output (for fallback on failure)
%     max_output       [-]  output saturation upper limit
%     min_output       [-]  output saturation lower limit
%     max_integral     [-]  anti-windup integral clamp upper
%     min_integral     [-]  anti-windup integral clamp lower
%     initialized      bool first-call flag
%     name             char controller label (for console diagnostics)
%
%   See also: WEC_Main_Optimizer, WEC_Driver, run_2d_optimizer

    properties
        Kp              % [-]  proportional gain
        Ki              % [-]  integral gain
        Kd              % [-]  derivative gain

        error_sum       % [-]  integral-term accumulator
        error_prev      % [-]  previous error for derivative calculation
        output_prev     % [-]  previous output (fallback if update fails)

        max_output      % [-]  output saturation upper limit
        min_output      % [-]  output saturation lower limit

        max_integral    % [-]  anti-windup integral upper limit
        min_integral    % [-]  anti-windup integral lower limit

        initialized     % bool  first-iteration flag
        name            % char  controller name for diagnostics
    end

    methods

        %% ─────────────────────────────────────────────────────────────
        %%  CONSTRUCTOR
        %% ─────────────────────────────────────────────────────────────

        function obj = WEC_PID_Controller(Kp, Ki, Kd, varargin)
        % WEC_PID_CONTROLLER  Create a PID controller with optional limits.
        %
        %   pid = WEC_PID_CONTROLLER(Kp, Ki, Kd)
        %   pid = WEC_PID_CONTROLLER(Kp, Ki, Kd, 'OutputLimits', [lo hi])
        %   pid = WEC_PID_CONTROLLER(Kp, Ki, Kd, 'IntegralLimits', [lo hi])
        %   pid = WEC_PID_CONTROLLER(Kp, Ki, Kd, 'Name', 'MassCorr')
        %
        %   INPUTS
        %     Kp, Ki, Kd       : PID gains (required, unitless)
        %     'OutputLimits'   : [min max] saturation on PID output
        %     'IntegralLimits' : [min max] anti-windup clamp on ∫e dt
        %     'Name'           : string label for fprintf diagnostics

            if nargin < 3
                error('WEC_PID_Controller:InvalidInput', ...
                      'Requires at least 3 inputs: Kp, Ki, Kd');
            end

            obj.Kp = Kp;
            obj.Ki = Ki;
            obj.Kd = Kd;

            p = inputParser;
            addParameter(p, 'OutputLimits',   [-inf, inf], ...
                         @(x) isnumeric(x) && length(x)==2);
            addParameter(p, 'IntegralLimits', [-inf, inf], ...
                         @(x) isnumeric(x) && length(x)==2);
            addParameter(p, 'Name', 'PID', @ischar);
            parse(p, varargin{:});

            obj.min_output   = p.Results.OutputLimits(1);
            obj.max_output   = p.Results.OutputLimits(2);
            obj.min_integral = p.Results.IntegralLimits(1);
            obj.max_integral = p.Results.IntegralLimits(2);
            obj.name         = p.Results.Name;

            obj.reset();
        end

        %% ─────────────────────────────────────────────────────────────
        %%  UPDATE — one PID time step
        %% ─────────────────────────────────────────────────────────────

        function output = update(obj, error_normalized, dt)
        % UPDATE  Compute one PID step from a UNITLESS normalised error.
        %
        %   output = pid.UPDATE(error_normalized, dt)
        %
        %   INPUTS
        %     error_normalized : [-]  dimensionless, = (target/actual) − 1
        %     dt               : [-]  time step (default 1.0 for discrete)
        %
        %   OUTPUT
        %     output : [-]  saturated control signal
        %
        %   ALGORITHM
        %     1. P = Kp × e
        %     2. I = Ki × ∫e dt   (clamped to IntegralLimits)
        %     3. D = Kd × de/dt
        %     4. u_raw = P + I + D
        %     5. u = clamp(u_raw, OutputLimits)
        %     6. Anti-windup back-calculation: if u ≠ u_raw, recompute
        %        the integral accumulator so that it would produce exactly
        %        the saturated output.  This prevents the integrator from
        %        winding up while the output is pinned at a limit.

            try
                if nargin < 3, dt = 1.0; end

                error = error_normalized;

                % First call: initialise derivative baseline
                if ~obj.initialized
                    obj.error_prev  = error;
                    obj.output_prev = 0;
                    obj.initialized = true;
                end

                % --- P term ---
                P_term = obj.Kp * error;

                % --- I term (with anti-windup clamping) ---
                obj.error_sum = obj.error_sum + error * dt;
                obj.error_sum = max(obj.min_integral, ...
                                min(obj.max_integral, obj.error_sum));
                I_term = obj.Ki * obj.error_sum;

                % --- D term ---
                error_diff = (error - obj.error_prev) / dt;
                D_term     = obj.Kd * error_diff;

                % --- Raw output ---
                output_raw = P_term + I_term + D_term;

                % --- Saturation ---
                output = max(obj.min_output, min(obj.max_output, output_raw));

                % --- Anti-windup back-calculation ---
                %  WHY back-calculate?
                %    When the output saturates, the integrator keeps
                %    accumulating error (windup).  On release, the stored
                %    integral causes a large overshoot.  Back-calculation
                %    resets the integral to the value that would produce
                %    exactly the saturated output, eliminating windup.
                if output ~= output_raw && obj.Ki ~= 0
                    obj.error_sum = (output - P_term - D_term) / obj.Ki;
                end

                % --- Store state for next iteration ---
                obj.error_prev  = error;
                obj.output_prev = output;

            catch ME
                warning('WEC_PID_Controller:UpdateFailed', ...
                        'PID update failed for %s: %s', obj.name, ME.message);
                if obj.initialized
                    output = obj.output_prev;
                else
                    output = (obj.min_output + obj.max_output) / 2;
                end
            end
        end

        %% ─────────────────────────────────────────────────────────────
        %%  RESET — clear controller state
        %% ─────────────────────────────────────────────────────────────

        function reset(obj)
        % RESET  Clear integrator and derivative state to initial conditions.

            obj.error_sum   = 0;
            obj.error_prev  = 0;
            obj.output_prev = 0;
            obj.initialized = false;
        end

        %% ─────────────────────────────────────────────────────────────
        %%  TUNE — update gains at runtime
        %% ─────────────────────────────────────────────────────────────

        function tune(obj, Kp, Ki, Kd)
        % TUNE  Update PID gains during runtime (for adaptive control).
        %
        %   Does NOT reset the integral accumulator — call reset()
        %   separately if a fresh start is needed after re-tuning.

            if nargin ~= 4
                error('WEC_PID_Controller:InvalidTune', 'Requires Kp, Ki, Kd');
            end

            obj.Kp = Kp;
            obj.Ki = Ki;
            obj.Kd = Kd;

            fprintf('  %s: Gains updated to Kp=%.2f, Ki=%.2f, Kd=%.2f\n', ...
                    obj.name, Kp, Ki, Kd);
        end

        %% ─────────────────────────────────────────────────────────────
        %%  STATE CHECKPOINT — serialise / restore
        %% ─────────────────────────────────────────────────────────────

        function state = get_state(obj)
        % GET_STATE  Export controller state for checkpointing / logging.

            state = struct( ...
                'error_sum',   obj.error_sum, ...
                'error_prev',  obj.error_prev, ...
                'output_prev', obj.output_prev, ...
                'gains',       [obj.Kp, obj.Ki, obj.Kd], ...
                'name',        obj.name, ...
                'initialized', obj.initialized);
        end

        function set_state(obj, state)
        % SET_STATE  Restore controller state from a checkpoint struct.

            if ~isstruct(state)
                warning('WEC_PID_Controller:InvalidState', 'State must be struct');
                return;
            end

            obj.error_sum   = state.error_sum;
            obj.error_prev  = state.error_prev;
            obj.output_prev = state.output_prev;
            obj.initialized = true;

            fprintf('  %s: State restored from checkpoint\n', obj.name);
        end

    end  % methods


    %% ═════════════════════════════════════════════════════════════════
    %%  STATIC METHODS — factory constructors + penalty shaping
    %% ═════════════════════════════════════════════════════════════════

    methods (Static)

        %% ─────────────────────────────────────────────────────────────
        %%  FACTORY: mass correction PID
        %% ─────────────────────────────────────────────────────────────

        function pid = mass_correction_pid()
        % MASS_CORRECTION_PID  Factory for the 2D→3D mass-bias PID.
        %
        %   pid = WEC_PID_Controller.mass_correction_pid()
        %
        %   TUNING RATIONALE (weak-start)
        %     Kp = 0.3  : moderate proportional — prevents initial overshoot
        %     Ki = 0.05 : slow integral — eliminates steady-state bias
        %     Kd = 0.02 : light derivative — damps oscillation
        %
        %   WHY weak start?
        %     Early iterations have high uncertainty in the error signal
        %     (the 2D surrogate has never been validated).  Aggressive
        %     gains would overshoot k_vol, causing the next 3D validation
        %     to swing the other way.  Weak gains converge in ~4–6 iters
        %     without oscillation.
        %
        %   NOTE: gains here are defaults for standalone use.  The driver
        %   overrides them via config.pid_mass_gains in the pipeline.

            pid = WEC_PID_Controller(0.3, 0.05, 0.02, ...
                'OutputLimits',   [0.2, 10.0], ...
                'IntegralLimits', [-2, 2], ...
                'Name',           'MassCorrection');
        end

        %% ─────────────────────────────────────────────────────────────
        %%  FACTORY: GM correction PID
        %% ─────────────────────────────────────────────────────────────

        function pid = gm_correction_pid()
        % GM_CORRECTION_PID  Factory for the 2D→3D GM-bias PID.
        %
        %   pid = WEC_PID_Controller.gm_correction_pid()
        %
        %   TUNING RATIONALE (weak-start)
        %     Kp = 0.2  : weaker than mass PID — GM is more sensitive
        %                  to CG shifts than mass is to volume shifts
        %     Ki = 0.03 : slow integral — GM error converges on a
        %                  slower time scale than mass error
        %     Kd = 0.01 : minimal derivative — GM signal is noisier
        %                  (CG is a ratio of two sums, amplifying noise)
        %
        %   NOTE: gains here are defaults for standalone use.  The driver
        %   overrides them via config.pid_gm_gains in the pipeline.

            pid = WEC_PID_Controller(0.2, 0.03, 0.01, ...
                'OutputLimits',   [0.2, 5.0], ...
                'IntegralLimits', [-0.3, 0.3], ...
                'Name',           'GM_Correction');
        end

        %% ─────────────────────────────────────────────────────────────
        %%  ZONE PENALTY — C1-continuous period penalty
        %% ─────────────────────────────────────────────────────────────

        function phi = zone_penalty(T, T_target, T_range, k_amp)
        % ZONE_PENALTY  C1-continuous piecewise-quadratic period penalty.
        %
        %   phi = WEC_PID_Controller.zone_penalty(T, T_target, T_range, k_amp)
        %
        %   Dimensionless penalty for natural-period optimisation.
        %   Zero at T = T_target, rising quadratically inside the
        %   acceptable band, and with amplified curvature (×k_amp)
        %   outside the band.  Value and slope are continuous at the
        %   band boundaries (C1 continuity).
        %
        %   WHY C1-continuous?
        %     fmincon (SQP) approximates the Hessian from finite-difference
        %     gradients.  A slope discontinuity at the band boundary would
        %     create a spurious curvature spike, causing the solver to
        %     oscillate or stall near the boundary.  Matching both value
        %     and slope at r_lo and r_hi eliminates this artefact.
        %
        %   WHY dimensionless r = T/T_target − 1?
        %     The relative error normalises the penalty so that a 10%
        %     deviation from a 7 s target is treated the same as a 10%
        %     deviation from a 2.5 s target.  Without this, heave (longer
        %     periods) would always dominate pitch in the objective.
        %
        %   INPUTS
        %     T        : [s]  actual natural period
        %     T_target : [s]  target natural period
        %     T_range  : [1×2 s]  [T_lo, T_hi] acceptable band
        %     k_amp    : [-]  curvature amplifier outside band (default 3)
        %
        %   OUTPUT
        %     phi : [-]  dimensionless penalty  (0 at target, ~0.08 at boundary)
        %
        %   PIECEWISE DEFINITION  (in r-space, r = T/T_target − 1)
        %
        %     r < r_lo :  phi = r_lo² + 2·r_lo·δ + k_amp·δ²
        %                 where δ = r − r_lo
        %                 (Taylor expansion around r_lo ensures C1 join)
        %
        %     r_lo ≤ r ≤ r_hi :  phi = r²
        %                        (simple quadratic centred on target)
        %
        %     r > r_hi :  phi = r_hi² + 2·r_hi·δ + k_amp·δ²
        %                 where δ = r − r_hi
        %
        %   EXAMPLE
        %     phi = WEC_PID_Controller.zone_penalty(8.0, 7.0, [5 9], 3)
        %     % → 0.0204  (8 s is inside the [5,9] band, small penalty)

            if nargin < 4
                k_amp = 3.0;
            end

            % Degenerate input guard
            if isinf(T) || isnan(T)
                phi = 1.0;          % fixed penalty, not catastrophic
                return;
            end

            % Relative error from target (dimensionless)
            r = T / T_target - 1.0;

            % Band boundaries in r-space
            r_lo = T_range(1) / T_target - 1.0;
            r_hi = T_range(2) / T_target - 1.0;

            if r < r_lo
                % Below band — amplified quadratic, C1-matched at r_lo
                delta = r - r_lo;
                phi   = r_lo^2 + 2*r_lo*delta + k_amp * delta^2;

            elseif r > r_hi
                % Above band — amplified quadratic, C1-matched at r_hi
                delta = r - r_hi;
                phi   = r_hi^2 + 2*r_hi*delta + k_amp * delta^2;

            else
                % Inside band — simple quadratic
                phi = r^2;
            end
        end

    end  % methods (Static)

end  % classdef WEC_PID_Controller
