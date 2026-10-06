classdef PID_Controller < handle
%PID_CONTROLLER Discrete PID controller for Stage-1 correction factors.
% Gains, limits, and state are dimensionless because the input error is normalized.
% Constructor accepts optional OutputLimits, IntegralLimits, and Name parameters.
% update(error,dt) applies proportional/integral/derivative control with saturation and
% anti-windup back-calculation; reset clears state and preserves configured limits.
% See docs/METHODS_ENGINE.md#optim-pid-correction
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

        %% CONSTRUCTOR

        function obj = PID_Controller(Kp, Ki, Kd, varargin)
        % PID_CONTROLLER  Create a PID controller with optional limits.
        % Kp, Ki, Kd are unitless gains. Optional name/value pairs are
        % OutputLimits [min max], IntegralLimits [min max], and Name.

            if nargin < 3
                error('mwecmass:optim:InvalidInput', ...
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

        %% UPDATE — one PID time step

        function output = update(obj, error_normalized, dt)
        % UPDATE  Compute one PID step from a UNITLESS normalised error.
        % error_normalized is dimensionless (typically target/actual−1), dt
        % defaults to 1.0, and output is saturated to OutputLimits. P/I/D use
        % the clamped integral; saturation applies anti-windup back-calculation.

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
                % Back-calculate the integral when saturation would cause windup.
                if output ~= output_raw && obj.Ki ~= 0
                    obj.error_sum = (output - P_term - D_term) / obj.Ki;
                end

                % --- Store state for next iteration ---
                obj.error_prev  = error;
                obj.output_prev = output;

            catch ME
                warning('mwecmass:optim:UpdateFailed', ...
                        'PID update failed for %s: %s', obj.name, ME.message);
                if obj.initialized
                    output = obj.output_prev;
                else
                    output = (obj.min_output + obj.max_output) / 2;
                end
            end
        end

        %% RESET — clear controller state

        function reset(obj)
        % RESET  Clear integrator and derivative state to initial conditions.

            obj.error_sum   = 0;
            obj.error_prev  = 0;
            obj.output_prev = 0;
            obj.initialized = false;
        end

    end  % methods

end  % classdef PID_Controller
