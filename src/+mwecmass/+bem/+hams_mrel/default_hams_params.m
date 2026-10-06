function params = default_hams_params(config)
%DEFAULT_HAMS_PARAMS  Return default HAMS solver parameters, with config overrides.
% Defaults use deep water, one head-seas heading (0 deg), irregular-frequency removal, and
% zero_inf_limits=1.  Frequencies are period-uniform input (2:0.5:20 s) and omega [rad/s]
% output for the .1/.3 parsers; config.period_min/max/step [s] replace this grid when present.
% Period spacing controls sampling resolution, while panel_size controls irregular-frequency
% ringing.  wave_diffrac_soln=1 requests total Froude-Krylov plus diffraction excitation for RAOs.
% The frequency convention follows docs/HAMS_MREL_ROUTE.md#function-map.

    if nargin < 1, config = struct(); end

    % Depth left NaN for placeholder call to avoid inventing a default.
    if isfield(config, 'water_depth') && ~isempty(config.water_depth)
        params.depth = config.water_depth;   % [m] site water depth from WEC_User_Input; -1 => deep water
    else
        params.depth = NaN;   % [m] unknown until a config with water_depth is attached
    end
    params.zero_inf_limits = 1;  % compute A(0) and A(inf)

    % Frequency grid: period-uniform.  HAMS reads (min, step)
    % as PERIOD in seconds when input_freq_type = 4.
    params.input_freq_type  = 4;   % period [s]
    params.output_freq_type = 3;   % omega [rad/s] (parser convention)
    T_min  = 2.0;    % [s]  shortest period
    T_max  = 20.0;   % [s]  longest period
    T_step = 0.5;    % [s]  period step

    % Apply config overrides if present
    if isfield(config, 'period_min')  && ~isempty(config.period_min)
        T_min  = config.period_min;
    end
    if isfield(config, 'period_max')  && ~isempty(config.period_max)
        T_max  = config.period_max;
    end
    if isfield(config, 'period_step') && ~isempty(config.period_step)
        T_step = config.period_step;
    end

    n_T = round((T_max - T_min) / T_step) + 1;
    assert(T_step > 0,     'period_step must be > 0 (got %g)', T_step);
    assert(T_max > T_min,  'period_max (%g) must be > period_min (%g)', T_max, T_min);
    assert(n_T >= 2,       'period grid must have >=2 entries (got %d)', n_T);

    params.min_frequency = T_min;    % T_min [s]  (label says "Wmin" but it's period when type=4)
    params.freq_step     = T_step;   % ΔT [s]
    params.n_frequencies = -n_T;     % negative → uniform stepping

    % Single heading (head seas) — auto-range mode (n_headings < 0).
    % write_control_file only supports the Minimum_heading+Heading_step
    % pair (verified Fortran format). n_headings > 0 writes a bare value
    % line that mis-parses in HAMS.
    params.n_headings  = -1;    % auto-range: 1 heading
    params.min_heading  = 0.0;  % heading start [deg]
    params.heading_step = 90.0; % heading step  [deg] (irrelevant for 1 heading)

    % Reference body center (rotation centre XR).
    % Always [0,0,0]: HAMS outputs A, B, Fe at the global origin.
    % The post-processor (rebuild_config_hydro) applies
    % the single congruence transform origin → CG.
    % DO NOT override this with CG — doing so causes a double-transform.
    params.ref_body_center = [0, 0, 0];
    params.ref_body_length = 1.0;

    % Solver settings
    %   wave_diffrac_soln = 1: total excitation (FK + diffraction)
    %   wave_diffrac_soln = 2: diffraction potential only (WRONG for RAO)
    params.wave_diffrac_soln = 1;   % total excitation force
    params.remove_irr_freq = 1;     % remove irregular frequencies
    params.n_threads = 4;

    % Minimal field points (required by format but not used)
    params.n_field_points = 1;
    params.field_points = [0, 0, 0];
end
