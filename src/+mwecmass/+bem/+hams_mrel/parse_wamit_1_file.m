function hams_data = parse_wamit_1_file(filepath, ref_body_length, output_freq_type)
%PARSE_WAMIT_1_FILE  Read and dimensionalise HAMS WAMIT-.1 output.
%   filepath names the solver output; ref_body_length L [m] defaults to 1.
%   output_freq_type 3/4 selects omega [rad/s] or period [s] in column 1.
%   Returns A_inf/A_zero [6x6], A/B [6x6xN], omega/periods [N x 1], and raw
%   values using HAMS density and reference-length scaling.
%   See docs/METHODS_ENGINE.md#bem-hams-normalization.
    if nargin < 2 || isempty(ref_body_length)
        ref_body_length = 1.0;
    end
    if nargin < 3 || isempty(output_freq_type)
        output_freq_type = 3;  % DEFAULT: omega in rad/s
    end

    L = ref_body_length;
    hc = mwecmass.bem.hams_mrel.hams_constants();
    rho = hc.RHO_WATER;

    assert(exist(filepath, 'file') == 2, ...
        'WAMIT .1 file not found: %s', filepath);
    assert(ismember(output_freq_type, [3, 4]), ...
        'output_freq_type must be 3 (rad/s) or 4 (period)');

    % -- Read line-by-line (handles mixed 4/5 column format) --
    fid = fopen(filepath, 'r');
    col1_raw = [];
    i_modes = [];
    j_modes = [];
    A_vals = [];
    B_vals = [];

    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line) || isempty(strtrim(line))
            continue;
        end
        vals = sscanf(line, '%f');
        if length(vals) >= 4
            col1_raw(end+1,1) = vals(1); %#ok<AGROW>
            i_modes(end+1,1) = round(vals(2)); %#ok<AGROW>
            j_modes(end+1,1) = round(vals(3)); %#ok<AGROW>
            A_vals(end+1,1) = vals(4); %#ok<AGROW>
            if length(vals) >= 5
                B_vals(end+1,1) = vals(5); %#ok<AGROW>
            else
                B_vals(end+1,1) = NaN; %#ok<AGROW>
            end
        end
    end
    fclose(fid);

    % Store raw values
    hams_data.raw_nondim.col1 = col1_raw;
    hams_data.raw_nondim.i_modes = i_modes;
    hams_data.raw_nondim.j_modes = j_modes;
    hams_data.raw_nondim.A = A_vals;
    hams_data.raw_nondim.B = B_vals;
    hams_data.L = L;
    hams_data.output_freq_type = output_freq_type;

    % Select the length exponent for each DOF pair.
    k_exp = zeros(6, 6);
    for ii = 1:6
        for jj = 1:6
            n_rot = (ii > 3) + (jj > 3);
            k_exp(ii, jj) = 3 + n_rot;
        end
    end

    % Separate special and regular frequencies.
    inf_mask  = col1_raw < 0;    % omega -> inf
    zero_mask = col1_raw == 0;   % omega -> 0
    reg_mask  = col1_raw > 0;

    unique_col1 = unique(col1_raw(reg_mask));

    % Convert the frequency column to omega.
    if output_freq_type == 3
        omega_all = unique_col1;          % already omega
    elseif output_freq_type == 4
        omega_all = 2 * pi ./ unique_col1; % period -> omega
    end

    % Sort omega ascending
    [omega_sorted, sort_idx] = sort(omega_all, 'ascend');
    col1_sorted = unique_col1(sort_idx);
    n_freq = length(omega_sorted);

    % Fail on empty frequency list before indexing omega_sorted(1)/(end).
    if n_freq == 0
        error('mwecmass:hams_mrel:NoRegularFrequencyData', ...
            'parse_wamit_1_file: %s has zero regular-frequency (col1 > 0) rows; the .1 file is truncated or the HAMS run aborted before regular frequencies were written.', ...
            filepath);
    end

    % Initialize dimensional matrices
    hams_data.A_inf  = zeros(6, 6);
    hams_data.A_zero = zeros(6, 6);
    hams_data.A = zeros(6, 6, n_freq);
    hams_data.B = zeros(6, 6, n_freq);

    % Populate A(inf).
    for idx = find(inf_mask)'
        i = i_modes(idx);  j = j_modes(idx);
        if i >= 1 && i <= 6 && j >= 1 && j <= 6
            hams_data.A_inf(i, j) = A_vals(idx) * rho * L^k_exp(i,j);
        end
    end

    % Populate A(0).
    for idx = find(zero_mask)'
        i = i_modes(idx);  j = j_modes(idx);
        if i >= 1 && i <= 6 && j >= 1 && j <= 6
            hams_data.A_zero(i, j) = A_vals(idx) * rho * L^k_exp(i,j);
        end
    end

    % Populate frequency-dependent A and B.
    for idx = find(reg_mask)'
        c1 = col1_raw(idx);
        i = i_modes(idx);  j = j_modes(idx);
        if i < 1 || i > 6 || j < 1 || j > 6; continue; end

        f_idx = find(col1_sorted == c1, 1);
        if isempty(f_idx); continue; end

        omega_k = omega_sorted(f_idx);

        scale_A = rho * L^k_exp(i,j);
        scale_B = rho * omega_k * L^k_exp(i,j);

        hams_data.A(i, j, f_idx) = A_vals(idx) * scale_A;
        if ~isnan(B_vals(idx))
            hams_data.B(i, j, f_idx) = B_vals(idx) * scale_B;
        end
    end

    hams_data.omega = omega_sorted;
    hams_data.periods = 2 * pi ./ omega_sorted;

    fprintf('  Parsed %s: %d frequencies, omega = [%.3f, %.3f] rad/s\n', ...
        filepath, n_freq, omega_sorted(1), omega_sorted(end));
    fprintf('    A33(inf) = %.1f kg,  A55(inf) = %.1f kg*m^2  [dimensional]\n', ...
        hams_data.A_inf(3,3), hams_data.A_inf(5,5));
    if output_freq_type == 3
        fprintf('    (first column = omega [rad/s])\n');
    else
        fprintf('    (first column = period [s])\n');
    end
end
