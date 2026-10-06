function fe_data = parse_wamit_3_file(filepath, ref_body_length, output_freq_type)
%PARSE_WAMIT_3_FILE  Read HAMS WAMIT-.3 output and dimensionalize exciting force.
% filepath names the .3 file; optional ref_body_length L is [m] (default 1), and
% output_freq_type is 3 for omega input [rad/s] or 4 for period input [s].  The parser
% auto-detects six- or seven-column records and uses Fe = Fe_bar*rho*g*L^m, with m=2 for
% translation modes 1:3 and m=3 for rotation modes 4:6.  Outputs are dimensional SI:
% omega [N x 1, rad/s], periods [N x 1, s], Fe/Fe_mag/Fe_phase [6 x N, N or N*m/deg],
% and headings [H x 1, deg], plus metadata.  See docs/HAMS_MREL_ROUTE.md#function-map.

    if nargin < 2 || isempty(ref_body_length)
        ref_body_length = 1.0;
    end
    if nargin < 3 || isempty(output_freq_type)
        output_freq_type = 3;
    end

    L = ref_body_length;
    hc = mwecmass.bem.hams_mrel.hams_constants();
    rho = hc.RHO_WATER;
    grav = hc.G;

    assert(exist(filepath, 'file') == 2, ...
        'WAMIT .3 file not found: %s', filepath);

    % Read numeric records.
    fid = fopen(filepath, 'r');
    raw_data = [];
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line) || isempty(strtrim(line))
            continue;
        end
        vals = sscanf(line, '%f');
        if length(vals) >= 6
            raw_data(end+1, 1:length(vals)) = vals(:)'; %#ok<AGROW>
        end
    end
    fclose(fid);

    if isempty(raw_data)
        warning('mwecmass:hams_mrel:Empty3File', 'No data in .3 file: %s', filepath);
        fe_data = struct('omega', [], 'Fe', []);
        return;
    end

    n_cols = size(raw_data, 2);

    % Detect the available output-column format.
    col2_vals = unique(raw_data(:, 2));
    if n_cols >= 7 && any(col2_vals > 6 | col2_vals < 1 | ...
            mod(col2_vals, 1) ~= 0)
        heading_col = 2; mode_col = 3;
        mag_col = 4; phase_col = 5; re_col = 6; im_col = 7;
        fmt = 7;
    else
        heading_col = 0; mode_col = 2;
        mag_col = 3; phase_col = 4; re_col = 5; im_col = 6;
        fmt = 6;
    end

    all_col1  = raw_data(:, 1);
    all_modes = round(raw_data(:, mode_col));
    if heading_col > 0
        all_headings = raw_data(:, heading_col);
    else
        all_headings = zeros(size(all_col1));
    end

    reg_mask = all_col1 > 0;
    unique_col1 = unique(all_col1(reg_mask));

    if output_freq_type == 3
        omega_all = unique_col1;
    elseif output_freq_type == 4
        omega_all = 2 * pi ./ unique_col1;
    end

    [omega_sorted, sort_idx] = sort(omega_all, 'ascend');
    col1_sorted = unique_col1(sort_idx);
    n_freq = length(omega_sorted);
    unique_headings = unique(all_headings(reg_mask));

    % Fail on empty frequency list before indexing unique_headings(1).
    if n_freq == 0
        error('mwecmass:hams_mrel:NoRegularFrequencyData', ...
            'parse_wamit_3_file: %s has zero regular-frequency (col1 > 0) rows; the .3 file is truncated or the HAMS run aborted before regular frequencies were written.', ...
            filepath);
    end

    m_exp = [2, 2, 2, 3, 3, 3];

    Fe_complex = zeros(6, n_freq);
    Fe_mag     = zeros(6, n_freq);
    Fe_phase   = zeros(6, n_freq);

    target_heading = unique_headings(1);

    for idx = 1:size(raw_data, 1)
        c1 = all_col1(idx);
        if c1 <= 0; continue; end
        heading = all_headings(idx);
        if abs(heading - target_heading) > 1e-6; continue; end
        i = all_modes(idx);
        if i < 1 || i > 6; continue; end
        f_idx = find(col1_sorted == c1, 1);
        if isempty(f_idx); continue; end

        scale = rho * grav * L^m_exp(i);
        Fe_mag(i, f_idx)   = raw_data(idx, mag_col) * scale;
        Fe_phase(i, f_idx) = raw_data(idx, phase_col);
        re_val = raw_data(idx, re_col) * scale;
        im_val = raw_data(idx, im_col) * scale;
        Fe_complex(i, f_idx) = re_val + 1i * im_val;
    end

    fe_data.omega    = omega_sorted;
    fe_data.periods  = 2 * pi ./ omega_sorted;
    fe_data.Fe       = Fe_complex;
    fe_data.Fe_mag   = Fe_mag;
    fe_data.Fe_phase = Fe_phase;
    fe_data.headings = unique_headings;
    fe_data.L        = L;
    fe_data.output_freq_type = output_freq_type;
    fe_data.format_detected = fmt;

    fprintf('  Parsed %s: %d frequencies, %d headings (format: %d-col)\n', ...
        filepath, n_freq, length(unique_headings), fmt);
    nonzero3 = Fe_mag(3, Fe_mag(3,:) > 0);
    if ~isempty(nonzero3)
        fprintf('    |Fe3| range: [%.1f, %.1f] N  [dimensional]\n', ...
            min(nonzero3), max(nonzero3));
    end
end
