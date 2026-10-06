function WAMIT_Cache_Builder(in, cache_out)
%WAMIT_CACHE_BUILDER Build a WAMIT hydrodynamic cache from per-shift .1 and .3 output.
% Inputs: in (WEC_User_Input, default WEC_User_Input()); cache_out (output file,
% default Input/<hull_tag>_wamit_cache.mat). Discovers all WAMIT output folders for the
% hull, parses .1 and .3 files at each draft, averages radiation damping over [4, 10] s
% periods, and stores the hydrodynamic tables including added mass, damping, and exciting forces.

    if nargin < 1 || isempty(in)
        in = WEC_User_Input();
    end
    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repo_root);
    addpath(fullfile(repo_root, 'src'));

    [~, hull_tag] = fileparts(in.files.ms2_file);
    if nargin < 2 || isempty(cache_out)
        cache_out = fullfile(repo_root, 'Input', [hull_tag '_wamit_cache.mat']);
    end

    % Parameter block: WEC_User_Input has no radiation-damping averaging-band field.
    period_band = [4.0, 10.0]; % [s], positive periods included in the damping average.
    output_root = fullfile(repo_root, 'Input', 'WAMIT', 'WAMIT_OUTPUT');
    [vertical_shifts, output_folders] = discover_drafts(output_root, hull_tag);
    geometry_config = mwecmass.driver.build_config(in, []);
    hydro_table = init_hydro_cache_struct(in, period_band, repo_root);
    omega_ref = [];

    for k = 1:numel(vertical_shifts)
        vertical_shift = vertical_shifts(k); % [m], body-frame shift, signed z positive up.
        output_folder = output_folders{k};
        file_one = fullfile(output_folder, [lower(hull_tag) '.1']);
        file_three = fullfile(output_folder, [lower(hull_tag) '.3']);
        data_one = parse_wamit_1_local(file_one, in.constants.rho_water, in.constants.wamit_L);
        data_three = parse_wamit_3_local(file_three, in.constants.rho_water, ...
            in.constants.g, in.constants.wamit_L);
        if isempty(omega_ref)
            omega_ref = data_one.omega(:);
        end

        [cg_world, ~, ~, displaced_mass, submerged] = mwecmass.hydrostatics.hydrostatic_inputs( ...
            geometry_config.ms2_model, -vertical_shift, geometry_config);
        periods_s = 2 * pi ./ omega_ref; % [s], positive regular wave periods.
        band_mask = periods_s >= period_band(1) & periods_s <= period_band(2);
        radiation_damping_band_avg = mean(data_one.radiation_damping_omega(:, :, band_mask), 3); % [N s/m and rotational equivalents], origin reference.

        hydro_table.drafts(k, 1) = vertical_shift; % [m], body frame, signed z positive up.
        hydro_table.z_cg(k, 1) = cg_world(3); % [m], reference-point z, positive up.
        hydro_table.submerged_volume(k, 1) = submerged.V_sub; % [m^3], positive displaced volume.
        hydro_table.displaced_mass(k, 1) = displaced_mass; % [kg], positive Archimedean displaced mass.
        hydro_table.added_mass_inf{k, 1} = data_one.added_mass_inf; % [kg, kg m, kg m^2], body-origin reference.
        hydro_table.radiation_damping_band_avg{k, 1} = radiation_damping_band_avg; % [N s/m and rotational equivalents], body-origin reference.
        hydro_table.added_mass_omega{k, 1} = data_one.added_mass_omega; % [kg, kg m, kg m^2], body-origin reference.
        hydro_table.radiation_damping_omega{k, 1} = data_one.radiation_damping_omega; % [N s/m and rotational equivalents], body-origin reference.
        hydro_table.exciting_force_omega{k, 1} = squeeze(data_three.exciting_force_omega(:, 1, :)); % [N and N m], body-origin reference.
    end
    hydro_table.omega = omega_ref; % [rad/s], positive circular wave frequency.
    save(cache_out, 'hydro_table', '-v7.3');
end

function [vertical_shifts, folders] = discover_drafts(output_root, hull_tag)
    entries = dir(output_root);
    entries = entries([entries.isdir] & ~startsWith({entries.name}, '.'));
    pattern = ['^' regexptranslate('escape', hull_tag) '_vs_([+-]?\d+\.?\d*)$'];
    vertical_shifts = [];
    folders = {};
    for k = 1:numel(entries)
        token = regexp(entries(k).name, pattern, 'tokens', 'once');
        if ~isempty(token)
            vertical_shifts(end + 1, 1) = str2double(token{1}); %#ok<AGROW>
            folders{end + 1, 1} = fullfile(output_root, entries(k).name); %#ok<AGROW>
        end
    end
    [vertical_shifts, order] = sort(vertical_shifts);
    folders = folders(order);
end

function data = parse_wamit_1_local(pathname, rho_water, ulen)
    fid = fopen(pathname, 'r');
    close_file = onCleanup(@() fclose(fid));
    column_one = []; mode_i = []; mode_j = []; added_mass = []; damping = [];
    while ~feof(fid)
        values = sscanf(fgetl(fid), '%f');
        if numel(values) >= 4
            column_one(end + 1, 1) = values(1); %#ok<AGROW>
            mode_i(end + 1, 1) = round(values(2)); %#ok<AGROW>
            mode_j(end + 1, 1) = round(values(3)); %#ok<AGROW>
            added_mass(end + 1, 1) = values(4); %#ok<AGROW>
            if numel(values) >= 5
                damping(end + 1, 1) = values(5); %#ok<AGROW>
            else
                damping(end + 1, 1) = NaN; %#ok<AGROW>
            end
        end
    end
    exponents = zeros(6, 6);
    for row = 1:6
        for column = 1:6
            exponents(row, column) = 3 + (row > 3) + (column > 3);
        end
    end
    regular_mask = column_one > 0;
    infinity_mask = column_one == 0;
    omega = sort(unique(column_one(regular_mask)), 'ascend');
    data.omega = omega;
    data.added_mass_inf = zeros(6, 6);
    data.added_mass_omega = zeros(6, 6, numel(omega));
    data.radiation_damping_omega = zeros(6, 6, numel(omega));
    for index = find(infinity_mask)'
        row = mode_i(index); column = mode_j(index);
        data.added_mass_inf(row, column) = added_mass(index) * rho_water * ulen^exponents(row, column);
    end
    for index = find(regular_mask)'
        row = mode_i(index); column = mode_j(index);
        frequency_index = find(omega == column_one(index), 1);
        exponent = exponents(row, column);
        data.added_mass_omega(row, column, frequency_index) = added_mass(index) * rho_water * ulen^exponent;
        if ~isnan(damping(index))
            data.radiation_damping_omega(row, column, frequency_index) = ...
                damping(index) * rho_water * omega(frequency_index) * ulen^exponent;
        end
    end
end

function data = parse_wamit_3_local(pathname, rho_water, gravity, ulen)
    fid = fopen(pathname, 'r');
    close_file = onCleanup(@() fclose(fid));
    frequency = []; heading = []; mode = []; real_part = []; imaginary_part = [];
    while ~feof(fid)
        values = sscanf(fgetl(fid), '%f');
        if numel(values) >= 7
            frequency(end + 1, 1) = values(1); %#ok<AGROW>
            heading(end + 1, 1) = values(2); %#ok<AGROW>
            mode(end + 1, 1) = round(values(3)); %#ok<AGROW>
            real_part(end + 1, 1) = values(6); %#ok<AGROW>
            imaginary_part(end + 1, 1) = values(7); %#ok<AGROW>
        end
    end
    omega = sort(unique(frequency(frequency > 0)), 'ascend');
    headings = sort(unique(heading));
    data.exciting_force_omega = zeros(6, numel(headings), numel(omega));
    for index = 1:numel(frequency)
        frequency_index = find(omega == frequency(index), 1);
        heading_index = find(headings == heading(index), 1);
        exponent = 2 + (mode(index) > 3);
        data.exciting_force_omega(mode(index), heading_index, frequency_index) = ...
            (real_part(index) + 1i * imaginary_part(index)) * rho_water * gravity * ulen^exponent;
    end
end

function hydro_table = init_hydro_cache_struct(in, period_band, repo_root)
    hydro_table.drafts = []; % [m], body-frame vertical shift, signed z positive up.
    hydro_table.z_cg = []; % [m], coefficient reference-point z, positive up.
    hydro_table.submerged_volume = []; % [m^3], positive displaced volume.
    hydro_table.displaced_mass = []; % [kg], positive Archimedean displaced mass.
    hydro_table.omega = []; % [rad/s], positive circular wave frequency.
    hydro_table.added_mass_inf = {}; % [kg, kg m, kg m^2], body-origin reference.
    hydro_table.radiation_damping_band_avg = {}; % [N s/m and rotational equivalents], body-origin reference.
    hydro_table.added_mass_omega = {}; % [kg, kg m, kg m^2], body-origin reference.
    hydro_table.radiation_damping_omega = {}; % [N s/m and rotational equivalents], body-origin reference.
    hydro_table.exciting_force_omega = {}; % [N and N m], body-origin reference.
    hydro_table.period_band = period_band; % [s], positive periods used for damping average.
    hydro_table.solver_params = struct('ref_body_center', [0 0 0], ... % [m], body-origin reference.
        'ref_body_length', in.constants.wamit_L, ... % [m], positive WAMIT reference length.
        'wave_diffrac_soln', 1, 'remove_irr_freq', 1, ...
        'depth', in.bem.water_depth); % [m], positive down.
    hydro_table.water_depth = in.bem.water_depth; % [m], positive down.
    hydro_table.ulen = in.constants.wamit_L; % [m], positive WAMIT reference length.
    hydro_table.mesh_Nu = NaN; % [-], not applicable to the higher-order geometry file.
    hydro_table.mesh_Nv = NaN; % [-], not applicable to the higher-order geometry file.
    hydro_table.panel_size = NaN; % [m], not applicable to the higher-order geometry file.
    hydro_table.wp_target_edge = NaN; % [m], not applicable to the higher-order geometry file.
    hydro_table.period_min = in.bem.T_min; % [s], positive analysis-period minimum.
    hydro_table.period_max = in.bem.T_max; % [s], positive analysis-period maximum.
    hydro_table.period_step = in.bem.T_step; % [s], positive analysis-period spacing.
    hydro_table.timestamp = datestr(now); %#ok<DATST,TNOW1>
    hydro_table.ms2_file = in.files.ms2_file;
    ms2_info = dir(fullfile(repo_root, 'Input', in.files.ms2_file));
    if isempty(ms2_info)
        hydro_table.ms2_date = '';
    else
        hydro_table.ms2_date = ms2_info.date;
    end
    hydro_table.solver = 'WAMIT-higher-order';
end
