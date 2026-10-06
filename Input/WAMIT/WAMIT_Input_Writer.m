function WAMIT_Input_Writer(in)
%WAMIT_INPUT_WRITER Write WAMIT .pot, .frc, and .cfg decks for each vertical shift.
% Input: in (WEC_User_Input, default WEC_User_Input()); reads periods in.bem.T_min:T_step:T_max,
% water depth, and geometric model. Discovers vertical shifts from existing WAMIT_OUTPUT
% directories. For each shift, computes displaced mass, center of gravity, and inertia,
% then writes three WAMIT input deck files with wave period, frequency-domain modes, and
% body mass/inertia. No return values; modifies Input/WAMIT/WAMIT_INPUT/ on disk.

    if nargin < 1 || isempty(in)
        in = WEC_User_Input();
    end

    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(repo_root);
    addpath(fullfile(repo_root, 'src'));

    % Parameter block: deck-format controls have no WEC_User_Input field.
    controls.append_infinite_frequency = true;
    controls.headings_deg = 0.0;             % [deg] wave heading, positive clockwise.
    controls.irad = 1;
    controls.idiff = 1;
    controls.modes = [1 1 1 1 1 1];
    controls.ioptn = [1 2 2 2 0 0 0 0 0];
    controls.cfg = struct('panel_size', 0.20, ... % [m] higher-order panel edge target.
                          'ksplin', 4, 'iquadi', 5, 'iquado', 4, ...
                          'irr', 3, 'isolve', 1, 'ncpu', 16, ...
                          'ramgbmax', 32.00, 'userid_path', 'C:\WAMITv7', ...
                          'xtrim_round', 3);

    [~, hull_tag] = fileparts(in.files.ms2_file);
    input_dir = fullfile(repo_root, 'Input', 'WAMIT', 'WAMIT_INPUT');
    output_dir = fullfile(repo_root, 'Input', 'WAMIT', 'WAMIT_OUTPUT');
    vertical_shifts = discover_vertical_shifts(output_dir, hull_tag);
    periods_s = in.bem.T_min:in.bem.T_step:in.bem.T_max; % [s], positive analysis periods.
    if controls.append_infinite_frequency
        periods_s = [0.0, periods_s]; % [s], WAMIT's infinite-frequency marker.
    end

    geometry_config = mwecmass.driver.build_config(in, []);
    total_volume = geometry_config.total_wec_volume; % [m^3], positive enclosed hull volume.
    gdf_file = [hull_tag '.gdf'];

    for k = 1:numel(vertical_shifts)
        vertical_shift = vertical_shifts(k); % [m], body-frame shift, signed z positive up.
        [cg_world, ~, ~, displaced_mass, ~] = mwecmass.hydrostatics.hydrostatic_inputs( ...
            geometry_config.ms2_model, -vertical_shift, geometry_config);
        cg_body = cg_world - [0, 0, vertical_shift]; % [m], world-to-body; z positive up.
        body_density = displaced_mass / total_volume; % [kg/m^3], positive uniform density.
        cz_body = cg_body(3); % [m], body frame, z positive up.
        inertia_cg = diag([ ...
            max(0, body_density * (geometry_config.hull_int_y2 + geometry_config.hull_int_z2 - total_volume * cz_body^2)), ... % [kg m^2]
            max(0, body_density * (geometry_config.hull_int_x2 + geometry_config.hull_int_z2 - total_volume * cz_body^2)), ... % [kg m^2]
            max(0, body_density * (geometry_config.hull_int_x2 + geometry_config.hull_int_y2))]); % [kg m^2]
        external_mass = mwecmass.hydrostatics.build_mass_matrix( ...
            displaced_mass, cg_body, inertia_cg); % [kg, kg m, kg m^2], body-origin convention.

        deck_name = sprintf('%s_vs%+0.4f', hull_tag, vertical_shift);
        write_wamit_pot(fullfile(input_dir, [deck_name '.pot']), gdf_file, periods_s, ...
            controls.headings_deg, controls.modes, in.bem.water_depth, ...
            controls.irad, controls.idiff, hull_tag, vertical_shift);
        write_wamit_frc(fullfile(input_dir, [deck_name '.frc']), controls.ioptn, ...
            in.constants.rho_water, cg_body, external_mass, hull_tag, vertical_shift);
        write_wamit_cfg(fullfile(input_dir, [deck_name '.cfg']), vertical_shift, ...
            controls.cfg, hull_tag);
    end
end

function vertical_shifts = discover_vertical_shifts(output_dir, hull_tag)
    entries = dir(output_dir);
    entries = entries([entries.isdir] & ~startsWith({entries.name}, '.'));
    pattern = ['^' regexptranslate('escape', hull_tag) '_vs_([+-]?\d+\.?\d*)$'];
    vertical_shifts = [];
    for k = 1:numel(entries)
        token = regexp(entries(k).name, pattern, 'tokens', 'once');
        if ~isempty(token)
            vertical_shifts(end + 1, 1) = str2double(token{1}); %#ok<AGROW>
        end
    end
    vertical_shifts = sort(vertical_shifts);
end

function write_wamit_pot(pathname, gdf_file, periods_s, headings_deg, modes, water_depth, irad, idiff, hull_tag, vertical_shift)
    fid = fopen(pathname, 'w');
    close_file = onCleanup(@() fclose(fid));
    fprintf(fid, '%s hull, vertical_shift = %+0.4f m\n', hull_tag, vertical_shift);
    fprintf(fid, '%-14g    HBOT     (water depth; <= 0 => infinite)\n', water_depth);
    fprintf(fid, '%-3d %-3d        IRAD IDIFF\n', irad, idiff);
    fprintf(fid, '%-5d            NPER  (periods in s)\n', numel(periods_s));
    write_value_block(fid, periods_s, 8);
    fprintf(fid, '%-5d            NBETA\n', numel(headings_deg));
    write_value_block(fid, headings_deg, 8);
    fprintf(fid, '%-5d            NBODY\n', 1);
    fprintf(fid, '%s\n', gdf_file);
    fprintf(fid, '0.0  0.0  0.0  0.0    XBODY YBODY ZBODY ALPHA\n');
    fprintf(fid, '%d %d %d %d %d %d    IMODE(1..6)\n', modes);
end

function write_wamit_frc(pathname, ioptn, rho_water, cg_body, external_mass, hull_tag, vertical_shift)
    fid = fopen(pathname, 'w');
    close_file = onCleanup(@() fclose(fid));
    fprintf(fid, '%s hull, vs=%+0.4f m | EXMASS about body origin\n', hull_tag, vertical_shift);
    fprintf(fid, '%d %d %d %d %d %d %d %d %d    IOPTN(1..9)\n', ioptn);
    fprintf(fid, '%-15.6e RHO  (kg/m^3)\n', rho_water);
    fprintf(fid, '%14.6e %14.6e %14.6e   XCG YCG ZCG  (body coords)\n', cg_body);
    fprintf(fid, '%-5d            IMASS\n', 1);
    write_6x6_block(fid, external_mass);
    fprintf(fid, '\n%-5d            IDAMP\n', 1);
    write_6x6_block(fid, zeros(6));
    fprintf(fid, '\n%-5d            ISTIFF\n', 1);
    write_6x6_block(fid, zeros(6));
    fprintf(fid, '\n%-5d            NBETAH\n', 0);
    fprintf(fid, '%-5d            NFIELD\n', 0);
end

function write_wamit_cfg(pathname, vertical_shift, cfg, hull_tag)
    fid = fopen(pathname, 'w');
    close_file = onCleanup(@() fclose(fid));
    xtrim = round(vertical_shift, cfg.xtrim_round); % [m], body heave, positive z up.
    fprintf(fid, '! MODAL ANALYSIS  (hull %s, vertical_shift = %+0.4f m)\n', hull_tag, vertical_shift);
    fprintf(fid, ' IPLTDAT=15\n ITRIMWL=1\n');
    fprintf(fid, ' XTRIM= %.*f 0.0 0.0\n', cfg.xtrim_round, xtrim);
    fprintf(fid, ' ILOWHI=1\n ISOLVE=%d\n ISOR = 0\n', cfg.isolve);
    fprintf(fid, ' PANEL_SIZE=%.2f\n KSPLIN=%d\n IQUADI=%d\n IQUADO=%d\n', ...
        cfg.panel_size, cfg.ksplin, cfg.iquadi, cfg.iquado);
    fprintf(fid, ' IRR=%d\n IPNLBPT=0\n IPERIN=1\n IPEROUT=2\n', cfg.irr);
    fprintf(fid, ' MONITR=0\n NUMHDR=1\n IALTFRC=2\n TOLGAPWL= 0.001\n VMAXOPT9=-1\n');
    fprintf(fid, ' NCPU=%d\n RAMGBMAX=%.2f\n USERID_PATH=%s\n', ...
        cfg.ncpu, cfg.ramgbmax, cfg.userid_path);
end

function write_value_block(fid, values, per_line)
    for start_index = 1:per_line:numel(values)
        stop_index = min(start_index + per_line - 1, numel(values));
        fprintf(fid, '  %.6f', values(start_index:stop_index));
        fprintf(fid, '\n');
    end
end

function write_6x6_block(fid, matrix)
    for row = 1:6
        fprintf(fid, '  %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e\n', matrix(row, :));
    end
end
