function rename_cache_fields(in)
%RENAME_CACHE_FIELDS Rename the shipped cache schema without changing its coefficient arrays.
% Accepts a WEC_User_Input (default WEC_User_Input()) specifying the cache location.
% Loads hydro_table from disk, renames its fields to match the current schema,
% adds depth and reference length, and saves the updated structure.
% Modifies the cache file in place; no return values.

    if nargin < 1 || isempty(in)
        in = WEC_User_Input();
    end
    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    cache_path = fullfile(repo_root, 'Input', in.bem.bem_cache_file);
    loaded = load(cache_path, 'hydro_table');
    hydro_table = loaded.hydro_table;

    legacy_period_min = ['hams' '_T_min'];
    legacy_period_max = ['hams' '_T_max'];
    legacy_period_step = ['hams' '_T_step'];
    legacy_solver_params = ['hams' '_params'];
    old_names = {['V' '_sub'], 'mass', 'A_inf', 'B_avg', 'A', 'B', 'Fe', 'T_band', ...
        legacy_period_min, legacy_period_max, legacy_period_step};
    new_names = {'submerged_volume', 'displaced_mass', 'added_mass_inf', ...
        'radiation_damping_band_avg', 'added_mass_omega', ...
        'radiation_damping_omega', 'exciting_force_omega', 'period_band', ...
        'period_min', 'period_max', 'period_step'};

    for k = 1:numel(old_names)
        value_before = hydro_table.(old_names{k});
        hydro_table.(new_names{k}) = value_before;
        assert(isequal(value_before, hydro_table.(new_names{k})), ...
            'rename_cache_fields:ChangedArray', ...
            'Cache array changed while renaming %s.', old_names{k});
    end
    hydro_table.solver_params = hydro_table.(legacy_solver_params);
    hydro_table.solver_params.depth = in.bem.water_depth; % [m], positive down.
    hydro_table.water_depth = in.bem.water_depth; % [m], positive down.
    hydro_table.ulen = in.constants.wamit_L; % [m], positive WAMIT reference length.
    hydro_table = rmfield(hydro_table, [old_names, {legacy_solver_params}]);
    save(cache_path, 'hydro_table', '-v7.3');
end
