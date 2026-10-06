function hydro_table = load_hydro_cache(path)
%LOAD_HYDRO_CACHE  Load, validate, and sort a WAMIT hydro_table cache.
% Input path is a .mat file containing hydro_table; output retains its 25-field BEM schema.
% drafts stores signed vertical_shift [m] in the body frame (z up), while water_depth is [m]
% positive down.  The eight per-node fields are reordered together by ascending drafts; scalar
% metadata is unchanged.  Each node must provide 6x6xN omega-dependent matrices, 6x6
% infinite/band-averaged matrices, and 6xN exciting forces.  Missing fields, duplicate drafts,
% or size violations raise bem:wamit:MissingField, bem:wamit:DuplicateDraft, or
% bem:wamit:SizeMismatch, respectively.  See docs/METHODS_ENGINE.md#73-hydrodynamic-cache-interpolation-and-cg-transfer and
% build_config for the cache contract.

    L = load(path, 'hydro_table');
    if ~isfield(L, 'hydro_table')
        error('bem:wamit:MissingField', ...
            'load_hydro_cache: %s does not contain a variable named ''hydro_table''.', path);
    end
    hydro_table = L.hydro_table;

    % Validate required fields.
    required_fields = {'drafts', 'z_cg', 'submerged_volume', 'displaced_mass', 'omega', ...
        'added_mass_inf', 'radiation_damping_band_avg', 'added_mass_omega', ...
        'radiation_damping_omega', 'exciting_force_omega', 'period_band', 'solver_params', ...
        'water_depth', 'ulen', 'mesh_Nu', 'mesh_Nv', 'panel_size', 'wp_target_edge', ...
        'period_min', 'period_max', 'period_step', 'timestamp', 'ms2_file', 'ms2_date', 'solver'};
    missing = required_fields(~isfield(hydro_table, required_fields));
    if ~isempty(missing)
        error('bem:wamit:MissingField', ...
            'load_hydro_cache: %s is missing hydro_table field(s): %s.', path, strjoin(missing, ', '));
    end

    % Sort draft nodes ascending, reordering each parallel per-node field.
    drafts = hydro_table.drafts(:);
    n = numel(drafts);
    [drafts_sorted, order] = sort(drafts);
    hydro_table.drafts = drafts_sorted;
    hydro_table.z_cg   = hydro_table.z_cg(order);
    hydro_table.submerged_volume = hydro_table.submerged_volume(order);
    hydro_table.displaced_mass = hydro_table.displaced_mass(order);
    hydro_table.added_mass_inf = hydro_table.added_mass_inf(order);
    hydro_table.radiation_damping_band_avg = hydro_table.radiation_damping_band_avg(order);
    hydro_table.added_mass_omega = hydro_table.added_mass_omega(order);
    hydro_table.radiation_damping_omega = hydro_table.radiation_damping_omega(order);
    hydro_table.exciting_force_omega = hydro_table.exciting_force_omega(order);

    % Reject duplicate draft nodes after sorting.
    if ~issorted(hydro_table.drafts, 'strictascend')
        error('bem:wamit:DuplicateDraft', ...
            'load_hydro_cache: %s has two draft nodes equal after sorting (interp1 would be ambiguous).', ...
            path);
    end

    % Validate per-node and frequency dimensions.
    n_omega = numel(hydro_table.omega);
    if numel(hydro_table.z_cg) ~= n || numel(hydro_table.submerged_volume) ~= n || ...
            numel(hydro_table.displaced_mass) ~= n || numel(hydro_table.added_mass_omega) ~= n || ...
            numel(hydro_table.radiation_damping_omega) ~= n || numel(hydro_table.added_mass_inf) ~= n || ...
            numel(hydro_table.radiation_damping_band_avg) ~= n || numel(hydro_table.exciting_force_omega) ~= n
        error('bem:wamit:SizeMismatch', ...
            'load_hydro_cache: %s has a per-node field whose length disagrees with numel(drafts)=%d.', ...
            path, n);
    end
    for k = 1:n
        if ~isequal(size(hydro_table.added_mass_omega{k}), [6, 6, n_omega]) || ...
                ~isequal(size(hydro_table.radiation_damping_omega{k}), [6, 6, n_omega])
            error('bem:wamit:SizeMismatch', ...
                'load_hydro_cache: %s node %d has added_mass_omega/radiation_damping_omega not sized [6 6 numel(omega)=%d].', path, k, n_omega);
        end
        if ~isequal(size(hydro_table.added_mass_inf{k}), [6, 6]) || ...
                ~isequal(size(hydro_table.radiation_damping_band_avg{k}), [6, 6])
            error('bem:wamit:SizeMismatch', ...
                'load_hydro_cache: %s node %d has added_mass_inf/radiation_damping_band_avg not sized [6 6].', path, k);
        end
        if ~isequal(size(hydro_table.exciting_force_omega{k}), [6, n_omega])
            error('bem:wamit:SizeMismatch', ...
                'load_hydro_cache: %s node %d has exciting_force_omega not sized [6 numel(omega)=%d].', path, k, n_omega);
        end
    end
end
