function [grids, n_zero_interior] = build_geometry_grid( ...
        config, t_steel, n_z, max_slope_factor)
%BUILD_GEOMETRY_GRID Build outer/inner section areas and moments on a z grid.
%   [grids,n_zero_interior] = build_geometry_grid(config,t_steel,n_z,max_slope_factor)
%   uses precomputed outer tables and contour offsets for a uniform jacket.
%   Areas are [m^2], second moments [m^4], z [m]; the count reports interior
%   samples whose tabulated outer area is zero (expected boundary degeneracy).

    persistent cache_key cache_grids

    % Hard error if the parser doesn't expose a filename — silently
    % using a constant key would let stale grids leak across runs.
    ms2_model = config.ms2_model;
    if ~isprop(ms2_model, 'filename') && ~isfield(ms2_model, 'filename')
        error('mwecmass:thin_shell:NoParserFilename', ...
              'ms2_model has no .filename property — cannot key cache.');
    end
    parser_key = ms2_model.filename;
    if isempty(parser_key)
        error('mwecmass:thin_shell:EmptyParserFilename', ...
              'ms2_model.filename is empty.  Cannot key cache.');
    end

    % Cache key includes the Aw_table identity (length + endpoints)
    % so a builder-side change to the table invalidates the cache.
    tab_z = config.Aw_table_z(:);
    this_key = sprintf( ...
        '%s|t=%.10g|n=%d|s=%.6g|zN=%d|z0=%.6g|zE=%.6g|hlo=%.6g|hhi=%.6g', ...
        parser_key, t_steel, n_z, max_slope_factor, ...
        length(tab_z), tab_z(1), tab_z(end), ...
        config.hull_z_min, config.hull_z_max);
    if ~isempty(cache_key) && strcmp(cache_key, this_key)
        grids = cache_grids;
        n_zero_interior = sum(grids.A_outer(2:end-1) <= 0);
        return;
    end

    % Validate boundary cache once (defensive — solve() already checks)
    if ~isfield(config, 'boundary_cache') || ...
            ~isstruct(config.boundary_cache) || ...
            ~isfield(config.boundary_cache, 'sources') || ...
            ~isfield(config.boundary_cache, 'u_samples')
        error('mwecmass:thin_shell:BadBoundaryCache', ...
              'config.boundary_cache is missing required fields ("sources"/"u_samples").');
    end
    n_u = length(config.boundary_cache.u_samples);

    %% --- Outer geometry: vectorised lookup from precomputed tables ---
    z = linspace(config.hull_z_min, config.hull_z_max, n_z)';
    A_outer   = max(0, interp1(config.Aw_table_z, config.Aw_table,      z, 'linear', 0));
    Iyy_outer = max(0, interp1(config.Aw_table_z, config.I_wp_yy_table, z, 'linear', 0));
    % Ixx_outer (∫∫ y² dA) needed for the realised hull's roll/yaw inertia.
    % I_wp_xx_table is the outer-section y-second-moment table on this z-grid.
    if isfield(config, 'I_wp_xx_table') && ~isempty(config.I_wp_xx_table)
        Ixx_outer = max(0, interp1(config.Aw_table_z, config.I_wp_xx_table, z, 'linear', 0));
    else
        % Symmetric-hull fallback: Ixx_cross == Iyy_cross when the cross-section
        % is symmetric in x↔y (true for revolution surfaces and for E1).
        Ixx_outer = Iyy_outer;
    end

    A_inner   = zeros(n_z, 1);
    Iyy_inner = zeros(n_z, 1);
    Ixx_inner = zeros(n_z, 1);

    n_zero_interior = 0;
    for k = 1:n_z
        if A_outer(k) <= 1e-10
            if k > 1 && k < n_z
                n_zero_interior = n_zero_interior + 1;
            end
            continue;
        end

        [A_inner(k), Iyy_inner(k), Ixx_inner(k)] = ...
            mwecmass.realise.thin_shell.inner_properties_at_z( ...
                config, z(k), t_steel, max_slope_factor, ...
                n_u, A_outer(k), Iyy_outer(k), Ixx_outer(k));
    end

    % Convex-polygon noise guards
    A_inner   = min(A_inner,   A_outer);
    A_inner   = max(A_inner,   0);
    Iyy_inner = min(Iyy_inner, Iyy_outer);
    Iyy_inner = max(Iyy_inner, 0);
    Ixx_inner = min(Ixx_inner, Ixx_outer);
    Ixx_inner = max(Ixx_inner, 0);

    grids = struct('z', z, ...
                   'A_outer',   A_outer, ...
                   'A_inner',   A_inner, ...
                   'Iyy_outer', Iyy_outer, ...
                   'Iyy_inner', Iyy_inner, ...
                   'Ixx_outer', Ixx_outer, ...
                   'Ixx_inner', Ixx_inner);

    cache_key   = this_key;
    cache_grids = grids;
end
