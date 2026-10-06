function [grids, n_zero_interior] = build_geometry_grid( ...
        config, strip_edges, t_offset_strip, is_solid_strip, ...
        n_z, max_slope_factor)
%BUILD_GEOMETRY_GRID Build z-gridded outer/inner strip areas and moments for UHPC realisation.
% Inputs define hull geometry, per-strip offsets, and solid-strip masks. Outputs hold integration arrays
% and the count of zero-interior nodes; thickness is clipped by configured geometric safeguards.
    persistent cache_key cache_grids

    ms2_model = config.ms2_model;
    if ~isprop(ms2_model, 'filename') && ~isfield(ms2_model, 'filename')
        error('mwecmass:modular_precast:NoParserFilename', ...
              'ms2_model has no .filename property — cannot key cache.');
    end
    parser_key = ms2_model.filename;
    if isempty(parser_key)
        error('mwecmass:modular_precast:EmptyParserFilename', ...
              'ms2_model.filename is empty.  Cannot key cache.');
    end

    tab_z = config.Aw_table_z(:);
    this_key = sprintf( ...
        '%s|sa|nstrip=%d|t=%s|sld=%s|n=%d|s=%.6g|zN=%d|hlo=%.6g|hhi=%.6g', ...
        parser_key, length(strip_edges)-1, ...
        num2str(t_offset_strip(:)', '%.6g,'), ...
        num2str(double(is_solid_strip(:)'), '%d,'), ...
        n_z, max_slope_factor, length(tab_z), ...
        config.hull_z_min, config.hull_z_max);
    if ~isempty(cache_key) && strcmp(cache_key, this_key)
        grids = cache_grids;
        n_zero_interior = sum(grids.A_outer(2:end-1) <= 0);
        return;
    end

    if ~isfield(config, 'boundary_cache') || ...
            ~isstruct(config.boundary_cache) || ...
            ~isfield(config.boundary_cache, 'sources') || ...
            ~isfield(config.boundary_cache, 'u_samples')
        error('mwecmass:modular_precast:BadBoundaryCache', ...
              'config.boundary_cache is missing required fields.');
    end
    n_u = length(config.boundary_cache.u_samples);

    z = linspace(config.hull_z_min, config.hull_z_max, n_z)';
    A_outer   = max(0, interp1(config.Aw_table_z, config.Aw_table,      z, 'linear', 0));
    Iyy_outer = max(0, interp1(config.Aw_table_z, config.I_wp_yy_table, z, 'linear', 0));
    if isfield(config, 'I_wp_xx_table') && ~isempty(config.I_wp_xx_table)
        Ixx_outer = max(0, interp1(config.Aw_table_z, config.I_wp_xx_table, z, 'linear', 0));
    else
        Ixx_outer = Iyy_outer;
    end

    A_inner   = zeros(n_z, 1);
    Iyy_inner = zeros(n_z, 1);
    Ixx_inner = zeros(n_z, 1);

    % Map each z-sample to the strip that contains it
    N_strips = length(t_offset_strip);
    strip_of_z = zeros(n_z, 1);
    for k = 1:n_z
        idx = find(z(k) >= strip_edges(1:end-1) - 1e-9 & ...
                   z(k) <= strip_edges(2:end)   + 1e-9, 1, 'first');
        if isempty(idx)
            if z(k) < strip_edges(1)
                idx = 1;
            else
                idx = N_strips;
            end
        end
        strip_of_z(k) = idx;
    end

    n_zero_interior = 0;
    for k = 1:n_z
        if A_outer(k) <= 1e-10
            if k > 1 && k < n_z
                n_zero_interior = n_zero_interior + 1;
            end
            continue;
        end
        i_strip = strip_of_z(k);
        if is_solid_strip(i_strip)
            % Solid strip → no inner void
            A_inner(k) = 0;  Iyy_inner(k) = 0;  Ixx_inner(k) = 0;
            continue;
        end
        t_k = t_offset_strip(i_strip);
        if ~isfinite(t_k) || t_k <= 0
            A_inner(k) = 0;  Iyy_inner(k) = 0;  Ixx_inner(k) = 0;
            continue;
        end
        [A_inner(k), Iyy_inner(k), Ixx_inner(k)] = ...
            mwecmass.realise.thin_shell.inner_properties_at_z( ...
                config, z(k), t_k, max_slope_factor, ...
                n_u, A_outer(k), Iyy_outer(k), Ixx_outer(k));
    end

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
                   'Ixx_inner', Ixx_inner, ...
                   'strip_of_z', strip_of_z);

    cache_key   = this_key;
    cache_grids = grids;
end
