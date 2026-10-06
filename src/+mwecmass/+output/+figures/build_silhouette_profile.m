function profile = build_silhouette_profile(config, n_levels)
%BUILD_SILHOUETTE_PROFILE Build a closed hull-silhouette polygon in the midplane (y=0) view.
% n_levels (default 300); z in [m, positive upward]. Returns [2*n_levels x 2]
% polygon of [x, z] vertices, CCW closed. Falls back to config.profile when
% required fields are missing.

    if nargin < 2, n_levels = 300; end

    has_model = isfield(config, 'ms2_model')      && ~isempty(config.ms2_model);
    has_cache = isfield(config, 'boundary_cache') && ~isempty(config.boundary_cache);
    has_zlims = isfield(config, 'hull_z_min')     && isfield(config, 'hull_z_max');

    if ~has_model || ~has_cache || ~has_zlims
        profile = config.profile;
        return;
    end

    z_lo   = config.hull_z_min;
    z_hi   = config.hull_z_max;
    z_vals = linspace(z_lo + 1e-4, z_hi - 1e-4, n_levels)';

    n_u_cache = length(config.boundary_cache.u_samples);
    x_right   = zeros(n_levels, 1);

    for k = 1:n_levels
        pts = mwecmass.geometry.extract_isocurve_at_z( ...
            config.ms2_model, z_vals(k), n_u_cache, config.boundary_cache);
        if ~isempty(pts) && size(pts, 1) >= 2
            x_right(k) = max(pts(:, 1));
        end
    end

    % Levels above or below the hull return an empty (zero-width) contour; drop them.
    valid = x_right > 1e-6;
    x_r   = x_right(valid);
    z_r   = z_vals(valid);

    if length(x_r) < 3
        profile = config.profile;
        return;
    end

    % The polygon closes at the keel, where x_r approaches zero at both ends.
    profile = [x_r, z_r; -flipud(x_r), flipud(z_r)];
end
