function [hams_data, cache] = get_or_run_hydro(vertical_shift, cache, config, hams_dir, hams_exe, tol, options)
%GET_OR_RUN_HYDRO  Return cached hydrodynamics or run HAMS at one vertical shift.
% vertical_shift and tol are scalars in [m]; config comes from build_config, and options
% is passed to mwecmass.bem.hams_mrel.run_at_draft on a miss.  A hit returns the existing
% cache entry without invoking HAMS; a miss appends one entry.  hams_data contains the
% cache schema's 6x6 frequency-dependent and infinite-frequency coefficients, scalar
% z_cg [m], displaced/submerged quantities, omega [rad/s], and status.  cache is updated
% in place and retains its existing field layout.  See docs/METHODS_ENGINE.md#73-hydrodynamic-cache-interpolation-and-cg-transfer;
% HAMS selection is controlled by run_HAMS_MREL.

    if nargin < 6 || isempty(tol); tol = 0.01; end
    if nargin < 7; options = struct(); end

    % Cache lookup
    if ~isempty(cache.drafts)
        [min_dist, idx] = min(abs(cache.drafts - vertical_shift));
        if min_dist < tol
            % Cache hit
            hams_data.added_mass_inf = cache.added_mass_inf{idx};
            hams_data.radiation_damping_band_avg = cache.radiation_damping_band_avg{idx};
            hams_data.added_mass_omega = cache.added_mass_omega{idx};
            hams_data.radiation_damping_omega = cache.radiation_damping_omega{idx};
            hams_data.exciting_force_omega = cache.exciting_force_omega{idx};
            hams_data.omega  = cache.omega;
            hams_data.z_cg   = cache.z_cg(idx);
            hams_data.submerged_volume = cache.submerged_volume(idx);
            hams_data.displaced_mass = cache.displaced_mass(idx);
            hams_data.status = 'cached';
            if isfield(options, 'verbose') && options.verbose
                fprintf('  Cache hit: vs=%+.4f (matched vs=%+.4f, dist=%.4f m)\n', ...
                        vertical_shift, cache.drafts(idx), min_dist);
            end
            return;
        end
    end

    % Cache miss — run HAMS
    hams_data = mwecmass.bem.hams_mrel.run_at_draft( ...
        config, vertical_shift, hams_dir, hams_exe, options);

    % Append to cache
    n = length(cache.drafts) + 1;
    cache.drafts(n,1)  = vertical_shift;
    cache.z_cg(n,1)    = hams_data.z_cg;
    cache.submerged_volume(n,1) = hams_data.submerged_volume;
    cache.displaced_mass(n,1) = hams_data.displaced_mass;
    cache.added_mass_inf{n,1} = hams_data.added_mass_inf;
    cache.radiation_damping_band_avg{n,1} = hams_data.radiation_damping_band_avg;
    cache.added_mass_omega{n,1} = hams_data.added_mass_omega;
    cache.radiation_damping_omega{n,1} = hams_data.radiation_damping_omega;
    cache.exciting_force_omega{n,1} = hams_data.exciting_force_omega;
    if isempty(cache.omega) && ~isempty(hams_data.omega)
        cache.omega = hams_data.omega;
    end
    cache.timestamp = datestr(now); %#ok<DATST,TNOW1> -- cache schema
end
