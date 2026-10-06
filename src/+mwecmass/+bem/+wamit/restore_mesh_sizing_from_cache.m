function [Nu, Nv, wp_edge] = restore_mesh_sizing_from_cache(hydro_table, in)
%RESTORE_MESH_SIZING_FROM_CACHE  Restore diagnostic mesh sizing from a hydro cache.
% Input hydro_table supplies mesh_Nu/mesh_Nv and optional wp_target_edge [m]; outputs Nu, Nv
% are panel counts (fallback 12,12), and wp_edge is [m] or [] when no cached edge exists.
% Optional in supplies bem.T_min/T_max/T_step [s] for a best-effort stale period-grid warning.
% This supports the WAMIT path in mwecmass.driver.run and preserves the cache metadata
% consumed by downstream diagnostic mesh generation.

    % Old caches fall back to (12, 12); optional in enables period-grid staleness checks.
    if isfield(hydro_table, 'mesh_Nu') && isfield(hydro_table, 'mesh_Nv')
        Nu = hydro_table.mesh_Nu;
        Nv = hydro_table.mesh_Nv;
        fprintf('  Mesh sizing : Nu=%d, Nv=%d (restored from cache)\n', Nu, Nv);
    else
        Nu = 12;
        Nv = 12;
        warning('WEC:OldCacheNoMeshSize', ...
            ['Hydro cache pre-dates the panel_size sizer. ' ...
             'Using fallback (Nu, Nv) = (12, 12) for the diagnostic mesh.\n' ...
             '         Re-running with [R]egenerate will adopt in.geometry.panel_size.']);
    end
    if isfield(hydro_table, 'wp_target_edge') && ~isempty(hydro_table.wp_target_edge)
        wp_edge = hydro_table.wp_target_edge;
        fprintf('  WP target edge : %.3f m (restored from cache)\n', wp_edge);
    else
        wp_edge = [];   % signals caller to leave its own wp_target_edge unchanged
    end

    % Period-grid sanity check (best-effort — old caches lack these fields)
    if nargin >= 2 && isstruct(in)
        cache_keys   = {'period_min', 'period_max', 'period_step'};
        input_keys   = {'T_min',      'T_max',      'T_step'};
        cache_labels = {'T_min', 'T_max', 'T_step'};
        mismatches = {};
        for k = 1:length(cache_keys)
            ck = cache_keys{k};
            ik = input_keys{k};
            if isfield(hydro_table, ck) && ~isempty(hydro_table.(ck)) && ...
                    isfield(in.bem, ik) && ~isempty(in.bem.(ik))
                if abs(hydro_table.(ck) - in.bem.(ik)) > 1e-6
                    mismatches{end+1} = sprintf( ...
                        '%s: cache=%.3f s, input=%.3f s', ...
                        cache_labels{k}, hydro_table.(ck), in.bem.(ik)); %#ok<AGROW>
                end
            end
        end
        if ~isempty(mismatches)
            warning('WEC:CachePeriodGridMismatch', ...
                ['Loaded hydro cache uses a DIFFERENT period grid than ' ...
                 'in.bem.T_min/T_max/T_step:\n         %s\n' ...
                 '         Cached hydrodynamic coefficients are sampled on the ' ...
                 'OLD grid; downstream consumers will interpolate but cannot ' ...
                 'recover data outside the cached omega range.\n' ...
                 '         Choose [R] at the next run to regenerate on the ' ...
                 'current grid.'], strjoin(mismatches, '; '));
        end
    end
end
