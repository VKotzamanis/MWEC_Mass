function config = rebuild_config_hydro(config, cache)
%REBUILD_CONFIG_HYDRO  Repopulate the config hydrodynamic-coefficient fields from a hydro_table/hydro_cache.
% Sorts by draft, extracts surge/heave/pitch indices [1,3,5] from each 6x6,
% transforms matrices from origin to CG, and repopulates common BEM fields used
% by downstream interpolation; an empty cache is left unchanged with warning.

    N = length(cache.drafts);
    if N == 0
        warning('mwecmass:bem:EmptyCache', 'Hydro cache is empty.');
        return;
    end

    % Sort by draft (ascending)
    [drafts_sorted, si] = sort(cache.drafts);
    z_cg_sorted  = cache.z_cg(si);
    added_mass_inf_sorted = cache.added_mass_inf(si);
    radiation_damping_band_avg_sorted = cache.radiation_damping_band_avg(si);

    config.hydro_drafts               = drafts_sorted;   % [m] draft (vertical-shift) nodes, signed z up
    config.hydro_z_cg                 = z_cg_sorted;     % [m] vertical CG each node's matrices are referred to
    config.added_mass_diagonal        = zeros(N, 3);     % [A11 kg, A33 kg, A55 kg m^2] per draft node
    config.radiation_damping_diagonal = zeros(N, 3);     % [B11 N s/m, B33 N s/m, B55 N m s/rad], band-averaged
    config.added_mass_full            = cell(N, 1);      % 3x3 [surge, heave, pitch] added mass per node; surge-pitch terms kg m
    config.radiation_damping_full     = cell(N, 1);      % 3x3 band-averaged radiation damping per node
    config.exciting_force_full        = cell(N, 1);      % [3xM] wave exciting force/moment at the node CG, N and N m per period

    idx_3dof = [1, 3, 5];  % surge, heave, pitch

    for i = 1:N
        A_6x6 = added_mass_inf_sorted{i};
        B_6x6 = radiation_damping_band_avg_sorted{i};
        z_cg_i = z_cg_sorted(i);

        % Extract 3×3 submatrix
        A_3x3 = A_6x6(idx_3dof, idx_3dof);
        B_3x3 = B_6x6(idx_3dof, idx_3dof);

        % Transform from origin to CG
        [A_cg, B_cg] = mwecmass.bem.transform_to_cg( ...
            A_3x3, B_3x3, z_cg_i);

        config.added_mass_diagonal(i, :)        = [A_cg(1,1), A_cg(2,2), A_cg(3,3)];
        config.radiation_damping_diagonal(i, :) = [B_cg(1,1), B_cg(2,2), B_cg(3,3)];
        config.added_mass_full{i}               = A_cg;
        config.radiation_damping_full{i}        = B_cg;

        % Fe vector transform: origin → CG
        %   Force vectors (surge, heave) are invariant under reference-point
        %   shift. Only the pitch moment changes:
        %     Fe_pitch_CG = Fe_pitch_origin − z_cg × Fe_surge_origin
        %   In matrix form: Fe_CG = T' * Fe_origin
        %   where T = [1 0 -z_cg; 0 1 0; 0 0 1]  (same T used for A/B).
        %   NOTE: this is the linear (vector) form, NOT the bilinear T'MT
        %   used for A and B.
        %
        %   The raw cache Fe is at the global origin (XR=[0,0,0]).
        %   This CG-shifted version is for downstream 2D/3D EOM consumers.
        %   mwecmass.output.figures.plot_hydrodynamics reads the raw cache
        %   directly (origin-frame EOM) and does NOT use this field.
        % Guard: config.hydro_cache is only populated on the HAMS path.
        % On the WAMIT .1 path (PATH B) rebuild_config_hydro is
        % not called, so this branch never executes there.  The isfield
        % guard makes the function safe even if called defensively.
        if isfield(config, 'hydro_cache') && ...
                ~isempty(config.hydro_cache) && ...
                isfield(config.hydro_cache, 'drafts') && ...
                ~isempty(config.hydro_cache.drafts)
            cache_local = config.hydro_cache;
            [~, cache_i] = min(abs(cache_local.drafts - drafts_sorted(i)));
            if isfield(cache_local, 'exciting_force_omega') && ...
                    ~isempty(cache_local.exciting_force_omega) && ...
                    cache_i <= length(cache_local.exciting_force_omega) && ...
                    ~isempty(cache_local.exciting_force_omega{cache_i})
                Fe_6xM = cache_local.exciting_force_omega{cache_i}; % [6×M] at origin
                Fe_3xM = Fe_6xM(idx_3dof, :);              % [3×M] surge/heave/pitch
                T_fe   = [1, 0, -z_cg_i; 0, 1, 0; 0, 0, 1];
                config.exciting_force_full{i} = T_fe' * Fe_3xM; % [3×M] at CG
            else
                config.exciting_force_full{i} = [];
            end
        else
            config.exciting_force_full{i} = [];
        end
    end

    config.hydro_ready = true;
end
