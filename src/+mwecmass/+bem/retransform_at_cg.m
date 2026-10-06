function config = retransform_at_cg(config, vertical_shift, actual_cg_z)
%RETRANSFORM_AT_CG  Re-do the origin-to-CG hydrodynamic-matrix transform for one draft at the actual CG.
%   Replaces the cached [surge,heave,pitch] coefficients at the nearest
%   draft with a transform to the supplied CG.
%   See docs/METHODS_ENGINE.md#bem-cg-congruence.

    if isempty(config.hydro_drafts); return; end

    % Find the closest draft in the interpolation table
    [~, idx] = min(abs(config.hydro_drafts - vertical_shift));

    % Retrieve A_origin and B_origin from the cache (6×6, at origin)
    cache = config.hydro_cache;
    [~, cache_idx] = min(abs(cache.drafts - vertical_shift));

    A_6x6 = cache.added_mass_inf{cache_idx};
    B_6x6 = cache.radiation_damping_band_avg{cache_idx};

    % Extract 3×3 [surge, heave, pitch]
    idx_3dof = [1, 3, 5];
    A_3x3_origin = A_6x6(idx_3dof, idx_3dof);
    B_3x3_origin = B_6x6(idx_3dof, idx_3dof);

    % Re-transform with ACTUAL CG from non-uniform density
    [A_cg, B_cg] = mwecmass.bem.transform_to_cg( ...
        A_3x3_origin, B_3x3_origin, actual_cg_z);

    % Overwrite the matching entry
    config.added_mass_diagonal(idx, :)        = [A_cg(1,1), A_cg(2,2), A_cg(3,3)];
    config.radiation_damping_diagonal(idx, :) = [B_cg(1,1), B_cg(2,2), B_cg(3,3)];
    config.added_mass_full{idx}               = A_cg;
    config.radiation_damping_full{idx}        = B_cg;
    % Record the new reference height so downstream interpolation does not
    % apply a second transform to these already-transformed values.
    config.hydro_z_cg(idx) = actual_cg_z;
end
