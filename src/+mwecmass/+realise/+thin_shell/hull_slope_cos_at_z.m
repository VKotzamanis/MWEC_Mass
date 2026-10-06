function cos_alpha = hull_slope_cos_at_z(config, z, A_here, P_here)
%HULL_SLOPE_COS_AT_Z Estimate the local hull-surface slope cosine at elevation z.
%   Uses dA/dz ≈ P*dr/dz from the waterplane-area table and the caller's
%   perimeter P_here [m]. Returns cos(alpha) in [0,1], with edge safeguards.
%   See docs/METHODS_ENGINE.md#thin-shell-slope-offset

    dz = 0.01;

    A_lo = max(0, interp1(config.Aw_table_z, config.Aw_table, z - dz, 'linear', 0));
    A_hi = max(0, interp1(config.Aw_table_z, config.Aw_table, z + dz, 'linear', 0));

    if A_here < 1e-10 || P_here < 1e-10
        cos_alpha = 1.0;
        return;
    end

    % If we're at the top of the hull (Aw vanishes one dz above), assume
    % a near-horizontal cap (the Steiner correction will saturate).
    if A_hi < 1e-10 && A_lo > 1e-10
        cos_alpha = 0.1;
        return;
    end

    % Centred FD where possible; fall back to one-sided at table edges.
    if A_hi > 0 && A_lo > 0
        dA_dz = (A_hi - A_lo) / (2 * dz);
    elseif A_hi > 0
        dA_dz = (A_hi - A_here) / dz;
    else
        dA_dz = (A_here - A_lo) / dz;
    end

    dr_dz     = abs(dA_dz) / max(P_here, 1e-10);
    cos_alpha = 1.0 / sqrt(1.0 + dr_dz^2);
    cos_alpha = max(cos_alpha, 0.01);
end
