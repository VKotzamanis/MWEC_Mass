function pts = isocurve_bloftsurf(z_wl, n_u, cache_entry)
%ISOCURVE_BLOFTSURF Compute a cached BLoftSurf iso-z contour.
%   z_wl and returned [N x 3] points are in m.

    d = cache_entry;
    degree  = d.degree;
    knots_v = d.knots_v;
    pts = [];

    % ── Detect azimuthal loft ────────────────────────────────
    flat_z_tol = 1e-4;
    i_mid = max(1, round(n_u / 2));
    spread = @(col) max(d.z_sections(:, col)) - min(d.z_sections(:, col));
    is_azimuthal = spread(1) < flat_z_tol && ...
                   spread(i_mid) < flat_z_tol && ...
                   spread(n_u) < flat_z_tol;

    if is_azimuthal
        % ── BRANCH A: azimuthal loft ─────────────────────────
        % z varies with u; use row 1 (any row — all equal) as
        % the z-profile across u.  Root-find u* by sign-change
        % scan, then output the full v-arc at u* from cached pts.
        z_at_u = d.z_sections(1, :)';   % [n_u × 1]
        n_v   = 200;
        v_arc = linspace(0, 1, n_v)';

        for k = 1:n_u - 1
            % Skip if no sign change (crossing) in this interval
            if (z_at_u(k) - z_wl) * (z_at_u(k+1) - z_wl) > 0
                continue;
            end

            % Linear interpolation for fractional position
            dz_k = z_at_u(k+1) - z_at_u(k);
            if abs(dz_k) < 1e-14
                t = 0;
            else
                t = (z_wl - z_at_u(k)) / dz_k;
            end
            t = max(0, min(1, t));

            % Interpolate section_pts [n_sec × 3] at u*
            % d.section_pts is [n_sec × n_u × 3]
            sp = (1 - t) * squeeze(d.section_pts(:, k,   :)) + ...
                      t  * squeeze(d.section_pts(:, k+1, :));

            % Output the full horizontal arc: sweep v in [0,1]
            for iv = 1:n_v
                pv = mwecmass.geometry.MS2Parser.bspline_curve_eval( ...
                         knots_v, sp, degree, v_arc(iv));
                pts = [pts; pv(1,:)]; %#ok<AGROW>
            end
        end
        return;
    end

    % ── BRANCH B: axial loft (original algorithm) ─────────────
    % z varies with v at each u.  For each u, convex-hull check
    % then bisect v* where z(v*) = z_wl.
    for i = 1:n_u
        z_ctrl = d.z_sections(:, i);

        if z_wl < min(z_ctrl) - 1e-10 || z_wl > max(z_ctrl) + 1e-10
            continue;
        end

        v_star = mwecmass.geometry.bspline_param_at_z( ...
                     knots_v, z_ctrl, degree, z_wl);
        if isnan(v_star)
            continue;
        end

        sec_pts_at_u = squeeze(d.section_pts(:, i, :));
        pt = mwecmass.geometry.MS2Parser.bspline_curve_eval( ...
                 knots_v, sec_pts_at_u, degree, v_star);
        pts = [pts; pt(1,:)]; %#ok<AGROW>
    end
end
