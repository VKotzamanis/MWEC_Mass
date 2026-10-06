function pts = isocurve_devsurf(z_wl, pts_boundary_1, pts_boundary_2)
%ISOCURVE_DEVSURF Compute a cached DevSurf/RuledSurf iso-z contour.
%   Boundary arrays and z_wl are in m; returned points are [N x 3].
%   See docs/METHODS_ENGINE.md#geometry-ruled-isocurve.

    z1 = pts_boundary_1(:, 3);
    z2 = pts_boundary_2(:, 3);
    n_u = length(z1);
    dz = z2 - z1;
    pts = [];

    for i = 1:n_u
        if abs(dz(i)) < 1e-14
            if abs(z1(i) - z_wl) < 1e-10
                v_star = 0.5;
            else
                continue;
            end
        else
            v_star = (z_wl - z1(i)) / dz(i);
        end

        if v_star < -1e-10 || v_star > 1 + 1e-10
            continue;
        end
        v_star = max(0, min(1, v_star));

        pt = (1 - v_star) * pts_boundary_1(i, :) + ...
             v_star * pts_boundary_2(i, :);
        pts = [pts; pt]; %#ok<AGROW>
    end
end
