function [r_min, r_pts] = compute_rmin_at_z(parser, z_target, n_samples, cache)
%COMPUTE_RMIN_AT_Z Compute minimum centroid-to-boundary radius.
%   z_target and r_min are in m; r_pts is an [N x 2] x-y contour.

    if nargin < 3 || isempty(n_samples), n_samples = 100; end

    r_pts = [];

    % Build or reuse boundary cache
    if nargin >= 4 && ~isempty(cache)
        b_cache = cache;
    else
        b_cache = mwecmass.geometry.precompute_boundary_cache( ...
                      parser, n_samples);
    end

    wl_pts = mwecmass.geometry.extract_isocurve_at_z( ...
                 parser, z_target, n_samples, b_cache);

    if isempty(wl_pts) || size(wl_pts, 1) < 3
        r_min = 0;
        return;
    end

    r_pts = wl_pts(:, 1:2);
    x_c = r_pts(:,1); y_c = r_pts(:,2);
    xp_c = circshift(x_c,-1); yp_c = circshift(y_c,-1);
    a_vec = x_c.*yp_c - xp_c.*y_c;
    A_poly = 0.5 * sum(a_vec);
    if abs(A_poly) > 1e-14
        centroid = [sum((x_c+xp_c).*a_vec)/(6*A_poly), ...
                    sum((y_c+yp_c).*a_vec)/(6*A_poly)];
    else
        centroid = mean(r_pts, 1);
    end
    dists = sqrt(sum((r_pts - centroid).^2, 2));
    r_min = min(dists);
end
