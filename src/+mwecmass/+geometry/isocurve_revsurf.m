function pts = isocurve_revsurf(parser, z_wl, n_arc, cache_entry)
%ISOCURVE_REVSURF Compute a RevSurf iso-z arc from cached geometry.
%   z_wl and returned [N x 3] points are in m; n_arc controls sampling.
%   See docs/METHODS_ENGINE.md#geometry-revolution-isocurve.

    d = cache_entry;
    u_grid = linspace(0, 1, length(d.z_profile))';

    u_star = mwecmass.geometry.profile_param_at_z( ...
                 d.z_profile, u_grid, z_wl, parser, d.profile_name);
    if isnan(u_star)
        pts = [];
        return;
    end

    % Exact profile point at u* (one eval_curve_or_snake call)
    profile_pt = parser.eval_curve_or_snake(d.profile_name, u_star);
    if size(profile_pt, 1) > 1, profile_pt = profile_pt(1,:); end

    % Axial + radial decomposition using cached axis
    v_rel = profile_pt - d.axis_start;
    z_along = dot(v_rel, d.axis_dir);
    proj = d.axis_start + z_along * d.axis_dir;
    radial = profile_pt - proj;
    r = norm(radial);

    if r < 1e-12
        pts = repmat(profile_pt, n_arc, 1);
        return;
    end

    e_r = radial / r;
    e_t = cross(d.axis_dir, e_r);
    e_t = e_t / norm(e_t);

    % Analytic arc generation
    phi_start = deg2rad(d.angle_start);
    phi_end   = deg2rad(d.angle_end);
    phi = linspace(phi_start, phi_end, n_arc)';

    pts = proj + r * cos(phi) .* e_r + r * sin(phi) .* e_t;
end
