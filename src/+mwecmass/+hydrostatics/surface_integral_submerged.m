function [V_s, Cxyz_s, wl_pts, int_x2_s, int_y2_s, int_z2_s] = ...
        surface_integral_submerged( ...
        parser, surf_name, n_quad, z_wl, interior_pt, n_bisect, orient_sign, z_ref)
%SURFACE_INTEGRAL_SUBMERGED Integrate a surface portion below z_wl [m].
% Per-u submerged limits come from find_submerged_v_limits. z_ref [m] shifts
% z in the moments for numerical conditioning; returned z moments remain
% shifted for the caller to unshift. A supplied orient_sign fixes orientation
% from the complete surface, avoiding unreliable partial-domain checks.

    if nargin < 6, n_bisect = 40; end
    if nargin < 7, orient_sign = []; end
    if nargin < 8, z_ref = 0; end

    [u_gl, w_u] = mwecmass.internal.gauss_legendre(n_quad);
    [v_gl, w_v] = mwecmass.internal.gauss_legendre(n_quad);

    V_s      = 0;
    Cxyz_s   = [0, 0, 0];
    int_x2_s = 0;
    int_y2_s = 0;
    int_z2_s = 0;
    wl_pts   = [];

    for i = 1:n_quad
        u = u_gl(i);

        % Get the submerged v-interval for this u.
        [v_lo, v_hi, v_wl] = mwecmass.geometry.find_submerged_v_limits( ...
            parser, surf_name, u, z_wl, n_bisect);

        dv = v_hi - v_lo;
        if dv < 1e-14, continue; end

        % Collect waterplane boundary point (if crossing)
        if ~isnan(v_wl)
            wl_pts = [wl_pts; parser.eval_surface(surf_name, u, v_wl)]; %#ok<AGROW>
        end

        % Apply Gauss-Legendre quadrature over [v_lo, v_hi].
        for j = 1:n_quad
            v = v_lo + v_gl(j) * dv;

            [S, Su, Sv] = parser.eval_surface_with_derivs(surf_name, u, v);

            n_vec = cross(Su, Sv);
            w = w_u(i) * w_v(j) * dv;

            % Shifted z for numerical conditioning
            z_s = S(3) - z_ref;

            % Volume: (1/3) ∮ S'·n dA
            V_s = V_s + w * (S(1)*n_vec(1) + S(2)*n_vec(2) + z_s*n_vec(3)) / 3;

            % Centroid: x,y unshifted; z shifted
            Cxyz_s = Cxyz_s + w * [ ...
                S(1)^2 * n_vec(1) / 2, ...
                S(2)^2 * n_vec(2) / 2, ...
                z_s^2  * n_vec(3) / 2];

            % Second moments: x,y unshifted; z shifted
            int_x2_s = int_x2_s + w * S(1)^3 * n_vec(1) / 3;
            int_y2_s = int_y2_s + w * S(2)^3 * n_vec(2) / 3;
            int_z2_s = int_z2_s + w * z_s^3  * n_vec(3) / 3;
        end
    end

    % Apply the orientation sign from the complete surface.
    %  The sign is pre-determined from the FULL surface integral
    %  orientation check (dot product of normal at (0.5,0.5)
    %  with outward direction).  This must be applied consistently
    %  to all partial-domain integrals on this surface.
    %
    %    The raw GL integral ∫∫ S·(S_u×S_v) du dv over a partial
    %    v-range can have a different sign than the full integral.
    %    The orientation (inward vs outward normal) is a property
    %    of the parameterization, not the integration domain.
    %    Determining it from the partial integral gives wrong
    %    signs for some waterlines, making V_sub non-monotonic.

    if ~isempty(orient_sign)
        % Use pre-determined sign
        if orient_sign < 0
            V_s      = -V_s;
            Cxyz_s   = -Cxyz_s;
            int_x2_s = -int_x2_s;
            int_y2_s = -int_y2_s;
            int_z2_s = -int_z2_s;
        end
    else
        % Fallback: local check at (0.5, 0.5)
        [S_mid, Su_mid, Sv_mid] = parser.eval_surface_with_derivs(surf_name, 0.5, 0.5);
        n_mid = cross(Su_mid, Sv_mid);
        if dot(n_mid, S_mid - interior_pt(:)') < 0
            V_s      = -V_s;
            Cxyz_s   = -Cxyz_s;
            int_x2_s = -int_x2_s;
            int_y2_s = -int_y2_s;
            int_z2_s = -int_z2_s;
        end
    end
end
