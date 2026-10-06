function [V_s, Cz_s, ix2_s, iy2_s, iz2_s] = surface_integral_strip( ...
        parser, surf_name, n_quad, z_lo, z_hi, interior_pt, orient_sign, z_ref)
%SURFACE_INTEGRAL_STRIP Integrate a surface over the horizontal band [z_lo,z_hi].
% Per-u limits come from find_strip_v_limits. z_ref [m] shifts z in the
% divergence-theorem moments for conditioning; orient_sign, when supplied,
% fixes the full-surface orientation. Returns volume and raw second moments.

    if nargin < 7, orient_sign = []; end
    if nargin < 8, z_ref = 0; end
    n_bisect = 40;
    [u_gl, w_u] = mwecmass.internal.gauss_legendre(n_quad);
    [v_gl, w_v] = mwecmass.internal.gauss_legendre(n_quad);

    V_s   = 0;
    Cz_s  = 0;
    ix2_s = 0;  iy2_s = 0;  iz2_s = 0;

    for i = 1:n_quad
        u = u_gl(i);

        [v_lo_i, v_hi_i] = mwecmass.geometry.find_strip_v_limits( ...
            parser, surf_name, u, z_lo, z_hi, n_bisect);

        dv = v_hi_i - v_lo_i;
        if dv < 1e-14, continue; end

        for j = 1:n_quad
            v = v_lo_i + v_gl(j) * dv;
            [S, Su, Sv] = parser.eval_surface_with_derivs(surf_name, u, v);
            n_vec = cross(Su, Sv);
            w = w_u(i) * w_v(j) * dv;

            z_s = S(3) - z_ref;

            V_s   = V_s + w * (S(1)*n_vec(1) + S(2)*n_vec(2) + z_s*n_vec(3)) / 3;
            Cz_s  = Cz_s + w * z_s^2  * n_vec(3) / 2;
            ix2_s = ix2_s + w * S(1)^3 * n_vec(1) / 3;
            iy2_s = iy2_s + w * S(2)^3 * n_vec(2) / 3;
            iz2_s = iz2_s + w * z_s^3  * n_vec(3) / 3;
        end
    end

    % Apply pre-determined orientation sign
    if ~isempty(orient_sign)
        if orient_sign < 0
            V_s   = -V_s;
            Cz_s  = -Cz_s;
            ix2_s = -ix2_s;
            iy2_s = -iy2_s;
            iz2_s = -iz2_s;
        end
    else
        [S_mid, Su_mid, Sv_mid] = parser.eval_surface_with_derivs( ...
            surf_name, 0.5, 0.5);
        n_mid = cross(Su_mid, Sv_mid);
        if dot(n_mid, S_mid - interior_pt(:)') < 0
            V_s   = -V_s;
            Cz_s  = -Cz_s;
            ix2_s = -ix2_s;
            iy2_s = -iy2_s;
            iz2_s = -iz2_s;
        end
    end
end
