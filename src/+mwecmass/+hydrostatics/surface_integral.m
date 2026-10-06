function [V_i, Cxyz_i, orient, int_x2_i, int_y2_i, int_z2_i] = ...
        surface_integral(parser, surf_name, n_quad, ~)
%SURFACE_INTEGRAL Integrate one parametric surface by the divergence theorem.
% Uses n_quad-by-n_quad Gauss-Legendre points on [0,1]^2 and analytic
% surface derivatives. Returns signed volume V [m^3], centroid numerators
% Cxyz, orientation sign, and raw second-volume moments [m^5]. A coarse
% signed-volume check selects the outward orientation.
% See docs/METHODS_ENGINE.md#divergence-theorem-surface-moments.

    [u_gl, w_u] = mwecmass.internal.gauss_legendre(n_quad);
    [v_gl, w_v] = mwecmass.internal.gauss_legendre(n_quad);

    V_i      = 0;
    Cxyz_i   = [0, 0, 0];
    int_x2_i = 0;
    int_y2_i = 0;
    int_z2_i = 0;

    for i = 1:n_quad
        u = u_gl(i);
        for j = 1:n_quad
            v = v_gl(j);

            % Surface point and analytical derivatives
            [S, Su, Sv] = parser.eval_surface_with_derivs(surf_name, u, v);

            % Oriented area element
            n_vec = cross(Su, Sv);

            % Quadrature weight
            w = w_u(i) * w_v(j);

            % Volume: (1/3) ∫ S · (S_u × S_v) du dv
            V_i = V_i + w * dot(S, n_vec) / 3;

            % Centroid numerators
            Cxyz_i = Cxyz_i + w * [ ...
                S(1)^2 * n_vec(1) / 2, ...
                S(2)^2 * n_vec(2) / 2, ...
                S(3)^2 * n_vec(3) / 2];

            % Second volume moments
            int_x2_i = int_x2_i + w * S(1)^3 * n_vec(1) / 3;
            int_y2_i = int_y2_i + w * S(2)^3 * n_vec(2) / 3;
            int_z2_i = int_z2_i + w * S(3)^3 * n_vec(3) / 3;
        end
    end

    % Determine orientation from a coarse signed-volume integral.
    % A coarse signed-volume integral is robust when a non-convex surface's
    % midpoint normal points inward despite the overall outward orientation.
    % A 5×5 signed-volume integral gives the correct sign of
    % ∫ S·(Su×Sv) du dv — it is always negative for an inward-
    % oriented parameterisation and positive for outward.
    n_orient_chk = 5;
    [u_chk, wu_chk] = mwecmass.internal.gauss_legendre(n_orient_chk);
    [v_chk, wv_chk] = mwecmass.internal.gauss_legendre(n_orient_chk);
    V_sign_test = 0;
    for ic = 1:n_orient_chk
        for jc = 1:n_orient_chk
            [S_c, Su_c, Sv_c] = parser.eval_surface_with_derivs( ...
                surf_name, u_chk(ic), v_chk(jc));
            V_sign_test = V_sign_test + wu_chk(ic) * wv_chk(jc) * ...
                dot(S_c, cross(Su_c, Sv_c)) / 3;
        end
    end

    if V_sign_test < 0
        V_i      = -V_i;
        Cxyz_i   = -Cxyz_i;
        int_x2_i = -int_x2_i;
        int_y2_i = -int_y2_i;
        int_z2_i = -int_z2_i;
        orient   = -1;
    else
        orient = +1;
    end
end
