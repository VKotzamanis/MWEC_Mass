function A_sub = compute_wetted_surface_area(parser, z_wl, n_quad)
%COMPUTE_WETTED_SURFACE_AREA Compute wetted hull area below z_wl [m].
% Integrates ||S_u x S_v|| over submerged source surfaces and adds reflected
% mirrors. The waterplane is excluded by definition. n_quad is the
% Gauss-Legendre order (default 20); parser is the MS2 geometry model.

    if nargin < 3 || isempty(n_quad), n_quad = 20; end

    [u_gl, w_u] = mwecmass.internal.gauss_legendre(n_quad);
    [v_gl, w_v] = mwecmass.internal.gauss_legendre(n_quad);

    topo      = parser.classify_visible_surfaces();
    n_bisect  = 40;
    A_sub     = 0;

    % Map to hold source wetted area (needed for mirror lookup)
    source_A_sub = containers.Map('KeyType','char','ValueType','double');

    % Integrate source surfaces.
    for s = 1:length(topo.sources)
        sname  = topo.sources{s};
        A_surf = 0;

        for i = 1:n_quad
            u = u_gl(i);

            [v_lo, v_hi, ~] = mwecmass.geometry.find_submerged_v_limits( ...
                parser, sname, u, z_wl, n_bisect);

            dv = v_hi - v_lo;
            if dv < 1e-14, continue; end

            for j = 1:n_quad
                v = v_lo + v_gl(j) * dv;
                [~, Su, Sv] = parser.eval_surface_with_derivs(sname, u, v);
                n_vec  = cross(Su, Sv);
                w      = w_u(i) * w_v(j) * dv;
                A_surf = A_surf + w * norm(n_vec);
            end
        end

        source_A_sub(sname) = A_surf;
        A_sub = A_sub + A_surf;
    end

    % Add mirrored surfaces.
    %   Reflection preserves distances → A_mirror = A_source.
    %   Iterates topo.mirrors exactly as compute_submerged does —
    %   handles 1, 2, or 3 mirrors per source (not a simple ×2).
    for m = 1:length(topo.mirrors)
        mirr = topo.mirrors(m);
        ult  = mirr.ultimate_source;

        if source_A_sub.isKey(ult)
            A_sub = A_sub + source_A_sub(ult);
        else
            A_mirr = 0;
            for i = 1:n_quad
                u = u_gl(i);
                [v_lo, v_hi, ~] = mwecmass.geometry.find_submerged_v_limits( ...
                    parser, mirr.name, u, z_wl, n_bisect);
                dv = v_hi - v_lo;
                if dv < 1e-14, continue; end
                for j = 1:n_quad
                    v = v_lo + v_gl(j) * dv;
                    [~, Su, Sv] = parser.eval_surface_with_derivs(mirr.name, u, v);
                    n_vec  = cross(Su, Sv);
                    w      = w_u(i) * w_v(j) * dv;
                    A_mirr = A_mirr + w * norm(n_vec);
                end
            end
            A_sub = A_sub + A_mirr;
        end
    end

    % Waterplane cap: deliberately NOT added.
end
