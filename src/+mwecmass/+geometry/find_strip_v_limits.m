function [v_lo, v_hi] = find_strip_v_limits( ...
        parser, surf_name, u, z_lo, z_hi, n_bisect)
%FIND_STRIP_V_LIMITS Find the v interval within an elevation strip.
%   Bounds are in m; returns [v_lo,v_hi] in parameter space.

    if nargin < 6, n_bisect = 40; end

    % ── RevSurf fast-path ─────────────────────────────────
    if parser.entities.isKey(surf_name)
        se = parser.entities(surf_name);
        if strcmp(se.type, 'RevSurf')
            pt = parser.eval_surface(surf_name, u, 0);
            z_u = pt(3);
            if z_u >= z_lo - 1e-10 && z_u <= z_hi + 1e-10
                v_lo = 0; v_hi = 1;
            else
                v_lo = 0; v_hi = 0;
            end
            return;
        end
        if strcmp(se.type, 'MirrSurf')
            src = parser.entities(se.params.source);
            if strcmp(src.type, 'RevSurf')
                pt = parser.eval_surface(surf_name, u, 0);
                z_u = pt(3);
                if z_u >= z_lo - 1e-10 && z_u <= z_hi + 1e-10
                    v_lo = 0; v_hi = 1;
                else
                    v_lo = 0; v_hi = 0;
                end
                return;
            end
        end
    end

    % ── Coarse sweep ──────────────────────────────────────
    n_sample = 21;
    v_sample = linspace(0, 1, n_sample);
    z_sample = zeros(n_sample, 1);
    for k = 1:n_sample
        pt = parser.eval_surface(surf_name, u, v_sample(k));
        z_sample(k) = pt(3);
    end

    % Classify: find first and last v where z is inside [z_lo, z_hi]
    inside = (z_sample >= z_lo - 1e-10) & (z_sample <= z_hi + 1e-10);

    if ~any(inside)
        v_lo = 0; v_hi = 0;
        return;
    end

    if all(inside)
        v_lo = 0; v_hi = 1;
        return;
    end

    % Find transitions: outside→inside and inside→outside
    first_in = find(inside, 1, 'first');
    last_in  = find(inside, 1, 'last');

    % Refine lower boundary
    if first_in == 1
        v_lo = 0;
    else
        % Bisect between v_sample(first_in-1) and v_sample(first_in)
        a = v_sample(first_in - 1);
        b = v_sample(first_in);
        for iter = 1:n_bisect
            m = (a + b) / 2;
            pt = parser.eval_surface(surf_name, u, m);
            if pt(3) >= z_lo - 1e-10 && pt(3) <= z_hi + 1e-10
                b = m;  % inside — tighten from right
            else
                a = m;  % outside — tighten from left
            end
        end
        v_lo = (a + b) / 2;
    end

    % Refine upper boundary
    if last_in == n_sample
        v_hi = 1;
    else
        a = v_sample(last_in);
        b = v_sample(last_in + 1);
        for iter = 1:n_bisect
            m = (a + b) / 2;
            pt = parser.eval_surface(surf_name, u, m);
            if pt(3) >= z_lo - 1e-10 && pt(3) <= z_hi + 1e-10
                a = m;  % inside — tighten from left
            else
                b = m;  % outside — tighten from right
            end
        end
        v_hi = (a + b) / 2;
    end
end
