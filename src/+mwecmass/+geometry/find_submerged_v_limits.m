function [v_lo, v_hi, v_wl] = find_submerged_v_limits( ...
        parser, surf_name, u, z_wl, n_bisect)
%FIND_SUBMERGED_V_LIMITS Find the v interval below a waterline.
%   z_wl is in m; returns [v_lo,v_hi] and the first crossing v_wl.

    if nargin < 5, n_bisect = 40; end

    % ── RevSurf fast-path ─────────────────────────────────
    %  For RevSurf, z depends on u only — revolution preserves z.
    %  A single z-check at any v suffices.  No v-sweep or bisection.
    %  For C0: surface2 accounts for ~33% of v-limit searches.
    if parser.entities.isKey(surf_name)
        se = parser.entities(surf_name);
        if strcmp(se.type, 'RevSurf')
            pt = parser.eval_surface(surf_name, u, 0);
            v_wl = NaN;
            if pt(3) <= z_wl + 1e-10
                v_lo = 0; v_hi = 1;
            else
                v_lo = 0; v_hi = 0;
            end
            return;
        end
        % MirrSurf of a RevSurf also preserves z under revolution
        if strcmp(se.type, 'MirrSurf')
            src = parser.entities(se.params.source);
            if strcmp(src.type, 'RevSurf')
                pt = parser.eval_surface(surf_name, u, 0);
                v_wl = NaN;
                if pt(3) <= z_wl + 1e-10
                    v_lo = 0; v_hi = 1;
                else
                    v_lo = 0; v_hi = 0;
                end
                return;
            end
        end
    end

    % ── Coarse sweep for sign change ────────────────────────
    n_sample = 21;
    v_sample = linspace(0, 1, n_sample);
    f_sample = zeros(n_sample, 1);
    for k = 1:n_sample
        pt = parser.eval_surface(surf_name, u, v_sample(k));
        f_sample(k) = pt(3) - z_wl;
    end

    % ── No crossing: all below or all above ─────────────────
    sign_changes = find(f_sample(1:end-1) .* f_sample(2:end) < 0, 1);

    if isempty(sign_changes)
        v_wl = NaN;
        if f_sample(1) <= 0
            v_lo = 0; v_hi = 1;   % all below
        else
            v_lo = 0; v_hi = 0;   % all above
        end
        return;
    end

    % ── Bisection to refine crossing ────────────────────────
    a = v_sample(sign_changes);
    b = v_sample(sign_changes + 1);
    fa = f_sample(sign_changes);

    for iter = 1:n_bisect
        m = (a + b) / 2;
        pt = parser.eval_surface(surf_name, u, m);
        fm = pt(3) - z_wl;
        if abs(fm) < 1e-10, break; end
        if fa * fm < 0
            b = m;
        else
            a = m; fa = fm;
        end
    end

    v_wl = (a + b) / 2;

    % ── Which side is submerged? ────────────────────────────
    if f_sample(1) <= 0
        v_lo = 0;     v_hi = v_wl;    % below is v < v_wl
    else
        v_lo = v_wl;  v_hi = 1;        % below is v > v_wl
    end
end
