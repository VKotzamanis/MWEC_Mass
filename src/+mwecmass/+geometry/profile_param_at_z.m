function u_star = profile_param_at_z(z_profile, u_grid, z_wl, parser, profile_name)
%PROFILE_PARAM_AT_Z Find a profile parameter at a target elevation.
%   z_wl and profile elevations are in m; returns the first crossing or NaN.

    f = z_profile - z_wl;
    u_star = NaN;

    for k = 1:length(f) - 1
        if f(k) * f(k + 1) <= 0
            a  = u_grid(k);
            b  = u_grid(k + 1);
            fa = f(k);

            if nargin >= 5 && ~isempty(parser) && ~isempty(profile_name)
                for iter = 1:40
                    m = (a + b) / 2;
                    pt_m = parser.eval_curve_or_snake(profile_name, m);
                    fm = pt_m(1, 3) - z_wl;
                    if abs(fm) < 1e-12, break; end
                    if fa * fm < 0
                        b = m;
                    else
                        a = m; fa = fm;
                    end
                end
                u_star = (a + b) / 2;
            else
                if abs(f(k) - f(k+1)) < 1e-14
                    u_star = a;
                else
                    u_star = a - fa * (b - a) / (f(k+1) - f(k));
                end
            end
            return;
        end
    end
end
