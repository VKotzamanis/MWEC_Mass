function inner_prof = compute_inner_profile(outer_profile, t_shell)
%COMPUTE_INNER_PROFILE Compute an inner profile offset inward from an outer profile by shell thickness.
% outer_profile is n-by-2 [x, z] vertices; t_shell is thickness [m]. Returns
% [n x 2] inner profile with vertices offset normal by t_shell; clamps miter
% at 3*t_shell where vertex angle exceeds 70 degrees.
    x = outer_profile(:,1);
    z = outer_profile(:,2);
    Nv = length(x);
    if abs(x(end)-x(1)) < 1e-12 && abs(z(end)-z(1)) < 1e-12
        x = x(1:end-1); z = z(1:end-1); Nv = Nv - 1;
    end
    if Nv < 3, inner_prof = outer_profile; return; end
    signed_area = 0.5 * sum(x .* circshift(z,-1) - circshift(x,-1) .* z);
    if signed_area < 0
        x = flipud(x); z = flipud(z);
    end
    x_off = zeros(Nv, 1);
    z_off = zeros(Nv, 1);
    for j = 1:Nv
        jm = mod(j-2, Nv) + 1;
        jp = mod(j,   Nv) + 1;
        e_prev = [x(j)-x(jm), z(j)-z(jm)];
        e_next = [x(jp)-x(j), z(jp)-z(j)];
        lp = norm(e_prev); ln = norm(e_next);
        if lp < 1e-12 || ln < 1e-12
            x_off(j) = x(j); z_off(j) = z(j); continue;
        end
        n_prev = [-e_prev(2),  e_prev(1)] / lp;
        n_next = [-e_next(2),  e_next(1)] / ln;
        n_avg    = n_prev + n_next;
        len_avg  = norm(n_avg);
        if len_avg < 1e-12, n_avg = n_prev; len_avg = 1; end
        n_avg    = n_avg / len_avg;
        cos_half = dot(n_avg, n_prev);
        if cos_half > 0.33
            miter = t_shell / cos_half;
        else
            miter = t_shell * 3.0;
        end
        x_off(j) = x(j) + n_avg(1) * miter;
        z_off(j) = z(j) + n_avg(2) * miter;
    end
    inner_prof = [x_off, z_off];
end
