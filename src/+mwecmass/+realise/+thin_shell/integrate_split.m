function R = integrate_split(grids, z_ballast)
%INTEGRATE_SPLIT Integrate geometry-grid areas and moments below/above z_ballast.
%   R contains trapezoidal volumes [m^3], z-area integrals [m^4/m^5], and
%   x/y second-moment volume integrals [m^5] for outer, inner, and jacket
%   regions. z_ballast is clamped to the grid and inserted by linear interpolation.
%   See docs/METHODS_ENGINE.md#thin-shell-split-integrals

    z      = grids.z(:);
    A_o    = grids.A_outer(:);
    A_i    = grids.A_inner(:);
    A_j    = A_o - A_i;
    Iyy_o  = grids.Iyy_outer(:);
    Iyy_i  = grids.Iyy_inner(:);
    Iyy_j  = Iyy_o - Iyy_i;
    Ixx_o  = grids.Ixx_outer(:);
    Ixx_i  = grids.Ixx_inner(:);
    Ixx_j  = Ixx_o - Ixx_i;

    % Clamp z_ballast to grid range
    z_ballast = max(z(1), min(z(end), z_ballast));

    below_mask = z <= z_ballast;
    above_mask = z >= z_ballast;
    z_below = z(below_mask);
    z_above = z(above_mask);

    tol = 1e-12;
    if isempty(z_below) || abs(z_below(end) - z_ballast) > tol
        A_o_zf   = interp1(z, A_o,   z_ballast, 'linear');
        A_i_zf   = interp1(z, A_i,   z_ballast, 'linear');
        A_j_zf   = interp1(z, A_j,   z_ballast, 'linear');
        Iyy_o_zf = interp1(z, Iyy_o, z_ballast, 'linear');
        Iyy_i_zf = interp1(z, Iyy_i, z_ballast, 'linear');
        Iyy_j_zf = interp1(z, Iyy_j, z_ballast, 'linear');
        Ixx_o_zf = interp1(z, Ixx_o, z_ballast, 'linear');
        Ixx_i_zf = interp1(z, Ixx_i, z_ballast, 'linear');
        Ixx_j_zf = interp1(z, Ixx_j, z_ballast, 'linear');

        z_below   = [z_below;  z_ballast];
        A_o_below = [A_o(below_mask);  A_o_zf];
        A_i_below = [A_i(below_mask);  A_i_zf];  %#ok<NASGU>
        A_j_below = [A_j(below_mask);  A_j_zf];  %#ok<NASGU>
        Iyy_o_below = [Iyy_o(below_mask); Iyy_o_zf];
        Iyy_i_below = [Iyy_i(below_mask); Iyy_i_zf];  %#ok<NASGU>
        Iyy_j_below = [Iyy_j(below_mask); Iyy_j_zf];  %#ok<NASGU>
        Ixx_o_below = [Ixx_o(below_mask); Ixx_o_zf];
        Ixx_i_below = [Ixx_i(below_mask); Ixx_i_zf];  %#ok<NASGU>
        Ixx_j_below = [Ixx_j(below_mask); Ixx_j_zf];  %#ok<NASGU>

        z_above   = [z_ballast;  z_above];
        A_o_above = [A_o_zf;     A_o(above_mask)];  %#ok<NASGU>
        A_i_above = [A_i_zf;     A_i(above_mask)];
        A_j_above = [A_j_zf;     A_j(above_mask)];
        Iyy_o_above = [Iyy_o_zf; Iyy_o(above_mask)];  %#ok<NASGU>
        Iyy_i_above = [Iyy_i_zf; Iyy_i(above_mask)];
        Iyy_j_above = [Iyy_j_zf; Iyy_j(above_mask)];
        Ixx_o_above = [Ixx_o_zf; Ixx_o(above_mask)];  %#ok<NASGU>
        Ixx_i_above = [Ixx_i_zf; Ixx_i(above_mask)];
        Ixx_j_above = [Ixx_j_zf; Ixx_j(above_mask)];
    else
        A_o_below = A_o(below_mask); A_i_below = A_i(below_mask); A_j_below = A_j(below_mask);  %#ok<NASGU>
        Iyy_o_below = Iyy_o(below_mask); Iyy_i_below = Iyy_i(below_mask); Iyy_j_below = Iyy_j(below_mask);  %#ok<NASGU>
        Ixx_o_below = Ixx_o(below_mask); Ixx_i_below = Ixx_i(below_mask); Ixx_j_below = Ixx_j(below_mask);  %#ok<NASGU>
        A_o_above = A_o(above_mask); A_i_above = A_i(above_mask); A_j_above = A_j(above_mask);  %#ok<NASGU>
        Iyy_o_above = Iyy_o(above_mask); Iyy_i_above = Iyy_i(above_mask); Iyy_j_above = Iyy_j(above_mask);  %#ok<NASGU>
        Ixx_o_above = Ixx_o(above_mask); Ixx_i_above = Ixx_i(above_mask); Ixx_j_above = Ixx_j(above_mask);  %#ok<NASGU>
    end

    R = struct();

    % --- BELOW z_ballast: cross-section is fully solid (use A_outer everywhere) ---
    R.V_below_outer       = tz(z_below, A_o_below);
    R.int_zA_below_outer  = tz(z_below, z_below .* A_o_below);
    R.int_z2A_below_outer = tz(z_below, z_below.^2 .* A_o_below);
    R.int_Ix_below_outer  = tz(z_below, Iyy_o_below);   % ∫∫∫ x² dV
    R.int_Iy_below_outer  = tz(z_below, Ixx_o_below);   % ∫∫∫ y² dV

    % --- ABOVE z_ballast: jacket annulus is steel; inner area is air ---
    R.V_above_jacket       = tz(z_above, A_j_above);
    R.int_zA_above_jacket  = tz(z_above, z_above .* A_j_above);
    R.int_z2A_above_jacket = tz(z_above, z_above.^2 .* A_j_above);
    R.int_Ix_above_jacket  = tz(z_above, Iyy_j_above);
    R.int_Iy_above_jacket  = tz(z_above, Ixx_j_above);

    R.V_above_inner       = tz(z_above, A_i_above);
    R.int_zA_above_inner  = tz(z_above, z_above .* A_i_above);
    R.int_z2A_above_inner = tz(z_above, z_above.^2 .* A_i_above);
    R.int_Ix_above_inner  = tz(z_above, Iyy_i_above);
    R.int_Iy_above_inner  = tz(z_above, Ixx_i_above);
end

function v = tz(x, y)
%TZ Safe trapezoidal integration: returns 0 for fewer than two sample points.
    if length(x) < 2
        v = 0;
    else
        v = trapz(x, y);
    end
end
