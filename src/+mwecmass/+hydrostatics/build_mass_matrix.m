function M = build_mass_matrix(m, r_G, J_CG)
%BUILD_MASS_MATRIX Build the 6-by-6 rigid-body mass matrix about the origin.
% Inputs: m [kg], r_G [3x1 m] CG position, and J_CG [3x3 kg m^2] inertia
% tensor using the negative-product convention. Output M uses [surge sway
% heave roll pitch yaw] order and includes the parallel-axis and coupling terms.
% For m <= 0, returns eye(6) and emits the existing warning.

    if m <= 0
        warning('mwecmass:hydrostatics:NonPositiveMass', ...
                ['Body mass is %g kg (<= 0), so the 6x6 mass matrix is undefined; ' ...
                 'returning the identity, which carries neither the mass nor the inertia.'], m);
        M = eye(6);
        return;
    end

    r_G = r_G(:);  % ensure column
    x = r_G(1); y = r_G(2); z = r_G(3);

    % Steiner correction: J^O = J_CG + m*((r.r)*I - r*r')
    d2 = x^2 + y^2 + z^2;
    J_O = J_CG + m * (d2 * eye(3) - r_G * r_G');

    % Skew-symmetric coupling: m * S(r_G)
    %   S(r) = [  0   z  -y ]
    %          [ -z   0   x ]
    %          [  y  -x   0 ]
    C_block = m * [  0,  z, -y;
                  -z,  0,  x;
                   y, -x,  0 ];

    % Assemble 6x6 (symmetric)
    M = [ m*eye(3),   C_block;
          C_block',   J_O     ];
end
