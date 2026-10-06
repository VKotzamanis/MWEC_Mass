function [V_cap, Cxyz_cap, A_cap, z_cap] = cap_contribution( ...
        open_edges, interior_pt, tol)
%CAP_CONTRIBUTION Compute volume and centroid numerators for open-edge caps.
% open_edges are chained into loops; each nearly flat loop is assigned area
% by the shoelace formula and signed cap contributions by its z relative to
% interior_pt. tol defaults to 1e-3 [m]. Outputs are per-loop V [m^3],
% centroid numerator [m^4], area [m^2], and cap elevation [m].

    if nargin < 3, tol = 1e-3; end

    % Chain open edges into loop(s)
    groups = mwecmass.geometry.chain_open_edges(open_edges, tol);
    n_groups = length(groups);

    V_cap    = zeros(n_groups, 1);
    Cxyz_cap = zeros(n_groups, 3);
    A_cap    = zeros(n_groups, 1);
    z_cap    = zeros(n_groups, 1);

    for g = 1:n_groups
        loop = groups{g}.loop_pts;
        z_g  = groups{g}.z_mean;

        if groups{g}.z_range > 0.05
            warning('mwecmass:hydrostatics:NonFlatCap', ...
                    'Cap %d z-range = %.4f m (not flat)', ...
                    g, groups{g}.z_range);
        end

        % Compute signed area with the shoelace formula.
        %  This IS Green's theorem evaluated discretely:
        %    A = (1/2) Σ (x_i y_{i+1} − x_{i+1} y_i)
        x = loop(:, 1);
        y = loop(:, 2);
        x_next = [x(2:end); x(1)];
        y_next = [y(2:end); y(1)];
        A_signed = sum(x .* y_next - x_next .* y) / 2;

        A_cap(g) = abs(A_signed);
        z_cap(g) = z_g;

        % Determine the outward normal direction.
        %  If z_cap > interior_z → top cap → outward = +z
        %  If z_cap < interior_z → bottom cap → outward = −z
        if z_g > interior_pt(3)
            sign_nz = +1;
        else
            sign_nz = -1;
        end

        % The cap area is always the unsigned polygon area.
        % The sign of the volume contribution comes from sign_nz.
        A_use = abs(A_signed);

        % V_cap = (sign_nz / 3) * z_cap * A_cap
        V_cap(g) = sign_nz * z_g * A_use / 3;

        % Centroid: only z-component (n_x = n_y = 0)
        % z̄·V contribution = (sign_nz / 2) * z_cap² * A_cap
        Cxyz_cap(g, :) = [0, 0, sign_nz * z_g^2 * A_use / 2];
    end
end
