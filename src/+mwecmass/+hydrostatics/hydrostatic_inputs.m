function [CG_global, M_6x6, C_6x6, mass, sub_props] = hydrostatic_inputs(parser, z_wl, config)
%HYDROSTATIC_INPUTS Build uniform-density Hydrostatic.in quantities.
% Inputs are parser geometry, waterline z_wl [m], and config constants/tables.
% Outputs: CG_global [1x3 m], M_6x6 and C_6x6, mass [kg], and submerged
% properties (V_sub, CB, Aw, waterplane moments). Uses configured area tables
% when available and the submerged-surface path otherwise.

    rho_w = config.RHO_WATER;
    g     = config.G;
    V_total  = config.total_wec_volume;
    centroid = config.hull_centroid;
    draft = -z_wl;   % m, positive down; waterline depth for the BEM deck, NOT the hull draft abs(hull_z_min+vs) used elsewhere.

    % 1. Submerged properties from precomputed tables (v9.0)
    %
    %  Uses the configured Aw-table trapezoidal integration path.
    %  Bypasses compute_submerged entirely — avoids the divergence-
    %  theorem partial-surface errors for offset-axis geometries.
    %
    %  Fallback: if V_sub_table is absent (geometry-only config or
    %  config struct), calls compute_submerged as before.
    if isfield(config, 'V_sub_table') && ~isempty(config.V_sub_table)
        z_wl_clamped = max(min(z_wl, config.Aw_table_z(end)), ...
                               config.Aw_table_z(1));
        sub_props = struct();
        sub_props.V_sub   = max(0, interp1(config.Aw_table_z, ...
                            config.V_sub_table,   z_wl_clamped, 'linear'));
        sub_props.CB      = [0, 0, interp1(config.Aw_table_z, ...
                            config.CB_z_table,    z_wl_clamped, 'linear')];
        sub_props.Aw      = max(0, interp1(config.Aw_table_z, ...
                            config.Aw_table,      z_wl_clamped, 'linear'));
        sub_props.I_wp_xx = max(0, interp1(config.Aw_table_z, ...
                            config.I_wp_xx_table, z_wl_clamped, 'linear'));
        sub_props.I_wp_yy = max(0, interp1(config.Aw_table_z, ...
                            config.I_wp_yy_table, z_wl_clamped, 'linear'));
    else
        % Fallback: compute_submerged (may have errors for
        % offset-axis geometries — retained for backward compatibility)
        sub_opts = struct('n_quad', 16, 'verbose', false);
        if isfield(config, 'Aw_table_z') && ~isempty(config.Aw_table_z)
            sub_opts.Aw_override = interp1( ...
                config.Aw_table_z, config.Aw_table, z_wl, 'linear', 0);
            sub_opts.I_wp_xx_override = interp1( ...
                config.Aw_table_z, config.I_wp_xx_table, z_wl, 'linear', 0);
            sub_opts.I_wp_yy_override = interp1( ...
                config.Aw_table_z, config.I_wp_yy_table, z_wl, 'linear', 0);
        end
        sub_props = mwecmass.hydrostatics.compute_submerged(parser, z_wl, sub_opts);
    end

    % 2. Mass (uniform density, Archimedes)
    mass = rho_w * sub_props.V_sub;
    rho_body = mass / V_total;

    % 3. CG in global frame
    %   Bilateral symmetry: CG must lie on the z-axis.
    %   The parametric centroid has x,y ≈ 0 (O(1e-17) from
    %   numerical noise in the divergence theorem). Force exact.
    CG_global = [0, 0, centroid(3) + draft];

    % 4. Inertia tensor about CG (frame-invariant)
    %   Use cx = cy = 0 (bilateral symmetry) to avoid noise.
    cx = 0; cy = 0; cz = centroid(3);
    Ixx_cg = max(0, rho_body * (config.hull_int_y2 + config.hull_int_z2 ...
             - V_total * (cy^2 + cz^2)));
    Iyy_cg = max(0, rho_body * (config.hull_int_x2 + config.hull_int_z2 ...
             - V_total * (cx^2 + cz^2)));
    Izz_cg = max(0, rho_body * (config.hull_int_x2 + config.hull_int_y2 ...
             - V_total * (cx^2 + cy^2)));
    J_CG = diag([Ixx_cg, Iyy_cg, Izz_cg]);

    % 5. 6x6 mass matrix about origin
    M_6x6 = mwecmass.hydrostatics.build_mass_matrix(mass, CG_global, J_CG);

    % 6. Restoring matrix about origin
    z_B = sub_props.CB(3) + draft;
    z_G = CG_global(3);
    C_6x6 = zeros(6);
    C_6x6(3,3) = rho_w * g * sub_props.Aw;
    C_6x6(4,4) = rho_w*g*sub_props.I_wp_xx + rho_w*g*sub_props.V_sub*z_B - mass*g*z_G;
    C_6x6(5,5) = rho_w*g*sub_props.I_wp_yy + rho_w*g*sub_props.V_sub*z_B - mass*g*z_G;
end
