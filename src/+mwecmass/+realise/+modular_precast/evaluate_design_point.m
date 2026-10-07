function out = evaluate_design_point(vs, t_offset_strip, is_solid_strip, ...
                                        z_ballast, grids, ctx)
%EVALUATE_DESIGN_POINT Evaluate UHPC/void hydrostatics, mass properties, GM, and periods.
% grids already encode per-strip thickness and solid masks; vs and z_ballast define the realised pose.
% Returns the standard realised-property struct, with infeasible/non-positive-mass points flagged.
    cfg = ctx.config;
    out = mwecmass.realise.empty_realised_properties();

    R = mwecmass.realise.thin_shell.integrate_split(grids, z_ballast);

    V_uhpc = R.V_below_outer + R.V_above_jacket;
    V_air   = R.V_above_inner;

    int_z_x_A_uhpc = R.int_zA_below_outer + R.int_zA_above_jacket;
    int_z_x_A_air   = R.int_zA_above_inner;
    int_x2_uhpc    = R.int_Ix_below_outer + R.int_Ix_above_jacket;
    int_x2_air      = R.int_Ix_above_inner;
    int_y2_uhpc    = R.int_Iy_below_outer + R.int_Iy_above_jacket;
    int_y2_air      = R.int_Iy_above_inner;
    int_z2_A_uhpc  = R.int_z2A_below_outer + R.int_z2A_above_jacket;
    int_z2_A_air    = R.int_z2A_above_inner;

    M_uhpc = ctx.rho_uhpc * V_uhpc;
    M_air   = ctx.rho_air   * V_air;
    M_total = M_uhpc + M_air;

    out.V_uhpc = V_uhpc;
    out.V_air   = V_air;
    out.M_uhpc = M_uhpc;
    out.M_air   = M_air;
    out.M_total = M_total;

    if M_total <= 0
        out.feasible = false;
        return;
    end

    if V_uhpc > 1e-12
        z_cg_uhpc = int_z_x_A_uhpc / V_uhpc;
    else
        z_cg_uhpc = 0.5 * (ctx.hull_z_min + z_ballast);
    end
    if V_air > 1e-12
        z_cg_air = int_z_x_A_air / V_air;
    else
        z_cg_air = 0.5 * (z_ballast + ctx.hull_z_max);
    end
    CG_z_body = (M_uhpc * z_cg_uhpc + M_air * z_cg_air) / M_total;

    Iyy_total_origin = ctx.rho_uhpc * (int_x2_uhpc + int_z2_A_uhpc) + ...
                       ctx.rho_air   * (int_x2_air   + int_z2_A_air);
    Ixx_total_origin = ctx.rho_uhpc * (int_y2_uhpc + int_z2_A_uhpc) + ...
                       ctx.rho_air   * (int_y2_air   + int_z2_A_air);
    Izz_total_origin = ctx.rho_uhpc * (int_x2_uhpc + int_y2_uhpc) + ...
                       ctx.rho_air   * (int_x2_air   + int_y2_air);
    Iyy_about_cg = max(0, Iyy_total_origin - M_total * CG_z_body^2);
    Ixx_about_cg = max(0, Ixx_total_origin - M_total * CG_z_body^2);
    Izz_about_cg = max(0, Izz_total_origin);

    out.z_cg_uhpc       = z_cg_uhpc;
    out.z_cg_air         = z_cg_air;
    out.CG_z_body        = CG_z_body;
    out.Iyy_total_origin = Iyy_total_origin;
    out.Iyy_about_cg     = Iyy_about_cg;
    out.Ixx_total_origin = Ixx_total_origin;
    out.Ixx_about_cg     = Ixx_about_cg;
    out.Izz_total_origin = Izz_total_origin;
    out.Izz_about_cg     = Izz_about_cg;

    out.draft          = abs(ctx.hull_z_min + vs);
    out.vertical_shift = vs;
    out.CG_z_world     = CG_z_body + vs;

    z_wl_body = -vs;
    z_sub_top = min(z_wl_body, ctx.hull_z_max);

    V_sub      = max(0, interp1(cfg.Aw_table_z, cfg.V_sub_table,    z_sub_top, 'linear', 0));
    Aw         = max(0, interp1(cfg.Aw_table_z, cfg.Aw_table,       z_sub_top, 'linear', 0));
    I_wp_yy    = max(0, interp1(cfg.Aw_table_z, cfg.I_wp_yy_table,  z_sub_top, 'linear', 0));
    CB_z_body  =        interp1(cfg.Aw_table_z, cfg.CB_z_table,     z_sub_top, 'linear', 'extrap');
    CB_z_world = CB_z_body + vs;

    out.V_sub              = V_sub;
    out.Aw                 = Aw;
    out.I_wp_yy            = I_wp_yy;
    out.CB_z_world         = CB_z_world;
    out.mass_buoyant_force = cfg.RHO_WATER * V_sub;
    out.mass_balance_error_abs = abs(M_total - out.mass_buoyant_force);

    if V_sub > 1e-10
        KM_world = CB_z_world + I_wp_yy / V_sub;
    else
        out.feasible = false;
        return;
    end
    out.KM_world = KM_world;

    GM = KM_world - out.CG_z_world;
    K33_hydro = cfg.RHO_WATER * cfg.G * Aw;
    % K55_hydro = M_total·g·GM (total mass at every site; equals rho_w·V_sub·g·GM at flotation balance).
    if GM > 0
        K55_hydro = M_total * cfg.G * GM;
    else
        K55_hydro = 0;
    end

    try
        % Pass the candidate's CG — same reasoning as evaluate_design_point (thin-shell).
        [A11, A33, A55, ~, ~] = ...
            mwecmass.bem.interpolate_at_draft( ...
                vs, cfg, out.CG_z_world);
    catch ME
        error('mwecmass:modular_precast:CoefficientInterpolationFailed', ...
              'mwecmass.bem.interpolate_at_draft failed at vs=%.4f m: %s.', ...
              vs, ME.message);
    end

    if K33_hydro > 1e-6
        T_heave = 2*pi * sqrt((M_total + A33) / K33_hydro);
    else
        T_heave = inf;
    end
    if K55_hydro > 1e-6 && (Iyy_about_cg + A55) > 0
        T_pitch = 2*pi * sqrt((Iyy_about_cg + A55) / K55_hydro);
    else
        T_pitch = inf;
    end

    out.feasible    = isfinite(GM) && isfinite(T_heave) && isfinite(T_pitch);
    out.GM          = GM;
    out.T_heave     = T_heave;
    out.T_pitch     = T_pitch;
    out.K33_hydro   = K33_hydro;
    out.K55_hydro   = K55_hydro;
    out.A11         = A11;
    out.A33         = A33;
    out.A55         = A55;
    % Suppress unused warnings (signature parallelism)
    t_offset_strip; is_solid_strip; %#ok<VUNUS>
end
