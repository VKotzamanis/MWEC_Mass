function out = evaluate_design_point(vs, t_steel, z_ballast, grids, ctx)  %#ok<INUSD>
%EVALUATE_DESIGN_POINT Evaluate thin-shell properties at (vs,t_steel,z_ballast).
%   Uses the supplied geometry grids without solving for draft. out contains
%   hydrostatics, mass/CG/inertia, GM, added mass, periods, and the absolute
%   buoyancy residual; SI units are used throughout. ctx supplies config and
%   hydrodynamic data needed by downstream property calculations.

    cfg = ctx.config;
    out = mwecmass.realise.empty_realised_properties();

    %% --- Split integrals (steel below z_ballast, jacket+air above) ---
    R = mwecmass.realise.thin_shell.integrate_split(grids, z_ballast);

    % V_ballast is the solid region below z_ballast; V_shell is the jacket above it.
    % V_steel and its integrals remain aggregate solid-region quantities.
    V_ballast  = R.V_below_outer;
    V_shell = R.V_above_jacket;
    V_steel = V_ballast + V_shell;
    V_air   = R.V_above_inner;

    int_z_x_A_ballast  = R.int_zA_below_outer;
    int_z_x_A_shell = R.int_zA_above_jacket;
    int_z_x_A_steel = int_z_x_A_ballast + int_z_x_A_shell;
    int_z_x_A_air   = R.int_zA_above_inner;

    int_x2_ballast  = R.int_Ix_below_outer;
    int_x2_shell = R.int_Ix_above_jacket;
    int_x2_steel = int_x2_ballast + int_x2_shell;
    int_x2_air   = R.int_Ix_above_inner;

    int_y2_ballast  = R.int_Iy_below_outer;
    int_y2_shell = R.int_Iy_above_jacket;
    int_y2_steel = int_y2_ballast + int_y2_shell;
    int_y2_air   = R.int_Iy_above_inner;

    int_z2_A_ballast  = R.int_z2A_below_outer;
    int_z2_A_shell = R.int_z2A_above_jacket;
    int_z2_A_steel = int_z2_A_ballast + int_z2_A_shell;
    int_z2_A_air   = R.int_z2A_above_inner;

    % Compute the per-region masses for the two-density material model.
    M_shell_region = ctx.rho_shell * V_shell;
    M_ballast_region  = ctx.rho_ballast  * V_ballast;

    % Keep the aggregate operation order for the equal-density default; the
    % difference term supplies the ballast-region correction for two densities.
    M_steel = ctx.rho_shell * V_steel + (ctx.rho_ballast - ctx.rho_shell) * V_ballast;
    M_air   = ctx.rho_air   * V_air;
    M_total = M_steel + M_air;

    out.V_steel = V_steel;
    out.V_air   = V_air;
    out.V_shell = V_shell;
    out.V_ballast  = V_ballast;
    out.M_steel = M_steel;
    out.M_air   = M_air;
    out.M_shell = M_shell_region;
    out.M_ballast  = M_ballast_region;
    out.M_total = M_total;

    if M_total <= 0
        out.feasible = false;
        return;
    end

    % z-centroids (body frame)
    if V_steel > 1e-12
        z_cg_steel = int_z_x_A_steel / V_steel;
    else
        z_cg_steel = 0.5 * (ctx.hull_z_min + z_ballast);
    end
    if V_air > 1e-12
        z_cg_air = int_z_x_A_air / V_air;
    else
        z_cg_air = 0.5 * (z_ballast + ctx.hull_z_max);
    end
    % Per-region centroids, required by the two-density model.
    % Degenerate-volume fallbacks
    % mirror the z_cg_steel/z_cg_air pattern immediately above (same body-frame ranges:
    % [hull_z_min, z_ballast] for ballast, [z_ballast, hull_z_max] for shell/air).
    if V_ballast > 1e-12
        z_cg_ballast = int_z_x_A_ballast / V_ballast;
    else
        z_cg_ballast = 0.5 * (ctx.hull_z_min + z_ballast);
    end
    if V_shell > 1e-12
        z_cg_shell = int_z_x_A_shell / V_shell;
    else
        % Degenerate-volume convention (fires only when V_shell <= 1e-12): a region's centroid
        % is undefined when its volume is zero, so a finite placeholder is assigned rather than
        % left NaN/Inf. The implemented finite value here is the midpoint of the shell's
        % nominal z-range, 0.5*(z_ballast + hull_z_max) -- not z_ballast alone. This is inert for
        % every mass total: int_z_x_A_shell is itself zero when V_shell is (the numerator this
        % centroid would otherwise divide), so CG_z_body's solid_first_moment term picks up no
        % contribution from this branch regardless of what z_cg_shell is set to here.
        z_cg_shell = 0.5 * (z_ballast + ctx.hull_z_max);
    end
    % Form the first moment directly by region; equal densities use the aggregate form.
    if ctx.rho_ballast == ctx.rho_shell
        solid_first_moment = M_steel * z_cg_steel;
    else
        solid_first_moment = ctx.rho_ballast * int_z_x_A_ballast + ctx.rho_shell * int_z_x_A_shell;
    end
    CG_z_body = (solid_first_moment + M_air * z_cg_air) / M_total;

    % Iyy (pitch axis, body origin):  ∫∫∫ ρ (x² + z²) dV
    % The single-density expression (ctx.rho_shell * combined-region sum) plus a
    % difference term over the ballast-only sub-integral, exactly zero at rho_ballast == rho_shell
    % (see the M_steel derivation above for the same discipline and its algebraic proof).
    Iyy_total_origin = ctx.rho_shell * (int_x2_steel + int_z2_A_steel) + ...
                       (ctx.rho_ballast - ctx.rho_shell) * (int_x2_ballast + int_z2_A_ballast) + ...
                       ctx.rho_air   * (int_x2_air   + int_z2_A_air);
    % Ixx (roll axis):                  ∫∫∫ ρ (y² + z²) dV
    Ixx_total_origin = ctx.rho_shell * (int_y2_steel + int_z2_A_steel) + ...
                       (ctx.rho_ballast - ctx.rho_shell) * (int_y2_ballast + int_z2_A_ballast) + ...
                       ctx.rho_air   * (int_y2_air   + int_z2_A_air);
    % Izz (yaw):                        ∫∫∫ ρ (x² + y²) dV
    Izz_total_origin = ctx.rho_shell * (int_x2_steel + int_y2_steel) + ...
                       (ctx.rho_ballast - ctx.rho_shell) * (int_x2_ballast + int_y2_ballast) + ...
                       ctx.rho_air   * (int_x2_air   + int_y2_air);

    Iyy_about_cg = max(0, Iyy_total_origin - M_total * CG_z_body^2);
    Ixx_about_cg = max(0, Ixx_total_origin - M_total * CG_z_body^2);
    Izz_about_cg = max(0, Izz_total_origin);

    % Export z_cg_steel as the geometric centroid for equal densities and as
    % the aggregate solid-mass centroid when the two densities differ.
    if ctx.rho_ballast == ctx.rho_shell
        out.z_cg_steel = z_cg_steel;
    else
        out.z_cg_steel = solid_first_moment / M_steel;
    end
    out.z_cg_ballast        = z_cg_ballast;
    out.z_cg_shell       = z_cg_shell;
    out.z_cg_air         = z_cg_air;
    out.CG_z_body        = CG_z_body;
    out.Iyy_total_origin = Iyy_total_origin;
    out.Iyy_about_cg     = Iyy_about_cg;
    out.Ixx_total_origin = Ixx_total_origin;
    out.Ixx_about_cg     = Ixx_about_cg;
    out.Izz_total_origin = Izz_total_origin;
    out.Izz_about_cg     = Izz_about_cg;

    % Draft / vs from input (vs is now a DV, not derived)
    out.draft          = abs(ctx.hull_z_min + vs);
    out.vertical_shift = vs;
    out.CG_z_world     = CG_z_body + vs;

    %% --- Waterplane / submerged hydrostatics at the given vs ---
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

    % KM (world frame)
    if V_sub > 1e-10
        KM_world = CB_z_world + I_wp_yy / V_sub;
    else
        out.feasible = false;
        return;
    end
    out.KM_world = KM_world;

    GM = KM_world - out.CG_z_world;

    %% --- Stiffness, added mass, natural periods ---
    K33_hydro = cfg.RHO_WATER * cfg.G * Aw;
    % K55_hydro = M_total·g·GM (total mass at every site; equals rho_w·V_sub·g·GM at flotation balance).
    if GM > 0
        K55_hydro = M_total * cfg.G * GM;
    else
        K55_hydro = 0;
    end

    try
        % Pass the candidate CG so A55 and surge-pitch cross terms use the
        % correct reference axis during optimisation.
        [A11, A33, A55, ~, ~] = ...
            mwecmass.bem.interpolate_at_draft( ...
                vs, cfg, out.CG_z_world);
    catch ME
        error('mwecmass:thin_shell:CoefficientInterpolationFailed', ...
              ['mwecmass.bem.interpolate_at_draft failed at vs=%.4f m: %s. ' ...
               'The hydro cache must cover the realisation draft range.'], ...
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
end
