function props = properties_3d(x, config)
%PROPERTIES_3D Full 3D WEC physical properties via strip integration on parametric B-spline cross-sections.
% See docs/METHODS_ENGINE.md#coupled-modal-energy-share.
% Inputs: x = [vertical_shift (m); ballast densities (kg/m^3)]; config from build_config,
% including precomputed strip geometry and WAMIT hydrodynamic data.
% Outputs include volume/waterplane geometry, mass and inertia, hydrostatic/added-mass matrices,
% and periods [s]. Coupled heave/pitch periods are named by kinetic-energy share; per-axis values
% are retained. Sign/frame conventions: world z=0 is the free surface, z up, and positive
% vertical_shift moves the hull up; draft is positive downward and GM_L = KM - CG_z.
% The free-floating model has K_pto = 0 and K_total = K_hydro.
    try
        %% ===== 1. UNPACK INPUTS =====
        props.vertical_shift = x(1);
        x = x(:);
        densities_at_nodes = x(2:end);

        parser     = config.ms2_model;
        hull_z_min = config.hull_z_min;
        hull_z_max = config.hull_z_max;

        % Waterline in body frame: z_wl = -vertical_shift
        %   vertical_shift > 0 -> hull moves up -> waterline is below hull centre
        %   vertical_shift < 0 -> hull moves down -> more submerged
        z_wl = -props.vertical_shift;

        % True draft: distance from waterline (z=0 world) to keel.
        % keel_z_world = hull_z_min + vertical_shift; draft = -keel_z_world (positive)
        props.draft = abs(hull_z_min + props.vertical_shift);

        %% ===== 2-3. SUBMERGED VOLUME, CB, WATERPLANE =====
        %  Fast path (preferred): interpolate V_sub, CB_z, Aw, I_wp from precomputed tables built
        %  at config time. Bypasses compute_submerged entirely -- correct because these quantities
        %  depend only on z_wl and hull geometry, not density, and the tables share config's z-grid.
        %  Slow path (fallback): calls compute_submerged directly when V_sub_table is absent
        %  (geometry-only config or a config without tables) -- same result.
        if z_wl > hull_z_min
            z_sub_top = min(z_wl, hull_z_max);

            if isfield(config, 'V_sub_table') && ~isempty(config.V_sub_table)
                % Fast path: all four quantities from precomputed tables
                props.Aw      = interp1(config.Aw_table_z, config.Aw_table,      z_sub_top, 'linear', 0);
                props.I_wp_xx = interp1(config.Aw_table_z, config.I_wp_xx_table, z_sub_top, 'linear', 0);
                props.I_wp_yy = interp1(config.Aw_table_z, config.I_wp_yy_table, z_sub_top, 'linear', 0);
                props.V_sub   = max(0, interp1(config.Aw_table_z, config.V_sub_table, z_sub_top, 'linear', 0));
                CB_z_body     = interp1(config.Aw_table_z, config.CB_z_table,    z_sub_top, 'linear', 'extrap');
                props.CB      = [0, 0, CB_z_body + props.vertical_shift];

                % Wetted surface area (hull sides only, waterplane excluded)
                if isfield(config, 'S_wet_table') && ~isempty(config.S_wet_table)
                    props.A_sub = max(0, interp1(config.Aw_table_z, config.S_wet_table, ...
                                                  z_sub_top, 'linear', 0));
                else
                    % Table absent: compute directly from the geometry.
                    props.A_sub = mwecmass.hydrostatics.compute_wetted_surface_area( ...
                                      parser, z_sub_top);
                end

            else
                % Fallback: full surface integral when tables are unavailable.
                sub_opts = struct('n_quad', 16, 'verbose', false);

                if isfield(config, 'Aw_table_z') && ~isempty(config.Aw_table_z)
                    sub_opts.Aw_override = interp1( ...
                        config.Aw_table_z, config.Aw_table, ...
                        z_sub_top, 'linear', 0);
                    sub_opts.I_wp_xx_override = interp1( ...
                        config.Aw_table_z, config.I_wp_xx_table, ...
                        z_sub_top, 'linear', 0);
                    sub_opts.I_wp_yy_override = interp1( ...
                        config.Aw_table_z, config.I_wp_yy_table, ...
                        z_sub_top, 'linear', 0);
                end

                hp_sub        = mwecmass.hydrostatics.compute_submerged(parser, z_sub_top, sub_opts);
                props.V_sub   = hp_sub.V_sub;
                props.Aw      = hp_sub.Aw;
                props.I_wp_xx = hp_sub.I_wp_xx;
                props.I_wp_yy = hp_sub.I_wp_yy;
                props.CB      = [0, 0, hp_sub.CB(3) + props.vertical_shift];
                % Wetted surface area: direct parametric integration (no table)
                props.A_sub   = mwecmass.hydrostatics.compute_wetted_surface_area( ...
                                    parser, z_sub_top);
            end
        else
            props.V_sub   = 0;
            props.Aw      = 0;
            props.I_wp_xx = 0;
            props.I_wp_yy = 0;
            props.CB      = [0, 0, 0];
            props.A_sub   = 0;
        end

        % Metacentric radius and KM
        if props.V_sub > 1e-6
            % I_wp_yy (waterplane_properties.m) is the second moment of area about the ORIGIN
            % (x = 0), not about the waterplane centre of flotation; no parallel-axis shift is
            % applied here. This equals the centre-of-flotation value only while the waterplane
            % centroid lies on x = 0 (true for C1 by axisymmetry); for a non-axisymmetric hull
            % it is in error by A_w * x_F^2.
            BM_L = props.I_wp_yy / props.V_sub;
        else
            BM_L = 0;
        end

        props.mass_buoyant_force = props.V_sub * config.RHO_WATER;
        props.KM = props.CB(3) + BM_L;

        %% ===== 4. MASS PROPERTIES - PRECOMPUTED STRIPS =====
        %  Strip volumes, centroids, and inertias were precomputed at config time using
        %  mwecmass.hydrostatics.compute_strip (parametric divergence theorem). Mass integration is now
        %  just dot products: mass_i = V_i * rho_i.
        N = length(densities_at_nodes);

        mass_total   = 0;
        cg_z_num     = 0;
        Iyy_total    = 0;
        Ixx_total    = 0;
        Izz_total    = 0;

        for i = 1:N
            V_strip = config.strip_V(i);
            if V_strip < 1e-12
                continue;
            end

            rho_i = max(config.ballast_density_bounds(1), ...
                    min(config.ballast_density_bounds(2), densities_at_nodes(i)));

            strip_mass = V_strip * rho_i;

            mass_total = mass_total + strip_mass;
            cg_z_num   = cg_z_num + strip_mass * config.strip_CB_z(i);
            Iyy_total  = Iyy_total + rho_i * config.strip_Iyy(i);
            Ixx_total  = Ixx_total + rho_i * config.strip_Ixx(i);
            Izz_total  = Izz_total + rho_i * config.strip_Izz(i);
        end

        props.mass_total = mass_total;

        % CG in body frame, then shift to world frame
        if mass_total > 1e-6
            cg_z_body = cg_z_num / mass_total;
        else
            cg_z_body = 0.5 * (hull_z_min + hull_z_max);
        end
        props.CG_total = [0, 0, cg_z_body + props.vertical_shift];

        %% ===== 5. INERTIA TENSOR =====
        %  Parallel axis theorem: I_cg = I_origin - M*d^2. For Iyy/Ixx, d = cg_z. For Izz, d = 0
        %  (CG is on the z-axis for a symmetric hull, so Izz already excludes a z^2 term).
        Iyy_about_cg = Iyy_total - mass_total * cg_z_body^2;
        Iyy_about_cg = max(0, Iyy_about_cg);  % guard against numerical noise

        Ixx_about_cg = Ixx_total - mass_total * cg_z_body^2;
        Ixx_about_cg = max(0, Ixx_about_cg);

        Izz_about_cg = max(0, Izz_total);

        % AutoCAD validation override (if provided)
        if ~isempty(config.autocad_Iyy)
            discrepancy_pct = abs(Iyy_about_cg - config.autocad_Iyy) / config.autocad_Iyy * 100;
            if discrepancy_pct > config.autocad_discrepancy_pct
                fprintf('  WARNING: Iyy discrepancy: %.1f%% (MATLAB vs AutoCAD)\n', discrepancy_pct);
                Iyy_about_cg = config.autocad_Iyy;
                if ~isempty(config.autocad_Ixx)
                    Ixx_about_cg = config.autocad_Ixx;
                end
                if ~isempty(config.autocad_Izz)
                    Izz_about_cg = config.autocad_Izz;
                end
            end
        end

        props.Inertia_Tensor = diag([Ixx_about_cg, Iyy_about_cg, Izz_about_cg]);
        props.Ixx = Ixx_about_cg;
        props.Iyy = Iyy_about_cg;
        props.Izz = Izz_about_cg;

        %% ===== 6. STABILITY =====
        props.GM_L = props.KM - props.CG_total(3);

        %% ===== 7. MASS BALANCE =====
        props.mass_discrepancy = props.mass_total - props.mass_buoyant_force;

        %% ===== 8. HYDROSTATIC STIFFNESS =====
        K33_hydro = config.RHO_WATER * config.G * props.Aw;

        % K55_hydro = M_total·g·GM (total mass at every site; equals rho_w·V_sub·g·GM at flotation balance).
        if props.GM_L > 0
            K55_hydro = props.mass_total * config.G * props.GM_L;   % N·m/rad, hydrostatic pitch restoring stiffness
        else
            K55_hydro = 0;
        end

        props.K_hydro = diag([0, K33_hydro, K55_hydro]);

        %% ===== 9. ADDED MASS =====
        %  Interpolate from the hydro table and apply the delta congruence transform from
        %  HAMS-CG (uniform-density CG stored in the cache) to the actual optimised CG in one
        %  call; props.A_full/A55 reference props.CG_total(3) directly, no further downstream
        %  retransform.
        [props.A11, props.A33, props.A55, props.A_full, props.B_full] = ...
            mwecmass.bem.interpolate_at_draft( ...
                props.vertical_shift, config, props.CG_total(3));

        %% ===== 10. TOTAL STIFFNESS (no PTO block -- see file header) =====
        %  Free-floating body: no PTO stiffness contribution. K_pto is the fixed zero matrix;
        %  K_total equals K_hydro (no surge-pitch cross term).
        K11_total = 0;
        K33_total = K33_hydro;
        K55_total = K55_hydro;
        props.K_pto = zeros(3, 3);
        props.K_total = [K11_total, 0,       0;
                         0,         K33_total, 0;
                         0,         0,       K55_total];

        %% ===== 11. NATURAL PERIODS (UNCOUPLED) =====
        M11_virtual = props.mass_total + props.A11;
        M33_virtual = props.mass_total + props.A33;
        M55_virtual = Iyy_about_cg + props.A55;

        if K11_total > 1e-6
            props.periods.surge = 2*pi * sqrt(M11_virtual / K11_total);
        else
            props.periods.surge = inf;
        end

        if K33_total > 1e-6
            props.periods.heave = 2*pi * sqrt(M33_virtual / K33_total);
        else
            props.periods.heave = inf;
        end

        if K55_total > 1e-6
            props.periods.pitch = 2*pi * sqrt(M55_virtual / K55_total);
        else
            props.periods.pitch = inf;
        end

        %% ===== 12. COUPLED EIGENPROBLEM; HEAVE AND PITCH NAMED BY SHARE =====
        %  props.periods.heave/.pitch are the eigen-periods of the coupled 3-DOF
        %  (1 surge, 2 heave, 3 pitch) undamped problem, each named by the DOF holding the
        %  largest share of that mode's kinetic energy. The per-axis values are
        %  kept beside them as .heave_uncoupled/.pitch_uncoupled; the coupled problem
        %  is solved on every call.
        M_phys = diag([props.mass_total, props.mass_total, Iyy_about_cg]);   % [kg, kg, kg m^2]
        M_total = M_phys + props.A_full;   % [kg, kg m, kg m^2] physical + omega->inf added mass
        K_total_full = props.K_total;      % [N/m, N/m, N m/rad] hydrostatic only, K_pto = 0

        props.periods.heave_uncoupled = props.periods.heave;   % [s] per-axis value
        props.periods.pitch_uncoupled = props.periods.pitch;   % [s] per-axis value

        [T_coupled, share, Phi] = ...
            mwecmass.hydrostatics.coupled_periods_by_share(M_total, K_total_full);

        props.periods.heave = T_coupled.heave;   % [s] heave-dominated coupled mode
        props.periods.pitch = T_coupled.pitch;   % [s] pitch-dominated coupled mode
        props.participation_factors = share;     % [%] kinetic-energy share, rows surge/heave/pitch
        props.coupled_periods = [T_coupled.surge; T_coupled.heave; T_coupled.pitch];   % [s]
        props.coupled_modes = Phi;               % columns ordered surge, heave, pitch
        props.surge_per_pitch = Phi(1,3) / Phi(3,3);   % [m/rad] surge per unit pitch in the pitch mode

        %% ===== 13. COMPONENT BREAKDOWN =====
        props.components = struct('density', num2cell(densities_at_nodes), ...
                                  'z_level', num2cell(config.density_nodes_z + props.vertical_shift));

        %% ===== 14. MASS MATRICES (6x6) =====
        props.MassMatrix_CG = mwecmass.hydrostatics.build_mass_matrix(...
            props.mass_total, [0, 0, 0], props.Inertia_Tensor);
        props.MassMatrix_Origin = mwecmass.hydrostatics.build_mass_matrix(...
            props.mass_total, props.CG_total, props.Inertia_Tensor);

        %% ===== 15. STORE DENSITY PROFILE =====
        props.densities_at_nodes = densities_at_nodes;

    catch ME
        warning('mwecmass:hydrostatics:ComputationFailed', ...
                'Property calculation failed: %s. Returning defaults.', ME.message);

        props.vertical_shift = x(1);
        props.draft = 0;
        props.mass_total = 0;
        props.CG_total = [0, 0, 0];
        props.CB = [0, 0, 0];
        props.GM_L = 0;
        props.V_sub = 0;
        props.Aw = 0;
        props.I_wp_xx = 0;
        props.I_wp_yy = 0;
        props.KM = 0;
        props.mass_buoyant_force = 0;
        props.Inertia_Tensor = zeros(3,3);
        props.Ixx = 0;
        props.Iyy = 0;
        props.Izz = 0;
        props.A11 = 0;
        props.A33 = 0;
        props.A55 = 0;
        props.A_full = zeros(3, 3);
        props.B_full = zeros(3, 3);
        props.K_hydro = zeros(3, 3);
        props.K_pto = zeros(3, 3);
        props.K_total = zeros(3, 3);
        props.periods.surge = inf;
        props.periods.heave = inf;
        props.periods.pitch = inf;
        props.coupled_periods = [inf; inf; inf];
        props.coupled_modes = eye(3);
        props.participation_factors = 100 * eye(3);
        props.mass_discrepancy = 0;
        props.components = struct('density', {}, 'z_level', {});
        props.MassMatrix_CG = eye(6);
        props.MassMatrix_Origin = eye(6);
        props.densities_at_nodes = x(2:end);
        props.A_sub = 0;
    end
end
