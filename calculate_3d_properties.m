function props = calculate_3d_properties(x, config)
% CALCULATE_3D_PROPERTIES  Compute full 3D WEC physical properties via
%   strip integration on parametric B-spline cross-sections.
%
%   REFACTORED (v4.0): No mesh objects.  All geometry queries go through
%   WEC_Core_Functions.evaluateCrossSectionMS2 which evaluates the parent
%   B-spline curves of each parametric surface.  The vertical_shift (x(1))
%   offsets the z-coordinate system — the waterline in body frame is at
%   z_wl = -vertical_shift, and "submerged" means z < z_wl.
%
%   COMPUTATION SEQUENCE
%     1. Unpack vertical_shift and density vector.
%     2. Waterplane properties: cross-section at z_wl → Aw, I_wp_yy, BM.
%     3. Submerged volume & CB: strip integration from z_min to z_wl.
%     4. Mass, CG, Iyy: strip integration over FULL hull with density.
%     5. GM = KM − CG_z  (KM = KB + BM, in world frame).
%     6. Stiffness, added mass, periods (unchanged from v3.0).
%     7. Coupled eigenvalue analysis (unchanged from v3.0).
%
%   INPUTS
%     x      : [1+N × 1] optimisation vector [vertical_shift; rho_1..rho_N]
%     config : struct from WEC_Configuration_Builder
%              Must contain: ms2_model, density_nodes_z, hull_z_min,
%              hull_z_max, strip_edges (or num_ballast_sections),
%              ballast_density_bounds, RHO_WATER, G, and WAMIT data.
%
%   OUTPUT
%     props : struct — see field list at end of function
%
%   See also: WEC_Core_Functions.compute_strip_bspline,
%             WEC_Core_Functions.evaluateCrossSectionMS2,
%             calculate_2d_properties
%
%   Author:  WEC Optimisation Team
%   Version: 4.0 — Spline-based strip integration (no mesh)

    try
        %% ===== 1. UNPACK INPUTS =====
        props.vertical_shift = x(1);
        x = x(:);
        densities_at_nodes = x(2:end);

        parser     = config.ms2_model;
        hull_z_min = config.hull_z_min;
        hull_z_max = config.hull_z_max;

        % Waterline in body frame: z_wl = -vertical_shift
        %   vertical_shift > 0 → hull moves up → waterline is below hull center
        %   vertical_shift < 0 → hull moves down → more submerged
        z_wl = -props.vertical_shift;

        % True draft: distance from waterline (z=0 world) to keel
        % In world frame: keel_z_world = hull_z_min + vertical_shift
        % Draft = distance from z=0 down to keel = -keel_z_world (positive)
        props.draft = abs(hull_z_min + props.vertical_shift);

        %% ===== 2–3. SUBMERGED VOLUME, CB, WATERPLANE (v4.4) =====
        %  Fast path (preferred): interpolate V_sub, CB_z, Aw, I_wp from
        %  precomputed tables built at config time (§2c of
        %  WEC_Configuration_Builder).  Bypasses compute_submerged entirely.
        %
        %  Slow path (fallback): calls compute_submerged when V_sub_table
        %  is absent (geometry-only config, or old config struct).  Identical
        %  behaviour to v4.3 — no regression for existing callers.
        %
        %  WHY the fast path is correct:
        %    V_sub, CB_z, Aw, I_wp depend only on z_wl and hull geometry,
        %    not on density.  Tables are built on the same z-grid as
        %    config.Aw_table_z using the same compute_submerged call.
        %    interp1('linear') matches the resolution of the adaptive table.

        if z_wl > hull_z_min
            z_sub_top = min(z_wl, hull_z_max);

            if isfield(config, 'V_sub_table') && ~isempty(config.V_sub_table)
                % ── Fast path: all four quantities from precomputed tables ──
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
                    % Table absent (old config or geometry-only): compute directly
                    props.A_sub = WEC_HydroProperties.compute_wetted_surface_area( ...
                                      parser, z_sub_top);
                end

            else
                % ── Slow path: full surface integral (backward-compatible) ──
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

                hp_sub        = WEC_HydroProperties.compute_submerged(parser, z_sub_top, sub_opts);
                props.V_sub   = hp_sub.V_sub;
                props.Aw      = hp_sub.Aw;
                props.I_wp_xx = hp_sub.I_wp_xx;
                props.I_wp_yy = hp_sub.I_wp_yy;
                props.CB      = [0, 0, hp_sub.CB(3) + props.vertical_shift];
                % Wetted surface area: direct parametric integration (no table)
                props.A_sub   = WEC_HydroProperties.compute_wetted_surface_area( ...
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
            BM_L = props.I_wp_yy / props.V_sub;
        else
            BM_L = 0;
        end

        props.mass_buoyant_force = props.V_sub * config.RHO_WATER;
        props.KM = props.CB(3) + BM_L;

        %% ===== 4. MASS PROPERTIES — PRECOMPUTED STRIPS (v4.3) =====
        %  Strip volumes, centroids, and inertias were precomputed at config
        %  time (§3f of WEC_Configuration_Builder) using
        %  WEC_HydroProperties.compute_strip (parametric divergence theorem).
        %  Mass integration is now just dot products: mass_i = V_i × ρ_i.

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
        %  Parallel axis theorem: I_cg = I_origin − M × d²
        %  For Iyy (pitch) and Ixx (roll): d = cg_z (distance along z)
        %  For Izz (yaw): d = 0 (CG is on the z-axis for symmetric hulls)

        Iyy_about_cg = Iyy_total - mass_total * cg_z_body^2;
        Iyy_about_cg = max(0, Iyy_about_cg);  % guard against numerical noise

        Ixx_about_cg = Ixx_total - mass_total * cg_z_body^2;
        Ixx_about_cg = max(0, Ixx_about_cg);

        % Izz_total is already ∫ ρ(z)[∫∫(x²+y²)dA]dz — no z² term,
        % so for a hull with x_cg = y_cg = 0: Izz_cg = Izz_origin.
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

        if props.GM_L > 0
            K55_hydro = props.mass_buoyant_force * config.G * props.GM_L;
            % equivalently: config.RHO_WATER * config.G * props.V_sub * props.GM_L
        else
            K55_hydro = 0;
        end

        props.K_hydro = diag([0, K33_hydro, K55_hydro]);

        %% ===== 9. ADDED MASS =====
        %  Interpolate from the hydro table and apply the delta congruence
        %  transform from HAMS-CG (= uniform-density CG stored in the cache)
        %  to the ACTUAL optimised CG in one call.  See
        %  WEC_Core_Functions.interpolate_wamit_added_mass for the math.
        %  Centralising the transform here replaces the former §9-9b two-step
        %  pattern; props.A_full and props.A55 now reference props.CG_total(3)
        %  directly with no further retransform downstream.
        [props.A11, props.A33, props.A55, ~, props.A_full, props.B_full] = ...
            WEC_Core_Functions.interpolate_wamit_added_mass( ...
                props.vertical_shift, config, props.CG_total(3));

        %% ===== 10. PTO STIFFNESS =====
        %  From CONVENTION_BRIDGE_MODAL_ANALYSIS.md §5.4:
        %    K_PTO = 2K × [cos²α, 0, −ℓcosα; 0, sin²α, 0; −ℓcosα, 0, ℓ²]
        %    ℓ = h_t cosα + b_t sinα   (tether eccentricity)
        %
        %  α = angle from horizontal. Vertical tether: α = 90°.
        %    cos²(90°) = 0 → no surge stiffness from vertical tether. ✓
        %    sin²(90°) = 1 → full heave stiffness from vertical tether. ✓
        if config.enable_pto_effects == 1
            K_pto_base = K33_hydro;  % [N/m] per-tether axial stiffness placeholder
            a_rad = deg2rad(config.pto_angle_deg);

            % Tether eccentricity (Convention Bridge §5.4)
            ell = config.ht * cos(a_rad) + config.bt * sin(a_rad);

            K11_pto = 2 * K_pto_base * cos(a_rad)^2;    % surge
            K33_pto = 2 * K_pto_base * sin(a_rad)^2;    % heave
            K55_pto = 2 * K_pto_base * ell^2;            % pitch
            K13_pto = -2 * K_pto_base * ell * cos(a_rad); % surge-pitch coupling
        else
            K11_pto = 0;
            K33_pto = 0;
            K55_pto = 0;
            K13_pto = 0;
        end

        props.K_pto = [K11_pto, 0, K13_pto;
                       0,       K33_pto, 0;
                       K13_pto, 0, K55_pto];

        %% ===== 11. TOTAL STIFFNESS =====
        K11_total = K11_pto;
        K33_total = K33_hydro + K33_pto;
        K55_total = K55_hydro + K55_pto;
        props.K_total = [K11_total, 0,       K13_pto;
                         0,         K33_total, 0;
                         K13_pto,   0,       K55_total];

        %% ===== 12. NATURAL PERIODS (UNCOUPLED) =====
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

        %% ===== 13. COUPLED EIGENVALUE ANALYSIS =====
        try
            M_phys = diag([props.mass_total, props.mass_total, Iyy_about_cg]);
            M_total = M_phys + props.A_full;
            K_total_full = props.K_total;

            if det(M_total) > 1e-12 && det(K_total_full) > 1e-12
                [V_eig, D_eig] = eig(K_total_full, M_total);
                omega_sq = diag(D_eig);

                valid_modes = omega_sq > 1e-6;
                if any(valid_modes)
                    omega_n = sqrt(omega_sq(valid_modes));
                    T_n = 2*pi ./ omega_n;

                    [T_n_sorted, sort_idx] = sort(T_n, 'descend');
                    props.coupled_periods = T_n_sorted;

                    V_valid = V_eig(:, valid_modes);
                    props.coupled_modes = V_valid(:, sort_idx);
                    props.participation_factors = compute_participation_factors(...
                        props.coupled_modes, M_total);
                end
            else
                props.coupled_periods = [props.periods.surge; props.periods.heave; props.periods.pitch];
                props.coupled_modes = eye(3);
                props.participation_factors = 100 * eye(3);
            end
        catch
            props.coupled_periods = [props.periods.surge; props.periods.heave; props.periods.pitch];
            props.coupled_modes = eye(3);
            props.participation_factors = 100 * eye(3);
        end

        %% ===== 14. COMPONENT BREAKDOWN =====
        props.components = struct('density', num2cell(densities_at_nodes), ...
                                  'z_level', num2cell(config.density_nodes_z + props.vertical_shift));

        %% ===== 15. MASS MATRICES (6x6) =====
        props.MassMatrix_CG = WEC_Core_Functions.calculate6x6MassMatrix(...
            props.mass_total, [0, 0, 0], props.Inertia_Tensor);
        props.MassMatrix_Origin = WEC_Core_Functions.calculate6x6MassMatrix(...
            props.mass_total, props.CG_total, props.Inertia_Tensor);

        %% ===== 16. STORE DENSITY PROFILE =====
        props.densities_at_nodes = densities_at_nodes;

    catch ME
        warning('calculate_3d_properties:ComputationFailed', ...
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


%% ===== HELPER FUNCTION =====
function PF = compute_participation_factors(modes, M)
% COMPUTE_PARTICIPATION_FACTORS  DOF contribution to each coupled mode.
%
%   PF(i,j) = φ_i² M_ii / (φ' M φ) × 100%

try
    n_modes = size(modes, 2);
    PF = zeros(3, n_modes);

    for j = 1:n_modes
        phi = modes(:, j);
        modal_mass = phi' * M * phi;

        if abs(modal_mass) > 1e-12
            for i = 1:3
                PF(i, j) = (phi(i)^2 * M(i,i)) / modal_mass * 100;
            end
        end
    end

    for j = 1:n_modes
        total = sum(PF(:, j));
        if abs(total - 100) > 5
            PF(:, j) = PF(:, j) * 100 / total;
        end
    end

catch
    PF = zeros(3, size(modes, 2));
end

end