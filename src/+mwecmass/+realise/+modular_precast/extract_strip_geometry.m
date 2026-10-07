function cstr = extract_strip_geometry(config, t_offset_arg, z_ballast, ...
                                       rho_UHPC, rho_air, solve_data)
%EXTRACT_STRIP_GEOMETRY Extract per-strip UHPC jacket/void contours and volumes.
% t_offset_arg may be scalar or per-strip; Inf and is_solid_strip mark fully solid strips. z_ballast separates
% solid lower regions from annular jacket/void regions. cstr inherits solve_data and adds visualisation data.
% See docs/METHODS_ENGINE.md#realise-modular-precast
    % Inherit all global solve fields (M_total, z_ballast, t_uhpc, etc.)
    cstr = solve_data;

    %% UNPACK CONFIG
    ms2_model   = config.ms2_model;
    n_sub       = config.constructability_n_sub;
    N           = length(config.density_nodes_z);
    Z_max       = config.hull_z_max;
    Z_min       = config.hull_z_min;
    t_min       = config.constructability_t_min;
    wall_height = config.constructability_wall_height;

    if isfield(config, 'wall_position') && strcmp(config.wall_position, 'bottom')
        wall_z_bottom = Z_min + wall_height;  %#ok<NASGU>
    else
        wall_z_bottom = Z_max - wall_height;  %#ok<NASGU>
    end

    if isfield(config, 'boundary_cache') && ~isempty(config.boundary_cache)
        b_cache = config.boundary_cache;
    else
        b_cache = mwecmass.geometry.precompute_boundary_cache(ms2_model, 100);
    end

    has_P_table = isfield(config, 'P_table') && ~isempty(config.P_table) && ...
                  length(config.P_table) == length(config.Aw_table_z);

    %% STRIP BOUNDARIES
    density_nodes_z = config.density_nodes_z(:);
    if ~isempty(config.strip_edges)
        strip_edges = config.strip_edges(:);
    else
        strip_edges = zeros(N + 1, 1);
        strip_edges(1) = Z_min;
        for ii = 2:N
            strip_edges(ii) = 0.5 * (density_nodes_z(ii-1) + density_nodes_z(ii));
        end
        strip_edges(N + 1) = Z_max;
    end

    %% PER-STRIP UNIFIED LOOP
    %
    % Classify each z-sample as solid (wall or below z_ballast) or annular; trapz
    % integrates strips that straddle z_ballast without a separate split.

    if isfield(config, 'wall_position') && strcmp(config.wall_position, 'bottom')
        wall_strip_idx = 1;    % _180: structural wall at hull bottom
    else
        wall_strip_idx = N;    % default: structural wall at hull top
    end

    % ── Resolve per-strip thickness vector (NEW Phase 1b support) ──
    if isscalar(t_offset_arg)
        t_offset_strip = t_offset_arg * ones(N, 1);
    else
        t_offset_strip = t_offset_arg(:);
        if length(t_offset_strip) ~= N
            error('mwecmass:modular_precast:BadTOffsetLength', ...
                  ['t_offset_arg length %d does not match strip count %d.'], ...
                  length(t_offset_strip), N);  %#ok<NBRAK2>
        end
    end
    % Strips with Inf t_offset are fully solid
    t_strip_is_solid_input = ~isfinite(t_offset_strip);

    % Phase 1b is_solid_strip flag (overrides for upgraded strips)
    if isfield(solve_data, 'is_solid_strip') && ~isempty(solve_data.is_solid_strip)
        t_strip_is_solid_input = t_strip_is_solid_input | solve_data.is_solid_strip(:);
    end
    % Wall is always solid (defensive)
    t_strip_is_solid_input(wall_strip_idx) = true;

    strip_z_lo        = zeros(N, 1);
    strip_z_hi        = zeros(N, 1);
    strip_V_total     = zeros(N, 1);
    strip_V_UHPC      = zeros(N, 1);
    strip_V_void      = zeros(N, 1);
    strip_mass_UHPC   = zeros(N, 1);
    strip_mass_void   = zeros(N, 1);
    strip_mass_total  = zeros(N, 1);
    strip_t_actual    = inf(N, 1);    % Inf = solid; finite = annular jacket
    strip_r_min       = inf(N, 1);
    strip_z_cg        = zeros(N, 1);
    strip_Iyy_UHPC    = zeros(N, 1);
    strip_Iyy_void    = zeros(N, 1);
    strip_is_wall     = false(N, 1);  % designated structural wall only
    strip_is_solid    = false(N, 1);  % wall OR entirely below z_ballast
    strip_is_feasible = true(N, 1);
    strip_scale       = zeros(N, 1);
    strip_rho_eff     = zeros(N, 1);

    contours_outer = cell(N, 1);
    contours_inner = cell(N, 1);

    for i = 1:N
        z_lo = strip_edges(i);
        z_hi = strip_edges(i + 1);
        strip_z_lo(i) = z_lo;
        strip_z_hi(i) = z_hi;

        n_sub_i = max(n_sub, ceil((z_hi - z_lo) / 0.01));
        z_i = linspace(z_lo, z_hi, n_sub_i)';
        z_i(end) = min(z_i(end), z_hi - 1e-4);

        % Use the same options as the driver strip integrator for consistent geometry.
        strip_phys = mwecmass.hydrostatics.compute_strip( ...
            ms2_model, z_lo, z_hi, ...
            struct('n_quad',        16, ...
                   'Aw_table_z',    config.Aw_table_z, ...
                   'Aw_table',      config.Aw_table, ...
                   'I_wp_xx_table', config.I_wp_xx_table, ...
                   'I_wp_yy_table', config.I_wp_yy_table));
        V_i   = strip_phys.V;      % [m^3] strip enclosed volume
        CBz_i = strip_phys.CB_z;   % [m]   strip centroid elevation (body frame)

        % ── Outer contours at each z-sample ─────────────────
        conts_i   = cell(n_sub_i, 1);
        A_outer_k = zeros(n_sub_i, 1);
        for k_c = 1:n_sub_i
            wl_k = mwecmass.geometry.extract_isocurve_at_z( ...
                       ms2_model, z_i(k_c), 100, b_cache);
            if ~isempty(wl_k) && size(wl_k, 1) >= 3
                cx = mean(wl_k(:,1));  cy = mean(wl_k(:,2));
                ang = atan2(wl_k(:,2)-cy, wl_k(:,1)-cx);
                [~, ord] = sort(ang);
                conts_i{k_c} = wl_k(ord, 1:2);
                A_outer_k(k_c) = polyarea(wl_k(ord,1), wl_k(ord,2));
            end
        end
        contours_outer{i} = conts_i;

        % ── L(z) slope-correction for non-vertical surfaces ──
        % Converts perpendicular thickness t_UHPC to in-plane
        % offset distance used by offset_polygon.
        dA_dz_k = gradient(A_outer_k, z_i);
        if has_P_table
            P_k = max(1e-6, interp1(config.Aw_table_z, ...
                      config.P_table, z_i, 'linear', 0));
            L_k = sqrt(1 + (dA_dz_k ./ P_k).^2);
        else
            R_eq_k = sqrt(max(A_outer_k, 0) / pi);
            L_k    = sqrt(1 + gradient(R_eq_k, z_i).^2);
        end
        L_k(A_outer_k < 1e-10) = 1;
        L_k = max(1, L_k);
        % CAP L_k to match the thin-shell inner_properties_at_z:
        %   cos_alpha = max(cos_alpha, 1/max_slope_factor)
        %   → L_k = 1/cos_alpha ≤ max_slope_factor
        % Without this cap, steep hull zones produce t_inplane >> t_UHPC,
        % collapsing inner polygons and over-assigning volume to UHPC (~10%
        % mass discrepancy vs Phase 1). With cap, strip integration is
        % consistent with Phase 1's integration model.
        % UHPC-specific cap must precede steel cap; otherwise this extraction diverges from the
        % solver when the two caps differ.
        if isfield(config, 'uhpc_max_slope_factor') && isfinite(config.uhpc_max_slope_factor)
            L_k_cap = config.uhpc_max_slope_factor;
        elseif isfield(config, 'steel_max_slope_factor') && isfinite(config.steel_max_slope_factor)
            L_k_cap = config.steel_max_slope_factor;
        else
            L_k_cap = 5.0;   % the thin-shell solve's own default
        end
        L_k = min(L_k, L_k_cap);

        % ── Per-sample arrays ────────────────────────────────
        inner_conts     = cell(n_sub_i, 1);
        A_UHPC_samples  = zeros(n_sub_i, 1);   % UHPC area per sample
        A_void_samples  = zeros(n_sub_i, 1);   % void area per sample
        Iyy_UHPC_samp   = zeros(n_sub_i, 1);   % Iyy of UHPC per sample
        Iyy_void_samp   = zeros(n_sub_i, 1);   % Iyy of void per sample
        r_min_strip     = Inf;
        t_perp_min      = Inf;

        % ── UNIFIED z-SAMPLE DECISION LOOP ──────────────────
        %
        %  is_solid_sample: true when the strip is the structural
        %  wall, when Phase 1b upgraded the strip to fully solid,
        %  OR when the sample elevation is at or below z_ballast.
        %  Every other sample gets a UHPC annular jacket + void.
        %
        %  Per-strip jacket thickness:
        %    t_UHPC_i = t_offset_strip(i)  (NEW Phase 1b: per-strip)
        t_UHPC_i = t_offset_strip(i);
        strip_solid = t_strip_is_solid_input(i);

        for k = 1:n_sub_i
            z_k     = z_i(k);
            outer_k = conts_i{k};
            A_out_k = A_outer_k(k);

            is_solid_sample = strip_solid || (z_k <= z_ballast);

            if is_solid_sample
                % ── SOLID UHPC: entire cross-section ────────
                % No offset, no inner contour, no void.
                A_UHPC_samples(k) = A_out_k;
                A_void_samples(k) = 0;
                if ~isempty(outer_k) && size(outer_k, 1) >= 3
                    [~, iner_out, ~] = mwecmass.hydrostatics.polygon_properties( ...
                        outer_k(:,1), outer_k(:,2));
                    Iyy_UHPC_samp(k) = abs(iner_out(2));
                end
                % inner_conts{k} stays empty → no void boundary

            else
                % ── ANNULAR UHPC JACKET + VOID INTERIOR ─────
                if isempty(outer_k) || size(outer_k, 1) < 3
                    % No valid contour — treat as solid UHPC
                    A_UHPC_samples(k) = A_out_k;
                    continue;
                end

                t_inplane_k = t_UHPC_i * L_k(k);
                [x_off, y_off] = mwecmass.internal.offset_polygon( ...
                    outer_k(:,1), outer_k(:,2), t_inplane_k);

                if length(x_off) < 3
                    % Offset polygon collapsed — treat as solid
                    [geom_out, iner_out, ~] = mwecmass.hydrostatics.polygon_properties( ...
                        outer_k(:,1), outer_k(:,2));
                    A_UHPC_samples(k)  = geom_out(1);
                    Iyy_UHPC_samp(k)   = abs(iner_out(2));
                    continue;
                end

                [geom_in,  iner_in,  ~] = mwecmass.hydrostatics.polygon_properties(x_off, y_off);
                [geom_out, iner_out, ~] = mwecmass.hydrostatics.polygon_properties( ...
                    outer_k(:,1), outer_k(:,2));
                A_in_k  = geom_in(1);
                A_out_recomputed = geom_out(1);

                if A_in_k <= 1e-10 || A_in_k >= A_out_recomputed
                    % Inner polygon collapsed or inverted — treat as solid
                    A_UHPC_samples(k)  = A_out_recomputed;
                    Iyy_UHPC_samp(k)   = abs(iner_out(2));
                else
                    % Valid annular cross-section (physics)
                    A_UHPC_samples(k) = A_out_recomputed - A_in_k;
                    A_void_samples(k) = A_in_k;
                    Iyy_UHPC_samp(k)  = max(0, abs(iner_out(2)) - abs(iner_in(2)));
                    Iyy_void_samp(k)  = abs(iner_in(2));

                    cx_in  = geom_in(2);
                    cy_in  = geom_in(3);
                    r_in_k = min(sqrt((x_off - cx_in).^2 + (y_off - cy_in).^2));
                    if r_in_k < r_min_strip
                        r_min_strip = r_in_k;
                    end
                    t_perp_min = min(t_perp_min, t_UHPC_i);

                    % VISUALIZATION inner contour: store the L_k polygon
                    % so the rendered void volume matches the integrated
                    % physics.  Earlier code stored a direct-t_UHPC viz
                    % polygon, but that systematically understates the
                    % UHPC ring on sloped surfaces (the visual void
                    % appears too large vs the reported %).  Using the
                    % L_k polygon — the same one used for A_UHPC_samples
                    % — keeps visualization and accounting consistent.
                    inner_conts{k} = [x_off, y_off];
                end
            end
        end  % k-sample loop

        % ── Integrate over strip ─────────────────────────────
        contours_inner{i} = inner_conts;
        strip_V_total(i)  = V_i;

        if n_sub_i >= 2
            V_UHPC_i = max(0, trapz(z_i, A_UHPC_samples));
            V_void_i = max(0, trapz(z_i, A_void_samples));
        else
            V_UHPC_i = A_UHPC_samples(1) * (z_hi - z_lo);
            V_void_i = A_void_samples(1)  * (z_hi - z_lo);
        end

        strip_V_UHPC(i)     = V_UHPC_i;
        strip_V_void(i)     = V_void_i;
        strip_mass_UHPC(i)  = V_UHPC_i * rho_UHPC;
        strip_mass_void(i)  = V_void_i * rho_air;
        strip_mass_total(i) = strip_mass_UHPC(i) + strip_mass_void(i);
        strip_z_cg(i)       = CBz_i;
        % Divide by the envelope built from the SAME integrator as the
        % parts (V_UHPC_i and V_void_i both come from trapz over the
        % polyarea contour samples).  Mixing integrators leaves a
        % 0.3-0.4% residual that pushes fully-solid strips above
        % rho_UHPC, which is physically impossible.  With this form
        % rho_eff is a mass-weighted average of two densities and is
        % therefore bounded by them by construction.
        % Units: V_env_i [m^3], strip_rho_eff [kg/m^3], strip_scale [-].
        V_env_i             = V_UHPC_i + V_void_i;
        strip_rho_eff(i)    = strip_mass_total(i) / max(V_env_i, eps);
        strip_scale(i)      = sqrt(max(V_void_i, 0) / max(V_env_i, eps));

        % Wall/solid classification
        strip_is_wall(i)  = (i == wall_strip_idx);
        strip_is_solid(i) = strip_is_wall(i) || (z_hi <= z_ballast) || ...
                            t_strip_is_solid_input(i);

        % Feasibility (only meaningful for strips with any void)
        any_void_samples = any(A_void_samples > 1e-10);
        if any_void_samples && isfinite(t_perp_min) && t_perp_min < 1e6
            strip_t_actual(i)    = t_perp_min;
            strip_r_min(i)       = max(0, r_min_strip);
            strip_is_feasible(i) = (t_perp_min >= t_min - 1e-6);
        end
        % Solid strips keep strip_t_actual = Inf, strip_is_feasible = true

        if n_sub_i >= 2
            strip_Iyy_UHPC(i) = rho_UHPC * ...
                mwecmass.internal.integrate_piecewise_cubic( ...
                    spline(z_i, Iyy_UHPC_samp), z_lo, z_hi);
            strip_Iyy_void(i) = rho_air * ...
                mwecmass.internal.integrate_piecewise_cubic( ...
                    spline(z_i, Iyy_void_samp), z_lo, z_hi);
        end
    end  % strip loop

    %% PACKAGE OUTPUT
    cstr.mode              = 'constructable_hull';
    cstr.fill_method       = 'uhpc_fill';   % 'mode' above is a distinct, unrelated, unread solver-identity tag -- not touched
    cstr.Z_max             = Z_max;
    cstr.Z_min             = Z_min;

    % Structural wall boundaries (from strip edges, not z_ballast)
    cstr.wall_strip_idx    = wall_strip_idx;
    cstr.wall_z_bottom     = strip_edges(wall_strip_idx);
    cstr.wall_z_top        = strip_edges(wall_strip_idx + 1);
    cstr.rho_hull          = rho_UHPC;
    cstr.rho_air           = rho_air;
    cstr.t_min             = t_min;
    cstr.wall_height       = strip_edges(wall_strip_idx + 1) - ...
                             strip_edges(wall_strip_idx);

    % Expose realised Phase 1+1b DOFs explicitly so downstream
    % consumers (visualize, verify) can use them directly.
    % t_UHPC remains scalar for compatibility with existing consumers.
    if isfield(solve_data, 't_uhpc')
        cstr.t_UHPC = solve_data.t_uhpc;   % global Phase 1a thickness
    else
        cstr.t_UHPC = max(t_offset_strip(isfinite(t_offset_strip)));
    end
    cstr.t_offset_strip   = t_offset_strip;   % per-strip thickness (Inf=solid)
    cstr.z_ballast           = z_ballast;
    cstr.vertical_shift   = solve_data.vertical_shift;
    cstr.draft            = solve_data.draft;

    cstr.strip_z_lo           = strip_z_lo;
    cstr.strip_z_hi           = strip_z_hi;
    cstr.strip_rho_eff        = strip_rho_eff;
    cstr.strip_scale_factor   = strip_scale;
    cstr.strip_V_total        = strip_V_total;
    cstr.strip_V_UHPC         = strip_V_UHPC;
    cstr.strip_V_void         = strip_V_void;
    cstr.strip_mass_UHPC      = strip_mass_UHPC;
    cstr.strip_mass_void      = strip_mass_void;
    cstr.strip_mass_total     = strip_mass_total;
    cstr.strip_t_min_actual   = strip_t_actual;
    cstr.strip_r_min          = strip_r_min;
    cstr.strip_z_cg           = strip_z_cg;
    cstr.strip_is_wall        = strip_is_wall;    % structural wall strip only
    cstr.strip_is_solid       = strip_is_solid;   % wall OR entirely below z_ballast
    cstr.strip_is_feasible    = strip_is_feasible;
    cstr.strip_Iyy_UHPC       = strip_Iyy_UHPC;
    cstr.strip_Iyy_void       = strip_Iyy_void;

    cstr.contours_outer       = contours_outer;
    cstr.contours_inner       = contours_inner;

    cstr.total_mass           = solve_data.M_total;
    cstr.total_V_UHPC         = solve_data.V_uhpc;
    cstr.total_V_void         = solve_data.V_air;
    cstr.total_V_hull         = solve_data.V_hull;
    cstr.UHPC_volume_fraction = solve_data.V_uhpc / max(solve_data.V_hull, eps);

    % Feasibility count excludes solid strips (they cannot violate t_min)
    n_infeas = sum(~strip_is_feasible & ~strip_is_solid);
    cstr.feasibility.n_infeasible       = n_infeas;
    cstr.feasibility.all_feasible       = (n_infeas == 0);
    cstr.feasibility.rho_min_achievable = zeros(N, 1);
    cstr.feasibility.s_max              = zeros(N, 1);

    % Console summary
    n_solid_ballast = sum(strip_is_solid & ~strip_is_wall);
    fprintf('\n      Per-strip geometry extracted (t_UHPC* = %.4f m, z_ballast = %.4f m)\n', ...
            cstr.t_UHPC, z_ballast);
    fprintf('        Per-strip t_offset (mm): ');
    for ii_p = 1:N
        if t_strip_is_solid_input(ii_p)
            fprintf('S ');
        else
            fprintf('%.0f ', t_offset_strip(ii_p)*1000);
        end
    end
    fprintf('\n');
    fprintf('        Wall strip: %d | Solid-ballast strips: %d | Annular strips: %d\n', ...
            1, n_solid_ballast, N - 1 - n_solid_ballast);
    fprintf('        V_hull = %.4f m³  V_UHPC = %.4f m³  V_void = %.4f m³\n', ...
            cstr.total_V_hull, cstr.total_V_UHPC, cstr.total_V_void);
    fprintf('        Mass = %.1f kg  GM = %.4f m  T_heave = %.3f s  T_pitch = %.3f s\n', ...
            cstr.total_mass, solve_data.GM_realised, ...
            solve_data.T_heave_realised, solve_data.T_pitch_realised);
    if n_infeas > 0
        fprintf('        WARNING: %d strip(s) have t_actual < t_min\n', n_infeas);
    else
        fprintf('        All annular strips feasible (t_actual >= %.1f mm)\n', t_min*1000);
    end

    % ── MASS BALANCE VERIFICATION ─────────────────────────────────
    %  Compare strip-level integrated mass against Phase 1 authoritative value.
    %  Discrepancy is EXPECTED and NORMAL — the two use different integration
    %  methods (strip trapz vs. precomputed Aw_table / V_sub_table fzero).
    %  The pipeline uses Phase 1 (solve_data) for final_props, NOT the strip sum.
    strip_mass_sum    = sum(strip_mass_total);
    mass_balance_err  = strip_mass_sum - solve_data.M_total;
    mass_balance_pct  = 100 * mass_balance_err / max(solve_data.M_total, eps);
    fprintf('\n      ── Mass Balance Verification ──────────────────────────\n');
    fprintf('        Phase 1 M_total (authoritative) : %.1f kg\n', solve_data.M_total);
    fprintf('        Strip integration sum           : %.1f kg  (%+.2f%%)\n', ...
            strip_mass_sum, mass_balance_pct);
    if abs(mass_balance_pct) < 2.0
        fprintf('        Status: OK (strip/Phase-1 discrepancy < 2%%)\n');
    else
        fprintf('        Status: LARGE DISCREPANCY (%.1f%%) — expected at transition zones\n', ...
                mass_balance_pct);
        fprintf('                Physics (GM, T_heave, draft) from Phase 1 is still correct.\n');
    end
    fprintf('        NOTE: final_props uses Phase 1 values.  Strip data is visual only.\n');
    fprintf('      ────────────────────────────────────────────────────────\n');
end
