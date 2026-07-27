classdef WEC_Constructable_Hull
    % WEC_CONSTRUCTABLE_HULL  Post-processing module: translate optimised
    %   effective densities into physically realisable UHPC + void geometry.
    %
    %   All methods are static (no instance state), following the pattern
    %   of WEC_Shell_Offset.
    %
    %   PHYSICAL MODEL
    %     The hull is split into two regions:
    %       (a) WALL REGION: top of hull down by wall_height [m].  Solid UHPC.
    %           Not optimised.  Single material, uniform density.
    %       (b) PLATFORM REGION: everything below the wall.  N strips, each
    %           with an optimiser-assigned effective density rho_eff_i.
    %           Realised as a UHPC annular shell (cross-section contour
    %           scaled inward about its centroid) with an air void interior.
    %
    %   SCALED-CONTOUR METHOD
    %     At each z-level within a strip, the hull cross-section contour
    %     is scaled by factor s_i about the contour centroid.  The region
    %     between the original and scaled contours is UHPC; the interior
    %     is void (air).
    %
    %     Scale factor:
    %       s_i = sqrt( (rho_UHPC - rho_eff_i) / (rho_UHPC - rho_air) )
    %
    %     This produces a UHPC wall whose thickness varies naturally with
    %     the hull shape — thicker where the hull is wider, thinner where
    %     it narrows — which is structurally rational.
    %
    %   CONSTRUCTABILITY
    %     Each strip is cast as a separate "cake layer" of UHPC with a
    %     void cavity shaped by the scaled contour.  Strips are connected
    %     by post-tensioned tendons.  Bilateral symmetry (xz and yz planes)
    %     is preserved by scaling about the centroid.
    %
    %   METHOD INVENTORY
    %     realize                 — Main entry point (post-processing)
    %     scale_contour           — Scale a 2D polygon about its centroid
    %     compute_strip_volume    — Integrate A(z) over a z-range via trapz
    %     compute_annular_props   — Area, centroid, Iyy of annular region
    %     verify_properties       — Compare actual vs idealised mass properties
    %     check_feasibility       — Per-strip minimum wall thickness check
    %     visualize               — Multi-panel figure of realised geometry
    %     validate_cylinder       — Analytical validation (closed-form cylinder)
    %
    %   DEPENDENCIES
    %     WEC_HydroProperties.extract_isocurve_at_z  (cross-section contour)
    %     WEC_Core_Functions.polygeom                     (shoelace formula)
    %
    %   See also: WEC_Shell_Offset, WEC_Configuration_Builder,
    %             calculate_3d_properties, WEC_Main_Optimizer
    %
    %   Author:  WEC Optimisation Team
    %   Version: 1.0 — Constructability post-processing module

    methods (Static)

        %% ═════════════════════════════════════════════════════════════
        %%  MAIN ENTRY POINT
        %% ═════════════════════════════════════════════════════════════

        function cstr = realize(config, x_opt, final_props)
        % REALIZE  Post-process converged optimiser output into a
        %   physically realisable UHPC + void geometry.
        %
        %   cstr = WEC_Constructable_Hull.realize(config, x_opt, final_props)
        %
        %   Called ONCE after Stage 2 converges.  Does not modify any
        %   optimiser output — pure post-processing.
        %
        %   INPUTS
        %     config      : struct from WEC_Configuration_Builder
        %     x_opt       : [1+N x 1]  converged design vector [draft; rho_1..rho_N]
        %     final_props : struct from calculate_3d_properties(x_opt, config)
        %
        %   OUTPUT
        %     cstr : struct with per-strip void geometry, verification,
        %            and feasibility data (see §7 PACKAGE OUTPUT below)
        %
        %   ALGORITHM
        %     1. Separate wall region (solid UHPC) from platform region
        %     2. Compute scale factor per platform strip
        %     3. Extract contours, compute volumes, generate void boundaries
        %     4. Verify mass properties against idealised (uniform rho_eff)
        %     5. Check structural feasibility (min wall thickness)

            fprintf('\n    WEC_Constructable_Hull.realize (UHPC two-phase global solve):\n');

            %% §1  RESOLVE UHPC MATERIAL PARAMETERS  ─────────────────
            rho_UHPC     = config.constructability_rho_hull;
            rho_void_mat = config.constructability_rho_fill;
            t_min_uhpc   = config.constructability_t_min;

            % Build opts struct overriding Shell_Offset defaults with UHPC params.
            % opt_or_cfg in Shell_Offset reads opts field first, then config field,
            % so overriding here does NOT require changing config.
            uhpc_opts = struct();
            uhpc_opts.rho_steel = rho_UHPC;
            uhpc_opts.rho_air   = rho_void_mat;
            uhpc_opts.t_min     = t_min_uhpc;

            % Warm-start: UHPC needs ~rho_steel/rho_UHPC times more wall volume
            % than steel for the same mass, so scale t_init proportionally.
            if isfield(config, 'uhpc_t_init') && ~isempty(config.uhpc_t_init)
                uhpc_opts.t_init = config.uhpc_t_init;
            else
                uhpc_opts.t_init = config.steel_t_init * ...
                                   (config.rho_steel / rho_UHPC);
            end
            if isfield(config, 'uhpc_max_slope_factor')
                uhpc_opts.max_slope_factor = config.uhpc_max_slope_factor;
            end
            if isfield(config, 'uhpc_n_z_grid')
                uhpc_opts.n_z_grid = config.uhpc_n_z_grid;
            end
            % NOTE: uhpc_w_GM/Th/Tp were retired in the fmincon redesign.
            % The new solve_constructable uses range-normalised
            % phi(r_heave) + phi(r_pitch); GM and mass are constraints,
            % not weighted objective terms.

            fprintf('      rho_UHPC = %.0f kg/m³,  rho_void = %.2f kg/m³\n', ...
                    rho_UHPC, rho_void_mat);
            fprintf('      t_min = %.4f m (%.1f mm),  t_init = %.4f m\n', ...
                    t_min_uhpc, t_min_uhpc*1000, uhpc_opts.t_init);

            %% §2  FEASIBILITY PRE-CHECK  ──────────────────────────────
            %
            %  At t = t_min, M_min = jacket-only mass (void inside).
            %  At t = inf (fully solid), M_max = rho_UHPC * V_hull.
            %  Target mass must lie within [M_min, M_max].

            n_z_pre   = 100;
            slope_f   = WEC_Shell_Offset.opt_or_cfg(uhpc_opts, 'max_slope_factor', ...
                            config.steel_max_slope_factor);
            [grids_pre, ~] = WEC_Shell_Offset.build_geometry_grid( ...
                config, t_min_uhpc, n_z_pre, slope_f);
            A_jacket_pre = grids_pre.A_outer - grids_pre.A_inner;
            M_min_uhpc   = rho_UHPC    * trapz(grids_pre.z, A_jacket_pre) + ...
                           rho_void_mat * trapz(grids_pre.z, grids_pre.A_inner);
            M_max_uhpc   = rho_UHPC    * trapz(grids_pre.z, grids_pre.A_outer);

            fprintf('      Achievable mass range: [%.0f, %.0f] kg  |  target: %.0f kg\n', ...
                    M_min_uhpc, M_max_uhpc, final_props.mass_total);

            if final_props.mass_total < M_min_uhpc * 0.95
                error('WEC_Constructable_Hull:MassTooLight', ...
                      ['Target mass %.1f kg < minimum UHPC+void mass %.1f kg ', ...
                       '(at t_min=%.4f m). Hull is too light for UHPC construction.'], ...
                      final_props.mass_total, M_min_uhpc, t_min_uhpc);
            end
            if final_props.mass_total > M_max_uhpc * 1.05
                error('WEC_Constructable_Hull:MassTooHeavy', ...
                      ['Target mass %.1f kg > fully-solid UHPC mass %.1f kg. ', ...
                       'Hull is too heavy for UHPC; use a denser material.'], ...
                      final_props.mass_total, M_max_uhpc);
            end

            %% §3  PHASE 1 — WALL-AWARE GLOBAL SOLVE + PHASE 1b  ───────
            %
            %  WEC_Shell_Offset.solve_constructable runs:
            %    Phase 1a — global (t_UHPC*, z_fill*) under wall=solid,
            %    Phase 1b — bottom-up greedy per-strip thickening to
            %               jointly close mass + GM residuals.
            %
            %  Returns the steel_data-shaped struct PLUS:
            %    .t_offset_strip   per-strip jacket thickness [m]  (Inf = solid)
            %    .is_solid_strip   per-strip solid flag
            %    .strip_edges      strip boundary z-coords (body frame)
            %    .wall_strip_idx   structural wall strip index
            %
            %  Wall info: derive strip_edges and wall_strip_idx the same
            %  way extract_strip_geometry does, then pass through opts.

            N_strips = length(config.density_nodes_z);
            density_nodes_z = config.density_nodes_z(:);
            if isfield(config, 'strip_edges') && ~isempty(config.strip_edges)
                strip_edges_realize = config.strip_edges(:);
            else
                strip_edges_realize = zeros(N_strips + 1, 1);
                strip_edges_realize(1) = config.hull_z_min;
                for ii = 2:N_strips
                    strip_edges_realize(ii) = 0.5 * ...
                        (density_nodes_z(ii-1) + density_nodes_z(ii));
                end
                strip_edges_realize(N_strips + 1) = config.hull_z_max;
            end
            if isfield(config, 'wall_position') && strcmp(config.wall_position, 'bottom')
                wall_strip_idx = 1;
            else
                wall_strip_idx = N_strips;
            end

            uhpc_opts.wall_strip_idx = wall_strip_idx;
            uhpc_opts.strip_edges    = strip_edges_realize;

            solve_data = WEC_Shell_Offset.solve_constructable( ...
                config, x_opt, final_props, uhpc_opts);

            if ~solve_data.feasible
                warning('WEC_Constructable_Hull:InfeasibleUHPCSolve', ...
                        ['UHPC constructable solve returned infeasible result. ' ...
                         'Targets may not be achievable under UHPC+void model.']);
            end
            solve_data.realisation_mode = 'uhpc_fill';

            %% §4  PHASE 2 — PER-STRIP GEOMETRY EXTRACTION  ────────────
            %
            %  Re-runs the strip loop with the per-strip t_offset vector
            %  AND is_solid flags from Phase 1b.  Produces contours_outer,
            %  contours_inner, per-strip volumes/masses for visualisation
            %  and the realised_strips substruct in build_realised_props.

            cstr = WEC_Constructable_Hull.extract_strip_geometry( ...
                       config, solve_data.t_offset_strip, solve_data.z_fill, ...
                       rho_UHPC, rho_void_mat, solve_data);

            fprintf('    WEC_Constructable_Hull.realize: done.\n');
        end


        %% ═══════════════════════════════════════════════════════════════
        %%  EXTRACT_STRIP_GEOMETRY — Phase 2 geometry extraction (UHPC)
        %% ═══════════════════════════════════════════════════════════════

        function cstr = extract_strip_geometry(config, t_offset_arg, z_fill, ...
                                               rho_UHPC, rho_void, solve_data)
        % EXTRACT_STRIP_GEOMETRY  Re-extract per-strip contours and volumes
        %   from per-strip (t_offset_strip, z_fill) for downstream visualisation.
        %
        %   t_offset_arg may be a SCALAR (legacy uniform-thickness mode) or a
        %   VECTOR (Phase 1b per-strip thickness — Inf entries mean fully solid).
        %   When solve_data carries .is_solid_strip, those strips are forced
        %   to solid regardless of t_offset_arg.
        %
        %   Called after Phase 1 global solve + Phase 1b thickening.  cstr
        %   inherits from solve_data (mass / GM / inertia) and adds per-strip
        %   visualisation fields.
        %
        %   TOPOLOGY
        %     z <= z_fill            → solid UHPC (sample-level)
        %     strip i is wall/solid  → strip is solid UHPC
        %     otherwise              → UHPC annular jacket + void interior

            % Inherit all global solve fields (M_total, z_fill, t_steel, etc.)
            cstr = solve_data;

            %% §1  UNPACK CONFIG  ────────────────────────────────────
            ms2_model   = config.ms2_model;
            n_sub       = config.constructability_n_sub;
            N           = length(config.density_nodes_z);
            Z_max       = config.hull_z_max;
            Z_min       = config.hull_z_min;
            t_min       = config.constructability_t_min;
            wall_height = config.constructability_wall_height;

            if isfield(config, 'wall_position') && strcmp(config.wall_position, 'bottom')
                wall_z_bottom = Z_min + wall_height;
            else
                wall_z_bottom = Z_max - wall_height;
            end

            if isfield(config, 'boundary_cache') && ~isempty(config.boundary_cache)
                b_cache = config.boundary_cache;
            else
                b_cache = WEC_HydroProperties.precompute_boundary_cache(ms2_model, 100);
            end

            has_P_table = isfield(config, 'P_table') && ~isempty(config.P_table) && ...
                          length(config.P_table) == length(config.Aw_table_z);

            %% §2  STRIP BOUNDARIES  ─────────────────────────────────
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

            %% §3  PER-STRIP UNIFIED LOOP  ───────────────────────────
            %
            %  UNIFIED PHYSICAL MODEL (replaces the old two-branch approach):
            %
            %    At every z-sample within every strip, the decision is made
            %    at the SAMPLE level — not the strip level:
            %
            %      is_solid_sample = (this is the designated wall strip)
            %                     OR (z_sample <= z_fill)
            %
            %    Solid samples → A_UHPC = A_outer,  A_void = 0  (no inner contour)
            %    Annular samples → A_UHPC = A_outer − A_inner(t_UHPC),  A_void = A_inner
            %
            %    Strips that straddle z_fill are handled automatically:
            %    their lower samples are solid and upper samples are annular.
            %    trapz integrates both contributions without any explicit split.
            %
            %  WALL POSITION:
            %    wall_position = 'top'    → wall_strip_idx = N (highest strip)
            %    wall_position = 'bottom' → wall_strip_idx = 1 (_180 flag)
            %
            %  strip_is_wall:  the DESIGNATED structural wall strip only
            %  strip_is_solid: wall strip OR strips whose ENTIRE z-range is ≤ z_fill
            %                  (used for rendering — partially-solid strips report false
            %                   because they have annular samples at their top)

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
                    error('WEC_Constructable_Hull:BadTOffsetLength', ...
                          ['t_offset_arg length %d does not match strip count %d.'], ...
                          length(t_offset_strip), N);
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
            strip_is_solid    = false(N, 1);  % wall OR entirely below z_fill
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

                strip_phys = WEC_HydroProperties.compute_strip(ms2_model, z_lo, z_hi);
                V_i   = strip_phys.V;
                CBz_i = strip_phys.CB_z;

                % ── Outer contours at each z-sample ─────────────────
                conts_i   = cell(n_sub_i, 1);
                A_outer_k = zeros(n_sub_i, 1);
                for k_c = 1:n_sub_i
                    wl_k = WEC_HydroProperties.extract_isocurve_at_z( ...
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
                % offset distance used by offset_vertices_raw.
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
                % CAP L_k to match WEC_Shell_Offset.inner_props_at_z line 1407:
                %   cos_alpha = max(cos_alpha, 1/max_slope_factor)
                %   → L_k = 1/cos_alpha ≤ max_slope_factor
                % Without this cap, steep hull zones produce t_inplane >> t_UHPC,
                % collapsing inner polygons and over-assigning volume to UHPC (~10%
                % mass discrepancy vs Phase 1). With cap, strip integration is
                % consistent with Phase 1's integration model.
                if isfield(config, 'steel_max_slope_factor') && isfinite(config.steel_max_slope_factor)
                    L_k_cap = config.steel_max_slope_factor;
                elseif isfield(config, 'max_slope_factor') && isfinite(config.max_slope_factor)
                    L_k_cap = config.max_slope_factor;
                else
                    L_k_cap = 5.0;   % WEC_Shell_Offset default
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
                %  OR when the sample elevation is at or below z_fill.
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

                    is_solid_sample = strip_solid || (z_k <= z_fill);

                    if is_solid_sample
                        % ── SOLID UHPC: entire cross-section ────────
                        % No offset, no inner contour, no void.
                        A_UHPC_samples(k) = A_out_k;
                        A_void_samples(k) = 0;
                        if ~isempty(outer_k) && size(outer_k, 1) >= 3
                            [~, iner_out, ~] = WEC_Core_Functions.polygeom( ...
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
                        [x_off, y_off] = WEC_Shell_Offset.offset_vertices_raw( ...
                            outer_k(:,1), outer_k(:,2), t_inplane_k);

                        if length(x_off) < 3
                            % Offset polygon collapsed — treat as solid
                            [geom_out, iner_out, ~] = WEC_Core_Functions.polygeom( ...
                                outer_k(:,1), outer_k(:,2));
                            A_UHPC_samples(k)  = geom_out(1);
                            Iyy_UHPC_samp(k)   = abs(iner_out(2));
                            continue;
                        end

                        [geom_in,  iner_in,  ~] = WEC_Core_Functions.polygeom(x_off, y_off);
                        [geom_out, iner_out, ~] = WEC_Core_Functions.polygeom( ...
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
                strip_mass_void(i)  = V_void_i * rho_void;
                strip_mass_total(i) = strip_mass_UHPC(i) + strip_mass_void(i);
                strip_z_cg(i)       = CBz_i;
                strip_rho_eff(i)    = strip_mass_total(i) / max(V_i, eps);
                strip_scale(i)      = sqrt(max(V_void_i, 0) / max(V_i, eps));

                % Wall/solid classification
                strip_is_wall(i)  = (i == wall_strip_idx);
                strip_is_solid(i) = strip_is_wall(i) || (z_hi <= z_fill) || ...
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
                        WEC_Core_Functions.integrate_pp( ...
                            spline(z_i, Iyy_UHPC_samp), z_lo, z_hi);
                    strip_Iyy_void(i) = rho_void * ...
                        WEC_Core_Functions.integrate_pp( ...
                            spline(z_i, Iyy_void_samp), z_lo, z_hi);
                end
            end  % strip loop

            %% §4  PACKAGE OUTPUT  ──────────────────────────────────
            cstr.mode              = 'constructable_hull';
            cstr.realisation_mode  = 'uhpc_fill';
            cstr.Z_max             = Z_max;
            cstr.Z_min             = Z_min;

            % Structural wall boundaries (from strip edges, not z_fill)
            cstr.wall_strip_idx    = wall_strip_idx;
            cstr.wall_z_bottom     = strip_edges(wall_strip_idx);
            cstr.wall_z_top        = strip_edges(wall_strip_idx + 1);
            cstr.rho_hull          = rho_UHPC;
            cstr.rho_fill          = rho_void;
            cstr.t_min             = t_min;
            cstr.wall_height       = strip_edges(wall_strip_idx + 1) - ...
                                     strip_edges(wall_strip_idx);

            % Expose realised Phase 1+1b DOFs explicitly so downstream
            % consumers (visualize, verify) can use them directly.
            % t_UHPC stays as a scalar for backward compatibility (it equals
            % the Phase 1a global thickness, used for legacy single-offset plots).
            if isfield(solve_data, 't_steel')
                cstr.t_UHPC = solve_data.t_steel;   % global Phase 1a thickness
            else
                cstr.t_UHPC = max(t_offset_strip(isfinite(t_offset_strip)));
            end
            cstr.t_offset_strip   = t_offset_strip;   % per-strip thickness (Inf=solid)
            cstr.z_fill           = z_fill;
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
            cstr.strip_is_solid       = strip_is_solid;   % wall OR entirely below z_fill
            cstr.strip_is_feasible    = strip_is_feasible;
            cstr.strip_Iyy_UHPC       = strip_Iyy_UHPC;
            cstr.strip_Iyy_void       = strip_Iyy_void;

            cstr.contours_outer       = contours_outer;
            cstr.contours_inner       = contours_inner;

            cstr.total_mass           = solve_data.M_total;
            cstr.total_V_UHPC         = solve_data.V_steel;
            cstr.total_V_void         = solve_data.V_air;
            cstr.total_V_hull         = solve_data.V_hull;
            cstr.UHPC_volume_fraction = solve_data.V_steel / max(solve_data.V_hull, eps);

            % Feasibility count excludes solid strips (they cannot violate t_min)
            n_infeas = sum(~strip_is_feasible & ~strip_is_solid);
            cstr.feasibility.n_infeasible       = n_infeas;
            cstr.feasibility.all_feasible       = (n_infeas == 0);
            cstr.feasibility.rho_min_achievable = zeros(N, 1);
            cstr.feasibility.s_max              = zeros(N, 1);

            % Console summary
            n_solid_fill = sum(strip_is_solid & ~strip_is_wall);
            fprintf('\n      Per-strip geometry extracted (t_UHPC* = %.4f m, z_fill = %.4f m)\n', ...
                    cstr.t_UHPC, z_fill);
            fprintf('        Per-strip t_offset (mm): ');
            for ii_p = 1:N
                if t_strip_is_solid_input(ii_p)
                    fprintf('S ');
                else
                    fprintf('%.0f ', t_offset_strip(ii_p)*1000);
                end
            end
            fprintf('\n');
            fprintf('        Wall strip: %d | Solid-fill strips: %d | Annular strips: %d\n', ...
                    1, n_solid_fill, N - 1 - n_solid_fill);
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


        %%  SCALE CONTOUR
        %% ═════════════════════════════════════════════════════════════

        function [inner_xy, centroid] = scale_contour(contour_xy, s)
        % SCALE_CONTOUR  Scale a 2D polygon about its centroid.
        %
        %   [inner_xy, centroid] = WEC_Constructable_Hull.scale_contour(contour_xy, s)
        %
        %   The void boundary in each strip is obtained by scaling the
        %   cross-section contour inward.  s = 0 collapses to the
        %   centroid; s = 1 returns the original contour.
        %
        %   WHY scale about the centroid (not the origin)?
        %     The hull cross-section centroid may not coincide with (0,0)
        %     at some z-levels.  Scaling about the centroid preserves
        %     bilateral symmetry of the annular region regardless of
        %     centroid position.  For a double-symmetric hull, centroid
        %     is near (0,0) and results are nearly identical.
        %
        %   INPUTS
        %     contour_xy : [K x 2]  ordered polygon vertices (x, y)
        %     s          : scalar   scale factor in [0, 1]
        %
        %   OUTPUTS
        %     inner_xy   : [K x 2]  scaled polygon (void boundary)
        %     centroid   : [1 x 2]  polygon centroid used for scaling

            x = contour_xy(:, 1);
            y = contour_xy(:, 2);
            K = length(x);

            % Close polygon temporarily for shoelace formula
            if K < 3
                inner_xy = contour_xy;
                centroid = mean(contour_xy, 1);
                return;
            end

            % Remove duplicate closing vertex if present
            if abs(x(end) - x(1)) < 1e-12 && abs(y(end) - y(1)) < 1e-12
                x = x(1:end-1);
                y = y(1:end-1);
                K = length(x);
            end

            % Shoelace centroid
            xp = circshift(x, -1);
            yp = circshift(y, -1);
            a  = x .* yp - xp .* y;
            A  = 0.5 * sum(a);

            if abs(A) < 1e-14
                % Degenerate polygon — fall back to arithmetic mean
                centroid = [mean(x), mean(y)];
            else
                cx = sum((x + xp) .* a) / (6 * A);
                cy = sum((y + yp) .* a) / (6 * A);
                centroid = [cx, cy];
            end

            % Scale each vertex toward/from centroid
            inner_xy = centroid + s * ([x, y] - centroid);

            % Restore closing vertex if original had one
            if size(contour_xy, 1) > K
                inner_xy(end+1, :) = inner_xy(1, :);
            end
        end


        %% ═════════════════════════════════════════════════════════════
        %%  STRIP VOLUME VIA CROSS-SECTION INTEGRATION
        %% ═════════════════════════════════════════════════════════════

        function [V_strip, contours, A_samples, z_samples] = ...
                compute_strip_volume(ms2_model, z_lo, z_hi, n_sub)
        % COMPUTE_STRIP_VOLUME  Strip volume via parametric divergence theorem
        %   with contour shapes from evaluateCrossSectionMS2.
        %
        %   Uses WEC_HydroProperties.compute_strip for correct volume,
        %   and evaluateCrossSectionMS2 for contour shapes (used by
        %   constructability post-processing only).
        %
        %   See also: WEC_HydroProperties.compute_strip

            if n_sub < 2, n_sub = 2; end

            % Correct volume from divergence theorem
            strip_phys = WEC_HydroProperties.compute_strip(ms2_model, z_lo, z_hi);
            V_strip = strip_phys.V;

            % Contour shapes for post-processing
            z_samples = linspace(z_lo, z_hi, n_sub)';
            contours  = cell(n_sub, 1);
            A_samples = zeros(n_sub, 1);
            for k = 1:n_sub
                wp = WEC_Core_Functions.evaluateCrossSectionMS2( ...
                         ms2_model, z_samples(k));
                if ~isempty(wp) && ~isempty(wp{1}) && size(wp{1},1) >= 3
                    contours{k}  = wp{1}(:, 1:2);
                    A_samples(k) = polyarea(wp{1}(:,1), wp{1}(:,2));
                end
            end
        end


        %% ═════════════════════════════════════════════════════════════
        %%  ANNULAR CROSS-SECTION PROPERTIES
        %% ═════════════════════════════════════════════════════════════

        function [A_ann, cx_ann, cy_ann, Iyy_ann, Iyy_inn] = ...
                compute_annular_props(contour_outer, contour_inner)
        % COMPUTE_ANNULAR_PROPS  Area, centroid, and second moments of
        %   the annular region between an outer and inner polygon.
        %
        %   [A_ann, cx, cy, Iyy_ann, Iyy_inn] = ...
        %       WEC_Constructable_Hull.compute_annular_props(outer, inner)
        %
        %   WHY subtraction is valid:
        %     The inner polygon is entirely contained within the outer
        %     (guaranteed by the scaling construction with s < 1).  Both
        %     have consistent winding.  The shoelace formula gives second
        %     moments about the coordinate origin, so direct subtraction
        %     gives the annular moment about the same origin.
        %
        %   NOTE: Iyy values are second moments of area about the
        %     COORDINATE ORIGIN (not the annular centroid).  The caller
        %     applies the parallel-axis theorem when computing mass
        %     moment of inertia about the body CG.
        %
        %   INPUTS
        %     contour_outer : [K x 2]  outer polygon vertices (x, y)
        %     contour_inner : [K x 2]  inner polygon vertices (x, y)
        %
        %   OUTPUTS
        %     A_ann   : [m²]  annular area = A_outer - A_inner
        %     cx_ann  : [m]   annular centroid x
        %     cy_ann  : [m]   annular centroid y
        %     Iyy_ann : [m⁴]  Iyy of annular region (about origin)
        %     Iyy_inn : [m⁴]  Iyy of inner region (about origin, for void)

            A_ann = 0; cx_ann = 0; cy_ann = 0; Iyy_ann = 0; Iyy_inn = 0;

            if isempty(contour_outer) || size(contour_outer, 1) < 3
                return;
            end
            if isempty(contour_inner) || size(contour_inner, 1) < 3
                % No inner polygon — entire section is solid
                [geom_out, iner_out, ~] = WEC_Core_Functions.polygeom( ...
                    contour_outer(:,1), contour_outer(:,2));
                A_ann   = geom_out(1);
                cx_ann  = geom_out(2);
                cy_ann  = geom_out(3);
                Iyy_ann = abs(iner_out(2));
                Iyy_inn = 0;
                return;
            end

            % Outer properties (about coordinate origin)
            [geom_out, iner_out, ~] = WEC_Core_Functions.polygeom( ...
                contour_outer(:,1), contour_outer(:,2));
            A_out   = geom_out(1);
            cx_out  = geom_out(2);
            cy_out  = geom_out(3);
            Iyy_out = abs(iner_out(2));

            % Inner properties (about coordinate origin)
            [geom_in, iner_in, ~] = WEC_Core_Functions.polygeom( ...
                contour_inner(:,1), contour_inner(:,2));
            A_in   = geom_in(1);
            cx_in  = geom_in(2);
            cy_in  = geom_in(3);
            Iyy_in = abs(iner_in(2));

            % Annular properties by subtraction
            A_ann = A_out - A_in;
            if A_ann > 1e-12
                cx_ann = (A_out * cx_out - A_in * cx_in) / A_ann;
                cy_ann = (A_out * cy_out - A_in * cy_in) / A_ann;
            else
                cx_ann = cx_out;
                cy_ann = cy_out;
            end
            Iyy_ann = Iyy_out - Iyy_in;
            Iyy_ann = max(Iyy_ann, 0);   % guard against numerical noise

            % Return inner Iyy separately for void contribution
            Iyy_inn = Iyy_in;
        end


        %% ═════════════════════════════════════════════════════════════
        %%  PROPERTY VERIFICATION
        %% ═════════════════════════════════════════════════════════════

        function verification = verify_properties( ...
                strip_mass, strip_z_cg, strip_Iyy_UHPC, strip_Iyy_void, ...
                strip_V_total, draft, final_props, config)
        % VERIFY_PROPERTIES  Recompute mass properties from the realised
        %   UHPC + void geometry and compare against the idealised
        %   (uniform rho_eff) values from the optimiser.
        %
        %   Reports absolute and percentage discrepancies for:
        %     mass, CG_z, Iyy, GM, total volume
        %
        %   INPUTS
        %     strip_mass     : [N x 1 kg]    total mass per strip
        %     strip_z_cg     : [N x 1 m]     z-centroid per strip (body frame)
        %     strip_Iyy_UHPC : [N x 1]       UHPC Iyy contribution per strip
        %     strip_Iyy_void : [N x 1]       void Iyy contribution per strip
        %     strip_V_total  : [N x 1 m³]    total volume per strip
        %     draft          : [m]            vertical shift
        %     final_props    : struct         from calculate_3d_properties
        %     config         : struct         from WEC_Configuration_Builder

            fprintf('\n      ═══ PROPERTY VERIFICATION ═══\n');

            % Derive local variables from inputs
            N           = length(strip_mass);
            ms2_model   = config.ms2_model;   %#ok<NASGU> retained for potential future use
            strip_edges = config.strip_edges;

            %% §1-§3  MASS, CG, Iyy FROM REALISED STRIP GEOMETRY
            %
            %  Uses strip_mass (= cstr.strip_mass_total) from the unified
            %  z-sample loop in extract_strip_geometry.  These already
            %  account for solid UHPC below z_fill and annular jacket+void
            %  above z_fill.  The old optimizer-density × V_strip approach
            %  is retired (it ignored z_fill entirely).
            mass_actual = 0;
            cg_z_num    = 0;
            Iyy_actual  = 0;

            for ii = 1:N
                m_ii       = strip_mass(ii);           % realised strip mass [kg]
                cg_z_ii    = strip_z_cg(ii);           % body-frame centroid [m]
                mass_actual = mass_actual + m_ii;
                cg_z_num    = cg_z_num    + m_ii * cg_z_ii;
                Iyy_actual  = Iyy_actual  + strip_Iyy_UHPC(ii) + strip_Iyy_void(ii);
            end

            mass_ideal  = final_props.mass_total;
            err_mass_pct = 100 * (mass_actual - mass_ideal) / max(abs(mass_ideal), eps);

            % CG_z
            if abs(mass_actual) > 1e-6
                CG_z_body = cg_z_num / mass_actual;
            else
                CG_z_body = 0;
            end
            CG_z_actual = CG_z_body + draft;
            CG_z_ideal  = final_props.CG_total(3);
            err_CG_pct  = 100 * (CG_z_actual - CG_z_ideal) / max(abs(CG_z_ideal), eps);

            % Iyy (already computed in strip loop — about body origin,
            % shift to CG using parallel axis theorem)
            cg_orig = [0, 0, CG_z_actual - draft];
            Iyy_actual = Iyy_actual - mass_actual * cg_orig(3)^2;
            % Add back the correct parallel axis term about CG
            % (strip_Iyy already includes origin-referenced terms)

            Iyy_ideal   = final_props.Iyy;
            err_Iyy_pct = 100 * (Iyy_actual - Iyy_ideal) / max(abs(Iyy_ideal), eps);

            %% §4  METACENTRIC HEIGHT
            %
            %  KM is geometry-dependent (CB, waterplane area, V_sub).
            %  The void realization does not change the external hull
            %  geometry, so KM is unchanged.  Only CG_z shifts.

            KM_actual = final_props.KM;   % unchanged
            GM_actual = KM_actual - CG_z_actual;
            GM_ideal  = final_props.GM_L;
            err_GM_pct = 100 * (GM_actual - GM_ideal) / max(abs(GM_ideal), eps);

            %% §5  TOTAL VOLUME CHECK
            %
            %  Compare sum of strip volumes against the divergence-theorem
            %  total volume from WEC_HydroProperties.compute.

            V_strips = sum(strip_V_total);
            V_divthm = config.total_wec_volume;
            err_V_pct = 100 * (V_strips - V_divthm) / max(abs(V_divthm), eps);

            %% §6  REPORT

            thresholds = struct('mass', 0.5, 'CG_z', 2.0, ...
                                'Iyy', 5.0, 'GM', 5.0, 'V', 1.0);

            pass_mass = abs(err_mass_pct) < thresholds.mass;
            pass_CG   = abs(err_CG_pct)  < thresholds.CG_z;
            pass_Iyy  = abs(err_Iyy_pct) < thresholds.Iyy;
            pass_GM   = abs(err_GM_pct)   < thresholds.GM;
            pass_V    = abs(err_V_pct)    < thresholds.V;
            all_passed = pass_mass && pass_CG && pass_Iyy && pass_GM && pass_V;

            pf = @(ok) ternary(ok, 'PASS', 'FAIL');

            fprintf('      %-16s %12s %12s %8s %6s  (%s)\n', ...
                    'Quantity', 'Idealised', 'Actual', 'Error%', '', 'Threshold');
            fprintf('      %s\n', repmat('─', 1, 66));
            fprintf('      %-16s %12.1f %12.1f %+7.2f%%  %s  (<%.1f%%)\n', ...
                    'Mass [kg]', mass_ideal, mass_actual, err_mass_pct, ...
                    pf(pass_mass), thresholds.mass);
            fprintf('      %-16s %12.4f %12.4f %+7.2f%%  %s  (<%.1f%%)\n', ...
                    'CG_z [m]', CG_z_ideal, CG_z_actual, err_CG_pct, ...
                    pf(pass_CG), thresholds.CG_z);
            fprintf('      %-16s %12.1f %12.1f %+7.2f%%  %s  (<%.1f%%)\n', ...
                    'Iyy [kg·m²]', Iyy_ideal, Iyy_actual, err_Iyy_pct, ...
                    pf(pass_Iyy), thresholds.Iyy);
            fprintf('      %-16s %12.4f %12.4f %+7.2f%%  %s  (<%.1f%%)\n', ...
                    'GM_L [m]', GM_ideal, GM_actual, err_GM_pct, ...
                    pf(pass_GM), thresholds.GM);
            fprintf('      %-16s %12.4f %12.4f %+7.2f%%  %s  (<%.1f%%)\n', ...
                    'V_total [m³]', V_divthm, V_strips, err_V_pct, ...
                    pf(pass_V), thresholds.V);
            fprintf('      %s\n', repmat('─', 1, 66));

            if all_passed
                fprintf('      ✓ ALL CHECKS PASSED\n');
            else
                fprintf('      ✗ SOME CHECKS FAILED — review discrepancies above\n');
            end
            fprintf('      ═══════════════════════════════\n');

            % Pack output
            verification.mass_ideal    = mass_ideal;
            verification.mass_actual   = mass_actual;
            verification.err_mass_pct  = err_mass_pct;
            verification.CG_z_ideal    = CG_z_ideal;
            verification.CG_z_actual   = CG_z_actual;
            verification.err_CG_z_pct  = err_CG_pct;
            verification.Iyy_ideal     = Iyy_ideal;
            verification.Iyy_actual    = Iyy_actual;
            verification.err_Iyy_pct   = err_Iyy_pct;
            verification.GM_ideal      = GM_ideal;
            verification.GM_actual     = GM_actual;
            verification.err_GM_pct    = err_GM_pct;
            verification.V_strips      = V_strips;
            verification.V_divthm      = V_divthm;
            verification.err_V_pct     = err_V_pct;
            verification.thresholds    = thresholds;
            verification.all_passed    = all_passed;
        end


        %% ═════════════════════════════════════════════════════════════
        %%  FEASIBILITY CHECK
        %% ═════════════════════════════════════════════════════════════

        function feasibility = check_feasibility( ...
                strip_scale, strip_r_min, strip_t_actual, ...
                strip_is_wall, strip_is_feasible, ...
                strip_z_lo, strip_z_hi, strip_rho_eff, ...
                t_min, rho_hull, rho_fill)
        % CHECK_FEASIBILITY  Per-strip minimum wall thickness check.
        %
        %   Reports which strips violate t_min and computes the minimum
        %   achievable density for each strip given its geometry.
        %
        %   The minimum achievable density is:
        %     s_max_i = 1 - t_min / r_min_i
        %     rho_min_i = rho_hull - s_max_i² × (rho_hull - rho_fill)
        %
        %   If the optimiser converged to rho_eff_i < rho_min_i, the
        %   strip is structurally infeasible.

            N = length(strip_scale);
            n_infeasible = 0;

            rho_min_achievable = zeros(N, 1);
            s_max              = zeros(N, 1);

            fprintf('\n      ═══ FEASIBILITY CHECK ═══\n');
            fprintf('      t_min = %.1f mm (%.2f in)\n', t_min*1000, t_min/0.0254);

            for i = 1:N
                if strip_is_wall(i)
                    rho_min_achievable(i) = rho_hull;
                    s_max(i)              = 0;
                    continue;
                end

                if strip_r_min(i) > 0 && isfinite(strip_r_min(i))
                    s_max(i) = max(0, 1 - t_min / strip_r_min(i));
                    rho_min_achievable(i) = rho_hull - s_max(i)^2 * (rho_hull - rho_fill);
                else
                    s_max(i) = 0;
                    rho_min_achievable(i) = rho_hull;
                end

                if ~strip_is_feasible(i)
                    n_infeasible = n_infeasible + 1;
                    fprintf('      WARNING: Strip %d  z=[%.3f, %.3f] m\n', ...
                            i, strip_z_lo(i), strip_z_hi(i));
                    fprintf('        t_actual = %.1f mm < t_min = %.1f mm\n', ...
                            strip_t_actual(i)*1000, t_min*1000);
                    fprintf('        s = %.4f  (s_max = %.4f)\n', ...
                            strip_scale(i), s_max(i));
                    fprintf('        rho_eff = %.0f kg/m³  (min achievable = %.0f kg/m³)\n', ...
                            strip_rho_eff(i), rho_min_achievable(i));
                    fprintf('        r_min = %.4f m\n', strip_r_min(i));
                end
            end

            if n_infeasible == 0
                fprintf('      ✓ ALL STRIPS FEASIBLE\n');
            else
                fprintf('      ✗ %d STRIP(S) INFEASIBLE\n', n_infeasible);
                fprintf('        Consider raising the optimizer lower density bound\n');
                fprintf('        to per-strip rho_min values above.\n');
            end
            fprintf('      ═══════════════════════════════\n');

            feasibility.n_infeasible      = n_infeasible;
            feasibility.all_feasible      = (n_infeasible == 0);
            feasibility.rho_min_achievable = rho_min_achievable;
            feasibility.s_max             = s_max;
        end


        %% ═════════════════════════════════════════════════════════════
        %%  VISUALIZATION
        %% ═════════════════════════════════════════════════════════════

        function visualize(cstr, config)
        % VISUALIZE  Plot 1 (XZ midplane elevation) and Figure 2 (per-strip
        %   plan-view cross-sections) of the as-built UHPC + void geometry.
        %
        %   TWO-STEP METHOD per strip:
        %     STEP A — perpendicular offset (preserved exactly):
        %       Compute the inner offset of the FULL world-frame hull
        %       profile (Plot 1) or the strip's plan-view contour
        %       (Figure 2) by perpendicular distance t_strip(i).  Plot 1
        %       uses offset_vertices_raw on the FULL profile then clips
        %       per strip — this avoids horizontal "shelves" at strip
        %       boundaries.  Figure 2 uses Minkowski erosion (polybuffer)
        %       with offset_vertices_raw as fallback.
        %
        %     STEP B — volume-based contraction:
        %       The 3" perpendicular offset alone is not enough volume to
        %       account for the strip's UHPC requirement on sloped walls
        %       (slope correction makes the integrated UHPC volume larger
        %       than a pure 2D offset would give).  After STEP A the
        %       resulting void polygon is contracted INWARD until its 2D
        %       area equals
        %           A_void_target = A_strip * V_void(i) / V_total(i).
        %       Plot 1 contracts horizontally about the void's x-centroid
        %       (z is preserved — no shelves added).  Figure 2 contracts
        %       uniformly about the void's centroid.  Contraction is
        %       symmetric so both sides move toward the centerline by
        %       proportional amounts; the 3" offset on the curved hull
        %       surface is untouched.
        %
        %     Strip colour:
        %       wall   → WALL_GRAY (no inner)
        %       solid  → FILL_SOLID_GRAY (no inner)
        %       annular → outer UHPC_GRAY + inner white + 45-degree hatch

            % ── Style constants (match plot_steel_solve.m) ────────────
            FONT_NAME       = 'Times New Roman';
            FONT_SIZE_AXIS  = 12;
            FONT_SIZE_TITLE = 14;
            FONT_SIZE_LABEL = 13;
            FONT_SIZE_ANNOT = 11;
            FONT_SIZE_LEGEND = 11;
            UHPC_GRAY       = [0.74, 0.76, 0.80];
            WALL_GRAY       = [0.45, 0.46, 0.50];
            FILL_SOLID_GRAY = [0.55, 0.55, 0.60];
            VOID_WHITE      = [1.00, 1.00, 1.00];
            HATCH_COLOR     = [0.50, 0.52, 0.58];
            HATCH_SPACING   = 0.06;
            BOUND_COLOR     = [0.08, 0.08, 0.08];
            INNER_COLOR     = [0.30, 0.30, 0.32];
            WATERLINE_COLOR = [0.15, 0.55, 0.95];
            ZFILL_COLOR     = [0.95, 0.55, 0.10];
            WALL_BOUND_COL  = [0.80, 0.15, 0.10];

            N = length(cstr.strip_z_lo);

            % Per-strip thickness (Phase 1b)
            if isfield(cstr, 't_offset_strip') && ~isempty(cstr.t_offset_strip)
                t_strip = cstr.t_offset_strip;
            else
                t_strip = cstr.t_UHPC * ones(N, 1);
            end
            if isfield(cstr, 'strip_is_solid') && ~isempty(cstr.strip_is_solid)
                is_solid = cstr.strip_is_solid(:);
            else
                is_solid = false(N, 1);
            end
            if isfield(cstr, 'strip_is_wall') && ~isempty(cstr.strip_is_wall)
                is_wall = cstr.strip_is_wall(:);
            else
                is_wall = false(N, 1);
            end

            vs       = cstr.vertical_shift;
            z_fill_w = cstr.z_fill + vs;
            wall_z_b = cstr.wall_z_bottom + vs;
            t_min_mm = cstr.t_min * 1000;

            % Smooth world-frame hull silhouette
            try
                prof_body = WEC_Visualization.build_smooth_viz_profile(config);
            catch
                if isfield(config, 'profile') && ~isempty(config.profile)
                    prof_body = config.profile;
                else
                    prof_body = [];
                end
            end
            if isempty(prof_body)
                warning('WEC_Constructable_Hull:visualize:NoProfile', ...
                        'No smooth profile available — falling back to outer contours.');
                prof_body = build_profile_from_outer_contours(cstr);
            end
            px_outer = prof_body(:, 1);
            pz_outer = prof_body(:, 2) + vs;

            %% ═════════════════════════════════════════════════════════
            %%  FIGURE 1 — XZ midplane elevation
            %% ═════════════════════════════════════════════════════════
            fig1 = figure('Name', 'Constructability: XZ midplane', ...
                          'Color', 'w', 'Position', [60 60 1500 800]);
            ax1 = axes('Parent', fig1, 'Position', [0.06 0.10 0.62 0.84]);
            hold(ax1, 'on');

            x_min_hull = min(px_outer);
            x_max_hull = max(px_outer);

            % ── Inner offset polygons computed on the FULL profile ────
            % CRITICAL: offset the entire hull profile (NOT a strip-clipped
            % piece).  Offsetting a clipped piece adds the strip's horizontal
            % top/bottom edges to the polygon, and the offset of those fake
            % edges produces the "shelves" we keep seeing at strip
            % boundaries.  Offsetting the full profile gives a uniform
            % perpendicular wall thickness on the actual hull surface.
            %
            % Cache one inner polygon per unique thickness.  Use
            % polyshape/polybuffer when available (handles convex AND
            % concave regions correctly via Minkowski erosion); fall back
            % to plot_steel_solve's offset_vertices_raw otherwise.
            inner_keys   = {};
            inner_polys  = {};
            for i = 1:N
                if is_wall(i) || is_solid(i), continue; end
                ti = t_strip(i);
                if ~isfinite(ti) || ti <= 0, continue; end
                key = sprintf('%.6f', ti);
                if any(strcmp(inner_keys, key)), continue; end
                inner_keys{end+1} = key;                                  %#ok<AGROW>
                inner_polys{end+1} = compute_inner_offset_local( ...
                                          px_outer, pz_outer, ti);        %#ok<AGROW>
            end
            inner_lookup = @(t) inner_polys{find(strcmp(inner_keys, ...
                                                       sprintf('%.6f', t)), 1)};

            % ── Per-strip render ─────────────────────────────────────
            for i = 1:N
                z_lo_w = cstr.strip_z_lo(i) + vs;
                z_hi_w = cstr.strip_z_hi(i) + vs;

                % Clip the FULL outer profile to this strip
                [x_top, z_top]   = clip_halfspace_local(px_outer, pz_outer, z_lo_w, false);
                [x_clip, z_clip] = clip_halfspace_local(x_top, z_top,        z_hi_w, true);
                if length(x_clip) < 3, continue; end

                if is_wall(i)
                    fc = WALL_GRAY;
                elseif is_solid(i)
                    fc = FILL_SOLID_GRAY;
                else
                    fc = UHPC_GRAY;
                end
                patch(ax1, x_clip, z_clip, fc, ...
                      'EdgeColor', 'none', 'FaceAlpha', 1.0, ...
                      'HandleVisibility', 'off');

                % Annular strips: clip the precomputed inner-offset polygon
                % (full hull, perpendicular t) to this strip's z-range,
                % then HORIZONTALLY contract the resulting void polygon so
                % its 2D area equals the volume target
                %     A_void_target = A_strip_XZ * V_void(i)/V_total(i).
                % This accounts for the 3" shell already taking some UHPC
                % volume; the void shrinks symmetrically inward in x to
                % make room for the remaining UHPC required by the strip.
                if ~is_wall(i) && ~is_solid(i) && isfinite(t_strip(i)) && t_strip(i) > 0
                    P_in = inner_lookup(t_strip(i));
                    if ~isempty(P_in) && size(P_in, 1) >= 3
                        [xi_top, zi_top]   = clip_halfspace_local( ...
                                                P_in(:,1), P_in(:,2), z_lo_w, false);
                        [xi_clip, zi_clip] = clip_halfspace_local( ...
                                                xi_top, zi_top, z_hi_w, true);
                        if length(xi_clip) >= 3
                            % Volume-based target void area in this strip
                            V_u = 0; V_v = 0;
                            if ~isempty(cstr.strip_V_UHPC), V_u = cstr.strip_V_UHPC(i); end
                            if ~isempty(cstr.strip_V_void), V_v = cstr.strip_V_void(i); end
                            V_t = V_u + V_v;
                            A_strip_XZ = polyarea(x_clip, z_clip);
                            A_void_offset = polyarea(xi_clip, zi_clip);
                            if V_t > 1e-12 && A_strip_XZ > 1e-9 && A_void_offset > 1e-9
                                A_void_target = A_strip_XZ * (V_v / V_t);
                                if A_void_offset > A_void_target
                                    % Horizontal scale toward x-centroid
                                    s_h = A_void_target / A_void_offset;
                                    cx  = mean(xi_clip);
                                    xi_clip = (xi_clip - cx) * s_h + cx;
                                end
                            end
                            patch(ax1, xi_clip, zi_clip, VOID_WHITE, ...
                                  'EdgeColor', INNER_COLOR, ...
                                  'LineStyle', '--', 'LineWidth', 1.0, ...
                                  'HandleVisibility', 'off');
                            draw_hatch_local(ax1, xi_clip, zi_clip, ...
                                             HATCH_SPACING, HATCH_COLOR);
                        end
                    end
                end
            end

            % Outer hull outline (over the patches)
            plot(ax1, [px_outer; px_outer(1)], [pz_outer; pz_outer(1)], ...
                 '-', 'Color', BOUND_COLOR, 'LineWidth', 1.8, ...
                 'HandleVisibility', 'off');

            % Reference lines and strip boundaries
            x_range = [x_min_hull - 0.30, x_max_hull + 0.30];
            for i = 1:N
                z_b = cstr.strip_z_lo(i) + vs;
                line(ax1, x_range, [z_b z_b], 'Color', [0.55 0.55 0.55], ...
                     'LineStyle', ':', 'LineWidth', 0.6, ...
                     'HandleVisibility', 'off');
            end
            line(ax1, x_range, ...
                 [cstr.strip_z_hi(N) + vs, cstr.strip_z_hi(N) + vs], ...
                 'Color', [0.55 0.55 0.55], 'LineStyle', ':', 'LineWidth', 0.6, ...
                 'HandleVisibility', 'off');

            h_wl = plot(ax1, x_range, [0 0], '--', ...
                        'Color', WATERLINE_COLOR, 'LineWidth', 2.0);
            h_wb = plot(ax1, x_range, [wall_z_b wall_z_b], '-', ...
                        'Color', WALL_BOUND_COL, 'LineWidth', 1.6);
            h_zf = plot(ax1, x_range, [z_fill_w z_fill_w], '-.', ...
                        'Color', ZFILL_COLOR, 'LineWidth', 1.8);

            % Legend proxies
            h_uhpc = patch(ax1, NaN, NaN, UHPC_GRAY,       'EdgeColor', 'none');
            h_fs   = patch(ax1, NaN, NaN, FILL_SOLID_GRAY, 'EdgeColor', 'none');
            h_wall = patch(ax1, NaN, NaN, WALL_GRAY,       'EdgeColor', 'none');
            h_void = patch(ax1, NaN, NaN, VOID_WHITE, ...
                           'EdgeColor', INNER_COLOR, 'LineStyle', '--', 'LineWidth', 1.0);

            % Per-strip annotation (right of the hull)
            x_annot = x_range(2) + 0.10;
            for i = 1:N
                z_mid = 0.5 * (cstr.strip_z_lo(i) + cstr.strip_z_hi(i)) + vs;
                if is_wall(i)
                    lbl = sprintf('S%d: 100%% UHPC (wall)', i);
                elseif is_solid(i)
                    lbl = sprintf('S%d: 100%% UHPC (solid fill)', i);
                else
                    V_u = 0; V_v = 0;
                    if ~isempty(cstr.strip_V_UHPC), V_u = cstr.strip_V_UHPC(i); end
                    if ~isempty(cstr.strip_V_void), V_v = cstr.strip_V_void(i); end
                    V_t = V_u + V_v;
                    if V_t > 1e-12
                        pct_u = 100 * V_u / V_t;
                        pct_v = 100 * V_v / V_t;
                    else
                        pct_u = 100; pct_v = 0;
                    end
                    if isfinite(t_strip(i))
                        lbl = sprintf('S%d: %.0f%% UHPC / %.0f%% Void   t = %.2f in', ...
                                      i, pct_u, pct_v, t_strip(i)/0.0254);
                    else
                        lbl = sprintf('S%d: %.0f%% UHPC / %.0f%% Void   t = solid', ...
                                      i, pct_u, pct_v);
                    end
                end
                text(ax1, x_annot, z_mid, lbl, ...
                     'FontSize', FONT_SIZE_ANNOT, 'FontName', FONT_NAME, ...
                     'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                     'Color', [0.15 0.15 0.15], 'Interpreter', 'none');
            end

            xlabel(ax1, 'x  [m]', 'FontSize', FONT_SIZE_LABEL, 'FontName', FONT_NAME);
            ylabel(ax1, 'z  [m]   (world frame, z = 0 = waterline)', ...
                   'FontSize', FONT_SIZE_LABEL, 'FontName', FONT_NAME);
            title(ax1, sprintf(['Constructable Hull (XZ midplane)  |  ' ...
                                't_min = %.0f mm  |  void = perpendicular ' ...
                                'offset, contracted horizontally to match V_void/V_total'], ...
                                t_min_mm), ...
                  'FontSize', FONT_SIZE_TITLE, 'FontName', FONT_NAME, ...
                  'Interpreter', 'none');

            legend(ax1, [h_uhpc, h_fs, h_wall, h_void, h_wl, h_zf, h_wb], ...
                   {'UHPC annular shell', 'UHPC solid fill', ...
                    'UHPC structural wall', 'Void (offset, x-contracted)', ...
                    'Waterline', 'z_{fill}', 'Wall boundary'}, ...
                   'Location', 'eastoutside', 'FontSize', FONT_SIZE_LEGEND, ...
                   'FontName', FONT_NAME, 'Interpreter', 'tex');
            axis(ax1, 'equal');
            grid(ax1, 'on');
            set(ax1, 'FontSize', FONT_SIZE_AXIS, 'FontName', FONT_NAME, ...
                     'LineWidth', 0.6, 'Box', 'on', 'TickDir', 'out', ...
                     'GridLineStyle', ':', 'GridAlpha', 0.3);
            xlim(ax1, [x_range(1), x_annot + 4.5]);

            % Routed through the shared exporter: Plots/ directory, 300 dpi.
            WEC_Visualization.save_figure(fig1, 'WEC_Constructability_XZ', ...
                                          struct('timestamp', false));

            %% ═════════════════════════════════════════════════════════
            %%  FIGURE 2 — Per-strip plan-view cross-sections
            %%  Plot BOTH the bottom (z = strip_z_lo) AND the top
            %%  (z = strip_z_hi) of every strip.  Strips whose end is
            %%  degenerate (e.g. the keel point at strip 1's bottom or
            %%  the hull peak at the wall strip's top) yield empty / tiny
            %%  contours and are skipped automatically.
            %%  Void is drawn as solid white (no hatching).
            %% ═════════════════════════════════════════════════════════

            % Build the (strip_idx, end_label, contour, z) panel list.
            panels = struct('i', {}, 'end_label', {}, 'contour', {}, 'z', {});
            min_area_keep = 1e-4;   % m^2 — drop degenerate end-cross-sections
            for i = 1:N
                conts = cstr.contours_outer{i};
                if isempty(conts), continue; end
                bot_c = conts{1};
                if ~isempty(bot_c) && size(bot_c, 1) >= 3 && ...
                        polyarea(bot_c(:,1), bot_c(:,2)) > min_area_keep
                    panels(end+1) = struct( ...
                        'i', i, 'end_label', 'Bottom', ...
                        'contour', bot_c, 'z', cstr.strip_z_lo(i)); %#ok<AGROW>
                end
                top_c = conts{end};
                if ~isempty(top_c) && size(top_c, 1) >= 3 && ...
                        polyarea(top_c(:,1), top_c(:,2)) > min_area_keep
                    panels(end+1) = struct( ...
                        'i', i, 'end_label', 'Top', ...
                        'contour', top_c, 'z', cstr.strip_z_hi(i)); %#ok<AGROW>
                end
            end

            n_panels = length(panels);
            if n_panels == 0
                fprintf('      Figure 2 skipped: no usable strip cross-sections.\n');
            else
                % Near-square grid.  The old fixed 2-row layout degenerated into
                % a very wide strip of postage-stamp axes once a design had more
                % than ~8 cross-sections.
                n_cols = max(1, ceil(sqrt(n_panels)));
                n_rows = ceil(n_panels / n_cols);
                fig2 = figure('Name', 'Constructability: Strip Plan-View', ...
                              'Color', 'w', ...
                              'Position', [80 80 ...
                                           min(max(n_cols*280, 600), 1800), ...
                                           min(max(n_rows*300, 400), 1000)]);

                for p = 1:n_panels
                    ax = subplot(n_rows, n_cols, p, 'Parent', fig2);
                    hold(ax, 'on');
                    info  = panels(p);
                    i     = info.i;
                    outer = info.contour;

                    xo = outer(:, 1);  yo = outer(:, 2);
                    if abs(xo(end)-xo(1)) > 1e-10 || abs(yo(end)-yo(1)) > 1e-10
                        xo = [xo; xo(1)];  yo = [yo; yo(1)];
                    end

                    if is_wall(i)
                        fc = WALL_GRAY;
                    elseif is_solid(i)
                        fc = FILL_SOLID_GRAY;
                    else
                        fc = UHPC_GRAY;
                    end
                    patch(ax, xo, yo, fc, ...
                          'EdgeColor', 'none', 'FaceAlpha', 1.0);

                    % Void overlay — perpendicular offset by t_strip(i),
                    % then 2D-contracted about its centroid until its
                    % area matches A_outer_plan * V_void(i)/V_total(i).
                    % Filled SOLID WHITE (no hatching).
                    if ~is_wall(i) && ~is_solid(i) && ...
                            isfinite(t_strip(i)) && t_strip(i) > 0
                        P_in = compute_inner_offset_local(xo, yo, t_strip(i));
                        if ~isempty(P_in) && size(P_in, 1) >= 3
                            xi = P_in(:,1);  yi = P_in(:,2);
                            A_outer_plan = polyarea(xo, yo);
                            V_u = 0; V_v = 0;
                            if ~isempty(cstr.strip_V_UHPC), V_u = cstr.strip_V_UHPC(i); end
                            if ~isempty(cstr.strip_V_void), V_v = cstr.strip_V_void(i); end
                            V_t = V_u + V_v;
                            A_in_now = polyarea(xi, yi);
                            if V_t > 1e-12 && A_outer_plan > 1e-9 && A_in_now > 1e-9
                                A_void_target = A_outer_plan * (V_v / V_t);
                                if A_in_now > A_void_target
                                    s_p = sqrt(A_void_target / A_in_now);
                                    cx = mean(xi);  cy = mean(yi);
                                    xi = (xi - cx) * s_p + cx;
                                    yi = (yi - cy) * s_p + cy;
                                end
                            end
                            patch(ax, xi, yi, VOID_WHITE, ...
                                  'EdgeColor', INNER_COLOR, ...
                                  'LineStyle', '--', 'LineWidth', 1.0);
                        end
                    end

                    plot(ax, xo, yo, '-', 'Color', BOUND_COLOR, 'LineWidth', 1.6);
                    plot(ax, 0, 0, '+', 'Color', [0.6 0.6 0.6], ...
                         'MarkerSize', 8, 'LineWidth', 0.8);

                    % Volume-based percentages (same as Plot 1 annotation)
                    V_u = 0; V_v = 0;
                    if ~isempty(cstr.strip_V_UHPC), V_u = cstr.strip_V_UHPC(i); end
                    if ~isempty(cstr.strip_V_void), V_v = cstr.strip_V_void(i); end
                    V_t = V_u + V_v;
                    if V_t > 1e-12
                        pct_u = 100 * V_u / V_t;
                        pct_v = 100 * V_v / V_t;
                    else
                        pct_u = 100; pct_v = 0;
                    end

                    if is_wall(i)
                        ttl = sprintf('Strip %d (%s, z=%.2f m): WALL', ...
                                      i, info.end_label, info.z);
                    elseif is_solid(i)
                        ttl = sprintf('Strip %d (%s, z=%.2f m): SOLID FILL', ...
                                      i, info.end_label, info.z);
                    else
                        if isfinite(t_strip(i))
                            t_str = sprintf('%.2f in', t_strip(i)/0.0254);
                        else
                            t_str = 'solid';
                        end
                        ttl = sprintf(['Strip %d (%s, z=%.2f m): ' ...
                                       '%.0f%% UHPC / %.0f%% Void   t = %s'], ...
                                       i, info.end_label, info.z, ...
                                       pct_u, pct_v, t_str);
                    end
                    title(ax, ttl, 'FontSize', 9, 'FontName', FONT_NAME, ...
                          'Interpreter', 'none');

                    axis(ax, 'equal');
                    grid(ax, 'on');
                    set(ax, 'FontSize', 8, 'FontName', FONT_NAME, ...
                            'LineWidth', 0.5, 'Box', 'on', 'TickDir', 'out');
                    xlabel(ax, 'x [m]', 'FontSize', 8, 'FontName', FONT_NAME);
                    ylabel(ax, 'y [m]', 'FontSize', 8, 'FontName', FONT_NAME);
                    hold(ax, 'off');
                end
            end

            if n_panels > 0
                annotation(fig2, 'textbox', [0.0 0.95 1.0 0.05], ...
                           'String', sprintf(['Plan-View Cross-Sections (bottom + top of each strip)  |  ' ...
                                              'void = inner offset by t_strip(i), ' ...
                                              'contracted to match V_void/V_total  |  ' ...
                                              '%d cross-sections'], n_panels), ...
                           'HorizontalAlignment', 'center', ...
                           'FontSize', FONT_SIZE_TITLE, 'FontName', FONT_NAME, ...
                           'FontWeight', 'bold', 'EdgeColor', 'none', ...
                           'FitBoxToText', 'off');

                WEC_Visualization.save_figure(fig2, 'WEC_Constructability_Strips', ...
                                              struct('timestamp', false));
            end
        end


        %% ═════════════════════════════════════════════════════════════
        %%  HATCH HELPER (kept for backward compatibility)
        %% ═════════════════════════════════════════════════════════════

        function draw_hatch(ax, xp, yp, spacing, color)
        % DRAW_HATCH  Public wrapper around the local hatcher.
            draw_hatch_local(ax, xp, yp, spacing, color);
        end




        %% ═════════════════════════════════════════════════════════════
        %%  CYLINDER VALIDATION
        %% ═════════════════════════════════════════════════════════════

        function report = validate_cylinder(R, H, rho_eff, rho_hull, rho_fill, N_strips, t_min)
        % VALIDATE_CYLINDER  Analytical validation of the constructability
        %   module against a closed-form cylindrical hull.
        %
        %   report = WEC_Constructable_Hull.validate_cylinder( ...
        %       R, H, rho_eff, rho_hull, rho_fill, N_strips, t_min)
        %
        %   Creates a cylindrical mesh with radius R and height H, assigns
        %   uniform rho_eff to all strips, and runs the realization.
        %   Compares against analytical:
        %
        %     s = sqrt( (rho_hull - rho_eff) / (rho_hull - rho_fill) )
        %     R_inner = s × R
        %     V_UHPC  = pi × (R² - R_inner²) × H
        %     V_void  = pi × R_inner² × H
        %     mass    = V_UHPC × rho_hull + V_void × rho_fill
        %     Iyy_area_annulus = pi/4 × (R⁴ - R_inner⁴)   (per unit height)
        %     t_min_actual = (1 - s) × R
        %
        %   Pass criterion: all quantities within 3%.
        %
        %   INPUTS
        %     R        : [m]      cylinder outer radius
        %     H        : [m]      cylinder height
        %     rho_eff  : [kg/m³]  uniform effective density for all strips
        %     rho_hull : [kg/m³]  UHPC density
        %     rho_fill : [kg/m³]  void fill density (air)
        %     N_strips : [-]      number of strips (default 5)
        %     t_min    : [m]      min wall thickness (default 0.0762)
        %
        %   OUTPUT
        %     report : struct with exact, computed, pct_errors, passed

            if nargin < 6 || isempty(N_strips), N_strips = 5; end
            if nargin < 7 || isempty(t_min), t_min = 0.0762; end

            % NOTE (v4.0): This validation method creates an internal
            %   triangulated cylinder mesh and builds a config struct with
            %   config.full_mesh.  Since realize() now expects config.ms2_model
            %   (a WEC_MS2_Parser object), this method is temporarily disabled.
            %   It will be re-implemented when the constructability module is
            %   validated against the C0 geometry using the MS2 parser.
            warning('WEC_Constructable_Hull:ValidateCylinderDisabled', ...
                    'validate_cylinder is disabled pending MS2 integration.');
            report = struct('passed', false, 'note', 'Disabled pending MS2 integration');
            return;

            fprintf('\n═══ CYLINDER VALIDATION (Constructable Hull) ═══\n');
            fprintf('  R = %.3f m, H = %.3f m, rho_eff = %.0f kg/m³\n', R, H, rho_eff);
            fprintf('  rho_hull = %.0f, rho_fill = %.1f, N_strips = %d\n', ...
                    rho_hull, rho_fill, N_strips);

            %% ── Create cylindrical mesh ──────────────────────────
            %  Same method as WEC_Shell_Offset.validate_cylinder.

            n_circ  = 60;
            n_axial = 40;
            [theta_grid, z_grid] = meshgrid( ...
                linspace(0, 2*pi, n_circ+1), ...
                linspace(0, H, n_axial+1));

            x_cyl = R * cos(theta_grid);
            y_cyl = R * sin(theta_grid);
            z_cyl = z_grid;

            fv = surf2patch(x_cyl, y_cyl, z_cyl, 'triangles');

            % Add top and bottom caps
            theta_cap = linspace(0, 2*pi, n_circ+1)';
            theta_cap = theta_cap(1:end-1);
            x_cap = R * cos(theta_cap);
            y_cap = R * sin(theta_cap);

            n_ex = size(fv.vertices, 1);

            % Bottom cap (z = 0)
            v_bot_center = [0, 0, 0];
            v_bot_ring   = [x_cap, y_cap, zeros(n_circ, 1)];
            fv.vertices  = [fv.vertices; v_bot_center; v_bot_ring];
            ic_bot = n_ex + 1;
            for k = 1:n_circ
                kn = mod(k, n_circ) + 1;
                fv.faces(end+1, :) = [ic_bot, ic_bot + k, ic_bot + kn]; %#ok<AGROW>
            end

            n_ex2 = size(fv.vertices, 1);

            % Top cap (z = H)
            v_top_center = [0, 0, H];
            v_top_ring   = [x_cap, y_cap, H * ones(n_circ, 1)];
            fv.vertices  = [fv.vertices; v_top_center; v_top_ring];
            ic_top = n_ex2 + 1;
            for k = 1:n_circ
                kn = mod(k, n_circ) + 1;
                fv.faces(end+1, :) = [ic_top, ic_top + kn, ic_top + k]; %#ok<AGROW>
            end

            cyl_mesh = triangulation(double(fv.faces), double(fv.vertices));

            %% ── Build minimal config for realize() ───────────────

            density_nodes_z = linspace(0, H, N_strips)';
            densities       = rho_eff * ones(N_strips, 1);

            % Topology for volume check
            F  = cyl_mesh.ConnectivityList;
            V  = cyl_mesh.Points;
            p1 = V(F(:,1), :);
            p2 = V(F(:,2), :);
            p3 = V(F(:,3), :);
            normals_mesh = cross(p2 - p1, p3 - p1, 2);
            face_vol_contribs = dot(p1, normals_mesh, 2) / 6;

            config_cyl = struct();
            config_cyl.full_mesh              = cyl_mesh;
            config_cyl.density_nodes_z        = density_nodes_z;
            config_cyl.topology.face_volume_contribs = face_vol_contribs;
            config_cyl.constructability_rho_hull     = rho_hull;
            config_cyl.constructability_rho_fill     = rho_fill;
            config_cyl.constructability_t_min        = t_min;
            config_cyl.constructability_wall_height  = 0;   % no wall for test
            config_cyl.constructability_n_sub        = 20;

            x_opt = [0; densities];   % draft = 0

            % Minimal final_props for verification
            final_props_cyl.mass_total    = rho_eff * pi * R^2 * H;
            final_props_cyl.CG_total      = [0, 0, H/2];
            final_props_cyl.Iyy           = rho_eff * pi * R^2 * H * (R^2/4 + H^2/12);
            final_props_cyl.KM            = H/2;  % placeholder
            final_props_cyl.GM_L          = 0;    % placeholder
            final_props_cyl.CB            = [0, 0, H/2];

            %% ── Run realize() ────────────────────────────────────

            cstr = WEC_Constructable_Hull.realize(config_cyl, x_opt, final_props_cyl);

            %% ── Analytical solution ──────────────────────────────

            s = sqrt((rho_hull - rho_eff) / (rho_hull - rho_fill));
            R_in = s * R;

            exact.s          = s;
            exact.R_inner    = R_in;
            exact.V_UHPC     = pi * (R^2 - R_in^2) * H;
            exact.V_void     = pi * R_in^2 * H;
            exact.V_total    = pi * R^2 * H;
            exact.mass       = exact.V_UHPC * rho_hull + exact.V_void * rho_fill;
            exact.Iyy_area_ann = pi/4 * (R^4 - R_in^4) * H;   % integrated over height
            exact.t_actual   = (1 - s) * R;

            %% ── Comparison ───────────────────────────────────────

            computed.V_UHPC  = cstr.total_V_UHPC;
            computed.V_void  = cstr.total_V_void;
            computed.V_total = cstr.total_V_hull;
            computed.mass    = cstr.total_mass;

            pct = @(c, e) 100 * (c - e) / max(abs(e), eps);

            pct_errors.V_UHPC  = pct(computed.V_UHPC,  exact.V_UHPC);
            pct_errors.V_void  = pct(computed.V_void,   exact.V_void);
            pct_errors.V_total = pct(computed.V_total,  exact.V_total);
            pct_errors.mass    = pct(computed.mass,     exact.mass);

            fprintf('\n  %-20s %12s %12s %8s\n', 'Quantity', 'Exact', 'Computed', 'Error%');
            fprintf('  %s\n', repmat('─', 1, 56));
            fprintf('  %-20s %12.6f %12.6f %+8.2f%%\n', ...
                    'V_UHPC [m³]', exact.V_UHPC, computed.V_UHPC, pct_errors.V_UHPC);
            fprintf('  %-20s %12.6f %12.6f %+8.2f%%\n', ...
                    'V_void [m³]', exact.V_void, computed.V_void, pct_errors.V_void);
            fprintf('  %-20s %12.6f %12.6f %+8.2f%%\n', ...
                    'V_total [m³]', exact.V_total, computed.V_total, pct_errors.V_total);
            fprintf('  %-20s %12.1f %12.1f %+8.2f%%\n', ...
                    'Mass [kg]', exact.mass, computed.mass, pct_errors.mass);
            fprintf('  %-20s %12.4f %12s %8s\n', ...
                    's_exact [-]', exact.s, '—', '—');
            fprintf('  %-20s %12.4f %12s %8s\n', ...
                    't_actual [m]', exact.t_actual, '—', '—');

            max_err = max(abs([pct_errors.V_UHPC, pct_errors.V_void, ...
                               pct_errors.V_total, pct_errors.mass]));

            if max_err < 3.0
                fprintf('\n  ✓ PASSED: max error = %.2f%% (< 3%%)\n', max_err);
            else
                fprintf('\n  ✗ FAILED: max error = %.2f%% (≥ 3%%)\n', max_err);
            end
            fprintf('═════════════════════════════════════════════════\n\n');

            report.exact      = exact;
            report.computed   = computed;
            report.pct_errors = pct_errors;
            report.max_error  = max_err;
            report.passed     = (max_err < 3.0);
            report.cstr       = cstr;
        end


        %% ═════════════════════════════════════════════════════════════
        %%  BUILD REALISED FINAL_PROPS (UHPC+void)
        %% ═════════════════════════════════════════════════════════════

        function realised = build_realised_props(final_props, cstr, config)
        % BUILD_REALISED_PROPS  Drop-in replacement for final_props that
        %   carries the AS-BUILT constructable-hull (UHPC + void) properties.
        %
        %   realised = WEC_Constructable_Hull.build_realised_props( ...
        %                  final_props, cstr, config)
        %
        %   PHYSICS:
        %     Phase 1a (global solve) finds (t_UHPC*, z_fill*) under a
        %     wall-aware mass model.  Phase 1b (per-strip thickening)
        %     locally raises t_offset_strip(i) until joint mass+GM residuals
        %     close.  The realised configuration therefore differs from
        %     the optimiser's per-strip rho_eff distribution; mass, CG, GM
        %     and pitch period are RE-SOLVED, NOT carried over.
        %
        %   FIELD MAP:
        %     - vertical_shift / draft        ← steel_data (re-solved via
        %                                       solve_draft_for_mass to
        %                                       enforce mass balance at
        %                                       the realised partition)
        %     - V_sub / Aw / I_wp_yy / I_wp_xx / CB / A_sub / KM /
        %       mass_buoyant_force            ← interpolated at re-solved draft
        %     - mass_total / CG_total         ← from steel_data (Phase 1+1b)
        %     - Inertia_Tensor (Ixx,Iyy,Izz)  ← realised radial distribution
        %     - GM_L                          ← KM − realised CG_z
        %     - K_hydro / K_pto / K_total     ← K_hydro recomputed; K_pto copied
        %     - A_full / B_full / A11/33/55   ← WAMIT @ re-solved draft, retransformed to realised CG
        %     - periods.heave/pitch/surge     ← recomputed with realised mass + Iyy + added mass
        %     - coupled_periods/modes/PF      ← re-eig with realised mass tensor
        %     - MassMatrix_CG/Origin          ← rebuilt
        %     - realised_strips substruct     ← per-strip realised partition
        %     - components                    ← REPLACED with realised
        %                                       per-strip UHPC/void breakdown
        %     - densities_at_nodes/cross_section ← carried from optimiser

            % cstr originates from WEC_Shell_Offset.solve_constructable() via
            % extract_strip_geometry(), so it carries all steel_data fields
            % required by Shell_Offset's version PLUS per-strip realisation
            % arrays consumed by build_realised_props.
            realised = WEC_Shell_Offset.build_realised_props(final_props, cstr, config);
            realised.realisation_mode = 'uhpc_fill';
            if isfield(realised, 'realised_strips') && isstruct(realised.realised_strips)
                realised.realised_strips.realisation_mode = 'uhpc_fill';
            end
        end


        function PF = compute_participation_factors(modes, M)
            try
                n_modes = size(modes, 2);
                PF = zeros(3, n_modes);
                for j = 1:n_modes
                    phi = modes(:, j);
                    modal_mass = phi' * M * phi;
                    if abs(modal_mass) > 1e-12
                        for i = 1:3
                            PF(i, j) = (phi(i)^2 * M(i, i)) / modal_mass * 100;
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

        function [xc, zc] = clip_z(x, z, z_cut, side)
        % CLIP_Z  Sutherland–Hodgman clip of a 2D polygon against a horizontal line.
        %
        %   [xc, zc] = WEC_Constructable_Hull.clip_z(x, z, z_cut, side)
        %
        %   side = 'below'  → keep vertices with z <= z_cut
        %   side = 'above'  → keep vertices with z >= z_cut
        %
        %   Identical algorithm used in plot_steel_solve.clip_halfspace.
        %   Returns empty arrays when the input has < 3 points.
            keep_below = strcmp(side, 'below');
            n = length(x);
            if n < 3,  xc = [];  zc = [];  return;  end
            xc = zeros(2*n, 1);
            zc = zeros(2*n, 1);
            cnt = 0;
            for i = 1:n
                j  = mod(i, n) + 1;
                zi = z(i);  zj = z(j);
                if keep_below
                    in_i = (zi <= z_cut);
                    in_j = (zj <= z_cut);
                else
                    in_i = (zi >= z_cut);
                    in_j = (zj >= z_cut);
                end
                if in_i
                    cnt = cnt + 1;
                    xc(cnt) = x(i);  zc(cnt) = zi;
                end
                if in_i ~= in_j && abs(zj - zi) > 1e-14
                    t    = (z_cut - zi) / (zj - zi);
                    cnt  = cnt + 1;
                    xc(cnt) = x(i) + t * (x(j) - x(i));
                    zc(cnt) = z_cut;
                end
            end
            xc = xc(1:cnt);
            zc = zc(1:cnt);
        end

    end  % methods (Static)
end  % classdef


%% ═════════════════════════════════════════════════════════════════════
%%  MODULE-LEVEL HELPER (outside classdef — accessible by all methods)
%% ═════════════════════════════════════════════════════════════════════

function r = ternary(cond, t, f)
% TERNARY  Inline if-else for fprintf convenience.
%   Matches the helper used in WEC_Driver, WEC_Configuration_Builder,
%   and WEC_Main_Optimizer.
    if cond, r = t; else, r = f; end
end


function R = find_runs(mask)
% FIND_RUNS  Return [start, stop] indices for each maximal run of true in mask.
    mask = logical(mask(:));
    if isempty(mask)
        R = zeros(0, 2);
        return;
    end
    d = diff([false; mask; false]);
    starts = find(d == 1);
    stops  = find(d == -1) - 1;
    n = min(length(starts), length(stops));
    R = [starts(1:n), stops(1:n)];
end


function [xc, zc] = clip_halfspace_local(x, z, z_cut, keep_below)
% CLIP_HALFSPACE_LOCAL  Sutherland-Hodgman polygon clip against z = z_cut.
%   keep_below = true  -> retain z <= z_cut
%   keep_below = false -> retain z >= z_cut
%   Mirrors plot_steel_solve.clip_halfspace exactly.
    n = length(x);
    if n < 3
        xc = [];  zc = [];
        return;
    end
    xc = zeros(2*n, 1);
    zc = zeros(2*n, 1);
    cnt = 0;
    for i = 1:n
        j = mod(i, n) + 1;
        zi = z(i);  zj = z(j);
        if keep_below
            in_i = (zi <= z_cut);
            in_j = (zj <= z_cut);
        else
            in_i = (zi >= z_cut);
            in_j = (zj >= z_cut);
        end
        if in_i
            cnt = cnt + 1;
            xc(cnt) = x(i);
            zc(cnt) = zi;
        end
        if in_i ~= in_j && abs(zj - zi) > 1e-14
            t = (z_cut - zi) / (zj - zi);
            cnt = cnt + 1;
            xc(cnt) = x(i) + t * (x(j) - x(i));
            zc(cnt) = z_cut;
        end
    end
    xc = xc(1:cnt);
    zc = zc(1:cnt);
end


function draw_hatch_local(ax, xp, yp, spacing, color)
% DRAW_HATCH_LOCAL  Fill polygon (xp,yp) with 45-degree hatching.
%   Mirrors plot_steel_solve.draw_hatch.
    if length(xp) < 3, return; end
    x_lo = min(xp);  x_hi = max(xp);
    y_lo = min(yp);  y_hi = max(yp);
    if (x_hi - x_lo) < 1e-9 || (y_hi - y_lo) < 1e-9, return; end
    c_min = x_lo - y_hi;
    c_max = x_hi - y_lo;
    c_vals = c_min:spacing:c_max;
    % One NaN-separated polyline for the whole hatch instead of one line object
    % per segment — the per-segment version created hundreds of objects per
    % strip, which is what made this figure slow to draw and heavy to export.
    hx = [];  hy = [];
    for ci = 1:length(c_vals)
        c = c_vals(ci);
        y_line = linspace(y_lo, y_hi, 200)';
        x_line = y_line + c;
        in = inpolygon(x_line, y_line, xp, yp);
        if ~any(in), continue; end
        d = diff([false; in; false]);
        starts = find(d == 1);
        stops  = find(d == -1) - 1;
        n_seg = min(length(starts), length(stops));
        for s = 1:n_seg
            idx = starts(s):stops(s);
            hx = [hx; x_line(idx); NaN];   %#ok<AGROW>
            hy = [hy; y_line(idx); NaN];   %#ok<AGROW>
        end
    end
    if ~isempty(hx)
        plot(ax, hx, hy, '-', 'Color', color, 'LineWidth', 0.4, ...
             'HandleVisibility', 'off');
    end
end


function P_in = compute_inner_offset_local(x_outer, y_outer, t)
% COMPUTE_INNER_OFFSET_LOCAL  Inner offset of polygon (x_outer, y_outer)
%   by perpendicular distance t.
%
%   Method order (best→worst, falls through on failure):
%     1. polyshape + polybuffer(-t)   — Minkowski erosion, the
%        literature-standard for offsetting polygons with mixed convex/
%        concave regions (handles self-intersection correctly).  Available
%        in MATLAB R2018b+.  This is the GEOMETRICALLY CORRECT choice.
%     2. WEC_Shell_Offset.offset_vertices_raw(x, y, t) — vertex-normal
%        offset with miter limit.  Same routine plot_steel_solve.m uses.
%        Robust on convex polygons; can produce small artefacts on sharp
%        concave corners but rarely catastrophic for hull-like shapes.
%
%   Returns Nx2 polygon (closed-implicitly, last vertex != first), or []
%   if the offset polygon collapsed.
    P_in = [];
    if length(x_outer) < 3 || ~isfinite(t) || t <= 0, return; end
    x_outer = x_outer(:);  y_outer = y_outer(:);

    % --- Method 1: polybuffer ---------------------------------------------
    try
        % polyshape silently warns on duplicate / collinear points; suppress.
        ws = warning('off', 'MATLAB:polyshape:repairedBySimplify');
        ps = polyshape(x_outer, y_outer);
        ps_in = polybuffer(ps, -t, 'JointType', 'miter', 'MiterLimit', 10);
        warning(ws);
        if ~isempty(ps_in.Vertices) && area(ps_in) > 1e-10
            % polybuffer may return multiple regions; keep the largest.
            R = regions(ps_in);
            if numel(R) > 1
                [~, kbig] = max(arrayfun(@area, R));
                ps_in = R(kbig);
            end
            V = ps_in.Vertices;
            % polyshape vertices may include NaN row separators; drop them.
            V = V(all(~isnan(V), 2), :);
            if size(V, 1) >= 3
                P_in = V;
                return;
            end
        end
    catch
        % polyshape/polybuffer not available or failed — fall through.
    end

    % --- Method 2: vertex-normal offset (steel's method) ------------------
    try
        [xi, yi] = WEC_Shell_Offset.offset_vertices_raw(x_outer, y_outer, t);
        if length(xi) >= 3
            A_outer = polyarea(x_outer, y_outer);
            A_inner = polyarea(xi, yi);
            if A_inner > 1e-10 && A_inner < A_outer
                P_in = [xi(:), yi(:)];
            end
        end
    catch
    end
end


function prof = build_profile_from_outer_contours(cstr)
% BUILD_PROFILE_FROM_OUTER_CONTOURS  Last-resort fallback when
%   WEC_Visualization.build_smooth_viz_profile is unavailable and
%   config.profile is empty.  Reconstructs a coarse XZ silhouette
%   from cstr.contours_outer by extracting (x_min, x_max) at each
%   z-sample.  Body frame.
    N = length(cstr.strip_z_lo);
    pts_R = [];  pts_L = [];
    for i = 1:N
        ci = cstr.contours_outer{i};
        if isempty(ci), continue; end
        n_z = length(ci);
        z_b = linspace(cstr.strip_z_lo(i), cstr.strip_z_hi(i), n_z)';
        for k = 1:n_z
            pk = ci{k};
            if ~isempty(pk) && size(pk, 1) >= 3
                pts_R(end+1, :) = [max(pk(:,1)), z_b(k)];   %#ok<AGROW>
                pts_L(end+1, :) = [min(pk(:,1)), z_b(k)];   %#ok<AGROW>
            end
        end
    end
    if isempty(pts_R)
        prof = [];
        return;
    end
    [~, ord_R] = sort(pts_R(:,2));
    [~, ord_L] = sort(pts_L(:,2), 'descend');
    prof = [pts_R(ord_R, :); pts_L(ord_L, :)];
end