classdef WEC_Shell_Offset
    % WEC_SHELL_OFFSET  Solid-steel fill-level solver (post-optimisation).
    %
    %   This class is invoked AFTER the Stage-2 fmincon converges, alongside
    %   WEC_Constructable_Hull.realize.  It does not modify the optimiser's
    %   converged solution — it solves an inverse problem on top of it.
    %
    %   PHYSICAL MODEL  (replaces the old composite-shell forward model)
    %     The hull is built as:
    %       (a) A uniform-thickness steel jacket (offset t_steel inward
    %           from the outer hull surface, same Steiner / vertex-normal
    %           machinery the old shell used).
    %       (b) A monolithic interior fill that switches at a single
    %           horizontal cut z_fill:
    %               z < z_fill  →  interior is solid steel
    %                              (cross-section is ALL steel: A_outer)
    %               z > z_fill  →  interior is air
    %                              (cross-section is steel jacket annulus
    %                               + air core: A_jacket + A_inner)
    %
    %   FREE DESIGN VARIABLES
    %     t_steel  [m]   uniform jacket thickness (constant across hull)
    %     z_fill   [m]   solid-steel/air transition height (continuous)
    %
    %   DERIVED VARIABLE  (enforces mass balance — hard constraint)
    %     draft d  [m]   re-solved each candidate so that
    %                    ρ_water · V_sub(d) = M_total(t_steel, z_fill)
    %
    %   SOFT TARGETS  (best-effort L² fit)
    %     final_props.GM_L
    %     final_props.periods.heave
    %     final_props.periods.pitch
    %
    %   METHOD INVENTORY
    %     solve                      — Public entry: post-optimiser inverse solve
    %     build_geometry_grid        — Sample A_outer/A_inner over a z-grid at fixed t_steel
    %     inner_props_at_z           — polygeom on miter polygon: A_inner + Iyy at one z (private)
    %     slope_cos_from_table       — Slope cos from FD on Aw_table (private)
    %     integrate_split            — Trapezoidal split at z_fill into steel/air regions
    %     evaluate_realised          — Hydrostatics, GM, natural periods at given (t,z_fill)
    %     solve_draft_for_mass       — Inner fzero on draft to enforce mass balance
    %     offset_vertices_raw        — Inward vertex-normal offset with miter (kept)
    %
    %   DEPENDENCIES
    %     WEC_HydroProperties.extract_isocurve_at_z    (parametric cross-section, robust)
    %     WEC_HydroProperties.waterplane_properties    (angular-sort + Green's theorem)
    %     WEC_Core_Functions.polygeom                  (offset polygon Iyy)
    %     WEC_Core_Functions.interpolate_wamit_added_mass  (A11/33/55 at draft)
    %     config.Aw_table*, I_wp_yy_table, V_sub_table, CB_z_table, boundary_cache
    %       (all precomputed once by WEC_Configuration_Builder)
    %
    %   NO SILENT FAILURES POLICY
    %     This module errors loudly when any required input is missing or
    %     malformed.  No try/catch zero-fills, no hard-coded fallbacks for
    %     missing config fields.  Hull-boundary degeneracies (cross-section
    %     evaluation returning empty above hull_z_max or below hull_z_min)
    %     ARE expected and handled — but with an explicit count and a
    %     fail-loud check that no more than 5% of interior z-grid points
    %     went degenerate.
    %
    %   See also: WEC_Constructable_Hull, WEC_Configuration_Builder,
    %             calculate_3d_properties, plot_steel_solve
    %
    %   Author:  WEC Optimisation Team
    %   Version: 2.0 — Post-optimisation steel-fill inverse solver

    methods (Static)

        %% ═════════════════════════════════════════════════════════════
        %%  MAIN ENTRY POINT
        %% ═════════════════════════════════════════════════════════════

        function steel_data = solve(config, x_opt_3d, final_props, opts)
        % SOLVE  Constrained fmincon SQP solve for (vs, t_steel, z_fill).
        %
        %   steel_data = WEC_Shell_Offset.solve(config, x_opt_3d, final_props)
        %   steel_data = WEC_Shell_Offset.solve(config, x_opt_3d, final_props, opts)
        %
        %   Called ONCE after Stage-2 fmincon converges (see WEC_Main_Optimizer §4c).
        %
        %   FORMULATION  (mirrors Stage 2 objective_function_3d / constraint_function_3d)
        %     DVs         x = [vs; t_steel; z_fill]
        %     Equality    ceq = M_total(t, z_fill) / (ρ_w · V_sub(-vs)) − 1 = 0
        %                       (buoyancy balance — same as Stage 2 §1179)
        %     Inequality  c  = 1 − GM(x) / gm_min ≤ 0
        %                       (stability floor — same as Stage 2 §1191)
        %     Objective   Φ  = phi(r_heave) + phi(r_pitch)
        %                       r_h = (T_heave - T_heave_goal) / heave_half
        %                       r_p = (T_pitch - T_pitch_goal) / pitch_half
        %                       phi = range_penalty (C1 piecewise quadratic)
        %     Solver      fmincon SQP
        %     Warm start  x0 = [x_opt_3d(1); t_init; analytic z_fill seed]
        %
        %   GM does NOT appear in the objective (user spec: inequality-only).
        %   Mass does NOT appear in the objective (it's the equality constraint).
        %   No per-target weights — range normalization already balances the terms.
        %
        %   INPUTS
        %     config       struct from WEC_Configuration_Builder.  Required:
        %                    ms2_model, hull_z_min, hull_z_max,
        %                    Aw_table_z, Aw_table, V_sub_table, CB_z_table,
        %                    I_wp_yy_table, RHO_WATER, G,
        %                    rho_steel, rho_air, steel_t_init, steel_t_min,
        %                    steel_max_slope_factor, steel_n_z_grid,
        %                    vertical_shift_bounds,
        %                    gm_min, T_heave_goal, T_heave_range,
        %                    T_pitch_goal, T_pitch_range, zone_k_amp.
        %     x_opt_3d     [1+N x 1] converged design vector [vertical_shift; rho_1..rho_N]
        %     final_props  struct from calculate_3d_properties(x_opt_3d, config).
        %                  Required: GM_L, mass_total, periods.heave, periods.pitch.
        %                  Used for warm-start info and pre-flight diagnostics only —
        %                  the realisation no longer tracks final_props.mass_total
        %                  as a soft target.
        %     opts         (optional) struct overriding rho_steel, rho_air,
        %                    t_init, t_min, max_slope_factor, n_z_grid.
        %
        %   OUTPUT
        %     steel_data   struct (see end of this function for full field list).

            t_solve_start = tic;
            fprintf('\n    WEC_Shell_Offset.solve (fmincon SQP, constrained):\n');

            %% §1  INPUT VALIDATION  ──────────────────────────────────
            if nargin < 3
                error('WEC_Shell_Offset:NotEnoughInputs', ...
                      'solve(config, x_opt_3d, final_props[, opts]) requires 3+ args.');
            end
            if nargin < 4 || isempty(opts), opts = struct(); end

            required_cfg = {'ms2_model','hull_z_min','hull_z_max', ...
                            'Aw_table_z','Aw_table','V_sub_table','CB_z_table', ...
                            'I_wp_yy_table','RHO_WATER','G', ...
                            'rho_steel','rho_air','steel_t_init','steel_t_min', ...
                            'steel_max_slope_factor','steel_n_z_grid', ...
                            'vertical_shift_bounds','gm_min', ...
                            'T_heave_goal','T_heave_range', ...
                            'T_pitch_goal','T_pitch_range','zone_k_amp'};
            for k = 1:length(required_cfg)
                if ~isfield(config, required_cfg{k}) || isempty(config.(required_cfg{k}))
                    error('WEC_Shell_Offset:MissingConfig', ...
                          'Required config field "%s" is missing or empty.', ...
                          required_cfg{k});
                end
            end

            x_opt_3d = x_opt_3d(:);
            if length(x_opt_3d) < 2
                error('WEC_Shell_Offset:BadOptVector', ...
                      'x_opt_3d must be [vertical_shift; rho_1; ...] (length >= 2).');
            end

            %% §2  RESOLVE OPTIONS
            rho_steel        = WEC_Shell_Offset.opt_or_cfg(opts, 'rho_steel',        config.rho_steel);
            rho_air          = WEC_Shell_Offset.opt_or_cfg(opts, 'rho_air',          config.rho_air);
            t_init           = WEC_Shell_Offset.opt_or_cfg(opts, 't_init',           config.steel_t_init);
            t_min            = WEC_Shell_Offset.opt_or_cfg(opts, 't_min',            config.steel_t_min);
            max_slope_factor = WEC_Shell_Offset.opt_or_cfg(opts, 'max_slope_factor', config.steel_max_slope_factor);
            n_z_grid         = WEC_Shell_Offset.opt_or_cfg(opts, 'n_z_grid',         config.steel_n_z_grid);

            assert(rho_steel > 0,        'rho_steel must be > 0 (got %g)',        rho_steel);
            assert(rho_air   > 0,        'rho_air must be > 0 (got %g)',          rho_air);
            assert(rho_air   < rho_steel,'rho_air (%g) must be < rho_steel (%g)', rho_air, rho_steel);
            assert(n_z_grid  >= 50,      'n_z_grid must be >= 50 (got %d)',       n_z_grid);
            assert(t_min     >= 0,       't_min must be >= 0 (got %g)',           t_min);

            %% §3  HULL BOUNDS AND THICKNESS BRACKET

            hull_z_min = config.hull_z_min;
            hull_z_max = config.hull_z_max;
            assert(hull_z_max > hull_z_min, 'hull extents inverted (%g >= %g)', ...
                   hull_z_min, hull_z_max);

            % Upper bound on t_steel:  half the minimum hull radius in the
            % middle 80% of the height (avoid the apex/keel-tip pinch points).
            %% §3  HULL BOUNDS AND THICKNESS BRACKET
            n_probe = 9;
            z_probes = linspace(0.1*hull_z_min + 0.9*hull_z_max, ...
                                0.9*hull_z_min + 0.1*hull_z_max, n_probe);
            r_probes = zeros(n_probe, 1);

            if ~isfield(config, 'boundary_cache') || ...
                    ~isstruct(config.boundary_cache) || ...
                    ~isfield(config.boundary_cache, 'sources')
                error('WEC_Shell_Offset:NoBoundaryCache', ...
                      ['config.boundary_cache is missing or malformed. ' ...
                       'Configuration_Builder must populate it via ' ...
                       'WEC_HydroProperties.precompute_boundary_cache.']);
            end
            cache = config.boundary_cache;
            for k = 1:n_probe
                r_probes(k) = WEC_HydroProperties.compute_rmin_at_z( ...
                                  config.ms2_model, z_probes(k), 60, cache);
            end
            r_min_global = min(r_probes(r_probes > 0));
            if isempty(r_min_global) || r_min_global <= 0
                error('WEC_Shell_Offset:NoFiniteRadius', ...
                      'compute_rmin_at_z returned no positive radius across %d probes.', ...
                      n_probe);
            end
            t_max = 0.5 * r_min_global;
            if t_min >= t_max
                error('WEC_Shell_Offset:TMinExceedsTMax', ...
                      ['t_min = %.5f m >= t_max = %.5f m — hull too narrow ', ...
                       'for requested fabrication floor.'], t_min, t_max);
            end

            t_init_requested = t_init;
            t_init = max(min(t_init, 0.9*t_max), max(t_min, 1e-4));
            if abs(t_init - t_init_requested) > 1e-6
                warning('WEC_Shell_Offset:TInitClamped', ...
                        'Initial t_steel guess clamped %.5f → %.5f m.', ...
                        t_init_requested, t_init);
            end

            vs_lb = config.vertical_shift_bounds(1);
            vs_ub = config.vertical_shift_bounds(2);
            vs_opt = x_opt_3d(1);
            draft_opt = abs(hull_z_min + vs_opt);

            fprintf('      Hull z=[%.4f, %.4f] m  t∈[%.5f, %.4f] m  vs∈[%.3f, %.3f] m\n', ...
                    hull_z_min, hull_z_max, t_min, 0.95*t_max, vs_lb, vs_ub);
            fprintf('      Targets (from optimiser): GM_min=%.3f m  T_h=%.3f s  T_p=%.3f s  M=%.1f kg\n', ...
                    config.gm_min, config.T_heave_goal, config.T_pitch_goal, ...
                    final_props.mass_total);
            fprintf('      Warm-start vs (from x_opt_3d): %.4f m  (draft=%.4f m)\n', ...
                    vs_opt, draft_opt);

            %% §4  ANALYTIC z_fill SEED at (vs=vs_opt, t=t_init)
            %  Picks the z_fill that exactly satisfies buoyancy balance for
            %  the warm-start hull mass at the optimiser's draft, so fmincon
            %  starts near feasibility on ceq.
            V_sub_at_vs_opt = max(0, interp1(config.Aw_table_z, ...
                                              config.V_sub_table, ...
                                              -vs_opt, 'linear', 0));
            M_buoy_target = config.RHO_WATER * V_sub_at_vs_opt;

            [grids_seed, n_zero_seed] = WEC_Shell_Offset.build_geometry_grid( ...
                config, t_init, n_z_grid, max_slope_factor);
            if n_zero_seed > 0.05 * n_z_grid
                error('WEC_Shell_Offset:DegenerateSeedGeometry', ...
                      'Seed geometry grid has %d/%d (>5%%) degenerate sections.', ...
                      n_zero_seed, n_z_grid);
            end
            z_grid_seed = grids_seed.z;
            A_o_seed    = grids_seed.A_outer;
            A_i_seed    = grids_seed.A_inner;
            V_jacket    = trapz(z_grid_seed, A_o_seed - A_i_seed);
            V_inner     = trapz(z_grid_seed, A_i_seed);
            V_inner_below_target = (M_buoy_target - rho_steel*V_jacket - ...
                                    rho_air*V_inner) / (rho_steel - rho_air);

            if V_inner_below_target < 0
                z_fill_seed = hull_z_min + 1e-3;
            elseif V_inner_below_target > V_inner
                z_fill_seed = hull_z_max - 1e-3;
            else
                cumV = cumtrapz(z_grid_seed, A_i_seed);
                interior = A_i_seed > 1e-10;
                if sum(interior) < 2
                    z_fill_seed = 0.5 * (hull_z_min + hull_z_max);
                else
                    [V_uniq, ia_uniq] = unique(cumV(interior), 'stable');
                    z_int = z_grid_seed(interior);
                    z_fill_seed = interp1(V_uniq, z_int(ia_uniq), ...
                                          V_inner_below_target, 'linear', 'extrap');
                end
                z_fill_seed = max(hull_z_min + 1e-3, ...
                                  min(hull_z_max - 1e-3, z_fill_seed));
            end
            fprintf('      Warm-start z_fill seed: %.4f m (analytic, M_buoy=%.0f kg)\n', ...
                    z_fill_seed, M_buoy_target);

            %% §5  CONTEXT (shared by nested objective + constraint)
            ctx = struct();
            ctx.config           = config;
            ctx.rho_steel        = rho_steel;
            ctx.rho_air          = rho_air;
            ctx.max_slope_factor = max_slope_factor;
            ctx.n_z_grid         = n_z_grid;
            ctx.hull_z_min       = hull_z_min;
            ctx.hull_z_max       = hull_z_max;

            %% §6  fmincon SQP
            %  DVs:        x = [vs; t_steel; z_fill]
            %  Bounds:     lb / ub on each
            %  Equality:   ceq = M_total / (ρ_w·V_sub(-vs)) − 1 = 0
            %  Inequality: c   = 1 − GM/gm_min ≤ 0
            %  Objective:  Φ = phi(r_heave) + phi(r_pitch)

            lb = [vs_lb;     t_min;          hull_z_min + 1e-3];
            ub = [vs_ub;     0.95*t_max;     hull_z_max - 1e-3];
            x0 = [vs_opt;    t_init;         z_fill_seed];

            % Clamp x0 into bounds (defensive)
            x0 = max(lb, min(ub, x0));

            fminopts = optimoptions('fmincon', ...
                'Algorithm',              'sqp', ...
                'Display',                'iter', ...
                'StepTolerance',          1e-8, ...
                'OptimalityTolerance',    1e-6, ...
                'ConstraintTolerance',    1e-6, ...
                'MaxIterations',          200, ...
                'MaxFunctionEvaluations', 1000);

            fprintf('      Solver: fmincon SQP  (DVs=3 ; ceq=buoyancy ; c=GM≥gm_min)\n');
            [x_star, phi_star, exitflag, output] = fmincon( ...
                @obj_fn, x0, [], [], [], [], lb, ub, @con_fn, fminopts);

            vs_star = x_star(1);
            t_star  = x_star(2);
            zf_star = x_star(3);

            t_min_active = ((t_star - t_min) / max(t_max - t_min, eps) < 0.01);

            %% §7  EVALUATE AT THE OPTIMUM
            [grids_star, n_zero_interior] = WEC_Shell_Offset.build_geometry_grid( ...
                config, t_star, n_z_grid, max_slope_factor);
            if n_zero_interior > 0.05 * n_z_grid
                error('WEC_Shell_Offset:TooManyDegenerateSections', ...
                      'Cross-section evaluation failed at %d/%d interior z-grid points (>5%%).', ...
                      n_zero_interior, n_z_grid);
            end
            realised = WEC_Shell_Offset.evaluate_at( ...
                vs_star, t_star, zf_star, grids_star, ctx);
            if ~realised.feasible
                warning('WEC_Shell_Offset:InfeasibleOptimum', ...
                        'fmincon optimum is infeasible — realised hydrostatics produced NaN/Inf.');
            end

            %% §8  PACKAGE OUTPUT  (shape preserved for downstream consumers)
            steel_data = struct();
            steel_data.t_steel          = t_star;
            steel_data.z_fill           = zf_star;
            steel_data.draft            = realised.draft;
            steel_data.vertical_shift   = realised.vertical_shift;
            steel_data.draft_optimiser  = draft_opt;
            steel_data.vs_optimiser     = vs_opt;
            steel_data.rho_steel        = rho_steel;
            steel_data.rho_air          = rho_air;
            steel_data.t_max            = t_max;
            steel_data.t_min            = t_min;
            steel_data.t_min_active     = t_min_active;

            steel_data.V_steel          = realised.V_steel;
            steel_data.V_air            = realised.V_air;
            steel_data.V_hull           = realised.V_steel + realised.V_air;
            steel_data.M_steel          = realised.M_steel;
            steel_data.M_air            = realised.M_air;
            steel_data.M_total          = realised.M_total;

            steel_data.z_cg_steel       = realised.z_cg_steel;
            steel_data.z_cg_air         = realised.z_cg_air;
            steel_data.CG_z_body        = realised.CG_z_body;
            steel_data.CG_z_world       = realised.CG_z_world;

            steel_data.Iyy_total_origin = realised.Iyy_total_origin;
            steel_data.Iyy_about_cg     = realised.Iyy_about_cg;
            steel_data.Ixx_total_origin = realised.Ixx_total_origin;
            steel_data.Ixx_about_cg     = realised.Ixx_about_cg;
            steel_data.Izz_total_origin = realised.Izz_total_origin;
            steel_data.Izz_about_cg     = realised.Izz_about_cg;

            steel_data.V_sub            = realised.V_sub;
            steel_data.Aw               = realised.Aw;
            steel_data.I_wp_yy          = realised.I_wp_yy;
            steel_data.CB_z_world       = realised.CB_z_world;
            steel_data.KM_world         = realised.KM_world;

            steel_data.GM_realised      = realised.GM;
            steel_data.T_heave_realised = realised.T_heave;
            steel_data.T_pitch_realised = realised.T_pitch;
            steel_data.K33_hydro        = realised.K33_hydro;
            steel_data.K55_hydro        = realised.K55_hydro;

            steel_data.A11              = realised.A11;
            steel_data.A33              = realised.A33;
            steel_data.A55              = realised.A55;

            % Targets struct (kept for downstream reporting + plot_steel_solve)
            tgt = struct( ...
                'GM',      final_props.GM_L, ...
                'T_heave', final_props.periods.heave, ...
                'T_pitch', final_props.periods.pitch, ...
                'mass',    final_props.mass_total);
            steel_data.targets = tgt;
            steel_data.residuals.dGM_pct      = 100 * (realised.GM      - tgt.GM)      / tgt.GM;
            steel_data.residuals.dT_heave_pct = 100 * (realised.T_heave - tgt.T_heave) / tgt.T_heave;
            steel_data.residuals.dT_pitch_pct = 100 * (realised.T_pitch - tgt.T_pitch) / tgt.T_pitch;
            steel_data.residuals.dmass_pct    = 100 * (realised.M_total - tgt.mass)    / tgt.mass;
            steel_data.mass_balance_error_pct = 100 * realised.mass_balance_error_abs / max(realised.M_total, eps);

            steel_data.phi_star = phi_star;
            steel_data.feasible = realised.feasible;
            steel_data.exitflag = exitflag;
            steel_data.solver   = 'fmincon-sqp';

            steel_data.z_grid           = grids_star.z;
            steel_data.A_outer_grid     = grids_star.A_outer;
            steel_data.A_inner_grid     = grids_star.A_inner;
            steel_data.A_jacket_grid    = grids_star.A_outer - grids_star.A_inner;
            steel_data.Iyy_outer_grid   = grids_star.Iyy_outer;
            steel_data.Iyy_inner_grid   = grids_star.Iyy_inner;

            %% §8b  PER-STRIP REALISED EQUIVALENT DENSITY
            %  Integrate the realised material distribution over each strip
            %  defined in config.strip_edges so the cake-layer visualisations
            %  (visualize_3d_cross_section, visualize_2d_equivalent,
            %  visualize_3D_equivalent) can colour the hull by the AS-BUILT
            %  effective density instead of the optimiser's continuous rho.
            %    Below z_fill: solid steel (entire cross-section A_outer)
            %    Above z_fill: steel jacket annulus + air interior
            if isfield(config, 'strip_edges') && ~isempty(config.strip_edges)
                strip_edges_v = config.strip_edges(:);
                N_strips_v    = length(strip_edges_v) - 1;
                rho_eff       = zeros(N_strips_v, 1);
                z_grid_v      = steel_data.z_grid;
                Ao_v          = steel_data.A_outer_grid;
                Ai_v          = steel_data.A_inner_grid;
                Aj_v          = steel_data.A_jacket_grid;
                for kk = 1:N_strips_v
                    zlo = strip_edges_v(kk);
                    zhi = strip_edges_v(kk+1);
                    % Breakpoints: strip edges + z_fill (if inside) + original grid samples
                    bp = unique(sort([zlo; zhi; zf_star; ...
                                      z_grid_v(z_grid_v > zlo & z_grid_v < zhi)]));
                    Ao_b = interp1(z_grid_v, Ao_v, bp, 'linear', 0);
                    Ai_b = interp1(z_grid_v, Ai_v, bp, 'linear', 0);
                    Aj_b = max(0, Ao_b - Ai_b);
                    V_outer = trapz(bp, Ao_b);
                    if V_outer <= 1e-12
                        rho_eff(kk) = NaN;
                        continue;
                    end
                    below = bp <= zf_star;
                    above = bp >= zf_star;
                    V_steel = 0;  V_air = 0;
                    if sum(below) >= 2
                        V_steel = V_steel + trapz(bp(below), Ao_b(below));
                    end
                    if sum(above) >= 2
                        V_steel = V_steel + trapz(bp(above), Aj_b(above));
                        V_air   = V_air   + trapz(bp(above), Ai_b(above));
                    end
                    M_strip = rho_steel * V_steel + rho_air * V_air;
                    rho_eff(kk) = M_strip / V_outer;
                end
                steel_data.strip_rho_eff = rho_eff;
                steel_data.strip_edges   = strip_edges_v;
            else
                steel_data.strip_rho_eff = [];
                steel_data.strip_edges   = [];
            end
            steel_data.realisation_mode = 'steel_fill';

            steel_data.elapsed_seconds  = toc(t_solve_start);

            %% §9  CONSOLE REPORT
            fprintf('      ──────────────── Steel-Solve Result (fmincon SQP) ────────────────\n');
            if t_min_active
                fprintf('      t_steel*       : %.5f m  (%.2f in)   *** at t_min bound ***\n', ...
                        t_star, t_star/0.0254);
            else
                fprintf('      t_steel*       : %.5f m  (%.2f in)\n', ...
                        t_star, t_star/0.0254);
            end
            fprintf('      z_fill*        : %.4f m  (hull range [%.3f, %.3f])\n', ...
                    zf_star, hull_z_min, hull_z_max);
            fprintf('      vs*            : %.4f m  (optimiser was %.4f m)\n', ...
                    vs_star, vs_opt);
            fprintf('      M_total        : %.1f kg   target  : %.1f kg   (%+.3f%%)\n', ...
                    realised.M_total, tgt.mass, steel_data.residuals.dmass_pct);
            fprintf('      Mass-balance   : %.3e kg  (%.4f%% of M_total)  ← ceq residual\n', ...
                    realised.mass_balance_error_abs, steel_data.mass_balance_error_pct);
            fprintf('      GM realised    : %.4f m   gm_min=%.3f m   gm_target=%.3f m\n', ...
                    realised.GM, config.gm_min, config.gm_target);
            fprintf('      T_heave        : %.3f s   target %.3f s   (%+.2f%%)\n', ...
                    realised.T_heave, tgt.T_heave, steel_data.residuals.dT_heave_pct);
            fprintf('      T_pitch        : %.3f s   target %.3f s   (%+.2f%%)\n', ...
                    realised.T_pitch, tgt.T_pitch, steel_data.residuals.dT_pitch_pct);
            fprintf('      Phi*           : %.6g  feasible=%d  exitflag=%d  iters=%d  elapsed=%.2f s\n', ...
                    phi_star, realised.feasible, exitflag, output.iterations, ...
                    steel_data.elapsed_seconds);
            fprintf('      ────────────────────────────────────────────────────────────────\n');

            %% --- Nested functions (close over ctx) -----------------------

            function f = obj_fn(x)
                % Objective: phi(r_heave) + phi(r_pitch)
                vs_ = x(1);  t_ = x(2);  zf_ = x(3);
                [grids_, n_zero_] = WEC_Shell_Offset.build_geometry_grid( ...
                    ctx.config, t_, ctx.n_z_grid, ctx.max_slope_factor);
                if n_zero_ > 0.05 * ctx.n_z_grid
                    f = 1e4; return;
                end
                r_ = WEC_Shell_Offset.evaluate_at(vs_, t_, zf_, grids_, ctx);
                if ~isfinite(r_.T_heave) || ~isfinite(r_.T_pitch)
                    f = 1e4; return;
                end
                cfg_ = ctx.config;
                heave_half = 0.5 * (cfg_.T_heave_range(2) - cfg_.T_heave_range(1));
                pitch_half = 0.5 * (cfg_.T_pitch_range(2) - cfg_.T_pitch_range(1));
                r_h = (r_.T_heave - cfg_.T_heave_goal) / max(heave_half, 1e-6);
                r_p = (r_.T_pitch - cfg_.T_pitch_goal) / max(pitch_half, 1e-6);
                f = WEC_Shell_Offset.range_penalty(r_h, cfg_.zone_k_amp) + ...
                    WEC_Shell_Offset.range_penalty(r_p, cfg_.zone_k_amp);
            end

            function [c, ceq] = con_fn(x)
                % c   = 1 − GM/gm_min ≤ 0   (stability floor)
                % ceq = M_total/(ρ_w·V_sub(-vs)) − 1 = 0   (buoyancy balance)
                vs_ = x(1);  t_ = x(2);  zf_ = x(3);
                [grids_, n_zero_] = WEC_Shell_Offset.build_geometry_grid( ...
                    ctx.config, t_, ctx.n_z_grid, ctx.max_slope_factor);
                if n_zero_ > 0.05 * ctx.n_z_grid
                    c = 1.0;  ceq = 1.0;  return;
                end
                r_ = WEC_Shell_Offset.evaluate_at(vs_, t_, zf_, grids_, ctx);
                if ~isfinite(r_.M_total) || ~isfinite(r_.mass_buoyant_force) || ...
                        r_.mass_buoyant_force <= 0
                    c = 1.0;  ceq = 1.0;  return;
                end
                ceq = r_.M_total / r_.mass_buoyant_force - 1.0;
                if isfinite(r_.GM)
                    c = 1.0 - r_.GM / ctx.config.gm_min;
                else
                    c = 1.0;
                end
            end

        end


        %% ═════════════════════════════════════════════════════════════
        %%  RANGE PENALTY  (shared by both realisation solvers)
        %% ═════════════════════════════════════════════════════════════

        function phi = range_penalty(r, k_amp)
        % RANGE_PENALTY  C1 piecewise-quadratic penalty (mirror of
        % WEC_Main_Optimizer's local helper).
        %   r = 0   → phi = 0       (at target)
        %   r = ±1  → phi = 1.0     (at range boundary)
        %   |r|>1   → phi grows k_amp× faster, C1 across the boundary
            if r < -1
                delta = -r - 1;
                phi = 1 + 2*delta + k_amp * delta^2;
            elseif r > 1
                delta = r - 1;
                phi = 1 + 2*delta + k_amp * delta^2;
            else
                phi = r^2;
            end
        end


        %% ═════════════════════════════════════════════════════════════
        %%  EVALUATE AT (vs, t_steel, z_fill)  — no inner fzero
        %% ═════════════════════════════════════════════════════════════

        function out = evaluate_at(vs, t_steel, z_fill, grids, ctx)
        % EVALUATE_AT  Hydrostatics + GM + periods at given (vs, t_steel, z_fill).
        %
        %   Unlike the previous evaluate_realised, this version does NOT
        %   solve an inner fzero for draft — vs is an input.  The mass
        %   balance |M_total − ρ_w·V_sub(-vs)| is exposed as
        %   `out.mass_balance_error_abs` and as `out.mass_buoyant_force`
        %   so the outer fmincon constraint can take it as ceq.
        %
        %   Returned struct field set matches the old evaluate_realised
        %   (downstream consumers unchanged).

            cfg = ctx.config;
            out = WEC_Shell_Offset.empty_realised();

            %% --- Split integrals (steel below z_fill, jacket+air above) ---
            R = WEC_Shell_Offset.integrate_split(grids, z_fill);

            V_steel = R.V_below_outer + R.V_above_jacket;
            V_air   = R.V_above_inner;

            int_z_x_A_steel = R.int_zA_below_outer + R.int_zA_above_jacket;
            int_z_x_A_air   = R.int_zA_above_inner;

            int_x2_steel = R.int_Ix_below_outer + R.int_Ix_above_jacket;
            int_x2_air   = R.int_Ix_above_inner;

            int_y2_steel = R.int_Iy_below_outer + R.int_Iy_above_jacket;
            int_y2_air   = R.int_Iy_above_inner;

            int_z2_A_steel = R.int_z2A_below_outer + R.int_z2A_above_jacket;
            int_z2_A_air   = R.int_z2A_above_inner;

            M_steel = ctx.rho_steel * V_steel;
            M_air   = ctx.rho_air   * V_air;
            M_total = M_steel + M_air;

            out.V_steel = V_steel;
            out.V_air   = V_air;
            out.M_steel = M_steel;
            out.M_air   = M_air;
            out.M_total = M_total;

            if M_total <= 0
                out.feasible = false;
                return;
            end

            % z-centroids (body frame)
            if V_steel > 1e-12
                z_cg_steel = int_z_x_A_steel / V_steel;
            else
                z_cg_steel = 0.5 * (ctx.hull_z_min + z_fill);
            end
            if V_air > 1e-12
                z_cg_air = int_z_x_A_air / V_air;
            else
                z_cg_air = 0.5 * (z_fill + ctx.hull_z_max);
            end
            CG_z_body = (M_steel * z_cg_steel + M_air * z_cg_air) / M_total;

            % Iyy (pitch axis, body origin):  ∫∫∫ ρ (x² + z²) dV
            Iyy_total_origin = ctx.rho_steel * (int_x2_steel + int_z2_A_steel) + ...
                               ctx.rho_air   * (int_x2_air   + int_z2_A_air);
            % Ixx (roll axis):                  ∫∫∫ ρ (y² + z²) dV
            Ixx_total_origin = ctx.rho_steel * (int_y2_steel + int_z2_A_steel) + ...
                               ctx.rho_air   * (int_y2_air   + int_z2_A_air);
            % Izz (yaw):                        ∫∫∫ ρ (x² + y²) dV
            Izz_total_origin = ctx.rho_steel * (int_x2_steel + int_y2_steel) + ...
                               ctx.rho_air   * (int_x2_air   + int_y2_air);

            Iyy_about_cg = max(0, Iyy_total_origin - M_total * CG_z_body^2);
            Ixx_about_cg = max(0, Ixx_total_origin - M_total * CG_z_body^2);
            Izz_about_cg = max(0, Izz_total_origin);

            out.z_cg_steel       = z_cg_steel;
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
            if GM > 0
                K55_hydro = M_total * cfg.G * GM;
            else
                K55_hydro = 0;
            end

            try
                % Pass the CANDIDATE's CG so A55 (and the surge-pitch
                % cross terms) reference the right axis during fmincon
                % iterations — previously this was at HAMS-CG and only
                % corrected post-solve in build_realised_props.
                [A11, A33, A55, ~, ~, ~] = ...
                    WEC_Core_Functions.interpolate_wamit_added_mass( ...
                        vs, cfg, out.CG_z_world);
            catch ME
                error('WEC_Shell_Offset:WAMITInterpolationFailed', ...
                      ['interpolate_wamit_added_mass failed at vs=%.4f m: %s. ' ...
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


        %% ═════════════════════════════════════════════════════════════
        %%  STRIP-AWARE GEOMETRY GRID  (wall + per-strip-thickness)
        %% ═════════════════════════════════════════════════════════════

        function [grids, n_zero_interior] = build_geometry_grid_strip_aware( ...
                config, strip_edges, t_offset_strip, is_solid_strip, ...
                n_z, max_slope_factor)
        % BUILD_GEOMETRY_GRID_STRIP_AWARE  Per-strip thickness + wall/solid mask.
        %
        %   Like build_geometry_grid, but A_inner / Iyy_inner / Ixx_inner are
        %   computed PER STRIP using a different t_offset for each strip, AND
        %   strips flagged is_solid_strip(i) are forced to A_inner = 0
        %   (entire cross-section is solid).
        %
        %   This is the unified geometry model used by the constructable
        %   (UHPC) path:
        %     - Wall strip (strip_is_wall=true) → solid (A_inner=0)
        %     - Below z_fill samples → handled by integrate_split (still
        %       treated as solid because A_outer flows through)
        %     - Phase-1b thickened strip → larger t_offset(i), smaller A_inner
        %     - Phase-1b fully-solid strip → A_inner=0
        %
        %   strip_edges  : (N+1)x1 body-frame z-boundaries
        %   t_offset_strip : Nx1 perpendicular jacket thickness per strip [m]
        %   is_solid_strip : Nx1 logical, true = strip has no inner void

            persistent cache_key cache_grids

            ms2_model = config.ms2_model;
            if ~isprop(ms2_model, 'filename') && ~isfield(ms2_model, 'filename')
                error('WEC_Shell_Offset:NoParserFilename', ...
                      'ms2_model has no .filename property — cannot key cache.');
            end
            parser_key = ms2_model.filename;
            if isempty(parser_key)
                error('WEC_Shell_Offset:EmptyParserFilename', ...
                      'ms2_model.filename is empty.  Cannot key cache.');
            end

            tab_z = config.Aw_table_z(:);
            this_key = sprintf( ...
                '%s|sa|nstrip=%d|t=%s|sld=%s|n=%d|s=%.6g|zN=%d|hlo=%.6g|hhi=%.6g', ...
                parser_key, length(strip_edges)-1, ...
                num2str(t_offset_strip(:)', '%.6g,'), ...
                num2str(double(is_solid_strip(:)'), '%d,'), ...
                n_z, max_slope_factor, length(tab_z), ...
                config.hull_z_min, config.hull_z_max);
            if ~isempty(cache_key) && strcmp(cache_key, this_key)
                grids = cache_grids;
                n_zero_interior = sum(grids.A_outer(2:end-1) <= 0);
                return;
            end

            if ~isfield(config, 'boundary_cache') || ...
                    ~isstruct(config.boundary_cache) || ...
                    ~isfield(config.boundary_cache, 'sources') || ...
                    ~isfield(config.boundary_cache, 'u_samples')
                error('WEC_Shell_Offset:BadBoundaryCache', ...
                      'config.boundary_cache is missing required fields.');
            end
            n_u = length(config.boundary_cache.u_samples);

            z = linspace(config.hull_z_min, config.hull_z_max, n_z)';
            A_outer   = max(0, interp1(config.Aw_table_z, config.Aw_table,      z, 'linear', 0));
            Iyy_outer = max(0, interp1(config.Aw_table_z, config.I_wp_yy_table, z, 'linear', 0));
            if isfield(config, 'I_wp_xx_table') && ~isempty(config.I_wp_xx_table)
                Ixx_outer = max(0, interp1(config.Aw_table_z, config.I_wp_xx_table, z, 'linear', 0));
            else
                Ixx_outer = Iyy_outer;
            end

            A_inner   = zeros(n_z, 1);
            Iyy_inner = zeros(n_z, 1);
            Ixx_inner = zeros(n_z, 1);

            % Map each z-sample to the strip that contains it
            N_strips = length(t_offset_strip);
            strip_of_z = zeros(n_z, 1);
            for k = 1:n_z
                idx = find(z(k) >= strip_edges(1:end-1) - 1e-9 & ...
                           z(k) <= strip_edges(2:end)   + 1e-9, 1, 'first');
                if isempty(idx)
                    if z(k) < strip_edges(1)
                        idx = 1;
                    else
                        idx = N_strips;
                    end
                end
                strip_of_z(k) = idx;
            end

            n_zero_interior = 0;
            for k = 1:n_z
                if A_outer(k) <= 1e-10
                    if k > 1 && k < n_z
                        n_zero_interior = n_zero_interior + 1;
                    end
                    continue;
                end
                i_strip = strip_of_z(k);
                if is_solid_strip(i_strip)
                    % Solid strip → no inner void
                    A_inner(k) = 0;  Iyy_inner(k) = 0;  Ixx_inner(k) = 0;
                    continue;
                end
                t_k = t_offset_strip(i_strip);
                if ~isfinite(t_k) || t_k <= 0
                    A_inner(k) = 0;  Iyy_inner(k) = 0;  Ixx_inner(k) = 0;
                    continue;
                end
                [A_inner(k), Iyy_inner(k), Ixx_inner(k)] = ...
                    WEC_Shell_Offset.inner_props_at_z( ...
                        config, z(k), t_k, max_slope_factor, ...
                        n_u, A_outer(k), Iyy_outer(k), Ixx_outer(k));
            end

            A_inner   = min(A_inner,   A_outer);
            A_inner   = max(A_inner,   0);
            Iyy_inner = min(Iyy_inner, Iyy_outer);
            Iyy_inner = max(Iyy_inner, 0);
            Ixx_inner = min(Ixx_inner, Ixx_outer);
            Ixx_inner = max(Ixx_inner, 0);

            grids = struct('z', z, ...
                           'A_outer',   A_outer, ...
                           'A_inner',   A_inner, ...
                           'Iyy_outer', Iyy_outer, ...
                           'Iyy_inner', Iyy_inner, ...
                           'Ixx_outer', Ixx_outer, ...
                           'Ixx_inner', Ixx_inner, ...
                           'strip_of_z', strip_of_z);

            cache_key   = this_key;
            cache_grids = grids;
        end


        %% ═════════════════════════════════════════════════════════════
        %%  EVALUATE REALISED — STRIP-AWARE  (per-strip t_offset + wall)
        %% ═════════════════════════════════════════════════════════════

        function out = evaluate_at_strip_aware(vs, t_offset_strip, is_solid_strip, ...
                                                z_fill, grids, ctx)
        % EVALUATE_AT_STRIP_AWARE  Same as evaluate_at, but the geometry grid
        %   was built with per-strip thickness via build_geometry_grid_strip_aware.
        %   t_offset_strip / is_solid_strip are not used directly in this body —
        %   the per-strip geometry is already baked into `grids`.  They are
        %   kept in the signature so the call site is self-documenting and
        %   parallel to the steel evaluator.
        %
        %   Returned struct has the same shape as evaluate_at.

            cfg = ctx.config;
            out = WEC_Shell_Offset.empty_realised();

            R = WEC_Shell_Offset.integrate_split(grids, z_fill);

            V_steel = R.V_below_outer + R.V_above_jacket;
            V_air   = R.V_above_inner;

            int_z_x_A_steel = R.int_zA_below_outer + R.int_zA_above_jacket;
            int_z_x_A_air   = R.int_zA_above_inner;
            int_x2_steel    = R.int_Ix_below_outer + R.int_Ix_above_jacket;
            int_x2_air      = R.int_Ix_above_inner;
            int_y2_steel    = R.int_Iy_below_outer + R.int_Iy_above_jacket;
            int_y2_air      = R.int_Iy_above_inner;
            int_z2_A_steel  = R.int_z2A_below_outer + R.int_z2A_above_jacket;
            int_z2_A_air    = R.int_z2A_above_inner;

            M_steel = ctx.rho_steel * V_steel;
            M_air   = ctx.rho_air   * V_air;
            M_total = M_steel + M_air;

            out.V_steel = V_steel;
            out.V_air   = V_air;
            out.M_steel = M_steel;
            out.M_air   = M_air;
            out.M_total = M_total;

            if M_total <= 0
                out.feasible = false;
                return;
            end

            if V_steel > 1e-12
                z_cg_steel = int_z_x_A_steel / V_steel;
            else
                z_cg_steel = 0.5 * (ctx.hull_z_min + z_fill);
            end
            if V_air > 1e-12
                z_cg_air = int_z_x_A_air / V_air;
            else
                z_cg_air = 0.5 * (z_fill + ctx.hull_z_max);
            end
            CG_z_body = (M_steel * z_cg_steel + M_air * z_cg_air) / M_total;

            Iyy_total_origin = ctx.rho_steel * (int_x2_steel + int_z2_A_steel) + ...
                               ctx.rho_air   * (int_x2_air   + int_z2_A_air);
            Ixx_total_origin = ctx.rho_steel * (int_y2_steel + int_z2_A_steel) + ...
                               ctx.rho_air   * (int_y2_air   + int_z2_A_air);
            Izz_total_origin = ctx.rho_steel * (int_x2_steel + int_y2_steel) + ...
                               ctx.rho_air   * (int_x2_air   + int_y2_air);
            Iyy_about_cg = max(0, Iyy_total_origin - M_total * CG_z_body^2);
            Ixx_about_cg = max(0, Ixx_total_origin - M_total * CG_z_body^2);
            Izz_about_cg = max(0, Izz_total_origin);

            out.z_cg_steel       = z_cg_steel;
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
            if GM > 0
                K55_hydro = M_total * cfg.G * GM;
            else
                K55_hydro = 0;
            end

            try
                % Pass the candidate's CG — same reasoning as evaluate_at.
                [A11, A33, A55, ~, ~, ~] = ...
                    WEC_Core_Functions.interpolate_wamit_added_mass( ...
                        vs, cfg, out.CG_z_world);
            catch ME
                error('WEC_Shell_Offset:WAMITInterpolationFailed', ...
                      'interpolate_wamit_added_mass failed at vs=%.4f m: %s.', ...
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


        %% ═════════════════════════════════════════════════════════════
        %%  CONSTRUCTABLE SOLVE — wall-aware Phase 1 + per-strip Phase 1b
        %% ═════════════════════════════════════════════════════════════

        function solve_data = solve_constructable(config, x_opt_3d, final_props, opts)
        % SOLVE_CONSTRUCTABLE  UHPC realisation via fmincon SQP.
        %
        %   Single-pass constrained optimisation that replaces the prior
        %   two-phase (global + greedy strip walk) implementation.
        %
        %   FORMULATION  (mirrors WEC_Shell_Offset.solve and Stage 2)
        %     DVs         x = [vs; z_fill; t_offset_strip(non-wall)]
        %                       N_dv = 2 + (N_strips − 1)
        %     Equality    ceq = M_total / (ρ_w · V_sub(-vs)) − 1 = 0
        %                       (buoyancy balance)
        %     Inequality  c   = 1 − GM/gm_min ≤ 0
        %                       (stability floor)
        %     Objective   Φ   = phi(r_heave) + phi(r_pitch)
        %                       range-normalised, same form as steel solve
        %     Solver      fmincon SQP
        %     Warm start  vs = x_opt_3d(1), z_fill = analytic seed,
        %                 t_offset uniform = t_init
        %
        %   The wall strip is pinned solid throughout.  Below-z_fill strips
        %   are integrated as solid via integrate_split (no explicit
        %   is_solid flip needed during the solve).  After convergence,
        %   any strip whose upper edge sits at or below z_fill* is flagged
        %   is_solid for downstream consumers (visualisation,
        %   extract_strip_geometry).
        %
        %   GM and mass are NOT in the objective — they are constraints.
        %   No per-target weights — range normalisation balances the
        %   T_heave and T_pitch residuals automatically.
        %
        %   solve_data is shape-compatible with WEC_Shell_Offset.solve's
        %   steel_data, with the per-strip extras (t_offset_strip,
        %   is_solid_strip, strip_edges, wall_strip_idx) appended.

            t_solve_start = tic;
            fprintf('\n    WEC_Shell_Offset.solve_constructable (fmincon SQP, constrained):\n');

            if nargin < 4 || isempty(opts), opts = struct(); end

            % ── Required wall + strip bookkeeping from caller ──────────
            if ~isfield(opts, 'wall_strip_idx') || isempty(opts.wall_strip_idx)
                error('WEC_Shell_Offset:NoWallStripIdx', ...
                      'opts.wall_strip_idx must be supplied to solve_constructable.');
            end
            if ~isfield(opts, 'strip_edges') || isempty(opts.strip_edges)
                error('WEC_Shell_Offset:NoStripEdges', ...
                      'opts.strip_edges must be supplied to solve_constructable.');
            end
            wall_strip_idx = opts.wall_strip_idx;
            strip_edges    = opts.strip_edges(:);
            N_strips       = length(strip_edges) - 1;
            nw_idx         = setdiff(1:N_strips, wall_strip_idx);
            N_nw           = length(nw_idx);

            % Required config fields (mirror steel solve)
            required_cfg = {'ms2_model','hull_z_min','hull_z_max', ...
                            'Aw_table_z','Aw_table','V_sub_table','CB_z_table', ...
                            'I_wp_yy_table','RHO_WATER','G', ...
                            'rho_steel','rho_air','steel_t_init','steel_t_min', ...
                            'steel_max_slope_factor','steel_n_z_grid', ...
                            'vertical_shift_bounds','gm_min', ...
                            'T_heave_goal','T_heave_range', ...
                            'T_pitch_goal','T_pitch_range','zone_k_amp'};
            for k = 1:length(required_cfg)
                if ~isfield(config, required_cfg{k}) || isempty(config.(required_cfg{k}))
                    error('WEC_Shell_Offset:MissingConfig', ...
                          'Required config field "%s" is missing or empty.', ...
                          required_cfg{k});
                end
            end

            % Resolve material + numerical options
            rho_uhpc         = WEC_Shell_Offset.opt_or_cfg(opts, 'rho_steel',        config.rho_steel);
            rho_air          = WEC_Shell_Offset.opt_or_cfg(opts, 'rho_air',          config.rho_air);
            t_init           = WEC_Shell_Offset.opt_or_cfg(opts, 't_init',           config.steel_t_init);
            t_min            = WEC_Shell_Offset.opt_or_cfg(opts, 't_min',            config.steel_t_min);
            max_slope_factor = WEC_Shell_Offset.opt_or_cfg(opts, 'max_slope_factor', config.steel_max_slope_factor);
            n_z_grid         = WEC_Shell_Offset.opt_or_cfg(opts, 'n_z_grid',         config.steel_n_z_grid);

            assert(rho_uhpc > 0, 'rho_uhpc must be > 0 (got %g)', rho_uhpc);
            assert(rho_air > 0 && rho_air < rho_uhpc, ...
                'rho_air (%g) must be in (0, rho_uhpc=%g)', rho_air, rho_uhpc);

            hull_z_min = config.hull_z_min;
            hull_z_max = config.hull_z_max;
            x_opt_3d   = x_opt_3d(:);
            vs_opt     = x_opt_3d(1);
            draft_opt  = abs(hull_z_min + vs_opt);

            % ── Probe maximum permissible thickness (same as solve) ────
            cache = config.boundary_cache;
            n_probe = 9;
            z_probes = linspace(0.1*hull_z_min + 0.9*hull_z_max, ...
                                0.9*hull_z_min + 0.1*hull_z_max, n_probe);
            r_probes = zeros(n_probe, 1);
            for k = 1:n_probe
                r_probes(k) = WEC_HydroProperties.compute_rmin_at_z( ...
                                  config.ms2_model, z_probes(k), 60, cache);
            end
            r_min_global = min(r_probes(r_probes > 0));
            if isempty(r_min_global) || r_min_global <= 0
                error('WEC_Shell_Offset:NoFiniteRadius', ...
                      'compute_rmin_at_z returned no positive radius.');
            end
            t_max = 0.5 * r_min_global;
            if t_min >= t_max
                error('WEC_Shell_Offset:TMinExceedsTMax', ...
                      't_min=%.5f >= t_max=%.5f m.', t_min, t_max);
            end
            t_init = max(min(t_init, 0.9*t_max), max(t_min, 1e-4));

            vs_lb = config.vertical_shift_bounds(1);
            vs_ub = config.vertical_shift_bounds(2);

            fprintf('      Hull z=[%.3f, %.3f] m  N_strips=%d  wall=%d  N_nw=%d\n', ...
                    hull_z_min, hull_z_max, N_strips, wall_strip_idx, N_nw);
            fprintf('      Targets (from optimiser): GM_min=%.3f m  T_h=%.3f s  T_p=%.3f s  M=%.1f kg\n', ...
                    config.gm_min, config.T_heave_goal, config.T_pitch_goal, ...
                    final_props.mass_total);
            fprintf('      Bounds : t∈[%.5f, %.4f] m  vs∈[%.3f, %.3f] m  ρ_UHPC=%.0f  ρ_void=%.2f\n', ...
                    t_min, 0.95*t_max, vs_lb, vs_ub, rho_uhpc, rho_air);

            %% Analytic z_fill seed at (vs=vs_opt, t=t_init, wall=solid)
            is_solid_seed = false(N_strips, 1);
            is_solid_seed(wall_strip_idx) = true;
            t_strip_seed = t_init * ones(N_strips, 1);
            t_strip_seed(wall_strip_idx) = inf;

            [grids_seed, n_zero_seed] = WEC_Shell_Offset.build_geometry_grid_strip_aware( ...
                config, strip_edges, t_strip_seed, is_solid_seed, n_z_grid, max_slope_factor);
            if n_zero_seed > 0.05 * n_z_grid
                error('WEC_Shell_Offset:DegenerateSeedGeometry', ...
                      'Seed geometry has %d/%d (>5%%) degenerate sections.', ...
                      n_zero_seed, n_z_grid);
            end

            V_sub_at_vs_opt = max(0, interp1(config.Aw_table_z, ...
                                              config.V_sub_table, ...
                                              -vs_opt, 'linear', 0));
            M_buoy_target = config.RHO_WATER * V_sub_at_vs_opt;
            V_jacket = trapz(grids_seed.z, grids_seed.A_outer - grids_seed.A_inner);
            V_inner  = trapz(grids_seed.z, grids_seed.A_inner);
            V_inner_below_target = (M_buoy_target - rho_uhpc*V_jacket - rho_air*V_inner) / ...
                                   (rho_uhpc - rho_air);

            if V_inner_below_target < 0
                z_fill_seed = hull_z_min + 1e-3;
            elseif V_inner_below_target > V_inner
                z_fill_seed = hull_z_max - 1e-3;
            else
                cumV = cumtrapz(grids_seed.z, grids_seed.A_inner);
                interior = grids_seed.A_inner > 1e-10;
                if sum(interior) < 2
                    z_fill_seed = 0.5 * (hull_z_min + hull_z_max);
                else
                    [V_uniq, ia_uniq] = unique(cumV(interior), 'stable');
                    z_int = grids_seed.z(interior);
                    z_fill_seed = interp1(V_uniq, z_int(ia_uniq), ...
                                          V_inner_below_target, 'linear', 'extrap');
                end
                z_fill_seed = max(hull_z_min + 1e-3, ...
                                  min(hull_z_max - 1e-3, z_fill_seed));
            end
            fprintf('      Warm-start z_fill seed: %.4f m (analytic, M_buoy=%.0f kg)\n', ...
                    z_fill_seed, M_buoy_target);

            %% ctx for nested closures
            ctx = struct();
            ctx.config           = config;
            ctx.rho_steel        = rho_uhpc;
            ctx.rho_air          = rho_air;
            ctx.max_slope_factor = max_slope_factor;
            ctx.n_z_grid         = n_z_grid;
            ctx.hull_z_min       = hull_z_min;
            ctx.hull_z_max       = hull_z_max;
            ctx.strip_edges      = strip_edges;
            ctx.N_strips         = N_strips;
            ctx.wall_strip_idx   = wall_strip_idx;
            ctx.nw_idx           = nw_idx;
            % is_solid_strip held fixed during the solve: wall pinned solid;
            % all others non-solid (z_fill controls the below-z_fill solid
            % region implicitly via integrate_split).  Promotion to
            % is_solid_strip(i)=true for below-z_fill strips happens AFTER
            % the solve, for visualisation/extract_strip_geometry consumers.
            ctx.is_solid_strip = is_solid_seed;

            %% fmincon SQP
            % DVs: x = [vs; z_fill; t_offset_strip(nw_idx)]
            lb = [vs_lb;     hull_z_min + 1e-3;     t_min        * ones(N_nw, 1)];
            ub = [vs_ub;     hull_z_max - 1e-3;     0.95 * t_max * ones(N_nw, 1)];
            x0 = [vs_opt;    z_fill_seed;           t_init       * ones(N_nw, 1)];
            x0 = max(lb, min(ub, x0));

            fminopts = optimoptions('fmincon', ...
                'Algorithm',              'sqp', ...
                'Display',                'iter', ...
                'StepTolerance',          1e-8, ...
                'OptimalityTolerance',    1e-6, ...
                'ConstraintTolerance',    1e-6, ...
                'MaxIterations',          300, ...
                'MaxFunctionEvaluations', 2000);

            fprintf('      Solver: fmincon SQP  (DVs=%d ; ceq=buoyancy ; c=GM≥gm_min)\n', 2 + N_nw);
            [x_star, phi_star, exitflag, output] = fmincon( ...
                @obj_fn, x0, [], [], [], [], lb, ub, @con_fn, fminopts);

            vs_star    = x_star(1);
            zf_star    = x_star(2);
            t_nw_star  = x_star(3:end);

            %% Rebuild per-strip arrays from converged x*
            t_offset_strip = zeros(N_strips, 1);
            is_solid_final = ctx.is_solid_strip;     % wall already solid
            t_offset_strip(wall_strip_idx) = inf;
            for kk = 1:N_nw
                t_offset_strip(nw_idx(kk)) = t_nw_star(kk);
            end
            % Promote strips fully below z_fill* to solid (consumer parity)
            for i = 1:N_strips
                if i == wall_strip_idx, continue; end
                if strip_edges(i+1) <= zf_star + 1e-9
                    is_solid_final(i) = true;
                    t_offset_strip(i) = inf;
                end
            end

            t_min_active = any(abs(t_nw_star - t_min) < 1e-5);

            %% Final evaluation at the converged x*
            [grids_final, n_zero_interior] = WEC_Shell_Offset.build_geometry_grid_strip_aware( ...
                config, strip_edges, t_offset_strip, is_solid_final, n_z_grid, max_slope_factor);
            if n_zero_interior > 0.05 * n_z_grid
                error('WEC_Shell_Offset:TooManyDegenerateSections', ...
                      'Final strip-aware grid has %d/%d (>5%%) degenerate sections.', ...
                      n_zero_interior, n_z_grid);
            end
            realised = WEC_Shell_Offset.evaluate_at_strip_aware( ...
                vs_star, t_offset_strip, is_solid_final, zf_star, grids_final, ctx);

            if ~realised.feasible
                warning('WEC_Shell_Offset:InfeasibleUHPCOptimum', ...
                        'fmincon UHPC optimum is infeasible — hydrostatics produced NaN/Inf.');
            end

            %% Targets struct (for downstream reporting)
            tgt = struct( ...
                'GM',      final_props.GM_L, ...
                'T_heave', final_props.periods.heave, ...
                'T_pitch', final_props.periods.pitch, ...
                'mass',    final_props.mass_total);

            %% Package output (steel_data shape + per-strip extras)
            solve_data = struct();
            solve_data.t_steel          = mean(t_nw_star);   % representative global value
            solve_data.z_fill           = zf_star;
            solve_data.draft            = realised.draft;
            solve_data.vertical_shift   = vs_star;
            solve_data.draft_optimiser  = draft_opt;
            solve_data.vs_optimiser     = vs_opt;
            solve_data.rho_steel        = rho_uhpc;
            solve_data.rho_air          = rho_air;
            solve_data.t_max            = t_max;
            solve_data.t_min            = t_min;
            solve_data.t_min_active     = t_min_active;

            solve_data.V_steel = realised.V_steel;
            solve_data.V_air   = realised.V_air;
            solve_data.V_hull  = realised.V_steel + realised.V_air;
            solve_data.M_steel = realised.M_steel;
            solve_data.M_air   = realised.M_air;
            solve_data.M_total = realised.M_total;

            solve_data.z_cg_steel       = realised.z_cg_steel;
            solve_data.z_cg_air         = realised.z_cg_air;
            solve_data.CG_z_body        = realised.CG_z_body;
            solve_data.CG_z_world       = realised.CG_z_world;

            solve_data.Iyy_total_origin = realised.Iyy_total_origin;
            solve_data.Iyy_about_cg     = realised.Iyy_about_cg;
            solve_data.Ixx_total_origin = realised.Ixx_total_origin;
            solve_data.Ixx_about_cg     = realised.Ixx_about_cg;
            solve_data.Izz_total_origin = realised.Izz_total_origin;
            solve_data.Izz_about_cg     = realised.Izz_about_cg;

            solve_data.V_sub      = realised.V_sub;
            solve_data.Aw         = realised.Aw;
            solve_data.I_wp_yy    = realised.I_wp_yy;
            solve_data.CB_z_world = realised.CB_z_world;
            solve_data.KM_world   = realised.KM_world;

            solve_data.GM_realised      = realised.GM;
            solve_data.T_heave_realised = realised.T_heave;
            solve_data.T_pitch_realised = realised.T_pitch;
            solve_data.K33_hydro        = realised.K33_hydro;
            solve_data.K55_hydro        = realised.K55_hydro;
            solve_data.A11              = realised.A11;
            solve_data.A33              = realised.A33;
            solve_data.A55              = realised.A55;

            solve_data.targets          = tgt;
            solve_data.residuals.dGM_pct      = 100 * (realised.GM      - tgt.GM)      / tgt.GM;
            solve_data.residuals.dT_heave_pct = 100 * (realised.T_heave - tgt.T_heave) / tgt.T_heave;
            solve_data.residuals.dT_pitch_pct = 100 * (realised.T_pitch - tgt.T_pitch) / tgt.T_pitch;
            solve_data.residuals.dmass_pct    = 100 * (realised.M_total - tgt.mass)    / tgt.mass;
            solve_data.mass_balance_error_pct = 100 * realised.mass_balance_error_abs / max(realised.M_total, eps);

            solve_data.phi_star = phi_star;
            solve_data.feasible = realised.feasible;
            solve_data.exitflag = exitflag;
            solve_data.solver   = 'fmincon-sqp';

            solve_data.z_grid           = grids_final.z;
            solve_data.A_outer_grid     = grids_final.A_outer;
            solve_data.A_inner_grid     = grids_final.A_inner;
            solve_data.A_jacket_grid    = grids_final.A_outer - grids_final.A_inner;
            solve_data.Iyy_outer_grid   = grids_final.Iyy_outer;
            solve_data.Iyy_inner_grid   = grids_final.Iyy_inner;

            % Per-strip realisation
            solve_data.t_offset_strip = t_offset_strip;
            solve_data.is_solid_strip = is_solid_final;
            solve_data.strip_edges    = strip_edges;
            solve_data.wall_strip_idx = wall_strip_idx;

            solve_data.elapsed_seconds = toc(t_solve_start);

            %% Console report
            fprintf('      ──────── Constructable Solve Result (fmincon SQP) ────────\n');
            fprintf('      vs*     : %.4f m   z_fill* : %.4f m   draft* : %.4f m\n', ...
                    vs_star, zf_star, realised.draft);
            fprintf('      M_total : %.1f kg   target  %.1f kg   (%+.2f%%)\n', ...
                    realised.M_total, tgt.mass, solve_data.residuals.dmass_pct);
            fprintf('      Mass-balance : %.3e kg  (%.4f%% of M_total)  ← ceq residual\n', ...
                    realised.mass_balance_error_abs, solve_data.mass_balance_error_pct);
            fprintf('      GM      : %.4f m   gm_min=%.3f m   gm_target=%.3f m\n', ...
                    realised.GM, config.gm_min, config.gm_target);
            fprintf('      T_heave : %.3f s   target %.3f s   (%+.2f%%)\n', ...
                    realised.T_heave, tgt.T_heave, solve_data.residuals.dT_heave_pct);
            fprintf('      T_pitch : %.3f s   target %.3f s   (%+.2f%%)\n', ...
                    realised.T_pitch, tgt.T_pitch, solve_data.residuals.dT_pitch_pct);
            fprintf('      Per-strip t_offset (mm): ');
            for ii = 1:N_strips
                if is_solid_final(ii)
                    fprintf('S ');
                else
                    fprintf('%.0f ', t_offset_strip(ii)*1000);
                end
            end
            fprintf('\n');
            fprintf('      Phi*    : %.6g  exitflag=%d  iters=%d  elapsed=%.2f s\n', ...
                    phi_star, exitflag, output.iterations, solve_data.elapsed_seconds);
            fprintf('      ──────────────────────────────────────────────────────────\n');

            %% --- Nested functions ------------------------------------------

            function [t_off, is_sol] = unpack_dvs(x)
                t_off = zeros(ctx.N_strips, 1);
                is_sol = ctx.is_solid_strip;
                t_off(ctx.wall_strip_idx) = inf;
                for kk2 = 1:length(ctx.nw_idx)
                    t_off(ctx.nw_idx(kk2)) = x(2 + kk2);
                end
            end

            function f = obj_fn(x)
                vs_ = x(1);  zf_ = x(2);
                [t_off, is_sol] = unpack_dvs(x);
                [g_, n_zero_] = WEC_Shell_Offset.build_geometry_grid_strip_aware( ...
                    ctx.config, ctx.strip_edges, t_off, is_sol, ...
                    ctx.n_z_grid, ctx.max_slope_factor);
                if n_zero_ > 0.05 * ctx.n_z_grid
                    f = 1e4; return;
                end
                r_ = WEC_Shell_Offset.evaluate_at_strip_aware( ...
                    vs_, t_off, is_sol, zf_, g_, ctx);
                if ~isfinite(r_.T_heave) || ~isfinite(r_.T_pitch)
                    f = 1e4; return;
                end
                cfg_ = ctx.config;
                h_half = 0.5 * (cfg_.T_heave_range(2) - cfg_.T_heave_range(1));
                p_half = 0.5 * (cfg_.T_pitch_range(2) - cfg_.T_pitch_range(1));
                r_h = (r_.T_heave - cfg_.T_heave_goal) / max(h_half, 1e-6);
                r_p = (r_.T_pitch - cfg_.T_pitch_goal) / max(p_half, 1e-6);
                f = WEC_Shell_Offset.range_penalty(r_h, cfg_.zone_k_amp) + ...
                    WEC_Shell_Offset.range_penalty(r_p, cfg_.zone_k_amp);
            end

            function [c, ceq] = con_fn(x)
                vs_ = x(1);  zf_ = x(2);
                [t_off, is_sol] = unpack_dvs(x);
                [g_, n_zero_] = WEC_Shell_Offset.build_geometry_grid_strip_aware( ...
                    ctx.config, ctx.strip_edges, t_off, is_sol, ...
                    ctx.n_z_grid, ctx.max_slope_factor);
                if n_zero_ > 0.05 * ctx.n_z_grid
                    c = 1.0; ceq = 1.0; return;
                end
                r_ = WEC_Shell_Offset.evaluate_at_strip_aware( ...
                    vs_, t_off, is_sol, zf_, g_, ctx);
                if ~isfinite(r_.M_total) || ~isfinite(r_.mass_buoyant_force) || ...
                        r_.mass_buoyant_force <= 0
                    c = 1.0; ceq = 1.0; return;
                end
                ceq = r_.M_total / r_.mass_buoyant_force - 1.0;
                if isfinite(r_.GM)
                    c = 1.0 - r_.GM / ctx.config.gm_min;
                else
                    c = 1.0;
                end
            end

        end


        %% ═════════════════════════════════════════════════════════════
        %%  EMPTY REALISED-PROPS PROTOTYPE  (always-populated struct)
        %% ═════════════════════════════════════════════════════════════

        function out = empty_realised()
        % EMPTY_REALISED  Pre-fills every field that evaluate_realised may
        %   produce with NaN (or false), so the caller can unpack the
        %   struct uniformly even when the candidate is infeasible.
            out = struct( ...
                'feasible',               false, ...
                'draft',                  NaN, ...
                'vertical_shift',         NaN, ...
                'V_steel',                NaN, ...
                'V_air',                  NaN, ...
                'M_steel',                NaN, ...
                'M_air',                  NaN, ...
                'M_total',                NaN, ...
                'mass_buoyant_force',     NaN, ...
                'mass_balance_error_abs', NaN, ...
                'z_cg_steel',             NaN, ...
                'z_cg_air',               NaN, ...
                'CG_z_body',              NaN, ...
                'CG_z_world',             NaN, ...
                'Iyy_total_origin',       NaN, ...
                'Iyy_about_cg',           NaN, ...
                'Ixx_total_origin',       NaN, ...
                'Ixx_about_cg',           NaN, ...
                'Izz_total_origin',       NaN, ...
                'Izz_about_cg',           NaN, ...
                'V_sub',                  NaN, ...
                'Aw',                     NaN, ...
                'I_wp_yy',                NaN, ...
                'CB_z_world',             NaN, ...
                'KM_world',               NaN, ...
                'GM',                     NaN, ...
                'T_heave',                NaN, ...
                'T_pitch',                NaN, ...
                'K33_hydro',              NaN, ...
                'K55_hydro',              NaN, ...
                'A11',                    NaN, ...
                'A33',                    NaN, ...
                'A55',                    NaN);
        end


        %% ═════════════════════════════════════════════════════════════
        %%  BUILD REALISED FINAL_PROPS  (replace optimiser values with
        %%  the as-built steel-fill hull's properties)
        %% ═════════════════════════════════════════════════════════════

        function realised = build_realised_props(final_props, steel_data, config)
        % BUILD_REALISED_PROPS  Drop-in replacement for final_props that
        %   carries the AS-BUILT steel-fill hull's mass / CG / inertia /
        %   GM / hydrostatics / added-mass / natural periods.
        %
        %   realised = WEC_Shell_Offset.build_realised_props(final_props, ...
        %                                                    steel_data, config)
        %
        %   Field-by-field map:
        %     - Mass / CG / V_sub / Aw / GM / KM / CB        ← steel_data (re-solved at the steel draft)
        %     - Inertia_Tensor (Ixx, Iyy, Izz)               ← steel_data parallel-axis to realised CG
        %     - A11 / A33 / A55 / A_full / B_full            ← interpolate WAMIT at the realised draft AT the realised CG via interpolate_wamit_added_mass(..., z_cg_target) — single centralised transform; no inline delta block here
        %     - K_hydro / K_pto / K_total                    ← K_pto is hull-independent (copy); K_hydro recomputed; K_total = K_hydro + K_pto
        %     - periods.heave / pitch                         ← from steel_data
        %     - periods.surge                                 ← recompute via realised mass + A11 + K11_total
        %     - coupled_periods / coupled_modes / participation_factors ← re-eig with realised mass tensor
        %     - MassMatrix_CG / MassMatrix_Origin            ← rebuild via WEC_Core_Functions.calculate6x6MassMatrix
        %     - components / densities_at_nodes / cross_section ← copied from final_props (the OPTIMISER's per-strip rho — kept for reference, not strictly meaningful for steel-fill)
        %
        %   INFEASIBLE-FALLBACK
        %     If steel_data.feasible == false, returns final_props unchanged
        %     (a loud warning is emitted by the caller).

            if ~isstruct(steel_data) || ~isfield(steel_data, 'feasible') || ~steel_data.feasible
                realised = final_props;
                return;
            end

            realised = struct();

            % ── Frame ──────────────────────────────────────────────────
            realised.vertical_shift = steel_data.vertical_shift;
            realised.draft          = steel_data.draft;

            % ── Hydrostatics from steel_data + interp1 from config tables ─
            z_wl_body = -steel_data.vertical_shift;
            z_sub_top = min(z_wl_body, config.hull_z_max);

            realised.Aw      = steel_data.Aw;
            realised.I_wp_yy = steel_data.I_wp_yy;
            if isfield(config, 'I_wp_xx_table') && ~isempty(config.I_wp_xx_table)
                realised.I_wp_xx = max(0, interp1(config.Aw_table_z, ...
                                                   config.I_wp_xx_table, ...
                                                   z_sub_top, 'linear', 0));
            else
                realised.I_wp_xx = realised.I_wp_yy;   % symmetric-hull fallback
            end
            realised.V_sub = steel_data.V_sub;
            realised.CB    = [0, 0, steel_data.CB_z_world];

            % Wetted surface area (handle missing S_wet_table)
            if isfield(config, 'S_wet_table') && ~isempty(config.S_wet_table)
                realised.A_sub = max(0, interp1(config.Aw_table_z, ...
                                                 config.S_wet_table, ...
                                                 z_sub_top, 'linear', 0));
            elseif isfield(final_props, 'A_sub')
                realised.A_sub = final_props.A_sub;
            else
                realised.A_sub = 0;
            end

            realised.mass_buoyant_force = config.RHO_WATER * realised.V_sub;
            realised.KM                 = steel_data.KM_world;

            % ── Mass / CG / Inertia tensor (full 3×3) ─────────────────
            realised.mass_total = steel_data.M_total;
            realised.CG_total   = [0, 0, steel_data.CG_z_world];

            realised.Ixx = steel_data.Ixx_about_cg;
            realised.Iyy = steel_data.Iyy_about_cg;
            realised.Izz = steel_data.Izz_about_cg;
            realised.Inertia_Tensor = diag([realised.Ixx, realised.Iyy, realised.Izz]);

            realised.GM_L              = steel_data.GM_realised;
            realised.mass_discrepancy  = realised.mass_total - realised.mass_buoyant_force;

            % Loud warnings on physically suspect realised values
            if realised.mass_total <= 0
                error('WEC_Shell_Offset:RealisedZeroMass', ...
                      'Realised mass_total = %g (<= 0).', realised.mass_total);
            end
            if realised.GM_L <= 0
                warning('WEC_Shell_Offset:RealisedNegativeGM', ...
                        ['Realised GM_L = %.4f m is non-positive — the steel-fill hull is ', ...
                         'STATICALLY UNSTABLE.  The optimiser-targeted GM may not have been ', ...
                         'achievable under the steel/air mass partition + t_min constraint.'], ...
                        realised.GM_L);
            end
            if realised.Iyy <= 0
                warning('WEC_Shell_Offset:RealisedZeroIyy', ...
                        'Realised Iyy_about_cg = %g — pitch period will be undefined.', ...
                        realised.Iyy);
            end

            % ── Hydrostatic stiffness ─────────────────────────────────
            K33_hydro = config.RHO_WATER * config.G * realised.Aw;
            if realised.GM_L > 0
                K55_hydro = realised.mass_total * config.G * realised.GM_L;
            else
                K55_hydro = 0;
            end
            realised.K_hydro = diag([0, K33_hydro, K55_hydro]);

            % ── Added mass at the realised draft, referenced to the
            %    REALISED CG in a single call.  The delta congruence
            %    transform (formerly inline here) now lives inside
            %    WEC_Core_Functions.interpolate_wamit_added_mass.  See
            %    that function for the math.
            try
                [A11, A33, A55, ~, A_full, B_full] = ...
                    WEC_Core_Functions.interpolate_wamit_added_mass( ...
                        steel_data.vertical_shift, config, realised.CG_total(3));
            catch ME
                error('WEC_Shell_Offset:WAMITRealisedFailed', ...
                      ['interpolate_wamit_added_mass failed for the realised draft ', ...
                       'vs=%.4f m: %s. The hydro cache must cover the steel draft range.'], ...
                       steel_data.vertical_shift, ME.message);
            end

            realised.A11    = A11;
            realised.A33    = A33;     % translation-invariant
            realised.A55    = A55;
            realised.A_full = A_full;
            realised.B_full = B_full;

            % ── PTO stiffness — hull-independent, copy from optimiser ─
            if isfield(final_props, 'K_pto')
                realised.K_pto = final_props.K_pto;
            else
                realised.K_pto = zeros(3, 3);
            end

            % ── Total stiffness (K_hydro + K_pto, including K13 cross-term) ─
            K11_pto = realised.K_pto(1, 1);
            K33_pto = realised.K_pto(2, 2);
            K55_pto = realised.K_pto(3, 3);
            K13_pto = realised.K_pto(1, 3);
            realised.K_total = [K11_pto,         0,            K13_pto; ...
                                0,               K33_hydro+K33_pto, 0; ...
                                K13_pto,         0,            K55_hydro+K55_pto];

            % ── Uncoupled natural periods ────────────────────────────
            M11_virtual = realised.mass_total + A11;
            M33_virtual = realised.mass_total + A33;
            M55_virtual = realised.Iyy        + A55;

            K11_total = K11_pto;
            K33_total = realised.K_total(2, 2);
            K55_total = realised.K_total(3, 3);

            if K11_total > 1e-6
                realised.periods.surge = 2 * pi * sqrt(M11_virtual / K11_total);
            else
                realised.periods.surge = inf;
            end
            % Heave/pitch already computed by the steel solver — use those
            % directly so we don't re-derive and risk a tiny numerical drift.
            realised.periods.heave = steel_data.T_heave_realised;
            realised.periods.pitch = steel_data.T_pitch_realised;

            % ── Coupled eigenvalue analysis (3-DOF: surge, heave, pitch) ─
            try
                M_phys  = diag([realised.mass_total, realised.mass_total, realised.Iyy]);
                M_total_3x3 = M_phys + A_full;
                if det(M_total_3x3) > 1e-12 && det(realised.K_total) > 1e-12
                    [V_eig, D_eig] = eig(realised.K_total, M_total_3x3);
                    omega_sq = diag(D_eig);
                    valid = omega_sq > 1e-6;
                    if any(valid)
                        omega_n = sqrt(omega_sq(valid));
                        T_n = 2 * pi ./ omega_n;
                        [T_n_sorted, sort_idx] = sort(T_n, 'descend');
                        realised.coupled_periods = T_n_sorted;
                        V_valid = V_eig(:, valid);
                        realised.coupled_modes = V_valid(:, sort_idx);
                        realised.participation_factors = ...
                            WEC_Shell_Offset.compute_participation_factors( ...
                                realised.coupled_modes, M_total_3x3);
                    else
                        realised.coupled_periods = [realised.periods.surge; ...
                                                     realised.periods.heave; ...
                                                     realised.periods.pitch];
                        realised.coupled_modes = eye(3);
                        realised.participation_factors = 100 * eye(3);
                    end
                else
                    realised.coupled_periods = [realised.periods.surge; ...
                                                 realised.periods.heave; ...
                                                 realised.periods.pitch];
                    realised.coupled_modes = eye(3);
                    realised.participation_factors = 100 * eye(3);
                end
            catch
                realised.coupled_periods = [realised.periods.surge; ...
                                             realised.periods.heave; ...
                                             realised.periods.pitch];
                realised.coupled_modes = eye(3);
                realised.participation_factors = 100 * eye(3);
            end

            % ── 6×6 mass matrices (centre and origin) ────────────────
            realised.MassMatrix_CG = WEC_Core_Functions.calculate6x6MassMatrix( ...
                realised.mass_total, [0, 0, 0], realised.Inertia_Tensor);
            realised.MassMatrix_Origin = WEC_Core_Functions.calculate6x6MassMatrix( ...
                realised.mass_total, realised.CG_total, realised.Inertia_Tensor);

            % ── Realised per-strip mass partition ─────────────────────
            % If steel_data carries per-strip realisation arrays (UHPC path
            % via solve_constructable, or Phase-2 extraction in the steel
            % path), expose them on final_props.realised_strips so
            % diagnostics, visualisation and IO can read the AS-BUILT state
            % directly without recomputing.

            realised.realisation_mode       = 'steel_fill';
            realised.density_profile_source = 'realised_partition';

            rs = struct();
            rs.realisation_mode = realised.realisation_mode;
            rs.t_offset = [];
            rs.is_solid = [];
            rs.is_wall  = [];
            rs.z_lo     = [];
            rs.z_hi     = [];
            rs.V_uhpc   = [];
            rs.V_void   = [];
            rs.mass_uhpc = [];
            rs.mass_void = [];
            rs.uhpc_volume_fraction = [];
            rs.contours_outer = {};
            rs.contours_inner = {};

            if isfield(steel_data, 't_offset_strip') && ~isempty(steel_data.t_offset_strip)
                rs.t_offset = steel_data.t_offset_strip(:);
            end
            if isfield(steel_data, 'is_solid_strip') && ~isempty(steel_data.is_solid_strip)
                rs.is_solid = steel_data.is_solid_strip(:);
            end
            if isfield(steel_data, 'wall_strip_idx') && ~isempty(steel_data.wall_strip_idx)
                N_rs = max([length(rs.t_offset), length(rs.is_solid), 0]);
                if N_rs > 0
                    rs.is_wall = false(N_rs, 1);
                    rs.is_wall(steel_data.wall_strip_idx) = true;
                end
            end
            if isfield(steel_data, 'strip_z_lo')
                rs.z_lo = steel_data.strip_z_lo(:);
            end
            if isfield(steel_data, 'strip_z_hi')
                rs.z_hi = steel_data.strip_z_hi(:);
            end
            if isfield(steel_data, 'strip_V_UHPC')
                rs.V_uhpc = steel_data.strip_V_UHPC(:);
            end
            if isfield(steel_data, 'strip_V_void')
                rs.V_void = steel_data.strip_V_void(:);
            end
            if isfield(steel_data, 'strip_mass_UHPC')
                rs.mass_uhpc = steel_data.strip_mass_UHPC(:);
            end
            if isfield(steel_data, 'strip_mass_void')
                rs.mass_void = steel_data.strip_mass_void(:);
            end
            if ~isempty(rs.V_uhpc) && ~isempty(rs.V_void)
                Vt = rs.V_uhpc + rs.V_void;
                Vt(Vt <= 0) = 1;
                rs.uhpc_volume_fraction = rs.V_uhpc ./ Vt;
            end
            if isfield(steel_data, 'contours_outer')
                rs.contours_outer = steel_data.contours_outer;
            end
            if isfield(steel_data, 'contours_inner')
                rs.contours_inner = steel_data.contours_inner;
            end
            realised.realised_strips = rs;

            % Backward-compat carry-overs (cross_section is the hull profile,
            % not material; densities_at_nodes is the optimiser's per-node rho —
            % we keep both for legacy plotting consumers).
            if isfield(final_props, 'cross_section')
                realised.cross_section = final_props.cross_section;
            else
                realised.cross_section = [];
            end
            if isfield(final_props, 'densities_at_nodes')
                realised.densities_at_nodes = final_props.densities_at_nodes;
            else
                realised.densities_at_nodes = [];
            end

            % Realised per-strip equivalent density (set by either pipeline):
            %   Steel:  steel_data.strip_rho_eff  (computed in solve §8b)
            %   UHPC :  cstr.strip_rho_eff        (computed in
            %           WEC_Constructable_Hull.extract_strip_geometry)
            % Exposed at the top level of realised so visualisation
            % functions can render the AS-BUILT cake layers without
            % digging into solver-specific sub-structs.
            if isfield(steel_data, 'strip_rho_eff') && ~isempty(steel_data.strip_rho_eff)
                realised.realised_strip_density = steel_data.strip_rho_eff(:);
            else
                realised.realised_strip_density = [];
            end
            if isfield(steel_data, 'strip_edges') && ~isempty(steel_data.strip_edges)
                realised.realised_strip_edges = steel_data.strip_edges(:);
            elseif isfield(config, 'strip_edges')
                realised.realised_strip_edges = config.strip_edges(:);
            else
                realised.realised_strip_edges = [];
            end
            if isfield(steel_data, 'realisation_mode') && ~isempty(steel_data.realisation_mode)
                realised.realisation_mode = steel_data.realisation_mode;
            else
                realised.realisation_mode = '';
            end

            % `components`: REPLACE the optimiser's per-strip rho with the
            % realised UHPC/void partition where available, so any consumer
            % that reads final_props.components sees as-built data.
            if ~isempty(rs.mass_uhpc)
                N = length(rs.mass_uhpc);
                comps = repmat(struct('density', 0, 'z_level', 0, ...
                                      'V_uhpc', 0, 'V_void', 0, ...
                                      'mass_uhpc', 0, 'mass_void', 0, ...
                                      'is_solid', false, 'is_wall', false, ...
                                      't_offset', NaN), N, 1);
                for ii = 1:N
                    if ~isempty(rs.z_lo) && ~isempty(rs.z_hi)
                        z_mid = 0.5 * (rs.z_lo(ii) + rs.z_hi(ii));
                    else
                        z_mid = 0;
                    end
                    Vti = 0;
                    if ~isempty(rs.V_uhpc), Vti = Vti + rs.V_uhpc(ii); end
                    if ~isempty(rs.V_void), Vti = Vti + rs.V_void(ii); end
                    Mti = rs.mass_uhpc(ii);
                    if ~isempty(rs.mass_void), Mti = Mti + rs.mass_void(ii); end
                    if Vti > 1e-12
                        comps(ii).density = Mti / Vti;
                    end
                    comps(ii).z_level = z_mid;
                    if ~isempty(rs.V_uhpc), comps(ii).V_uhpc = rs.V_uhpc(ii); end
                    if ~isempty(rs.V_void), comps(ii).V_void = rs.V_void(ii); end
                    comps(ii).mass_uhpc = rs.mass_uhpc(ii);
                    if ~isempty(rs.mass_void), comps(ii).mass_void = rs.mass_void(ii); end
                    if ~isempty(rs.is_solid), comps(ii).is_solid = rs.is_solid(ii); end
                    if ~isempty(rs.is_wall),  comps(ii).is_wall  = rs.is_wall(ii);  end
                    if ~isempty(rs.t_offset), comps(ii).t_offset = rs.t_offset(ii); end
                end
                realised.components = comps;
            elseif isfield(final_props, 'components')
                realised.components = final_props.components;
            else
                realised.components = struct('density', {}, 'z_level', {});
            end
        end


        function PF = compute_participation_factors(modes, M)
            % Local helper — mirrors calculate_3d_properties.compute_participation_factors.
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


        %% ═════════════════════════════════════════════════════════════
        %%  INNER DRAFT SOLVE  (mass-balance enforcement)
        %% ═════════════════════════════════════════════════════════════

        function [d_star, vs_star, mass_err_abs] = solve_draft_for_mass(M_total, cfg, draft_warm)
        % SOLVE_DRAFT_FOR_MASS  Find draft d so that ρ_w · V_sub(d) = M_total.
        %
        %   Returns NaN, NaN, NaN if no draft in the feasible bracket
        %   satisfies the equation (M_total exceeds full-displacement
        %   buoyancy or falls below freeboard buoyancy).

            rho_w      = cfg.RHO_WATER;
            hull_z_min = cfg.hull_z_min;
            hull_z_max = cfg.hull_z_max;

            % g(d) = ρ_w·V_sub(z_wl(d)) − M_total      where z_wl = hull_z_min + draft (world frame waterline at z=0,
            %                                          but z_wl_body uses the body-frame draft def from existing code).
            %
            % From calculate_3d_properties.m:55:
            %   props.draft = abs(hull_z_min + props.vertical_shift)
            %   z_wl_body   = -vertical_shift
            % so for draft d > 0 (with hull_z_min < 0):
            %   vertical_shift = -(d + hull_z_min)        --> note hull_z_min < 0
            %   z_wl_body      = d + hull_z_min
            %
            % Equivalently z_wl_body ranges from hull_z_min (d=0, hull lifted out
            % of water) to hull_z_max (d = hull_z_max - hull_z_min, fully submerged).
            % V_sub(z_wl_body) is monotone non-decreasing on this range
            % (per the precomputed config.V_sub_table).

            d_min = 0;
            d_max = hull_z_max - hull_z_min;        % fully submerged
            d_lo  = max(d_min,  draft_warm - 0.5 * (hull_z_max - hull_z_min));
            d_hi  = min(d_max,  draft_warm + 0.5 * (hull_z_max - hull_z_min));
            if d_hi <= d_lo
                d_lo = d_min;
                d_hi = d_max;
            end

            g = @(d) rho_w * WEC_Shell_Offset.V_sub_at_draft(d, cfg) - M_total;

            g_lo = g(d_lo);
            g_hi = g(d_hi);

            % Expand bracket if the warm-start window doesn't straddle 0
            if g_lo * g_hi > 0
                d_lo = d_min;
                d_hi = d_max;
                g_lo = g(d_lo);
                g_hi = g(d_hi);
            end

            if g_lo > 0
                % Even at zero draft (hull lifted out of water) buoyancy already
                % exceeds the steel-air mass — impossible (the hull is too light).
                % g(d_lo) = ρ_w·V_sub(0) - M_total = 0 - M_total = -M_total <= 0,
                % so this branch normally cannot trigger; it would only fire if
                % the table extends V_sub > 0 at d = 0 (curtain effect). Treat
                % as infeasible.
                d_star = NaN;  vs_star = NaN;  mass_err_abs = NaN;
                return;
            end
            if g_hi < 0
                % Fully submerged the buoyancy is still less than M_total → infeasible.
                d_star = NaN;  vs_star = NaN;  mass_err_abs = NaN;
                return;
            end

            opts_fz = optimset('Display','off','TolX',1e-7);
            try
                d_star = fzero(g, [d_lo, d_hi], opts_fz);
            catch ME
                error('WEC_Shell_Offset:DraftSolveFailed', ...
                      'fzero failed in inner draft solve: %s', ME.message);
            end

            vs_star = -(d_star + hull_z_min);
            mass_err_abs = abs(g(d_star));
        end


        %% ═════════════════════════════════════════════════════════════
        %%  SUBMERGED VOLUME LOOK-UP (config tables)
        %% ═════════════════════════════════════════════════════════════

        function V = V_sub_at_draft(d, cfg)
        % V_SUB_AT_DRAFT  Submerged volume from precomputed table at draft d.
        %
        %   Body-frame waterline z_wl = d + hull_z_min  (hull_z_min < 0 typical).
        %   Outside the table range, interp1(...,'linear',0) returns 0.

            z_wl_body = d + cfg.hull_z_min;
            z_top     = min(z_wl_body, cfg.hull_z_max);
            V = max(0, interp1(cfg.Aw_table_z, cfg.V_sub_table, z_top, 'linear', 0));
        end


        %% ═════════════════════════════════════════════════════════════
        %%  GEOMETRY GRID  (build A_outer, A_inner, Iyy_outer, Iyy_inner)
        %% ═════════════════════════════════════════════════════════════

        function [grids, n_zero_interior] = build_geometry_grid( ...
                config, t_steel, n_z, max_slope_factor)
        % BUILD_GEOMETRY_GRID  Sample (A_outer, A_inner, Iyy_outer, Iyy_inner) on a z-grid.
        %
        %   Robust path: A_outer and Iyy_outer come from the precomputed
        %   config.Aw_table / I_wp_yy_table (both built by the same
        %   WEC_HydroProperties.extract_isocurve_at_z + waterplane_properties
        %   pipeline used everywhere else in the suite).  A_inner / Iyy_inner
        %   are computed per-z by:
        %     1. extract_isocurve_at_z  (same boundary cache as the builder)
        %     2. waterplane_properties  (angular-sort the contour)
        %     3. offset_vertices_raw → inner miter polygon
        %     4. polygeom on inner polygon for A_inner, Iyy_inner, Ixx_inner
        %
        %   The "evaluateCrossSectionMS2" path used by the original Shell_Offset
        %   was retired in v2.0 — it failed silently at majority of z-grid
        %   points on hulls with multi-surface topology (column joining
        %   platform to keel), even where Aw_table reported non-zero area.
        %
        %   n_zero_interior counts INTERIOR z's where A_outer is genuinely zero
        %   per the Aw table (i.e. above hull_z_max or below hull_z_min after
        %   surface-range trimming).  These are EXPECTED boundary degeneracies,
        %   not failures.

            persistent cache_key cache_grids

            % Hard error if the parser doesn't expose a filename — silently
            % using a constant key would let stale grids leak across runs.
            ms2_model = config.ms2_model;
            if ~isprop(ms2_model, 'filename') && ~isfield(ms2_model, 'filename')
                error('WEC_Shell_Offset:NoParserFilename', ...
                      'ms2_model has no .filename property — cannot key cache.');
            end
            parser_key = ms2_model.filename;
            if isempty(parser_key)
                error('WEC_Shell_Offset:EmptyParserFilename', ...
                      'ms2_model.filename is empty.  Cannot key cache.');
            end

            % Cache key includes the Aw_table identity (length + endpoints)
            % so a builder-side change to the table invalidates the cache.
            tab_z = config.Aw_table_z(:);
            this_key = sprintf( ...
                '%s|t=%.10g|n=%d|s=%.6g|zN=%d|z0=%.6g|zE=%.6g|hlo=%.6g|hhi=%.6g', ...
                parser_key, t_steel, n_z, max_slope_factor, ...
                length(tab_z), tab_z(1), tab_z(end), ...
                config.hull_z_min, config.hull_z_max);
            if ~isempty(cache_key) && strcmp(cache_key, this_key)
                grids = cache_grids;
                n_zero_interior = sum(grids.A_outer(2:end-1) <= 0);
                return;
            end

            % Validate boundary cache once (defensive — solve() already checks)
            if ~isfield(config, 'boundary_cache') || ...
                    ~isstruct(config.boundary_cache) || ...
                    ~isfield(config.boundary_cache, 'sources') || ...
                    ~isfield(config.boundary_cache, 'u_samples')
                error('WEC_Shell_Offset:BadBoundaryCache', ...
                      'config.boundary_cache is missing required fields ("sources"/"u_samples").');
            end
            n_u = length(config.boundary_cache.u_samples);

            %% --- Outer geometry: vectorised lookup from precomputed tables ---
            z = linspace(config.hull_z_min, config.hull_z_max, n_z)';
            A_outer   = max(0, interp1(config.Aw_table_z, config.Aw_table,      z, 'linear', 0));
            Iyy_outer = max(0, interp1(config.Aw_table_z, config.I_wp_yy_table, z, 'linear', 0));
            % Ixx_outer (∫∫ y² dA) needed for the realised hull's roll/yaw inertia.
            % I_wp_xx_table is built by Configuration_Builder §2c on the same z-grid.
            if isfield(config, 'I_wp_xx_table') && ~isempty(config.I_wp_xx_table)
                Ixx_outer = max(0, interp1(config.Aw_table_z, config.I_wp_xx_table, z, 'linear', 0));
            else
                % Symmetric-hull fallback: Ixx_cross == Iyy_cross when the cross-section
                % is symmetric in x↔y (true for revolution surfaces and for E1).
                Ixx_outer = Iyy_outer;
            end

            A_inner   = zeros(n_z, 1);
            Iyy_inner = zeros(n_z, 1);
            Ixx_inner = zeros(n_z, 1);

            n_zero_interior = 0;
            for k = 1:n_z
                if A_outer(k) <= 1e-10
                    if k > 1 && k < n_z
                        n_zero_interior = n_zero_interior + 1;
                    end
                    continue;
                end

                [A_inner(k), Iyy_inner(k), Ixx_inner(k)] = ...
                    WEC_Shell_Offset.inner_props_at_z( ...
                        config, z(k), t_steel, max_slope_factor, ...
                        n_u, A_outer(k), Iyy_outer(k), Ixx_outer(k));
            end

            % Convex-polygon noise guards
            A_inner   = min(A_inner,   A_outer);
            A_inner   = max(A_inner,   0);
            Iyy_inner = min(Iyy_inner, Iyy_outer);
            Iyy_inner = max(Iyy_inner, 0);
            Ixx_inner = min(Ixx_inner, Ixx_outer);
            Ixx_inner = max(Ixx_inner, 0);

            grids = struct('z', z, ...
                           'A_outer',   A_outer, ...
                           'A_inner',   A_inner, ...
                           'Iyy_outer', Iyy_outer, ...
                           'Iyy_inner', Iyy_inner, ...
                           'Ixx_outer', Ixx_outer, ...
                           'Ixx_inner', Ixx_inner);

            cache_key   = this_key;
            cache_grids = grids;
        end


        %% ═════════════════════════════════════════════════════════════
        %%  INNER PROPS AT ONE Z  (polygeom on miter polygon — robust path)
        %% ═════════════════════════════════════════════════════════════

        function [A_inner, Iyy_inner, Ixx_inner] = inner_props_at_z( ...
                config, z_level, t_steel, max_slope_factor, n_u, ...
                A_outer, Iyy_outer, Ixx_outer)
        % INNER_PROPS_AT_Z  Inner-area + Iyy + Ixx at one z, given non-empty A_outer.
        %
        %   Uses extract_isocurve_at_z (the same robust path the builder
        %   uses for the Aw_table) to get the cross-section polygon, then
        %   offsets it inward via offset_vertices_raw and calls polygeom.
        %   Returns BOTH cross-section second moments:
        %     Iyy_inner = ∫∫ x² dA   (pitch axis)
        %     Ixx_inner = ∫∫ y² dA   (roll axis)
        %   so the caller can build a full inertia tensor for the realised hull.

            wl_pts = WEC_HydroProperties.extract_isocurve_at_z( ...
                         config.ms2_model, z_level, n_u, config.boundary_cache);

            if isempty(wl_pts) || size(wl_pts, 1) < 3
                error('WEC_Shell_Offset:NoIsocurveAtZ', ...
                      ['extract_isocurve_at_z returned <3 contour points at ', ...
                       'z=%g where Aw_table reports A_outer=%g (>0).  ', ...
                       'Check that config.boundary_cache and config.Aw_table* ', ...
                       'were built for the same ms2_model.'], z_level, A_outer);
            end

            [~, ~, ~, pts_ord] = WEC_HydroProperties.waterplane_properties(wl_pts);
            x_poly = pts_ord(:, 1);
            y_poly = pts_ord(:, 2);

            x_cl = [x_poly; x_poly(1)];
            y_cl = [y_poly; y_poly(1)];
            P_outer = sum(sqrt(diff(x_cl).^2 + diff(y_cl).^2));

            cos_alpha = WEC_Shell_Offset.slope_cos_from_table( ...
                            config, z_level, A_outer, P_outer);
            cos_alpha = max(cos_alpha, 1.0 / max_slope_factor);
            offset_dist = t_steel / cos_alpha;

            % Use polygeom on the miter polygon for BOTH area and second moments
            % (consistent: same boundary for mass and inertia integrals).
            [x_off, y_off] = WEC_Shell_Offset.offset_vertices_raw( ...
                                  x_poly, y_poly, offset_dist);
            if length(x_off) >= 3
                [geom_in, iner_in, ~] = WEC_Core_Functions.polygeom(x_off, y_off);
                A_miter = geom_in(1);
                % Validity guard: collapsed polygon (A<=0) or acute-corner overshoot
                % (A_miter >= A_outer occurs when interior angles are acute and the
                % miter vertex moves past the outer boundary — physically invalid).
                if A_miter <= 1e-10 || A_miter >= A_outer
                    A_inner = 0;  Iyy_inner = 0;  Ixx_inner = 0;
                    return;
                end
                A_inner   = A_miter;
                Iyy_inner = abs(iner_in(2));   % ∫∫ x² dA about origin
                Ixx_inner = abs(iner_in(1));   % ∫∫ y² dA about origin
                if ~isfinite(Iyy_inner)
                    Iyy_inner = Iyy_outer * (A_inner / A_outer)^2;
                end
                if ~isfinite(Ixx_inner)
                    Ixx_inner = Ixx_outer * (A_inner / A_outer)^2;
                end
            else
                A_inner = 0;  Iyy_inner = 0;  Ixx_inner = 0;
                return;
            end

            Iyy_inner = min(Iyy_inner, Iyy_outer);
            Ixx_inner = min(Ixx_inner, Ixx_outer);
        end


        %% ═════════════════════════════════════════════════════════════
        %%  SLOPE COS FROM Aw_TABLE  (mesh-free, table-based finite difference)
        %% ═════════════════════════════════════════════════════════════

        function cos_alpha = slope_cos_from_table(config, z, A_here, P_here)
        % SLOPE_COS_FROM_TABLE  cos(alpha) of hull surface at z, from Aw_table.
        %
        %   Uses dA/dz ≈ P · dr/dz   (uniform radial expansion model)
        %   →  tan(alpha) = |dA/dz| / P,   cos(alpha) = 1 / sqrt(1 + tan²)
        %
        %   dA/dz is a centred finite difference on config.Aw_table at z;
        %   P_here is the actual cross-section perimeter from the polygon
        %   the caller has already extracted.

            dz = 0.01;

            A_lo = max(0, interp1(config.Aw_table_z, config.Aw_table, z - dz, 'linear', 0));
            A_hi = max(0, interp1(config.Aw_table_z, config.Aw_table, z + dz, 'linear', 0));

            if A_here < 1e-10 || P_here < 1e-10
                cos_alpha = 1.0;
                return;
            end

            % If we're at the top of the hull (Aw vanishes one dz above), assume
            % a near-horizontal cap (the Steiner correction will saturate).
            if A_hi < 1e-10 && A_lo > 1e-10
                cos_alpha = 0.1;
                return;
            end

            % Centred FD where possible; fall back to one-sided at table edges.
            if A_hi > 0 && A_lo > 0
                dA_dz = (A_hi - A_lo) / (2 * dz);
            elseif A_hi > 0
                dA_dz = (A_hi - A_here) / dz;
            else
                dA_dz = (A_here - A_lo) / dz;
            end

            dr_dz     = abs(dA_dz) / max(P_here, 1e-10);
            cos_alpha = 1.0 / sqrt(1.0 + dr_dz^2);
            cos_alpha = max(cos_alpha, 0.01);
        end


        %% ═════════════════════════════════════════════════════════════
        %%  INTEGRATE SPLIT  (trapezoidal with z_fill as extra knot)
        %% ═════════════════════════════════════════════════════════════

        function R = integrate_split(grids, z_fill)
        % INTEGRATE_SPLIT  Trapezoidal volume / 1st-moment / 2nd-moment splits
        %   below and above z_fill, separately for A_outer, A_inner, A_jacket.
        %
        %   Inserts z_fill as an additional integration knot via linear
        %   interpolation, so the split is exact to leading order in dz.

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

            % Clamp z_fill to grid range
            z_fill = max(z(1), min(z(end), z_fill));

            below_mask = z <= z_fill;
            above_mask = z >= z_fill;
            z_below = z(below_mask);
            z_above = z(above_mask);

            tol = 1e-12;
            if isempty(z_below) || abs(z_below(end) - z_fill) > tol
                A_o_zf   = interp1(z, A_o,   z_fill, 'linear');
                A_i_zf   = interp1(z, A_i,   z_fill, 'linear');
                A_j_zf   = interp1(z, A_j,   z_fill, 'linear');
                Iyy_o_zf = interp1(z, Iyy_o, z_fill, 'linear');
                Iyy_i_zf = interp1(z, Iyy_i, z_fill, 'linear');
                Iyy_j_zf = interp1(z, Iyy_j, z_fill, 'linear');
                Ixx_o_zf = interp1(z, Ixx_o, z_fill, 'linear');
                Ixx_i_zf = interp1(z, Ixx_i, z_fill, 'linear');
                Ixx_j_zf = interp1(z, Ixx_j, z_fill, 'linear');

                z_below   = [z_below;  z_fill];
                A_o_below = [A_o(below_mask);  A_o_zf];
                A_i_below = [A_i(below_mask);  A_i_zf];
                A_j_below = [A_j(below_mask);  A_j_zf];
                Iyy_o_below = [Iyy_o(below_mask); Iyy_o_zf];
                Iyy_i_below = [Iyy_i(below_mask); Iyy_i_zf];
                Iyy_j_below = [Iyy_j(below_mask); Iyy_j_zf];
                Ixx_o_below = [Ixx_o(below_mask); Ixx_o_zf];
                Ixx_i_below = [Ixx_i(below_mask); Ixx_i_zf];
                Ixx_j_below = [Ixx_j(below_mask); Ixx_j_zf];

                z_above   = [z_fill;     z_above];
                A_o_above = [A_o_zf;     A_o(above_mask)];
                A_i_above = [A_i_zf;     A_i(above_mask)];
                A_j_above = [A_j_zf;     A_j(above_mask)];
                Iyy_o_above = [Iyy_o_zf; Iyy_o(above_mask)];
                Iyy_i_above = [Iyy_i_zf; Iyy_i(above_mask)];
                Iyy_j_above = [Iyy_j_zf; Iyy_j(above_mask)];
                Ixx_o_above = [Ixx_o_zf; Ixx_o(above_mask)];
                Ixx_i_above = [Ixx_i_zf; Ixx_i(above_mask)];
                Ixx_j_above = [Ixx_j_zf; Ixx_j(above_mask)];
            else
                A_o_below = A_o(below_mask); A_i_below = A_i(below_mask); A_j_below = A_j(below_mask);
                Iyy_o_below = Iyy_o(below_mask); Iyy_i_below = Iyy_i(below_mask); Iyy_j_below = Iyy_j(below_mask);
                Ixx_o_below = Ixx_o(below_mask); Ixx_i_below = Ixx_i(below_mask); Ixx_j_below = Ixx_j(below_mask);
                A_o_above = A_o(above_mask); A_i_above = A_i(above_mask); A_j_above = A_j(above_mask);
                Iyy_o_above = Iyy_o(above_mask); Iyy_i_above = Iyy_i(above_mask); Iyy_j_above = Iyy_j(above_mask);
                Ixx_o_above = Ixx_o(above_mask); Ixx_i_above = Ixx_i(above_mask); Ixx_j_above = Ixx_j(above_mask);
            end

            R = struct();

            % --- BELOW z_fill: cross-section is fully solid (use A_outer everywhere) ---
            R.V_below_outer       = WEC_Shell_Offset.tz(z_below, A_o_below);
            R.int_zA_below_outer  = WEC_Shell_Offset.tz(z_below, z_below .* A_o_below);
            R.int_z2A_below_outer = WEC_Shell_Offset.tz(z_below, z_below.^2 .* A_o_below);
            R.int_Ix_below_outer  = WEC_Shell_Offset.tz(z_below, Iyy_o_below);   % ∫∫∫ x² dV
            R.int_Iy_below_outer  = WEC_Shell_Offset.tz(z_below, Ixx_o_below);   % ∫∫∫ y² dV

            % --- ABOVE z_fill: jacket annulus is steel; inner area is air ---
            R.V_above_jacket       = WEC_Shell_Offset.tz(z_above, A_j_above);
            R.int_zA_above_jacket  = WEC_Shell_Offset.tz(z_above, z_above .* A_j_above);
            R.int_z2A_above_jacket = WEC_Shell_Offset.tz(z_above, z_above.^2 .* A_j_above);
            R.int_Ix_above_jacket  = WEC_Shell_Offset.tz(z_above, Iyy_j_above);
            R.int_Iy_above_jacket  = WEC_Shell_Offset.tz(z_above, Ixx_j_above);

            R.V_above_inner       = WEC_Shell_Offset.tz(z_above, A_i_above);
            R.int_zA_above_inner  = WEC_Shell_Offset.tz(z_above, z_above .* A_i_above);
            R.int_z2A_above_inner = WEC_Shell_Offset.tz(z_above, z_above.^2 .* A_i_above);
            R.int_Ix_above_inner  = WEC_Shell_Offset.tz(z_above, Iyy_i_above);
            R.int_Iy_above_inner  = WEC_Shell_Offset.tz(z_above, Ixx_i_above);
        end


        function v = tz(x, y)
            % Safe trapz: returns 0 when fewer than 2 sample points.
            if length(x) < 2
                v = 0;
            else
                v = trapz(x, y);
            end
        end


        %% ═════════════════════════════════════════════════════════════
        %%  SIGMOID (and its inverse) for bound-encoded parameters
        %% ═════════════════════════════════════════════════════════════

        function s = sigmoid(x)
            s = 1 ./ (1 + exp(-x));
        end

        function x = inv_sigmoid(s)
            s = max(min(s, 1 - 1e-12), 1e-12);
            x = log(s ./ (1 - s));
        end


        %% ═════════════════════════════════════════════════════════════
        %%  OPTIONS HELPER
        %% ═════════════════════════════════════════════════════════════

        function v = opt_or_cfg(opts, name, cfg_value)
            % Use opts.(name) if present and non-empty; otherwise use cfg_value.
            if isfield(opts, name) && ~isempty(opts.(name))
                v = opts.(name);
            else
                v = cfg_value;
            end
        end


        %% ═════════════════════════════════════════════════════════════
        %%  RAW VERTEX-NORMAL OFFSET (KEPT from v1.0)
        %% ═════════════════════════════════════════════════════════════

        function [x_off, y_off] = offset_vertices_raw(x, y, dist)
        % OFFSET_VERTICES_RAW  Shrink convex polygon by constant distance.
        %
        %   For convex polygons, vertex-normal offset with miter correction
        %   never self-intersects; raw offset vertices in input order are
        %   exactly what polygeom needs for an Iyy estimate.

            x = x(:);  y = y(:);
            Nv = length(x);

            if Nv < 3 || dist <= 0
                x_off = x;  y_off = y;
                return;
            end

            % Remove duplicate closing vertex if present
            if abs(x(end)-x(1)) < 1e-12 && abs(y(end)-y(1)) < 1e-12
                x = x(1:end-1);
                y = y(1:end-1);
                Nv = length(x);
            end
            if Nv < 3
                x_off = [];  y_off = [];
                return;
            end

            % Ensure CCW winding
            signed_area = 0.5 * sum(x .* circshift(y,-1) - circshift(x,-1) .* y);
            if signed_area < 0
                x = flipud(x);
                y = flipud(y);
            end

            x_off = zeros(Nv, 1);
            y_off = zeros(Nv, 1);

            for j = 1:Nv
                jm = mod(j-2, Nv) + 1;
                jp = mod(j,   Nv) + 1;

                e_prev = [x(j)-x(jm), y(j)-y(jm)];
                e_next = [x(jp)-x(j), y(jp)-y(j)];

                len_prev = norm(e_prev);
                len_next = norm(e_next);

                if len_prev < 1e-12 || len_next < 1e-12
                    x_off(j) = x(j);
                    y_off(j) = y(j);
                    continue;
                end

                % Inward normals (CCW: left of edge)
                n_prev = [-e_prev(2), e_prev(1)] / len_prev;
                n_next = [-e_next(2), e_next(1)] / len_next;

                n_avg = n_prev + n_next;
                len_avg = norm(n_avg);
                if len_avg < 1e-12
                    n_avg = n_prev;
                    len_avg = 1.0;
                end
                n_avg = n_avg / len_avg;

                % Miter correction
                cos_half = dot(n_avg, n_prev);
                if cos_half > 0.33
                    miter_dist = dist / cos_half;
                else
                    miter_dist = dist * 3.0;
                end

                x_off(j) = x(j) + n_avg(1) * miter_dist;
                y_off(j) = y(j) + n_avg(2) * miter_dist;
            end
        end


    end  % methods (Static)
end  % classdef
