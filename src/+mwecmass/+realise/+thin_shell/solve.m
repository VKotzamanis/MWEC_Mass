function steel_data = solve(config, x_opt_3d, final_props, opts, fids)
%SOLVE Optimise thin-shell vertical shift, thickness, and fill elevation.
% steel_data = solve(config,x_opt_3d,final_props,opts,fids) runs constrained
% fmincon SQP for x=[vs;t_steel;z_fill], enforcing buoyancy and GM constraints
% with a heave/pitch range objective. config supplies geometry/material data;
% opts overrides solver defaults and fids receives logs. Numeric results use SI.
% See docs/METHODS_ENGINE.md#thin-shell-realisation-optimisation

    if nargin < 5 || isempty(fids), fids = 1; end
    t_solve_start = tic;
    mwecmass.output.emit(fids, '\n    mwecmass.realise.thin_shell.solve (fmincon SQP, constrained):\n');

    %% Input validation
    if nargin < 3
        error('mwecmass:thin_shell:NotEnoughInputs', ...
              'solve(config, x_opt_3d, final_props[, opts]) requires 3+ args.');
    end
    if nargin < 4 || isempty(opts), opts = struct(); end

    required_cfg = {'ms2_model','hull_z_min','hull_z_max', ...
                    'Aw_table_z','Aw_table','V_sub_table','CB_z_table', ...
                    'I_wp_yy_table','RHO_WATER','G', ...
                    'rho_shell','rho_fill','rho_air','steel_t_init','steel_t_min', ...
                    'steel_max_slope_factor','steel_n_z_grid', ...
                    'vertical_shift_bounds','gm_min', ...
                    'T_heave_goal','T_heave_range', ...
                    'T_pitch_goal','T_pitch_range','zone_k_amp'};
    for k = 1:length(required_cfg)
        if ~isfield(config, required_cfg{k}) || isempty(config.(required_cfg{k}))
            error('mwecmass:thin_shell:MissingConfig', ...
                  'Required config field "%s" is missing or empty.', ...
                  required_cfg{k});
        end
    end

    x_opt_3d = x_opt_3d(:);
    if length(x_opt_3d) < 2
        error('mwecmass:thin_shell:BadOptVector', ...
              'x_opt_3d must be [vertical_shift; rho_1; ...] (length >= 2).');
    end

    %% Resolve options
    % rho_shell is wall density; rho_fill is the density below z_fill and
    % defaults to rho_shell. opts can override these and solver settings.
    rho_shell        = mwecmass.internal.option_or_config(opts, 'rho_shell',        config.rho_shell);
    rho_fill         = mwecmass.internal.option_or_config(opts, 'rho_fill',         config.rho_fill);
    rho_air          = mwecmass.internal.option_or_config(opts, 'rho_air',          config.rho_air);
    t_init           = mwecmass.internal.option_or_config(opts, 't_init',           config.steel_t_init);
    t_min            = mwecmass.internal.option_or_config(opts, 't_min',            config.steel_t_min);
    max_slope_factor = mwecmass.internal.option_or_config(opts, 'max_slope_factor', config.steel_max_slope_factor);
    n_z_grid         = mwecmass.internal.option_or_config(opts, 'n_z_grid',         config.steel_n_z_grid);

    assert(rho_shell > 0,        'rho_shell must be > 0 (got %g)',        rho_shell);
    assert(rho_fill  > 0,        'rho_fill must be > 0 (got %g)',         rho_fill);
    assert(rho_air   > 0,        'rho_air must be > 0 (got %g)',          rho_air);
    % Air must be lighter than both solid regions; the fill ordering below
    % also ensures a monotone cumulative integrand for the warm-start seed.
    assert(rho_air   < min(rho_shell, rho_fill), ...
           'rho_air (%g) must be < min(rho_shell, rho_fill) = %g', rho_air, min(rho_shell, rho_fill));
    assert(n_z_grid  >= 50,      'n_z_grid must be >= 50 (got %d)',       n_z_grid);
    assert(t_min     >= 0,       't_min must be >= 0 (got %g)',           t_min);
    % Enforce the fill ordering at the solver boundary because opts may
    % override the input values; this keeps the seed cumulative integral monotone.
    assert(rho_fill  >= rho_shell, ...
           ['rho_fill (%g) must be >= rho_shell (%g): the analytic z_fill seed''s cumulative ' ...
            'integrand is only guaranteed monotone under this ordering (Amendment A2)'], ...
           rho_fill, rho_shell);

    %% Hull bounds and thickness bracket

    hull_z_min = config.hull_z_min;
    hull_z_max = config.hull_z_max;
    assert(hull_z_max > hull_z_min, 'hull extents inverted (%g >= %g)', ...
           hull_z_min, hull_z_max);

    % Upper bound on t_steel:  half the minimum hull radius in the
    % middle 80% of the height (avoid the apex/keel-tip pinch points).
    n_probe = 9;
    z_probes = linspace(0.1*hull_z_min + 0.9*hull_z_max, ...
                        0.9*hull_z_min + 0.1*hull_z_max, n_probe);
    r_probes = zeros(n_probe, 1);

    if ~isfield(config, 'boundary_cache') || ...
            ~isstruct(config.boundary_cache) || ...
            ~isfield(config.boundary_cache, 'sources')
        error('mwecmass:thin_shell:NoBoundaryCache', ...
              ['config.boundary_cache is missing or malformed. ' ...
               'Configuration_Builder must populate it via ' ...
               'mwecmass.geometry.precompute_boundary_cache.']);
    end
    cache = config.boundary_cache;
    for k = 1:n_probe
        r_probes(k) = mwecmass.geometry.compute_rmin_at_z( ...
                          config.ms2_model, z_probes(k), 60, cache);
    end
    r_min_global = min(r_probes(r_probes > 0));
    if isempty(r_min_global) || r_min_global <= 0
        error('mwecmass:thin_shell:NoFiniteRadius', ...
              'compute_rmin_at_z returned no positive radius across %d probes.', ...
              n_probe);
    end
    t_max = 0.5 * r_min_global;
    if t_min >= t_max
        error('mwecmass:thin_shell:TMinExceedsTMax', ...
              ['t_min = %.5f m >= t_max = %.5f m — hull too narrow ', ...
               'for requested fabrication floor.'], t_min, t_max);
    end

    t_init_requested = t_init;
    t_init = max(min(t_init, 0.9*t_max), max(t_min, 1e-4));
    if abs(t_init - t_init_requested) > 1e-6
        warning('mwecmass:thin_shell:TInitClamped', ...
                'Initial t_steel guess clamped %.5f → %.5f m.', ...
                t_init_requested, t_init);
    end

    vs_lb = config.vertical_shift_bounds(1);
    vs_ub = config.vertical_shift_bounds(2);
    vs_opt = x_opt_3d(1);
    draft_opt = abs(hull_z_min + vs_opt);

    mwecmass.output.emit(fids, '      Hull z=[%.4f, %.4f] m  t∈[%.5f, %.4f] m  vs∈[%.3f, %.3f] m\n', ...
            hull_z_min, hull_z_max, t_min, 0.95*t_max, vs_lb, vs_ub);
    mwecmass.output.emit(fids, '      Targets (from optimiser): GM_min=%.3f m  T_h=%.3f s  T_p=%.3f s  M=%.1f kg\n', ...
            config.gm_min, config.T_heave_goal, config.T_pitch_goal, ...
            final_props.mass_total);
    mwecmass.output.emit(fids, '      Warm-start vs (from x_opt_3d): %.4f m  (draft=%.4f m)\n', ...
            vs_opt, draft_opt);

    %% Analytic z_fill seed at the optimiser draft and initial thickness
    %  Picks the z_fill that satisfies buoyancy balance for the warm-start
    %  hull mass at the optimiser's draft, so fmincon starts near
    %  feasibility on ceq. z_fill remains a free fmincon design variable
    %  z_fill remains a free design variable; this seed affects convergence only.
    %  The default (rho_fill == rho_shell) branch below inverts the
    %  single-density closed form exactly. The general
    %  branch's cumulative integrand g(z) is piecewise LINEAR between
    %  z_grid_seed nodes, so its cumulative integral is piecewise
    %  QUADRATIC; inverting it via interp1 on the node values of
    %  cumtrapz linearly interpolates a function that is
    %  actually quadratic between those nodes -- not exact. The general
    %  branch now inverts the exact quadratic on the bracketing segment
    %  (invert_piecewise_linear_cumulative, local function below). The
    %  z_fill_seed = hull_z_min+1e-3 / hull_z_max-1e-3 clamps at the two
    %  ends of the C_target range are BOUNDED FALLBACKS (nearest feasible
    %  seed to an out-of-range target), not exact solutions.
    V_sub_at_vs_opt = max(0, interp1(config.Aw_table_z, ...
                                      config.V_sub_table, ...
                                      -vs_opt, 'linear', 0));
    M_buoy_target = config.RHO_WATER * V_sub_at_vs_opt;

    [grids_seed, n_zero_seed] = mwecmass.realise.thin_shell.build_geometry_grid( ...
        config, t_init, n_z_grid, max_slope_factor);
    if n_zero_seed > 0.05 * n_z_grid
        error('mwecmass:thin_shell:DegenerateSeedGeometry', ...
              'Seed geometry grid has %d/%d (>5%%) degenerate sections.', ...
              n_zero_seed, n_z_grid);
    end
    z_grid_seed = grids_seed.z;
    A_o_seed    = grids_seed.A_outer;
    A_i_seed    = grids_seed.A_inner;
    V_jacket    = trapz(z_grid_seed, A_o_seed - A_i_seed);
    V_inner     = trapz(z_grid_seed, A_i_seed);

    % two-density seed (no root-find): the single-solid-density closed form
    % below generalises to a weighted-cumulative inversion. At rho_fill == rho_shell the general
    % form is mathematically identical to the single-density inversion but has
    % a different accumulation order in the general branch
    % (cumtrapz(c*f) accumulates rounding differently than c*cumtrapz(f)) -- short-circuited to
    % the single-density expressions in that case so the default stays bit-for-bit,
    % rho_fill >= rho_shell for C(z) to be monotone, as enforced above.
    if rho_fill == rho_shell
        V_inner_below_target = (M_buoy_target - rho_shell*V_jacket - ...
                                rho_air*V_inner) / (rho_shell - rho_air);

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
    else
        % C(z) = integral_{hull_z_min}^{z} [(rho_fill-rho_shell)*(A_o-A_i) +
        %        (rho_fill-rho_air)*A_i] dz'; solve C(z_fill) = M_buoy_target -
        %        rho_shell*V_jacket - rho_air*V_inner. Reduces to the single-density inversion (scaled)
        % when rho_fill == rho_shell -- see the M-1 derivation cited above.
        C_target = M_buoy_target - rho_shell*V_jacket - rho_air*V_inner;
        g_seed   = (rho_fill - rho_shell)*(A_o_seed - A_i_seed) + (rho_fill - rho_air)*A_i_seed;
        C_full   = trapz(z_grid_seed, g_seed);

        if C_target < 0
            z_fill_seed = hull_z_min + 1e-3;   % bounded fallback, not exact
        elseif C_target > C_full
            z_fill_seed = hull_z_max - 1e-3;   % bounded fallback, not exact
        else
            % Exact quadratic-segment inversion, not a linear interp1 of the
            % cumulative table.
            z_fill_seed = invert_piecewise_linear_cumulative(z_grid_seed, g_seed, C_target);
            z_fill_seed = max(hull_z_min + 1e-3, ...
                              min(hull_z_max - 1e-3, z_fill_seed));
        end
    end
    mwecmass.output.emit(fids, '      Warm-start z_fill seed: %.4f m (analytic, M_buoy=%.0f kg)\n', ...
            z_fill_seed, M_buoy_target);

    %% Context shared by nested objective and constraint
    ctx = struct();
    ctx.config           = config;
    ctx.rho_shell        = rho_shell;   % was ctx.rho_steel
    ctx.rho_fill         = rho_fill;    % new in the two-density model
    ctx.rho_air          = rho_air;
    ctx.max_slope_factor = max_slope_factor;
    ctx.n_z_grid         = n_z_grid;
    ctx.hull_z_min       = hull_z_min;
    ctx.hull_z_max       = hull_z_max;

    %% fmincon SQP
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

    mwecmass.output.emit(fids, '      Solver: fmincon SQP  (DVs=3 ; ceq=buoyancy ; c=GM≥gm_min)\n');
    [x_star, phi_star, exitflag, output] = fmincon( ...
        @obj_fn, x0, [], [], [], [], lb, ub, @con_fn, fminopts);

    vs_star = x_star(1);
    t_star  = x_star(2);
    zf_star = x_star(3);

    t_min_active = ((t_star - t_min) / max(t_max - t_min, eps) < 0.01);

    %% Evaluate at the optimum
    [grids_star, n_zero_interior] = mwecmass.realise.thin_shell.build_geometry_grid( ...
        config, t_star, n_z_grid, max_slope_factor);
    if n_zero_interior > 0.05 * n_z_grid
        error('mwecmass:thin_shell:TooManyDegenerateSections', ...
              'Cross-section evaluation failed at %d/%d interior z-grid points (>5%%).', ...
              n_zero_interior, n_z_grid);
    end
    realised = mwecmass.realise.thin_shell.evaluate_design_point( ...
        vs_star, t_star, zf_star, grids_star, ctx);
    if ~realised.feasible
        warning('mwecmass:thin_shell:InfeasibleOptimum', ...
                'fmincon optimum is infeasible — realised hydrostatics produced NaN/Inf.');
    end

    %% Package output
    steel_data = struct();
    steel_data.t_steel          = t_star;
    steel_data.z_fill           = zf_star;
    steel_data.draft            = realised.draft;
    steel_data.vertical_shift   = realised.vertical_shift;
    steel_data.draft_optimiser  = draft_opt;
    steel_data.vs_optimiser     = vs_opt;
    % steel_data.rho_steel stays populated (not NaN) and now carries the FILL
    % density, not the shell density -- plot_steel_solve.m's single annotated "SOLID STEEL"
    % region is drawn below z_fill (clip_polygon_below_z), which is the fill region, so this is
    % the physically correct aggregate alias. steel_data.rho_shell/.rho_fill
    % provide the explicit per-region densities.
    steel_data.rho_steel        = rho_fill;
    steel_data.rho_shell        = rho_shell;
    steel_data.rho_fill         = rho_fill;
    steel_data.rho_air          = rho_air;
    steel_data.t_max            = t_max;
    steel_data.t_min            = t_min;
    steel_data.t_min_active     = t_min_active;

    % .V_steel/.M_steel are kept as aggregate aliases (= shell+fill) so
    % plot_steel_solve.m keeps running unmodified, correct in aggregate; .V_shell/.V_fill/
    % .M_shell/.M_fill are the per-region fields of the two-density model.
    steel_data.V_steel          = realised.V_steel;
    steel_data.V_air            = realised.V_air;
    steel_data.V_hull           = realised.V_steel + realised.V_air;
    steel_data.V_shell          = realised.V_shell;
    steel_data.V_fill           = realised.V_fill;
    steel_data.M_steel          = realised.M_steel;
    steel_data.M_air            = realised.M_air;
    steel_data.M_shell          = realised.M_shell;
    steel_data.M_fill           = realised.M_fill;
    steel_data.M_total          = realised.M_total;

    steel_data.z_cg_steel       = realised.z_cg_steel;
    steel_data.z_cg_air         = realised.z_cg_air;
    % The per-region
    % centroids evaluate_design_point.m now computes (out.z_cg_fill/out.z_cg_shell) were missing from this
    % field copy -- added here, matching the existing V_shell/V_fill/M_shell/M_fill pattern
    % immediately above, so they actually reach results.steel_data instead of being silently
    % dropped at packaging.
    steel_data.z_cg_fill        = realised.z_cg_fill;
    steel_data.z_cg_shell       = realised.z_cg_shell;
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

    %% Per-strip realised equivalent density
    %  Integrate the realised material distribution over each strip
    %  defined in config.strip_edges so the cake-layer visualisations
    %  (visualize_3d_cross_section, visualize_2d_equivalent,
    %  visualize_3D_equivalent) can colour the hull by the AS-BUILT
    %  effective density instead of the optimiser's continuous rho.
    %    Below z_fill: solid steel (entire cross-section A_outer)
    %    Above z_fill: steel jacket annulus + air interior
    %  The helper inserts z_fill only for strips that contain it, preserving
    %  volume closure at strip boundaries.
    if isfield(config, 'strip_edges') && ~isempty(config.strip_edges)
        [V_env_s, V_sol_s, V_vd_s, V_fill_s, V_shell_s] = ...
            mwecmass.realise.thin_shell.strip_partition_volumes( ...
                steel_data.z_grid, steel_data.A_outer_grid, ...
                steel_data.A_inner_grid, zf_star, config.strip_edges);

        % Preserve the single-density operation order when densities match while
        % retaining the physically correct two-region mass expression.
        M_strip_s = rho_shell * V_sol_s + (rho_fill - rho_shell) * V_fill_s + rho_air * V_vd_s;   % [kg]
        rho_eff   = M_strip_s ./ max(V_env_s, eps);            % [kg/m^3]
        rho_eff(V_env_s <= 1e-12) = NaN;   % degenerate strip: no volume

        steel_data.strip_rho_eff = rho_eff;
        steel_data.strip_edges   = config.strip_edges(:);
        steel_data.strip_V_env   = V_env_s;    % [m^3]
        steel_data.strip_V_solid = V_sol_s;    % [m^3]
        steel_data.strip_V_void  = V_vd_s;     % [m^3]
        steel_data.strip_V_fill  = V_fill_s;   % [m^3] two-density model
        steel_data.strip_V_shell = V_shell_s;  % [m^3] two-density model
    else
        steel_data.strip_rho_eff = [];
        steel_data.strip_edges   = [];
        steel_data.strip_V_env   = [];
        steel_data.strip_V_solid = [];
        steel_data.strip_V_void  = [];
        steel_data.strip_V_fill  = [];
        steel_data.strip_V_shell = [];
    end
        steel_data.fill_method = 'steel_fill';

    steel_data.elapsed_seconds  = toc(t_solve_start);

    %% Console report
    mwecmass.output.emit(fids, '      ──────────────── Steel-Solve Result (fmincon SQP) ────────────────\n');
    if t_min_active
        mwecmass.output.emit(fids, '      t_steel*       : %.5f m  (%.2f in)   *** at t_min bound ***\n', ...
                t_star, t_star/0.0254);
    else
        mwecmass.output.emit(fids, '      t_steel*       : %.5f m  (%.2f in)\n', ...
                t_star, t_star/0.0254);
    end
    mwecmass.output.emit(fids, '      z_fill*        : %.4f m  (hull range [%.3f, %.3f])\n', ...
            zf_star, hull_z_min, hull_z_max);
    mwecmass.output.emit(fids, '      vs*            : %.4f m  (optimiser was %.4f m)\n', ...
            vs_star, vs_opt);
    mwecmass.output.emit(fids, '      M_total        : %.1f kg   target  : %.1f kg   (%+.3f%%)\n', ...
            realised.M_total, tgt.mass, steel_data.residuals.dmass_pct);
    mwecmass.output.emit(fids, '      Mass-balance   : %.3e kg  (%.4f%% of M_total)  ← ceq residual\n', ...
            realised.mass_balance_error_abs, steel_data.mass_balance_error_pct);
    mwecmass.output.emit(fids, '      GM realised    : %.4f m   gm_min=%.3f m   gm_target=%.3f m\n', ...
            realised.GM, config.gm_min, config.gm_target);
    mwecmass.output.emit(fids, '      T_heave        : %.3f s   target %.3f s   (%+.2f%%)\n', ...
            realised.T_heave, tgt.T_heave, steel_data.residuals.dT_heave_pct);
    mwecmass.output.emit(fids, '      T_pitch        : %.3f s   target %.3f s   (%+.2f%%)\n', ...
            realised.T_pitch, tgt.T_pitch, steel_data.residuals.dT_pitch_pct);
    mwecmass.output.emit(fids, '      Phi*           : %.6g  feasible=%d  exitflag=%d  iters=%d  elapsed=%.2f s\n', ...
            phi_star, realised.feasible, exitflag, output.iterations, ...
            steel_data.elapsed_seconds);
    mwecmass.output.emit(fids, '      ────────────────────────────────────────────────────────────────\n');

    %% --- Nested functions (close over ctx) -----------------------

    function f = obj_fn(x)
        % Objective: phi(r_heave) + phi(r_pitch)
        vs_ = x(1);  t_ = x(2);  zf_ = x(3);
        [grids_, n_zero_] = mwecmass.realise.thin_shell.build_geometry_grid( ...
            ctx.config, t_, ctx.n_z_grid, ctx.max_slope_factor);
        if n_zero_ > 0.05 * ctx.n_z_grid
            f = 1e4; return;
        end
        r_ = mwecmass.realise.thin_shell.evaluate_design_point(vs_, t_, zf_, grids_, ctx);
        if ~isfinite(r_.T_heave) || ~isfinite(r_.T_pitch)
            f = 1e4; return;
        end
        cfg_ = ctx.config;
        heave_half = 0.5 * (cfg_.T_heave_range(2) - cfg_.T_heave_range(1));
        pitch_half = 0.5 * (cfg_.T_pitch_range(2) - cfg_.T_pitch_range(1));
        r_h = (r_.T_heave - cfg_.T_heave_goal) / max(heave_half, 1e-6);
        r_p = (r_.T_pitch - cfg_.T_pitch_goal) / max(pitch_half, 1e-6);
        f = mwecmass.optim.range_penalty(r_h, cfg_.zone_k_amp) + ...
            mwecmass.optim.range_penalty(r_p, cfg_.zone_k_amp);
    end

    function [c, ceq] = con_fn(x)
        % c   = 1 − GM/gm_min ≤ 0   (stability floor)
        % ceq = M_total/(ρ_w·V_sub(-vs)) − 1 = 0   (buoyancy balance)
        vs_ = x(1);  t_ = x(2);  zf_ = x(3);
        [grids_, n_zero_] = mwecmass.realise.thin_shell.build_geometry_grid( ...
            ctx.config, t_, ctx.n_z_grid, ctx.max_slope_factor);
        if n_zero_ > 0.05 * ctx.n_z_grid
            c = 1.0;  ceq = 1.0;  return;
        end
        r_ = mwecmass.realise.thin_shell.evaluate_design_point(vs_, t_, zf_, grids_, ctx);
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


function z_star = invert_piecewise_linear_cumulative(z_grid, g, target)
%INVERT_PIECEWISE_LINEAR_CUMULATIVE Invert the cumulative trapezoidal integral
% of piecewise-linear g(z) at target. The integral is quadratic within each
% grid segment, so this solves that local quadratic rather than interpolating
% cumulative node values. Requires nonnegative g (monotone cumulative).
% z_grid and g are [N×1] vectors; z_star lies in the grid range.

    cumG = cumtrapz(z_grid, g);
    idx = find(cumG >= target, 1, 'first');
    if isempty(idx) || idx <= 1
        z_star = z_grid(1);
        return;
    end
    i = idx - 1;
    z_i  = z_grid(i);   dz  = z_grid(i+1) - z_i;
    g_i  = g(i);        g_ip1 = g(i+1);
    dC   = target - cumG(i);
    slope = (g_ip1 - g_i) / dz;

    if abs(slope) < 1e-12 * max(abs(g_i), 1)
        % g effectively constant on this segment: C(z) is linear, not quadratic.
        if g_i > 1e-12
            delta = dC / g_i;
        else
            delta = 0.5 * dz;   % zero-integrand segment: C(z) does not move; any point in it
                                 % is an equally valid inverse -- midpoint is the neutral choice.
        end
    else
        % 0.5*slope*delta^2 + g_i*delta - dC = 0; take the root with delta in [0, dz].
        disc  = max(g_i^2 + 2*slope*dC, 0);   % nonneg by construction when g >= 0 (rho_fill>=rho_shell)
        delta = (-g_i + sqrt(disc)) / slope;
    end
    delta  = max(0, min(dz, delta));
    z_star = z_i + delta;
end
