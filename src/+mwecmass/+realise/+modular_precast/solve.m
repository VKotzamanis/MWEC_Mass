function solve_data = solve(config, x_opt_3d, final_props, opts)
%SOLVE Constrained UHPC realisation solve for vertical shift, fill level, and strip jackets.
% DVs are [vertical_shift; z_fill; non-wall thicknesses]. Flotation equality and GM floor are constraints;
% heave/pitch range penalties form the objective. The wall is pinned solid and below-fill regions are solid.
% Returns steel_data-compatible fields plus per-strip thickness, solid mask, edges, and wall index.
% See docs/METHODS_ENGINE.md#realise-modular-precast
    t_solve_start = tic;
    fprintf('\n    mwecmass.realise.modular_precast.solve (fmincon SQP, constrained):\n');

    if nargin < 4 || isempty(opts), opts = struct(); end

    % ── Required wall + strip bookkeeping from caller ──────────
    if ~isfield(opts, 'wall_strip_idx') || isempty(opts.wall_strip_idx)
        error('mwecmass:modular_precast:NoWallStripIdx', ...
              'opts.wall_strip_idx must be supplied to solve.');
    end
    if ~isfield(opts, 'strip_edges') || isempty(opts.strip_edges)
        error('mwecmass:modular_precast:NoStripEdges', ...
              'opts.strip_edges must be supplied to solve.');
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
            error('mwecmass:modular_precast:MissingConfig', ...
                  'Required config field "%s" is missing or empty.', ...
                  required_cfg{k});
        end
    end

    % Resolve material + numerical options
    % config.rho_steel (7500 kg/m^3 steel, WEC_User_Input.m in.materials.thin_shell.rho_shell)
    % is a UHPC-density fallback here only because +modular_precast/solve_and_extract.m:38 always sets
    % opts.rho_steel = rho_UHPC (2500 kg/m^3) before calling this function. If that assignment
    % were ever dropped, opt_or_cfg would silently substitute the wrong material's density with
    % no warning. Assert the caller-supplied option is present; no value change.
    assert(isfield(opts, 'rho_steel') && ~isempty(opts.rho_steel), ...
           'mwecmass:modular_precast:MissingUHPCDensity', ...
           ['modular_precast.solve requires opts.rho_steel (the UHPC/hull ', ...
            'density) to be set by the caller; falling back to config.rho_steel would ', ...
            'silently substitute the thin-shell steel density (%.0f kg/m^3) for UHPC.'], ...
           config.rho_steel);
    rho_uhpc         = mwecmass.internal.option_or_config(opts, 'rho_steel',        config.rho_steel);
    rho_air          = mwecmass.internal.option_or_config(opts, 'rho_air',          config.rho_air);
    t_init           = mwecmass.internal.option_or_config(opts, 't_init',           config.steel_t_init);
    t_min            = mwecmass.internal.option_or_config(opts, 't_min',            config.steel_t_min);
    max_slope_factor = mwecmass.internal.option_or_config(opts, 'max_slope_factor', config.steel_max_slope_factor);
    n_z_grid         = mwecmass.internal.option_or_config(opts, 'n_z_grid',         config.steel_n_z_grid);

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
        r_probes(k) = mwecmass.geometry.compute_rmin_at_z( ...
                          config.ms2_model, z_probes(k), 60, cache);
    end
    r_min_global = min(r_probes(r_probes > 0));
    if isempty(r_min_global) || r_min_global <= 0
        error('mwecmass:modular_precast:NoFiniteRadius', ...
              'compute_rmin_at_z returned no positive radius.');
    end
    t_max = 0.5 * r_min_global;
    if t_min >= t_max
        error('mwecmass:modular_precast:TMinExceedsTMax', ...
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

    %% Analytic z_fill seed
    is_solid_seed = false(N_strips, 1);
    is_solid_seed(wall_strip_idx) = true;
    t_strip_seed = t_init * ones(N_strips, 1);
    t_strip_seed(wall_strip_idx) = inf;

    [grids_seed, n_zero_seed] = mwecmass.realise.modular_precast.build_geometry_grid( ...
        config, strip_edges, t_strip_seed, is_solid_seed, n_z_grid, max_slope_factor);
    if n_zero_seed > 0.05 * n_z_grid
        error('mwecmass:modular_precast:DegenerateSeedGeometry', ...
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

    %% Nested-solver context
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
    % is_solid_strip held constant during the solve: wall pinned solid;
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

    %% DESIGN-TRAJECTORY LOGGING
    % Record accepted solver iterates and realised scalars for diagnostics.
    iter_hist = struct( ...
        'iter',        [], ...   % fmincon iteration index
        'vs',          [], ...   % vertical shift            [m]
        'z_fill',      [], ...   % UHPC/void transition      [m]
        'draft',       [], ...   % |hull_z_min + vs|         [m]
        't_strip',     [], ...   % N_strips x n_iter, Inf = solid [m]
        'M_total',     [], ...   % realised mass             [kg]
        'M_buoy',      [], ...   % rho_w * V_sub(-vs)        [kg]
        'GM',          [], ...   % metacentric height        [m]
        'T_heave',     [], ...   % uncoupled heave period    [s]
        'T_pitch',     [], ...   % uncoupled pitch period    [s]
        'V_steel',     [], ...   % UHPC volume               [m^3]
        'V_air',       [], ...   % void volume               [m^3]
        'CG_z_world',  [], ...   % CG elevation, world frame [m]
        'phi',         [], ...   % objective value           [-]
        'ceq',         [], ...   % M_total/M_buoy - 1        [-]
        'c_gm',        []);      % 1 - GM/gm_min             [-]

    fminopts = optimoptions('fmincon', ...
        'Algorithm',              'sqp', ...
        'Display',                'iter', ...
        'StepTolerance',          1e-8, ...
        'OptimalityTolerance',    1e-6, ...
        'ConstraintTolerance',    1e-6, ...
        'MaxIterations',          300, ...
        'MaxFunctionEvaluations', 2000, ...
        'OutputFcn',              @out_fn);

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
    [grids_final, n_zero_interior] = mwecmass.realise.modular_precast.build_geometry_grid( ...
        config, strip_edges, t_offset_strip, is_solid_final, n_z_grid, max_slope_factor);
    if n_zero_interior > 0.05 * n_z_grid
        error('mwecmass:modular_precast:TooManyDegenerateSections', ...
              'Final strip-aware grid has %d/%d (>5%%) degenerate sections.', ...
              n_zero_interior, n_z_grid);
    end
    realised = mwecmass.realise.modular_precast.evaluate_design_point( ...
        vs_star, t_offset_strip, is_solid_final, zf_star, grids_final, ctx);

    if ~realised.feasible
        warning('mwecmass:modular_precast:InfeasibleUHPCOptimum', ...
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

    % Design trajectory (one record per accepted SQP iterate).
    % NOTE: t_strip here is the DV state during the solve, where only
    % the wall is Inf.  The post-solve promotion of below-z_fill
    % strips to solid is NOT applied retroactively to the history —
    % t_offset_strip above is the final, promoted vector.
    solve_data.iter_history = iter_hist;

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
        [g_, n_zero_] = mwecmass.realise.modular_precast.build_geometry_grid( ...
            ctx.config, ctx.strip_edges, t_off, is_sol, ...
            ctx.n_z_grid, ctx.max_slope_factor);
        if n_zero_ > 0.05 * ctx.n_z_grid
            f = 1e4; return;
        end
        r_ = mwecmass.realise.modular_precast.evaluate_design_point( ...
            vs_, t_off, is_sol, zf_, g_, ctx);
        if ~isfinite(r_.T_heave) || ~isfinite(r_.T_pitch)
            f = 1e4; return;
        end
        cfg_ = ctx.config;
        h_half = 0.5 * (cfg_.T_heave_range(2) - cfg_.T_heave_range(1));
        p_half = 0.5 * (cfg_.T_pitch_range(2) - cfg_.T_pitch_range(1));
        r_h = (r_.T_heave - cfg_.T_heave_goal) / max(h_half, 1e-6);
        r_p = (r_.T_pitch - cfg_.T_pitch_goal) / max(p_half, 1e-6);
        f = mwecmass.optim.range_penalty(r_h, cfg_.zone_k_amp) + ...
            mwecmass.optim.range_penalty(r_p, cfg_.zone_k_amp);
    end

    function [c, ceq] = con_fn(x)
        vs_ = x(1);  zf_ = x(2);
        [t_off, is_sol] = unpack_dvs(x);
        [g_, n_zero_] = mwecmass.realise.modular_precast.build_geometry_grid( ...
            ctx.config, ctx.strip_edges, t_off, is_sol, ...
            ctx.n_z_grid, ctx.max_slope_factor);
        if n_zero_ > 0.05 * ctx.n_z_grid
            c = 1.0; ceq = 1.0; return;
        end
        r_ = mwecmass.realise.modular_precast.evaluate_design_point( ...
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

    function stop = out_fn(x, optimValues, state)
        % OUT_FN  fmincon OutputFcn — append one record per accepted
        %   iterate to iter_hist (nested-scope variable).  Never
        %   alters the solve: always returns stop = false.
        %
        %   'init' is logged as iteration 0 so the warm start appears
        %   as the first record.  A degenerate geometry grid is
        %   skipped rather than logged, so the history never contains
        %   points the forward model could not evaluate.
        stop = false;
        if ~(strcmp(state,'init') || strcmp(state,'iter') || strcmp(state,'done'))
            return;
        end

        vs_l = x(1);
        zf_l = x(2);
        [t_l, is_l] = unpack_dvs(x);

        [g_l, nz_l] = mwecmass.realise.modular_precast.build_geometry_grid( ...
            ctx.config, ctx.strip_edges, t_l, is_l, ...
            ctx.n_z_grid, ctx.max_slope_factor);
        if nz_l > 0.05 * ctx.n_z_grid
            return;
        end

        r_l = mwecmass.realise.modular_precast.evaluate_design_point( ...
            vs_l, t_l, is_l, zf_l, g_l, ctx);

        iter_hist.iter(end+1)       = optimValues.iteration;
        iter_hist.vs(end+1)         = vs_l;
        iter_hist.z_fill(end+1)     = zf_l;
        iter_hist.draft(end+1)      = r_l.draft;
        iter_hist.t_strip           = [iter_hist.t_strip, t_l(:)];
        iter_hist.M_total(end+1)    = r_l.M_total;
        iter_hist.M_buoy(end+1)     = r_l.mass_buoyant_force;
        iter_hist.GM(end+1)         = r_l.GM;
        iter_hist.T_heave(end+1)    = r_l.T_heave;
        iter_hist.T_pitch(end+1)    = r_l.T_pitch;
        iter_hist.V_steel(end+1)    = r_l.V_steel;
        iter_hist.V_air(end+1)      = r_l.V_air;
        iter_hist.CG_z_world(end+1) = r_l.CG_z_world;
        iter_hist.phi(end+1)        = optimValues.fval;

        if isfinite(r_l.M_total) && isfinite(r_l.mass_buoyant_force) && ...
                r_l.mass_buoyant_force > 0
            iter_hist.ceq(end+1) = r_l.M_total / r_l.mass_buoyant_force - 1.0;
        else
            iter_hist.ceq(end+1) = NaN;
        end
        if isfinite(r_l.GM)
            iter_hist.c_gm(end+1) = 1.0 - r_l.GM / ctx.config.gm_min;
        else
            iter_hist.c_gm(end+1) = NaN;
        end
    end

end
