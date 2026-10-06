function [x_opt, props_2d, conv_data, hist] = ...
        stage1_sweep(x0, config, hist)
%STAGE1_SWEEP Rank a draft grid, then refine the best candidates with constrained fmincon.
% Tier 1 uses one mass-balanced probe per draft; Tier 2 applies the full Stage-2 objective and constraints.
% A missing or singleton hydro_drafts grid falls back to x0. Outputs include the selected design,
% properties, convergence data, and history.
% See docs/METHODS_ENGINE.md#optim-stage1-two-tier
    % ── Guard: fall back to passthrough if no draft grid ──────────
    if ~isfield(config, 'hydro_drafts') || length(config.hydro_drafts) <= 1
        fprintf(' Skipping Stage 1 (no draft grid). Using default x0.\n');
        x_opt    = x0;
        props_2d = mwecmass.hydrostatics.properties_2d(x0, config);
        conv_data = struct('final_objective', NaN, 'final_exitflag', 0);
        return;
    end

    vs_grid = config.hydro_drafts(:)';
    N       = length(vs_grid);

    fprintf(' Stage 1: two-tier draft-landscape sweep (%d drafts)\n', N);
    fprintf('   vs range: [%.3f, %.3f] m\n', vs_grid(1), vs_grid(end));

    % Stage-2 bounds reused by each refined candidate
    [lb_base, ub_base] = mwecmass.optim.stage2_bounds(config);

    % Sweep results
    sweep.vs       = vs_grid;
    sweep.fval     = inf(1, N);
    sweep.exitflag = zeros(1, N);
    sweep.x        = cell(1, N);
    sweep.props    = cell(1, N);
    sweep.feasible = false(1, N);

    obj_fun = @(x) mwecmass.optim.stage2_objective(x, config);
    con_fun = @(x) mwecmass.optim.stage2_constraints(x, config);

    % ── TIER 1: physics-based single-point screen ─────────────────
    fprintf('\n Tier 1: physics-based screen...\n');

    for k = 1:N
        vs_k   = vs_grid(k);
        screen = mwecmass.optim.stage1_screen_draft(vs_k, config);

        sweep.fval(k)     = screen.fval;
        sweep.exitflag(k) = 1;
        sweep.x{k}        = screen.x;
        sweep.props{k}    = screen.props;
        sweep.feasible(k) = screen.feasible;

        fprintf('   [%2d/%d] vs=%+.3f  f=%.4f  GM=%.3f  T_h=%.2f  T_p=%.2f  %s\n', ...
            k, N, vs_k, screen.fval, screen.props.GM_L, ...
            screen.props.periods.heave, screen.props.periods.pitch, ...
            mwecmass.internal.ternary(screen.feasible, 'FEAS', 'infeas'));
    end

    % ── TIER 2: targeted fmincon on top-K candidates ──────────────
    % Rank by heave-period proximity; Tier 2 still enforces GM and all density constraints.
    K_refine   = config.n_sweep_refine;
    T_h_vals   = arrayfun(@(k) sweep.props{k}.periods.heave, 1:N);
    T_h_vals(~isfinite(T_h_vals)) = Inf;   % dry / degenerate drafts rank last
    [~, rank]  = sort(abs(T_h_vals - config.T_heave_goal));
    top_idx    = rank(1:min(K_refine, N));

    fprintf('\n Tier 1 → Tier 2 candidates ranked by |T_heave − %.1f s|:\n', ...
            config.T_heave_goal);
    for ki = 1:length(top_idx)
        k_ki = top_idx(ki);
        fprintf('   [%d] vs=%+.3f m  T_h=%.2f s  |err|=%.2f s  GM=%.3f m  %s\n', ...
                ki, vs_grid(k_ki), T_h_vals(k_ki), ...
                abs(T_h_vals(k_ki) - config.T_heave_goal), ...
                sweep.props{k_ki}.GM_L, ...
                mwecmass.internal.ternary(sweep.feasible(k_ki), 'GM-feas', 'GM-infeas'));
    end

    % Lightweight fmincon — looser than Stage 2 (warm-start only)
    opts_refine = optimoptions('fmincon', ...
        'Algorithm',              'sqp', ...
        'Display',                'none', ...
        'MaxFunctionEvaluations', 150, ...
        'MaxIterations',          25, ...
        'ConstraintTolerance',    1e-4, ...
        'OptimalityTolerance',    1e-4, ...
        'StepTolerance',          1e-6, ...
        'ScaleProblem',           true);

    fprintf('\n Tier 2: refining %d candidates...\n', length(top_idx));

    for j = 1:length(top_idx)
        k    = top_idx(j);
        vs_k = vs_grid(k);

        lb_k = lb_base;  lb_k(1) = vs_k;
        ub_k = ub_base;  ub_k(1) = vs_k;
        x0_k = max(lb_k, min(ub_k, sweep.x{k}));

        try
            [x_k, fval_k, ef_k] = fmincon(obj_fun, x0_k, ...
                [], [], [], [], lb_k, ub_k, con_fun, opts_refine);

            props_k  = mwecmass.hydrostatics.properties_3d(x_k, config);
            mass_err = abs(props_k.mass_total - props_k.mass_buoyant_force) ...
                       / max(props_k.mass_total, 1);
            is_feas  = ef_k > 0 && mass_err < 0.01 && props_k.GM_L >= config.gm_min;

            % Only update if Tier 2 improved on Tier 1
            if fval_k < sweep.fval(k)
                sweep.fval(k)     = fval_k;
                sweep.exitflag(k) = ef_k;
                sweep.x{k}        = x_k;
                sweep.props{k}    = props_k;
                sweep.feasible(k) = is_feas;
            end

            % Print stored candidate, not raw Tier-2; Tier 2 is only used if it improves Tier 1.
            fprintf('   [T2 %d/%d] vs=%+.3f  f=%.4f  GM=%.3f  %s\n', ...
                j, length(top_idx), vs_k, sweep.fval(k), sweep.props{k}.GM_L, ...
                mwecmass.internal.ternary(sweep.feasible(k), 'FEAS', 'infeas'));

        catch ME_t2
            fprintf('   [T2 %d/%d] vs=%+.3f  FAILED: %s\n', ...
                j, length(top_idx), vs_k, ME_t2.message);
        end
    end

    % ── Pick best point for Stage-2 warm-start ────────────────────
    feas_idx_final = find(sweep.feasible);
    if ~isempty(feas_idx_final)
        [~, best_in_feas] = min(sweep.fval(feas_idx_final));
        best = feas_idx_final(best_in_feas);
        fprintf('\n   Best feasible: draft %d (vs=%+.3f m, f=%.4f)\n', ...
                best, vs_grid(best), sweep.fval(best));
    else
        [~, best] = min(sweep.fval);
        warning('WEC:NoFeasibleDraft', ...
                'No feasible draft after Tier-2 refinement. Using lowest objective (vs=%+.3f m).', ...
                vs_grid(best));
    end

    x_opt    = sweep.x{best};
    props_2d = mwecmass.hydrostatics.properties_2d(x_opt, config);
    conv_data = struct('final_objective',  sweep.fval(best), ...
                       'final_exitflag',   sweep.exitflag(best), ...
                       'sweep',            sweep);

    fprintf('   Warm-start for Stage 2: vs=%+.3f m, rho=[', x_opt(1));
    fprintf('%.0f ', x_opt(2:end));
    fprintf(']\n');
end
