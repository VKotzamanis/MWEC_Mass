function [c, ceq] = stage2_constraints(x, config)
%STAGE2_CONSTRAINTS Nonlinear constraints for Stage-2 fmincon.
% x=[vertical_shift, rho_1..rho_N]. c<=0 enforces GM >= gm_min and the adjacent density ratio;
% ceq enforces mass/buoyant-force-1=0. A pinned wall node is excluded from adjacent pair checks so
% its material density does not constrain neighboring platform nodes.
% See docs/METHODS_ENGINE.md#optim-stage2-formulation
try
    props     = mwecmass.hydrostatics.properties_3d(x, config);
    densities = x(2:end);
    N         = length(densities);

    % props.GM_L is unbiased (no k_gm correction); config.gm_min is the same floor Stage-1's 2D
    % constraint applies to the k_gm-biased props.GM (mwecmass.optim.solve_2d_surrogate.m). Documentation-only
    % note, no value change.
    c_gm = 1.0 - props.GM_L / config.gm_min;

    % Identify which adjacent pairs to constrain.
    % Skip any pair bridging the wall-platform boundary.
    w_idx = config.wall_strip_index;   % [] if not constructability
    constrained_pairs = [];
    for i = 1:(N-1)
        if ~isempty(w_idx) && (i == w_idx || i + 1 == w_idx)
            continue;   % skip wall-platform boundary
        end
        constrained_pairs(end+1) = i; %#ok<AGROW>
    end
    n_pairs = length(constrained_pairs);

    c_ratio = zeros(n_pairs, 1);
    for j = 1:n_pairs
        i = constrained_pairs(j);
        c_ratio(j) = densities(i) / (densities(i+1) + 1) - config.max_density_ratio;
    end

    c   = [c_gm; c_ratio];
    ceq = props.mass_total / props.mass_buoyant_force - 1.0;
catch
    % Return O(1) violations for diverging solutions.
    % Size must be consistent — use safe fallback.
    n_pairs_fallback = max(0, length(x) - 2);
    if config.enable_constructability && ~isempty(config.wall_strip_index)
        n_pairs_fallback = max(0, n_pairs_fallback - 1);
    end
    c   = ones(1 + n_pairs_fallback, 1);
    ceq = 1.0;
end
end
