function [c, ceq] = stage2_constraints(x, config)
%STAGE2_CONSTRAINTS Nonlinear constraints for Stage-2 fmincon.
% x=[vertical_shift, rho_1..rho_N]. c<=0 enforces GM, adjacent density ratio, and heavy-bottom
% monotonicity; ceq enforces mass/buoyant-force-1=0. A pinned wall node is excluded from adjacent
% pair checks so its material density does not constrain neighboring platform nodes.
% See docs/METHODS_ENGINE.md#optim-stage2-formulation
try
    props     = mwecmass.hydrostatics.properties_3d(x, config);
    densities = x(2:end);
    % When shell is enabled, densities are CORE densities (not bulk).
    % The ratio and monotonic constraints apply to the core only,
    % because the shell density is uniform across all strips and
    % poses no manufacturability concern.
    N         = length(densities);
    rho_max   = config.ballast_density_bounds(2);

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
    c_mono  = zeros(n_pairs, 1);
    for j = 1:n_pairs
        i = constrained_pairs(j);
        c_ratio(j) = densities(i) / (densities(i+1) + 1) - config.max_density_ratio;
        c_mono(j)  = (densities(i+1) - densities(i)) / rho_max;
    end

    % Minimum mass constraint (constructability realisation type).
    % Prevents the optimizer from requesting density distributions
    % that produce less mass than is physically achievable.
    if config.m_min_constructability > 0
        c_mass_min = 1.0 - props.mass_total / config.m_min_constructability;
    else
        c_mass_min = [];
    end

    c   = [c_gm; c_ratio; c_mono; c_mass_min];
    ceq = props.mass_total / props.mass_buoyant_force - 1.0;
catch
    % Return O(1) violations for diverging solutions.
    % Size must be consistent — use safe fallback.
    n_pairs_fallback = max(0, length(x) - 2);
    if config.enable_constructability && ~isempty(config.wall_strip_index)
        n_pairs_fallback = max(0, n_pairs_fallback - 1);
    end
    n_extra = 0;
    if isfield(config, 'm_min_constructability') && config.m_min_constructability > 0
        n_extra = 1;
    end
    c   = ones(1 + 2 * n_pairs_fallback + n_extra, 1);
    ceq = 1.0;
end
end
