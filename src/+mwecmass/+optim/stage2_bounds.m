function [lb, ub] = stage2_bounds(config)
%STAGE2_BOUNDS Bounds for x=[vertical_shift, rho_1..rho_N].
% Uses vertical_shift_bounds and ballast_density_bounds; the density floors of the shell
% realisations (config.per_strip_density_lb, both modes) raise the per-strip lower bounds, and the
% modular-precast wall strip is pinned to constructability_rho_hull.
    N = config.num_ballast_sections;

    lb = [config.vertical_shift_bounds(1), ...
          ones(1, N) * config.ballast_density_bounds(1)];
    ub = [config.vertical_shift_bounds(2), ...
          ones(1, N) * config.ballast_density_bounds(2)];

    if ~isempty(config.per_strip_density_lb)
        for i = 1:N
            lb(1 + i) = max(lb(1 + i), config.per_strip_density_lb(i));
        end
    end
    if ~isempty(config.wall_strip_index)
        w_idx = config.wall_strip_index;
        lb(1 + w_idx) = config.constructability_rho_hull;
        ub(1 + w_idx) = config.constructability_rho_hull;
    end
end
