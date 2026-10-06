function draw_composite_strips(ax, outer_profile, inner_profile, ...
                               node_z, densities, ~, style)
%DRAW_COMPOSITE_STRIPS Draw colour-filled density strips between an outer and an inner (shell-offset) profile.
% node_z: vertical position of density nodes (m, z positive upward).
% densities: one scalar per node, used as patch face colors.
    if nargin < 7 || isempty(style)
        style = mwecmass.output.figures.presentation_style();
    end
    N      = length(densities);
    node_z = node_z(:);
    dz_half = (node_z(2) - node_z(1)) / 2;
    bounds  = [ node_z(1) - dz_half; ...
               (node_z(1:end-1) + node_z(2:end)) / 2; ...
                node_z(end) + dz_half ];
    outer_slabs = cell(N, 1);
    inner_slabs = cell(N, 1);
    for k = 1:N
        outer_slabs{k} = mwecmass.output.figures.clip_profile_to_z_range( ...
            outer_profile, bounds(k), bounds(k+1));
        inner_slabs{k} = mwecmass.output.figures.clip_profile_to_z_range( ...
            inner_profile, bounds(k), bounds(k+1));
    end
    for k = 1:N
        outer_slab = outer_slabs{k};
        if isempty(outer_slab) || size(outer_slab,1) < 3, continue; end
        patch(ax, outer_slab(:,1), outer_slab(:,2), style.fill_palette.shell, ...
              'EdgeColor', 'none', 'FaceAlpha', 1.0);
    end
    for k = 1:N
        inner_slab = inner_slabs{k};
        if isempty(inner_slab) || size(inner_slab,1) < 3, continue; end
        patch(ax, inner_slab(:,1), inner_slab(:,2), densities(k), ...
              'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 1.0);
    end
    for k = 1:N
        z_lo = bounds(k);
        z_hi = bounds(k+1);
        outer_slab = outer_slabs{k};
        if isempty(outer_slab) || size(outer_slab,1) < 3, continue; end
        for bz = [z_lo, z_hi]
            x_cross = [];
            n_out = size(outer_profile,1);
            for j = 1:n_out
                j2 = mod(j, n_out) + 1;
                z1 = outer_profile(j,2);  z2 = outer_profile(j2,2);
                if (z1 - bz) * (z2 - bz) <= 0 && abs(z2-z1) > 1e-12
                    t_edge = (bz - z1) / (z2 - z1);
                    x_cross(end+1) = outer_profile(j,1) + ...
                        t_edge*(outer_profile(j2,1) - outer_profile(j,1)); %#ok<AGROW>
                end
            end
            if length(x_cross) >= 2
                line(ax, [min(x_cross), max(x_cross)], [bz, bz], ...
                     'Color', style.fill_palette.boundary, 'LineWidth', 0.8);
            end
        end
    end
end
