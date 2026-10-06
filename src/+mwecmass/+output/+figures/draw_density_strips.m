function draw_density_strips(ax, profile, node_z, densities, ~, ~, strip_bounds, style)
%DRAW_DENSITY_STRIPS Draw N horizontal slab polygons, one flat colour per density layer, on the given axes.
    if nargin < 8 || isempty(style)
        style = mwecmass.output.figures.presentation_style();
    end
    N      = length(densities);
    node_z = node_z(:);
    if nargin >= 7 && ~isempty(strip_bounds) && length(strip_bounds) == N + 1
        bounds = strip_bounds(:);
    elseif N > 1
        dz_half  = (node_z(2) - node_z(1)) / 2;
        bounds   = [ node_z(1) - dz_half; ...
                    (node_z(1:end-1) + node_z(2:end)) / 2; ...
                     node_z(end) + dz_half ];
    else
        % At least ±0.5 m half-height when profile z-range is small.
        z_rng  = max(abs(node_z), 0.5);
        bounds = [ node_z(1) - z_rng; node_z(1) + z_rng ];
    end

    for k = 1:N
        z_lo      = bounds(k);
        z_hi      = bounds(k+1);
        rho_strip = densities(k);

        slab = mwecmass.output.figures.clip_profile_to_z_range(profile, z_lo, z_hi);
        if isempty(slab) || size(slab,1) < 3
            continue;
        end

        patch(ax, slab(:,1), slab(:,2), rho_strip, ...
              'EdgeColor', 'none', 'FaceAlpha', 0.85);
    end
    for k = 1:length(bounds)
        x_cross = [];
        for j = 1:size(profile,1)
            j2 = mod(j, size(profile,1)) + 1;
            z1 = profile(j,2);  z2 = profile(j2,2);
            if (z1 - bounds(k)) * (z2 - bounds(k)) <= 0 && abs(z2-z1) > 1e-12
                t = (bounds(k) - z1) / (z2 - z1);
                x_cross(end+1) = profile(j,1) + t*(profile(j2,1) - profile(j,1)); %#ok<AGROW>
            end
        end
        if length(x_cross) >= 2
            plot(ax, [min(x_cross), max(x_cross)], [bounds(k), bounds(k)], ...
                 '-', 'Color', style.fill_palette.boundary, 'LineWidth', 0.8);
        end
    end
end
