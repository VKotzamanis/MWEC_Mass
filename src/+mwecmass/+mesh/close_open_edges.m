function [verts, panels, surf_ids, is_cap, n_cap] = ...
        close_open_edges(verts, panels, surf_ids, is_cap)
%CLOSE_OPEN_EDGES Fill horizontal open mesh boundaries with fan-triangulated caps.
% Syntax: [verts,panels,surf_ids,is_cap,n_cap] = close_open_edges(verts,panels,surf_ids,is_cap).
% Inputs: verts [Nv x 3] m; panels [Np x 4] indices (triangle uses v3=v4); surf_ids/is_cap [Np x 1].
% Outputs append cap vertices/panels and report n_cap; boundaries within 0.01 m in z are capped.

    n_cap = 0;

    edge_map = containers.Map('KeyType', 'char', 'ValueType', 'int32');
    for p = 1:size(panels, 1)
        v = panels(p, :);
        if v(3) == v(4)
            ee = [v(1) v(2); v(2) v(3); v(3) v(1)];
        else
            ee = [v(1) v(2); v(2) v(3); v(3) v(4); v(4) v(1)];
        end
        for e = 1:size(ee, 1)
            key = sprintf('%d_%d', min(ee(e,:)), max(ee(e,:)));
            if edge_map.isKey(key)
                edge_map(key) = edge_map(key) + 1;
            else
                edge_map(key) = 1;
            end
        end
    end

    boundary_edges = [];
    keys = edge_map.keys();
    for i = 1:length(keys)
        if edge_map(keys{i}) == 1
            v = sscanf(keys{i}, '%d_%d');
            boundary_edges = [boundary_edges; v']; %#ok<AGROW>
        end
    end

    if isempty(boundary_edges), return; end

    boundary_vidx = unique(boundary_edges(:));
    boundary_z    = verts(boundary_vidx, 3);
    z_tol         = 0.01;  % [m]
    z_levels      = unique(round(boundary_z / z_tol) * z_tol);

    for iz = 1:length(z_levels)
        z_lvl    = z_levels(iz);
        at_level = boundary_vidx(abs(boundary_z - z_lvl) < z_tol);

        if length(at_level) < 3, continue; end

        cap_center = mean(verts(at_level, :), 1);
        center_idx = size(verts, 1) + 1;
        verts      = [verts; cap_center]; %#ok<AGROW>

        rel    = verts(at_level, :) - cap_center;
        angles = atan2(rel(:, 2), rel(:, 1));
        [~, order]    = sort(angles);
        ordered_verts = at_level(order);

        n_fan = length(ordered_verts);
        for k = 1:n_fan
            k_next = mod(k, n_fan) + 1;
            new_tri = [center_idx, ordered_verts(k), ...
                       ordered_verts(k_next), ordered_verts(k_next)];
            panels   = [panels; new_tri]; %#ok<AGROW>
            surf_ids = [surf_ids; 0]; %#ok<AGROW>
            is_cap   = [is_cap; true]; %#ok<AGROW>
            n_cap    = n_cap + 1;
        end
    end
end
