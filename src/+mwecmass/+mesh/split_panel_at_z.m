function [sub_panels, sub_verts, wl_idx] = split_panel_at_z(verts, vidx, z_wl)
%SPLIT_PANEL_AT_Z Clip one panel to the submerged side of a horizontal waterline.
% Syntax: [sub_panels,sub_verts,wl_idx] = split_panel_at_z(verts,vidx,z_wl).
% Inputs: verts [Nv x 3] m; vidx [1 x 4] panel indices; z_wl [m].
% Outputs: local clipped vertices/panels and intersection indices; triangles use v4=v3.
% Perimeter traversal preserves the original panel winding.
% See docs/METHODS_ENGINE.md#waterline-clipping

    pts   = verts(vidx, :);
    z     = pts(:, 3);
    below = z <= z_wl + 1e-10;

    edges   = [1 2; 2 3; 3 4; 4 1];
    sub_poly = zeros(0, 3);
    is_wl    = false(0, 1);

    for e = 1:4
        i1 = edges(e, 1);
        i2 = edges(e, 2);

        % Add start vertex if below waterline
        if below(i1)
            sub_poly(end+1, :) = pts(i1, :); %#ok<AGROW>
            is_wl(end+1, 1)    = false; %#ok<AGROW>
        end

        % Add intersection point if edge crosses waterline
        if below(i1) ~= below(i2)
            dz = z(i2) - z(i1);
            if abs(dz) > 1e-14
                t = (z_wl - z(i1)) / dz;
                t = max(0, min(1, t));
                new_pt    = pts(i1,:) + t * (pts(i2,:) - pts(i1,:));
                new_pt(3) = z_wl;   % snap to exact waterline
                sub_poly(end+1, :) = new_pt; %#ok<AGROW>
                is_wl(end+1, 1)    = true; %#ok<AGROW>
            end
        end
    end

    n = size(sub_poly, 1);
    if n < 3
        sub_panels = [];
        sub_verts  = [];
        wl_idx     = [];
        return;
    end

    sub_verts = sub_poly;

    if n == 3
        sub_panels = [1, 2, 3, 3];           % triangle
    elseif n == 4
        sub_panels = [1, 2, 3, 4];           % quad
    else
        % Fan triangulation from first vertex
        sub_panels = zeros(n - 2, 4);
        for k = 1:n-2
            sub_panels(k, :) = [1, k+1, k+2, k+2];
        end
    end

    wl_idx = find(is_wl);
end
