function normals = compute_panel_normals(verts, panels)
%COMPUTE_PANEL_NORMALS Normalize panel normals using HAMS CalTransNormals winding rules.
% Syntax: normals = mwecmass.mesh.compute_panel_normals(verts,panels).
% Inputs: verts [Nv x 3] m; panels [Np x 4] indices (v3=v4 denotes a triangle).
% Output: normals [Np x 3], unit vectors; degenerate panels receive [0 0 1].
% Quad uses cross(V3-V1,V4-V2); triangle uses cross(V2-V1,V3-V2).
% See docs/METHODS_ENGINE.md#panel-normals

    n_p = size(panels, 1);
    normals = zeros(n_p, 3);

    for p = 1:n_p
        v = panels(p,:);
        if v(3) == v(4)
            % Triangle: cross(V2-V1, V3-V2)
            d1 = verts(v(2),:) - verts(v(1),:);
            d2 = verts(v(3),:) - verts(v(2),:);
        else
            % Quad: cross(V3-V1, V4-V2)
            d1 = verts(v(3),:) - verts(v(1),:);
            d2 = verts(v(4),:) - verts(v(2),:);
        end
        n_vec = cross(d1, d2);
        n_len = norm(n_vec);

        if n_len > 1e-14
            normals(p,:) = n_vec / n_len;
        else
            normals(p,:) = [0 0 1];
        end
    end
end
