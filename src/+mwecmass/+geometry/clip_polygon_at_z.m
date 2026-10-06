function poly_out = clip_polygon_at_z(poly, z_plane, keep_side)
%CLIP_POLYGON_AT_Z Clip an [x,z] polygon against a horizontal plane.
%   z_plane is in m; keep_side selects below or above with a 1e-9 m tolerance.

    try
        if isempty(poly) || size(poly, 1) < 2
            poly_out = [];
            return;
        end

        n = size(poly, 1);
        output = zeros(2*n, 2);
        out_count = 0;
        epsilon = 1e-9;

        for i = 1:n
            p1 = poly(i, :);
            next_idx = mod(i, n) + 1;
            p2 = poly(next_idx, :);

            if strcmp(keep_side, 'below')
                p1_in = p1(2) <= z_plane + epsilon;
                p2_in = p2(2) <= z_plane + epsilon;
            else
                p1_in = p1(2) >= z_plane - epsilon;
                p2_in = p2(2) >= z_plane - epsilon;
            end

            if p1_in
                out_count = out_count + 1;
                output(out_count, :) = p1;
            end

            if xor(p1_in, p2_in)
                dz = p2(2) - p1(2);
                if abs(dz) > epsilon
                    t = (z_plane - p1(2)) / dz;
                    t = max(0, min(1, t));
                    intersection_point = p1 + t * (p2 - p1);
                    out_count = out_count + 1;
                    output(out_count, :) = intersection_point;
                end
            end
        end

        if out_count > 2
            poly_out_temp = output(1:out_count, :);
            [~, idx] = unique(round(poly_out_temp*1e6)/1e6, 'rows', 'stable');
            poly_out = poly_out_temp(sort(idx), :);
        else
            poly_out = [];
        end

    catch
        poly_out = [];
    end
end
