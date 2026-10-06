function out = clip_to_halfplane(poly, z_bound, sign_dir)
%CLIP_TO_HALFPLANE Sutherland-Hodgman clip of a polygon against one half-plane bound.
% poly is n-by-2 [x, z] vertices; z_bound is the z-coordinate of the half-plane
% boundary; sign_dir is +1 to keep z >= z_bound, -1 to keep z <= z_bound.
    if isempty(poly) || size(poly, 1) < 2
        out = [];
        return;
    end
    n   = size(poly, 1);
    out = zeros(2*n, 2);
    cnt = 0;
    for i = 1:n
        j  = mod(i, n) + 1;
        zi = poly(i, 2);
        zj = poly(j, 2);
        inside_i = sign_dir * (zi - z_bound) >= -1e-12;
        inside_j = sign_dir * (zj - z_bound) >= -1e-12;
        if inside_i && inside_j
            cnt = cnt + 1;
            out(cnt, :) = poly(j, :);
        elseif inside_i && ~inside_j
            t = (z_bound - zi) / (zj - zi);
            cnt = cnt + 1;
            out(cnt, :) = poly(i, :) + t * (poly(j, :) - poly(i, :));
        elseif ~inside_i && inside_j
            t = (z_bound - zi) / (zj - zi);
            cnt = cnt + 1;
            out(cnt, :) = poly(i, :) + t * (poly(j, :) - poly(i, :));
            cnt = cnt + 1;
            out(cnt, :) = poly(j, :);
        end
    end
    out = out(1:cnt, :);
end
