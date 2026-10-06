function intersections = find_waterline_intersections(profile, z_level)
%FIND_WATERLINE_INTERSECTIONS Intersect a closed profile with a horizontal waterline.
%   profile and returned intersections are [x,z] coordinates in m.

    try
        intersections = [];

        if isempty(profile) || size(profile, 1) < 2
            return;
        end

        if ~isscalar(z_level) || ~isnumeric(z_level)
            return;
        end

        n_vertices = size(profile, 1);
        epsilon = 1e-10;

        for i = 1:n_vertices
            p1 = profile(i, :);
            p2_idx = mod(i, n_vertices) + 1;
            p2 = profile(p2_idx, :);

            z1 = p1(2);
            z2 = p2(2);

            if (z1 - z_level) * (z2 - z_level) < -epsilon
                dz = z2 - z1;
                if abs(dz) > epsilon
                    t = (z_level - z1) / dz;
                    t = max(0, min(1, t));
                    x_intersect = p1(1) + t * (p2(1) - p1(1));
                    intersections = [intersections; x_intersect, z_level]; %#ok<AGROW>
                end
            end
        end

    catch
        intersections = [];
    end
end
