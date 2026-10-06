function [geom, iner, cpmo] = polygon_properties(x, y)
%POLYGON_PROPERTIES Compute signed area, centroid, and planar moments.
% x,y are polygon vertex vectors; an unclosed polygon is closed internally.
% Outputs geom=[A xc yc 0], iner=[Ixx Iyy Ixy 0 0 0] about the origin,
% and cpmo=[Iuu Ivv Iuv 0 0 0] about the centroid. Degenerate input returns zeros.

    try
        if nargin < 2
            error('polygon_properties requires x and y coordinates');
        end

        x = x(:);
        y = y(:);

        % Close polygon if not already closed
        if (x(1) ~= x(end)) || (y(1) ~= y(end))
            x(end+1) = x(1);
            y(end+1) = y(1);
        end

        n = length(x);

        % Require at least 3 unique vertices
        if n < 4
            geom = [0, 0, 0, 0];
            iner = [0, 0, 0, 0, 0, 0];
            cpmo = [0, 0, 0, 0, 0, 0];
            return;
        end

        xi = x(1:n-1);
        yi = y(1:n-1);
        xip1 = x(2:n);
        yip1 = y(2:n);

        % Shoelace formula
        a = xi .* yip1 - xip1 .* yi;
        Area = 0.5 * sum(a);

        if abs(Area) < 1e-12
            xc = mean(xi);
            yc = mean(yi);
            Ixx = 0; Iyy = 0; Ixy = 0;
            Iuu = 0; Ivv = 0; Iuv = 0;
        else
            xc = sum((xi + xip1) .* a) / (6 * Area);
            yc = sum((yi + yip1) .* a) / (6 * Area);

            Ixx = sum((yi.^2 + yi.*yip1 + yip1.^2) .* a) / 12;
            Iyy = sum((xi.^2 + xi.*xip1 + xip1.^2) .* a) / 12;
            Ixy = sum((xi.*yip1 + 2*xi.*yi + 2*xip1.*yip1 + xip1.*yi) .* a) / 24;

            Iuu = Ixx - Area * yc^2;
            Ivv = Iyy - Area * xc^2;
            Iuv = Ixy - Area * xc * yc;
        end

        geom = [abs(Area), xc, yc, 0];
        iner = [Ixx, Iyy, Ixy, 0, 0, 0];
        cpmo = [abs(Iuu), abs(Ivv), Iuv, 0, 0, 0];

    catch ME
        warning('mwecmass:hydrostatics:PolygeomFailed', ...
                'Polygon geometry calculation failed: %s', ME.message);
        geom = [0, 0, 0, 0];
        iner = [0, 0, 0, 0, 0, 0];
        cpmo = [0, 0, 0, 0, 0, 0];
    end
end
