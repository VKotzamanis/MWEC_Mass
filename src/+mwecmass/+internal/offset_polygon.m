function [x_off, y_off] = offset_polygon(x, y, dist)
%OFFSET_POLYGON Offset a polygon inward by a constant distance.
% Vertices are treated as a planar polygon; winding is normalized to CCW,
% then adjacent inward normals are averaged with a bounded miter correction.
% Inputs x,y are vertex vectors and dist is the offset distance; outputs are
% offset vertices, unchanged for degenerate or nonpositive-distance inputs.

    x = x(:);  y = y(:);
    Nv = length(x);

    if Nv < 3 || dist <= 0
        x_off = x;  y_off = y;
        return;
    end

    % Remove duplicate closing vertex if present
    if abs(x(end)-x(1)) < 1e-12 && abs(y(end)-y(1)) < 1e-12
        x = x(1:end-1);
        y = y(1:end-1);
        Nv = length(x);
    end
    if Nv < 3
        x_off = [];  y_off = [];
        return;
    end

    % Ensure CCW winding
    signed_area = 0.5 * sum(x .* circshift(y,-1) - circshift(x,-1) .* y);
    if signed_area < 0
        x = flipud(x);
        y = flipud(y);
    end

    x_off = zeros(Nv, 1);
    y_off = zeros(Nv, 1);

    for j = 1:Nv
        jm = mod(j-2, Nv) + 1;
        jp = mod(j,   Nv) + 1;

        e_prev = [x(j)-x(jm), y(j)-y(jm)];
        e_next = [x(jp)-x(j), y(jp)-y(j)];

        len_prev = norm(e_prev);
        len_next = norm(e_next);

        if len_prev < 1e-12 || len_next < 1e-12
            x_off(j) = x(j);
            y_off(j) = y(j);
            continue;
        end

        % Inward normals (CCW: left of edge)
        n_prev = [-e_prev(2), e_prev(1)] / len_prev;
        n_next = [-e_next(2), e_next(1)] / len_next;

        n_avg = n_prev + n_next;
        len_avg = norm(n_avg);
        if len_avg < 1e-12
            n_avg = n_prev;
            len_avg = 1.0;
        end
        n_avg = n_avg / len_avg;

        % Miter correction
        cos_half = dot(n_avg, n_prev);
        if cos_half > 0.33
            miter_dist = dist / cos_half;
        else
            miter_dist = dist * 3.0;
        end

        x_off(j) = x(j) + n_avg(1) * miter_dist;
        y_off(j) = y(j) + n_avg(2) * miter_dist;
    end
end
