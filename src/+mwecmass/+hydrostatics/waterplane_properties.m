function [Aw, I_xx, I_yy, pts_ordered] = waterplane_properties(wl_pts)
%WATERPLANE_PROPERTIES Compute polygon waterplane area and second moments.
% wl_pts [N x 3] are boundary samples; points are angularly sorted about
% their arithmetic-mean centroid before Green's-theorem shoelace sums.
% Outputs Aw [m^2], I_xx=∫y^2 dA and I_yy=∫x^2 dA [m^4], plus ordered points.
% Fewer than three points returns zero properties.
% See docs/METHODS_ENGINE.md#waterplane-green-shoelace.

    if isempty(wl_pts) || size(wl_pts, 1) < 3
        Aw = 0; I_xx = 0; I_yy = 0; pts_ordered = wl_pts;
        return;
    end

    % Sort boundary points by angle around the centroid.
    x = wl_pts(:, 1);
    y = wl_pts(:, 2);
    cx = mean(x);
    cy = mean(y);
    angles = atan2(y - cy, x - cx);
    [~, order] = sort(angles);
    x = x(order);
    y = y(order);
    pts_ordered = wl_pts(order, :);

    x_next = [x(2:end); x(1)];
    y_next = [y(2:end); y(1)];

    % Compute area with the shoelace form of Green's theorem.
    cross_terms = x .* y_next - x_next .* y;
    Aw = abs(sum(cross_terms)) / 2;

    % Compute second moments with Green's theorem extensions.
    %  These are the standard formulas for polygon moments:
    %    I_xx = ∫∫ y² dA    (about x-axis)
    %    I_yy = ∫∫ x² dA    (about y-axis)
    I_xx = abs(sum(cross_terms .* (y.^2 + y .* y_next + y_next.^2))) / 12;
    I_yy = abs(sum(cross_terms .* (x.^2 + x .* x_next + x_next.^2))) / 12;
end
