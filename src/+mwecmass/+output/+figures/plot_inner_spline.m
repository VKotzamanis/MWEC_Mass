function plot_inner_spline(ax, inner_prof, style)
%PLOT_INNER_SPLINE Plot the inner profile as a smooth spline curve.
% Inputs: ax (axes handle); inner_prof (Nx2, [x z] in body frame); style (optional, defaults to presentation_style). Curves through >= 3 points; returns silently if profile is empty or degenerate.

    if isempty(inner_prof) || size(inner_prof,1) < 3, return; end
    if nargin < 3 || isempty(style)
        style = mwecmass.output.figures.presentation_style();
    end
    x = inner_prof(:,1);
    z = inner_prof(:,2);
    x_cl = [x; x(1)];
    z_cl = [z; z(1)];
    segs = sqrt(diff(x_cl).^2 + diff(z_cl).^2);
    segs = max(segs, 1e-12);
    t_param = [0; cumsum(segs)];
    t_fine  = linspace(0, t_param(end), 400);
    x_sp = interp1(t_param, x_cl, t_fine, 'spline');
    z_sp = interp1(t_param, z_cl, t_fine, 'spline');
    plot(ax, x_sp, z_sp, '--', ...
         'Color', style.fill_palette.inner_boundary, ...
         'LineWidth', style.line_width.boundary, ...
         'DisplayName', 'Shell boundary');
end
