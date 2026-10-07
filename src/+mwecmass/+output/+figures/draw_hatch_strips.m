function draw_hatch_strips(ax, xp, yp, spacing, color, line_width, holes)
%DRAW_HATCH_STRIPS Fills a polygon with 45-degree hatch lines for the modular-precast strip plots.
% holes (optional): cell of [n x 2] loops inside the polygon that are left unhatched.
    if nargin < 7
        holes = {};
    end
    draw_hatch_local(ax, xp, yp, spacing, color, line_width, holes);
end

function draw_hatch_local(ax, xp, yp, spacing, color, line_width, holes)
%DRAW_HATCH_LOCAL Fill polygon (xp,yp), less the holes, with 45-degree hatching.
    if length(xp) < 3, return; end
    x_lo = min(xp);  x_hi = max(xp);
    y_lo = min(yp);  y_hi = max(yp);
    if (x_hi - x_lo) < 1e-9 || (y_hi - y_lo) < 1e-9, return; end
    c_min = x_lo - y_hi;
    c_max = x_hi - y_lo;
    c_vals = c_min:spacing:c_max;
    for ci = 1:length(c_vals)
        c = c_vals(ci);
        y_line = linspace(y_lo, y_hi, 200)';
        x_line = y_line + c;
        in = inpolygon(x_line, y_line, xp, yp);
        for h = 1:numel(holes)
            in = in & ~inpolygon(x_line, y_line, holes{h}(:, 1), holes{h}(:, 2));
        end
        if ~any(in), continue; end
        d = diff([false; in; false]);
        starts = find(d == 1);
        stops  = find(d == -1) - 1;
        n_seg = min(length(starts), length(stops));
        for s = 1:n_seg
            plot(ax, x_line(starts(s):stops(s)), ...
                     y_line(starts(s):stops(s)), '-', ...
                 'Color', color, 'LineWidth', line_width, ...
                 'HandleVisibility', 'off');
        end
    end
end
