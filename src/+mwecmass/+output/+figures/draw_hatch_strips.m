function draw_hatch_strips(ax, xp, yp, spacing, color, line_width)
%DRAW_HATCH_STRIPS Fills a polygon with 45-degree hatch lines for the modular-precast strip plots.
    draw_hatch_local(ax, xp, yp, spacing, color, line_width);
end

function draw_hatch_local(ax, xp, yp, spacing, color, line_width)
%DRAW_HATCH_LOCAL Fill polygon (xp,yp) with 45-degree hatching.
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
