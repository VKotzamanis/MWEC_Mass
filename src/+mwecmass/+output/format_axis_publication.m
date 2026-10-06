function format_axis_publication(ax, style)
%FORMAT_AXIS_PUBLICATION Apply presentation style to axes, colorbar, and legend.
% Sets grid, box, font, size [pt], line width [pt], and interpreter.
% Font size cascade: title and labels scale from axes size via TitleFontSizeMultiplier
% and LabelFontSizeMultiplier (1.0, so labels match tick numbers).
% Colorbar and legend sizes are decoupled from axes (MATLAB defaults to 90%).
% style struct optional; defaults to presentation_style() when omitted.
    if nargin < 2 || isempty(style)
        style = mwecmass.output.figures.presentation_style();
    end
    grid(ax, 'on'); box(ax, 'on');
    set(ax, 'GridLineStyle', style.grid.line_style, ...
            'GridAlpha', style.grid.alpha, ...
            'LineWidth', style.line_width.axes, ...
            'FontName', style.font_name, ...
            'FontSize', style.font_size.axes, ...
            'LabelFontSizeMultiplier', 1, ...
            'TitleFontSizeMultiplier', style.font_size.title / style.font_size.axes, ...
            'TickLabelInterpreter', style.tick_label_interpreter);
    if ~isempty(ax.Colorbar)
        ax.Colorbar.FontSize = style.font_size.colorbar;
    end
    if ~isempty(ax.Legend)
        ax.Legend.FontSize = style.font_size.legend;
    end
end
