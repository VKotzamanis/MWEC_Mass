function style_legend(lg, style)
%STYLE_LEGEND  Apply the shared type settings to one legend.
% Inputs: lg is a legend handle; style is out.style.
% Outputs: none; legend text is updated in points (pt).
    set(lg, 'FontName', style.font_name, ...
            'FontSize', style.font_size.legend, ...
            'Interpreter', style.tick_label_interpreter);
end
