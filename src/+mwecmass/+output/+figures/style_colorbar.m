function style_colorbar(cb, style)
%STYLE_COLORBAR  Apply the shared type settings to one colorbar and its Label.
% Inputs: cb is a colorbar handle; style is out.style.
% Outputs: none; colorbar tick labels are updated in points (pt); the colorbar's Label text
%          object (always present, even with an empty String) is styled through style_text's
%          'label' role (style.font_size.annotation, no dedicated font_size.label field exists).
    set(cb, 'FontName', style.font_name, ...
            'FontSize', style.font_size.colorbar, ...
            'TickLabelInterpreter', style.tick_label_interpreter);
    mwecmass.output.figures.style_text(cb.Label, style, 'label');
end
