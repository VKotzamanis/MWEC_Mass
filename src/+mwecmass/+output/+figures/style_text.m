function style_text(h, style, role, font_name)
%STYLE_TEXT  Apply the shared font to a free-standing text, annotation, xline/yline label, or
%tiledlayout title/subtitle handle.
% Inputs: h is a text/annotation/ConstantLine/layout-title handle (or array); style is out.style;
%         role is 'annotation', 'title', or 'label'; font_name is an optional override (used for
%         a monospace numeric card, out.style.mono_font_name), defaulting to style.font_name.
% Outputs: none; FontName, FontSize and Interpreter are updated in place. FontSize is read from
%          style.font_size.(role), falling back to style.font_size.annotation when the role has
%          no dedicated font-size field (role 'label': no font_size.label exists).
    valid_roles = {'annotation', 'title', 'label'};
    if ~any(strcmp(role, valid_roles))
        error('mwecmass:figures:UnknownTextRole', ...
              'role must be annotation, title, or label.');
    end
    if nargin < 4 || isempty(font_name)
        font_name = style.font_name;
    end
    if isfield(style.font_size, role)
        font_size = style.font_size.(role);
    else
        font_size = style.font_size.annotation;
    end
    set(h, 'FontName', font_name, 'FontSize', font_size, ...
           'Interpreter', style.tick_label_interpreter);
end
