function fig = new_figure(style, size_name)
%NEW_FIGURE  Create a white publication figure at one named print size.
% Inputs: style is out.style; size_name is 'single_column' or 'double_column'.
% Outputs: fig is a figure handle sized in centimetres (cm).
    if ~isfield(style, 'figure_size') || ~isfield(style.figure_size, size_name)
        error('mwecmass:figures:UnknownSizeClass', ...
              'style.figure_size.%s is required.', size_name);
    end
    size_cm = style.figure_size.(size_name);
    fig = figure('Color', 'white', 'Units', 'centimeters', ...
                 'Position', [2, 2, size_cm]);
end
