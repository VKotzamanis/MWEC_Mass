function style = presentation_style(config)
%PRESENTATION_STYLE Return the presentation style struct: fonts, sizes, colours, line widths.
% Every figure function calls this at the top to source fonts [pt], line widths [pt],
% RGB colours [0,1], and hatch spacing [m] from config.output.style or repository defaults.
% Inputs: config (struct or empty, see output_options). Outputs: style (struct from opts.style).
    if nargin < 1
        config = [];
    end
    opts  = mwecmass.output.figures.output_options(config);
    style = opts.style;
end
