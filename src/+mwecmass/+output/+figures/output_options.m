function opts = output_options(config)
%OUTPUT_OPTIONS Return output options for figures: save flags, formats, style, timestamp convention.
% config: run config (carries .output), output-options struct, or empty (loads repository defaults).
% Returns struct with save, format, style, timestamp_filenames, output_dir. Units: font sizes in pt;
% line widths in pt; image resolution in dpi; colours as RGB triplets [0,1]; hatch spacing in m.
    opts = [];
    if nargin >= 1 && isstruct(config) && isscalar(config)
        if isfield(config, 'output') && isstruct(config.output) && isscalar(config.output)
            opts = config.output;
        elseif isfield(config, 'save') || isfield(config, 'style')
            opts = config;
        end
    end

    defaults = root_defaults();
    if isempty(opts)
        opts = defaults;
    else
        names = fieldnames(defaults);
        for k = 1:numel(names)
            if ~isfield(opts, names{k})
                opts.(names{k}) = defaults.(names{k});
            end
        end
    end
end


function out = root_defaults()
%ROOT_DEFAULTS Fetch the repository's default output options from the root options file.
% Resolves the repository root from this file's location (the driver adds only src/).
    if exist('WEC_Output_Options', 'file') ~= 2
        here = fileparts(mfilename('fullpath'));   % src/+mwecmass/+output/+figures
        addpath(fileparts(fileparts(fileparts(fileparts(here)))));
    end
    out = WEC_Output_Options();
end
