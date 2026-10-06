function v = option_or_config(opts, name, cfg_value)
%OPTION_OR_CONFIG Select a nonempty option or its configuration fallback.
% Inputs: opts struct, name field-name character vector, cfg_value fallback.
% Output v is opts.(name) when present and nonempty; otherwise cfg_value.
    if isfield(opts, name) && ~isempty(opts.(name))
        v = opts.(name);
    else
        v = cfg_value;
    end
end
