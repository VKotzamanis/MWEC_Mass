function cmap = figure_colormap(style, n)
%FIGURE_COLORMAP Return the n-row colour matrix of the named colormap (cividis, parula, or jet).
% style.colormap (default: cividis); n is the sample count. Returns RGB triplets in [0,1].
    name = 'cividis';
    if isstruct(style) && isfield(style, 'colormap') && ~isempty(style.colormap)
        name = char(style.colormap);
    end
    switch lower(name)
        case 'cividis'
            cmap = mwecmass.output.cividis_map(n);
        case 'parula'
            cmap = parula(n);
        case 'jet'
            cmap = jet(n);
        otherwise
            error('mwecmass:output:figures:unknownColormap', ...
                  'style.colormap = ''%s'' is not one of cividis/parula/jet.', name);
    end
end
