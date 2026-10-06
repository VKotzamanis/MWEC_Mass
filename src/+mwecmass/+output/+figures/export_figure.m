function [files, base_name] = export_figure(fig, name_stem, config, out_dir)
%EXPORT_FIGURE Write a figure to each configured format under a repeatable file name.
% Writes one file per format from config (default: current folder). PNG resolution in dpi,
% PDF content type, FIG as MATLAB object. Returns file paths and the base name for companions.
    if nargin < 4
        out_dir = '';
    end
    opts = mwecmass.output.figures.output_options(config);

    base_name = name_stem;
    if opts.export.timestamp_filenames
        base_name = sprintf('%s_%s', name_stem, ...
                            char(datetime('now', 'Format', 'yyyyMMdd_HHmmss')));
    end

    exts = opts.export.formats;
    if ~iscell(exts)
        exts = cellstr(exts);
    end

    files = cell(1, numel(exts));
    for k = 1:numel(exts)
        if isempty(out_dir)
            files{k} = sprintf('%s.%s', base_name, char(exts{k}));
        else
            files{k} = fullfile(out_dir, sprintf('%s.%s', base_name, char(exts{k})));
        end
        switch lower(char(exts{k}))
            case 'png'
                exportgraphics(fig, files{k}, 'Resolution', opts.export.dpi);
            case 'pdf'
                exportgraphics(fig, files{k}, 'ContentType', opts.export.pdf_content);
            case 'fig'
                savefig(fig, files{k});
            otherwise
                error('mwecmass:figures:UnsupportedExportFormat', ...
                      'Unsupported out.export format: %s.', char(exts{k}));
        end
    end
end
