function d = output_dir(type)
%OUTPUT_DIR Resolve and create Output/<type>/ under the repository root.
    repo_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    d = fullfile(repo_root, 'Output', type);
    if ~exist(d, 'dir')
        [ok, msg] = mkdir(d);
        % Fail loudly: a silent repo_root fallback here would misdirect every artefact this
        % run's two callers write, with no signal at the write site that it happened.
        if ~ok
            error('mwecmass:io:output_dir:MkdirFailed', ...
                    'Could not create %s: %s', d, msg);
        end
    end
end
