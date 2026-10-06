function fid = open_log(out, type, name)
%OPEN_LOG Create Output/<type>/ as needed and open one UTF-8 report log for writing.
    repo_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    log_dir = fullfile(repo_root, out.output_dir, type);
    if ~exist(log_dir, 'dir')
        % Capture mkdir's own [ok,msg] rather than relying on it to throw: an identified,
        % path-carrying error is more actionable than MATLAB's own uncaptured-mkdir error.
        [ok, msg] = mkdir(log_dir);
        if ~ok
            error('mwecmass:output:OpenLogFailed', ...
                'Could not create log directory %s: %s', log_dir, msg);
        end
    end

    log_file = fullfile(log_dir, [name, '.log']);
    fid = fopen(log_file, 'w', 'n', 'UTF-8');
    if fid < 0
        error('mwecmass:output:OpenLogFailed', 'Could not open log file: %s', log_file);
    end
end
