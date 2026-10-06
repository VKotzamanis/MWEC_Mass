function [status, result] = run_solver(hams_exe, run_dir)
%RUN_SOLVER  Execute HAMS-MREL with explicit input/output directories.
% hams_exe and run_dir are paths; the command is exe run_dir/Input run_dir/Output.
% Outputs are status (0 normally means success) and captured console result.  On non-Windows
% systems, Intel oneAPI must be present in LD_LIBRARY_PATH (source setvars.sh before MATLAB).
% Because HAMS can return status 0 after input errors, known abort markers in stdout also
% convert the result to failure.  See docs/HAMS_MREL_ROUTE.md for the calling convention.

    assert(exist(hams_exe, 'file') == 2, ...
        'HAMS executable not found: %s', hams_exe);

    % LINUX/MAC PRE-FLIGHT: verify Intel oneAPI is on LD_LIBRARY_PATH.
    % Bail out with an actionable message before invoking system();
    % otherwise the binary aborts with a cryptic "libmkl_intel_lp64.so.3:
    % cannot open shared object file" and the only signal upstream is
    % a non-zero status code.
    if ~ispc
        ld_path = getenv('LD_LIBRARY_PATH');
        if ~contains(ld_path, 'intel/oneapi') && ...
                ~contains(ld_path, 'intel\oneapi')
            error('mwecmass:hams_mrel:OneAPIEnvMissing', ...
                ['Intel oneAPI runtime not on LD_LIBRARY_PATH.\n' ...
                 'HAMS-MREL is linked against Intel MKL + iomp5 and ' ...
                 'cannot load them without setvars.sh.\n\n' ...
                 'Fix: quit MATLAB, then from a terminal run\n' ...
                 '    source /opt/intel/oneapi/setvars.sh\n' ...
                 '    matlab\n' ...
                 '(adjust the setvars path if your oneAPI install is elsewhere).']);
        end
    end

    input_dir  = fullfile(run_dir, 'Input');
    output_dir = fullfile(run_dir, 'Output');

    cmd = sprintf('"%s" "%s" "%s"', hams_exe, input_dir, output_dir);

    fprintf('  Running HAMS: %s\n', cmd);
    tic;
    [status, result] = system(cmd);
    elapsed = toc;

    % HAMS-MREL returns exit code 0 even on input-file errors (it
    % prints the diagnostic, then `stop` falls through with status
    % 0).  Scan stdout for the known abort markers so consumers
    % don't see "success" when no .1 file was produced.
    failure_markers = { ...
        'Terminating application', ...
        'Input files missing', ...
        'input file missing', ...
        'Error opening', ...
        'Error encountered reading'};
    stdout_failed = false;
    for k = 1:numel(failure_markers)
        if contains(result, failure_markers{k})
            stdout_failed = true;
            break;
        end
    end

    if status == 0 && ~stdout_failed
        fprintf('  HAMS completed in %.1f seconds\n', elapsed);
    else
        if status == 0 && stdout_failed
            status = 1;   % surface stdout failure to caller
        end
        warning('mwecmass:hams_mrel:HAMSFailed', ...
            'HAMS aborted (status=%d).\nOutput:\n%s', status, result);
    end
end
