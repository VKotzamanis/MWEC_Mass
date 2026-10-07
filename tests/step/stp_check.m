function r = stp_check(file, varargin)
%STP_CHECK  Run tests/step_check.py on a STEP file and return its JSON output as a struct.

here = fileparts(mfilename('fullpath'));
script = fullfile(here, '..', 'step_check.py');
[status, out] = system(sprintf('python3 -I "%s" "%s" %s', script, file, strjoin(varargin, ' ')));
if status ~= 0
    error('stp_check: step_check.py failed:\n%s', out);
end
r = jsondecode(out);
end
