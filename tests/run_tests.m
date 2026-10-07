% Run every tests/**/test_*.m (except tests/octave_shims) and exit non-zero if any fails.
%   octave --no-gui --quiet tests/run_tests.m
% Environment variables:
%   TESTS_FILTER (optional)  run only files whose path contains that text.
%   MWEC_REGRESSION=1        also run tests/regression/test_*.m (they run the whole pipeline, 24 to 94
%                            minutes per mode); without it they are listed as skipped.
% Each test is a function file with no inputs or outputs that calls error() on failure.
% No local functions here: baseline_run clears functions, which would remove functions defined in
% this script file.
tests_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(tests_dir);
warning('off', 'Octave:shadowed-function');
% One call so that the shims come first, then src, then the repo root (addpath prepends).
addpath(fullfile(tests_dir, 'octave_shims'), fullfile(repo_root, 'src'), repo_root);

files = {};
pending = {tests_dir};
while ~isempty(pending)
  folder = pending{end};
  pending(end) = [];
  entries = dir(folder);
  for k = 1:numel(entries)
    name = entries(k).name;
    full = fullfile(folder, name);
    if entries(k).isdir
      if ~any(strcmp(name, {'.', '..', 'octave_shims'}))
        pending{end+1} = full;
      end
    elseif ~isempty(regexp(name, '^test_.*\.m$', 'once'))
      files{end+1} = full;
    end
  end
end
files = sort(files);
run_regression = strcmp(getenv('MWEC_REGRESSION'), '1');
regression_dir = [fullfile(tests_dir, 'regression') filesep];
is_regression = strncmp(files, regression_dir, numel(regression_dir));
filter = getenv('TESTS_FILTER');
if ~isempty(filter)
  keep = ~cellfun(@isempty, strfind(files, filter));
  files = files(keep);
  is_regression = is_regression(keep);
end
if isempty(files)
  fprintf('run_tests: no test files found.\n');
  exit(1);
end
[~, names] = cellfun(@fileparts, files, 'UniformOutput', false);
if numel(unique(names)) ~= numel(names)
  fprintf('run_tests: test file names must be unique across folders.\n');
  exit(1);
end

n_fail = 0;
n_skip = 0;
failed = {};
t_all = tic;
for k = 1:numel(files)
  rel = files{k}(numel(repo_root) + 2:end);
  if is_regression(k) && ~run_regression
    n_skip = n_skip + 1;
    fprintf('SKIP  %s  (pipeline regression; set MWEC_REGRESSION=1 to run it)\n\n', rel);
    continue;
  end
  fprintf('---- %s\n', rel);
  t0 = tic;
  start_dir = pwd;
  addpath(fileparts(files{k}));
  ok = true;
  message = '';
  try
    feval(names{k});
  catch err
    ok = false;
    message = err.message;
  end
  rmpath(fileparts(files{k}));
  cd(start_dir);
  close all;
  dt = toc(t0);
  if ok
    fprintf('PASS  %s  (%.1f s)\n\n', rel, dt);
  else
    n_fail = n_fail + 1;
    failed{end+1} = rel;
    fprintf('FAIL  %s  (%.1f s)\n      %s\n\n', rel, dt, strrep(message, sprintf('\n'), sprintf('\n      ')));
  end
end
fprintf('%d passed, %d failed, %d skipped, %d total, %.1f s\n', ...
        numel(files) - n_fail - n_skip, n_fail, n_skip, numel(files), toc(t_all));
if n_fail > 0
  fprintf('Failed:\n');
  fprintf('  %s\n', failed{:});
  exit(1);
end
exit(0);
