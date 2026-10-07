% Run every tests/**/test_*.m (except tests/octave_shims) and exit non-zero if any fails.
%   octave --no-gui --quiet tests/run_tests.m
% Environment variable TESTS_FILTER (optional): run only files whose path contains that text.
% Each test is a function file with no inputs or outputs that calls error() on failure.
% No local functions here: tests such as test_pipeline_baseline end with "clear functions", which
% would remove functions defined in this script file.
tests_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(tests_dir);
warning('off', 'Octave:shadowed-function');
addpath(fullfile(tests_dir, 'octave_shims'));
addpath(fullfile(repo_root, 'src'));
addpath(repo_root);

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
filter = getenv('TESTS_FILTER');
if ~isempty(filter)
  files = files(~cellfun(@isempty, strfind(files, filter)));
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
failed = {};
t_all = tic;
for k = 1:numel(files)
  rel = files{k}(numel(repo_root) + 2:end);
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
fprintf('%d passed, %d failed, %d total, %.1f s\n', numel(files) - n_fail, n_fail, numel(files), toc(t_all));
if n_fail > 0
  fprintf('Failed:\n');
  fprintf('  %s\n', failed{:});
  exit(1);
end
exit(0);
