function write_octave_baseline(modes)
%WRITE_OCTAVE_BASELINE Run the pipeline and write tests/baseline/octave_v1_baseline.json.
%   write_octave_baseline() runs both realisation types; write_octave_baseline({'thin_shell'})
%   runs the listed types only and keeps the other entries of an existing file, so the two
%   types can be produced by two Octave processes started side by side. From the repository root:
%     octave --no-gui --quiet --eval "addpath('tools'); write_octave_baseline({'thin_shell'})"
%   File layout: one JSON member per line, "info", the modes, then "seconds_<mode>". The text of a
%   mode line is jsonencode(summary), the shortest text that reads back to the same double, and
%   tests/regression/test_pipeline_baseline.m compares that text exactly.
  repo = fileparts(fileparts(mfilename('fullpath')));
  addpath(fullfile(repo, 'tests', 'octave_shims'), fullfile(repo, 'src'), repo);
  warning('off', 'Octave:shadowed-function');
  all_modes = {'modular_precast', 'thin_shell'};
  if nargin < 1, modes = all_modes; end
  path = fullfile(repo, 'tests', 'baseline', 'octave_v1_baseline.json');

  for k = 1:numel(modes)
    t0 = tic;
    summary = baseline_run(modes{k});
    seconds = toc(t0);
    fprintf('%s: %.0f s\n', modes{k}, seconds);
    members = read_members(path);
    members(modes{k}) = jsonencode(summary);
    members(['seconds_' modes{k}]) = sprintf('%.1f', seconds);
    members('info') = jsonencode(struct('octave', version(), 'source', ...
                      'mwecmass.driver.run under Octave via tools/baseline_run.m'));
    write_members(path, members, all_modes);
    fprintf('wrote %s\n', path);
  end
end

function members = read_members(path)
  members = containers.Map();
  if ~exist(path, 'file'), return; end
  lines = strsplit(fileread(path), sprintf('\n'));
  for k = 1:numel(lines)
    tok = regexp(lines{k}, '^"([^"]+)": (.*?),?$', 'tokens', 'once');
    if ~isempty(tok), members(tok{1}) = tok{2}; end
  end
end

function write_members(path, members, all_modes)
  order = ['info', all_modes, strcat('seconds_', all_modes)];
  order = order(cellfun(@(n) isKey(members, n), order));
  lines = cell(1, numel(order));
  for k = 1:numel(order)
    lines{k} = sprintf('"%s": %s', order{k}, members(order{k}));
  end
  fid = fopen(path, 'w');
  fprintf(fid, '{\n%s\n}\n', strjoin(lines, sprintf(',\n')));
  fclose(fid);
end
