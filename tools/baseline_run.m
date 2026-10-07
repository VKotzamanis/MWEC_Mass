function summary = baseline_run(mode, keep_log, overrides)
%BASELINE_RUN Run mwecmass.driver.run headless under Octave and return the key scalars.
%   summary = baseline_run(mode) for mode 'modular_precast' or 'thin_shell' on Input/C1.ms2 and
%   the Octave-readable copy of its BEM cache (tests/fixtures/C1_wamit_cache_v5.mat, made by
%   tools/convert_v73_to_v5.py because Octave cannot read the cell arrays of the v7.3 file).
%
%   The pipeline derives its output folder from the location of src/, so the run uses a copy of
%   src/, the deck, the cache and validation/ in a temporary folder; Output/ of the repository is
%   never touched. Every figure, log and diagnostic switch in out.save is false; only the results
%   MAT-file is written, into the temporary folder, and read back for the numbers. The temporary
%   folder is removed afterwards. keep_log (default false) prints the pipeline console text.
%   overrides is a cell array of {'dotted.path.in.in', value; ...} applied to the input struct, or
%   the name of a preset: 'full' (the default) keeps the author inputs, 'fast' also coarsens the
%   z-grids (n_z_levels 40, n_z_grid 60) so that a rerun takes less time. The applied overrides
%   are returned in the summary.
%
%   Test-only: requires Octave with tests/octave_shims on the path.
  if nargin < 2 || isempty(keep_log), keep_log = false; end
  if ~any(strcmp(mode, {'modular_precast', 'thin_shell'}))
    error('baseline_run:mode', 'mode must be modular_precast or thin_shell.');
  end
  repo = fileparts(fileparts(mfilename('fullpath')));
  cache_v5 = fullfile(repo, 'tests', 'fixtures', 'C1_wamit_cache_v5.mat');
  if ~exist(cache_v5, 'file')
    error('baseline_run:cache', ['Missing %s. Make it with: python3 tools/convert_v73_to_v5.py ' ...
          'Input/C1_wamit_cache.mat tests/fixtures/C1_wamit_cache_v5.mat hydro_table'], cache_v5);
  end

  tmp = tempname();
  mkdir(tmp);
  mkdir(fullfile(tmp, 'Input'));
  mkdir(fullfile(tmp, 'Output'));
  copyfile(fullfile(repo, 'src'), fullfile(tmp, 'src'));
  copyfile(fullfile(repo, 'validation'), fullfile(tmp, 'validation'));
  copyfile(fullfile(repo, 'Input', 'C1.ms2'), fullfile(tmp, 'Input', 'C1.ms2'));
  copyfile(cache_v5, fullfile(tmp, 'Input', 'C1_wamit_cache.mat'));

  old_path = path();
  start_dir = pwd();
  cleanup = onCleanup(@() restore(old_path, start_dir, tmp));
  rmpath(fullfile(repo, 'src'));
  addpath(fullfile(tmp, 'src'));
  addpath(repo);
  clear functions;

  in = WEC_User_Input();
  if nargin < 3, overrides = 'full'; end
  if ischar(overrides), overrides = preset_overrides(mode, overrides); end
  in.materials.realisation_type = mode;
  for k = 1:size(overrides, 1)
    parts = strsplit(overrides{k, 1}, '.');
    in = setfield(in, parts{:}, overrides{k, 2}); %#ok<SFLD>
  end
  out = WEC_Output_Options();
  out.save = set_all(out.save, false);
  out.save.results_mat = true;
  out.console_echo = false;

  log_file = fullfile(tmp, 'console.txt');
  if ~keep_log, diary(log_file); diary on; end
  try
    mwecmass.driver.run(in, out);
  catch err
    if ~keep_log
      diary off;
      fprintf('---- last pipeline console lines before the error ----\n');
      fprintf('%s\n', tail_of(log_file, 40));
    end
    rethrow(err);
  end
  if ~keep_log, diary off; end

  loaded = load(fullfile(tmp, 'Output', ['C1_' mode '_results.mat']));
  summary = summarise(mode, loaded.results, loaded.final_props);
  summary.input_overrides = struct();
  for k = 1:size(overrides, 1)
    summary.input_overrides.(strrep(overrides{k, 1}, '.', '__')) = overrides{k, 2};
  end
end

function overrides = preset_overrides(mode, preset)
% The default n_sub = 100 makes build_config's circle-based floor code (compute_rmin_at_z.m:30-35)
% read a centroid from a polygon area of -1.1e-14 m^2 (pure roundoff, at z = -1.756 m in strip 3)
% and return a floor of 2500 kg/m^3 for strip 3, so the run stops with 'cannot float'.
% n_sub = 101 moves the sample heights off that point; this code is deleted by task T4.
  overrides = cell(0, 2);
  if strcmp(mode, 'modular_precast')
    overrides = {'materials.modular_precast.n_sub', 101};
  end
  switch preset
    case 'full'
    case 'fast'
      overrides = [overrides; {'geometry.n_z_levels', 40; ...
                               'materials.thin_shell.n_z_grid', 60; ...
                               'materials.modular_precast.n_z_grid', 60}];
    otherwise
      error('baseline_run:preset', 'unknown preset %s (full or fast).', preset);
  end
end

function restore(old_path, start_dir, tmp)
  diary off;
  path(old_path);
  cd(start_dir);
  clear functions;
  confirm_recursive_rmdir(false, 'local');
  rmdir(tmp, 's');
end

function s = set_all(s, value)
  names = fieldnames(s);
  for k = 1:numel(names)
    if isstruct(s.(names{k}))
      s.(names{k}) = set_all(s.(names{k}), value);
    else
      s.(names{k}) = value;
    end
  end
end

function txt = tail_of(file, n)
  txt = fileread(file);
  lines = strsplit(txt, sprintf('\n'));
  txt = strjoin(lines(max(1, end - n + 1):end), sprintf('\n'));
end

function summary = summarise(mode, results, final_props)
  s2 = results.stage2_3d;
  if strcmp(mode, 'thin_shell')
    c = results.steel_data;
  else
    c = results.constructability;
  end
  summary = struct();
  summary.mode = mode;
  summary.stage1_x = row(results.stage1_2d.x_optimal);
  summary.stage2_x = row(s2.x_optimal);
  summary.stage2_exitflag = s2.exitflag;
  summary.stage2_fval = s2.fval;
  summary.stage2_iterations = s2.output.iterations;
  summary.mass_total = final_props.mass_total;
  summary.mass_buoyant_force = final_props.mass_buoyant_force;
  summary.CG_total_z = final_props.CG_total(3);
  summary.GM_L = final_props.GM_L;
  summary.T_heave_coupled = final_props.coupled_periods(2);
  summary.T_pitch_coupled = final_props.coupled_periods(3);
  summary.vertical_shift = final_props.vertical_shift;
  summary.draft = final_props.draft;
  summary.realised_strip_density = row(final_props.realised_strip_density);
  stage3 = struct();
  names = {'t_steel', 'vertical_shift', 'draft', 't_max', 't_min_active', 'V_steel', ...
           'V_air', 'M_total', 'CG_z_world', 'GM_realised', 'T_heave_realised', ...
           'T_pitch_realised', 'mass_balance_error_pct', 'phi_star', 'feasible', 'exitflag'};
  for k = 1:numel(names)
    if isfield(c, names{k})
      stage3.(names{k}) = double(c.(names{k}));
    end
  end
  % The ballast level field is read under either of its two names.
  for name = {'z_ballast', 'z_fill'}
    if isfield(c, name{1})
      stage3.z_ballast = double(c.(name{1}));
      break;
    end
  end
  if isfield(c, 't_offset_strip')
    t = row(c.t_offset_strip);
    stage3.t_offset_strip_finite = t(isfinite(t));
  end
  summary.stage3 = stage3;
end

function r = row(v)
  r = double(v(:)');
end
