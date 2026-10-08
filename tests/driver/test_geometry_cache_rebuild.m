function test_geometry_cache_rebuild()
%TEST_GEOMETRY_CACHE_REBUILD Changing the deck, a keyed input or a keyed source file rebuilds.
%   build_config on the stand-in cylinder deck (thin shell): one scenario each for the deck, a
%   keyed input and a keyed source file (build_config run from a copy of src/ whose file is then
%   edited); after each change the next call recomputes and the one after reloads. Inputs the
%   cached steps do not read, and source files outside the keyed set, leave the cache valid.
  repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  deck = fullfile(repo_root, 'tests', 'standins', 'fixtures', 'cylinder.ms2');
  work = tempname();
  mkdir(work);
  saved_path = path();
  addpath(fullfile(repo_root, 'tests', 'standins'), '-end');   % kernel stand-ins until the real one is merged
  cleanup = onCleanup(@() restore(saved_path, work));
  miss = 'deck, inputs or source changed, computing';
  hit = 'Geometry products: reloaded from';

  % Deck bytes.
  cache_dir = fullfile(work, 'cache_deck');
  in = deck_input(deck, 'thin_shell', cache_dir);
  [log, config] = build(in);
  expect(log, 'no cache file, computing', 'first build');
  volume_before = config.total_wec_volume;
  [log, config] = build(in);
  expect(log, hit, 'second build');
  deck_copy = fullfile(work, 'cylinder.ms2');
  text = fileread(deck);
  changed = strrep(strrep(text, '1.5 0.0 -3.0', '1.6 0.0 -3.0'), '1.5 0.0 1.0', '1.6 0.0 1.0');
  if strcmp(changed, text), error('the stand-in deck no longer holds the text this test edits'); end
  fid = fopen(deck_copy, 'w'); fwrite(fid, changed); fclose(fid);
  in_deck = deck_input(deck_copy, 'thin_shell', cache_dir);
  [log, config] = build(in_deck);
  expect(log, miss, 'deck changed');
  volume_after = config.total_wec_volume;
  [log, config] = build(in_deck);
  expect(log, hit, 'changed deck repeated');
  [log, config] = build(in);
  expect(log, miss, 'original deck restored');
  fprintf('deck: changed points P1.x and P2.x 1.5 -> 1.6 m, hull volume %.4f -> %.4f m^3, rebuilt\n', ...
          volume_before, volume_after);

  % Keyed input.
  in_keyed = in;
  in_keyed.geometry.n_z_levels = in.geometry.n_z_levels + 1;
  [log, config] = build(in_keyed);
  expect(log, miss, 'keyed input changed');
  if numel(config.y_span_table) ~= in_keyed.geometry.n_z_levels
    error('the rebuild did not use the changed input (y_span_table has %d rows)', numel(config.y_span_table));
  end
  [log, config] = build(in_keyed);
  expect(log, hit, 'keyed input repeated');
  in_other = in_keyed;
  in_other.targets.T_heave_goal = in_keyed.targets.T_heave_goal + 1;
  in_other.materials.modular_precast.wall_height = 2.0;   % precast-only: not read in thin-shell mode
  [log, config] = build(in_other);
  expect(log, hit, 'inputs the cached steps do not read changed');
  if config.T_heave_goal ~= in_other.targets.T_heave_goal
    error('a reloaded config does not carry the current non-geometry inputs');
  end
  fprintf('input: n_z_levels %d -> %d rebuilt (y_span_table %d rows); unread inputs reloaded\n', ...
          in.geometry.n_z_levels, in_keyed.geometry.n_z_levels, numel(config.y_span_table));

  % Keyed source file: run build_config from a copy of src whose file is edited afterwards.
  copy_root = fullfile(work, 'copy');
  mkdir(copy_root);
  mkdir(fullfile(copy_root, 'Input'));
  copyfile(fullfile(repo_root, 'src'), fullfile(copy_root, 'src'));
  addpath(fullfile(copy_root, 'src'));
  clear config
  clear functions
  in_copy = deck_input(deck, 'thin_shell', fullfile(work, 'cache_source'), copy_root);
  [log, config] = build(in_copy);
  expect(log, 'no cache file, computing', 'build from the copy');
  [log, config] = build(in_copy);
  expect(log, hit, 'build from the copy repeated');
  append_text(fullfile(copy_root, 'src', '+mwecmass', '+optim', 'stage2_bounds.m'), sprintf('\n%% edit\n'));
  [log, config] = build(in_copy);
  expect(log, hit, 'file outside the keyed set edited');
  append_text(fullfile(copy_root, 'src', '+mwecmass', '+hydrostatics', 'compute_strip.m'), sprintf('\n%% edit\n'));
  clear config
  clear functions
  [log, config] = build(in_copy);
  expect(log, miss, 'keyed source file edited');
  [log, config] = build(in_copy);
  expect(log, hit, 'keyed source file edit repeated');
  fprintf('source: +hydrostatics/compute_strip.m edited, rebuilt; +optim/stage2_bounds.m edited, reloaded\n');
end

function [log, config] = build(in)
  log = evalc('config = mwecmass.driver.build_config(in, [], struct());');
end

function expect(log, text, step)
  if isempty(strfind(log, text))
    error('%s: expected "%s" in the console output, got:\n%s', step, text, log);
  end
end

function append_text(file, text)
  fid = fopen(file, 'a');
  fwrite(fid, text);
  fclose(fid);
end

function restore(saved_path, work)
  path(saved_path);
  clear functions
  if exist(work, 'dir')
    confirm_recursive_rmdir(false, 'local');
    rmdir(work, 's');
  end
end
