function test_geometry_cache_reload()
%TEST_GEOMETRY_CACHE_RELOAD A config reloaded from the geometry cache equals a freshly built one.
%   Both realisation paths of build_config (modular precast with the wall-pinned layout and the
%   per-strip bounds; thin shell with the uniform layout) on the stand-in cylinder deck: first
%   call computes and saves, second reloads. Equality is isequaln (NaN equal to NaN; z_cg_target is NaN), exact by identity of two code paths.
%   config.hydro_cache is left out: empty_hydro_cache stamps the wall-clock time into it, outside
%   the cached block. The saved file must hold plain data only.
  repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  deck = fullfile(repo_root, 'tests', 'standins', 'fixtures', 'cylinder.ms2');
  cache_dir = tempname();
  saved_path = path();
  addpath(fullfile(repo_root, 'tests', 'standins'), '-end');   % kernel stand-ins until the real one is merged
  cleanup = onCleanup(@() restore(saved_path, cache_dir));

  modes = {'modular_precast', 'thin_shell'};
  for m = 1:numel(modes)
    in = deck_input(deck, modes{m}, cache_dir);
    cache_file = fullfile(cache_dir, ['cylinder_' modes{m} '_geometry.mat']);
    if exist(cache_file, 'file'), error('cache file exists before the first build: %s', cache_file); end

    t = tic;
    log_fresh = evalc('fresh = mwecmass.driver.build_config(in, [], struct());');
    t_fresh = toc(t);
    if isempty(strfind(log_fresh, 'no cache file, computing')) || ~exist(cache_file, 'file')
      error('%s: first build did not compute and save', modes{m});
    end

    t = tic;
    log_load = evalc('loaded = mwecmass.driver.build_config(in, [], struct());');
    t_load = toc(t);
    if isempty(strfind(log_load, 'Geometry products: reloaded from'))
      error('%s: second build did not reload', modes{m});
    end
    if ~isempty(strfind(log_load, 'Computing per-strip density floors'))
      error('%s: reload still ran the per-strip bounds', modes{m});
    end

    if ~isequaln(rmfield(fresh, 'hydro_cache'), rmfield(loaded, 'hydro_cache'))
      where = first_difference(rmfield(fresh, 'hydro_cache'), rmfield(loaded, 'hydro_cache'));
      if isempty(where), where = '(no field named; isequaln still false)'; end
      error('%s: reloaded config differs from the fresh one at %s', modes{m}, where);
    end

    stored = load(cache_file);
    objects = find_objects(stored, 'file');
    if ~isempty(objects)
      error('%s: cache file holds a non-plain value at %s', modes{m}, strjoin(objects, ', '));
    end
    info = dir(cache_file);
    fprintf('%s: fresh build %.2f s, reload %.2f s, %d config fields identical, cache %.1f kB\n', ...
            modes{m}, t_fresh, t_load, numel(fieldnames(fresh)) - 1, info.bytes / 1024);
  end
end

function restore(saved_path, cache_dir)
  path(saved_path);
  remove_folder(cache_dir);
end

function paths = find_objects(v, path)
  paths = {};
  if isstruct(v)
    names = fieldnames(v);
    for e = 1:numel(v)
      for k = 1:numel(names)
        paths = [paths, find_objects(v(e).(names{k}), [path '.' names{k}])]; %#ok<AGROW>
      end
    end
  elseif iscell(v)
    for k = 1:numel(v)
      paths = [paths, find_objects(v{k}, sprintf('%s{%d}', path, k))]; %#ok<AGROW>
    end
  elseif isobject(v) || isa(v, 'function_handle')
    paths = {path};
  end
end

function remove_folder(folder)
  if exist(folder, 'dir')
    confirm_recursive_rmdir(false, 'local');
    rmdir(folder, 's');
  end
end
