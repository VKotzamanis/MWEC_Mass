function test_geometry_cache_key()
%TEST_GEOMETRY_CACHE_KEY What the geometry-cache key covers.
%   Key sensitivity: the deck bytes, every field of the keyed inputs and every source file in the
%   keyed set change the key; a file outside the set does not; the key does not depend on where
%   the sources sit; a missing source file is an error. Coverage: the cached function receives
%   author inputs only through its argument g, every g field it reads is built by
%   geometry_inputs, and the files its call graph reaches (found by reading the sources) are all
%   in the keyed set. The oracle for the call graph is the text of the sources.
  repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  src_root = fullfile(repo_root, 'src');
  deck = fullfile(repo_root, 'tests', 'standins', 'fixtures', 'cylinder.ms2');
  work = tempname();
  mkdir(work);
  cleanup = onCleanup(@() remove_folder(work));
  cache = @(varargin) mwecmass.internal.geometry_cache(varargin{:});

  inputs = struct('n', 20, 'bounds', [20 2500], 'flag', true, 'name', 'top', 'empty', [], ...
                  'nested', struct('a', 1, 'b', {{'x', 2}}));
  key0 = cache('key', deck, inputs, src_root);
  if numel(key0) ~= 64 || ~strcmp(key0, cache('key', deck, inputs, src_root))
    error('the key is not a stable 64-character hex digest');
  end

  % Inputs: change each field (and a nested one) by the smallest step of its type.
  changed = {'n', inputs.n + eps(inputs.n); 'bounds', [20 2500 + eps(2500)]; 'flag', false; ...
             'name', 'bottom'; 'empty', 0; 'nested', struct('a', 2, 'b', {{'x', 2}})};
  for k = 1:size(changed, 1)
    other = inputs;
    other.(changed{k, 1}) = changed{k, 2};
    if strcmp(cache('key', deck, other, src_root), key0)
      error('changing input %s did not change the key', changed{k, 1});
    end
  end
  other = inputs;
  other.nested.b{2} = 3;
  if strcmp(cache('key', deck, other, src_root), key0)
    error('changing a nested input did not change the key');
  end
  other = rmfield(inputs, 'flag');
  if strcmp(cache('key', deck, other, src_root), key0)
    error('removing an input did not change the key');
  end

  % Field names are framed: names that differ only where a name ends and the value begins.
  if strcmp(cache('key', deck, struct('xu', int8(5)), src_root), ...
            cache('key', deck, struct('x', uint8(5)), src_root))
    error('two different inputs serialise to the same key');
  end

  % Deck bytes.
  deck_copy = fullfile(work, 'cylinder.ms2');
  copyfile(deck, deck_copy);
  if ~strcmp(cache('key', deck_copy, inputs, src_root), key0)
    error('the key depends on the deck path, not only its bytes');
  end
  append_text(deck_copy, ' ');
  if strcmp(cache('key', deck_copy, inputs, src_root), key0)
    error('changing the deck bytes did not change the key');
  end

  % Sources: a copy of src gives the same key; each keyed file changes it; others do not.
  copy_root = fullfile(work, 'src');
  copyfile(src_root, copy_root);
  if ~strcmp(cache('key', deck, inputs, copy_root), key0)
    error('the key depends on where the sources sit');
  end
  files = cache('sources', copy_root);
  pkg = fullfile(copy_root, '+mwecmass');
  relative = cellfun(@(f) strrep(f(numel(pkg) + 2:end), '\', '/'), files, 'UniformOutput', false);
  required = {'+driver/build_config.m', '+driver/build_hydrostatic_tables.m', ...
              '+driver/build_strip_geometry_tables.m', '+driver/parse_hull_deck.m', ...
              '+geometry/MS2Parser.m', '+hydrostatics/compute_strip.m', '+internal/geometry_cache.m'};
  for k = 1:numel(required)
    if ~any(strcmp(relative, required{k}))
      error('keyed set lacks %s', required{k});
    end
  end
  for k = 1:numel(files)
    append_text(files{k}, sprintf('\n%% key test\n'));
    if strcmp(cache('key', deck, inputs, copy_root), key0)
      error('changing %s did not change the key', relative{k});
    end
    restore_file(files{k}, fullfile(src_root, '+mwecmass', relative{k}));
  end
  if ~strcmp(cache('key', deck, inputs, copy_root), key0)
    error('the key did not return after restoring the sources');
  end
  unkeyed = fullfile(pkg, '+optim', 'stage2_bounds.m');
  append_text(unkeyed, sprintf('\n%% key test\n'));
  if ~strcmp(cache('key', deck, inputs, copy_root), key0)
    error('a file outside the keyed set changed the key');
  end
  delete(fullfile(pkg, '+driver', 'parse_hull_deck.m'));
  try
    cache('key', deck, inputs, copy_root);
    error('a missing source file did not raise an error');
  catch err
    if ~strcmp(err.identifier, 'mwecmass:internal:GeometryCacheSourceMissing')
      rethrow(err);
    end
  end
  fprintf('key: %d input changes and a deck change detected, each of %d keyed source files detected, unkeyed file ignored\n', ...
          size(changed, 1) + 2, numel(files));

  % Coverage by reading the sources.
  build_text = strip_comments(fileread(fullfile(src_root, '+mwecmass', '+driver', 'build_config.m')));
  at = strfind(build_text, 'function [products, ms2_model] = compute_geometry_products');
  if numel(at) ~= 1, error('compute_geometry_products not found in build_config.m'); end
  compute_text = build_text(at:end);
  if ~isempty(regexp(compute_text, '\<in\s*[.(]', 'once'))
    error('compute_geometry_products reads the author inputs `in` directly');
  end
  inputs_at = strfind(build_text, 'function g = geometry_inputs');
  if numel(inputs_at) ~= 1, error('geometry_inputs not found in build_config.m'); end
  built = unique(cellfun(@(c) c{1}, regexp(build_text(inputs_at:at - 1), '\<g\.(\w+)\s*=', 'tokens'), ...
                         'UniformOutput', false));
  used = unique(cellfun(@(c) c{1}, regexp(compute_text, '\<g\.(\w+)', 'tokens'), 'UniformOutput', false));
  missing = setdiff(used, built);
  if ~isempty(missing)
    error('compute_geometry_products reads g.%s, which geometry_inputs does not build', strjoin(missing, ', g.'));
  end
  unused = setdiff(built, used);
  if isempty(unused), unused = {'none'}; end
  fprintf('inputs: %d keyed values read, none outside g; keyed but unread: %s\n', ...
          numel(used), strjoin(unused, ', '));

  roots = {compute_text};
  for name = {'build_hydrostatic_tables', 'build_strip_geometry_tables', 'parse_hull_deck'}
    roots{end+1} = strip_comments(fileread(fullfile(src_root, '+mwecmass', '+driver', [name{1} '.m']))); %#ok<AGROW>
  end
  reached = {};
  pending = roots;
  while ~isempty(pending)
    text = pending{end};
    pending(end) = [];
    if ~isempty(regexp(text, '\<(feval|str2func|eval|evalc|evalin|assignin|import|run)\s*[(\s]', 'once'))
      error('a file in the call graph calls by name or imports; the text search cannot follow it');
    end
    refs = regexp(text, 'mwecmass(\.\+?\w+)+', 'match');
    for k = 1:numel(refs)
      file = resolve_reference(src_root, refs{k});
      if isempty(file)
        error('cannot resolve %s to a source file', refs{k});
      end
      if ~any(strcmp(reached, file))
        reached{end+1} = file; %#ok<AGROW>
        pending{end+1} = strip_comments(fileread(file)); %#ok<AGROW>
      end
    end
  end
  keyed = cache('sources', src_root);
  outside = setdiff(reached, keyed);
  if ~isempty(outside)
    error('the cached steps can call files outside the keyed set: %s', strjoin(outside, ', '));
  end
  fprintf('call graph: %d files reachable from the cached steps, all in the keyed set of %d files\n', ...
          numel(reached), numel(keyed));
end

function file = resolve_reference(src_root, ref)
  parts = strsplit(ref, '.');
  parts(1) = [];
  file = '';
  for n = numel(parts):-1:1
    folders = cellfun(@(p) ['+' p], parts(1:n-1), 'UniformOutput', false);
    candidate = fullfile(src_root, '+mwecmass', folders{:}, [parts{n} '.m']);
    if exist(candidate, 'file')
      file = candidate;
      return;
    end
  end
end

function text = strip_comments(text)
% Remove % comments and ... continuations outside character vectors and strings; a quote after
% an identifier character, closing bracket or dot is a transpose.
  lines = strsplit(text, sprintf('\n'), 'CollapseDelimiters', false);
  for k = 1:numel(lines)
    line = lines{k};
    keep = true(1, numel(line));
    quote = '';
    previous = ' ';
    for c = 1:numel(line)
      ch = line(c);
      if ~isempty(quote)
        if ch == quote
          if c < numel(line) && line(c + 1) == quote
            % doubled quote stays inside the string
          else
            quote = '';
          end
        end
      elseif ch == '%' || (ch == '.' && c + 2 <= numel(line) && strcmp(line(c:c+2), '...'))
        keep(c:end) = false;
        break;
      elseif ch == '"'
        quote = '"';
      elseif ch == '''' && ~(isstrprop(previous, 'alphanum') || any(previous == '_)]}''.'))
        quote = '''';
      end
      if ch ~= ' ' && ch ~= sprintf('\t'), previous = ch; end
    end
    lines{k} = line(keep);
  end
  text = strjoin(lines, sprintf('\n'));
end

function append_text(file, text)
  fid = fopen(file, 'a');
  fwrite(fid, text);
  fclose(fid);
end

function restore_file(copy, original)
  copyfile(original, copy);
end

function remove_folder(folder)
  if exist(folder, 'dir')
    confirm_recursive_rmdir(false, 'local');
    rmdir(folder, 's');
  end
end
