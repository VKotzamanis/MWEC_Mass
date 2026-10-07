function varargout = geometry_cache(action, varargin)
%GEOMETRY_CACHE Key, save and reload the geometry products of mwecmass.driver.build_config.
% Syntax:
%   files = mwecmass.internal.geometry_cache('sources', src_root)
%   key   = mwecmass.internal.geometry_cache('key', ms2_file, inputs, src_root)
%   [hit, products, reason] = mwecmass.internal.geometry_cache('load', cache_file, key)
%   ok    = mwecmass.internal.geometry_cache('save', cache_file, key, products)
% src_root is the folder that holds +mwecmass; inputs is the struct of every author input the
% cached steps read. The key is the SHA-256 (hex) of the deck bytes, the serialised inputs, the
% runtime version and the bytes of every source file in 'sources': the whole of +geometry,
% +hydrostatics and +internal, and the driver files that build the products. A missing source
% file is an error. A saved file holds plain arrays and structs only; a containers.Map is stored
% as its keys and values and rebuilt on load, any other object is an error.

    switch action
        case 'sources'
            varargout{1} = source_files(varargin{1});
        case 'key'
            varargout{1} = make_key(varargin{:});
        case 'load'
            [varargout{1:3}] = load_products(varargin{:});
        case 'save'
            varargout{1} = save_products(varargin{:});
        otherwise
            error('mwecmass:internal:GeometryCacheAction', 'Unknown action ''%s''.', action);
    end
end

function files = source_files(src_root)
    pkg = fullfile(src_root, '+mwecmass');
    folders = {'+geometry', '+hydrostatics', '+internal'};
    driver_files = {'build_config.m', 'build_hydrostatic_tables.m', ...
                    'build_strip_geometry_tables.m', 'parse_hull_deck.m'};
    files = {};
    for k = 1:numel(folders)
        folder = fullfile(pkg, folders{k});
        found = list_m_files(folder);
        if isempty(found)
            error('mwecmass:internal:GeometryCacheSourceMissing', ...
                  'No .m source files in %s.', folder);
        end
        files = [files; found]; %#ok<AGROW>
    end
    for k = 1:numel(driver_files)
        file = fullfile(pkg, '+driver', driver_files{k});
        if ~exist(file, 'file')
            error('mwecmass:internal:GeometryCacheSourceMissing', ...
                  'Source file %s is missing.', file);
        end
        files{end+1, 1} = file; %#ok<AGROW>
    end
    files = sort(files);
end

function files = list_m_files(folder)
    files = {};
    if ~exist(folder, 'dir'), return; end
    entries = dir(folder);
    for k = 1:numel(entries)
        name = entries(k).name;
        full = fullfile(folder, name);
        if entries(k).isdir
            if ~any(strcmp(name, {'.', '..'}))
                files = [files; list_m_files(full)]; %#ok<AGROW>
            end
        elseif ~isempty(regexp(name, '\.m$', 'once'))
            files{end+1, 1} = full; %#ok<AGROW>
        end
    end
end

function key = make_key(ms2_file, inputs, src_root)
    files = source_files(src_root);
    pkg_root = fullfile(src_root, '+mwecmass');
    parts = cell(3 + numel(files), 1);
    parts{1} = frame('ms2', read_bytes(ms2_file));
    parts{2} = frame('inputs', serialise(inputs));
    parts{3} = frame('runtime', uint8(version()));
    for k = 1:numel(files)
        rel = strrep(files{k}(numel(pkg_root) + 2:end), '\', '/');
        parts{3 + k} = frame(rel, read_bytes(files{k}));
    end
    stream = vertcat(parts{:});
    digest = java.security.MessageDigest.getInstance('SHA-256');
    digest.update(stream);
    key = lower(sprintf('%02x', typecast(digest.digest(), 'uint8')));
end

function b = frame(name, bytes)
    b = [uint8(name(:)); uint8(0); typecast(uint64(numel(bytes)), 'uint8').'; bytes(:)];
end

function bytes = read_bytes(file)
    fid = fopen(file, 'r');
    if fid < 0
        error('mwecmass:internal:GeometryCacheSourceMissing', 'Cannot read %s.', file);
    end
    bytes = fread(fid, Inf, '*uint8');
    fclose(fid);
end

function b = serialise(v)
% Deterministic bytes of a value made of structs, cells, char, logical and real numeric arrays.
    if isstruct(v)
        names = sort(fieldnames(v));
        parts = {tag('struct', size(v))};
        for k = 1:numel(v)
            for f = 1:numel(names)
                parts{end+1, 1} = uint8(names{f}(:)); %#ok<AGROW>
                parts{end+1, 1} = serialise(v(k).(names{f})); %#ok<AGROW>
            end
        end
        b = vertcat(parts{:});
    elseif iscell(v)
        parts = {tag('cell', size(v))};
        for k = 1:numel(v)
            parts{end+1, 1} = serialise(v{k}); %#ok<AGROW>
        end
        b = vertcat(parts{:});
    elseif ischar(v)
        b = [tag('char', size(v)); numeric_bytes(uint16(v))];
    elseif islogical(v)
        b = [tag('logical', size(v)); numeric_bytes(double(v))];
    elseif isnumeric(v) && isreal(v)
        b = [tag(class(v), size(v)); numeric_bytes(v)];
    else
        error('mwecmass:internal:GeometryCacheKeyType', ...
              'Cannot serialise a value of class %s into the key.', class(v));
    end
end

function b = tag(name, sz)
    b = [uint8(name(:)); uint8(0); typecast(uint64(sz), 'uint8').'];
end

function b = numeric_bytes(v)
    if isempty(v)
        b = zeros(0, 1, 'uint8');
    else
        b = typecast(v(:).', 'uint8').';
    end
end

function [hit, products, reason] = load_products(cache_file, key)
    hit = false;
    products = [];
    if ~exist(cache_file, 'file')
        reason = 'no cache file';
        return;
    end
    try
        stored = load(cache_file);
        valid = isfield(stored, 'key') && isfield(stored, 'plain');
    catch
        valid = false;
    end
    if ~valid
        reason = 'cache file unreadable';
        return;
    end
    if ~strcmp(stored.key, key)
        reason = 'deck, inputs or source changed';
        return;
    end
    products = from_plain(stored.plain);
    hit = true;
    reason = '';
end

function ok = save_products(cache_file, key, products)
    plain = to_plain(products); %#ok<NASGU>
    folder = fileparts(cache_file);
    partial = [cache_file(1:end-4) '_writing.mat'];
    try
        if ~isempty(folder) && ~exist(folder, 'dir'), mkdir(folder); end
        save(partial, '-v7', 'key', 'plain');
        movefile(partial, cache_file, 'f');
        ok = true;
    catch err
        if exist(partial, 'file'), delete(partial); end
        warning('mwecmass:internal:GeometryCacheWrite', ...
                'Geometry cache not saved to %s: %s', cache_file, err.message);
        ok = false;
    end
end

function v = to_plain(v)
    if isa(v, 'containers.Map')
        if ~strcmp(v.KeyType, 'char')
            error('mwecmass:internal:GeometryCacheNotPlain', 'Only char-keyed containers.Map is stored.');
        end
        vals = values(v);
        out.mwecmass_map_keys = keys(v);
        out.mwecmass_map_values = cellfun(@to_plain, vals, 'UniformOutput', false);
        v = out;
    elseif isstruct(v)
        names = fieldnames(v);
        for k = 1:numel(v)
            for f = 1:numel(names)
                v(k).(names{f}) = to_plain(v(k).(names{f}));
            end
        end
    elseif iscell(v)
        v = cellfun(@to_plain, v, 'UniformOutput', false);
    elseif ~(isnumeric(v) || ischar(v) || islogical(v))
        error('mwecmass:internal:GeometryCacheNotPlain', ...
              'Cannot store a value of class %s; products hold plain arrays and structs only.', class(v));
    end
end

function v = from_plain(v)
    if isstruct(v) && isscalar(v) && isequal(sort(fieldnames(v)), ...
            {'mwecmass_map_keys'; 'mwecmass_map_values'})
        m = containers.Map();
        for k = 1:numel(v.mwecmass_map_keys)
            m(v.mwecmass_map_keys{k}) = from_plain(v.mwecmass_map_values{k});
        end
        v = m;
    elseif isstruct(v)
        names = fieldnames(v);
        for k = 1:numel(v)
            for f = 1:numel(names)
                v(k).(names{f}) = from_plain(v(k).(names{f}));
            end
        end
    elseif iscell(v)
        v = cellfun(@from_plain, v, 'UniformOutput', false);
    end
end
