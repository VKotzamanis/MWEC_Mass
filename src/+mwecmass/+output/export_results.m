function export_results(results, final_props, in, path)
%EXPORT_RESULTS Write a data-only results/final_props MAT export.
%   Adds schema metadata and typed-empty fields, converts config objects to
%   plain data, rejects object leaves, and leaves caller inputs unchanged.
    schema = mwecmass.output.export_schema();

    type = in.materials.realisation_type;
    if ~ismember(type, schema.types)
        error('mwec:massSchema:badMode', ...
            'in.materials.realisation_type = ''%s'' is not one of preliminary/thin_shell/modular_precast.', type);
    end

    results.schema_version = schema.schema_version;
    results.realisation_type = type;
    results.context = struct( ...
        'T_heave_goal',  in.context.T_heave_goal, ...
        'T_pitch_goal',  in.context.T_pitch_goal, ...
        'T_surge_goal',  in.context.T_surge_goal, ...
        'T_heave_range', in.context.T_heave_range, ...
        'T_pitch_range', in.context.T_pitch_range, ...
        'T_surge_range', in.context.T_surge_range, ...
        'gm_target',     in.context.gm_target, ...
        'WIS_station',   in.context.WIS_station, ...
        'data_year',     in.context.data_year, ...
        'water_depth',   in.bem.water_depth);

    for k = 1:numel(schema.final_props_fields)
        f = schema.final_props_fields(k);
        if ~isfield(final_props, f.name)
            final_props.(f.name) = mass_typed_empty(f.class);
        end
    end

    if ~isfield(results, 'stage3')
        results.stage3 = [];
    end

    results_out = results;
    if isfield(results_out, 'config') && isstruct(results_out.config) && isscalar(results_out.config)
        results_out.config = local_config_as_data(results_out.config);
    end
    local_assert_no_object_leaf('results', results_out);
    local_assert_no_object_leaf('final_props', final_props);

    results = results_out;
    save(path, 'results', 'final_props');
end

function config = local_config_as_data(config)
    if isfield(config, 'ms2_model')
        config = rmfield(config, 'ms2_model');
    end

    deck = '';
    if isfield(config, 'ms2_file')
        deck = config.ms2_file;
    end
    config.ms2_deck_sha256 = local_file_sha256(deck);
    if isempty(config.ms2_deck_sha256)
        warning('mwec:massExport:deckDigestUnavailable', ...
            ['Could not read the hull-deck file ''%s'' to compute its SHA-256 digest; ' ...
             'results.config.ms2_deck_sha256 is written empty.'], deck);
    end

    if isfield(config, 'boundary_cache') && isstruct(config.boundary_cache) ...
            && isscalar(config.boundary_cache) && isfield(config.boundary_cache, 'data') ...
            && isa(config.boundary_cache.data, 'containers.Map')
        map = config.boundary_cache.data;
        map_keys = keys(map);
        map_values = values(map);
        if ~isempty(map_keys) && all(cellfun(@isvarname, map_keys))
            data = struct();
            for i = 1:numel(map_keys)
                data.(map_keys{i}) = map_values{i};
            end
        else
            data = struct('keys', {{map_keys}}, 'values', {{map_values}});
        end
        config.boundary_cache.data = data;
    end
end

function hex = local_file_sha256(file)
    hex = '';
    if isempty(file) || ~(ischar(file) || isstring(file)) || ~isfile(char(file))
        return;
    end
    fid = fopen(char(file), 'r');
    if fid < 0
        return;
    end
    bytes = fread(fid, Inf, '*uint8');
    fclose(fid);
    digest = java.security.MessageDigest.getInstance('SHA-256');
    if ~isempty(bytes)
        digest.update(bytes);
    end
    hex = lower(sprintf('%02x', typecast(digest.digest(), 'uint8')));
end

function local_assert_no_object_leaf(pth, value)
    if isobject(value) && ~isstring(value)
        error('mwec:massExport:objectLeaf', ...
            ['Export leaf %s is an instance of class %s. The export stores data only, so that ' ...
             'moving a class between namespaces can never make a saved result file unreadable.'], ...
            pth, class(value));
    end
    if isstruct(value)
        names = fieldnames(value);
        for e = 1:numel(value)
            if isscalar(value)
                prefix = pth;
            else
                prefix = sprintf('%s(%d)', pth, e);
            end
            for i = 1:numel(names)
                local_assert_no_object_leaf([prefix '.' names{i}], value(e).(names{i}));
            end
        end
    elseif iscell(value)
        for i = 1:numel(value)
            local_assert_no_object_leaf(sprintf('%s{%d}', pth, i), value{i});
        end
    end
end

function v = mass_typed_empty(class_name)
    switch class_name
        case 'double'
            v = double.empty(0, 0);
        case 'logical'
            v = logical.empty(0, 0);
        case 'char'
            v = '';
        case 'cell'
            v = {};
        case 'struct'
            v = struct();
        otherwise
            error('mwec:massSchema:badClass', 'mass_typed_empty: unsupported class "%s".', class_name);
    end
end
