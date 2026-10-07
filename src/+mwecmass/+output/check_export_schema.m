function check_export_schema(file, expected_type)
%CHECK_EXPORT_SCHEMA Validate a results/final_props MAT export.
%   Optional expected_type selects type-scoped invariants; otherwise the recorded type is used.
%   Throws on missing fields, bad classes/shapes, inconsistent typed-empty fields, or object leaves.
    if nargin < 2
        expected_type = '';
    end

    schema = mwecmass.output.export_schema();
    L = load(file, 'results', 'final_props');
    if isfield(L, 'results') && isfield(L.results, 'mode') && ~isfield(L.results, 'realisation_type')
        L.results.realisation_type = L.results.mode;
    end
    if isfield(L, 'final_props') && isfield(L.final_props, 'realisation_mode') && ~isfield(L.final_props, 'fill_method')
        L.final_props.fill_method = L.final_props.realisation_mode;
    end

    if ~isfield(L, 'final_props')
        error('mwec:massSchema:missingField', 'Mass export is missing required field %s.', 'final_props');
    end
    if ~isfield(L, 'results')
        error('mwec:massSchema:missingField', 'Mass export is missing required field %s.', 'results');
    end

    for k = 1:numel(schema.results_top_fields)
        name = schema.results_top_fields{k};
        if ~isfield(L.results, name)
            error('mwec:massSchema:missingField', 'Mass export is missing required field %s.', ['results.' name]);
        end
    end

    for k = 1:numel(schema.final_props_fields)
        name = schema.final_props_fields(k).name;
        if ~isfield(L.final_props, name)
            error('mwec:massSchema:missingField', 'Mass export is missing required field %s.', ['final_props.' name]);
        end
    end

    type = L.results.realisation_type;
    if ~ismember(type, schema.types)
        error('mwec:massSchema:badMode', ...
            'results.realisation_type = ''%s'' is not one of preliminary/thin_shell/modular_precast.', type);
    end

    if ~strcmp(L.results.schema_version, schema.schema_version)
        error('mwec:massSchema:badVersion', ...
            'results.schema_version = ''%s'', expected ''%s''.', L.results.schema_version, schema.schema_version);
    end

    for k = 1:numel(schema.fixed_shape_checks)
        f = schema.fixed_shape_checks(k);
        container = L.final_props;   % all 8 spot-check fields live on final_props
        if ~isfield(container, f.name)
            error('mwec:massSchema:missingField', 'Mass export is missing required field %s.', ['final_props.' f.name]);
        end
        val = container.(f.name);
        if ~isa(val, f.class)
            error('mwec:massSchema:badClass', 'Field %s has class %s, expected %s.', ...
                ['final_props.' f.name], class(val), f.class);
        end
        if ~isequal(size(val), f.size)
            error('mwec:massSchema:badSize', 'Field %s has size %s, expected %s.', ...
                ['final_props.' f.name], mat2str(size(val)), mat2str(f.size));
        end
    end

    check_type = type;
    if ~isempty(expected_type)
        check_type = expected_type;
    end
    is_preliminary = strcmp(check_type, 'preliminary');

    for k = 1:numel(schema.final_props_fields)
        f = schema.final_props_fields(k);
        if ~f.realisation_only
            continue;
        end
        is_empty_val = local_is_typed_empty(L.final_props.(f.name));
        if is_preliminary && ~is_empty_val
            error('mwec:massSchema:notEmptyOfType', ...
                'Field %s is not empty-of-type under realisation type ''%s''.', ['final_props.' f.name], check_type);
        elseif ~is_preliminary && is_empty_val
            error('mwec:massSchema:unexpectedlyEmpty', ...
                'Field %s is empty under realisation type ''%s'', where it should be populated.', ['final_props.' f.name], check_type);
        end
    end

    st = L.results.stage3;
    if is_preliminary
        if ~local_is_typed_empty(st)
            error('mwec:massSchema:notEmptyOfType', ...
                'Field %s is not empty-of-type under realisation type ''%s''.', 'results.stage3', check_type);
        end
    else
        if ~isstruct(st) || ~isscalar(st)
            error('mwec:massSchema:unexpectedlyEmpty', ...
                'Field %s is empty under realisation type ''%s'', where it should be populated.', 'results.stage3', check_type);
        end
        for k = 1:numel(schema.stage3_fields)
            f = schema.stage3_fields(k);
            if ~isfield(st, f.name)
                error('mwec:massSchema:missingField', 'Mass export is missing required field %s.', ['results.stage3.' f.name]);
            end
            if ~isempty(f.class) && ~isa(st.(f.name), f.class)
                error('mwec:massSchema:badClass', 'Field %s has class %s, expected %s.', ...
                    ['results.stage3.' f.name], class(st.(f.name)), f.class);
            end
        end
        if ~strcmp(st.mode, check_type) || ~any(strcmp(st.status, {'accepted', 'failed'}))
            error('mwec:massSchema:badStage3', ...
                'results.stage3 has mode ''%s'' and status ''%s''; expected mode ''%s'' and status accepted or failed.', ...
                st.mode, st.status, check_type);
        end
    end

    if ~isfield(L.results, 'config')
        error('mwec:massSchema:missingField', 'Mass export is missing required field %s.', 'results.config');
    end
    for k = 1:numel(schema.config_required_subfields)
        name = schema.config_required_subfields{k};
        if ~isfield(L.results.config, name)
            error('mwec:massSchema:missingField', 'Mass export is missing required field %s.', ['results.config.' name]);
        end
    end

    for k = 1:numel(schema.config_required_char_subfields)
        name = schema.config_required_char_subfields{k};
        if ~isfield(L.results.config, name)
            error('mwec:massSchema:missingField', 'Mass export is missing required field %s.', ['results.config.' name]);
        end
        val = L.results.config.(name);
        if ~ischar(val)
            error('mwec:massSchema:badClass', 'Field %s has class %s, expected %s.', ...
                ['results.config.' name], class(val), 'char');
        end
    end

    local_assert_no_object_leaf('results', L.results);
    local_assert_no_object_leaf('final_props', L.final_props);
end

function local_assert_no_object_leaf(pth, value)
    if isobject(value) && ~isstring(value)
        error('mwec:massSchema:objectLeaf', ...
            'Field %s is an instance of class %s; the export stores data only.', pth, class(value));
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

function tf = local_is_typed_empty(val)
    if isstruct(val)
        tf = isempty(val) || isempty(fieldnames(val));
    else
        tf = isempty(val);
    end
end
