classdef MS2Parser
%MS2PARSER Parse MultiSurf .ms2 geometry into evaluatable entities.
% Geometry remains parametric: points are [1x3], curve evaluations are [Nx3],
% and surface evaluations use normalized parameters u,v in [0,1].
% Each entity's references are resolved to entity numbers once, at parse time, together with the
% quantities that do not depend on the parameters; B-spline evaluation is self-contained.

    properties
        entities        % containers.Map: name → entity struct
        visible_surfs   % cell array of visible surface names
        units           % string: 'm', 'ft'
        extents         % [xmin ymin zmin xmax ymax zmax]
        filename        % source .ms2 file path
        file_symmetry   % cell array: e.g. {'x','y'} from header
    end

    properties (Access = private)
        names           % cell array: entity names; the position is the entity number
        ent             % cell array: entity records (see resolve_entities)
        ids             % containers.Map: entity name → entity number
        id_empty        % number of the placeholder that stands for an empty name in a list; 0 if none
    end


    methods (Static)

        function model = parse(filename)
        % PARSE: Read a MultiSurf file and build its entity graph.

            if ~exist(filename, 'file')
                error('mwecmass:geometry:FileNotFound', ...
                       'MS2 file not found: %s', filename);
            end

            model = mwecmass.geometry.MS2Parser();
            model.filename = filename;
            model.entities = containers.Map();

            fid = fopen(filename, 'r');
            raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
            fclose(fid);
            lines = raw{1};

            model.units   = 'm';
            model.extents = zeros(1, 6);
            model.file_symmetry = {};    % from header 'Symmetry: x y'
            for i = 1:length(lines)
                L = strtrim(lines{i});
                if startsWith(L, 'Units:')
                    tokens = strsplit(L);
                    if length(tokens) >= 2, model.units = tokens{2}; end
                elseif startsWith(L, 'Extents:')
                    nums = sscanf(L, 'Extents: %f %f %f %f %f %f');
                    if length(nums) == 6, model.extents = nums'; end
                elseif startsWith(L, 'Symmetry:')
                    sym_tokens = strsplit(strtrim(L(10:end)));
                    model.file_symmetry = lower(sym_tokens);
                elseif strcmp(L, 'BeginModel;')
                    break;
                end
            end

            in_model  = false;
            line_buf  = '';

            for i = 1:length(lines)
                L = strtrim(lines{i});

                if strcmp(L, 'BeginModel;')
                    in_model = true;
                    continue;
                end
                if strcmp(L, 'EndModel;')
                    break;
                end
                if ~in_model, continue; end

                if isempty(L), continue; end
                if startsWith(L, 'Attribute:'), continue; end

                if ~isempty(line_buf)
                    if startsWith(lines{i}, ' ') || startsWith(lines{i}, char(9)) ...
                            || startsWith(L, '{') || startsWith(L, 'A:')
                        line_buf = [line_buf ' ' L]; %#ok<AGROW> -- continuation lines
                        if endsWith(L, ';')
                            model.parse_entity_line(line_buf);
                            line_buf = '';
                        end
                        continue;
                    else
                        if ~isempty(line_buf)
                            model.parse_entity_line(line_buf);
                        end
                        line_buf = '';
                    end
                end

                if endsWith(L, ';')
                    model.parse_entity_line(L);
                else
                    line_buf = L;
                end
            end

            if ~isempty(line_buf)
                model.parse_entity_line(line_buf);
            end

            model.visible_surfs = {};
            keys = model.entities.keys();
            for i = 1:length(keys)
                e = model.entities(keys{i});
                if e.is_surface && e.visibility > 0
                    model.visible_surfs{end+1} = keys{i};
                end
            end


            has_explicit_mirrors = false;
            for i = 1:length(model.visible_surfs)
                e_chk = model.entities(model.visible_surfs{i});
                if strcmp(e_chk.type, 'MirrSurf')
                    has_explicit_mirrors = true;
                    break;
                end
            end

            if ~has_explicit_mirrors && ~isempty(model.file_symmetry)
                orig_sources = model.visible_surfs;
                n_synth = 0;

                if any(strcmp(model.file_symmetry, 'x'))
                    for s = 1:length(orig_sources)
                        mname = [orig_sources{s} '_mirrX'];
                        em = struct();
                        em.type       = 'MirrSurf';
                        em.visibility = 1;
                        em.is_surface = true;
                        em.params     = struct('source', orig_sources{s}, ...
                                               'mirror_plane', 'X');
                        em.parents    = orig_sources(s);
                        model.entities(mname) = em;
                        model.visible_surfs{end+1} = mname;
                        n_synth = n_synth + 1;
                    end
                end

                if any(strcmp(model.file_symmetry, 'y'))
                    surfs_to_mirror = model.visible_surfs;
                    for s = 1:length(surfs_to_mirror)
                        mname = [surfs_to_mirror{s} '_mirrY'];
                        em = struct();
                        em.type       = 'MirrSurf';
                        em.visibility = 1;
                        em.is_surface = true;
                        em.params     = struct('source', surfs_to_mirror{s}, ...
                                               'mirror_plane', 'Y');
                        em.parents    = surfs_to_mirror(s);
                        model.entities(mname) = em;
                        model.visible_surfs{end+1} = mname;
                        n_synth = n_synth + 1;
                    end
                end

                if n_synth > 0
                    fprintf('    File symmetry [%s]: created %d synthetic mirrors\n', ...
                            strjoin(model.file_symmetry, ','), n_synth);
                end
            end

            model = model.resolve_entities();

            fprintf('  MS2 Parser: %s\n', filename);
            fprintf('    Entities: %d total, %d visible surfaces\n', ...
                    model.entities.Count, length(model.visible_surfs));
            for i = 1:length(model.visible_surfs)
                e = model.entities(model.visible_surfs{i});
                fprintf('      %s (%s)\n', model.visible_surfs{i}, e.type);
            end
        end

    end % methods (Static)



    methods (Access = private)

        function parse_entity_line(obj, line)
        % PARSE_ENTITY_LINE: Parse one entity definition into a normalized struct.

            line = strtrim(regexprep(line, ';$', ''));
            if isempty(line), return; end

            tokens = strsplit(line);
            if isempty(tokens), return; end

            etype = tokens{1};

            skip_types = {'RealList', 'Variable', 'Pathname'};
            if any(strcmp(etype, skip_types)), return; end

            try
                switch etype
                    case 'FramePoint'
                        obj.parse_frame_point(tokens, line);
                    case 'BCurve'
                        obj.parse_bcurve(tokens, line);
                    case 'Conic'
                        obj.parse_conic(tokens, line);
                    case 'CopyCurve'
                        obj.parse_copy_curve(tokens, line);
                    case 'Line'
                        obj.parse_line_entity(tokens, line);
                    case 'BSubCurve'
                        obj.parse_bsub_curve(tokens, line);
                    case 'MirrPoint'
                        obj.parse_mirr_point(tokens, line);
                    case 'Arc'
                        obj.parse_arc(tokens, line);
                    case 'PolyCurve2'
                        obj.parse_polycurve2(tokens, line);
                    case 'ProjCurve'
                        obj.parse_proj_curve(tokens, line);
                    case 'RuledSurf'
                        obj.parse_ruled_surf(tokens, line);
                    case 'RevSurf'
                        obj.parse_rev_surf(tokens, line);
                    case 'BLoftSurf'
                        obj.parse_bloft_surf(tokens, line);
                    case 'DevSurf'
                        obj.parse_dev_surf(tokens, line);
                    case 'MirrSurf'
                        obj.parse_mirr_surf(tokens, line);
                    case 'EdgeSnake'
                        obj.parse_edge_snake(tokens, line);
                    case 'AbsBead'
                        obj.parse_abs_bead(tokens, line);
                    case 'AbsRing'
                        obj.parse_abs_ring(tokens, line);
                    case 'BSubSnake'
                        obj.parse_bsub_snake(tokens, line);
                end
            catch ME
                warning('mwecmass:geometry:ParseError', ...
                        'Failed to parse: %s\n  Error: %s', line, ME.message);
            end
        end


        function parse_frame_point(obj, tokens, ~)
        % PARSE_FRAME_POINT: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            if isempty(slash_idx), return; end

            after_slash = tokens(slash_idx+1 : end);
            if length(after_slash) < 6, return; end

            parent1 = after_slash{2};
            parent2 = after_slash{3};
            dx = str2double(after_slash{4});
            dy = str2double(after_slash{5});
            dz = str2double(after_slash{6});

            if strcmp(parent1, '*'), parent1 = ''; end
            if strcmp(parent2, '*'), parent2 = ''; end

            e = struct();
            e.type       = 'FramePoint';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('parent1', parent1, 'parent2', parent2, ...
                                  'offset', [dx, dy, dz]);
            e.parents    = {};
            if ~isempty(parent1), e.parents{end+1} = parent1; end
            if ~isempty(parent2), e.parents{end+1} = parent2; end

            obj.entities(name) = e;
        end

        function parse_mirr_point(obj, tokens, ~)
        % PARSE_MIRR_POINT: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            source_pt  = after_slash{1};
            plane_str  = after_slash{2};   % e.g. '*X=0'
            mirror_axis = upper(plane_str(2));  % 'X', 'Y', or 'Z'

            e = struct();
            e.type       = 'MirrPoint';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('source', source_pt, 'plane', mirror_axis);
            e.parents    = {source_pt};

            obj.entities(name) = e;
        end


        function parse_bcurve(obj, tokens, line)
        % PARSE_BCURVE: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);
            degree = str2double(after_slash{2});

            brace_start = strfind(line, '{');
            brace_end   = strfind(line, '}');
            if ~isempty(brace_start) && ~isempty(brace_end)
                cp_str  = line(brace_start(1)+1 : brace_end(1)-1);
                cp_names = strsplit(strtrim(cp_str));
            else
                cp_names = {};
            end

            e = struct();
            e.type       = 'BCurve';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('degree', degree, 'ctrl_pt_names', {cp_names});
            e.parents    = cp_names;

            obj.entities(name) = e;
        end

        function parse_conic(obj, tokens, ~)
        % PARSE_CONIC: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            conic_type = str2double(after_slash{2});
            center_name = after_slash{3};
            radius_name = after_slash{4};
            apex_name   = after_slash{5};
            angle_start = str2double(after_slash{6});
            angle_end   = str2double(after_slash{7});

            e = struct();
            e.type       = 'Conic';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('conic_type', conic_type, ...
                                  'center', center_name, ...
                                  'radius_pt', radius_name, ...
                                  'apex_pt', apex_name, ...
                                  'angle_start', angle_start, ...
                                  'angle_end', angle_end);
            e.parents    = {center_name, radius_name, apex_name};

            obj.entities(name) = e;
        end

        function parse_copy_curve(obj, tokens, ~)
        % PARSE_COPY_CURVE: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            source_curve = after_slash{2};
            src_pt       = after_slash{3};
            dst_pt       = after_slash{4};
            sx = str2double(after_slash{5});
            sy = str2double(after_slash{6});
            sz = str2double(after_slash{7});

            e = struct();
            e.type       = 'CopyCurve';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('source', source_curve, ...
                                  'src_pt', src_pt, 'dst_pt', dst_pt, ...
                                  'scale', [sx, sy, sz]);
            e.parents    = {source_curve, src_pt, dst_pt};

            obj.entities(name) = e;
        end

        function parse_line_entity(obj, tokens, ~)
        % PARSE_LINE_ENTITY: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            pt_start = after_slash{2};
            pt_end   = after_slash{3};

            e = struct();
            e.type       = 'Line';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('pt_start', pt_start, 'pt_end', pt_end);
            e.parents    = {pt_start, pt_end};

            obj.entities(name) = e;
        end

        function parse_bsub_curve(obj, tokens, line)
        % PARSE_BSUB_CURVE: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);
            degree = str2double(after_slash{2});

            brace_start = strfind(line, '{');
            brace_end   = strfind(line, '}');
            if ~isempty(brace_start) && ~isempty(brace_end)
                bp_str  = line(brace_start(1)+1 : brace_end(1)-1);
                bead_names = strsplit(strtrim(bp_str));
            else
                bead_names = {};
            end

            e = struct();
            e.type       = 'BSubCurve';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('degree', degree, 'bead_names', {bead_names});
            e.parents    = bead_names;

            obj.entities(name) = e;
        end

        function parse_arc(obj, tokens, ~)
        % PARSE_ARC: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            arc_type = str2double(after_slash{2});
            pt1_name = after_slash{3};
            pt2_name = after_slash{4};
            pt3_name = after_slash{5};

            e = struct();
            e.type       = 'Arc';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('arc_type', arc_type, ...
                                  'pt_start', pt1_name, ...
                                  'pt_center', pt2_name, ...
                                  'pt_end', pt3_name);
            e.parents    = {pt1_name, pt2_name, pt3_name};

            obj.entities(name) = e;
        end

        function parse_polycurve2(obj, tokens, line)
        % PARSE_POLYCURVE2: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            brace_start = strfind(line, '{');
            brace_end   = strfind(line, '}');
            if ~isempty(brace_start) && ~isempty(brace_end)
                inner  = line(brace_start(1)+1 : brace_end(1)-1);
                curve_names = strsplit(strtrim(inner));
            else
                curve_names = {};
            end

            e = struct();
            e.type       = 'PolyCurve2';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('curve_names', {curve_names});
            e.parents    = curve_names;

            obj.entities(name) = e;
        end

        function parse_proj_curve(obj, tokens, ~)
        % PARSE_PROJ_CURVE: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            source_name = after_slash{2};
            plane_str   = after_slash{3};  % e.g. '*Y=0'
            proj_axis   = upper(plane_str(2));  % 'X', 'Y', or 'Z'

            e = struct();
            e.type       = 'ProjCurve';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('source', source_name, 'proj_plane', proj_axis);
            e.parents    = {source_name};

            obj.entities(name) = e;
        end


        function parse_edge_snake(obj, tokens, ~)
        % PARSE_EDGE_SNAKE: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            edge_index   = str2double(after_slash{2});
            surface_name = after_slash{3};

            e = struct();
            e.type       = 'EdgeSnake';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('edge_index', edge_index, ...
                                  'surface_name', surface_name);
            e.parents    = {surface_name};

            obj.entities(name) = e;
        end

        function parse_abs_bead(obj, tokens, ~)
        % PARSE_ABS_BEAD: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            parent_name = after_slash{1};
            param_val   = str2double(after_slash{2});

            e = struct();
            e.type       = 'AbsBead';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('parent_curve', parent_name, ...
                                  'parameter', param_val);
            e.parents    = {parent_name};

            obj.entities(name) = e;
        end

        function parse_abs_ring(obj, tokens, ~)
        % PARSE_ABS_RING: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            parent_name = after_slash{1};
            param_val   = str2double(after_slash{2});

            e = struct();
            e.type       = 'AbsRing';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('parent_snake', parent_name, ...
                                  'parameter', param_val);
            e.parents    = {parent_name};

            obj.entities(name) = e;
        end

        function parse_bsub_snake(obj, tokens, line)
        % PARSE_BSUB_SNAKE: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);
            degree = str2double(after_slash{2});

            brace_start = strfind(line, '{');
            brace_end   = strfind(line, '}');
            if ~isempty(brace_start) && ~isempty(brace_end)
                bp_str  = line(brace_start(1)+1 : brace_end(1)-1);
                bead_names = strsplit(strtrim(bp_str));
            else
                bead_names = {};
            end

            e = struct();
            e.type       = 'BSubSnake';
            e.visibility = vis;
            e.is_surface = false;
            e.params     = struct('degree', degree, 'bead_names', {bead_names});
            e.parents    = bead_names;

            obj.entities(name) = e;
        end


        function parse_ruled_surf(obj, tokens, ~)
        % PARSE_RULED_SURF: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            curve1_name = after_slash{2};
            curve2_name = after_slash{3};

            e = struct();
            e.type       = 'RuledSurf';
            e.visibility = vis;
            e.is_surface = true;
            e.params     = struct('curve1', curve1_name, 'curve2', curve2_name);
            e.parents    = {curve1_name, curve2_name};

            obj.entities(name) = e;
        end

        function parse_rev_surf(obj, tokens, ~)
        % PARSE_REV_SURF: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            profile_name = after_slash{2};
            axis_name    = after_slash{3};
            angle_start  = str2double(after_slash{4});
            angle_end    = str2double(after_slash{5});

            e = struct();
            e.type       = 'RevSurf';
            e.visibility = vis;
            e.is_surface = true;
            e.params     = struct('profile', profile_name, ...
                                  'axis', axis_name, ...
                                  'angle_start', angle_start, ...
                                  'angle_end', angle_end);
            e.parents    = {profile_name, axis_name};

            obj.entities(name) = e;
        end

        function parse_bloft_surf(obj, tokens, line)
        % PARSE_BLOFT_SURF: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);
            degree = str2double(after_slash{2});

            brace_start = strfind(line, '{');
            brace_end   = strfind(line, '}');
            if ~isempty(brace_start) && ~isempty(brace_end)
                sec_str = line(brace_start(1)+1 : brace_end(1)-1);
                section_names = strsplit(strtrim(sec_str));
            else
                section_names = {};
            end

            e = struct();
            e.type       = 'BLoftSurf';
            e.visibility = vis;
            e.is_surface = true;
            e.params     = struct('degree', degree, ...
                                  'section_names', {section_names});
            e.parents    = section_names;

            obj.entities(name) = e;
        end

        function parse_dev_surf(obj, tokens, ~)
        % PARSE_DEV_SURF: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            snake_name = after_slash{1};
            curve_name = after_slash{2};

            e = struct();
            e.type       = 'DevSurf';
            e.visibility = vis;
            e.is_surface = true;
            e.params     = struct('snake', snake_name, 'curve', curve_name);
            e.parents    = {snake_name, curve_name};

            obj.entities(name) = e;
        end

        function parse_mirr_surf(obj, tokens, ~)
        % PARSE_MIRR_SURF: Parse one entity definition into a normalized struct.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            source_name = after_slash{1};
            plane_token = after_slash{2};  % e.g. '*Y=0' or '*X=0'

            if contains(plane_token, 'Y=0')
                mirror_plane = 'Y';
            elseif contains(plane_token, 'X=0')
                mirror_plane = 'X';
            else
                mirror_plane = 'Y';  % default
                warning('mwecmass:geometry:UnknownPlane', ...
                        'Unknown mirror plane: %s, defaulting to Y=0', plane_token);
            end

            e = struct();
            e.type       = 'MirrSurf';
            e.visibility = vis;
            e.is_surface = true;
            e.params     = struct('source', source_name, ...
                                  'mirror_plane', mirror_plane);
            e.parents    = {source_name};

            obj.entities(name) = e;
        end

    end % methods (Access = private)



    methods (Access = private)

        function obj = resolve_entities(obj)
        % RESOLVE_ENTITIES: Number the entities, replace every name reference by its number and
        % evaluate every parameter-independent quantity (point coordinates, control points,
        % knots, revolution axes) once. A name that is referenced but never defined gets a
        % placeholder entry of type 'Missing', so evaluating it raises the error the lookup raised.

            names = obj.entities.keys();
            n_def = numel(names);
            ids = containers.Map();
            for i = 1:n_def
                ids(names{i}) = i;
            end

            ent = cell(1, n_def);
            id_empty = 0;
            for i = 1:n_def
                r = obj.entities(names{i});
                r.name = names{i};
                r.k_ok = false;
                for k = 1:numel(r.parents)
                    if isempty(r.parents{k})
                        if id_empty == 0
                            names{end+1} = ''; %#ok<AGROW>
                            id_empty = numel(names);
                            ent{end+1} = struct('type', 'Missing', 'name', '', ...
                                                'is_surface', false, 'k_ok', true); %#ok<AGROW>
                        end
                    elseif ~ids.isKey(r.parents{k})
                        names{end+1} = r.parents{k}; %#ok<AGROW>
                        ids(r.parents{k}) = numel(names);
                        ent{end+1} = struct('type', 'Missing', 'name', r.parents{k}, ...
                                            'is_surface', false, 'k_ok', true); %#ok<AGROW>
                    end
                end
                ent{i} = r;
            end

            is_point = false(1, n_def);
            for i = 1:n_def
                ent{i} = obj.resolve_refs(ent{i}, ids, id_empty);
                is_point(i) = any(strcmp(ent{i}.type, {'FramePoint', 'MirrPoint', 'AbsBead', 'AbsRing'}));
            end

            obj.names = names;
            obj.ent = ent;
            obj.ids = ids;
            obj.id_empty = id_empty;

            % Points first, so that the curves built on them find their coordinates ready. A
            % constant that cannot be evaluated stays unmarked; the evaluator then repeats the
            % evaluation and raises the same error.
            for pass = 1:2
                for i = 1:n_def
                    if is_point(i) ~= (pass == 1), continue; end
                    try
                        r = obj.compute_constants(obj.ent{i});
                        r.k_ok = true;
                        obj.ent{i} = r;
                    catch
                    end
                end
            end
        end


        function r = resolve_refs(obj, r, ids, id_empty)
        % RESOLVE_REFS: Replace the name references of entity record r by entity numbers. An
        % empty name inside a list stands for the 'Missing' placeholder id_empty, so evaluation
        % raises the error a lookup of that name raised; the optional parents of a FramePoint
        % keep 0 for an empty name.

            p = r.params;
            switch r.type
                case 'FramePoint'
                    r.parent1_id = obj.ref_ids(ids, p.parent1, 0);
                    r.parent2_id = obj.ref_ids(ids, p.parent2, 0);
                case 'MirrPoint'
                    r.source_id = obj.ref_ids(ids, p.source, id_empty);
                case 'AbsBead'
                    r.parent_id = obj.ref_ids(ids, p.parent_curve, id_empty);
                case 'AbsRing'
                    r.parent_id = obj.ref_ids(ids, p.parent_snake, id_empty);
                case 'BCurve'
                    r.cp_ids = obj.ref_ids(ids, p.ctrl_pt_names, id_empty);
                case 'Conic'
                    r.center_id = obj.ref_ids(ids, p.center, id_empty);
                    r.radius_id = obj.ref_ids(ids, p.radius_pt, id_empty);
                    r.apex_id   = obj.ref_ids(ids, p.apex_pt, id_empty);
                case 'CopyCurve'
                    r.source_id = obj.ref_ids(ids, p.source, id_empty);
                    r.src_pt_id = obj.ref_ids(ids, p.src_pt, id_empty);
                    r.dst_pt_id = obj.ref_ids(ids, p.dst_pt, id_empty);
                case 'Line'
                    r.start_id = obj.ref_ids(ids, p.pt_start, id_empty);
                    r.end_id   = obj.ref_ids(ids, p.pt_end, id_empty);
                case 'Arc'
                    r.start_id  = obj.ref_ids(ids, p.pt_start, id_empty);
                    r.centre_id = obj.ref_ids(ids, p.pt_center, id_empty);
                    r.end_id    = obj.ref_ids(ids, p.pt_end, id_empty);
                case 'PolyCurve2'
                    r.curve_ids = obj.ref_ids(ids, p.curve_names, id_empty);
                case 'ProjCurve'
                    r.source_id = obj.ref_ids(ids, p.source, id_empty);
                case 'EdgeSnake'
                    r.surface_id = obj.ref_ids(ids, p.surface_name, id_empty);
                case 'RuledSurf'
                    r.curve1_id = obj.ref_ids(ids, p.curve1, id_empty);
                    r.curve2_id = obj.ref_ids(ids, p.curve2, id_empty);
                case 'RevSurf'
                    r.profile_id = obj.ref_ids(ids, p.profile, id_empty);
                    r.axis_id    = obj.ref_ids(ids, p.axis, id_empty);
                case 'BLoftSurf'
                    r.section_ids = obj.ref_ids(ids, p.section_names, id_empty);
                    r.surf_key    = strjoin(p.section_names, '|');
                case 'DevSurf'
                    r.snake_id = obj.ref_ids(ids, p.snake, id_empty);
                    r.curve_id = obj.ref_ids(ids, p.curve, id_empty);
                case 'MirrSurf'
                    r.source_id = obj.ref_ids(ids, p.source, id_empty);
            end
        end


        function r = record_of(obj, e)
        % RECORD_OF: The record of an entity struct as model.entities holds it, resolved on entry
        % for the public evaluators named after an entity type.

            r = e;
            r.name = '';
            r.k_ok = false;
            r = obj.resolve_refs(r, obj.ids, obj.id_empty);
        end


        function r = compute_constants(obj, r)
        % COMPUTE_CONSTANTS: Evaluate what entity r needs that does not depend on the parameters.

            p = r.params;
            switch r.type
                case 'FramePoint'
                    if r.parent1_id > 0
                        base = obj.any_point_at(r.parent1_id);
                    elseif r.parent2_id > 0
                        base = obj.any_point_at(r.parent2_id);
                    else
                        base = [0, 0, 0];  % absolute position
                    end
                    r.pt = base + p.offset;

                case 'MirrPoint'
                    pt = obj.any_point_at(r.source_id);
                    switch p.plane
                        case 'X', pt(1) = -pt(1);
                        case 'Y', pt(2) = -pt(2);
                        case 'Z', pt(3) = -pt(3);
                    end
                    r.pt = pt;

                case {'AbsBead', 'AbsRing'}
                    try
                        pt_arr = obj.snake_at(r.parent_id, p.parameter);
                    catch
                        pt_arr = obj.curve_at(r.parent_id, p.parameter);
                    end
                    r.pt = pt_arr(1, :);  % ensure [1×3]

                case 'BCurve'
                    n_cp = numel(r.cp_ids);
                    ctrl_pts = zeros(n_cp, 3);
                    for i = 1:n_cp
                        ctrl_pts(i, :) = obj.point_at(r.cp_ids(i));
                    end
                    r.ctrl_pts = ctrl_pts;
                    r.knots = mwecmass.geometry.MS2Parser.make_clamped_knots(n_cp, p.degree);

                case 'Conic'
                    r.center = obj.any_point_at(r.center_id);
                    r.rad_pt = obj.any_point_at(r.radius_id);
                    r.apex   = obj.any_point_at(r.apex_id);

                case 'CopyCurve'
                    r.src_pos = obj.point_at(r.src_pt_id);
                    r.dst_pos = obj.point_at(r.dst_pt_id);

                case 'Line'
                    r.p_start = obj.any_point_at(r.start_id);
                    r.p_end   = obj.any_point_at(r.end_id);

                case 'Arc'
                    r.p_start  = obj.any_point_at(r.start_id);
                    r.p_centre = obj.any_point_at(r.centre_id);
                    r.p_end    = obj.any_point_at(r.end_id);

                case {'BSubCurve', 'BSubSnake'}
                    bead_names = p.bead_names;
                    bead1 = obj.entities(bead_names{1});
                    bead2 = obj.entities(bead_names{end});
                    r.t_start = bead1.params.parameter;
                    r.t_end   = bead2.params.parameter;
                    r.parent_id = obj.id_of(bead1.params.parent_curve);

                case 'RevSurf'
                    axis_ent = obj.ent_unchecked(r.axis_id);
                    r.axis_start = obj.any_point_at(obj.id_of(axis_ent.params.pt_start));
                    axis_end     = obj.any_point_at(obj.id_of(axis_ent.params.pt_end));
                    axis_vec = axis_end - r.axis_start;
                    r.axis_len = norm(axis_vec);
                    r.axis_dir = axis_vec / r.axis_len;

                case 'BLoftSurf'
                    r.knots_v = mwecmass.geometry.MS2Parser.make_clamped_knots( ...
                                    numel(r.section_ids), p.degree);
            end
        end


        function id = id_of(obj, name)
        % ID_OF: Number of the entity called name.

            id = find(strcmp(obj.names, name), 1);
            if isempty(id)
                error('mwecmass:geometry:EntityNotFound', ...
                       'Entity not found: %s', name);
            end
        end


        function id = id_unchecked(obj, name)
        % ID_UNCHECKED: Number of the entity called name; an unknown name raises the error of
        % an unchecked entities(name) access.

            id = find(strcmp(obj.names, name), 1);
            if isempty(id)
                obj.entities(name);
            end
        end


        function r = ent_unchecked(obj, id)
        % ENT_UNCHECKED: Entity record number id; a placeholder raises the error of an unchecked
        % entities(name) access.

            r = obj.ent{id};
            if strcmp(r.type, 'Missing')
                r = obj.entities(r.name);
            end
        end


        function pt = point_at(obj, id)
        % POINT_AT: Coordinates [1x3] of a FramePoint or MirrPoint.

            r = obj.ent{id};
            switch r.type
                case {'FramePoint', 'MirrPoint'}
                    if ~r.k_ok, r = obj.compute_constants(r); end
                    pt = r.pt;
                case 'Missing'
                    error('mwecmass:geometry:EntityNotFound', ...
                           'Entity not found: %s', r.name);
                otherwise
                    error('mwecmass:geometry:WrongType', ...
                           '%s is %s, not FramePoint/MirrPoint', r.name, r.type);
            end
        end


        function pt = any_point_at(obj, id)
        % ANY_POINT_AT: Coordinates [1x3] of a point, bead or ring.

            r = obj.ent{id};
            switch r.type
                case {'FramePoint', 'MirrPoint'}
                    pt = obj.point_at(id);

                case {'AbsBead', 'AbsRing'}
                    if ~r.k_ok, r = obj.compute_constants(r); end
                    pt = r.pt;

                case 'Missing'
                    error('mwecmass:geometry:EntityNotFound', ...
                           'Entity not found: %s', r.name);

                otherwise
                    error('mwecmass:geometry:CannotResolvePoint', ...
                           'Cannot resolve %s (type: %s) to a point', ...
                           r.name, r.type);
            end
        end


        function pts = curve_or_snake_at(obj, id, t)
        % CURVE_OR_SNAKE_AT: Evaluate curve or snake number id at parameters t; points are [Nx3].

            r = obj.ent{id};
            if strcmp(r.type, 'Missing')
                error('mwecmass:geometry:EntityNotFound', ...
                       'Entity not found: %s', r.name);
            end

            if any(strcmp(r.type, {'EdgeSnake', 'BSubSnake', 'AbsBead'}))
                pts = obj.snake_at(id, t);
            else
                pts = obj.curve_at(id, t);
            end
        end


        function pts = curve_at(obj, id, t)
        % CURVE_AT: Evaluate curve number id at parameters t; points are [Nx3].

            t = t(:);  % force column

            r = obj.ent{id};
            switch r.type

                case 'BCurve'
                    pts = obj.eval_bcurve(r, t);

                case 'Conic'
                    pts = obj.eval_conic(r, t);

                case 'CopyCurve'
                    pts = obj.eval_copy_curve(r, t);

                case 'Line'
                    pts = obj.eval_line(r, t);

                case 'BSubCurve'
                    pts = obj.eval_bsub_curve(r, t);

                case 'Arc'
                    pts = obj.eval_arc(r, t);

                case 'PolyCurve2'
                    pts = obj.eval_polycurve2(r, t);

                case 'ProjCurve'
                    pts = obj.eval_proj_curve(r, t);

                case 'Missing'
                    error('mwecmass:geometry:EntityNotFound', ...
                           'Entity not found: %s', r.name);

                otherwise
                    if any(strcmp(r.type, {'EdgeSnake', 'BSubSnake', 'AbsBead'}))
                        pts = obj.snake_at(id, t);
                    else
                        error('mwecmass:geometry:UnsupportedCurveType', ...
                               'Cannot evaluate %s as curve or snake (type: %s)', ...
                               r.name, r.type);
                    end
            end
        end


        function pts = snake_at(obj, id, t)
        % SNAKE_AT: Evaluate snake number id at parameters t; points are [Nx3].

            t = t(:);

            r = obj.ent{id};
            switch r.type

                case 'EdgeSnake'
                    pts = obj.eval_edge_snake(r, t);

                case 'BSubSnake'
                    pts = obj.eval_bsub_snake(r, t);

                case 'AbsBead'
                    pt = obj.snake_at(r.parent_id, r.params.parameter);
                    pts = repmat(pt, length(t), 1);

                case 'Missing'
                    error('mwecmass:geometry:EntityNotFound', ...
                           'Entity not found: %s', r.name);

                otherwise
                    pts = obj.curve_at(id, t);
            end
        end


        function pt = surface_at(obj, id, u, v)
        % SURFACE_AT: Evaluate surface number id at normalized parameters u and v.

            r = obj.ent{id};
            switch r.type

                case 'RuledSurf'
                    pt = obj.eval_ruled_surf(r, u, v);

                case 'RevSurf'
                    pt = obj.eval_rev_surf(r, u, v);

                case 'BLoftSurf'
                    pt = obj.eval_bloft_surf(r, u, v);

                case 'DevSurf'
                    pt = obj.eval_dev_surf(r, u, v);

                case 'MirrSurf'
                    pt = obj.eval_mirr_surf(r, u, v);

                case 'Missing'
                    error('mwecmass:geometry:EntityNotFound', ...
                           'Entity not found: %s', r.name);

                otherwise
                    error('mwecmass:geometry:UnsupportedSurface', ...
                           'Cannot evaluate %s as surface (type: %s)', ...
                           r.name, r.type);
            end
        end


        function [pts, dpts] = curve_deriv_at(obj, id, t)
        % CURVE_DERIV_AT: Evaluate curve number id and its parameter derivative at t.

            t = t(:);

            r = obj.ent{id};
            switch r.type
                case 'BCurve'
                    [pts, dpts] = obj.eval_bcurve_deriv(r, t);

                case 'Conic'
                    [pts, dpts] = obj.eval_conic_deriv(r, t);

                case 'CopyCurve'
                    [pts, dpts] = obj.eval_copy_curve_deriv(r, t);

                case 'Line'
                    [pts, dpts] = obj.eval_line_deriv(r, t);

                case 'BSubCurve'
                    [pts, dpts] = obj.eval_bsub_curve_deriv(r, t);

                case 'Arc'
                    [pts, dpts] = obj.eval_arc_deriv(r, t);

                case 'PolyCurve2'
                    [pts, dpts] = obj.eval_polycurve2_deriv(r, t);

                case 'ProjCurve'
                    [pts, dpts] = obj.eval_proj_curve_deriv(r, t);

                case 'Missing'
                    error('mwecmass:geometry:EntityNotFound', ...
                           'Entity not found: %s', r.name);

                otherwise
                    if any(strcmp(r.type, {'EdgeSnake', 'BSubSnake', 'AbsBead'}))
                        [pts, dpts] = obj.snake_deriv_at(id, t);
                    else
                        h = 1e-7;
                        pts = obj.curve_at(id, t);
                        pts_h = obj.curve_at(id, min(t+h, 1));
                        pts_l = obj.curve_at(id, max(t-h, 0));
                        dpts = (pts_h - pts_l) ./ (min(t+h,1) - max(t-h,0));
                    end
            end
        end


        function [pts, dpts] = snake_deriv_at(obj, id, t)
        % SNAKE_DERIV_AT: Evaluate snake number id and its parameter derivative at t.

            t = t(:);
            r = obj.ent_unchecked(id);

            switch r.type
                case 'EdgeSnake'
                    [pts, dpts] = obj.eval_edge_snake_deriv(r, t);

                case 'BSubSnake'
                    [pts, dpts] = obj.eval_bsub_snake_deriv(r, t);

                otherwise
                    h = 1e-7;
                    pts = obj.snake_at(id, t);
                    pts_h = obj.snake_at(id, min(t+h, 1));
                    pts_l = obj.snake_at(id, max(t-h, 0));
                    dpts = (pts_h - pts_l) ./ (min(t+h,1) - max(t-h,0));
            end
        end


        function [S, Su, Sv] = surface_deriv_at(obj, id, u, v)
        % SURFACE_DERIV_AT: Evaluate surface number id with its parameter derivatives.

            r = obj.ent{id};
            switch r.type
                case 'RuledSurf'
                    [S, Su, Sv] = obj.eval_ruled_surf_derivs(r, u, v);
                case 'RevSurf'
                    [S, Su, Sv] = obj.eval_rev_surf_derivs(r, u, v);
                case 'BLoftSurf'
                    [S, Su, Sv] = obj.eval_bloft_surf_derivs(r, u, v);
                case 'DevSurf'
                    [S, Su, Sv] = obj.eval_dev_surf_derivs(r, u, v);
                case 'MirrSurf'
                    [S, Su, Sv] = obj.eval_mirr_surf_derivs(r, u, v);
                case 'Missing'
                    error('mwecmass:geometry:EntityNotFound', ...
                           'Entity not found: %s', r.name);
                otherwise
                    h = 1e-6;
                    S = obj.surface_at(id, u, v);
                    Su = (obj.surface_at(id, min(u+h,1), v) - ...
                          obj.surface_at(id, max(u-h,0), v)) / ...
                         (min(u+h,1) - max(u-h,0));
                    Sv = (obj.surface_at(id, u, min(v+h,1)) - ...
                          obj.surface_at(id, u, max(v-h,0))) / ...
                         (min(v+h,1) - max(v-h,0));
            end
        end

        function v = ref_ids(~, ids, names, id_empty)
        % REF_IDS: Entity numbers for references by name; an empty name gives id_empty.

            if ischar(names)
                names = {names};
            end
            v = zeros(1, numel(names));
            for k = 1:numel(names)
                if isempty(names{k})
                    v(k) = id_empty;
                elseif ids.isKey(names{k})
                    v(k) = ids(names{k});
                else
                    error('mwecmass:geometry:EntityNotFound', ...
                           'Entity not found: %s', names{k});
                end
            end
        end


    end % methods (Access = private)



    methods
        % The eval_* methods that take a name are the entry points. Those named after an entity type
        % (eval_bcurve, eval_ruled_surf, ...) take the entity struct of model.entities; the
        % internal calls pass the resolved record, which they recognise by its field k_ok.

        function pt = eval_point(obj, name)
        % EVAL_POINT: Evaluate the entity at normalized parameter values; points are [Nx3].

            pt = obj.point_at(obj.id_of(name));
        end

        function pt = eval_any_point(obj, name)
        % EVAL_ANY_POINT: Evaluate the entity at normalized parameter values; points are [Nx3].

            pt = obj.any_point_at(obj.id_of(name));
        end


        function pts = eval_curve_or_snake(obj, name, t)
        % EVAL_CURVE_OR_SNAKE: Evaluate the entity at normalized parameter values; points are [Nx3].

            pts = obj.curve_or_snake_at(obj.id_of(name), t);
        end

        function pts = eval_curve(obj, name, t)
        % EVAL_CURVE: Evaluate the entity at normalized parameter values; points are [Nx3].

            pts = obj.curve_at(obj.id_of(name), t);
        end

        function pts = eval_bcurve(obj, e, t)
        % EVAL_BCURVE: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            if ~e.k_ok, e = obj.compute_constants(e); end
            pts = mwecmass.geometry.MS2Parser.bspline_curve_eval( ...
                      e.knots, e.ctrl_pts, e.params.degree, t);
        end

        function pts = eval_conic(obj, e, t)
        % EVAL_CONIC: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;
            if ~e.k_ok, e = obj.compute_constants(e); end
            center = e.center;
            rad_pt = e.rad_pt;
            apex   = e.apex;

            vec_a = rad_pt - center;
            vec_b = apex - center;
            a = norm(vec_a);
            b = norm(vec_b);
            e_a = vec_a / a;
            e_b = vec_b / b;

            theta = deg2rad(p.angle_start + t * (p.angle_end - p.angle_start));

            pts = center + a * cos(theta) .* e_a + b * sin(theta) .* e_b;
        end

        function pts = eval_copy_curve(obj, e, t)
        % EVAL_COPY_CURVE: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;
            if ~e.k_ok, e = obj.compute_constants(e); end
            src_pos = e.src_pos;
            dst_pos = e.dst_pos;
            offset  = dst_pos - src_pos;

            src_pts = obj.curve_at(e.source_id, t);

            pts = zeros(size(src_pts));
            for i = 1:size(src_pts, 1)
                rel = src_pts(i, :) - src_pos;
                pts(i, :) = src_pos + rel .* p.scale + offset;
            end
        end

        function pts = eval_line(obj, e, t)
        % EVAL_LINE: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            if ~e.k_ok, e = obj.compute_constants(e); end
            p_start = e.p_start;
            p_end   = e.p_end;

            pts = (1 - t) .* p_start + t .* p_end;
        end

        function pts = eval_bsub_curve(obj, e, t)
        % EVAL_BSUB_CURVE: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            if ~e.k_ok, e = obj.compute_constants(e); end
            t_start = e.t_start;
            t_end   = e.t_end;

            t_parent = t_start + t * (t_end - t_start);

            pts = obj.curve_at(e.parent_id, t_parent);
        end


        function pts = eval_arc(obj, e, t)
        % EVAL_ARC: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            if ~e.k_ok, e = obj.compute_constants(e); end

            [pts, ~] = mwecmass.geometry.MS2Parser.arc_evaluate(e.p_start, e.p_centre, e.p_end, t);
        end


        function pts = eval_polycurve2(obj, e, t)
        % EVAL_POLYCURVE2: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            cids = e.curve_ids;
            nc = length(cids);
            pts = zeros(length(t), 3);

            for k = 1:length(t)
                if t(k) >= 1 - 1e-10
                    p = obj.curve_at(cids(nc), 1);
                    pts(k, :) = p(1, :);
                    continue;
                end
                tk = max(t(k), 0);
                seg = min(floor(tk * nc) + 1, nc);     % segment index 1..nc
                t_local = tk * nc - (seg - 1);          % local parameter [0, 1)
                t_local = min(max(t_local, 0), 1);
                p = obj.curve_at(cids(seg), t_local);
                pts(k, :) = p(1, :);
            end
        end


        function pts = eval_proj_curve(obj, e, t)
        % EVAL_PROJ_CURVE: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;

            se = obj.ent{e.source_id};
            if any(strcmp(se.type, {'EdgeSnake', 'BSubSnake', 'AbsBead'}))
                pts = obj.snake_at(e.source_id, t);
            else
                pts = obj.curve_at(e.source_id, t);
            end

            switch p.proj_plane
                case 'X', pts(:,1) = 0;
                case 'Y', pts(:,2) = 0;
                case 'Z', pts(:,3) = 0;
            end
        end


        function [pts, dpts] = eval_arc_deriv(obj, e, t)
        % EVAL_ARC_DERIV: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            if ~e.k_ok, e = obj.compute_constants(e); end

            [pts, dpts] = mwecmass.geometry.MS2Parser.arc_evaluate(e.p_start, e.p_centre, e.p_end, t);
        end


        function [pts, dpts] = eval_polycurve2_deriv(obj, e, t)
        % EVAL_POLYCURVE2_DERIV: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            cids = e.curve_ids;
            nc = length(cids);
            pts  = zeros(length(t), 3);
            dpts = zeros(length(t), 3);

            for k = 1:length(t)
                if t(k) >= 1 - 1e-10
                    [p, dp] = obj.curve_deriv_at(cids(nc), 1);
                    pts(k, :)  = p(1, :);
                    dpts(k, :) = dp(1, :) * nc;
                    continue;
                end
                tk = max(t(k), 0);
                seg = min(floor(tk * nc) + 1, nc);
                t_local = tk * nc - (seg - 1);
                t_local = min(max(t_local, 0), 1);

                [p, dp] = obj.curve_deriv_at(cids(seg), t_local);
                pts(k, :)  = p(1, :);
                dpts(k, :) = dp(1, :) * nc;   % chain rule: dt_local/dt = nc
            end
        end


        function [pts, dpts] = eval_proj_curve_deriv(obj, e, t)
        % EVAL_PROJ_CURVE_DERIV: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;

            se = obj.ent{e.source_id};
            if any(strcmp(se.type, {'EdgeSnake', 'BSubSnake', 'AbsBead'}))
                [pts, dpts] = obj.snake_deriv_at(e.source_id, t);
            else
                [pts, dpts] = obj.curve_deriv_at(e.source_id, t);
            end

            switch p.proj_plane
                case 'X', pts(:,1) = 0; dpts(:,1) = 0;
                case 'Y', pts(:,2) = 0; dpts(:,2) = 0;
                case 'Z', pts(:,3) = 0; dpts(:,3) = 0;
            end
        end



        function pts = eval_snake(obj, name, t)
        % EVAL_SNAKE: Evaluate the entity at normalized parameter values; points are [Nx3].

            pts = obj.snake_at(obj.id_of(name), t);
        end

        function pts = eval_edge_snake(obj, e, t)
        % EVAL_EDGE_SNAKE: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;
            edge_idx = p.edge_index;
            sid = e.surface_id;

            N = length(t);
            pts = zeros(N, 3);

            switch edge_idx
                case 1  % v = 0, u = t
                    for i = 1:N
                        pts(i,:) = obj.surface_at(sid, t(i), 0);
                    end
                case 2  % u = 1, v = t
                    for i = 1:N
                        pts(i,:) = obj.surface_at(sid, 1, t(i));
                    end
                case 3  % v = 1, u = t (natural direction)
                    for i = 1:N
                        pts(i,:) = obj.surface_at(sid, t(i), 1);
                    end
                case 4  % u = 0, v = t (natural direction)
                    for i = 1:N
                        pts(i,:) = obj.surface_at(sid, 0, t(i));
                    end
                otherwise
                    error('mwecmass:geometry:BadEdge', ...
                           'Invalid edge index %d for %s', edge_idx, p.surface_name);
            end
        end

        function pts = eval_bsub_snake(obj, e, t)
        % EVAL_BSUB_SNAKE: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            if ~e.k_ok, e = obj.compute_constants(e); end
            t_start = e.t_start;
            t_end   = e.t_end;

            t_parent = t_start + t * (t_end - t_start);

            pts = obj.snake_at(e.parent_id, t_parent);
        end



        function pt = eval_surface(obj, name, u, v)
        % EVAL_SURFACE: Evaluate the entity at normalized parameter values; points are [Nx3].

            pt = obj.surface_at(obj.id_of(name), u, v);
        end

        function S = eval_surface_grid(obj, name, u_grid, v_grid)
        % EVAL_SURFACE_GRID: Evaluate the entity at normalized parameter values; points are [Nx3].

            Nu = length(u_grid);
            Nv = length(v_grid);
            S = zeros(Nu, Nv, 3);
            if Nu * Nv == 0, return; end
            id = obj.id_of(name);
            for i = 1:Nu
                for j = 1:Nv
                    S(i, j, :) = obj.surface_at(id, u_grid(i), v_grid(j));
                end
            end
        end



        function [pts, dpts] = eval_curve_with_deriv(obj, name, t)
        % EVAL_CURVE_WITH_DERIV: Evaluate the entity at normalized parameter values; points are [Nx3].

            [pts, dpts] = obj.curve_deriv_at(obj.id_of(name), t);
        end


        function [pts, dpts] = eval_bcurve_deriv(obj, e, t)
        % EVAL_BCURVE_DERIV: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            if ~e.k_ok, e = obj.compute_constants(e); end
            [pts, dpts] = mwecmass.geometry.MS2Parser.bspline_curve_eval_with_deriv( ...
                e.knots, e.ctrl_pts, e.params.degree, t);
        end


        function [pts, dpts] = eval_conic_deriv(obj, e, t)
        % EVAL_CONIC_DERIV: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;
            if ~e.k_ok, e = obj.compute_constants(e); end
            center = e.center;
            rad_pt = e.rad_pt;
            apex   = e.apex;

            vec_a = rad_pt - center;
            vec_b = apex - center;
            a = norm(vec_a);  b = norm(vec_b);
            e_a = vec_a / a;  e_b = vec_b / b;

            theta_start = deg2rad(p.angle_start);
            theta_end   = deg2rad(p.angle_end);
            dtheta_dt = theta_end - theta_start;

            theta = theta_start + t * dtheta_dt;

            pts  = center + a * cos(theta) .* e_a + b * sin(theta) .* e_b;
            dpts = (-a * sin(theta) .* e_a + b * cos(theta) .* e_b) * dtheta_dt;
        end


        function [pts, dpts] = eval_copy_curve_deriv(obj, e, t)
        % EVAL_COPY_CURVE_DERIV: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;
            if ~e.k_ok, e = obj.compute_constants(e); end
            src_pos = e.src_pos;
            dst_pos = e.dst_pos;
            offset  = dst_pos - src_pos;

            [src_pts, src_dpts] = obj.curve_deriv_at(e.source_id, t);

            pts = zeros(size(src_pts));
            dpts = zeros(size(src_dpts));
            for i = 1:size(src_pts, 1)
                rel = src_pts(i, :) - src_pos;
                pts(i, :)  = src_pos + rel .* p.scale + offset;
                dpts(i, :) = src_dpts(i, :) .* p.scale;
            end
        end


        function [pts, dpts] = eval_line_deriv(obj, e, t)
        % EVAL_LINE_DERIV: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            if ~e.k_ok, e = obj.compute_constants(e); end
            p_start = e.p_start;
            p_end   = e.p_end;

            pts  = (1 - t) .* p_start + t .* p_end;
            dpts = repmat(p_end - p_start, length(t), 1);
        end


        function [pts, dpts] = eval_bsub_curve_deriv(obj, e, t)
        % EVAL_BSUB_CURVE_DERIV: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            if ~e.k_ok, e = obj.compute_constants(e); end
            t_start = e.t_start;
            t_end   = e.t_end;

            t_parent = t_start + t * (t_end - t_start);

            [pts, dpts_parent] = obj.curve_deriv_at(e.parent_id, t_parent);
            dpts = dpts_parent * (t_end - t_start);
        end


        function [pts, dpts] = eval_snake_with_deriv(obj, name, t)
        % EVAL_SNAKE_WITH_DERIV: Evaluate the entity at normalized parameter values; points are [Nx3].

            [pts, dpts] = obj.snake_deriv_at(obj.id_unchecked(name), t);
        end


        function [pts, dpts] = eval_edge_snake_deriv(obj, e, t)
        % EVAL_EDGE_SNAKE_DERIV: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;
            sid = e.surface_id;
            N = length(t);
            pts  = zeros(N, 3);
            dpts = zeros(N, 3);

            for i = 1:N
                switch p.edge_index
                    case 1
                        [S, Su, ~] = obj.surface_deriv_at(sid, t(i), 0);
                        pts(i,:) = S; dpts(i,:) = Su;
                    case 2
                        [S, ~, Sv] = obj.surface_deriv_at(sid, 1, t(i));
                        pts(i,:) = S; dpts(i,:) = Sv;
                    case 3
                        [S, Su, ~] = obj.surface_deriv_at(sid, t(i), 1);
                        pts(i,:) = S; dpts(i,:) = Su;
                    case 4
                        [S, ~, Sv] = obj.surface_deriv_at(sid, 0, t(i));
                        pts(i,:) = S; dpts(i,:) = Sv;
                end
            end
        end


        function [pts, dpts] = eval_bsub_snake_deriv(obj, e, t)
        % EVAL_BSUB_SNAKE_DERIV: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            if ~e.k_ok, e = obj.compute_constants(e); end
            t_start = e.t_start;
            t_end   = e.t_end;

            t_parent = t_start + t * (t_end - t_start);

            [pts, dpts_parent] = obj.snake_deriv_at(e.parent_id, t_parent);
            dpts = dpts_parent * (t_end - t_start);
        end


        function [S, Su, Sv] = eval_surface_with_derivs(obj, name, u, v)
        % EVAL_SURFACE_WITH_DERIVS: Evaluate the entity at normalized parameter values; points are [Nx3].

            [S, Su, Sv] = obj.surface_deriv_at(obj.id_of(name), u, v);
        end


        function [S, Su, Sv] = eval_ruled_surf_derivs(obj, e, u, v)
        % EVAL_RULED_SURF_DERIVS: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            [c1, dc1] = obj.curve_deriv_at(e.curve1_id, u);
            [c2, dc2] = obj.curve_deriv_at(e.curve2_id, u);

            S  = (1 - v) * c1 + v * c2;
            Su = (1 - v) * dc1 + v * dc2;
            Sv = c2 - c1;
        end


        function [S, Su, Sv] = eval_rev_surf_derivs(obj, e, u, v)
        % EVAL_REV_SURF_DERIVS: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;

            [profile_pt, profile_dpdt] = obj.curve_deriv_at(e.profile_id, u);
            profile_pt = profile_pt(1,:);
            profile_dpdt = profile_dpdt(1,:);

            if ~e.k_ok, e = obj.compute_constants(e); end
            axis_start = e.axis_start;
            axis_dir   = e.axis_dir;

            v_rel = profile_pt - axis_start;
            z_along = dot(v_rel, axis_dir);
            proj = axis_start + z_along * axis_dir;
            radial = profile_pt - proj;
            r = norm(radial);

            phi_start = deg2rad(p.angle_start);
            phi_end   = deg2rad(p.angle_end);
            dphi_dv = phi_end - phi_start;
            phi = phi_start + v * dphi_dv;

            if r < 1e-12
                S = profile_pt;
                Su = profile_dpdt;
                Sv = [0, 0, 0];
                return;
            end

            e_r = radial / r;
            e_t = cross(axis_dir, e_r);
            e_t = e_t / norm(e_t);

            S = proj + r * cos(phi) * e_r + r * sin(phi) * e_t;

            Sv = r * (-sin(phi) * e_r + cos(phi) * e_t) * dphi_dv;

            dz_du = dot(profile_dpdt, axis_dir);
            dproj_du = dz_du * axis_dir;
            dradial_du = profile_dpdt - dproj_du;
            dr_du = dot(dradial_du, e_r);

            if r > 1e-10
                de_r_du = (dradial_du - dr_du * e_r) / r;
                de_t_du = cross(axis_dir, de_r_du);
            else
                de_r_du = [0, 0, 0];
                de_t_du = [0, 0, 0];
            end

            Su = dproj_du + dr_du * cos(phi) * e_r + r * cos(phi) * de_r_du ...
                          + dr_du * sin(phi) * e_t + r * sin(phi) * de_t_du;
        end


        function [S, Su, Sv] = eval_bloft_surf_derivs(obj, e, u, v)
        % EVAL_BLOFT_SURF_DERIVS: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;
            n_sec = numel(e.section_ids);
            degree = p.degree;

            persistent bloft_cache_key bloft_cache_u bloft_cache_pts ...
                       bloft_cache_dpts bloft_cache_knots;

            surf_key = e.surf_key;

            if ~isempty(bloft_cache_key) && ...
                    strcmp(bloft_cache_key, surf_key) && ...
                    abs(bloft_cache_u - u) < 1e-14
                sec_pts  = bloft_cache_pts;
                sec_dpts = bloft_cache_dpts;
                knots_v  = bloft_cache_knots;
            else
                sec_pts  = zeros(n_sec, 3);
                sec_dpts = zeros(n_sec, 3);
                for k = 1:n_sec
                    sid = e.section_ids(k);
                    se = obj.ent{sid};
                    if strcmp(se.type, 'Missing')
                        sec_pts(k,:)  = [0 0 0];
                        sec_dpts(k,:) = [0 0 0];
                        continue;
                    end
                    if any(strcmp(se.type, {'EdgeSnake', 'BSubSnake', 'AbsBead'}))
                        [pp, dd] = obj.snake_deriv_at(sid, u);
                    else
                        [pp, dd] = obj.curve_deriv_at(sid, u);
                    end
                    sec_pts(k,:)  = pp(1,:);
                    sec_dpts(k,:) = dd(1,:);
                end
                if ~e.k_ok, e = obj.compute_constants(e); end
                knots_v = e.knots_v;

                bloft_cache_key   = surf_key;
                bloft_cache_u     = u;
                bloft_cache_pts   = sec_pts;
                bloft_cache_dpts  = sec_dpts;
                bloft_cache_knots = knots_v;
            end

            [S_arr, Sv_arr] = mwecmass.geometry.MS2Parser.bspline_curve_eval_with_deriv( ...
                knots_v, sec_pts, degree, v);
            S  = S_arr(1,:);
            Sv = Sv_arr(1,:);

            Su_arr = mwecmass.geometry.MS2Parser.bspline_curve_eval(knots_v, sec_dpts, degree, v);
            Su = Su_arr(1,:);
        end


        function [S, Su, Sv] = eval_dev_surf_derivs(obj, e, u, v)
        % EVAL_DEV_SURF_DERIVS: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            [s_pt, s_dpt] = obj.snake_deriv_at(e.snake_id, u);
            [c_pt, c_dpt] = obj.curve_deriv_at(e.curve_id, u);

            s_pt = s_pt(1,:); s_dpt = s_dpt(1,:);
            c_pt = c_pt(1,:); c_dpt = c_dpt(1,:);

            S  = (1 - v) * s_pt + v * c_pt;
            Su = (1 - v) * s_dpt + v * c_dpt;
            Sv = c_pt - s_pt;
        end


        function [S, Su, Sv] = eval_mirr_surf_derivs(obj, e, u, v)
        % EVAL_MIRR_SURF_DERIVS: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;
            [S, Su, Sv] = obj.surface_deriv_at(e.source_id, u, v);

            switch p.mirror_plane
                case 'Y'
                    S(2) = -S(2);
                    Su(2) = -Su(2);
                    Sv(2) = -Sv(2);
                case 'X'
                    S(1) = -S(1);
                    Su(1) = -Su(1);
                    Sv(1) = -Sv(1);
            end
        end


        function pt = eval_ruled_surf(obj, e, u, v)
        % EVAL_RULED_SURF: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            c1 = obj.curve_at(e.curve1_id, u);
            c2 = obj.curve_at(e.curve2_id, u);
            pt = (1 - v) * c1 + v * c2;
        end


        function pt = eval_rev_surf(obj, e, u, v)
        % EVAL_REV_SURF: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;

            profile_pt = obj.curve_at(e.profile_id, u);

            if ~e.k_ok, e = obj.compute_constants(e); end
            axis_start = e.axis_start;
            axis_dir   = e.axis_dir;
            if e.axis_len < 1e-12
                pt = profile_pt;
                return;
            end

            v_rel = profile_pt - axis_start;
            z_along = dot(v_rel, axis_dir);
            proj = axis_start + z_along * axis_dir;

            radial = profile_pt - proj;
            r = norm(radial);

            if r < 1e-12
                pt = profile_pt;
                return;
            end

            e_r = radial / r;  % initial radial direction
            e_t = cross(axis_dir, e_r);
            e_t = e_t / norm(e_t);  % tangential direction

            phi = deg2rad(p.angle_start + v * (p.angle_end - p.angle_start));

            pt = proj + r * cos(phi) * e_r + r * sin(phi) * e_t;
        end


        function pt = eval_bloft_surf(obj, e, u, v)
        % EVAL_BLOFT_SURF: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;
            n_sections = numel(e.section_ids);
            degree = p.degree;

            section_pts = zeros(n_sections, 3);
            for k = 1:n_sections
                pts_k = obj.curve_or_snake_at(e.section_ids(k), u);
                section_pts(k, :) = pts_k(1, :);  % ensure [1×3]
            end

            if ~e.k_ok, e = obj.compute_constants(e); end

            pt = mwecmass.geometry.MS2Parser.bspline_curve_eval(e.knots_v, section_pts, degree, v);
        end

        function pts = eval_bloft_surf_at_u(obj, e, u, v_array)
        % EVAL_BLOFT_SURF_AT_U: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;
            n_sections = numel(e.section_ids);
            degree = p.degree;

            section_pts = zeros(n_sections, 3);
            for k = 1:n_sections
                pts_k = obj.curve_or_snake_at(e.section_ids(k), u);
                section_pts(k, :) = pts_k(1, :);
            end

            if ~e.k_ok, e = obj.compute_constants(e); end

            v_array = v_array(:);
            pts = mwecmass.geometry.MS2Parser.bspline_curve_eval(e.knots_v, section_pts, degree, v_array);
        end


        function pt = eval_dev_surf(obj, e, u, v)
        % EVAL_DEV_SURF: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            pt_snake = obj.snake_at(e.snake_id, u);
            pt_curve = obj.curve_at(e.curve_id, u);

            pt = (1 - v) * pt_snake + v * pt_curve;
        end


        function pt = eval_mirr_surf(obj, e, u, v)
        % EVAL_MIRR_SURF: Evaluate the entity at normalized parameter values; points are [Nx3].

            if ~isfield(e, 'k_ok'), e = obj.record_of(e); end

            p = e.params;
            pt = obj.surface_at(e.source_id, u, v);

            switch p.mirror_plane
                case 'Y'
                    pt(2) = -pt(2);
                case 'X'
                    pt(1) = -pt(1);
            end
        end

    end % methods



    methods (Static)

        function knots = make_clamped_knots(n_ctrl, degree)
        % MAKE_CLAMPED_KNOTS: Build an open knot vector on [0,1].

            n_knots   = n_ctrl + degree + 1;
            n_interior = n_knots - 2 * (degree + 1);

            if n_interior <= 0
                knots = [zeros(1, degree+1), ones(1, degree+1)];
            else
                interior = linspace(0, 1, n_interior + 2);
                interior = interior(2:end-1);
                knots = [zeros(1, degree+1), interior, ones(1, degree+1)];
            end
        end

        function pts = bspline_curve_eval(knots, ctrl_pts, degree, t)
        % BSPLINE_CURVE_EVAL: Evaluate a B-spline using the supplied knots and control points.

            t = t(:);
            N = length(t);
            n_ctrl = size(ctrl_pts, 1);
            pts = zeros(N, 3);

            for k = 1:N
                tk = min(max(t(k), 0), 1 - 1e-12);

                basis = mwecmass.geometry.MS2Parser.bspline_basis_all(knots, degree, tk, n_ctrl);

                pts(k, :) = basis * ctrl_pts;
            end

            endpoint_mask = (t >= 1 - 1e-10);
            if any(endpoint_mask)
                pts(endpoint_mask, :) = repmat(ctrl_pts(end, :), sum(endpoint_mask), 1);
            end
        end

        function N = bspline_basis_all(knots, degree, t, n_ctrl)
        % BSPLINE_BASIS_ALL: Evaluate a B-spline using the supplied knots and control points.

            n = length(knots) - 1;  % number of basis functions at degree 0

            N0 = zeros(1, n);
            for i = 1:n
                if knots(i) <= t && t < knots(i+1)
                    N0(i) = 1;
                end
            end

            N_prev = N0;
            for p = 1:degree
                N_cur = zeros(1, n - p);
                for i = 1:(n - p)
                    left_denom  = knots(i + p) - knots(i);
                    right_denom = knots(i + p + 1) - knots(i + 1);

                    left_term  = 0;
                    right_term = 0;

                    if abs(left_denom) > 1e-14
                        left_term = (t - knots(i)) / left_denom * N_prev(i);
                    end
                    if abs(right_denom) > 1e-14
                        right_term = (knots(i + p + 1) - t) / right_denom * N_prev(i + 1);
                    end

                    N_cur(i) = left_term + right_term;
                end
                N_prev = N_cur;
            end

            N = N_prev(1:n_ctrl);
        end


        function [pts, dpts] = bspline_curve_eval_with_deriv(knots, ctrl_pts, degree, t)
        % BSPLINE_CURVE_EVAL_WITH_DERIV: Evaluate a B-spline using the supplied knots and control points.

            pts = mwecmass.geometry.MS2Parser.bspline_curve_eval(knots, ctrl_pts, degree, t);

            if degree < 1
                dpts = zeros(size(pts));
                return;
            end

            n = size(ctrl_pts, 1);
            Q = zeros(n - 1, 3);
            for i = 1:n-1
                denom = knots(i + degree + 1) - knots(i + 1);
                if abs(denom) > 1e-14
                    Q(i, :) = degree * (ctrl_pts(i+1, :) - ctrl_pts(i, :)) / denom;
                end
            end

            d_knots = knots(2:end-1);
            d_degree = degree - 1;

            dpts = mwecmass.geometry.MS2Parser.bspline_curve_eval(d_knots, Q, d_degree, t);
        end


        function [pts, dpts] = arc_evaluate(p_start, p_centre, p_end, t)
        % ARC_EVALUATE: Evaluate a circular arc and its parameter derivative.

            t = t(:);
            N = length(t);

            v_s = p_start  - p_centre;
            v_e = p_end    - p_centre;
            r   = norm(v_s);

            if r < 1e-14
                pts  = repmat(p_centre, N, 1);
                dpts = zeros(N, 3);
                return;
            end

            e_r = v_s / r;
            cp  = cross(v_s, v_e);
            cp_len = norm(cp);

            if cp_len < 1e-14
                pts  = (1 - t) .* p_start + t .* p_end;
                dpts = repmat(p_end - p_start, N, 1);
                return;
            end

            n_hat = cp / cp_len;
            e_t   = cross(n_hat, e_r);

            theta_total = atan2(dot(v_e, e_t), dot(v_e, e_r));


            theta = t * theta_total;

            pts  = p_centre + r * (cos(theta) .* e_r + sin(theta) .* e_t);
            dpts = r * theta_total * (-sin(theta) .* e_r + cos(theta) .* e_t);
        end

    end % methods (Static)



    methods

        function info = classify_visible_surfaces(obj)
        % CLASSIFY_VISIBLE_SURFACES: Classify visible sources, mirrors, and shared edges.

            sources = {};
            mirrors = struct('name', {}, 'source', {}, 'plane', {});

            for i = 1:length(obj.visible_surfs)
                sname = obj.visible_surfs{i};
                e = obj.entities(sname);

                if strcmp(e.type, 'MirrSurf')
                    m = struct();
                    m.name   = sname;
                    m.source = e.params.source;
                    m.plane  = e.params.mirror_plane;
                    mirrors(end+1) = m; %#ok<AGROW>
                else
                    sources{end+1} = sname; %#ok<AGROW>
                end
            end

            for i = 1:length(mirrors)
                [ult_source, ult_flips] = obj.resolve_mirror_chain(mirrors(i).name);
                mirrors(i).ultimate_source = ult_source;
                mirrors(i).effective_flips = ult_flips;
            end

            adjacency = struct('surf_a', {}, 'edge_a', {}, ...
                               'surf_b', {}, 'edge_b', {}, ...
                               'via_entity', {});

            for i = 1:length(sources)
                sname = sources{i};
                edges_i = obj.get_boundary_entities(sname);

                edge_fields = {'edge1', 'edge3'};  % v=0 and v=1 (named boundaries)
                edge_indices = [1, 3];

                for ei = 1:length(edge_fields)
                    ent_name = edges_i.(edge_fields{ei});
                    if isempty(ent_name), continue; end

                    [found, parent_surf, parent_edge] = ...
                        obj.trace_to_edge_snake_of(ent_name, sources);

                    if found
                        adj = struct();
                        adj.surf_a      = sname;
                        adj.edge_a      = edge_indices(ei);
                        adj.surf_b      = parent_surf;
                        adj.edge_b      = parent_edge;
                        adj.via_entity  = ent_name;
                        adjacency(end+1) = adj; %#ok<AGROW>
                    end
                end
            end

            info = struct();
            info.sources   = sources;
            info.mirrors   = mirrors;
            info.adjacency = adjacency;
        end


        function [ult_source, flips] = resolve_mirror_chain(obj, name)
        % RESOLVE_MIRROR_CHAIN: Resolve a mirror to its source and accumulated coordinate flips.

            flips = {};
            current = name;

            for safety = 1:20  % guard against cycles
                if ~obj.entities.isKey(current), break; end
                e = obj.entities(current);
                if ~strcmp(e.type, 'MirrSurf'), break; end
                flips{end+1} = e.params.mirror_plane; %#ok<AGROW>
                current = e.params.source;
            end

            ult_source = current;
        end


        function edges = get_boundary_entities(obj, surf_name)
        % GET_BOUNDARY_ENTITIES: Return named entities defining each surface edge.

            edges = struct('edge1', '', 'edge2', '', 'edge3', '', 'edge4', '');

            if ~obj.entities.isKey(surf_name), return; end
            e = obj.entities(surf_name);

            switch e.type
                case 'RuledSurf'
                    edges.edge1 = e.params.curve1;
                    edges.edge3 = e.params.curve2;

                case 'DevSurf'
                    edges.edge1 = e.params.snake;
                    edges.edge3 = e.params.curve;

                case 'BLoftSurf'
                    sections = e.params.section_names;
                    if ~isempty(sections)
                        edges.edge1 = sections{1};
                        edges.edge3 = sections{end};
                    end

                case 'RevSurf'

                case 'MirrSurf'
                    src_edges = obj.get_boundary_entities(e.params.source);
                    edges = src_edges;  % same entity names — the flip is applied during eval
            end
        end


        function [found, parent_surf, parent_edge] = trace_to_edge_snake_of(obj, entity_name, source_list)
        % TRACE_TO_EDGE_SNAKE_OF: Trace an entity to a source surface edge.

            found = false;
            parent_surf = '';
            parent_edge = 0;

            if ~obj.entities.isKey(entity_name), return; end
            e = obj.entities(entity_name);

            if strcmp(e.type, 'EdgeSnake')
                parent_surf = e.params.surface_name;
                parent_edge = e.params.edge_index;
                found = any(strcmp(parent_surf, source_list));
                return;
            end

            if any(strcmp(e.type, {'BSubSnake', 'BSubCurve'}))
                if isfield(e.params, 'bead_names') && ~isempty(e.params.bead_names)
                    bead_name = e.params.bead_names{1};
                    if obj.entities.isKey(bead_name)
                        bead = obj.entities(bead_name);
                        if isfield(bead.params, 'parent_curve')
                            [found, parent_surf, parent_edge] = ...
                                obj.trace_to_edge_snake_of(bead.params.parent_curve, source_list);
                        end
                    end
                end
                return;
            end

            if strcmp(e.type, 'AbsBead')
                if isfield(e.params, 'parent_curve')
                    [found, parent_surf, parent_edge] = ...
                        obj.trace_to_edge_snake_of(e.params.parent_curve, source_list);
                end
                return;
            end

        end


        function names = get_required_entities(obj, surface_name)
        % GET_REQUIRED_ENTITIES: Return required dependencies in parent-before-child order.

            visited = containers.Map();
            order   = {};

            function dfs(name)
            % DFS: Visit an entity dependency and append it after its parents.
                if visited.isKey(name), return; end
                visited(name) = true;
                if obj.entities.isKey(name)
                    e = obj.entities(name);
                    for p = 1:length(e.parents)
                        dfs(e.parents{p});
                    end
                end
                order{end+1} = name;
            end

            dfs(surface_name);
            names = order;
        end

        function info = geometry_summary(obj)
        % GEOMETRY_SUMMARY: Print and return parsed-model counts and extents.

            info = struct();
            info.filename      = obj.filename;
            info.n_entities    = obj.entities.Count;
            info.n_visible     = length(obj.visible_surfs);
            info.visible_names = obj.visible_surfs;
            info.extents       = obj.extents;

            fprintf('\n  ┌─────────────────────────────────────────┐\n');
            fprintf('  │  MS2 GEOMETRY SUMMARY                    │\n');
            fprintf('  ├─────────────────────────────────────────┤\n');
            fprintf('  │  File:     %s\n', obj.filename);
            fprintf('  │  Entities: %d total\n', info.n_entities);
            fprintf('  │  Visible surfaces: %d\n', info.n_visible);
            for i = 1:info.n_visible
                e = obj.entities(info.visible_names{i});
                deps = obj.get_required_entities(info.visible_names{i});
                fprintf('  │    %s (%s) — %d dependencies\n', ...
                        info.visible_names{i}, e.type, length(deps));
            end
            fprintf('  │  Extents: [%.1f %.1f %.1f] to [%.1f %.1f %.1f] %s\n', ...
                    obj.extents(1:3), obj.extents(4:6), obj.units);
            fprintf('  └─────────────────────────────────────────┘\n');
        end


        function clear_cache(~)
        % CLEAR_CACHE: Nothing to invalidate: points and revolution axes are evaluated once at parse time.
        end

    end % methods

end % classdef
