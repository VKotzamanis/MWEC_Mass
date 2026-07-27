classdef WEC_MS2_Parser
% WEC_MS2_PARSER  Parse MultiSurf .ms2 files into evaluatable geometry.
%
%   model = WEC_MS2_PARSER.parse(filename)
%
%   Reads a MultiSurf .ms2 ASCII model file and builds an entity graph
%   where each entity (point, curve, surface, snake, bead) can be
%   evaluated at arbitrary parameter values.  The parser resolves the
%   full dependency chain so that evaluating a visible surface
%   automatically triggers evaluation of all its parent entities.
%
%   ARCHITECTURE
%   ────────────────────────────────────────────────────────────────────
%   WEC_MS2_Parser (this file)
%     ├─ §1  File reader — tokenises .ms2 lines into entity structs
%     ├─ §2  Dependency resolver — topological sort of entity graph
%     ├─ §3  Point evaluator — FramePoint recursive resolution
%     ├─ §4  Curve evaluators — BCurve, Conic, CopyCurve, Line
%     ├─ §5  Snake evaluators — EdgeSnake, BSubSnake, AbsBead
%     ├─ §6  Surface evaluators — RuledSurf, RevSurf, BLoftSurf,
%     │                           DevSurf, MirrSurf
%     ├─ §7  B-spline engine — basis functions, curve/surface eval
%     └─ §8  Utility methods — knot vector generation, normal computation
%   ────────────────────────────────────────────────────────────────────
%
%   DESIGN PRINCIPLES
%     1. Geometry is stored as EVALUATABLE OBJECTS, not sampled points.
%        Each curve has an eval(t) method, each surface has eval(u,v).
%        Panelisation is a separate downstream step (WEC_Panelizer).
%
%     2. Only entities REQUIRED by visible surfaces are evaluated.
%        Hidden construction surfaces are traced only when a visible
%        surface depends on them (e.g. EdgeSnake on a hidden RevSurf).
%
%     3. B-spline evaluation is implemented from scratch — no Curve
%        Fitting Toolbox dependency.  This ensures portability across
%        MATLAB installations.
%
%     4. Degree-p B-splines APPROXIMATE control points (do NOT
%        interpolate through them).  This matches MultiSurf's semantics.
%
%   ENTITY TYPES SUPPORTED
%     Points:   FramePoint
%     Curves:   BCurve, Conic, CopyCurve, Line, BSubCurve
%     Snakes:   EdgeSnake, AbsBead, BSubSnake
%     Surfaces: RuledSurf, RevSurf, BLoftSurf, DevSurf, MirrSurf
%
%   EDGE NUMBERING CONVENTION (MultiSurf EdgeSnake)
%     For a surface S(u,v) with u,v ∈ [0,1]:
%       Edge 1:  u varies 0→1, v = 0
%       Edge 2:  v varies 0→1, u = 1
%       Edge 3:  u varies 0→1, v = 1
%       Edge 4:  v varies 0→1, u = 0
%     All edges use natural (non-reversed) parameterization.
%     This was verified empirically against E1.ms2 — reversal on
%     edges 3/4 causes BLoftSurf self-intersection.
%
%   USAGE
%     model = WEC_MS2_Parser.parse('C0.ms2');
%     pt    = model.eval_point('pt3');
%     crv   = model.eval_curve('curve1', linspace(0,1,100));
%     S     = model.eval_surface('surface3', 0.5, 0.3);
%
%   See also: WEC_Panelizer, WEC_HAMS_Interface
%
%   Author:  WEC Optimisation Team
%   Version: 1.0 — Phase 1 initial implementation

    properties
        entities        % containers.Map: name → entity struct
        visible_surfs   % cell array of visible surface names
        units           % string: 'm', 'ft'
        extents         % [xmin ymin zmin xmax ymax zmax]
        filename        % source .ms2 file path
        file_symmetry   % cell array: e.g. {'x','y'} from header
    end

    properties (Access = private)
        point_cache     % containers.Map: FramePoint name → [1×3] coordinates
        rev_axis_cache  % containers.Map: axis entity name → struct(start, end, dir, len)
    end

    %% ═══════════════════════════════════════════════════════════════
    %%  §1  FILE READER — PUBLIC STATIC CONSTRUCTOR
    %% ═══════════════════════════════════════════════════════════════

    methods (Static)

        function model = parse(filename)
        % PARSE  Read a MultiSurf .ms2 file and build the entity graph.
        %
        %   model = WEC_MS2_Parser.parse('hull.ms2')
        %
        %   Reads every entity between BeginModel and EndModel,
        %   identifies visible surfaces, and resolves dependencies.
        %
        %   The returned model object can then evaluate any entity
        %   via eval_point, eval_curve, eval_surface, eval_snake.

            if ~exist(filename, 'file')
                error('WEC_MS2_Parser:FileNotFound', ...
                       'MS2 file not found: %s', filename);
            end

            model = WEC_MS2_Parser();
            model.filename = filename;
            model.entities = containers.Map();
            model.point_cache    = containers.Map();
            model.rev_axis_cache = containers.Map();

            % Read entire file
            fid = fopen(filename, 'r');
            raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
            fclose(fid);
            lines = raw{1};

            % ── Parse header ──────────────────────────────────────
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

            % ── Parse entities between BeginModel and EndModel ────
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

                % Skip empty lines and attribute lines
                if isempty(L), continue; end
                if startsWith(L, 'Attribute:'), continue; end

                % Accumulate multi-line definitions (continuation lines
                % start with whitespace or '{')
                if ~isempty(line_buf)
                    % Check if this is a continuation line
                    if startsWith(lines{i}, ' ') || startsWith(lines{i}, char(9)) ...
                            || startsWith(L, '{') || startsWith(L, 'A:')
                        line_buf = [line_buf ' ' L];
                        if endsWith(L, ';')
                            model.parse_entity_line(line_buf);
                            line_buf = '';
                        end
                        continue;
                    else
                        % Previous line was complete, parse it
                        if ~isempty(line_buf)
                            model.parse_entity_line(line_buf);
                        end
                        line_buf = '';
                    end
                end

                % Start new entity
                if endsWith(L, ';')
                    model.parse_entity_line(L);
                else
                    line_buf = L;
                end
            end

            % Parse any remaining buffer
            if ~isempty(line_buf)
                model.parse_entity_line(line_buf);
            end

            % ── Identify visible surfaces ─────────────────────────
            %  VISIBILITY RULE (MultiSurf .ms2 format):
            %    tokens{3} = color code (display only, NOT visibility)
            %    tokens{4} = layer number
            %      Positive layer → VISIBLE (part of the hull)
            %      Negative layer → HIDDEN (construction geometry)
            %
            %  The parser stores the layer number in e.visibility.
            %  Example from E1.ms2:
            %    surface1 layer=-11 → hidden (full-ellipsoid construction)
            %    surface3 layer=3   → visible (lower half-ellipsoid)
            model.visible_surfs = {};
            keys = model.entities.keys();
            for i = 1:length(keys)
                e = model.entities(keys{i});
                if e.is_surface && e.visibility > 0
                    model.visible_surfs{end+1} = keys{i};
                end
            end

            % ── Synthetic mirrors from file-level symmetry ─────────
            %  MultiSurf can declare bilateral symmetry in the file
            %  header ('Symmetry: x y') instead of using explicit
            %  MirrSurf entities.  When this is the case, the .ms2
            %  file contains only source surfaces for one quadrant.
            %
            %  We create synthetic MirrSurf entities so that the rest
            %  of the pipeline (classify_visible_surfaces, HydroProperties,
            %  Panelizer) works identically to files with explicit mirrors.
            %
            %  For 'Symmetry: x y':
            %    source1, source2                    (original quadrant)
            %    source1_mirrX, source2_mirrX        (X-reflected)
            %    source1_mirrY, source2_mirrY        (Y-reflected)
            %    source1_mirrX_mirrY, ...            (both reflected)
            %  Total: 4× the number of sources = full body.

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
                        em.parents    = {orig_sources{s}};
                        model.entities(mname) = em;
                        model.visible_surfs{end+1} = mname;
                        n_synth = n_synth + 1;
                    end
                end

                if any(strcmp(model.file_symmetry, 'y'))
                    % Mirror ALL current visible surfaces (originals + MirrX if created)
                    surfs_to_mirror = model.visible_surfs;
                    for s = 1:length(surfs_to_mirror)
                        mname = [surfs_to_mirror{s} '_mirrY'];
                        em = struct();
                        em.type       = 'MirrSurf';
                        em.visibility = 1;
                        em.is_surface = true;
                        em.params     = struct('source', surfs_to_mirror{s}, ...
                                               'mirror_plane', 'Y');
                        em.parents    = {surfs_to_mirror{s}};
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

            fprintf('  MS2 Parser: %s\n', filename);
            fprintf('    Entities: %d total, %d visible surfaces\n', ...
                    model.entities.Count, length(model.visible_surfs));
            for i = 1:length(model.visible_surfs)
                e = model.entities(model.visible_surfs{i});
                fprintf('      %s (%s)\n', model.visible_surfs{i}, e.type);
            end
        end

    end % methods (Static)


    %% ═══════════════════════════════════════════════════════════════
    %%  §1b  ENTITY LINE PARSER (private)
    %% ═══════════════════════════════════════════════════════════════

    methods (Access = private)

        function parse_entity_line(obj, line)
        % PARSE_ENTITY_LINE  Parse a single .ms2 entity definition.
        %
        %   Dispatches to type-specific parsers based on the first token.
        %   Each parser creates a struct with fields:
        %     .type       — string: 'FramePoint', 'BCurve', etc.
        %     .visibility — integer: positive = visible, negative = hidden
        %     .is_surface — logical: true for surface entities
        %     .params     — struct: type-specific parameters
        %     .parents    — cell array of parent entity names

            % Remove trailing semicolon and clean
            line = strtrim(regexprep(line, ';$', ''));
            if isempty(line), return; end

            tokens = strsplit(line);
            if isempty(tokens), return; end

            etype = tokens{1};

            % Skip non-geometry entities
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
                    % Silently skip unknown types
                end
            catch ME
                warning('WEC_MS2_Parser:ParseError', ...
                        'Failed to parse: %s\n  Error: %s', line, ME.message);
            end
        end

        %% ─── POINT PARSER ────────────────────────────────────────

        function parse_frame_point(obj, tokens, ~)
        % PARSE_FRAME_POINT  FramePoint name vis layer / flags parent1 parent2 dx dy dz
        %
        %   FramePoint definition in .ms2:
        %     FramePoint name 14 -1 / 0 parent1 parent2 dx dy dz
        %
        %   '*' in parent slot means 'no parent' (absolute position or
        %   offset from the other parent).

            name = tokens{2};
            vis  = str2double(tokens{4});

            % Find '/' separator
            slash_idx = find(strcmp(tokens, '/'));
            if isempty(slash_idx), return; end

            % Tokens after '/': flags parent1 parent2 dx dy dz
            after_slash = tokens(slash_idx+1 : end);
            if length(after_slash) < 6, return; end

            flags   = after_slash{1};
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
        % PARSE_MIRR_POINT  MirrPoint name vis layer / source_pt *plane ;
        %
        %   Mirrors a source point about a coordinate plane.
        %   *X=0 → negate x,  *Y=0 → negate y,  *Z=0 → negate z.

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

        %% ─── CURVE PARSERS ──────────────────────────────────────

        function parse_bcurve(obj, tokens, line)
        % PARSE_BCURVE  BCurve name vis layer divs A:... / * degree { cp1 cp2 ... cpN }

            name = tokens{2};
            vis  = str2double(tokens{4});

            % Extract degree: after "/ *" there's the degree integer
            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);
            % after_slash: * degree
            degree = str2double(after_slash{2});

            % Extract control point names from { ... }
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
        % PARSE_CONIC  Conic name vis layer divs / * type center radius_pt apex_pt angle_start angle_end
        %
        %   Type 1 = elliptical arc.  Center, radius-point (defines a semi-axis),
        %   apex-point (defines the other semi-axis), and angular sweep.

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);
            % after_slash: * type center radius_pt apex_pt angle_start angle_end

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
        % PARSE_COPY_CURVE  CopyCurve name vis layer divs / * source src_pt dst_pt sx sy sz

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);
            % after_slash: * source src_pt dst_pt sx sy sz

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
        % PARSE_LINE_ENTITY  Line name vis layer divs / * pt_start pt_end

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
        % PARSE_BSUB_CURVE  BSubCurve name vis layer divs A:... / * degree { bead1 bead2 }

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
        % PARSE_ARC  Arc name vis layer divs / * type p1 p2 p3
        %
        %   Type 2: p1 = start point, p2 = centre, p3 = end point.
        %   The arc sweeps from p1 to p3 around p2 (shortest path).

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
        % PARSE_POLYCURVE2  PolyCurve2 name vis layer divs ... / * { c1 c2 ... cN }
        %
        %   Joins N sub-curves end-to-end.  Each sub-curve occupies an
        %   equal fraction of the parameter space [0, 1].

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
        % PARSE_PROJ_CURVE  ProjCurve name vis layer divs / * source_curve *plane
        %
        %   Projects source_curve onto a coordinate plane.
        %   *Y=0 → set y=0,  *X=0 → set x=0,  *Z=0 → set z=0.

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

        %% ─── SNAKE / BEAD PARSERS ───────────────────────────────

        function parse_edge_snake(obj, tokens, ~)
        % PARSE_EDGE_SNAKE  EdgeSnake name vis layer divs / * edge_index surface_name

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
        % PARSE_ABS_BEAD  AbsBead name vis layer / snake_or_curve parameter

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
        % PARSE_ABS_RING  AbsRing name vis layer / parent_snake parameter
        %
        %   An AbsRing is a parametric point placed at a fixed parameter value
        %   on a snake (typically an EdgeSnake).  It is identical in structure
        %   to AbsBead but its parent is always a snake, never a bare curve.
        %
        %   MultiSurf format:
        %     AbsRing ring1 9 -1 / snake1 0.0 ;
        %
        %   The entity resolves to a 3D point via eval_any_point, which calls
        %   eval_snake(parent_snake, parameter).

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
        % PARSE_BSUB_SNAKE  BSubSnake name vis layer divs A:... / * degree { bead1 bead2 }

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

        %% ─── SURFACE PARSERS ────────────────────────────────────

        function parse_ruled_surf(obj, tokens, ~)
        % PARSE_RULED_SURF  RuledSurf name vis layer divsU divsV sym / * curve1 curve2

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
        % PARSE_REV_SURF  RevSurf name vis layer divsU divsV sym / * profile axis angle_start angle_end

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
        % PARSE_BLOFT_SURF  BLoftSurf name vis layer divsU divsV sym A:... / * degree { section1 ... sectionN }

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
        % PARSE_DEV_SURF  DevSurf name vis layer divsU divsV sym / snake curve

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
        % PARSE_MIRR_SURF  MirrSurf name vis layer divsU divsV flag / source *PLANE

            name = tokens{2};
            vis  = str2double(tokens{4});

            slash_idx = find(strcmp(tokens, '/'));
            after_slash = tokens(slash_idx+1 : end);

            source_name = after_slash{1};
            plane_token = after_slash{2};  % e.g. '*Y=0' or '*X=0'

            % Parse mirror plane
            if contains(plane_token, 'Y=0')
                mirror_plane = 'Y';
            elseif contains(plane_token, 'X=0')
                mirror_plane = 'X';
            else
                mirror_plane = 'Y';  % default
                warning('WEC_MS2_Parser:UnknownPlane', ...
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


    %% ═══════════════════════════════════════════════════════════════
    %%  §3  POINT EVALUATOR
    %% ═══════════════════════════════════════════════════════════════

    methods

        function pt = eval_point(obj, name)
        % EVAL_POINT  Evaluate a FramePoint by resolving its parent chain.
        %
        %   pt = model.eval_point('pt3')   → [x, y, z]
        %
        %   FramePoints define position as parent + offset.  When the
        %   parent is another FramePoint, recursion resolves the chain
        %   down to a root point (both parents = '*').
        %
        %   CACHING: Results are memoized in point_cache.  FramePoint
        %   positions never change unless control points are modified
        %   (geometry scaling).  Call clear_cache() after any modification.

            % ── Cache lookup ──────────────────────────────────────
            if obj.point_cache.isKey(name)
                pt = obj.point_cache(name);
                return;
            end

            if ~obj.entities.isKey(name)
                error('WEC_MS2_Parser:EntityNotFound', ...
                       'Entity not found: %s', name);
            end

            e = obj.entities(name);

            if strcmp(e.type, 'FramePoint')
                p = e.params;

                % Resolve parent position (recursive)
                if ~isempty(p.parent1)
                    base = obj.eval_any_point(p.parent1);
                elseif ~isempty(p.parent2)
                    base = obj.eval_any_point(p.parent2);
                else
                    base = [0, 0, 0];  % absolute position
                end

                pt = base + p.offset;

            elseif strcmp(e.type, 'MirrPoint')
                pt = obj.eval_any_point(e.params.source);
                switch e.params.plane
                    case 'X', pt(1) = -pt(1);
                    case 'Y', pt(2) = -pt(2);
                    case 'Z', pt(3) = -pt(3);
                end

            else
                error('WEC_MS2_Parser:WrongType', ...
                       '%s is %s, not FramePoint/MirrPoint', name, e.type);
            end

            % ── Store in cache ────────────────────────────────────
            obj.point_cache(name) = pt;
        end

        function pt = eval_any_point(obj, name)
        % EVAL_ANY_POINT  Resolve any entity to a 3D point.
        %
        %   pt = model.eval_any_point('bead1')  → [x, y, z]
        %
        %   Handles:
        %     FramePoint → recursive parent resolution
        %     AbsBead    → evaluate bead on its parent curve/snake
        %     Other      → error
        %
        %   WHY this exists:
        %     Line entities can reference AbsBeads as endpoints (e.g.
        %     C0.ms2 line1 goes from bead1 to pt12).  eval_point only
        %     handles FramePoints.  This method dispatches correctly.

            if ~obj.entities.isKey(name)
                error('WEC_MS2_Parser:EntityNotFound', ...
                       'Entity not found: %s', name);
            end

            e = obj.entities(name);

            switch e.type
                case 'FramePoint'
                    pt = obj.eval_point(name);

                case 'MirrPoint'
                    pt = obj.eval_point(name);

                case 'AbsBead'
                    parent_name = e.params.parent_curve;
                    t_val = e.params.parameter;
                    % Try snake first, then curve
                    try
                        pt_arr = obj.eval_snake(parent_name, t_val);
                    catch
                        pt_arr = obj.eval_curve(parent_name, t_val);
                    end
                    pt = pt_arr(1, :);  % ensure [1×3]

                case 'AbsRing'
                    % AbsRing: fixed-parameter point on a snake.
                    % Structurally identical to AbsBead evaluation;
                    % parent is always a snake entity (e.g. EdgeSnake).
                    parent_name = e.params.parent_snake;
                    t_val = e.params.parameter;
                    try
                        pt_arr = obj.eval_snake(parent_name, t_val);
                    catch
                        pt_arr = obj.eval_curve(parent_name, t_val);
                    end
                    pt = pt_arr(1, :);

                otherwise
                    error('WEC_MS2_Parser:CannotResolvePoint', ...
                           'Cannot resolve %s (type: %s) to a point', ...
                           name, e.type);
            end
        end

    %% ═══════════════════════════════════════════════════════════════
    %%  §4  CURVE EVALUATORS
    %% ═══════════════════════════════════════════════════════════════

        function pts = eval_curve_or_snake(obj, name, t)
        % EVAL_CURVE_OR_SNAKE  Evaluate any curve-like entity at parameter t.
        %
        %   pts = model.eval_curve_or_snake('snake4', 0.5)  → [1 × 3]
        %
        %   WHY this exists:
        %     BLoftSurf sections can be either curves (BCurve, CopyCurve)
        %     or snakes (EdgeSnake, BSubSnake).  Rather than relying on a
        %     try-catch fallback, this method dispatches cleanly based on
        %     entity type.

            if ~obj.entities.isKey(name)
                error('WEC_MS2_Parser:EntityNotFound', ...
                       'Entity not found: %s', name);
            end

            e = obj.entities(name);

            % Snake types → eval_snake
            if any(strcmp(e.type, {'EdgeSnake', 'BSubSnake', 'AbsBead'}))
                pts = obj.eval_snake(name, t);
            else
                pts = obj.eval_curve(name, t);
            end
        end

        function pts = eval_curve(obj, name, t)
        % EVAL_CURVE  Evaluate a curve entity at parameter values t.
        %
        %   pts = model.eval_curve('curve1', linspace(0,1,100))
        %
        %   Returns [N × 3] array of (x,y,z) points.
        %   t must be in [0, 1].

            t = t(:);  % force column
            N = length(t);

            if ~obj.entities.isKey(name)
                error('WEC_MS2_Parser:EntityNotFound', ...
                       'Entity not found: %s', name);
            end

            e = obj.entities(name);

            switch e.type

                case 'BCurve'
                    pts = obj.eval_bcurve(e, t);

                case 'Conic'
                    pts = obj.eval_conic(e, t);

                case 'CopyCurve'
                    pts = obj.eval_copy_curve(e, t);

                case 'Line'
                    pts = obj.eval_line(e, t);

                case 'BSubCurve'
                    pts = obj.eval_bsub_curve(e, t);

                case 'Arc'
                    pts = obj.eval_arc(e, t);

                case 'PolyCurve2'
                    pts = obj.eval_polycurve2(e, t);

                case 'ProjCurve'
                    pts = obj.eval_proj_curve(e, t);

                otherwise
                    % In MultiSurf, snakes (EdgeSnake, BSubSnake) are
                    % curve-on-surface entities used interchangeably with
                    % curves in surface definitions.  Dispatch transparently.
                    if any(strcmp(e.type, {'EdgeSnake', 'BSubSnake', 'AbsBead'}))
                        pts = obj.eval_snake(name, t);
                    else
                        error('WEC_MS2_Parser:UnsupportedCurveType', ...
                               'Cannot evaluate %s as curve or snake (type: %s)', ...
                               name, e.type);
                    end
            end
        end

        function pts = eval_bcurve(obj, e, t)
        % EVAL_BCURVE  Evaluate a B-spline curve at parameters t.
        %
        %   Uses degree from file (typically 2) with clamped knot vector.
        %   Control points are resolved by evaluating their FramePoint
        %   parents.
        %
        %   WHY degree-2 B-spline approximates (not interpolates)?
        %     With a clamped knot vector, the curve passes through the
        %     first and last control points but only APPROXIMATES the
        %     interior points.  This matches MultiSurf's semantics.

            cp_names = e.params.ctrl_pt_names;
            n_cp = length(cp_names);
            degree = e.params.degree;

            % Resolve control points
            ctrl_pts = zeros(n_cp, 3);
            for i = 1:n_cp
                ctrl_pts(i, :) = obj.eval_point(cp_names{i});
            end

            % Generate clamped knot vector
            knots = WEC_MS2_Parser.make_clamped_knots(n_cp, degree);

            % Evaluate B-spline
            pts = WEC_MS2_Parser.bspline_curve_eval(knots, ctrl_pts, degree, t);
        end

        function pts = eval_conic(obj, e, t)
        % EVAL_CONIC  Evaluate an elliptical arc at parameters t.
        %
        %   The conic is parameterized as:
        %     P(t) = center + a·cos(θ)·ê_a + b·sin(θ)·ê_b
        %   where:
        %     a = |radius_pt − center|, ê_a = unit(radius_pt − center)
        %     b = |apex_pt − center|,   ê_b = unit(apex_pt − center)
        %     θ = angle_start + t·(angle_end − angle_start)
        %
        %   eval_any_point is used (not eval_point) so that control points
        %   can be AbsRing or AbsBead entities, not just FramePoints.

            p = e.params;
            center = obj.eval_any_point(p.center);
            rad_pt = obj.eval_any_point(p.radius_pt);
            apex   = obj.eval_any_point(p.apex_pt);

            % Semi-axes
            vec_a = rad_pt - center;
            vec_b = apex - center;
            a = norm(vec_a);
            b = norm(vec_b);
            e_a = vec_a / a;
            e_b = vec_b / b;

            % Map t ∈ [0,1] to angle ∈ [angle_start, angle_end] (degrees)
            theta = deg2rad(p.angle_start + t * (p.angle_end - p.angle_start));

            pts = center + a * cos(theta) .* e_a + b * sin(theta) .* e_b;
        end

        function pts = eval_copy_curve(obj, e, t)
        % EVAL_COPY_CURVE  Evaluate a translated copy of a source curve.
        %
        %   offset = dst_pt − src_pt.  Each source point is shifted by offset.
        %   Scale factors are applied (typically all 1.0 = pure translation).

            p = e.params;
            src_pos = obj.eval_point(p.src_pt);
            dst_pos = obj.eval_point(p.dst_pt);
            offset  = dst_pos - src_pos;

            % Evaluate source curve
            src_pts = obj.eval_curve(p.source, t);

            % Apply translation (scale is relative to src_pt)
            pts = zeros(size(src_pts));
            for i = 1:size(src_pts, 1)
                rel = src_pts(i, :) - src_pos;
                pts(i, :) = src_pos + rel .* p.scale + offset;
            end
        end

        function pts = eval_line(obj, e, t)
        % EVAL_LINE  Evaluate a line segment: P(t) = (1−t)·start + t·end.
        %
        %   Uses eval_any_point for endpoints because Line entities can
        %   reference AbsBeads (e.g. C0.ms2: line1 from bead1 to pt12).

            p       = e.params;
            p_start = obj.eval_any_point(p.pt_start);
            p_end   = obj.eval_any_point(p.pt_end);

            pts = (1 - t) .* p_start + t .* p_end;
        end

        function pts = eval_bsub_curve(obj, e, t)
        % EVAL_BSUB_CURVE  Evaluate a sub-section of a parent curve.
        %
        %   The sub-curve is defined by bead positions on the parent.
        %   t ∈ [0,1] maps to [bead1_param, bead2_param] on the parent.

            p = e.params;
            bead_names = p.bead_names;

            % Get bead parameter values on the parent curve
            bead1 = obj.entities(bead_names{1});
            bead2 = obj.entities(bead_names{end});
            t_start = bead1.params.parameter;
            t_end   = bead2.params.parameter;

            % Map t to parent parameter space
            t_parent = t_start + t * (t_end - t_start);

            % Evaluate parent curve
            parent_name = bead1.params.parent_curve;
            pts = obj.eval_curve(parent_name, t_parent);
        end


        function pts = eval_arc(obj, e, t)
        % EVAL_ARC  Circular arc: start → end around centre.
        %
        %   P(t) = centre + r (cos(θ) e_r + sin(θ) e_t)
        %   where θ = t * θ_total and e_r points to the start.

            p = e.params;
            p_start  = obj.eval_any_point(p.pt_start);
            p_centre = obj.eval_any_point(p.pt_center);
            p_end    = obj.eval_any_point(p.pt_end);

            [pts, ~] = WEC_MS2_Parser.arc_evaluate(p_start, p_centre, p_end, t);
        end


        function pts = eval_polycurve2(obj, e, t)
        % EVAL_POLYCURVE2  Concatenation of N curves, equal parameter shares.
        %
        %   t ∈ [0, 1/N) → curve1(N*t)
        %   t ∈ [1/N, 2/N) → curve2(N*t − 1)
        %   etc.

            cnames = e.params.curve_names;
            nc = length(cnames);
            pts = zeros(length(t), 3);

            for k = 1:length(t)
                % FIX (L5): Handle t=1 exactly — return last curve endpoint.
                % The old clamp (1-1e-15) relied on the last sub-curve
                % implicitly handling near-1 values; bspline_curve_eval
                % handles t>=1-1e-10 explicitly, but PolyCurve2 did not.
                if t(k) >= 1 - 1e-10
                    p = obj.eval_curve(cnames{nc}, 1);
                    pts(k, :) = p(1, :);
                    continue;
                end
                tk = max(t(k), 0);
                seg = min(floor(tk * nc) + 1, nc);     % segment index 1..nc
                t_local = tk * nc - (seg - 1);          % local parameter [0, 1)
                t_local = min(max(t_local, 0), 1);
                p = obj.eval_curve(cnames{seg}, t_local);
                pts(k, :) = p(1, :);
            end
        end


        function pts = eval_proj_curve(obj, e, t)
        % EVAL_PROJ_CURVE  Project source curve onto a coordinate plane.
        %
        %   *Y=0 → set y=0,  *X=0 → set x=0,  *Z=0 → set z=0.

            p = e.params;

            % Source may be a curve or a snake (EdgeSnake)
            if obj.entities.isKey(p.source)
                se = obj.entities(p.source);
                if any(strcmp(se.type, {'EdgeSnake', 'BSubSnake', 'AbsBead'}))
                    pts = obj.eval_snake(p.source, t);
                else
                    pts = obj.eval_curve(p.source, t);
                end
            else
                pts = obj.eval_curve(p.source, t);
            end

            switch p.proj_plane
                case 'X', pts(:,1) = 0;
                case 'Y', pts(:,2) = 0;
                case 'Z', pts(:,3) = 0;
            end
        end


        function [pts, dpts] = eval_arc_deriv(obj, e, t)
        % EVAL_ARC_DERIV  Arc point and tangent.

            p = e.params;
            p_start  = obj.eval_any_point(p.pt_start);
            p_centre = obj.eval_any_point(p.pt_center);
            p_end    = obj.eval_any_point(p.pt_end);

            [pts, dpts] = WEC_MS2_Parser.arc_evaluate(p_start, p_centre, p_end, t);
        end


        function [pts, dpts] = eval_polycurve2_deriv(obj, e, t)
        % EVAL_POLYCURVE2_DERIV  Concatenated curve derivative.
        %   dP/dt_global = dP/dt_local × (dt_local/dt_global) = dP/dt_local × N

            cnames = e.params.curve_names;
            nc = length(cnames);
            pts  = zeros(length(t), 3);
            dpts = zeros(length(t), 3);

            for k = 1:length(t)
                % FIX (L5 companion): same endpoint handling as eval_polycurve2
                if t(k) >= 1 - 1e-10
                    [p, dp] = obj.eval_curve_with_deriv(cnames{nc}, 1);
                    pts(k, :)  = p(1, :);
                    dpts(k, :) = dp(1, :) * nc;
                    continue;
                end
                tk = max(t(k), 0);
                seg = min(floor(tk * nc) + 1, nc);
                t_local = tk * nc - (seg - 1);
                t_local = min(max(t_local, 0), 1);

                [p, dp] = obj.eval_curve_with_deriv(cnames{seg}, t_local);
                pts(k, :)  = p(1, :);
                dpts(k, :) = dp(1, :) * nc;   % chain rule: dt_local/dt = nc
            end
        end


        function [pts, dpts] = eval_proj_curve_deriv(obj, e, t)
        % EVAL_PROJ_CURVE_DERIV  Projection zeroes one component of the derivative too.

            p = e.params;

            if obj.entities.isKey(p.source)
                se = obj.entities(p.source);
                if any(strcmp(se.type, {'EdgeSnake', 'BSubSnake', 'AbsBead'}))
                    [pts, dpts] = obj.eval_snake_with_deriv(p.source, t);
                else
                    [pts, dpts] = obj.eval_curve_with_deriv(p.source, t);
                end
            else
                [pts, dpts] = obj.eval_curve_with_deriv(p.source, t);
            end

            switch p.proj_plane
                case 'X', pts(:,1) = 0; dpts(:,1) = 0;
                case 'Y', pts(:,2) = 0; dpts(:,2) = 0;
                case 'Z', pts(:,3) = 0; dpts(:,3) = 0;
            end
        end


    %% ═══════════════════════════════════════════════════════════════
    %%  §5  SNAKE EVALUATORS
    %% ═══════════════════════════════════════════════════════════════

        function pts = eval_snake(obj, name, t)
        % EVAL_SNAKE  Evaluate a snake (curve-on-surface) at parameters t.
        %
        %   pts = model.eval_snake('snake1', linspace(0,1,50))
        %
        %   Returns [N × 3] array.

            t = t(:);

            if ~obj.entities.isKey(name)
                error('WEC_MS2_Parser:EntityNotFound', ...
                       'Entity not found: %s', name);
            end

            e = obj.entities(name);

            switch e.type

                case 'EdgeSnake'
                    pts = obj.eval_edge_snake(e, t);

                case 'BSubSnake'
                    pts = obj.eval_bsub_snake(e, t);

                case 'AbsBead'
                    % A bead is a single point — evaluate its parent at the bead parameter
                    pt = obj.eval_snake(e.params.parent_curve, e.params.parameter);
                    pts = repmat(pt, length(t), 1);

                otherwise
                    % May be a curve — try curve evaluation
                    pts = obj.eval_curve(name, t);
            end
        end

        function pts = eval_edge_snake(obj, e, t)
        % EVAL_EDGE_SNAKE  Evaluate an edge of a parent surface.
        %
        %   EDGE NUMBERING (MultiSurf EdgeSnake convention):
        %     Edge 1: v = 0,   u varies 0→1
        %     Edge 2: u = 1,   v varies 0→1
        %     Edge 3: v = 1,   u varies 0→1
        %     Edge 4: u = 0,   v varies 0→1
        %
        %   WHY no reversal on edges 3 and 4?
        %     The boundary-loop convention (edges 3/4 reversed so that
        %     consecutive edges form a CCW loop) applies to topological
        %     analysis, NOT to EdgeSnake parameterization.  MultiSurf's
        %     EdgeSnake always parameterizes the edge in the NATURAL
        %     direction of the underlying surface parameter (u or v
        %     increasing from 0 to 1).
        %
        %     This was verified empirically on E1.ms2: snake4 = edge 4
        %     of surface3 (RevSurf 0→180°).  With reversal, snake4
        %     starts at φ=180° and the BLoftSurf self-intersects.
        %     Without reversal, snake4 starts at φ=0° and the loft
        %     produces the correct smooth transition.

            p = e.params;
            edge_idx = p.edge_index;
            surf_name = p.surface_name;

            N = length(t);
            pts = zeros(N, 3);

            switch edge_idx
                case 1  % v = 0, u = t
                    for i = 1:N
                        pts(i,:) = obj.eval_surface(surf_name, t(i), 0);
                    end
                case 2  % u = 1, v = t
                    for i = 1:N
                        pts(i,:) = obj.eval_surface(surf_name, 1, t(i));
                    end
                case 3  % v = 1, u = t (natural direction)
                    for i = 1:N
                        pts(i,:) = obj.eval_surface(surf_name, t(i), 1);
                    end
                case 4  % u = 0, v = t (natural direction)
                    for i = 1:N
                        pts(i,:) = obj.eval_surface(surf_name, 0, t(i));
                    end
                otherwise
                    error('WEC_MS2_Parser:BadEdge', ...
                           'Invalid edge index %d for %s', edge_idx, surf_name);
            end
        end

        function pts = eval_bsub_snake(obj, e, t)
        % EVAL_BSUB_SNAKE  Evaluate a sub-section of a parent snake.
        %
        %   The sub-snake is defined between two beads on the parent.
        %   t ∈ [0,1] maps linearly from bead_start to bead_end.
        %
        %   WHY no wrapping?
        %     When bead_end < bead_start, the snake traverses the parent
        %     curve BACKWARD (decreasing parameter), not forward-wrapping
        %     through the periodic boundary.  Example: on E1.ms2,
        %     snake2 = {bead2(0.75), bead1(0.25)} goes from the north
        %     pole backward through the x>0 equator to the south pole.
        %     Wrapping would go forward through the x<0 equator — wrong.
        %
        %     Verified: wrapping produces snake3 starting at (-1.5,0,-1)
        %     which cascades into the BLoftSurf connecting opposite sides.
        %     Linear backward gives (+1.5,0,-1) — correct.

            p = e.params;
            bead_names = p.bead_names;

            % Resolve bead parameters
            bead_start = obj.entities(bead_names{1});
            bead_end   = obj.entities(bead_names{end});
            t_start = bead_start.params.parameter;
            t_end   = bead_end.params.parameter;

            % Get parent curve/snake name from the first bead
            parent_name = bead_start.params.parent_curve;

            % Simple linear mapping — works for both forward and backward
            t_parent = t_start + t * (t_end - t_start);

            pts = obj.eval_snake(parent_name, t_parent);
        end


    %% ═══════════════════════════════════════════════════════════════
    %%  §6  SURFACE EVALUATORS
    %% ═══════════════════════════════════════════════════════════════

        function pt = eval_surface(obj, name, u, v)
        % EVAL_SURFACE  Evaluate a surface at parameter (u, v).
        %
        %   pt = model.eval_surface('surface3', 0.5, 0.3)  → [x, y, z]
        %
        %   u, v ∈ [0, 1].  Returns [1 × 3] point.
        %
        %   For grid evaluation, call in a loop or use eval_surface_grid.

            if ~obj.entities.isKey(name)
                error('WEC_MS2_Parser:EntityNotFound', ...
                       'Entity not found: %s', name);
            end

            e = obj.entities(name);

            switch e.type

                case 'RuledSurf'
                    pt = obj.eval_ruled_surf(e, u, v);

                case 'RevSurf'
                    pt = obj.eval_rev_surf(e, u, v);

                case 'BLoftSurf'
                    pt = obj.eval_bloft_surf(e, u, v);

                case 'DevSurf'
                    pt = obj.eval_dev_surf(e, u, v);

                case 'MirrSurf'
                    pt = obj.eval_mirr_surf(e, u, v);

                otherwise
                    error('WEC_MS2_Parser:UnsupportedSurface', ...
                           'Cannot evaluate %s as surface (type: %s)', ...
                           name, e.type);
            end
        end

        function S = eval_surface_grid(obj, name, u_grid, v_grid)
        % EVAL_SURFACE_GRID  Evaluate a surface on a (Nu × Nv) grid.
        %
        %   S = model.eval_surface_grid('surface3', ...
        %           linspace(0,1,20), linspace(0,1,20))
        %
        %   Returns [Nu × Nv × 3] array.

            Nu = length(u_grid);
            Nv = length(v_grid);
            S = zeros(Nu, Nv, 3);
            for i = 1:Nu
                for j = 1:Nv
                    S(i, j, :) = obj.eval_surface(name, u_grid(i), v_grid(j));
                end
            end
        end


    %% ═══════════════════════════════════════════════════════════════
    %%  §6b  ANALYTICAL DERIVATIVES
    %% ═══════════════════════════════════════════════════════════════
    %
    %  WHY analytical derivatives?
    %    The GL quadrature in WEC_HydroProperties computes
    %      n̂ dA = (∂S/∂u × ∂S/∂v) du dv
    %    Previously this used 5-point central finite differences:
    %      S_u ≈ (S(u+h,v) − S(u−h,v)) / 2h
    %    requiring 5 eval_surface calls per GL point (S, u±h, v±h).
    %    Analytical derivatives require 1 call to eval_surface_with_derivs,
    %    giving a 3-5× speedup on the GL quadrature.
    %
    %  CURVE DERIVATIVES
    %    BCurve:    B-spline derivative (degree reduction)
    %    Conic:     −a sin(θ) ê_a + b cos(θ) ê_b scaled by dθ/dt
    %    Line:      constant (end − start)
    %    CopyCurve: scale .* source derivative
    %    BSubCurve: parent derivative × dt_parent/dt
    %    Snakes:    chain rule through parent snake/surface
    %
    %  SURFACE DERIVATIVES
    %    RuledSurf:  S_u from curve derivs, S_v = C2(u) − C1(u)
    %    RevSurf:    S_u from profile deriv, S_v closed-form rotation
    %    BLoftSurf:  S_u from section derivs, S_v from loft B-spline deriv
    %    DevSurf:    same as RuledSurf
    %    MirrSurf:   flip coordinate of source derivs

        function [pts, dpts] = eval_curve_with_deriv(obj, name, t)
        % EVAL_CURVE_WITH_DERIV  Evaluate curve and tangent at parameter t.
        %
        %   [pts, dpts] = model.eval_curve_with_deriv('curve1', 0.3)
        %
        %   Returns [N × 3] pts and [N × 3] dpts = dP/dt.

            t = t(:);

            if ~obj.entities.isKey(name)
                error('WEC_MS2_Parser:EntityNotFound', ...
                       'Entity not found: %s', name);
            end

            e = obj.entities(name);

            switch e.type
                case 'BCurve'
                    [pts, dpts] = obj.eval_bcurve_deriv(e, t);

                case 'Conic'
                    [pts, dpts] = obj.eval_conic_deriv(e, t);

                case 'CopyCurve'
                    [pts, dpts] = obj.eval_copy_curve_deriv(e, t);

                case 'Line'
                    [pts, dpts] = obj.eval_line_deriv(e, t);

                case 'BSubCurve'
                    [pts, dpts] = obj.eval_bsub_curve_deriv(e, t);

                case 'Arc'
                    [pts, dpts] = obj.eval_arc_deriv(e, t);

                case 'PolyCurve2'
                    [pts, dpts] = obj.eval_polycurve2_deriv(e, t);

                case 'ProjCurve'
                    [pts, dpts] = obj.eval_proj_curve_deriv(e, t);

                otherwise
                    % Snake types — use eval_snake_with_deriv
                    if any(strcmp(e.type, {'EdgeSnake', 'BSubSnake', 'AbsBead'}))
                        [pts, dpts] = obj.eval_snake_with_deriv(name, t);
                    else
                        % Fallback to finite differences
                        h = 1e-7;
                        pts = obj.eval_curve(name, t);
                        pts_h = obj.eval_curve(name, min(t+h, 1));
                        pts_l = obj.eval_curve(name, max(t-h, 0));
                        dpts = (pts_h - pts_l) ./ (min(t+h,1) - max(t-h,0));
                    end
            end
        end


        function [pts, dpts] = eval_bcurve_deriv(obj, e, t)
        % EVAL_BCURVE_DERIV  B-spline curve and derivative.

            cp_names = e.params.ctrl_pt_names;
            n_cp = length(cp_names);
            degree = e.params.degree;

            ctrl_pts = zeros(n_cp, 3);
            for i = 1:n_cp
                ctrl_pts(i, :) = obj.eval_point(cp_names{i});
            end

            knots = WEC_MS2_Parser.make_clamped_knots(n_cp, degree);
            [pts, dpts] = WEC_MS2_Parser.bspline_curve_eval_with_deriv( ...
                knots, ctrl_pts, degree, t);
        end


        function [pts, dpts] = eval_conic_deriv(obj, e, t)
        % EVAL_CONIC_DERIV  Elliptical arc and derivative.
        %
        %   P(t) = center + a cos(θ) ê_a + b sin(θ) ê_b
        %   dP/dt = (−a sin(θ) ê_a + b cos(θ) ê_b) × dθ/dt
        %
        %   eval_any_point is used so AbsRing/AbsBead control points resolve.

            p = e.params;
            center = obj.eval_any_point(p.center);
            rad_pt = obj.eval_any_point(p.radius_pt);
            apex   = obj.eval_any_point(p.apex_pt);

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
        % EVAL_COPY_CURVE_DERIV  Translated copy: derivative = scale .* source derivative.

            p = e.params;
            src_pos = obj.eval_point(p.src_pt);
            dst_pos = obj.eval_point(p.dst_pt);
            offset  = dst_pos - src_pos;

            [src_pts, src_dpts] = obj.eval_curve_with_deriv(p.source, t);

            pts = zeros(size(src_pts));
            dpts = zeros(size(src_dpts));
            for i = 1:size(src_pts, 1)
                rel = src_pts(i, :) - src_pos;
                pts(i, :)  = src_pos + rel .* p.scale + offset;
                dpts(i, :) = src_dpts(i, :) .* p.scale;
            end
        end


        function [pts, dpts] = eval_line_deriv(obj, e, t)
        % EVAL_LINE_DERIV  Line segment: dP/dt = end − start (constant).

            p_start = obj.eval_any_point(e.params.pt_start);
            p_end   = obj.eval_any_point(e.params.pt_end);

            pts  = (1 - t) .* p_start + t .* p_end;
            dpts = repmat(p_end - p_start, length(t), 1);
        end


        function [pts, dpts] = eval_bsub_curve_deriv(obj, e, t)
        % EVAL_BSUB_CURVE_DERIV  Sub-curve derivative via chain rule.
        %   dP/dt = dP/dt_parent × dt_parent/dt
        %   where dt_parent/dt = (t_end − t_start)

            p = e.params;
            bead1 = obj.entities(p.bead_names{1});
            bead2 = obj.entities(p.bead_names{end});
            t_start = bead1.params.parameter;
            t_end   = bead2.params.parameter;

            t_parent = t_start + t * (t_end - t_start);
            parent_name = bead1.params.parent_curve;

            [pts, dpts_parent] = obj.eval_curve_with_deriv(parent_name, t_parent);
            dpts = dpts_parent * (t_end - t_start);
        end


        function [pts, dpts] = eval_snake_with_deriv(obj, name, t)
        % EVAL_SNAKE_WITH_DERIV  Snake derivative via chain rule.

            t = t(:);
            e = obj.entities(name);

            switch e.type
                case 'EdgeSnake'
                    [pts, dpts] = obj.eval_edge_snake_deriv(e, t);

                case 'BSubSnake'
                    [pts, dpts] = obj.eval_bsub_snake_deriv(e, t);

                otherwise
                    % Fallback to finite differences
                    h = 1e-7;
                    pts = obj.eval_snake(name, t);
                    pts_h = obj.eval_snake(name, min(t+h, 1));
                    pts_l = obj.eval_snake(name, max(t-h, 0));
                    dpts = (pts_h - pts_l) ./ (min(t+h,1) - max(t-h,0));
            end
        end


        function [pts, dpts] = eval_edge_snake_deriv(obj, e, t)
        % EVAL_EDGE_SNAKE_DERIV  Edge snake derivative.
        %
        %   Edge 1 (v=0, u=t): dP/dt = S_u(t, 0)
        %   Edge 2 (u=1, v=t): dP/dt = S_v(1, t)
        %   Edge 3 (v=1, u=t): dP/dt = S_u(t, 1)
        %   Edge 4 (u=0, v=t): dP/dt = S_v(0, t)

            p = e.params;
            N = length(t);
            pts  = zeros(N, 3);
            dpts = zeros(N, 3);

            for i = 1:N
                switch p.edge_index
                    case 1
                        [S, Su, ~] = obj.eval_surface_with_derivs(p.surface_name, t(i), 0);
                        pts(i,:) = S; dpts(i,:) = Su;
                    case 2
                        [S, ~, Sv] = obj.eval_surface_with_derivs(p.surface_name, 1, t(i));
                        pts(i,:) = S; dpts(i,:) = Sv;
                    case 3
                        [S, Su, ~] = obj.eval_surface_with_derivs(p.surface_name, t(i), 1);
                        pts(i,:) = S; dpts(i,:) = Su;
                    case 4
                        [S, ~, Sv] = obj.eval_surface_with_derivs(p.surface_name, 0, t(i));
                        pts(i,:) = S; dpts(i,:) = Sv;
                end
            end
        end


        function [pts, dpts] = eval_bsub_snake_deriv(obj, e, t)
        % EVAL_BSUB_SNAKE_DERIV  Sub-snake derivative via chain rule.

            p = e.params;
            bead_start = obj.entities(p.bead_names{1});
            bead_end   = obj.entities(p.bead_names{end});
            t_start = bead_start.params.parameter;
            t_end   = bead_end.params.parameter;

            t_parent = t_start + t * (t_end - t_start);
            parent_name = bead_start.params.parent_curve;

            [pts, dpts_parent] = obj.eval_snake_with_deriv(parent_name, t_parent);
            dpts = dpts_parent * (t_end - t_start);
        end


        function [S, Su, Sv] = eval_surface_with_derivs(obj, name, u, v)
        % EVAL_SURFACE_WITH_DERIVS  Surface point and partial derivatives.
        %
        %   [S, Su, Sv] = model.eval_surface_with_derivs('surface2', 0.5, 0.3)
        %
        %   Returns [1 × 3] arrays:
        %     S  — surface point
        %     Su — ∂S/∂u
        %     Sv — ∂S/∂v

            if ~obj.entities.isKey(name)
                error('WEC_MS2_Parser:EntityNotFound', ...
                       'Entity not found: %s', name);
            end

            e = obj.entities(name);

            switch e.type
                case 'RuledSurf'
                    [S, Su, Sv] = obj.eval_ruled_surf_derivs(e, u, v);
                case 'RevSurf'
                    [S, Su, Sv] = obj.eval_rev_surf_derivs(e, u, v);
                case 'BLoftSurf'
                    [S, Su, Sv] = obj.eval_bloft_surf_derivs(e, u, v);
                case 'DevSurf'
                    [S, Su, Sv] = obj.eval_dev_surf_derivs(e, u, v);
                case 'MirrSurf'
                    [S, Su, Sv] = obj.eval_mirr_surf_derivs(e, u, v);
                otherwise
                    % Fallback to finite differences
                    h = 1e-6;
                    S = obj.eval_surface(name, u, v);
                    Su = (obj.eval_surface(name, min(u+h,1), v) - ...
                          obj.eval_surface(name, max(u-h,0), v)) / ...
                         (min(u+h,1) - max(u-h,0));
                    Sv = (obj.eval_surface(name, u, min(v+h,1)) - ...
                          obj.eval_surface(name, u, max(v-h,0))) / ...
                         (min(v+h,1) - max(v-h,0));
            end
        end


        function [S, Su, Sv] = eval_ruled_surf_derivs(obj, e, u, v)
        % S(u,v) = (1−v)·C1(u) + v·C2(u)
        % Su = (1−v)·C1'(u) + v·C2'(u)
        % Sv = C2(u) − C1(u)

            p = e.params;
            [c1, dc1] = obj.eval_curve_with_deriv(p.curve1, u);
            [c2, dc2] = obj.eval_curve_with_deriv(p.curve2, u);

            S  = (1 - v) * c1 + v * c2;
            Su = (1 - v) * dc1 + v * dc2;
            Sv = c2 - c1;
        end


        function [S, Su, Sv] = eval_rev_surf_derivs(obj, e, u, v)
        % RevSurf analytical derivatives.
        %
        % S(u,v) = proj + r cos(φ) ê_r + r sin(φ) ê_t
        % where proj = axis_start + z_along * axis_dir
        %
        % Su = dproj/du + dr/du cos(φ) ê_r + dr/du sin(φ) ê_t
        %      + r cos(φ) dê_r/du + r sin(φ) dê_t/du
        %    (simplifies because ê_r, ê_t change direction with the profile)
        %
        % Sv = r (−sin(φ) ê_r + cos(φ) ê_t) dφ/dv

            p = e.params;

            % Profile point and derivative
            [profile_pt, profile_dpdt] = obj.eval_curve_with_deriv(p.profile, u);
            profile_pt = profile_pt(1,:);
            profile_dpdt = profile_dpdt(1,:);

            % Cached axis
            if obj.rev_axis_cache.isKey(p.axis)
                ax = obj.rev_axis_cache(p.axis);
                axis_start = ax.start;
                axis_dir   = ax.dir;
            else
                axis_ent = obj.entities(p.axis);
                axis_start = obj.eval_any_point(axis_ent.params.pt_start);
                axis_end   = obj.eval_any_point(axis_ent.params.pt_end);
                axis_vec = axis_end - axis_start;
                axis_dir = axis_vec / norm(axis_vec);
            end

            % Decompose profile point into axial + radial
            v_rel = profile_pt - axis_start;
            z_along = dot(v_rel, axis_dir);
            proj = axis_start + z_along * axis_dir;
            radial = profile_pt - proj;
            r = norm(radial);

            % Revolution angle
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

            % S
            S = proj + r * cos(phi) * e_r + r * sin(phi) * e_t;

            % Sv = r (-sin(φ) ê_r + cos(φ) ê_t) dφ/dv
            Sv = r * (-sin(phi) * e_r + cos(phi) * e_t) * dphi_dv;

            % Su: derivative of profile point projected into the rotated frame
            % dproj/du = (d(profile)/du · axis_dir) * axis_dir
            % dr/du and dê_r/du from the radial component of d(profile)/du
            dz_du = dot(profile_dpdt, axis_dir);
            dproj_du = dz_du * axis_dir;
            dradial_du = profile_dpdt - dproj_du;
            dr_du = dot(dradial_du, e_r);

            % The radial direction also rotates as the profile changes
            % dê_r/du = (dradial_du − dr_du ê_r) / r
            if r > 1e-10
                de_r_du = (dradial_du - dr_du * e_r) / r;
                de_t_du = cross(axis_dir, de_r_du);
                % FIX (L4): Removed rescaling de_t_du/|de_t_du|*|de_r_du|.
                % Correct: d/du[cross(axis,e_r)] = cross(axis, de_r_du) — no scaling.
                % Was a no-op for vertical axes; wrong for non-vertical.
            else
                de_r_du = [0, 0, 0];
                de_t_du = [0, 0, 0];
            end

            Su = dproj_du + dr_du * cos(phi) * e_r + r * cos(phi) * de_r_du ...
                          + dr_du * sin(phi) * e_t + r * sin(phi) * de_t_du;
        end


        function [S, Su, Sv] = eval_bloft_surf_derivs(obj, e, u, v)
        % BLoftSurf analytical derivatives.
        %
        % Su: evaluate section curve derivatives at u, B-spline loft in v
        % Sv: evaluate section points at u, B-spline loft derivative in v
        %
        % FIX (D2): Per-u section cache.
        %   In the GL quadrature loop the outer iteration is over u-nodes and
        %   the inner iteration is over v-nodes. This function is called once
        %   per (u,v) pair. Without caching, n_sec eval_curve_with_deriv calls
        %   fire for every v-node at a fixed u — i.e. n_quad × n_sec times per
        %   u-node. The cache stores the last (surface_key, u) evaluation so
        %   section evaluations only run on the FIRST v-call at each u.
        %   For n_sec=5, n_quad=20: reduces 100 section evals to 5 per u-node.
        %   Uses MATLAB persistent variables (function-scoped, survives calls).

            p = e.params;
            section_names = p.section_names;
            n_sec = length(section_names);
            degree = p.degree;

            % ── Per-u section cache ──────────────────────────────────────
            persistent bloft_cache_key bloft_cache_u bloft_cache_pts ...
                       bloft_cache_dpts bloft_cache_knots;

            surf_key = strjoin(section_names, '|');

            if ~isempty(bloft_cache_key) && ...
                    strcmp(bloft_cache_key, surf_key) && ...
                    abs(bloft_cache_u - u) < 1e-14
                % Cache hit — reuse section evaluations from previous v-call
                sec_pts  = bloft_cache_pts;
                sec_dpts = bloft_cache_dpts;
                knots_v  = bloft_cache_knots;
            else
                % Cache miss — evaluate all sections at this u
                sec_pts  = zeros(n_sec, 3);
                sec_dpts = zeros(n_sec, 3);
                for k = 1:n_sec
                    sname = section_names{k};
                    if ~obj.entities.isKey(sname)
                        sec_pts(k,:)  = [0 0 0];
                        sec_dpts(k,:) = [0 0 0];
                        continue;
                    end
                    se = obj.entities(sname);
                    if any(strcmp(se.type, {'EdgeSnake', 'BSubSnake', 'AbsBead'}))
                        [pp, dd] = obj.eval_snake_with_deriv(sname, u);
                    else
                        [pp, dd] = obj.eval_curve_with_deriv(sname, u);
                    end
                    sec_pts(k,:)  = pp(1,:);
                    sec_dpts(k,:) = dd(1,:);
                end
                knots_v = WEC_MS2_Parser.make_clamped_knots(n_sec, degree);

                % Store in cache
                bloft_cache_key   = surf_key;
                bloft_cache_u     = u;
                bloft_cache_pts   = sec_pts;
                bloft_cache_dpts  = sec_dpts;
                bloft_cache_knots = knots_v;
            end

            % S and Sv: loft through section_pts, derivative in v
            [S_arr, Sv_arr] = WEC_MS2_Parser.bspline_curve_eval_with_deriv( ...
                knots_v, sec_pts, degree, v);
            S  = S_arr(1,:);
            Sv = Sv_arr(1,:);

            % Su: loft through section u-derivatives at the same v
            Su_arr = WEC_MS2_Parser.bspline_curve_eval(knots_v, sec_dpts, degree, v);
            Su = Su_arr(1,:);
        end


        function [S, Su, Sv] = eval_dev_surf_derivs(obj, e, u, v)
        % DevSurf = ruled surface between snake and curve.
        % S(u,v) = (1−v)·snake(u) + v·curve(u)

            p = e.params;
            [s_pt, s_dpt] = obj.eval_snake_with_deriv(p.snake, u);
            [c_pt, c_dpt] = obj.eval_curve_with_deriv(p.curve, u);

            s_pt = s_pt(1,:); s_dpt = s_dpt(1,:);
            c_pt = c_pt(1,:); c_dpt = c_dpt(1,:);

            S  = (1 - v) * s_pt + v * c_pt;
            Su = (1 - v) * s_dpt + v * c_dpt;
            Sv = c_pt - s_pt;
        end


        function [S, Su, Sv] = eval_mirr_surf_derivs(obj, e, u, v)
        % MirrSurf: flip one coordinate of source surface derivatives.

            p = e.params;
            [S, Su, Sv] = obj.eval_surface_with_derivs(p.source, u, v);

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

        %% ─── RULED SURFACE ──────────────────────────────────────

        function pt = eval_ruled_surf(obj, e, u, v)
        % EVAL_RULED_SURF  S(u,v) = (1−v)·C₁(u) + v·C₂(u)
        %
        %   The simplest surface: linear interpolation between two
        %   boundary curves at each u-parameter.

            p = e.params;
            c1 = obj.eval_curve(p.curve1, u);
            c2 = obj.eval_curve(p.curve2, u);
            pt = (1 - v) * c1 + v * c2;
        end

        %% ─── REVOLUTION SURFACE ─────────────────────────────────

        function pt = eval_rev_surf(obj, e, u, v)
        % EVAL_REV_SURF  Revolve a profile curve around an axis.
        %
        %   u parameterizes the profile curve (generatrix).
        %   v parameterizes the revolution angle from angle_start to angle_end.
        %
        %   For a point P on the profile at distance r from the axis:
        %     S(u,v) = axis_origin + z·ê_axis + r·cos(φ)·ê_r + r·sin(φ)·ê_t
        %   where φ = angle_start + v·(angle_end − angle_start).
        %
        %   AXIS CACHING:
        %     The axis (a Line entity with two fixed endpoints) is evaluated
        %     once and cached in rev_axis_cache.  For C0, line1 is the only
        %     axis, called ~34,000 times per compute_submerged.  Caching
        %     eliminates all redundant eval_any_point calls.

            p = e.params;

            % Evaluate profile point at parameter u
            profile_pt = obj.eval_curve(p.profile, u);

            % ── Cached axis evaluation ────────────────────────────
            if obj.rev_axis_cache.isKey(p.axis)
                ax = obj.rev_axis_cache(p.axis);
                axis_start = ax.start;
                axis_dir   = ax.dir;
            else
                axis_ent = obj.entities(p.axis);
                axis_start = obj.eval_any_point(axis_ent.params.pt_start);
                axis_end   = obj.eval_any_point(axis_ent.params.pt_end);
                axis_vec  = axis_end - axis_start;
                axis_len  = norm(axis_vec);
                if axis_len < 1e-12
                    pt = profile_pt;
                    return;
                end
                axis_dir = axis_vec / axis_len;

                % Store in cache
                ax = struct('start', axis_start, 'end_pt', axis_end, ...
                            'dir', axis_dir, 'len', axis_len);
                obj.rev_axis_cache(p.axis) = ax;
            end

            % Project profile point onto axis
            v_rel = profile_pt - axis_start;
            z_along = dot(v_rel, axis_dir);
            proj = axis_start + z_along * axis_dir;

            % Radial vector (perpendicular to axis)
            radial = profile_pt - proj;
            r = norm(radial);

            if r < 1e-12
                % Point is on the axis — no rotation effect
                pt = profile_pt;
                return;
            end

            % Build local coordinate frame perpendicular to axis
            e_r = radial / r;  % initial radial direction
            e_t = cross(axis_dir, e_r);
            e_t = e_t / norm(e_t);  % tangential direction

            % Revolution angle
            phi = deg2rad(p.angle_start + v * (p.angle_end - p.angle_start));

            % Rotated point
            pt = proj + r * cos(phi) * e_r + r * sin(phi) * e_t;
        end

        %% ─── B-SPLINE LOFT SURFACE ─────────────────────────────

        function pt = eval_bloft_surf(obj, e, u, v)
        % EVAL_BLOFT_SURF  Degree-p B-spline loft through control sections.
        %
        %   u parameterizes along each section curve.
        %   v parameterizes across sections (the loft direction).
        %
        %   The sections are CONTROL SECTIONS for a degree-p B-spline
        %   in the v-direction.  The surface APPROXIMATES the sections
        %   (does not pass through them).  This is critical for smooth
        %   transitions — see the ellipsoid-to-column loft discussion.
        %
        %   WHY B-spline loft (not PCHIP, not CubicSpline)?
        %     PCHIP interpolates through sections → flat steps at
        %     transitions where two sections share the same z.
        %     CubicSpline overshoots at large shape changes.
        %     Degree-2 B-spline pulls the surface toward later sections,
        %     creating the smooth fillet observed in MultiSurf.

            p = e.params;
            section_names = p.section_names;
            n_sections = length(section_names);
            degree = p.degree;

            % Evaluate each section at parameter u → one 3D point per section
            section_pts = zeros(n_sections, 3);
            for k = 1:n_sections
                pts_k = obj.eval_curve_or_snake(section_names{k}, u);
                section_pts(k, :) = pts_k(1, :);  % ensure [1×3]
            end

            % Generate clamped knot vector for the loft direction
            knots_v = WEC_MS2_Parser.make_clamped_knots(n_sections, degree);

            % Evaluate B-spline at parameter v using section_pts as control points
            pt = WEC_MS2_Parser.bspline_curve_eval(knots_v, section_pts, degree, v);
        end

        function pts = eval_bloft_surf_at_u(obj, e, u, v_array)
        % EVAL_BLOFT_SURF_AT_U  Evaluate a BLoftSurf at fixed u, multiple v.
        %
        %   pts = model.eval_bloft_surf_at_u(e, 0.3, [0.1; 0.2; 0.5])
        %
        % ── NOTE (D2) ─────────────────────────────────────────────────────
        % The per-u redundancy this method aimed to solve is now handled by
        % the persistent cache inside eval_bloft_surf_derivs, which fires
        % transparently from the GL quadrature path. This method is still
        % available for POINT-ONLY (no-derivative) callers.
        % ──────────────────────────────────────────────────────────────────
        %
        %   Returns [length(v_array) × 3] array.
        %
        %   WHY this exists:
        %     In the GL quadrature loop, eval_bloft_surf(u, v) is called
        %     n_quad times at the same u (once per v-node).  Each call
        %     re-evaluates all section curves at parameter u.  This method
        %     evaluates sections ONCE at u, then runs the v-direction
        %     B-spline loft for each v value.
        %
        %     For E1 with 5 sections and n_quad=20: eliminates 19×5 = 95
        %     redundant section evaluations per u-node.

            p = e.params;
            section_names = p.section_names;
            n_sections = length(section_names);
            degree = p.degree;

            % ── Evaluate sections ONCE at u ───────────────────────
            section_pts = zeros(n_sections, 3);
            for k = 1:n_sections
                pts_k = obj.eval_curve_or_snake(section_names{k}, u);
                section_pts(k, :) = pts_k(1, :);
            end

            % ── Generate knot vector ONCE ─────────────────────────
            knots_v = WEC_MS2_Parser.make_clamped_knots(n_sections, degree);

            % ── Evaluate B-spline at each v ───────────────────────
            v_array = v_array(:);
            pts = WEC_MS2_Parser.bspline_curve_eval(knots_v, section_pts, degree, v_array);
        end

        %% ─── DEVELOPED SURFACE ──────────────────────────────────

        function pt = eval_dev_surf(obj, e, u, v)
        % EVAL_DEV_SURF  Developed (ruled) surface with arc-length matching.
        %
        %   S(u,v) = (1−v)·snake(u) + v·curve(u)
        %
        %   This is essentially a ruled surface between two boundary
        %   curves (a snake and a curve), where the parameter
        %   correspondence is by arc-length matching.
        %
        %   [SIMPLIFICATION] Currently uses direct parameter matching
        %   (same u on both curves).  Arc-length reparameterization can
        %   be added if validation shows discrepancies.

            p = e.params;
            pt_snake = obj.eval_snake(p.snake, u);
            pt_curve = obj.eval_curve(p.curve, u);

            pt = (1 - v) * pt_snake + v * pt_curve;
        end

        %% ─── MIRROR SURFACE ─────────────────────────────────────

        function pt = eval_mirr_surf(obj, e, u, v)
        % EVAL_MIRR_SURF  Mirror a source surface about X=0 or Y=0.
        %
        %   Simply evaluates the source surface and negates one coordinate.
        %   Normal orientation is automatically reversed by the negation
        %   (the winding order flips).

            p = e.params;
            pt = obj.eval_surface(p.source, u, v);

            switch p.mirror_plane
                case 'Y'
                    pt(2) = -pt(2);
                case 'X'
                    pt(1) = -pt(1);
            end
        end

    end % methods


    %% ═══════════════════════════════════════════════════════════════
    %%  §7  B-SPLINE ENGINE (Static)
    %% ═══════════════════════════════════════════════════════════════
    %
    %  Self-contained B-spline evaluation — no Curve Fitting Toolbox.
    %
    %  WHY implement from scratch?
    %    1. Portability: not all MATLAB installations have the toolbox.
    %    2. Transparency: the evaluation is 20 lines of Cox-de Boor,
    %       easy to verify and debug.
    %    3. Control: we need clamped knot vectors with specific degree,
    %       which toolbox functions sometimes obscure.

    methods (Static)

        function knots = make_clamped_knots(n_ctrl, degree)
        % MAKE_CLAMPED_KNOTS  Generate a clamped (open) knot vector.
        %
        %   knots = make_clamped_knots(n_ctrl, degree)
        %
        %   For n_ctrl control points and given degree:
        %     n_knots = n_ctrl + degree + 1
        %     First (degree+1) knots = 0
        %     Last  (degree+1) knots = 1
        %     Interior knots uniformly spaced in (0, 1)
        %
        %   The clamped knot vector ensures the curve passes through
        %   the first and last control points.

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
        % BSPLINE_CURVE_EVAL  Evaluate a B-spline curve at parameters t.
        %
        %   pts = bspline_curve_eval(knots, ctrl_pts, degree, t)
        %
        %   INPUTS
        %     knots    : [1 × (n_ctrl + degree + 1)]  knot vector
        %     ctrl_pts : [n_ctrl × 3]  control point coordinates
        %     degree   : integer  polynomial degree
        %     t        : [N × 1]  parameter values in [0, 1]
        %
        %   OUTPUT
        %     pts : [N × 3]  evaluated points
        %
        %   Uses the Cox-de Boor recursion for basis function evaluation.

            t = t(:);
            N = length(t);
            n_ctrl = size(ctrl_pts, 1);
            pts = zeros(N, 3);

            for k = 1:N
                % Clamp t to [0, 1−eps] to handle right endpoint
                tk = min(max(t(k), 0), 1 - 1e-12);

                % Evaluate all basis functions at tk
                basis = WEC_MS2_Parser.bspline_basis_all(knots, degree, tk, n_ctrl);

                % Weighted sum of control points
                pts(k, :) = basis * ctrl_pts;
            end

            % Handle exact t = 1 (endpoint interpolation for clamped spline)
            endpoint_mask = (t >= 1 - 1e-10);
            if any(endpoint_mask)
                pts(endpoint_mask, :) = repmat(ctrl_pts(end, :), sum(endpoint_mask), 1);
            end
        end

        function N = bspline_basis_all(knots, degree, t, n_ctrl)
        % BSPLINE_BASIS_ALL  Evaluate all B-spline basis functions at t.
        %
        %   N = bspline_basis_all(knots, degree, t, n_ctrl)
        %
        %   Returns [1 × n_ctrl] vector of basis function values.
        %   Uses the Cox-de Boor recursion formula:
        %
        %     N_{i,0}(t) = 1  if knots(i) ≤ t < knots(i+1), else 0
        %     N_{i,p}(t) = (t − knots(i)) / (knots(i+p) − knots(i)) · N_{i,p-1}(t)
        %                + (knots(i+p+1) − t) / (knots(i+p+1) − knots(i+1)) · N_{i+1,p-1}(t)
        %
        %   Convention: 0/0 = 0 (degenerate knot spans).

            n = length(knots) - 1;  % number of basis functions at degree 0

            % Degree 0
            N0 = zeros(1, n);
            for i = 1:n
                if knots(i) <= t && t < knots(i+1)
                    N0(i) = 1;
                end
            end

            % Recursion up to target degree
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
        % BSPLINE_CURVE_EVAL_WITH_DERIV  Evaluate B-spline curve AND its derivative.
        %
        %   [pts, dpts] = bspline_curve_eval_with_deriv(knots, ctrl_pts, degree, t)
        %
        %   The derivative of a degree-p B-spline with control points P_i
        %   and knot vector T is a degree-(p-1) B-spline with:
        %     Q_i = p * (P_{i+1} - P_i) / (T_{i+p+1} - T_{i+1})
        %   and knots T(2:end-1)  (remove first and last).
        %
        %   RETURNS
        %     pts  : [N × 3]  curve points
        %     dpts : [N × 3]  tangent vectors dP/dt

            % Evaluate the curve itself
            pts = WEC_MS2_Parser.bspline_curve_eval(knots, ctrl_pts, degree, t);

            if degree < 1
                dpts = zeros(size(pts));
                return;
            end

            % Derivative control points
            n = size(ctrl_pts, 1);
            Q = zeros(n - 1, 3);
            for i = 1:n-1
                denom = knots(i + degree + 1) - knots(i + 1);
                if abs(denom) > 1e-14
                    Q(i, :) = degree * (ctrl_pts(i+1, :) - ctrl_pts(i, :)) / denom;
                end
            end

            % Derivative knot vector (remove first and last)
            d_knots = knots(2:end-1);
            d_degree = degree - 1;

            % Evaluate derivative B-spline
            dpts = WEC_MS2_Parser.bspline_curve_eval(d_knots, Q, d_degree, t);
        end


        function [pts, dpts] = arc_evaluate(p_start, p_centre, p_end, t)
        % ARC_EVALUATE  Circular arc from start to end around centre.
        %
        %   [pts, dpts] = arc_evaluate(start, centre, end, t)
        %
        %   t ∈ [0,1]:  t=0 → start,  t=1 → end.
        %
        %   DERIVATION
        %     v_s = start − centre,  v_e = end − centre,  r = |v_s|.
        %     e_r = v_s / r  (initial radial),
        %     n   = cross(v_s, v_e) / |cross(v_s, v_e)|  (arc-plane normal),
        %     e_t = cross(n, e_r)  (tangential).
        %     θ_total = atan2(dot(v_e, e_t), dot(v_e, e_r)).
        %     P(t) = centre + r (cos(θ) e_r + sin(θ) e_t),  θ = t·θ_total.
        %     dP/dt = r (−sin(θ) e_r + cos(θ) e_t) · θ_total.

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
                % Degenerate: collinear points — treat as line
                pts  = (1 - t) .* p_start + t .* p_end;
                dpts = repmat(p_end - p_start, N, 1);
                return;
            end

            n_hat = cp / cp_len;
            e_t   = cross(n_hat, e_r);

            theta_total = atan2(dot(v_e, e_t), dot(v_e, e_r));

            % Ensure we take the short arc (|θ| ≤ π)
            % For C0.ms2 the arcs are 90° — always short.

            theta = t * theta_total;

            pts  = p_centre + r * (cos(theta) .* e_r + sin(theta) .* e_t);
            dpts = r * theta_total * (-sin(theta) .* e_r + cos(theta) .* e_t);
        end

    end % methods (Static)


    %% ═══════════════════════════════════════════════════════════════
    %%  §8  SURFACE TOPOLOGY — Source/Mirror Classification
    %% ═══════════════════════════════════════════════════════════════
    %
    %  WHY classify surfaces?
    %    The .ms2 file defines a hull as a patchwork of parametric
    %    surfaces.  Many visible surfaces are MirrSurfs — exact copies
    %    of a source surface with one coordinate flipped.  For both
    %    panelization and property computation, we only need to
    %    EVALUATE the source surfaces.  Mirror contributions are
    %    derived by coordinate flip (panels) or analytical symmetry
    %    (volume integrals).
    %
    %    Shared edges between surfaces are identified by tracing
    %    parent references in the entity DAG — not by sampling
    %    geometry and comparing point-by-point.  This is exact and
    %    costs zero eval_surface calls.
    %
    %  C0 EXAMPLE
    %    Sources:   surface2 (RevSurf),  surface4 (DevSurf)
    %    Mirrors:   surface8 = MirrY(surface2)
    %               surface5 = MirrX(surface4)
    %               surface6 = MirrY(surface4)
    %               surface7 = MirrY(surface5) → effective MirrXY(surface4)
    %    Adjacency: surface4.edge1 ↔ surface2.edge1 via snake3
    %               (snake3 = EdgeSnake(edge 1, surface2))

    methods

        function info = classify_visible_surfaces(obj)
        % CLASSIFY_VISIBLE_SURFACES  Identify sources, mirrors, and shared edges.
        %
        %   info = model.classify_visible_surfaces()
        %
        %   OUTPUT  struct with fields:
        %     .sources    — cell array of source surface names (non-MirrSurf)
        %     .mirrors    — struct array with fields:
        %                     .name   — mirror surface name
        %                     .source — source surface name
        %                     .plane  — 'X' or 'Y'
        %     .adjacency  — struct array with fields:
        %                     .surf_a, .edge_a — first surface and its edge index
        %                     .surf_b, .edge_b — second surface and its edge index
        %                     .via_entity      — the EdgeSnake entity connecting them
        %
        %   The classification uses ONLY the entity DAG (parent references).
        %   Zero geometry evaluations are performed.

            sources = {};
            mirrors = struct('name', {}, 'source', {}, 'plane', {});

            % ── Step 1: Classify each visible surface ─────────────
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

            % ── Step 2: Resolve chained mirrors ───────────────────
            %  A mirror may reference another mirror (e.g. surface7 =
            %  MirrY(surface5), surface5 = MirrX(surface4)).  Trace
            %  each mirror back to its ultimate non-mirror source and
            %  accumulate the effective flip.
            for i = 1:length(mirrors)
                [ult_source, ult_flips] = obj.resolve_mirror_chain(mirrors(i).name);
                mirrors(i).ultimate_source = ult_source;
                mirrors(i).effective_flips = ult_flips;
            end

            % ── Step 3: Detect shared edges between sources ───────
            adjacency = struct('surf_a', {}, 'edge_a', {}, ...
                               'surf_b', {}, 'edge_b', {}, ...
                               'via_entity', {});

            for i = 1:length(sources)
                sname = sources{i};
                edges_i = obj.get_boundary_entities(sname);

                % For each named boundary entity, trace to see if it's
                % an EdgeSnake of another source surface
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

            % ── Assemble output ───────────────────────────────────
            info = struct();
            info.sources   = sources;
            info.mirrors   = mirrors;
            info.adjacency = adjacency;
        end


        function [ult_source, flips] = resolve_mirror_chain(obj, name)
        % RESOLVE_MIRROR_CHAIN  Follow MirrSurf references to the ultimate source.
        %
        %   [source, flips] = model.resolve_mirror_chain('surface7')
        %
        %   Returns the non-MirrSurf ancestor and a cell array of flip
        %   planes applied in order.  For surface7 = MirrY(MirrX(surface4)):
        %     source = 'surface4', flips = {'X', 'Y'}

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
        % GET_BOUNDARY_ENTITIES  Identify the entity defining each parametric edge.
        %
        %   edges = model.get_boundary_entities('surface4')
        %
        %   Returns struct with fields edge1..edge4, each containing the
        %   entity name that analytically defines that boundary curve.
        %   Empty string means the edge is derived (not a named entity).
        %
        %   EDGE NUMBERING (same as EdgeSnake convention):
        %     edge1: v=0, u varies    edge3: v=1, u varies
        %     edge2: u=1, v varies    edge4: u=0, v varies
        %
        %   SURFACE TYPE MAPPING:
        %     RuledSurf(C1, C2):     edge1 = C1, edge3 = C2
        %     DevSurf(snake, curve): edge1 = snake, edge3 = curve
        %     RevSurf(profile, ...): edges are derived (arcs and profile
        %                            at specific angles — not named entities)
        %     BLoftSurf({sections}): edge1 ≈ first section, edge3 ≈ last
        %                            (approximate for degree > 1)

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
                    % RevSurf edges are arcs and profile slices — not
                    % independently named entities in the DAG.  The profile
                    % curve (e.params.profile) parameterizes the u-direction.
                    % v=0 and v=1 edges are the profile at angle_start and
                    % angle_end respectively — same 3D curve, different
                    % revolution angles.  We leave these empty.
                    % (Adjacency through RevSurf is detected when ANOTHER
                    % surface references an EdgeSnake of this RevSurf.)

                case 'MirrSurf'
                    % Boundary entities of a MirrSurf are the mirrors of
                    % the source surface's boundary entities.  Resolve by
                    % delegating to the source.
                    src_edges = obj.get_boundary_entities(e.params.source);
                    edges = src_edges;  % same entity names — the flip is applied during eval
            end
        end


        function [found, parent_surf, parent_edge] = trace_to_edge_snake_of(obj, entity_name, source_list)
        % TRACE_TO_EDGE_SNAKE_OF  Follow an entity to an EdgeSnake of a source surface.
        %
        %   [found, surf, edge] = model.trace_to_edge_snake_of('snake3', {'surface2','surface4'})
        %
        %   Traces through BSubSnake → parent → ... → EdgeSnake.
        %   If the EdgeSnake references a surface in source_list,
        %   returns found=true with the surface name and edge index.

            found = false;
            parent_surf = '';
            parent_edge = 0;

            if ~obj.entities.isKey(entity_name), return; end
            e = obj.entities(entity_name);

            % ── Direct EdgeSnake ──────────────────────────────────
            if strcmp(e.type, 'EdgeSnake')
                parent_surf = e.params.surface_name;
                parent_edge = e.params.edge_index;
                found = any(strcmp(parent_surf, source_list));
                return;
            end

            % ── BSubSnake or BSubCurve: trace through bead's parent ─
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

            % ── AbsBead: trace through parent curve ───────────────
            if strcmp(e.type, 'AbsBead')
                if isfield(e.params, 'parent_curve')
                    [found, parent_surf, parent_edge] = ...
                        obj.trace_to_edge_snake_of(e.params.parent_curve, source_list);
                end
                return;
            end

            % ── Other entity types: not traceable to a surface edge
        end


        function names = get_required_entities(obj, surface_name)
        % GET_REQUIRED_ENTITIES  Trace all entities needed by a surface.
        %
        %   names = model.get_required_entities('surface3')
        %
        %   Returns a cell array of all entity names in the dependency
        %   chain of the given surface, in topological order (parents
        %   before children).

            visited = containers.Map();
            order   = {};

            function dfs(name)
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
        % GEOMETRY_SUMMARY  Print and return a summary of the parsed model.

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


        function clear_cache(obj)
        % CLEAR_CACHE  Invalidate all cached geometry evaluations.
        %
        %   model.clear_cache()
        %
        %   MUST be called after modifying control point positions
        %   (e.g. geometry scaling via entity.params.offset changes).
        %   Without clearing, eval_point returns stale cached values.
        %
        %   IMPLEMENTATION NOTE:
        %     WEC_MS2_Parser is a value class, but containers.Map is a
        %     handle class.  We must clear the EXISTING Maps (which the
        %     caller still references) rather than replacing them with
        %     new empty Maps (which would only affect the local copy).
        %
        %   USAGE
        %     e = model.entities('pt2');
        %     e.params.offset(1) = e.params.offset(1) * 1.5;
        %     model.entities('pt2') = e;
        %     model.clear_cache();   % ← required

            k = keys(obj.point_cache);
            if ~isempty(k)
                remove(obj.point_cache, k);
            end

            k = keys(obj.rev_axis_cache);
            if ~isempty(k)
                remove(obj.rev_axis_cache, k);
            end
        end

    end % methods

end % classdef