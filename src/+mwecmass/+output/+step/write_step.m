function write_step(brep, filename)
%WRITE_STEP  Write a B-rep to an ISO 10303-21 file with AP214 (AUTOMOTIVE_DESIGN) entities.
%
%   mwecmass.output.step.write_step(brep, filename)
%
%   All lengths are metres and the file declares SI_UNIT($,.METRE.) with no prefix. Reals are
%   written with %.17g. Each body becomes one product with its own shape representation: solids
%   as ADVANCED_BREP_SHAPE_REPRESENTATION (MANIFOLD_SOLID_BREP, or BREP_WITH_VOIDS when the body
%   has void shells), sheets as MANIFOLD_SURFACE_SHAPE_REPRESENTATION (SHELL_BASED_SURFACE_MODEL
%   of OPEN_SHELLs). Only entities reachable from a body are written. Each body gets the layers
%   <name>_faces, <name>_edges and <name>_vertices; an edge or vertex shared with an earlier body
%   stays on that body's layer.
%
%   Input struct brep (all indices are 1-based):
%     uncertainty  optional scalar, stated distance accuracy in metres (default 1e-7)
%     vertices     [nv x 3] coordinates
%     curves       struct array with fields
%                    degree   spline degree
%                    ctrl     [n x 3] control points
%                    knots    clamped knot vector of length n+degree+1
%                    weights  [] (non-rational) or [n x 1] positive weights (rational)
%     edges        struct array with fields
%                    vertices [start end] vertex indices
%                    curve    curve index; the curve runs from the start vertex to the end vertex
%     surfaces     cell array of structs, each with field type
%                    'plane'   origin [1x3], normal [1x3], xdir [1x3] (optional, any direction
%                              not parallel to normal)
%                    'bspline' degree [du dv], ctrl [nu x nv x 3], knots {ku, kv} (clamped),
%                              weights [] or [nu x nv]
%     faces        struct array with fields
%                    surface    surface index
%                    same_sense true if the face normal equals the surface normal
%                    loops      cell array of signed edge index vectors; loops{1} is the outer
%                               bound, counter-clockwise seen from the face normal; a negative
%                               index traverses that edge from its end vertex to its start
%     bodies       struct array with fields
%                    name    printable ASCII, unique
%                    kind    'solid' or 'sheet'
%                    shells  cell array of signed face index vectors; a negative index uses the
%                            face with reversed normal, so neighbouring bodies can share a face.
%                            Solid: shells{1} is the outer closed shell (normals out of the
%                            solid), further shells are voids, each oriented as a standalone
%                            solid (normals out of the cavity). Sheet: each shell is an open
%                            shell.
%
%   Edges are shared by index, so faces that meet along an edge reference one EDGE_CURVE.

mwecmass.output.step.validate_brep(brep);
unc = mwecmass.output.step.declared_uncertainty(brep);
[~, base, ext] = fileparts(filename);

lines = cell(1024, 1);
nl = 0;
vid = zeros(size(brep.vertices, 1), 1);
cid = zeros(numel(brep.curves), 1);
eid = zeros(numel(brep.edges), 1);
sid = zeros(numel(brep.surfaces), 1);
fid = zeros(numel(brep.faces), 2);
used = struct('faces', [], 'edges', [], 'vertices', []);

ac = emit('APPLICATION_CONTEXT(''core data for automotive mechanical design processes'')');
emit(sprintf('APPLICATION_PROTOCOL_DEFINITION(''international standard'',''automotive_design'',2000,#%d)', ac));
mc = emit(sprintf('MECHANICAL_CONTEXT('''',#%d,''mechanical'')', ac));
pdc = emit(sprintf('PRODUCT_DEFINITION_CONTEXT(''part definition'',#%d,''design'')', ac));
len = emit('(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.))');
ang = emit('(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.))');
ste = emit('(NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT())');
um = emit(sprintf('UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(%s),#%d,''distance_accuracy_value'',''stated accuracy of the model'')', reals(unc), len));
ctx = emit(sprintf(['(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#%d)) ' ...
    'GLOBAL_UNIT_ASSIGNED_CONTEXT((#%d,#%d,#%d)) REPRESENTATION_CONTEXT(''MWEC_Mass'',''3D''))'], um, len, ang, ste));

layers = cell(1, 0);
for k = 1:numel(brep.bodies)
    used = struct('faces', [], 'edges', [], 'vertices', []);
    bd = brep.bodies(k);
    [item, reptype] = body_item(bd);
    frame = emit(sprintf('AXIS2_PLACEMENT_3D('''',#%d,#%d,#%d)', point([0 0 0]), direction([0 0 1]), direction([1 0 0])));
    rep = emit(sprintf('%s(%s,(#%d,#%d),#%d)', reptype, quote(bd.name), frame, item, ctx));
    prod = emit(sprintf('PRODUCT(%s,%s,'''',(#%d))', quote(bd.name), quote(bd.name), mc));
    emit(sprintf('PRODUCT_RELATED_PRODUCT_CATEGORY(''part'',$,(#%d))', prod));
    pdf = emit(sprintf('PRODUCT_DEFINITION_FORMATION('''','''',#%d)', prod));
    pd = emit(sprintf('PRODUCT_DEFINITION(''design'','''',#%d,#%d)', pdf, pdc));
    pds = emit(sprintf('PRODUCT_DEFINITION_SHAPE('''','''',#%d)', pd));
    emit(sprintf('SHAPE_DEFINITION_REPRESENTATION(#%d,#%d)', pds, rep));
    layers{end + 1} = layer(bd.name, 'faces', used.faces); %#ok<AGROW>
    layers{end + 1} = layer(bd.name, 'edges', used.edges); %#ok<AGROW>
    layers{end + 1} = layer(bd.name, 'vertices', used.vertices); %#ok<AGROW>
end
for k = 1:numel(layers)
    if ~isempty(layers{k})
        emit(layers{k});
    end
end

fid_out = fopen(filename, 'w');
if fid_out < 0
    error('mwecmass:step:io', 'write_step: cannot open %s for writing', filename);
end
fprintf(fid_out, 'ISO-10303-21;\nHEADER;\n');
fprintf(fid_out, 'FILE_DESCRIPTION((''MWEC_Mass B-rep''),''2;1'');\n');
fprintf(fid_out, 'FILE_NAME(%s,''%s'',(''MWEC_Mass''),(''''),''mwecmass.output.step.write_step'','''','''');\n', ...
    quote([base ext]), datestr(now, 'yyyy-mm-ddTHH:MM:SS'));
fprintf(fid_out, 'FILE_SCHEMA((''AUTOMOTIVE_DESIGN''));\nENDSEC;\nDATA;\n');
fprintf(fid_out, '%s\n', strjoin(lines(1:nl), newline));
fprintf(fid_out, 'ENDSEC;\nEND-ISO-10303-21;\n');
fclose(fid_out);

    function id = emit(body)
        nl = nl + 1;
        if nl > numel(lines)
            lines(2 * numel(lines), 1) = {[]};
        end
        lines{nl} = sprintf('#%d=%s;', nl, body);
        id = nl;
    end

    function id = point(xyz)
        id = emit(sprintf('CARTESIAN_POINT('''',(%s))', reals(xyz)));
    end

    function id = direction(d)
        id = emit(sprintf('DIRECTION('''',(%s))', reals(d / norm(d))));
    end

    function ids = points(P)
        ids = zeros(1, size(P, 1));
        for i = 1:size(P, 1)
            ids(i) = point(P(i, :));
        end
    end

    function id = vertex(i)
        if vid(i) == 0
            vid(i) = emit(sprintf('VERTEX_POINT('''',#%d)', point(brep.vertices(i, :))));
            used.vertices(end + 1) = vid(i);
        end
        id = vid(i);
    end

    function id = curve(i)
        if cid(i) == 0
            c = brep.curves(i);
            [mult, kn] = distinct_knots(c.knots);
            head = sprintf('%d,(%s),.UNSPECIFIED.,.F.,.F.', c.degree, refs(points(c.ctrl)));
            knots = sprintf('(%s),(%s),.UNSPECIFIED.', ints(mult), reals(kn));
            if isempty(c.weights)
                cid(i) = emit(sprintf('B_SPLINE_CURVE_WITH_KNOTS('''',%s,%s)', head, knots));
            else
                cid(i) = emit(sprintf(['(BOUNDED_CURVE() B_SPLINE_CURVE(%s) B_SPLINE_CURVE_WITH_KNOTS(%s) CURVE() ' ...
                    'GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((%s)) REPRESENTATION_ITEM(''''))'], ...
                    head, knots, reals(c.weights)));
            end
        end
        id = cid(i);
    end

    function id = edge(i)
        if eid(i) == 0
            ed = brep.edges(i);
            eid(i) = emit(sprintf('EDGE_CURVE('''',#%d,#%d,#%d,.T.)', ...
                vertex(ed.vertices(1)), vertex(ed.vertices(2)), curve(ed.curve)));
            used.edges(end + 1) = eid(i);
        end
        id = eid(i);
    end

    function id = surface(i)
        if sid(i) == 0
            s = brep.surfaces{i};
            if strcmp(s.type, 'plane')
                n = s.normal(:)' / norm(s.normal);
                if isfield(s, 'xdir') && ~isempty(s.xdir)
                    x = s.xdir(:)';
                else
                    [~, m] = min(abs(n));
                    x = circshift([1 0 0], m - 1);
                end
                x = x - (x * n') * n;
                axis2 = emit(sprintf('AXIS2_PLACEMENT_3D('''',#%d,#%d,#%d)', ...
                    point(s.origin(:)'), direction(n), direction(x)));
                sid(i) = emit(sprintf('PLANE('''',#%d)', axis2));
            else
                sid(i) = bspline_surface(s);
            end
        end
        id = sid(i);
    end

    function id = bspline_surface(s)
        [nu, nv, ~] = size(s.ctrl);
        rows = cell(1, nu);
        for a = 1:nu
            rows{a} = ['(' refs(points(reshape(s.ctrl(a, :, :), nv, 3))) ')'];
        end
        [mu, ku] = distinct_knots(s.knots{1});
        [mv, kv] = distinct_knots(s.knots{2});
        head = sprintf('%d,%d,(%s),.UNSPECIFIED.,.F.,.F.,.F.', s.degree(1), s.degree(2), strjoin(rows, ','));
        knots = sprintf('(%s),(%s),(%s),(%s),.UNSPECIFIED.', ints(mu), ints(mv), reals(ku), reals(kv));
        if isempty(s.weights)
            id = emit(sprintf('B_SPLINE_SURFACE_WITH_KNOTS('''',%s,%s)', head, knots));
        else
            wrows = cell(1, nu);
            for a = 1:nu
                wrows{a} = ['(' reals(s.weights(a, :)) ')'];
            end
            id = emit(sprintf(['(BOUNDED_SURFACE() B_SPLINE_SURFACE(%s) B_SPLINE_SURFACE_WITH_KNOTS(%s) ' ...
                'GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE((%s)) REPRESENTATION_ITEM('''') SURFACE())'], ...
                head, knots, strjoin(wrows, ',')));
        end
    end

    function id = face(signed)
        i = abs(signed);
        col = 1 + (signed < 0);
        if fid(i, col) == 0
            f = brep.faces(i);
            bounds = zeros(1, numel(f.loops));
            for j = 1:numel(f.loops)
                L = f.loops{j}(:)';
                if signed < 0
                    L = -fliplr(L);
                end
                oe = zeros(1, numel(L));
                for q = 1:numel(L)
                    oe(q) = emit(sprintf('ORIENTED_EDGE('''',*,*,#%d,%s)', edge(abs(L(q))), tf(L(q) > 0)));
                end
                lp = emit(sprintf('EDGE_LOOP('''',(%s))', refs(oe)));
                if j == 1
                    bounds(j) = emit(sprintf('FACE_OUTER_BOUND('''',#%d,.T.)', lp));
                else
                    bounds(j) = emit(sprintf('FACE_BOUND('''',#%d,.T.)', lp));
                end
            end
            fid(i, col) = emit(sprintf('ADVANCED_FACE('''',(%s),#%d,%s)', refs(bounds), ...
                surface(f.surface), tf(f.same_sense == (signed > 0))));
            used.faces(end + 1) = fid(i, col);
        end
        id = fid(i, col);
    end

    function id = shell(signed_faces, type)
        ids = zeros(1, numel(signed_faces));
        for q = 1:numel(signed_faces)
            ids(q) = face(signed_faces(q));
        end
        id = emit(sprintf('%s('''',(%s))', type, refs(ids)));
    end

    function [item, reptype] = body_item(bd)
        if strcmp(bd.kind, 'solid')
            reptype = 'ADVANCED_BREP_SHAPE_REPRESENTATION';
            outer = shell(bd.shells{1}, 'CLOSED_SHELL');
            if numel(bd.shells) == 1
                item = emit(sprintf('MANIFOLD_SOLID_BREP(%s,#%d)', quote(bd.name), outer));
            else
                voids = zeros(1, numel(bd.shells) - 1);
                for q = 2:numel(bd.shells)
                    cs = shell(bd.shells{q}, 'CLOSED_SHELL');
                    voids(q - 1) = emit(sprintf('ORIENTED_CLOSED_SHELL('''',*,#%d,.F.)', cs));
                end
                item = emit(sprintf('BREP_WITH_VOIDS(%s,#%d,(%s))', quote(bd.name), outer, refs(voids)));
            end
        else
            reptype = 'MANIFOLD_SURFACE_SHAPE_REPRESENTATION';
            shells = zeros(1, numel(bd.shells));
            for q = 1:numel(bd.shells)
                shells(q) = shell(bd.shells{q}, 'OPEN_SHELL');
            end
            item = emit(sprintf('SHELL_BASED_SURFACE_MODEL(%s,(%s))', quote(bd.name), refs(shells)));
        end
    end

    function text = layer(name, what, ids)
        if isempty(ids)
            text = '';
        else
            text = sprintf('PRESENTATION_LAYER_ASSIGNMENT(%s,%s,(%s))', quote([name '_' what]), ...
                quote(['Layer of all ' what ' of body ' name]), refs(ids));
        end
    end
end

function s = quote(txt)
s = ['''' strrep(strrep(txt, '\', '\\'), '''', '''''') ''''];
end

function s = tf(flag)
if flag
    s = '.T.';
else
    s = '.F.';
end
end

function s = refs(ids)
s = sprintf('#%d,', ids);
s(end) = [];
end

function s = ints(v)
s = sprintf('%d,', v);
s(end) = [];
end

function s = reals(v)
t = strsplit(sprintf('%.17g,', v), ',');
t(end) = [];
for k = 1:numel(t)
    if ~any(t{k} == '.')
        e = find(t{k} == 'e', 1);
        if isempty(e)
            t{k} = [t{k} '.'];
        else
            t{k} = [t{k}(1:e - 1) '.' t{k}(e:end)];
        end
    end
    t{k} = upper(t{k});
end
s = strjoin(t, ',');
end

function [mult, kn] = distinct_knots(t)
t = t(:)';
last = [find(diff(t) ~= 0) numel(t)];
kn = t(last);
mult = diff([0 last]);
end
