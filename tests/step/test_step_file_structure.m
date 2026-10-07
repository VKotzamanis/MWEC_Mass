function test_step_file_structure()
%TEST_STEP_FILE_STRUCTURE  Schema line, shared edges, layers and orphan points in a two-solid file.

addpath(fileparts(mfilename('fullpath')));
b = stp_new();
[b, lower] = stp_box(b, [0 0 0], [1 1 1]);
[b, upper] = stp_box(b, [0 0 1], [1 1 2]);
b.bodies(1) = struct('name', 'lower', 'kind', 'solid', 'shells', {{lower}});
b.bodies(2) = struct('name', 'upper', 'kind', 'solid', 'shells', {{upper}});
file = [tempname() '.step'];
cleanup = onCleanup(@() delete(file));
mwecmass.output.step.write_step(b, file);
text = fileread(file);

if isempty(strfind(text, 'FILE_SCHEMA((''AUTOMOTIVE_DESIGN''));'))
    error('test_step_file_structure: FILE_SCHEMA(AUTOMOTIVE_DESIGN) missing');
end
n_edges = numel(regexp(text, '=EDGE_CURVE\(', 'match'));
n_vertices = numel(regexp(text, '=VERTEX_POINT\(', 'match'));
fprintf('structure: EDGE_CURVE %d (brep edges %d), VERTEX_POINT %d (brep vertices %d)\n', ...
    n_edges, numel(b.edges), n_vertices, size(b.vertices, 1));
if n_edges ~= numel(b.edges) || n_vertices ~= size(b.vertices, 1)
    error('test_step_file_structure: edges and vertices must be written once each');
end

layered = regexp(text, '=PRESENTATION_LAYER_ASSIGNMENT\([^;]*?\(([^;]*)\)\);', 'tokens');
layer_ids = sort(cell2mat(cellfun(@(c) str2double(regexp(c{1}, '(?<=#)\d+', 'match')), layered, 'UniformOutput', false)));
topo = [];
for kind = {'ADVANCED_FACE', 'EDGE_CURVE', 'VERTEX_POINT'}
    ids = regexp(text, ['#(\d+)=' kind{1} '\('], 'tokens');
    topo = [topo cellfun(@(c) str2double(c{1}), ids)]; %#ok<AGROW>
end
fprintf('structure: %d layer assignments, %d faces, edges and vertices\n', numel(layered), numel(topo));
if ~isequal(layer_ids, sort(topo))
    error('test_step_file_structure: every face, edge and vertex must sit on exactly one layer');
end

points = regexp(text, '#(\d+)=CARTESIAN_POINT\(', 'tokens');
point_ids = cellfun(@(c) str2double(c{1}), points);
referenced = cellfun(@(id) ~isempty(regexp(text, sprintf('#%d[,)]', id), 'once')), num2cell(point_ids));
fprintf('structure: %d CARTESIAN_POINT, %d unreferenced\n', numel(point_ids), nnz(~referenced));
if any(~referenced)
    error('test_step_file_structure: free CARTESIAN_POINT found');
end
end
