function test_step_solid_and_sheet()
%TEST_STEP_SOLID_AND_SHEET  One file with a solid and a sheet sharing junction edges (steel combined file).

addpath(fileparts(mfilename('fullpath')));
b = stp_new();
[b, lateral] = stp_bulged_box(b);
[b, ballast] = stp_box(b, [0 0 -1], [1 1 0]);
b.bodies(1) = struct('name', 'ballast', 'kind', 'solid', 'shells', {{ballast}});
b.bodies(2) = struct('name', 'shell', 'kind', 'sheet', 'shells', {{lateral}});

% The planar faces are exact, so the ballast volume is 1; the junction is the z = 0 square.
junction = intersect(edge_set(b, ballast), edge_set(b, lateral));
fprintf('junction edges shared by the solid and the sheet: %d\n', numel(junction));
if numel(junction) ~= 4
    error('test_step_solid_and_sheet: the solid and the sheet must share the four junction edges');
end

file = [tempname() '.step'];
cleanup = onCleanup(@() delete(file));
mwecmass.output.step.write_step(b, file);
r = stp_check(file);

fprintf('solid+sheet: volumes %d, surfaces %d, open edges volumes %d, sheets %d, OCC volume %.17g, mesh volume %.17g\n', ...
    r.n_volumes, r.n_surfaces, r.open_edges_volumes, r.open_edges_sheets, r.occ_volumes(1), r.mesh_volumes(1));
if r.n_volumes ~= 1 || r.n_surfaces ~= 10
    error('test_step_solid_and_sheet: expected 1 volume and 10 surfaces');
end
if r.open_edges_volumes ~= 0
    error('test_step_solid_and_sheet: the solid must have no open mesh edge');
end
if r.open_edges_sheets == 0
    error('test_step_solid_and_sheet: the bounded sheet must report boundary edges');
end
if abs(r.occ_volumes(1) - 1) > 1e-9
    error('test_step_solid_and_sheet: ballast volume %.17g differs from 1', r.occ_volumes(1));
end

text = fileread(file);
n_edge_curves = numel(regexp(text, '=EDGE_CURVE\(', 'match'));
fprintf('EDGE_CURVE entities in the file: %d, edges in the B-rep reachable from the bodies: %d\n', ...
    n_edge_curves, numel(union(edge_set(b, ballast), edge_set(b, lateral))));
if n_edge_curves ~= numel(union(edge_set(b, ballast), edge_set(b, lateral)))
    error('test_step_solid_and_sheet: shared edges must be written once');
end
end

function e = edge_set(b, shell)
e = [];
for f = abs(shell(:)')
    for j = 1:numel(b.faces(f).loops)
        e = union(e, abs(b.faces(f).loops{j}(:)'));
    end
end
end
