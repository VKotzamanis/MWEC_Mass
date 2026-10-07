function test_step_shared_face()
%TEST_STEP_SHARED_FACE  Two solids sharing a face in one file.

addpath(fileparts(mfilename('fullpath')));
b = stp_new();
[b, lower] = stp_box(b, [0 0 0], [1 1 1]);
[b, upper] = stp_box(b, [0 0 1], [1 1 2]);
b.bodies(1) = struct('name', 'lower', 'kind', 'solid', 'shells', {{lower}});
b.bodies(2) = struct('name', 'upper', 'kind', 'solid', 'shells', {{upper}});
shared = intersect(abs(lower), abs(upper));
fprintf('shared faces: %d, faces %d, edges %d\n', numel(shared), numel(b.faces), numel(b.edges));
if numel(shared) ~= 1 || numel(b.faces) ~= 11
    error('test_step_shared_face: the two boxes must share exactly one face');
end
file = [tempname() '.step'];
cleanup = onCleanup(@() delete(file));
mwecmass.output.step.write_step(b, file);
r = stp_check(file);

fprintf('shared face: volumes %d, surfaces %d, open edges %d, OCC volumes %s, mesh volumes %s\n', ...
    r.n_volumes, r.n_surfaces, r.open_edges, mat2str(r.occ_volumes(:)', 17), mat2str(r.mesh_volumes(:)', 17));
if r.n_volumes ~= 2
    error('test_step_shared_face: expected 2 volumes');
end
end
