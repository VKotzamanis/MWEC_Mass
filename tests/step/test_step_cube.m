function test_step_cube()
%TEST_STEP_CUBE  Unit cube from PLANE faces: one volume, six surfaces, closed, volume 1.

addpath(fileparts(mfilename('fullpath')));
b = stp_new();
[b, shell] = stp_box(b, [0 0 0], [1 1 1]);
b.bodies(1) = struct('name', 'cube', 'kind', 'solid', 'shells', {{shell}});
file = [tempname() '.step'];
cleanup = onCleanup(@() delete(file));
mwecmass.output.step.write_step(b, file);
r = stp_check(file);

fprintf('cube: volumes %d, surfaces %d, open edges %d, OCC volume %.17g, mesh volume %.17g\n', ...
    r.n_volumes, r.n_surfaces, r.open_edges, r.occ_volumes, r.mesh_volumes);
assert_true(r.n_volumes == 1 && r.n_surfaces == 6, 'cube must import as 1 volume and 6 surfaces');
assert_true(r.open_edges == 0, 'cube mesh must have no open edges');
% 1e-9: planar faces are exact, only OCC rounding remains
assert_true(abs(r.occ_volumes - 1) < 1e-9, 'OCC cube volume must be 1');
assert_true(abs(r.mesh_volumes - 1) < 1e-9, 'mesh cube volume must be 1');
end

function assert_true(cond, msg)
if ~cond
    error('test_step_cube: %s', msg);
end
end
