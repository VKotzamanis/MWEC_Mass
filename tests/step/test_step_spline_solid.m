function test_step_spline_solid()
%TEST_STEP_SPLINE_SOLID  Closed solid with bicubic B-spline side faces and planar caps.

addpath(fileparts(mfilename('fullpath')));
b = stp_new();
[b, lateral, caps, exact] = stp_bulged_box(b);
b.bodies(1) = struct('name', 'bulged_box', 'kind', 'solid', 'shells', {{[lateral caps]}});
file = [tempname() '.step'];
cleanup = onCleanup(@() delete(file));
mwecmass.output.step.write_step(b, file);
r = stp_check(file);

fprintf('spline solid: volumes %d, surfaces %d, open edges %d\n', r.n_volumes, r.n_surfaces, r.open_edges);
fprintf('spline solid: exact volume %.17g, OCC %.17g (diff %.3g), mesh %.17g (diff %.3g, mesh size %.3g)\n', ...
    exact, r.occ_volumes, r.occ_volumes - exact, r.mesh_volumes, r.mesh_volumes - exact, r.mesh_size);
if r.n_volumes ~= 1 || r.n_surfaces ~= 6
    error('test_step_spline_solid: expected 1 volume and 6 surfaces');
end
if r.open_edges ~= 0
    error('test_step_spline_solid: the surface mesh has open edges');
end
end
