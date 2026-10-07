function test_step_void()
%TEST_STEP_VOID  BREP_WITH_VOIDS: unit cube with an inner cube void, volume 1 - a^3.

addpath(fileparts(mfilename('fullpath')));
a = 0.4;
lo = (1 - a) / 2;
b = stp_new();
[b, outer] = stp_box(b, [0 0 0], [1 1 1]);
[b, cavity] = stp_box(b, [lo lo lo], [lo lo lo] + a);
b.bodies(1) = struct('name', 'hollow_cube', 'kind', 'solid', 'shells', {{outer, cavity}});
file = [tempname() '.step'];
cleanup = onCleanup(@() delete(file));
mwecmass.output.step.write_step(b, file);
r = stp_check(file);

exact = 1 - a^3;
fprintf('void: volumes %d, surfaces %d, open edges %d, exact %.17g, OCC %.17g, mesh %.17g\n', ...
    r.n_volumes, r.n_surfaces, r.open_edges, exact, r.occ_volumes, r.mesh_volumes);
if r.n_volumes ~= 1 || r.n_surfaces ~= 12
    error('test_step_void: expected 1 volume and 12 surfaces');
end
if r.open_edges ~= 0
    error('test_step_void: the surface mesh has open edges');
end
% 1e-9: all faces planar, so the volume is exact up to rounding
if abs(r.occ_volumes - exact) > 1e-9 || abs(r.mesh_volumes - exact) > 1e-9
    error('test_step_void: volume differs from 1 - a^3');
end
end
