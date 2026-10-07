function test_step_rational_cylinder()
%TEST_STEP_RATIONAL_CYLINDER  Quarter circles as rational quadratic B-splines form a cylinder.
%
%   The cylinder is the independent oracle here (pi r^2 h); the writer only sees B-spline data.

addpath(fileparts(mfilename('fullpath')));
r = 0.75;
h = 1.5;
w = [1 sqrt(2) / 2 1];
q = [r 0; 0 r; -r 0; 0 -r];

b = stp_new();
for k = 1:4
    a = q(k, :);
    c = q(mod(k, 4) + 1, :);
    mid = a + c;
    ctrl = zeros(3, 2, 3);
    for j = 1:2
        z = (j - 1) * h;
        ctrl(:, j, :) = [a z; mid z; c z];
    end
    [b, lateral(k)] = stp_patch(b, [2 1], ctrl, {[0 0 0 1 1 1], [0 0 1 1]}, [w' w']); %#ok<AGROW>
end
v = zeros(1, 4);
for k = 1:4
    [b, v(k)] = stp_vertex(b, [q(k, :) 0]);
end
[b, bottom] = stp_cap(b, v([1 4 3 2]), [0 0 -1]);
for k = 1:4
    [b, v(k)] = stp_vertex(b, [q(k, :) h]);
end
[b, top] = stp_cap(b, v, [0 0 1]);
b.bodies(1) = struct('name', 'cylinder', 'kind', 'solid', 'shells', {{[lateral bottom top]}});
file = [tempname() '.step'];
cleanup = onCleanup(@() delete(file));
mwecmass.output.step.write_step(b, file);
res = stp_check(file);

exact = pi * r^2 * h;
fprintf('cylinder: volumes %d, surfaces %d, open edges %d\n', res.n_volumes, res.n_surfaces, res.open_edges);
fprintf('cylinder: pi r^2 h %.17g, OCC %.17g (diff %.3g), mesh %.17g (diff %.3g, mesh size %.3g)\n', ...
    exact, res.occ_volumes, res.occ_volumes - exact, res.mesh_volumes, res.mesh_volumes - exact, res.mesh_size);
if res.n_volumes ~= 1 || res.n_surfaces ~= 6
    error('test_step_rational_cylinder: expected 1 volume and 6 surfaces');
end
if res.open_edges ~= 0
    error('test_step_rational_cylinder: the surface mesh has open edges');
end
end
