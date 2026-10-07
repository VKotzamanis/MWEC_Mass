function test_realised_section_gap()
%TEST_REALISED_SECTION_GAP  F12: no gap in the elevation and no false boundary where the void changes
%   its interval count inside a module, and a material island inside the air is a hole of the outline.
%   mock_section/+mwecmass/+solid/body_section.m is a one-module body on [-1, 1] whose section along
%   y = 0 changes type at given heights (0 solid, 1 one void interval |x| <= 0.5, 2 two intervals
%   0.4 <= |x| <= 0.8 of one U-shaped void loop). Independent oracle: the closed forms of that mock.
%   Cases: one interval to two; two changes (solid -> one -> two) inside a single sampling cell of
%   n_z = 11 heights (cell width 0.2) at three height pairs; two intervals to one; one to two and back
%   (a material island). Checked: the polygons tile the strip between the module-edge margins, the void
%   area is the closed form, consecutive groups meet at adjacent floating-point heights, there is one
%   void region whose outer loop less its holes has the void polygons' area, the island is one hole
%   [-0.4, 0.4] x [z1, z2] to adjacent floats, and the horizontal boundary at every change between one
%   and two intervals is [-0.8, -0.5], [-0.4, 0.4] and [0.5, 0.8]. A one-ulp bound on a height is the
%   step between adjacent floats that bisection reaches; x bounds are 8 eps, the rounding of the y = 0
%   crossing on a straight piece; area bounds are 64 eps relative, the rounding of polyarea over a few
%   hundred vertices. Printed only: a feature inside one sampling cell
%   (types [1 2 1 2] at 0.02, 0.05, 0.08), which n_z does not resolve (F12 docstring).

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'standins'), '-end');
addpath(fullfile(root, 'tests', 'standins', 'fixtures'), '-end');
mock_dir = fullfile(root, 'tests', 'figures', 'mock_section');
addpath(mock_dir);
cleanup = onCleanup(@() remove_mock(mock_dir));

n_z = 11;
fprintf('one interval to two at 0.05\n');
run_case([-Inf 0.05], [1 2], n_z);
for zz = [0.05 0.07; 0.05 0.15; 0.25 0.30]'
    fprintf('two changes, solid -> 1 at %g -> 2 at %g\n', zz(1), zz(2));
    run_case([-Inf zz(1) zz(2)], [0 1 2], n_z);
end
fprintf('two intervals to one at 0.3\n');
run_case([-Inf 0.3], [2 1], n_z);
fprintf('one interval to two at 0.05 and back to one at 0.5 (material island)\n');
run_case([-Inf 0.05 0.5], [1 2 1], n_z);

fprintf('unresolved: types [1 2 1 2] from 0.02, 0.05, 0.08 with n_z = %d\n', n_z);
data = mock_data([-Inf 0.02 0.05 0.08], [1 2 1 2], n_z);
voids = data.polygons(strcmp({data.polygons.role}, 'void'));
A_void = sum(arrayfun(@(p) polyarea(p.xz(:, 1), p.xz(:, 2)), voids));
A_exp = void_area_closed_form([-Inf 0.02 0.05 0.08], [1 2 1 2]);
fprintf('  void area %.15f, closed form %.15f, difference %.3e (printed, not gated)\n', A_void, A_exp, A_void - A_exp);
fprintf('all F12 gap tests passed\n');
end

function data = mock_data(z_from, type, n_z)
design = struct('mode', 'modular_precast', 'edges', [-1; 1], 'vs', 5, 't', 0.2, 'z_ballast', -2, 'solid_modules', []);
realised = struct('hull_name', 'mock', 'design', design, ...
    'body', struct('design', design, 'mock', struct('z_from', z_from, 'type', type)), ...
    'props', struct('CG_total', [0 0 0], 'CB', [0 0 0]), 'status', 'accepted', 'reason', '', 'check', struct(), ...
    'vs', 5, 'mode', 'modular_precast');
data = mwecmass.output.figures.realised_section_data(realised, [], n_z);
end

function A = void_area_closed_form(z_from, type)
margin = 8 * eps(1);
width = [0 1 0.8];
lo = max(z_from, -1 + margin);
hi = [z_from(2:end), 1 - margin];
A = sum(width(type + 1) .* (hi - lo));
end

function run_case(z_from, type, n_z)
margin = 8 * eps(1);
A_strip = 2 * (2 - 2 * margin);
top = 1 - margin;
data = mock_data(z_from, type, n_z);

P = data.polygons;
roles = {P.role};
n_void = sum(strcmp(roles, 'void'));
n_shell = sum(strcmp(roles, 'shell'));
n_void_exp = sum(type == 1) + 2 * sum(type == 2);
n_shell_exp = sum(type == 0) + 2 * sum(type == 1) + 3 * sum(type == 2);
check(isequal(sort(unique(roles)), unique([{'shell', 'void'}])), 'roles drawn are %s', strjoin(unique(roles), ','));
check(n_void == n_void_exp && n_shell == n_shell_exp, '%d void and %d shell polygons, expected %d and %d', ...
    n_void, n_shell, n_void_exp, n_shell_exp);

A = sum(arrayfun(@(p) polyarea(p.xz(:, 1), p.xz(:, 2)), P));
fprintf('  polygon area %.15f, strip %.15f, difference %.2e\n', A, A_strip, A - A_strip);
check(abs(A - A_strip) <= 64 * eps * A_strip, 'polygons leave %.3e m^2 of the strip uncovered', A_strip - A);

voids = P(strcmp(roles, 'void'));
A_void = sum(arrayfun(@(p) polyarea(p.xz(:, 1), p.xz(:, 2)), voids));
A_exp = void_area_closed_form(z_from, type);
fprintf('  void area %.15f, closed form %.15f, difference %.2e\n', A_void, A_exp, A_void - A_exp);
check(abs(A_void - A_exp) <= 64 * eps * A_exp, 'void area %.15f, closed form %.15f', A_void, A_exp);

% consecutive groups meet at adjacent floating-point heights around each change
starts = z_from(2:end);
zs = unique([P.z_lo, P.z_hi]);
zs = zs(zs > -1 + margin & zs < top);
for c = 1:numel(starts)
    below = zs(zs < starts(c));
    above = zs(zs >= starts(c));
    lo = below(end);
    hi = above(1);
    check(hi == lo + eps(lo) || hi == lo + eps(hi), 'change at %g: groups end at %.17g and start at %.17g', starts(c), lo, hi);
    fprintf('  change %g: groups meet at %.17g and %.17g\n', starts(c), lo, hi);
end

vo = data.void_outlines;
check(numel(vo) == 1, '%d void regions for one connected air region', numel(vo));
check(signed_area(vo.xz) > 0, 'outer loop orientation');
islands = find(type(1:end - 2) == 1 & type(2:end - 1) == 2 & type(3:end) == 1);
check(numel(vo.holes) == numel(islands), '%d holes, expected %d', numel(vo.holes), numel(islands));
A_holes = 0;
for h = 1:numel(vo.holes)
    hz = vo.holes{h};
    a_h = signed_area(hz);
    check(a_h < 0, 'hole %d orientation: signed area %.3e', h, a_h);
    z1 = starts(islands(h));
    z2 = starts(islands(h) + 1);
    a_exp = 0.8 * (z2 - z1);
    fprintf('  hole %d: x [%.17g, %.17g], z [%.17g, %.17g], area %.15f, closed form %.15f\n', h, ...
        min(hz(:, 1)), max(hz(:, 1)), min(hz(:, 2)), max(hz(:, 2)), -a_h, a_exp);
    check(all(abs([min(hz(:, 1)), max(hz(:, 1))] - [-0.4, 0.4]) <= 8 * eps), 'hole %d x range', h);
    check(abs(min(hz(:, 2)) - z1) <= eps(z1) && abs(max(hz(:, 2)) - z2) <= eps(z2), 'hole %d z range', h);
    check(abs(-a_h - a_exp) <= 64 * eps * a_exp, 'hole %d area %.15f, closed form %.15f', h, -a_h, a_exp);
    A_holes = A_holes - a_h;
end
A_out = signed_area(vo.xz);
fprintf('  outer loop area %.15f, void polygons + holes %.15f, difference %.2e\n', A_out, A_void + A_holes, ...
    A_out - A_void - A_holes);
check(abs(A_out - A_void - A_holes) <= 64 * eps * A_out, 'outer loop area %.15f, void polygons %.15f + holes %.15f', ...
    A_out, A_void, A_holes);

% horizontal boundary at every change between one and two intervals, over all loops of the region
loops = [{vo.xz}, vo.holes];
for c = find(type(1:end - 1) + type(2:end) == 3)
    z_a = zs(find(zs < starts(c), 1, 'last'));
    z_b = zs(find(zs >= starts(c), 1, 'first'));
    seg = zeros(0, 2);
    for L = 1:numel(loops)
        xz = loops{L};
        nxt = circshift(xz, -1);
        at = xz(:, 2) >= z_a & xz(:, 2) <= z_b & nxt(:, 2) >= z_a & nxt(:, 2) <= z_b & abs(xz(:, 1) - nxt(:, 1)) > 8 * eps;
        seg = [seg; sort([xz(at, 1), nxt(at, 1)], 2)]; %#ok<AGROW>
    end
    seg = sortrows(seg);
    want = [-0.8 -0.5; -0.4 0.4; 0.5 0.8];
    check(isequal(size(seg), size(want)) && all(abs(seg(:) - want(:)) <= 8 * eps), ...
        'horizontal boundary at %g: %s, expected [-0.8 -0.5; -0.4 0.4; 0.5 0.8]', starts(c), mat2str(seg, 17));
    fprintf('  boundary at %g: %s\n', starts(c), mat2str(seg, 6));
end
end

function a = signed_area(xz)
x = xz(:, 1);
y = xz(:, 2);
a = 0.5 * sum(x .* circshift(y, -1) - circshift(x, -1) .* y);
end

function remove_mock(mock_dir)
rmpath(mock_dir);
clear functions
end

function check(cond, varargin)
if ~cond
    error('test_realised_section_gap:fail', varargin{:});
end
end
