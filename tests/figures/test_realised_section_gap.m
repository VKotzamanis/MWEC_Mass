function test_realised_section_gap()
%TEST_REALISED_SECTION_GAP  F12: no gap in the elevation and no false boundary where the void changes
%   its interval count inside a module.
%   mock_section/+mwecmass/+solid/body_section.m is a one-module body on [-1, 1] whose section along
%   y = 0 changes type at given heights (solid, one void interval |x| <= 0.5, two intervals
%   0.4 <= |x| <= 0.8). Cases: one change, and two changes (solid -> one -> two intervals) inside a
%   single sampling cell of n_z = 11 heights (cell width 0.2), at three height pairs. Checked against the
%   closed forms: the polygons tile the strip between the module-edge margins, the void area is the
%   closed form, consecutive groups meet at adjacent floating-point heights, there is one void outline
%   (one connected air region), and the horizontal boundary of that outline at the split of the void
%   is exactly the part of each interval that the other side does not cover: [-0.8, -0.5], [-0.4, 0.4]
%   and [0.5, 0.8].

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

margin = 8 * eps(1);
A_strip = 2 * (2 - 2 * margin);
n_z = 11;

fprintf('one change, type 1 -> 2 at 0.05\n');
run_case(-Inf, 0.05, [1 2], margin, A_strip, n_z);
for zz = [0.05 0.07; 0.05 0.15; 0.25 0.30]'
    fprintf('two changes, solid -> 1 at %g -> 2 at %g\n', zz(1), zz(2));
    run_case([-Inf zz(1)], zz(2), [0 1 2], margin, A_strip, n_z);
end
fprintf('all F12 gap tests passed\n');
end

function run_case(z_from_first, z_last, type, margin, A_strip, n_z)
% z_from_first: the starts of the types before the last; z_last: the start of the last type.
z_from = [z_from_first(:)', z_last];
if numel(type) == 2
    z_from = [-Inf, z_last];
end
design = struct('mode', 'modular_precast', 'edges', [-1; 1], 'vs', 5, 't', 0.2, 'z_ballast', -2, 'solid_modules', []);
realised = struct('hull_name', 'mock', 'design', design, ...
    'body', struct('design', design, 'mock', struct('z_from', z_from, 'type', type)), ...
    'props', struct('CG_total', [0 0 0], 'CB', [0 0 0]), 'status', 'accepted', 'reason', '', 'check', struct(), ...
    'vs', 5, 'mode', 'modular_precast');
data = mwecmass.output.figures.realised_section_data(realised, [], n_z);

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
top = 1 - margin;
starts = z_from(2:end);
if numel(type) == 3
    A_exp = (starts(2) - starts(1)) * 1 + 0.8 * (top - starts(2));
else
    A_exp = 1 * (starts(1) - (-1 + margin)) + 0.8 * (top - starts(1));
end
fprintf('  void area %.15f, closed form %.15f, difference %.2e\n', A_void, A_exp, A_void - A_exp);
check(abs(A_void - A_exp) <= 64 * eps * A_exp, 'void area %.15f, closed form %.15f', A_void, A_exp);

% consecutive groups meet at adjacent floating-point heights around each change
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
check(numel(vo) == 1, '%d void outlines for one connected air region', numel(vo));
A_out = polyarea(vo.xz(:, 1), vo.xz(:, 2));
check(abs(A_out - A_void) <= 64 * eps * A_void, 'outline area %.15f, void polygons %.15f', A_out, A_void);
check(signed_area(vo.xz) > 0, 'outline orientation');

% horizontal boundary at the split of the void (the last change)
z_a = zs(find(zs < starts(end), 1, 'last'));
z_b = zs(find(zs >= starts(end), 1, 'first'));
xz = vo.xz;
nxt = circshift(xz, -1);
at_split = xz(:, 2) >= z_a & xz(:, 2) <= z_b & nxt(:, 2) >= z_a & nxt(:, 2) <= z_b & abs(xz(:, 1) - nxt(:, 1)) > 8 * eps;
seg = sort([xz(at_split, 1), nxt(at_split, 1)], 2);
seg = sortrows(seg);
want = [-0.8 -0.5; -0.4 0.4; 0.5 0.8];
check(isequal(size(seg), size(want)) && all(abs(seg(:) - want(:)) <= 8 * eps), ...
    'horizontal boundary at the split: %s, expected [-0.8 -0.5; -0.4 0.4; 0.5 0.8]', mat2str(seg, 17));
fprintf('  boundary at the split: %s\n', mat2str(seg, 6));
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
