function test_realised_section_gap()
%TEST_REALISED_SECTION_GAP  F12: no gap in the elevation where the void changes its interval count.
%   mock_section/+mwecmass/+solid/body_section.m is a one-module body whose void is one interval at
%   y = 0 below z = 0.05 and two from there up. The change falls between two sampled heights; the
%   polygons must still tile the strip between the module-edge margins (area closure) and the two
%   groups must meet at adjacent floating-point heights.

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

z_split = 0.05;
design = struct('mode', 'modular_precast', 'edges', [-1; 1], 'vs', 5, 't', 0.2, 'z_ballast', -2, 'solid_modules', []);
realised = struct('hull_name', 'mock', 'design', design, 'body', struct('design', design, 'mock', struct('z_split', z_split)), ...
    'props', struct('CG_total', [0 0 0], 'CB', [0 0 0]), 'status', 'accepted', 'reason', '', 'check', struct(), ...
    'vs', 5, 'mode', 'modular_precast');
n_z = 11;
data = mwecmass.output.figures.realised_section_data(realised, [], n_z);

margin = 8 * eps(1);
P = data.polygons;
roles = {P.role};
check(isequal(sort(unique(roles)), {'shell', 'void'}), 'roles drawn are %s', strjoin(unique(roles), ','));
check(sum(strcmp(roles, 'void')) == 3 && sum(strcmp(roles, 'shell')) == 5, ...
    '%d void and %d shell polygons, expected 3 and 5', sum(strcmp(roles, 'void')), sum(strcmp(roles, 'shell')));

A = sum(arrayfun(@(p) polyarea(p.xz(:, 1), p.xz(:, 2)), P));
A_strip = 2 * (2 - 2 * margin);
fprintf('polygon areas %.15f, strip %.15f, difference %.2e\n', A, A_strip, A - A_strip);
check(abs(A - A_strip) <= 64 * eps * A_strip, 'polygons leave %.3e m^2 of the strip uncovered', A_strip - A);

% the two groups meet at adjacent floating-point heights around the change
one = P(strcmp(roles, 'void') & arrayfun(@(p) max(p.xz(:, 1)) == 0.5, P));
two = P(strcmp(roles, 'void') & arrayfun(@(p) min(p.xz(:, 1)) == 0.4, P));
check(numel(one) == 1 && numel(two) == 1, 'void groups: %d and %d', numel(one), numel(two));
check(one.z_hi < z_split && two.z_lo >= z_split && two.z_lo == nextfloat(one.z_hi), ...
    'groups end at %.17g and start at %.17g, not adjacent floats around %g', one.z_hi, two.z_lo, z_split);
fprintf('void groups meet at %.17g and %.17g\n', one.z_hi, two.z_lo);
fprintf('all F12 gap tests passed\n');
end

function z = nextfloat(a)
z = a + eps(a);
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
