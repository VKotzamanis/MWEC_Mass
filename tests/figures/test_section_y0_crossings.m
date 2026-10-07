function test_section_y0_crossings()
%TEST_SECTION_Y0_CROSSINGS  y = 0 crossings of closed section loops, on the exact curves.
%   The loops are built here from known closed forms (circle of rational quadratic arcs, polygons of
%   degree-1 pieces) as an independent oracle. No kernel function is used: section_y0_crossings
%   works on the control points of the curves.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);

% Circle of radius R about (cx, cy), four rational quadratic arcs starting at angle a0: the crossings of
% y = 0 are inside the pieces (a0 = 30 degrees), so they are found by the root solve. The roots are
% cx -+ sqrt(R^2 - cy^2). Evaluation of a NURBS point and the Illinois solve each carry a few ulp
% of the coordinates (R, cx, cy are of order 1), so 16 ulp of the largest coordinate bounds the error.
R = 1.7;
cx = 0.3;
cy = 0.2;
for a0 = [pi / 6, 0, pi / 2, 7 * pi / 4 + 0.1]
    loop = circle_loop(R, cx, cy, a0);
    x = mwecmass.output.figures.section_y0_crossings(loop);
    expected = [cx - sqrt(R^2 - cy^2); cx + sqrt(R^2 - cy^2)];
    err = max(abs(x - expected));
    fprintf('circle R %.2f centre (%.1f, %.1f), first arc at %.3f rad: crossings %.15f %.15f, error %.2e\n', ...
        R, cx, cy, a0, x, err);
    check(numel(x) == 2 && err <= 16 * eps(2), 'circle crossings off by %.3e', err);
end

% A circle centred on y = 0, first arc at angle 0: the vertex at angle 0 has y = 0 bitwise (found once,
% half-open rule) and the one at angle pi carries sin(pi) rounding (found inside a piece); the
% crossings are x = cx -+ R.
loop = circle_loop(R, cx, 0, 0);
x = mwecmass.output.figures.section_y0_crossings(loop);
check(numel(x) == 2 && abs(x(1) - (cx - R)) <= 16 * eps(2) && abs(x(2) - (cx + R)) <= 16 * eps(2), ...
    'circle on the axis: %s', mat2str(x));
fprintf('circle centred on y = 0: crossings %.15f %.15f (exact %.15f, %.15f)\n', x, cx - R, cx + R);

% Diamond with vertices (-1,0), (0,-1), (1,0), (0,1): vertices exactly on y = 0 are crossings, once each.
loop = polygon_loop([1 0; 0 1; -1 0; 0 -1]);
x = mwecmass.output.figures.section_y0_crossings(loop);
check(isequal(x, [-1; 1]), 'diamond: %s', mat2str(x));

% A vertex that only touches y = 0 from above is a double root: two crossings at the same x.
loop = polygon_loop([1 1; 0 0; -1 1; 0 2]);
x = mwecmass.output.figures.section_y0_crossings(loop);
check(isequal(x, [0; 0]), 'touch: %s', mat2str(x));

% A loop on one side of y = 0 has no crossing.
loop = polygon_loop([1 1; -1 1; -1 2; 1 2]);
x = mwecmass.output.figures.section_y0_crossings(loop);
check(isempty(x), 'loop above y = 0 has %d crossings', numel(x));
loop = polygon_loop([1 -1; 1 -2; -1 -2; -1 -1]);
check(isempty(mwecmass.output.figures.section_y0_crossings(loop)), 'loop below y = 0');

% Edge lying along y = 0 (a degenerate slit-like loop through the axis): [y > 0] counts it consistently.
loop = polygon_loop([2 0; 0 0; 0 1; 2 1]);
x = mwecmass.output.figures.section_y0_crossings(loop);
check(numel(x) == 2, 'edge on the axis: %d crossings', numel(x));
fprintf('loop with an edge on y = 0: crossings %s\n', mat2str(x'));

% Non-convex loop crossing the axis four times (a U shape): sorted, even.
loop = polygon_loop([0 -1; 3 -1; 3 1; 2 1; 2 -0.5; 1 -0.5; 1 1; 0 1]);
x = mwecmass.output.figures.section_y0_crossings(loop);
check(isequal(x, [0; 1; 2; 3]), 'U shape: %s', mat2str(x'));
% Two crossings inside a short parameter interval of one piece: x(s) = 2s - 1, y(s) = (s - 0.515)^2 - 1e-4,
% a quadratic Bezier piece closed by three lines. The exact roots are s = 0.505 and 0.525, so x = 0.01 and
% 0.05 (closed form, independent oracle). The control values of y are rounded to a few ulp of 0.27 and
% |dy/ds| = 0.02 at the roots, dx = 2 ds, so 16 * eps(0.27) / 0.02 bounds the error.
B = -1.03;
C = 0.515^2 - 1e-4;
ctrl = [-1, C, 0; 0, C + B / 2, 0; 1, 1 + B + C, 0];
loop = struct('z', 0, 'pieces', struct('patch', 1, 'u', 0, 'dir', 1, ...
    'curve', struct('degree', 2, 'ctrl', ctrl, 'knots', [0 0 0 1 1 1], 'weights', [])));
corners = [ctrl(3, :); 1 2 0; -1 2 0; ctrl(1, :)];
for q = 1:3
    loop.pieces(end + 1) = struct('patch', q + 1, 'u', 0, 'dir', 1, 'curve', struct('degree', 1, ...
        'ctrl', corners(q:q + 1, :), 'knots', [0 0 1 1], 'weights', []));
end
x = mwecmass.output.figures.section_y0_crossings(loop);
err = max(abs(x - [0.01; 0.05]));
fprintf('two crossings in one sample interval: %.15f %.15f, error %.2e\n', x, err);
check(numel(x) == 2 && err <= 16 * eps(0.27) / 0.02, 'double crossing: %s', mat2str(x'));

% A cubic piece with three simple roots at s = 0.2, 0.5, 0.8 (x = s), closed by one line from its end: the
% three roots of the cubic plus one crossing of the closing line make an even count; the closing line
% runs from (1, y(1)) = (1, 0.08) to (0, y(0)) = (0, -0.08) and crosses y = 0 at x = 0.5.
f = @(s) (s - 0.2) .* (s - 0.5) .* (s - 0.8);
s4 = [0 1/3 2/3 1];
V = f(s4);
% Bezier ordinates of the cubic through 4 points at s4 (interpolation of equispaced values)
M = zeros(4);
for i = 0:3
    for j = 0:3
        M(i + 1, j + 1) = nchoosek(3, j) * s4(i + 1)^j * (1 - s4(i + 1))^(3 - j);
    end
end
yb = M \ V';
ctrl = [(0:3)' / 3, yb, zeros(4, 1)];
loop = struct('z', 0, 'pieces', struct('patch', 1, 'u', 0, 'dir', 1, ...
    'curve', struct('degree', 3, 'ctrl', ctrl, 'knots', [0 0 0 0 1 1 1 1], 'weights', [])));
loop.pieces(2) = struct('patch', 2, 'u', 0, 'dir', 1, 'curve', struct('degree', 1, ...
    'ctrl', [ctrl(4, :); ctrl(1, :)], 'knots', [0 0 1 1], 'weights', []));
x = mwecmass.output.figures.section_y0_crossings(loop);
err = max(abs(x - [0.2; 0.5; 0.5; 0.8]));
fprintf('cubic with three roots plus closing line: %s, error %.2e\n', mat2str(x', 15), err);
check(numel(x) == 4 && err <= 64 * eps, 'cubic crossings: %s', mat2str(x'));

% An open chain (one crossing) is not a closed loop and is reported, not returned.
loop = polygon_loop([1 1; -1 1; -1 -1]);
loop.pieces(end) = [];
msg = '';
try
    mwecmass.output.figures.section_y0_crossings(loop);
catch err
    msg = err.identifier;
end
check(strcmp(msg, 'mwecmass:figures:SectionTopology'), 'open chain gave "%s"', msg);
fprintf('all crossing tests passed\n');
end

function loop = circle_loop(R, cx, cy, a0)
pieces = struct('patch', {}, 'u', {}, 'curve', {}, 'dir', {});
for q = 0:3
    a = a0 + q * pi / 2;
    p0 = [cx + R * cos(a), cy + R * sin(a), 0];
    p2 = [cx + R * cos(a + pi / 2), cy + R * sin(a + pi / 2), 0];
    p1 = [cx + R * sqrt(2) * cos(a + pi / 4), cy + R * sqrt(2) * sin(a + pi / 4), 0];
    curve = struct('degree', 2, 'ctrl', [p0; p1; p2], 'knots', [0 0 0 1 1 1], 'weights', [1; sqrt(2) / 2; 1]);
    pieces(end + 1) = struct('patch', q + 1, 'u', 0, 'curve', curve, 'dir', 1); %#ok<AGROW>
end
loop = struct('z', 0, 'pieces', pieces);
end

function loop = polygon_loop(P)
n = size(P, 1);
pieces = struct('patch', {}, 'u', {}, 'curve', {}, 'dir', {});
for k = 1:n
    a = P(k, :);
    b = P(mod(k, n) + 1, :);
    curve = struct('degree', 1, 'ctrl', [a 0; b 0], 'knots', [0 0 1 1], 'weights', []);
    pieces(end + 1) = struct('patch', k, 'u', 0, 'curve', curve, 'dir', 1); %#ok<AGROW>
end
loop = struct('z', 0, 'pieces', pieces);
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'standins'), '-end');
addpath(fullfile(root, 'tests', 'standins', 'fixtures'), '-end');
end

function check(cond, varargin)
if ~cond
    error('test_section_y0_crossings:fail', varargin{:});
end
end
