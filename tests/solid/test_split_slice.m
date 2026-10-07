function test_split_slice()
%TEST_SPLIT_SLICE  F3b split_bspline_surface and F4 slice_bspline_surface on fixture patches.
%   The cylinder and box patches of sti_closed_form (exact NURBS of the SK fixture decks) and their
%   closed-form sections serve as the independent oracle; a revolved cubic profile (rational,
%   degree [3 2]) checks the general case.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);

% F4 on the fixture patches against the closed-form sections
worst = struct('cylinder', 0, 'box', 0);
for name = {'cylinder', 'box'}
    fx = sti_closed_form('fixture', name{1});
    P = sti_closed_form('patches', fx, 0);
    for z = [fx.z(1) + 0.25, mean(fx.z), fx.z(2) - 0.5]
        L = mwecmass.solid.slice_bspline_surface(P, z);
        sc = sti_closed_form('section', fx, 0);
        ref = [sc.A, sc.Sx / sc.A, sc.Sy / sc.A, sc.Ixx, sc.Iyy, sc.Ixy];
        got = [L.area, L.centroid, L.I];
        worst.(name{1}) = max(worst.(name{1}), max(abs(got - ref) ./ max(abs(ref), 1)));
        check(L.simple && numel(L.pieces) == 4 && all(L.pts(:, 3) == z), '%s: loop at z = %g', name{1}, z);
        check(all(L.pts(:, 1) .* L.pts([2:end 1], 2) - L.pts([2:end 1], 1) .* L.pts(:, 2) >= -eps), ...
            '%s: loop not counter-clockwise', name{1});
    end
end
fprintf('F4 vs closed-form sections, largest relative difference: box %.3e, cylinder %.3e\n', worst.box, worst.cylinder);
% box rows are lines: Gauss-Legendre integrates them exactly, so only rounding of a few dozen
% operations on values <= 2.25 remains; the cylinder's rational arcs are a quadrature result (printed)
check(worst.box <= 64 * eps, 'F4 box sections differ from the closed form beyond rounding');

% F3b on the cylinder quarter (u-degree 1): parent and pieces agree at the same (u, v)
fx = sti_closed_form('fixture', 'cylinder');
P = sti_closed_form('patches', fx, 0);
[lo, hi, us] = mwecmass.solid.split_bspline_surface(P(1), -1.25);
check(isequal(lo.surf.ctrl(end, :, :), hi.surf.ctrl(1, :, :)) && all(all(hi.surf.ctrl(1, :, 3) == -1.25)), 'cut row shared');
dev = piece_deviation(P(1), lo, hi);
fprintf('cylinder quarter cut at z = -1.25: u* = %.17g, largest |parent - piece| %.3e m\n', us, dev);
check(dev <= 16 * eps * 4, 'cylinder pieces deviate from the parent');

% general case: a quarter revolution (rational arc in v) of a cubic profile with three spans
prof = [0 0 -2; 0.8 0 -2; 1.4 0 -1.6; 1.5 0 -0.7; 1.2 0 0.2; 0.6 0 0.9; 0 0 1];
ku = [0 0 0 0 0.3 0.55 0.8 1 1 1 1];
w = sqrt(2) / 2;
ctrl = zeros(7, 3, 3);
for i = 1:7
    r = prof(i, 1);
    ctrl(i, :, :) = reshape([r 0 prof(i, 3); r r prof(i, 3); 0 r prof(i, 3)], 1, 3, 3);
end
q = P(1);
q.name = 'cubic_rev';
q.surf = struct('type', 'bspline', 'degree', [3 2], 'ctrl', ctrl, 'knots', {{ku, [0 0 0 1 1 1]}}, ...
    'weights', repmat([1 w 1], 7, 1));
q.z_range = [-2 1];
q.c0_u = [];
q.pole = [true true];
Q4 = repmat(q, 1, 4);
for k = 2:4
    for c = find([k == 2 || k == 4, k == 3 || k == 4])
        Q4(k).surf.ctrl(:, :, c) = -Q4(k).surf.ctrl(:, :, c);
    end
end
for z = [-1.9, -0.7, 0.35, 0.95]
    [lo, hi, us] = mwecmass.solid.split_bspline_surface(q, z);
    dev = piece_deviation(q, lo, hi);
    zz = mwecmass.solid.eval_bspline_surface(q.surf, us, 0.3);
    fprintf('cubic revolution cut at z = %5.2f: u* = %.17g, z(u*) - z = %.3e, largest |parent - piece| %.3e m\n', ...
        z, us, zz(3) - z, dev);
    % knot insertion moves no point: rounding of up to 3 insertions on coordinates <= 2
    check(dev <= 64 * eps * 2, 'cubic pieces deviate from the parent');
    check(isequal(lo.surf.ctrl(end, :, :), hi.surf.ctrl(1, :, :)) && all(lo.surf.ctrl(end, :, 3) == z), 'cut row');
    for pc = {lo, hi}
        Z = pc{1}.surf.ctrl(:, :, 3);
        check(all(all(Z == Z(:, 1))), 'piece rows keep one z');
    end
    check(isequal(lo.z_range, [-2 z]) && isequal(hi.z_range, [z 1]) && lo.u_range(2) == us && hi.u_range(1) == us, 'ranges');
    L = mwecmass.solid.slice_bspline_surface(Q4, z);
    k1 = find([L.pieces.patch] == 1);
    check(isequal(L.pieces(k1).curve.ctrl, reshape(hi.surf.ctrl(1, :, :), [], 3)) && L.simple && ...
        numel(L.pieces) == 4, 'F4 row equals the F3b cut row');
end
% q alone is a quarter: its section is open
expect_error(@() mwecmass.solid.slice_bspline_surface(q, 0), 'mwecmass:solid:SectionNotClosed');

% errors
expect_error(@() mwecmass.solid.split_bspline_surface(q, 1), 'mwecmass:solid:ZOutside');
expect_error(@() mwecmass.solid.split_bspline_surface(q, -2.5), 'mwecmass:solid:ZOutside');
bad = q;
bad.z_of_u = false;
expect_error(@() mwecmass.solid.split_bspline_surface(bad, 0), 'mwecmass:solid:ZNotOneParameter');
bad = q;
bad.surf.ctrl(4, 2, 3) = 0;
expect_error(@() mwecmass.solid.split_bspline_surface(bad, 0), 'mwecmass:solid:ZNotOneParameter');
bump = q;
bump.surf.ctrl(:, :, 3) = repmat([-2; -1; 0.5; 0.5; -0.5; 0.6; 1], 1, 3);
expect_error(@() mwecmass.solid.split_bspline_surface(bump, 0), 'mwecmass:solid:ZNotMonotonic');
flat = q;
flat.surf.ctrl(:, :, 3) = repmat([-2; -1.5; -1; -1; -1; -1; 1], 1, 3);
expect_error(@() mwecmass.solid.split_bspline_surface(flat, -1), 'mwecmass:solid:ZNotMonotonic');
[~, ~, us] = mwecmass.solid.split_bspline_surface(flat, -1.2);
fprintf('row z [-2 -1.5 -1 -1 -1 -1 1]: z = -1 raises ZNotMonotonic (interval), z = -1.2 cut at u* = %.17g\n', us);
end

function dev = piece_deviation(parent, lo, hi)
dev = 0;
for pc = {lo, hi}
    p = pc{1};
    u = p.u_range(1) + (0:10)' / 10 * diff(p.u_range);
    [uu, vv] = meshgrid(u, linspace(0, 1, 5));
    A = mwecmass.solid.eval_bspline_surface(parent.surf, uu(:), vv(:));
    B = mwecmass.solid.eval_bspline_surface(p.surf, uu(:), vv(:));
    dev = max(dev, max(abs(A(:) - B(:))));
end
end

function expect_error(f, id)
try
    f();
catch err
    if ~strcmp(err.identifier, id)
        error('test_split_slice:fail', 'expected %s, got %s (%s)', id, err.identifier, err.message);
    end
    return
end
error('test_split_slice:fail', 'expected %s, got no error', id);
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'standins', 'fixtures'), '-end');
end

function check(cond, varargin)
if ~cond
    error('test_split_slice:fail', varargin{:});
end
end
