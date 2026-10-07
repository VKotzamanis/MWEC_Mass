function test_sk_offset_slice()
%TEST_SK_OFFSET_SLICE  Stand-ins F2, F2b, F3b, F4 and sti_inner_box against the fixture closed forms.
%   Reads geo.analytic (stand-in marker of F1) and asserts stand-in exactness: J1, which merges the
%   real F1-F4, deletes or rewrites this test.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
evalc('mc = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''tests'', ''standins'', ''fixtures'', ''cylinder.ms2''));');
evalc('mb = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''tests'', ''standins'', ''fixtures'', ''box.ms2''));');
gc = mwecmass.solid.outer_nurbs(mc);
gb = mwecmass.solid.outer_nurbs(mb);
t_min = 0.1;
opts = struct('t_min', t_min);

% F4 on outer and inner patches against the closed-form sections
worst = 0;
for g = {gc, gb}
    geo = g{1};
    fx = geo.analytic;
    for z = [fx.z(1) + 0.25, mean(fx.z), fx.z(2) - 0.5]
        L = mwecmass.solid.slice_bspline_surface(geo.outer, z);
        worst = max(worst, compare(L, sti_closed_form('section', fx, 0)));
        check(numel(L.pieces) == 4 && L.simple && size(L.pts, 1) == 32, '%s: loop at z = %g', fx.name, z);
        check(all(abs(L.pts(:, 3) - z) <= 4 * eps(abs(z) + 1)), '%s: points off z', fx.name);
    end
end
fprintf('F4 vs closed-form sections (area, centroid, I): largest relative difference %.3e\n', worst);
% both are closed forms of the same section; their difference is rounding of a few dozen operations
check(worst <= 64 * eps, 'F4 differs from the closed form beyond rounding');

% F2 on the cylinder
t = 0.1;
[inner, rep] = mwecmass.solid.offset_surface(mc, [], gc, t, [-3 1], opts);
d = t + 0.01 * t_min / 2;
check(inner.d == d && inner.eps_fit == 0.01 * t_min && numel(inner.patches) == 12, 'F2 set fields');
check(isequal(sort(fieldnames(inner))', sort({'t', 'd', 'eps_fit', 'z_range', 'z_lo', 'refit', 'patches', 'report'})), 'S2 fields');
check(isequal(rep, inner.report) && rep.ok && ~rep.cap_reached && all([rep.patches.t_local_min] == d), 'S2r report');
check(isequal(inner.z_range, [-3 + d, 1 - d]) && inner.z_lo == -3 + d, 'S2 z range');
Li = mwecmass.solid.slice_bspline_surface(inner.patches, -1);
fprintf('F2 cylinder at t = %.3f m: d = %.17g, inner area %.17g, closed form %.17g\n', t, d, Li.area, ...
    sti_closed_form('section', gc.analytic, d).A);
check(compare(Li, sti_closed_form('section', gc.analytic, d)) <= 64 * eps, 'inner section');
for k = 1:numel(inner.patches)
    q = inner.patches(k);
    [S, Su, Sv] = mwecmass.solid.eval_bspline_surface(q.surf, mean(q.u_range), 0.5);
    toward_axis = -[S(1:2) 0];
    if q.z_range(1) == q.z_range(2)
        toward_axis = [0 0 sign(mean(gc.analytic.z) - S(3))];
    end
    check(sign(dot(cross(Su, Sv), toward_axis)) == 2 * q.outward - 1, '%s: normal not into the void', q.name);
    fields = {'seam_v0', 'seam_u1', 'seam_v1', 'seam_u0'};
    for b = 1:4
        nb = q.(fields{b});
        if isempty(nb)
            check(b == 4 && q.pole(1) || b == 2 && q.pole(2), '%s: empty seam %d is not a pole', q.name, b);
            continue
        end
        check(isequal(inner.patches(nb(1)).(fields{nb(2)}), [k b]), '%s: seam %d not mutual', q.name, b);
        [ra, wa] = boundary(q.surf, b);
        [rb, wb] = boundary(inner.patches(nb(1)).surf, nb(2));
        check(isequal(ra, rb) && isequal(wa, wb), '%s: seam %d rows differ', q.name, b);
    end
end
[~, rep2] = mwecmass.solid.offset_surface(mc, [], gc, 0.2, [-2 0], setfield(opts, 'knots_from', inner));
check(rep2.ok, 'knots_from run');
d_close = mwecmass.solid.void_closing_distance(mc, [], gc, [-3 1]);
fprintf('F2b: cylinder d_close %.17g (R = 1.5, H/2 = 2), box %.17g\n', d_close, ...
    mwecmass.solid.void_closing_distance(mb, [], gb, [-2.5 0.5]));
check(d_close == 1.5 && mwecmass.solid.void_closing_distance(mb, [], gb, [-2.5 0.5]) == 0.75, 'F2b closed form');
expect_error(@() mwecmass.solid.offset_surface(mc, [], gc, 1.5, [-3 1], opts), 'mwecmass:solid:VoidClosed');
expect_error(@() mwecmass.solid.offset_surface(mb, [], gb, 0.1, [-2.5 0.5], opts), 'mwecmass:solid:FitNotConverged');

% sti_inner_box
ib = sti_inner_box(gb, t, t_min);
Lb = mwecmass.solid.slice_bspline_surface(ib.patches, -1);
check(compare(Lb, sti_closed_form('section', gb.analytic, d)) <= 64 * eps, 'inner box section');
check(numel(ib.patches) == 6 && all([ib.patches.outward] == ~[gb.outer.outward]), 'inner box patches');
expect_error(@() sti_inner_box(gb, 0.75, t_min), 'mwecmass:solid:VoidClosed');
fprintf('sti_inner_box at t = %.3f m: area %.17g = (2 - 2d)(1.5 - 2d) = %.17g\n', t, Lb.area, (2 - 2 * d) * (1.5 - 2 * d));

% F3b: cut rows are shared bitwise, the pieces reproduce the parent at the same (u, v)
worst = 0;
for p = [gc.outer(1), gb.outer(2)]
    for z = [-2.2, -0.3]
        [lo, hi, us] = mwecmass.solid.split_bspline_surface(p, z);
        check(isequal(lo.surf.ctrl(end, :, :), hi.surf.ctrl(1, :, :)), '%s: cut row not shared', p.name);
        check(isempty(lo.surf.weights) && isempty(hi.surf.weights) || isequal(lo.surf.weights(end, :), hi.surf.weights(1, :)), '%s: cut weights not shared', p.name);
        check(lo.u_range(2) == us && hi.u_range(1) == us && lo.z_range(2) == z && hi.z_range(1) == z, 'F3b ranges');
        check(all(abs(lo.surf.ctrl(end, :, 3) - z) <= eps(abs(z) + 4)), 'cut row z');
        v = linspace(0, 1, 5)';
        for part = {lo, hi}
            q = part{1};
            u = q.u_range(1) + [0.1; 0.5; 0.9] * diff(q.u_range);
            [U, Vv] = meshgrid(u, v);
            A = mwecmass.solid.eval_bspline_surface(p.surf, U(:), Vv(:));
            B = mwecmass.solid.eval_bspline_surface(q.surf, U(:), Vv(:));
            worst = max(worst, max(abs(A(:) - B(:))));
        end
    end
end
fprintf('F3b: largest |parent - piece| at the same (u, v): %.3e m\n', worst);
% knot insertion does not move the surface; the difference is rounding on coordinates <= 3 m
check(worst <= 64 * eps * 3, 'F3b pieces differ from the parent beyond rounding');
expect_error(@() mwecmass.solid.split_bspline_surface(gc.outer(1), 1), 'mwecmass:solid:ZOutside');
expect_error(@() mwecmass.solid.split_bspline_surface(gb.outer(1), -2.5), 'mwecmass:solid:ZOutside');
bad = gb.outer(2);
bad.surf.ctrl = cat(1, bad.surf.ctrl, bad.surf.ctrl(1, :, :));
bad.surf.ctrl(3, :, 3) = 0;
bad.surf.knots{1} = [0 0 0.5 1 1];
expect_error(@() mwecmass.solid.split_bspline_surface(bad, -1), 'mwecmass:solid:ZNotMonotonic');
tilt = gb.outer(2);
tilt.surf.ctrl(2, 2, 3) = 0;
expect_error(@() mwecmass.solid.split_bspline_surface(tilt, -1), 'mwecmass:solid:ZNotOneParameter');
fprintf('F3b errors ZOutside, ZNotMonotonic, ZNotOneParameter raised\n');
end

function r = compare(L, s)
ref = [s.A, s.Sx / s.A, s.Sy / s.A, s.Ixx, s.Iyy, s.Ixy];
got = [L.area, L.centroid, L.I];
scale = [s.A, sqrt(s.A), sqrt(s.A), s.Ixx + s.Iyy, s.Ixx + s.Iyy, s.Ixx + s.Iyy];
r = max(abs(got - ref) ./ scale);
end

function [r, w] = boundary(s, b)
W = s.weights;
switch b
    case 1, r = squeeze(s.ctrl(:, 1, :)); w = W(:, 1);
    case 2, r = squeeze(s.ctrl(end, :, :)); w = W(end, :)';
    case 3, r = squeeze(s.ctrl(:, end, :)); w = W(:, end);
    case 4, r = squeeze(s.ctrl(1, :, :)); w = W(1, :)';
end
end

function expect_error(f, id)
try
    f();
catch err
    check(strcmp(err.identifier, id), 'expected %s, got %s: %s', id, err.identifier, err.message);
    return
end
error('test_sk_offset_slice:fail', 'expected error %s', id);
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
    error('test_sk_offset_slice:fail', varargin{:});
end
end
