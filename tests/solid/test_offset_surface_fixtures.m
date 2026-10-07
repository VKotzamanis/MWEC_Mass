function test_offset_surface_fixtures()
%TEST_OFFSET_SURFACE_FIXTURES  F2 offset_surface and F2b void_closing_distance on the SK cylinder, the SK box and stepped_spar.
%   Independent references: the closed forms of the SK kit (sti_closed_form, sti_inner_box: the
%   normal offset of a plane is the parallel plane, of a cylinder the coaxial cylinder) and, for
%   stepped_spar, the distance in the meridian plane to the deck's straight profile.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
sk = fullfile(root, 'tests', 'standins', 'fixtures');
fix = fullfile(root, 'tests', 'solid', 'fixtures');
t_min = 0.0254;

% ---------------------------------------------------------------- cylinder (SK join fixture)
model = parse(fullfile(sk, 'cylinder.ms2'));
geo = mwecmass.solid.outer_nurbs(model);
fx = sti_closed_form('fixture', 'cylinder');
cyl_dist = @(X) min([fx.R - hypot(X(:, 1), X(:, 2)), X(:, 3) - fx.z(1), fx.z(2) - X(:, 3)], [], 2);
for t = [t_min 0.0762]
    tic;
    inner = mwecmass.solid.offset_surface(model, [], geo, t, geo.z_range, struct('t_min', t_min));
    t1 = toc;
    d = inner.d;
    % the closed form (SK stand-in recipe): each quarter of the coaxial cylinder of radius R - d
    % between z0 + d and z1 - d, cut at its rims into bottom disk, side and top disk
    full = sti_closed_form('patches', fx, d);
    rows = {[1 2], [2 3], [3 4]};
    u = {[0 1/3], [1/3 2/3], [2/3 1]};
    check(numel(inner.patches) == 3 * numel(full), 'cylinder: three crease pieces per quarter');
    for k = 1:numel(full)
        f = full(k);
        for p = 1:3
            q = inner.patches(3 * (k - 1) + p);
            check(strcmp(q.name, sprintf('%s_inner%d', f.name, p)) && q.visible == k, 'cylinder: piece name and visible');
            check(isequal(q.surf.ctrl, f.surf.ctrl(rows{p}, :, :)) && isequal(q.surf.weights, f.surf.weights(rows{p}, :)) && ...
                isequal(q.surf.knots{1}, u{p}([1 1 2 2])) && isequal(q.surf.knots{2}, f.surf.knots{2}), ...
                'cylinder: %s differs from the closed form', q.name);
            check(isequal(q.seam_v0, [3 * (f.seam_v0(1) - 1) + p, 1]) && isequal(q.seam_v1, [3 * (f.seam_v1(1) - 1) + p, 3]), ...
                'cylinder: %s v-seams', q.name);
            check(q.outward == ~f.outward && isequal(q.pole, [p == 1, p == 3]) && ~q.exact, 'cylinder: %s flags', q.name);
        end
    end
    check(isequal(inner.z_range, [fx.z(1) + d, fx.z(2) - d]) && inner.z_lo == fx.z(1) + d && ~inner.refit, ...
        'cylinder: z_range of the void');
    fprintf('cylinder t = %.4f: inner set equals the closed form bitwise (R - d = %.17g, rims trimmed at z = %.17g, %.17g), %.1f s\n', ...
        t, fx.R - d, fx.z(1) + d, fx.z(2) - d, t1);
    inner_set_checks(geo, inner, t_min, cyl_dist, sprintf('cylinder t = %.4f', t));
end
% cached per thickness: the same set again, without rebuilding
tic;
again = mwecmass.solid.offset_surface(model, [], geo, 0.0762, geo.z_range, struct('t_min', t_min));
t2 = toc;
check(isequal(again, inner), 'cylinder: cached set differs');
fprintf('cylinder: second call at t = 0.0762 returned the cached set in %.3f s (first %.1f s)\n', t2, t1);
% a partial range (precast: the void ends at the bottom of a solid module): open at its top
zc = [fx.z(1), 0];
part = mwecmass.solid.offset_surface(model, [], geo, t_min, zc, struct('t_min', t_min));
zz = reshape([part.patches.z_range], 2, [])';
check(max(zz(:)) >= zc(2) && part.z_lo == fx.z(1) + part.d, 'cylinder: partial set does not cover z_range');
check(~any(zz(:, 1) == zz(:, 2) & zz(:, 1) > 0), 'cylinder: partial set keeps the top disk');
fprintf('cylinder z_range [%g %g]: %d patches up to z = %.17g (open end above %g)\n', zc, numel(part.patches), max(zz(:)), zc(2));
inner_set_checks(geo, part, t_min, cyl_dist, 'cylinder partial', struct('z_clip', zc));
% F2b and the bound it sets
dc = mwecmass.solid.void_closing_distance(model, [], geo, geo.z_range);
dref = sti_closed_form('d_close', fx);
fprintf('cylinder: d_close %.17g m, closed form min(R, H/2) = %.17g m, t_max = %.17g m\n', dc, dref, dc - 0.01 * t_min / 2);
% two rational evaluations of the chord ends and one norm, each a few correctly rounded operations
% of values below 4 m: within 8 ulp of 1.5
check(abs(dc - dref) <= 8 * eps(dref), 'cylinder: d_close');
expect_error(@() mwecmass.solid.offset_surface(model, [], geo, dc, geo.z_range, struct('t_min', t_min)), ...
    'mwecmass:solid:VoidClosed');
expect_error(@() mwecmass.solid.offset_surface(model, [], geo, t_min, geo.z_range, struct()), ...
    'mwecmass:solid:FitInputMissing');
fprintf('cylinder: VoidClosed at t = d_close, FitInputMissing without t_min\n');

% ---------------------------------------------------------------- box: equals sti_inner_box
model = parse(fullfile(sk, 'box.ms2'));
geo = mwecmass.solid.outer_nurbs(model);
fb = sti_closed_form('fixture', 'box');
box_dist = @(X) min([X(:, 1) - fb.x(1), fb.x(2) - X(:, 1), X(:, 2) - fb.y(1), fb.y(2) - X(:, 2), ...
    X(:, 3) - fb.z(1), fb.z(2) - X(:, 3)], [], 2);
for t = [t_min 0.0762]
    inner = mwecmass.solid.offset_surface(model, [], geo, t, geo.z_range, struct('t_min', t_min));
    ref = sti_inner_box(geo, t, t_min);
    check(numel(inner.patches) == numel(ref.patches), 'box: patch count');
    dmax = 0;
    for k = 1:numel(ref.patches)
        a = inner.patches(k);
        b = ref.patches(k);
        dmax = max(dmax, max(abs(a.surf.ctrl(:) - b.surf.ctrl(:))));
        check(isequal(a.surf.ctrl, b.surf.ctrl) && isequal(a.surf.knots, b.surf.knots) && strcmp(a.name, b.name) && ...
            a.outward == b.outward, 'box: %s differs from sti_inner_box', a.name);
        for f = {'seam_u0', 'seam_u1', 'seam_v0', 'seam_v1'}
            check(isequal(a.(f{1}), b.(f{1})), 'box: %s %s differs from sti_inner_box', a.name, f{1});
        end
    end
    check(isequal(inner.z_range, ref.z_range) && inner.z_lo == ref.z_lo, 'box: z_range');
    fprintf('box t = %.4f: real F2 equals sti_inner_box (control points bitwise, max difference %g m; seams, names, orientation)\n', t, dmax);
    inner_set_checks(geo, inner, t_min, box_dist, sprintf('box t = %.4f', t));
end
dc = mwecmass.solid.void_closing_distance(model, [], geo, geo.z_range);
dref = sti_closed_form('d_close', fb);
fprintf('box: d_close %.17g m, closed form (half the smallest dimension) %.17g m\n', dc, dref);
% planar faces: the critical chord is found exactly
check(dc == dref, 'box: d_close');

% ---------------------------------------------------------------- stepped_spar: convex rim, concave corner, shelf
model = parse(fullfile(fix, 'stepped_spar.ms2'));
geo = mwecmass.solid.outer_nurbs(model);
prof = [0 -3; 1.5 -3; 1.5 -1; 0.75 -1; 0.75 1; 0 1];
spar_dist = @(X) polyline_distance([hypot(X(:, 1), X(:, 2)), X(:, 3)], prof);
t = t_min;
inner = mwecmass.solid.offset_surface(model, [], geo, t, geo.z_range, struct('t_min', t_min));
d = inner.d;
P = inner.patches;
check(numel(P) == 24, 'stepped_spar: six inner pieces per quarter');
expect = [-3 + d, -3 + d; -3 + d, -1 - d; -1 - d, -1 - d; -1 - d, -1; -1, 1 - d; 1 - d, 1 - d];
for q = 1:4
    i = 6 * (q - 1) + (1:6);
    check(isequal(reshape([P(i).z_range], 2, [])', expect), 'stepped_spar quarter %d: z ranges of the pieces', q);
end
% convex rim P2: the base wall and the shelf offsets are trimmed where they cross, (1.5 - d, -1 - d)
base = P(2).surf.ctrl;
shelf = P(3).surf.ctrl;
check(all(all(base(end, :, 3) == -1 - d)) && isequal(base(end, :, :), shelf(1, :, :)), ...
    'stepped_spar: base wall and shelf meet at the trimmed rim');
r_rim = hypot(base(end, 1, 1), base(end, 1, 2));
% concave corner P3: a face of its own, the crease circle offset by d along its fan of normals
fan = P(4);
check(fan.surf.degree(1) == 2 && ~isempty(fan.surf.weights), 'stepped_spar: fan face is a rational arc in u');
[UU, VV] = ndgrid(linspace(fan.u_range(1), fan.u_range(2), 41), linspace(0, 1, 21));
X = mwecmass.solid.eval_bspline_surface(fan.surf, UU(:), VV(:) * fan.surf.knots{2}(end));
dev = max(abs(hypot(hypot(X(:, 1), X(:, 2)) - 0.75, X(:, 3) + 1) - d));
fprintf(['stepped_spar t = %.4f (d = %.17g): rim P2 trimmed at r = %.17g (1.5 - d = %.17g), shelf at z = -1 - d, ' ...
    'fan face about P3: max | |X - crease| - d | = %.3g m\n'], t, d, r_rim, 1.5 - d, dev);
% a rational quadratic evaluation and two hypot calls on values below 2 m: a few ulp of 2
check(abs(r_rim - (1.5 - d)) <= 4 * eps(2) && dev <= 16 * eps(2), 'stepped_spar: rim and fan face');
st = inner_set_checks(geo, inner, t_min, spar_dist, 'stepped_spar');
% F4 inside the fan band and the inner volume (Pappus on the profile of the written faces)
zf = -1 - d / 2;
L = mwecmass.solid.slice_bspline_surface(P(arrayfun(@(e) e.z_range(1) ~= e.z_range(2), P)), zf);
rf = hypot(L.pts(:, 1), L.pts(:, 2));
fprintf('stepped_spar: section at z = -1 - d/2: radius %.17g .. %.17g, closed form 0.75 - sqrt(3)/2 d = %.17g\n', ...
    min(rf), max(rf), 0.75 - sqrt(3) / 2 * d);
V = pappus(P(1:6));
a = 0.75;
Vref = pi * (1.5 - d)^2 * (2 - 2 * d) + pi * (a^2 * d - a * pi * d^2 / 2 + 2 * d^3 / 3) + pi * (a - d)^2 * (2 - d);
fprintf('stepped_spar: inner volume %.15f m3 (Gauss on the written faces), closed form of the offset profile %.15f m3, difference %.3g m3\n', ...
    V, Vref, V - Vref);
dc = mwecmass.solid.void_closing_distance(model, [], geo, geo.z_range);
fprintf('stepped_spar: d_close %.17g m (column radius 0.75), %d dense points checked\n', dc, st.n_dense);
check(abs(dc - 0.75) <= 8 * eps(0.75), 'stepped_spar: d_close');
end

function V = pappus(P)
% volume of revolution of the profile of one quarter's faces: |sum of the integrals of pi r^2 dz| along the first column
[xg, wg] = gauss(12);
V = 0;
for k = 1:numel(P)
    s = P(k).surf;
    ku = unique(s.knots{1});
    for j = 1:numel(ku) - 1
        uu = ku(j) + (ku(j + 1) - ku(j)) * (xg + 1) / 2;
        [X, Xu] = mwecmass.solid.eval_bspline_surface(s, uu, repmat(s.knots{2}(1), numel(uu), 1));
        V = V + (ku(j + 1) - ku(j)) / 2 * sum(wg .* pi .* (X(:, 1).^2 + X(:, 2).^2) .* Xu(:, 3));
    end
end
V = abs(V);
end

function dist = polyline_distance(Q, prof)
% distance from points (r, z) to the profile polyline (segment by segment, exact)
dist = Inf(size(Q, 1), 1);
for i = 1:size(prof, 1) - 1
    A = prof(i, :);
    B = prof(i + 1, :);
    if A(1) == 0 && B(1) == 0
        continue
    end
    AB = B - A;
    s = min(max(((Q - A) * AB') / (AB * AB'), 0), 1);
    dist = min(dist, hypot(Q(:, 1) - A(1) - s * AB(1), Q(:, 2) - A(2) - s * AB(2)));
end
end

function [x, w] = gauss(n)
b = (1:n - 1) ./ sqrt(4 * (1:n - 1).^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D));
w = 2 * V(1, i)'.^2;
end

function m = parse(f)
evalc('m = mwecmass.geometry.MS2Parser.parse(f);');
end

function expect_error(f, id)
try
    f();
catch err
    if ~strcmp(err.identifier, id)
        error('test_offset_surface_fixtures:fail', 'expected %s, got %s (%s)', id, err.identifier, err.message);
    end
    return
end
error('test_offset_surface_fixtures:fail', 'expected %s, got no error', id);
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'solid'));
addpath(fullfile(root, 'tests', 'standins', 'fixtures'), '-end');
end

function check(cond, varargin)
if ~cond
    error('test_offset_surface_fixtures:fail', varargin{:});
end
end
