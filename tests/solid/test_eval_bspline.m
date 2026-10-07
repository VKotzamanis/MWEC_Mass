function test_eval_bspline()
%TEST_EVAL_BSPLINE  F3: eval_bspline_curve / eval_bspline_surface against exact polynomial and circle identities.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);

% cubic Bezier of x = s^3, y = s, z = s^2 (power basis converted to Bezier control points)
cub = struct('degree', 3, 'ctrl', [0 0 0; 0 1/3 0; 0 2/3 1/3; 1 1 1], 'knots', [0 0 0 0 1 1 1 1], 'weights', []);
s = linspace(0, 1, 41)';
[C, Cs, Css] = mwecmass.solid.eval_bspline_curve(cub, s);
err = max(max(abs([C - [s.^3 s s.^2], Cs - [3 * s.^2 ones(size(s)) 2 * s], Css - [6 * s zeros(size(s)) 2 * ones(size(s))]])));
fprintf('cubic Bezier: largest deviation from s^3, s, s^2 and derivatives %.3e\n', err);
% a few rounded operations on values <= 6 per term
check(err <= 64 * eps * 6, 'cubic Bezier derivatives');

% the same cubic after inserting knots 0.3 and 0.7 by hand (Boehm), now a 3-span B-spline
k2 = [0 0 0 0 0.25 0.5 0.75 1 1 1 1];
% control points of s^3, s, s^2 on k2 by interpolation at Greville points is not exact; use
% the polar forms instead: for f(s) = s^k the control point i is the blossom at (t_{i+1}, t_{i+2}, t_{i+3})
P = zeros(7, 3);
for i = 1:7
    a = k2(i + 1:i + 3);
    P(i, :) = [prod(a), mean(a), (a(1) * a(2) + a(1) * a(3) + a(2) * a(3)) / 3];
end
spl = struct('degree', 3, 'ctrl', P, 'knots', k2, 'weights', []);
[C, Cs, Css] = mwecmass.solid.eval_bspline_curve(spl, s);
err = max(max(abs([C - [s.^3 s s.^2], Cs - [3 * s.^2 ones(size(s)) 2 * s], Css - [6 * s zeros(size(s)) 2 * ones(size(s))]])));
fprintf('cubic B-spline (blossom control points): largest deviation %.3e\n', err);
check(err <= 64 * eps * 6 / 0.25^2, 'cubic B-spline derivatives');

% rational quadratic 90-degree arc of radius 2 about (1, -1): |C - c| = r, (C - c).C' = 0 and
% C'.C' + (C - c).C'' = 0 hold exactly for a circle (derivatives of |C - c|^2 = r^2)
c = [1 -1 0.5];
r = 2;
arc = struct('degree', 2, 'ctrl', [c + [r 0 0]; c + [r r 0]; c + [0 r 0]], 'knots', [0 0 0 1 1 1], ...
    'weights', [1; sqrt(2) / 2; 1]);
[C, Cs, Css] = mwecmass.solid.eval_bspline_curve(arc, s);
e1 = max(abs(sqrt(sum((C - c).^2, 2)) - r));
e2 = max(abs(sum((C - c) .* Cs, 2)));
e3 = max(abs(sum(Cs .* Cs, 2) + sum((C - c) .* Css, 2)));
fprintf('rational arc: radius %.3e, (C-c).C'' %.3e, C''.C'' + (C-c).C'''' %.3e\n', e1, e2, e3);
check(e1 <= 64 * eps * 4 && e2 <= 64 * eps * 16 && e3 <= 256 * eps * 64, 'rational arc identities');

% surface: tensor product of the cubic in u and the arc in v (rational), against products
[uu, vv] = meshgrid(linspace(0, 1, 9), linspace(0, 1, 7));
uu = uu(:);
vv = vv(:);
ctrl = zeros(4, 3, 3);
W = zeros(4, 3);
for i = 1:4
    for j = 1:3
        ctrl(i, j, :) = [arc.ctrl(j, 1:2), cub.ctrl(i, 1)];
        W(i, j) = arc.weights(j);
    end
end
srf = struct('type', 'bspline', 'degree', [3 2], 'ctrl', ctrl, 'knots', {{cub.knots, arc.knots}}, 'weights', W);
[S, Su, Sv, Suu, Suv, Svv] = mwecmass.solid.eval_bspline_surface(srf, uu, vv);
[A, As, Ass] = mwecmass.solid.eval_bspline_curve(arc, vv);
z = uu.^3;
ref = {[A(:, 1:2) z], [0 * A(:, 1:2) 3 * uu.^2], [As(:, 1:2) 0 * z], [0 * A(:, 1:2) 6 * uu], ...
    zeros(numel(uu), 3), [Ass(:, 1:2) 0 * z]};
got = {S, Su, Sv, Suu, Suv, Svv};
worst = 0;
for k = 1:6
    worst = max(worst, max(max(abs(got{k} - ref{k}))));
end
fprintf('surface (cubic x arc): largest deviation of S and derivatives from the curve products %.3e\n', worst);
check(worst <= 256 * eps * 64, 'surface derivatives');

% non-rational surface: bilinear patch evaluates to the bilinear form
bil = struct('type', 'bspline', 'degree', [1 1], 'ctrl', cat(3, [0 1; 2 3], [0 0; 1 1], [0 1; 0 1]), ...
    'knots', {{[0 0 1 1], [0 0 1 1]}}, 'weights', []);
[S, Su, Sv] = mwecmass.solid.eval_bspline_surface(bil, uu, vv);
ref = [(1 - uu) .* vv + uu .* (2 + vv), uu, vv];
err = max(max(abs(S - ref)));
check(err <= 16 * eps * 4 && all(abs(Su(:, 2) - 1) <= 4 * eps) && all(abs(Sv(:, 3) - 1) <= 4 * eps), 'bilinear patch');
fprintf('bilinear patch: largest deviation %.3e\n', err);
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
end

function check(cond, msg)
if ~cond
    error('test_eval_bspline:fail', '%s', msg);
end
end
