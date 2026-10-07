function test_trim_fold()
%TEST_TRIM_FOLD  trim_fold on a smooth fold: a V of two lines joined by a fillet arc of radius R, offset into the V by d.
%   Independent reference (closed form): for d > R the arc's offset is a reversed arc (a swallowtail)
%   and the two line offsets y = |x| + d sqrt(2) cross at (0, d sqrt(2)); for d < R nothing is trimmed.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
R = 0.05;
for d = [0.03 0.08]
    [F, rng] = chain(R, d);
    keep = mwecmass.solid.trim_fold(F, rng);
    if d < R
        check(numel(keep) == 3 && isequal([keep.k], 1:3) && ~any([keep.trim0]) && ~any([keep.trim1]), ...
            'd = %g < R: nothing trimmed', d);
        fprintf('trim_fold: d = %g < R = %g: three parts kept, no trim\n', d, R);
    else
        X = keep(1).X1;
        Xref = [0, d * sqrt(2)];
        check(numel(keep) == 2 && isequal([keep.k], [1 3]) && keep(1).trim1 && keep(2).trim0 && ...
            isequal(keep(1).X1, keep(2).X0), 'd = %g > R: the fold is not trimmed at one crossing', d);
        % Newton on two exact lines ends at adjacent doubles; the mean of the two evaluations adds one rounding
        check(max(abs(X - Xref)) <= 4 * eps(1), 'd = %g: crossing %s, closed form %s', d, mat2str(X, 17), mat2str(Xref, 17));
        fprintf('trim_fold: d = %g > R = %g: fillet offset removed, crossing (%.3g, %.17g), closed form (0, %.17g)\n', ...
            d, R, X(1), X(2), Xref(2));
    end
end
end

function [F, rng] = chain(R, d)
% left line from (-1, 1) to the fillet, fillet arc about (0, R sqrt(2)), right line to (1, 1);
% each offset by d along the normal into the V
s2 = sqrt(2);
A = [-1 1];
T1 = [-R / s2, R / s2];
T2 = [R / s2, R / s2];
B = [1 1];
Cc = [0, R * s2];
n1 = [1 1] / s2;
n3 = [-1 1] / s2;
F = {@(s) line_off(A, T1, n1, d, s), @(th) arc_off(Cc, R - d, th), @(s) line_off(T2, B, n3, d, s)};
rng = [0 1; -3 * pi / 4, -pi / 4; 0 1];
end

function [Q, dQ] = line_off(P0, P1, n, d, s)
s = s(:);
Q = P0 + s * (P1 - P0) + d * n;
dQ = repmat(P1 - P0, numel(s), 1);
end

function [Q, dQ] = arc_off(Cc, r, th)
th = th(:);
Q = Cc + r * [cos(th), sin(th)];
dQ = r * [-sin(th), cos(th)];
end

function check(cond, varargin)
if ~cond
    error('test_trim_fold:fail', varargin{:});
end
end
