function test_offset_surface_rules()
%TEST_OFFSET_SURFACE_RULES  F2 on decks written here: a rounded rim that folds (fillet radius
%   0.0625 m < d), the S2r section count, the knots_from cache key, and the hand-off of concave
%   creases and a concave vertex to fit_z_faces (test-only recording stub, fit_z_faces_recorder).
%   Independent references: the deck's profile (distance in the meridian plane to its lines and
%   arc, exact for a surface of revolution) and the deck coordinates of the planar faces.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
tmp = tempname();
mkdir(tmp);
stub = fit_z_faces_recorder();
cleanup = onCleanup(@() leave(stub, tmp));
t_min = 0.0254;
o = struct('t_min', t_min);

% ---------------------------------------------------------------- fillet that folds
deck = {
    'MultiSurf 1.44'
    'Units: m kg'
    'Symmetry:  x y'
    'BeginModel;'
    'FramePoint K 14 -1 / 0 * * 0.0 0.0 -3.0 ;'
    'FramePoint F1 14 -1 / 0 * * 1.4375 0.0 -3.0 ;'
    'FramePoint CF 14 -1 / 0 * * 1.4375 0.0 -2.9375 ;'
    'FramePoint F2 14 -1 / 0 * * 1.5 0.0 -2.9375 ;'
    'FramePoint F3 14 -1 / 0 * * 1.5 0.0 1.0 ;'
    'FramePoint T 14 -1 / 0 * * 0.0 0.0 1.0 ;'
    'Line keel_line 11 -1 8x4 / * K F1 ;'
    'Arc fillet 11 -1 8x4 / * 2 F1 CF F2 ;'
    'Line wall_line 11 -1 8x4 / * F2 F3 ;'
    'Line top_line 11 -1 8x4 / * F3 T ;'
    'PolyCurve2 profile 11 -1 8x4 / * { keel_line fillet wall_line top_line } ;'
    'Line axis 6 -1 8x4 / * K T ;'
    'RevSurf hull 2 11 36x4 12x4 0 / * profile axis 0.0 90.0 ;'
    'EndModel;'};
f = fullfile(tmp, 'fillet.ms2');
write_text(f, sprintf('%s\n', deck{:}));
model = parse(f);
geo = mwecmass.solid.outer_nurbs(model);
check(all([geo.outer.exact]), 'fillet: outer patches exact');
rc = 0.0625;
dist = @(X) profile_distance([hypot(X(:, 1), X(:, 2)), X(:, 3)], rc);

t = 0.0762;
tic;
[inner, rep] = mwecmass.solid.offset_surface(model, [], geo, t, geo.z_range, o);
el = toc;
d = inner.d;
check(d > rc, 'fillet: d exceeds the fillet radius');
check(rep.ok, 'fillet t = %.4f: report not ok', t);
q1 = inner.patches([inner.patches.visible] == 1);
check(numel(q1) == 3, 'fillet t = %.4f: three pieces per quarter (bottom disk, wall, top disk), got %d', t, numel(q1));
poly = [0, -3 + d; 1.5 - d, -3 + d; 1.5 - d, 1 - d; 0, 1 - d];
n_on = 0;
for k = 1:3
    col = reshape(q1(k).surf.ctrl(:, 1, :), [], 3);
    check(all(col(:, 2) == 0), 'fillet: first column in y = 0');
    A = poly(k, :);
    B = poly(k + 1, :);
    % on the segment, bitwise: the constant coordinate equal, the other between the ends
    c = find(A == B, 1);
    other = 3 - c;
    rz = col(:, [1 3]);
    check(all(rz(:, c) == A(c)) && isequal(rz(1, :), A) && isequal(rz(end, :), B) && ...
        all(rz(:, other) >= min(A(other), B(other)) & rz(:, other) <= max(A(other), B(other))), ...
        'fillet t = %.4f: piece %s does not lie on the offset polyline', t, q1(k).name);
    n_on = n_on + size(col, 1);
end
fprintf(['fillet t = %.4f (d = %.17g > r = %.4g): the fold is trimmed, quarter 1 has 3 pieces; their %d first-column ' ...
    'control points lie bitwise on (0, -3 + d), (1.5 - d, -3 + d), (1.5 - d, 1 - d), (0, 1 - d); %.1f s\n'], ...
    t, d, rc, n_on, el);
st = inner_set_checks(geo, inner, t_min, dist, sprintf('fillet t = %.4f', t));
% the revolved points are rational evaluations (a few roundings of coordinates below 3 m) and the
% oracle subtracts coordinates below 3 m: 8 eps(3)
check(abs(st.tl_min - d) <= 8 * eps(3) && abs(st.tl_max - d) <= 8 * eps(3), ...
    'fillet t = %.4f: dense t_local %.17g .. %.17g, d = %.17g', t, st.tl_min, st.tl_max, d);
n_sec = check_heights(inner, 6);
check(rep.n_sections == n_sec, 'fillet: S2r n_sections %d, distinct check heights %d', rep.n_sections, n_sec);
fprintf('fillet t = %.4f: S2r judged M3 on %d sections, the number of distinct check heights strictly inside the set\n', ...
    t, rep.n_sections);

t = t_min;
[inner2, rep2] = mwecmass.solid.offset_surface(model, [], geo, t, geo.z_range, o);
check(rep2.ok, 'fillet t = %.4f: report not ok', t);
st2 = inner_set_checks(geo, inner2, t_min, dist, sprintf('fillet t = %.4f', t));
n_sec = check_heights(inner2, 6);
check(rep2.n_sections == n_sec, 'fillet: S2r n_sections %d, distinct check heights %d', rep2.n_sections, n_sec);
fprintf('fillet t = %.4f (no fold, d < r): dense t_local %.12f .. %.12f m, %d M3 sections = distinct check heights\n', ...
    t, st2.tl_min, st2.tl_max, rep2.n_sections);

% ---------------------------------------------------------------- knots_from: keyed by its knots
kf = inner2;
i = find(arrayfun(@(p) numel(unique(p.surf.knots{1})) > 2 && p.visible == 1, kf.patches), 1);
check(~isempty(i), 'fillet: a fitted piece');
ku = kf.patches(i).surf.knots{1};
u = unique(ku);
kf.patches(i).surf.knots{1} = sort([ku, (u(1) + u(2)) / 2]);
t2 = 0.03;
ra = mwecmass.solid.offset_surface(model, [], geo, t2, geo.z_range, struct('t_min', t_min, 'knots_from', inner2));
rb = mwecmass.solid.offset_surface(model, [], geo, t2, geo.z_range, struct('t_min', t_min, 'knots_from', kf));
check(kf.t == inner2.t && numel(kf.patches) == numel(inner2.patches), 'knots_from: same t and count');
ka = ra.patches(i).surf.knots{1};
kb = rb.patches(i).surf.knots{1};
check(ra.refit && rb.refit && numel(kb) == numel(ka) + 1, ...
    'knots_from: two sets with equal t and count but different knots give the same refit (%d and %d knots)', numel(ka), numel(kb));
ri = mwecmass.solid.offset_surface(model, [], geo, t_min, geo.z_range, struct('t_min', t_min, 'knots_from', inner2));
check(isequal(ri.patches(i).surf.knots, inner2.patches(i).surf.knots), 'knots_from: identity keeps the knots');
fprintf('knots_from at t = %.3f: %s refitted on %d and on %d knots (sets with equal t = %.4f and %d pieces, different knots)\n', ...
    t2, ra.patches(i).name, numel(ka), numel(kb), kf.t, numel(kf.patches));

% ---------------------------------------------------------------- hand-off of concave creases and a vertex
% quarter of a box [-2, 2]^2 x [0, 2] with a cube [1, 2]^3 cut from each top corner (Symmetry x y);
% the notch has three concave creases meeting at (1, 1, 1), whose cone of normals is an octant
model = parse(write_notch(tmp));
geo = mwecmass.solid.outer_nurbs(model);
check(all([geo.outer.exact]), 'notch: outer patches exact');
nm = @(k) geo.outer(k).name;
addpath(stub);
rec = run_to_stub(@() mwecmass.solid.offset_surface(model, [], geo, t_min, geo.z_range, o));
rmpath(stub);
F = rec.faces;
g = rec.general;
check(isequal(rec.opts.kind, 'inner') && rec.opts.d == t_min + 0.01 * t_min / 2 && isequal(rec.opts.geo, geo), ...
    'notch: fit_z_faces options');
H = F(g);
kinds = arrayfun(@(e) fieldnames(e.offset_of), H, 'UniformOutput', false);
kinds = cellfun(@(c) c{1}, kinds, 'UniformOutput', false);
cr = H(strcmp(kinds, 'crease'));
vx = H(strcmp(kinds, 'vertex'));
check(numel(H) == 4 && numel(cr) == 3 && numel(vx) == 1, 'notch: %d entries handed (%s)', numel(H), strjoin(kinds, ' '));
pairs = sort(cellfun(@(c) strjoin(sort({nm(c(1)), nm(c(2))}), '+'), arrayfun(@(e) e.offset_of.crease.outer, cr, ...
    'UniformOutput', false), 'UniformOutput', false));
want = sort({'floor+wall_x1', 'floor+wall_y1', 'wall_x1+wall_y1'});
check(isequal(pairs, want), 'notch: creases %s', strjoin(pairs, ', '));
for e = cr
    c = e.offset_of.crease.curve;
    k = e.offset_of.crease.outer;
    check(isempty(e.surf) && ~e.exact && e.visible == min([geo.outer(k).visible]) && ...
        any(all(c.ctrl == [1 1 1], 2)), 'notch: crease entry %s', e.name);
end
check(isequal(vx.offset_of.vertex.point, [1 1 1]) && isempty(vx.surf) && ...
    isequal(sort({geo.outer(vx.offset_of.vertex.outer).name}), sort({'floor', 'wall_x1', 'wall_y1'})), ...
    'notch: vertex entry');
prim = arrayfun(@(e) isempty(e.flips), geo.outer);
for e = H
    fld = fieldnames(e.offset_of);
    outer = e.offset_of.(fld{1}).outer;
    check(any(prim(outer)), 'notch: %s offsets only mirror entries (fit_z_faces makes the mirrors)', e.name);
end
for e = F(~g)
    check(~isempty(e.surf), 'notch: structured piece %s without a surface', e.name);
end
% the planar faces keep their structure: the notch floor (z = 1) moves down by d, bounded by the
% offset planes of the outer walls (convex) and by the planes over its concave creases
fl = F(strcmp({F.name}, 'floor_inner'));
check(numel(fl) == 1, 'notch: floor_inner');
dd = rec.opts.d;
Q = sortrows(reshape(fl.surf.ctrl, [], 3));
want = sortrows([1 1 1 - dd; 2 - dd 1 1 - dd; 1 2 - dd 1 - dd; 2 - dd 2 - dd 1 - dd]);
check(isequal(Q, want), 'notch: floor_inner corners %s', mat2str(Q, 17));
fprintf(['notch: fit_z_faces receives %d structured pieces and %d entries (primaries only): creases %s with ' ...
    'offset_of.crease, vertex (%g, %g, %g) with offset_of.vertex; floor_inner corners bitwise at z = 1 - d\n'], ...
    nnz(~g), nnz(g), strjoin(pairs, ', '), vx.offset_of.vertex.point);
end

function n = check_heights(inner, n_gauss)
% the distinct heights of the S2r check points (Gauss points and span midpoints in u) strictly
% inside the set and off its constant-z pieces, from the written patches
[xg, ~] = gauss_legendre(n_gauss);
P = inner.patches;
zz = reshape([P.z_range], 2, [])';
flatz = zz(zz(:, 1) == zz(:, 2), 1);
z = zeros(0, 1);
for k = 1:numel(P)
    if P(k).z_range(1) == P(k).z_range(2)
        continue
    end
    ku = unique(P(k).surf.knots{1});
    us = zeros(0, 1);
    for j = 1:numel(ku) - 1
        us = [us; ku(j) + (ku(j + 1) - ku(j)) * (xg + 1) / 2; (ku(j) + ku(j + 1)) / 2]; %#ok<AGROW>
    end
    X = mwecmass.solid.eval_bspline_surface(P(k).surf, us, repmat(P(k).surf.knots{2}(1), numel(us), 1));
    z = [z; X(:, 3)]; %#ok<AGROW>
end
z = unique(z(z > min(zz(:)) & z < max(zz(:)) & ~ismember(z, flatz)));
n = numel(z);
end

function f = write_notch(tmp)
P = struct('A00', [0 0 0], 'A10', [1 0 0], 'A20', [2 0 0], 'A01', [0 1 0], 'A11', [1 1 0], 'A21', [2 1 0], ...
    'A02', [0 2 0], 'A12', [1 2 0], 'A22', [2 2 0], 'B11', [1 1 1], 'B21', [2 1 1], 'B12', [1 2 1], ...
    'B22', [2 2 1], 'C00', [0 0 2], 'C10', [1 0 2], 'C20', [2 0 2], 'C01', [0 1 2], 'C11', [1 1 2], ...
    'C21', [2 1 2], 'C02', [0 2 2], 'C12', [1 2 2]);
L = {};
for n = fieldnames(P)'
    L{end + 1} = sprintf('FramePoint %s 14 -1 / 0 * * %.1f %.1f %.1f ;', n{1}, P.(n{1})); %#ok<AGROW>
end
% faces: name, first line (two points), second line; rulings join the first points and the last points
faces = {
    'bottom_a', 'A00', 'A10', 'A01', 'A11'
    'bottom_b', 'A10', 'A20', 'A11', 'A21'
    'bottom_c', 'A01', 'A11', 'A02', 'A12'
    'bottom_d', 'A11', 'A21', 'A12', 'A22'
    'top_a', 'C00', 'C10', 'C01', 'C11'
    'top_b', 'C10', 'C20', 'C11', 'C21'
    'top_c', 'C01', 'C11', 'C02', 'C12'
    'floor', 'B11', 'B21', 'B12', 'B22'
    'side_xa', 'A20', 'C20', 'A21', 'C21'
    'side_xb', 'A21', 'B21', 'A22', 'B22'
    'side_ya', 'A02', 'C02', 'A12', 'C12'
    'side_yb', 'A12', 'B12', 'A22', 'B22'
    'wall_x1', 'B11', 'C11', 'B12', 'C12'
    'wall_y1', 'B11', 'C11', 'B21', 'C21'};
for i = 1:size(faces, 1)
    L{end + 1} = sprintf('Line %s_1 11 -1 8x4 / * %s %s ;', faces{i, 1}, faces{i, 2}, faces{i, 3}); %#ok<AGROW>
    L{end + 1} = sprintf('Line %s_2 11 -1 8x4 / * %s %s ;', faces{i, 1}, faces{i, 4}, faces{i, 5}); %#ok<AGROW>
end
for i = 1:size(faces, 1)
    L{end + 1} = sprintf('RuledSurf %s 2 11 36x4 12x4 0 / * %s_1 %s_2 ;', faces{i, 1}, faces{i, 1}, faces{i, 1}); %#ok<AGROW>
end
f = fullfile(tmp, 'notch.ms2');
write_text(f, sprintf('%s\n', 'MultiSurf 1.44', 'Units: m kg', 'Symmetry:  x y', 'BeginModel;', L{:}, 'EndModel;'));
end

function dist = profile_distance(Q, rc)
% distance from points (r, z) to the deck profile: the Lines K-F1, F2-F3, F3-T and the fillet arc
% of radius rc about (1.5 - rc, -3 + rc) from F1 to F2 (exact, segment by segment)
segs = [0 -3 1.5 - rc -3; 1.5 -3 + rc 1.5 1; 1.5 1 0 1];
dist = Inf(size(Q, 1), 1);
for i = 1:size(segs, 1)
    A = segs(i, 1:2);
    AB = segs(i, 3:4) - A;
    s = min(max(((Q - A) * AB') / (AB * AB'), 0), 1);
    dist = min(dist, hypot(Q(:, 1) - A(1) - s * AB(1), Q(:, 2) - A(2) - s * AB(2)));
end
c = [1.5 - rc, -3 + rc];
w = Q - c;
ang = atan2(w(:, 2), w(:, 1));
on_arc = ang >= -pi / 2 & ang <= 0;
da = abs(hypot(w(:, 1), w(:, 2)) - rc);
dist(on_arc) = min(dist(on_arc), da(on_arc));
end

function rec = run_to_stub(f)
global FIT_Z_FACES_RECORD
FIT_Z_FACES_RECORD = [];
try
    f();
catch err
    if ~strcmp(err.identifier, 'fit_z_faces_stub:reached')
        error('test_offset_surface_rules:fail', 'expected the fit_z_faces stub, got %s (%s)', err.identifier, err.message);
    end
    rec = FIT_Z_FACES_RECORD;
    return
end
error('test_offset_surface_rules:fail', 'fit_z_faces was not called');
end

function [x, w] = gauss_legendre(n)
b = (1:n - 1) ./ sqrt(4 * (1:n - 1).^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D));
w = 2 * V(1, i)'.^2;
end

function leave(stub, tmp)
if any(strcmp(strsplit(path(), pathsep()), stub))
    rmpath(stub);
end
rmdir(stub, 's');
if exist(tmp, 'dir')
    rmdir(tmp, 's');
end
end

function m = parse(f)
evalc('m = mwecmass.geometry.MS2Parser.parse(f);');
end

function write_text(f, txt)
fid = fopen(f, 'w');
fprintf(fid, '%s', txt);
fclose(fid);
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'solid'));
end

function check(cond, varargin)
if ~cond
    error('test_offset_surface_rules:fail', varargin{:});
end
end
