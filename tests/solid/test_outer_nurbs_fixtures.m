function test_outer_nurbs_fixtures()
%TEST_OUTER_NURBS_FIXTURES  F1 on the SK fixtures, box_swapped, split_wall_box, stepped_spar and the deck check.
%   The SK closed forms (sti_closed_form 'patches': the exact NURBS of cylinder.ms2 and box.ms2) and
%   the deck coordinates of the T2a fixtures serve as the independent reference.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
fix = fullfile(root, 'tests', 'solid', 'fixtures');
sk = fullfile(root, 'tests', 'standins', 'fixtures');

% SK fixtures: the real F1 equals the closed-form S1 patches field by field (bitwise)
common = {'name', 'source', 'type', 'flips', 'surf', 'outward', 'exact', 'z_of_u', 'u_range', 'z_range', ...
    'offset_kind', 'pole', 'c0_u', 'c0_v', 'seam_u0', 'seam_u1', 'seam_v0', 'seam_v1'};
for name = {'cylinder', 'box'}
    model = parse(fullfile(sk, [name{1} '.ms2']));
    geo = mwecmass.solid.outer_nurbs(model);
    ref = sti_closed_form('patches', sti_closed_form('fixture', name{1}), 0);
    check(numel(geo.outer) == numel(ref), '%s: entry count', name{1});
    for k = 1:numel(ref)
        for f = common
            a = geo.outer(k).(f{1});
            b = ref(k).(f{1});
            if strcmp(f{1}, 'surf')
                a = rmfield(a, 'type');
                b = rmfield(b, 'type');
            end
            check(isequal(a, b) || (isempty(a) && isempty(b)), '%s: %s.%s differs from the closed form', ...
                name{1}, ref(k).name, f{1});
        end
        check(geo.outer(k).visible == k && ~geo.outer(k).swap_uv && isempty(geo.outer(k).fit), '%s: added fields', name{1});
    end
    fprintf('%s: %d entries equal the closed-form S1 patches bitwise\n', name{1}, numel(ref));
end

% box_swapped: z of the side faces depends on v alone; F1 exchanges u and v
mb = parse(fullfile(sk, 'box.ms2'));
ms = parse(fullfile(fix, 'box_swapped.ms2'));
gb = mwecmass.solid.outer_nurbs(mb);
gs = mwecmass.solid.outer_nurbs(ms);
nsw = 0;
for k = 1:numel(gs.outer)
    e = gs.outer(k);
    ent = ms.entities(e.name);
    c1 = ms.eval_curve(ent.params.curve1, [0; 1]);
    c2 = ms.eval_curve(ent.params.curve2, [0; 1]);
    side = ~strcmp(e.name, 'box_bottom') && ~strcmp(e.name, 'box_top');
    if side
        check(all(c1(:, 3) == c1(1, 3)) && all(c2(:, 3) == c2(1, 3)) && c1(1, 3) ~= c2(1, 3), ...
            '%s: the parser''s z depends on v alone', e.name);
    end
    check(e.swap_uv == side && e.exact, '%s: swap_uv and exact', e.name);
    j = find(strcmp({gb.outer.name}, e.name));
    check(isequal(e.surf, gb.outer(j).surf), '%s: transposed control net equals box.ms2', e.name);
    for b = 1:4
        nb = e.([seam_field(b)]);
        nbb = gb.outer(j).([seam_field(b)]);
        check(isequal(nb, nbb), '%s: seams differ from box.ms2', e.name);
    end
    if side
        % the parser's edge e is S1 boundary 5 - e of a swapped patch
        t = [0; 0.3; 1];
        for edge = 1:4
            switch edge
                case 1, Pp = ms.eval_surface_grid(e.name, t, 0);
                case 2, Pp = ms.eval_surface_grid(e.name, 1, t);
                case 3, Pp = ms.eval_surface_grid(e.name, t, 1);
                case 4, Pp = ms.eval_surface_grid(e.name, 0, t);
            end
            c = boundary_curve(e.surf, 5 - edge);
            Q = mwecmass.solid.eval_bspline_curve(c, t);
            % lines evaluated two ways: rounding of one product and one sum
            check(max(max(abs(reshape(Pp, [], 3) - Q))) <= 8 * eps(3), '%s: parser edge %d is not S1 boundary %d', ...
                e.name, edge, 5 - edge);
        end
        nsw = nsw + 1;
    end
end
check_seams(gs.outer, 'box_swapped');
Lb = mwecmass.solid.slice_bspline_surface(gb.outer, -1);
Ls = mwecmass.solid.slice_bspline_surface(gs.outer, -1);
check(isequal(Lb, Ls), 'box_swapped: F4 section at z = -1 differs from box.ms2');
fprintf('box_swapped: %d side faces swapped, control nets and seams equal box.ms2, F4 at z = -1 bitwise equal (area %.17g)\n', ...
    nsw, Ls.area);

% split_wall_box: F1 decisions only (column corners at z = -1)
m = parse(fullfile(fix, 'split_wall_box.ms2'));
g = mwecmass.solid.outer_nurbs(m);
P = g.outer;
check(all([P.exact]), 'split_wall_box: every patch exact');
names = {P.name};
cnt = @(nm) sum(strcmp(names, nm));
for nm = {'box_side1', 'box_side3', 'box_side4'}
    i = find(strcmp(names, nm{1}));
    check(numel(i) == 2 && isequal(P(i(1)).z_range, [-2.5 -1]) && isequal(P(i(2)).z_range, [-1 0.5]), ...
        '%s: split at z = -1 into two pieces', nm{1});
    check(isequal(P(i(1)).surf.ctrl(end, :, :), P(i(2)).surf.ctrl(1, :, :)) && all(P(i(1)).surf.ctrl(end, :, 3) == -1), ...
        '%s: cut row shared bitwise at z = -1', nm{1});
end
for nm = {'box_side2lo', 'box_side2hi', 'box_bottom', 'box_top'}
    check(cnt(nm{1}) == 1, '%s: not split', nm{1});
end
M2 = m.eval_point('M2');
M3 = m.eval_point('M3');
s1 = find(strcmp(names, 'box_side1'));
s3 = find(strcmp(names, 'box_side3'));
check(isequal(reshape(P(s1(1)).surf.ctrl(end, end, :), 1, 3), M2), 'box_side1 cut row ends at M2');
check(isequal(reshape(P(s3(1)).surf.ctrl(end, 1, :), 1, 3), M3), 'box_side3 cut row ends at M3');
lo2 = find(strcmp(names, 'box_side2lo'));
bot = find(strcmp(names, 'box_bottom'));
check(isequal(P(s1(1)).seam_v1, [lo2 1]), 'lower piece of box_side1 meets box_side2lo along its v1 column');
ca = boundary_curve(P(s1(1)).surf, 3);
cb = boundary_curve(P(lo2).surf, 1);
check(isequal(ca.knots, [0 0 0.5 0.5]) && isequal(cb.knots, [0 0 1 1]) && isequal(ca.ctrl, cb.ctrl), ...
    'the column [0 0 u* u*] and the Line [0 0 1 1] are the same curve');
check(isequal(P(lo2).seam_u0, [bot 3]), 'box_side2lo meets box_bottom');
ca = boundary_curve(P(lo2).surf, 4);
cb = boundary_curve(P(bot).surf, 3);
check(isequal(ca.ctrl, flipud(cb.ctrl)), 'box_side2lo and box_bottom share the Line run the other way');
n_seam = check_seams(P, 'split_wall_box');
fprintf('split_wall_box: %d entries, %d seam boundaries mutual and the same curve\n', numel(P), n_seam);

% stepped_spar: rev_z, split at the step (existing C0 knots), F3b and F4 around the step
m = parse(fullfile(fix, 'stepped_spar.ms2'));
g = mwecmass.solid.outer_nurbs(m);
P = g.outer;
check(numel(P) == 12 && all([P.exact]) && all(strcmp({P.offset_kind}, 'rev_z')), 'stepped_spar: 12 exact rev_z entries');
prof = [0 0 -3; 1.5 0 -3; 1.5 0 -1; 0.75 0 -1; 0.75 0 1; 0 0 1];
for q = 1:4
    i = 3 * (q - 1) + (1:3);
    check(isequal(reshape([P(i).z_range], 2, [])', [-3 -1; -1 -1; -1 1]), 'quarter %d: z ranges', q);
    check(isequal(P(i(1)).surf.ctrl(end, :, :), P(i(2)).surf.ctrl(1, :, :)) && ...
        isequal(P(i(2)).surf.ctrl(end, :, :), P(i(3)).surf.ctrl(1, :, :)), 'quarter %d: cut rows shared', q);
    fl = [any(strcmp(P(i(1)).flips, 'X')) any(strcmp(P(i(1)).flips, 'Y'))];
    rows = [1 2 3; 3 4 0; 4 5 6];
    for j = 1:3
        rr = rows(j, rows(j, :) > 0);
        exp = zeros(numel(rr), 3, 3);
        for a = 1:numel(rr)
            r = prof(rr(a), 1);
            exp(a, :, :) = reshape([r 0 prof(rr(a), 3); r r prof(rr(a), 3); 0 r prof(rr(a), 3)], 1, 3, 3);
        end
        for c = 1:2
            if fl(c)
                exp(:, :, c) = -exp(:, :, c);
            end
        end
        check(isequal(P(i(j)).surf.ctrl, exp), 'quarter %d piece %d: control points are the deck''s rows', q, j);
    end
end
parent = P(1);
parent.surf.ctrl = cat(1, P(1).surf.ctrl, P(2).surf.ctrl(2, :, :), P(3).surf.ctrl(2:end, :, :));
parent.surf.weights = [P(1).surf.weights; P(2).surf.weights(2, :); P(3).surf.weights(2:end, :)];
parent.surf.knots{1} = [0 0 0.2 0.4 0.6 0.8 1 1];
parent.z_range = [-3 1];
parent.u_range = [0 1];
expect_error(@() mwecmass.solid.split_bspline_surface(parent, -1), 'mwecmass:solid:ZNotMonotonic');
[lo, hi] = mwecmass.solid.split_bspline_surface(P(1), -2);
check(all(lo.surf.ctrl(end, :, 3) == -2) && isequal(lo.surf.ctrl(end, :, :), hi.surf.ctrl(1, :, :)), 'F3b at z = -2');
[lo, hi] = mwecmass.solid.split_bspline_surface(P(3), 0);
check(all(lo.surf.ctrl(end, :, 3) == 0) && isequal(lo.surf.ctrl(end, :, :), hi.surf.ctrl(1, :, :)), 'F3b at z = 0');
zr = reshape([P.z_range], 2, [])';
below = P(max(zr, [], 2) <= -1 & zr(:, 1) ~= zr(:, 2));
above = P(min(zr, [], 2) >= -1 & zr(:, 1) ~= zr(:, 2));
Lb = mwecmass.solid.slice_bspline_surface(below, -1);
La = mwecmass.solid.slice_bspline_surface(above, -1);
rb = sqrt(sum(Lb.pts(:, 1:2).^2, 2));
ra = sqrt(sum(La.pts(:, 1:2).^2, 2));
fprintf(['stepped_spar at z = -1: faces below radius %.17g..%.17g (area %.17g, pi 1.5^2 = %.17g), ' ...
    'faces above radius %.17g..%.17g (area %.17g, pi 0.75^2 = %.17g)\n'], min(rb), max(rb), Lb.area, ...
    pi * 2.25, min(ra), max(ra), La.area, pi * 0.5625);
% points of a rational arc of radius r: rounding of one evaluation
check(max(abs(rb - 1.5)) <= 8 * eps(1.5) && max(abs(ra - 0.75)) <= 8 * eps(1.5), 'stepped_spar: step circles');
shelf = P(2:3:end);
check(all(arrayfun(@(e) ~e.outward, shelf([1 4]))) && all(arrayfun(@(e) e.outward, shelf([2 3]))), ...
    'stepped_spar: shelf normals');
fprintf('stepped_spar: 12 entries (3 per quarter), cut rows shared, deck rows unchanged, F3b at -2 and 0\n');

% deck check, on copies written here
tmp = tempname();
mkdir(tmp);
cleanup = onCleanup(@() rmdir_quiet(tmp));
txt = fileread(fullfile(sk, 'box.ms2'));
deck = strrep(txt, 'EndModel;', sprintf('Variable vv 0 -1 / 1.5 ;\nFooCurve junk 11 -1 / * B1 T1 ;\nEndModel;'));
f1 = fullfile(tmp, 'box.ms2');
write_text(f1, deck);
g1 = mwecmass.solid.outer_nurbs(parse(f1));
check(isequal(g1.outer, gb.outer), 'unreferenced Variable and unknown lines change geo');
write_text(f1, strrep(txt, 'FramePoint T1 ', 'FooPoint T1 '));
expect_error(@() mwecmass.solid.outer_nurbs(parse(f1)), 'mwecmass:solid:UnsupportedEntity');
write_text(f1, strrep(txt, 'RuledSurf box_top ', 'FooSurf box_top '));
msg = expect_error(@() mwecmass.solid.outer_nurbs(parse(f1)), 'mwecmass:solid:HullNotClosed');
check(~isempty(strfind(msg, 'FooSurf box_top')), 'HullNotClosed lists the unparsed line'); %#ok<STREMP>
txt = fileread(fullfile(sk, 'cylinder.ms2'));
f2 = fullfile(tmp, 'cylinder.ms2');
write_text(f2, strrep(txt, 'Symmetry:  x y', 'Symmetry:  z'));
expect_error(@() mwecmass.solid.outer_nurbs(parse(f2)), 'mwecmass:solid:UnsupportedMirror');
fprintf('deck check: unreferenced lines ignored; UnsupportedEntity, HullNotClosed, UnsupportedMirror raised\n');

% the general path needs the boundary cache and t_min
expect_error(@() mwecmass.solid.outer_nurbs(mb, [], struct('force_general', true)), 'mwecmass:solid:FitInputMissing');
end

function f = seam_field(b)
fields = {'seam_v0', 'seam_u1', 'seam_v1', 'seam_u0'};
f = fields{b};
end

function n_seam = check_seams(P, name)
n_seam = 0;
for k = 1:numel(P)
    for b = 1:4
        nb = P(k).(seam_field(b));
        if isempty(nb)
            check(b == 2 && P(k).pole(2) || b == 4 && P(k).pole(1), '%s: %s boundary %d open', name, P(k).name, b);
            continue
        end
        check(isequal(P(nb(1)).(seam_field(nb(2))), [k b]), '%s: %s seam %d not mutual', name, P(k).name, b);
        check(same_curve(boundary_curve(P(k).surf, b), boundary_curve(P(nb(1)).surf, nb(2))), ...
            '%s: %s seam %d is not the same curve', name, P(k).name, b);
        n_seam = n_seam + 1;
    end
end
end

function tf = same_curve(a, b)
% contract S1: degree, bitwise control points and weights, knots equal after the affine map to [0, 1]
% (4 * 2^-52: three correctly rounded operations); either direction
wa = a.weights;
wb = b.weights;
if isempty(wa), wa = ones(size(a.ctrl, 1), 1); end
if isempty(wb), wb = ones(size(b.ctrl, 1), 1); end
ka = (a.knots - a.knots(1)) / (a.knots(end) - a.knots(1));
kb = (b.knots - b.knots(1)) / (b.knots(end) - b.knots(1));
kr = fliplr(1 - kb);
tf = a.degree == b.degree && numel(ka) == numel(kb) && ...
    ((isequal(a.ctrl, b.ctrl) && isequal(wa, wb) && all(abs(ka - kb) <= 4 * 2^-52)) || ...
    (isequal(a.ctrl, flipud(b.ctrl)) && isequal(wa, flipud(wb)) && all(abs(ka - kr) <= 4 * 2^-52)));
end

function c = boundary_curve(s, b)
W = s.weights;
if isempty(W)
    W = ones(size(s.ctrl, 1), size(s.ctrl, 2));
end
switch b
    case 1, c = struct('degree', s.degree(1), 'ctrl', reshape(s.ctrl(:, 1, :), [], 3), 'knots', s.knots{1}, 'weights', W(:, 1));
    case 2, c = struct('degree', s.degree(2), 'ctrl', reshape(s.ctrl(end, :, :), [], 3), 'knots', s.knots{2}, 'weights', W(end, :)');
    case 3, c = struct('degree', s.degree(1), 'ctrl', reshape(s.ctrl(:, end, :), [], 3), 'knots', s.knots{1}, 'weights', W(:, end));
    case 4, c = struct('degree', s.degree(2), 'ctrl', reshape(s.ctrl(1, :, :), [], 3), 'knots', s.knots{2}, 'weights', W(1, :)');
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

function rmdir_quiet(d)
if exist(d, 'dir')
    rmdir(d, 's');
end
end

function msg = expect_error(f, id)
try
    f();
catch err
    if ~strcmp(err.identifier, id)
        error('test_outer_nurbs_fixtures:fail', 'expected %s, got %s (%s)', id, err.identifier, err.message);
    end
    msg = err.message;
    return
end
error('test_outer_nurbs_fixtures:fail', 'expected %s, got no error', id);
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
    error('test_outer_nurbs_fixtures:fail', varargin{:});
end
end
