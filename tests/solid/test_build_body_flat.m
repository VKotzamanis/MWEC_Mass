function test_build_body_flat()
%TEST_BUILD_BODY_FLAT  F5, F6, F6b with a module edge at the height of a flat part of the outer
%   surface (plan T3; contract F5, F6b): stepped_spar.ms2 (T2a; revolved profile with a shelf at
%   z = -1) and stepped_box.ms2 (a box whose deck steps from z = 0.5 to 1.5 over x > 0).
%   Asserted: body_checks (validate_brep, I2, I1, I7, STEP import); the plane faces at the edge (a
%   cap over the area of the loop from above, less the void where both modules are hollow) by their
%   loops and areas against F4 of the faces on each side; the flat part written once, as outer
%   faces of the module below; F6b on both sides of the edge equal to F4 of that side; for
%   stepped_box (bilinear faces, so F6 is exact up to rounding) the module volumes 9 and 1.5 m3.
%   Printed: stepped_spar's module volumes next to 4.5 pi and 1.125 pi m3 (rational faces).

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
fix = fullfile(root, 'tests', 'solid', 'fixtures');

% ---------------------------------------------------------------- stepped_spar
evalc('model = mwecmass.geometry.MS2Parser.parse(fullfile(fix, ''stepped_spar.ms2''));');
geo = mwecmass.solid.outer_nurbs(model);
zs = -1;
e = [geo.z_range(1); zs; geo.z_range(2)];
t_min = 0.1;
inner = mwecmass.solid.offset_surface(model, [], geo, t_min, geo.z_range, struct('t_min', t_min));
below = arrayfun(@(p) p.z_range(1) ~= p.z_range(2) && max(p.z_range) == zs, geo.outer);
above = arrayfun(@(p) p.z_range(1) ~= p.z_range(2) && min(p.z_range) == zs, geo.outer);
Lb = mwecmass.solid.slice_bspline_surface(geo.outer(below), zs);
La = mwecmass.solid.slice_bspline_surface(geo.outer(above), zs);
fprintf('stepped_spar at z = -1: F4 from below %.15f m2 (pi 1.5^2 = %.15f), from above %.15f m2 (pi 0.75^2 = %.15f)\n', ...
    Lb.area, pi * 1.5^2, La.area, pi * 0.75^2);
cases = {'modular_precast', [NaN; NaN], [1 2], e(1), 'stepped_spar precast, both modules solid'
         'modular_precast', [t_min; t_min], [], e(1), 'stepped_spar precast, both modules hollow at t_min'
         'thin_shell', [t_min; t_min], [], zs, 'stepped_spar thin shell, ballast up to the shelf'};
for c = 1:size(cases, 1)
    [mode, t, solid, zb, label] = cases{c, :};
    design = struct('mode', mode, 'edges', e, 'vs', 0, 't', t, 'z_ballast', zb, 'solid_modules', solid);
    body = mwecmass.solid.build_body(geo, design, inner);
    if strcmp(mode, 'modular_precast')
        rho = struct('uhpc', 2500, 'air', 1.2);
    else
        rho = struct('ballast', 7500, 'shell', 7850, 'air', 1.2);
    end
    st = body_checks(geo, body, rho, label);
    B = body.brep;
    flat = flat_faces(B, zs);
    shelf = flat(strcmp({B.faces(flat).role}, 'outer'));
    check(numel(shelf) == 4 && all(arrayfun(@(k) isequal(B.faces(k).module, 1), shelf)), ...
        '%s: the shelf is not written once, as outer faces of module 1', label);
    planes = flat(~strcmp({B.faces(flat).role}, 'outer'));
    hollow = all(isfinite(t)) && isempty(solid);
    if strcmp(mode, 'modular_precast')
        check(numel(planes) == 1 && strcmp(B.faces(planes).role, 'cap'), '%s: one cap at the edge', label);
        A = face_area(B, planes);
        if hollow
            sec = mwecmass.solid.body_section(body, zs);
            check(numel(B.faces(planes).loops) == 2, '%s: the cap is an annulus', label);
            ref = La.area - sec.inner.area;
        else
            ref = La.area;
        end
        % both from 16-point Gauss on the same rows: rounding of the two sums
        check(abs(A - ref) <= 64 * eps * La.area, '%s: cap area %.17g, loop from above %.17g', label, A, ref);
        fprintf('  cap at z = -1: %.15f m2 (F4 from above%s %.15f)\n', A, repmat(' less the void', 1, hollow), ref);
    else
        check(numel(planes) == 2 && all(strcmp({B.faces(planes).role}, 'ballast_top')), '%s: ballast_top annulus and disk', label);
        A = face_area(B, planes(1)) + face_area(B, planes(2));
        check(abs(A - La.area) <= 64 * eps * La.area, '%s: ballast_top area %.17g', label, A);
    end
    fprintf('  module volumes %.12f, %.12f m3 (4.5 pi = %.12f, 1.125 pi = %.12f)\n', ...
        st.bp.modules(1).V, st.bp.modules(2).V, 4.5 * pi, 1.125 * pi);
    s1 = mwecmass.solid.body_section(body, zs);
    s2 = mwecmass.solid.body_section(body, zs, 'below');
    check(s1.module == 2 && abs(s1.outer.area - La.area) <= 64 * eps * La.area, '%s: F6b above the shelf', label);
    check(s2.module == 1 && abs(s2.outer.area - Lb.area) <= 64 * eps * Lb.area, '%s: F6b below the shelf', label);
    check(isempty(s1.inner) == ~hollow, '%s: F6b void above the shelf', label);
end

% ---------------------------------------------------------------- stepped_box
evalc('model = mwecmass.geometry.MS2Parser.parse(fullfile(fix, ''stepped_box.ms2''));');
geo = mwecmass.solid.outer_nurbs(model);
check(numel(geo.outer) == 22 && all([geo.outer.exact]) && isempty(geo.flat), 'stepped_box: 22 exact patches');
zs = 0.5;
e = [geo.z_range(1); zs; geo.z_range(2)];
check(isequal(geo.z_range, [-2.5 1.5]), 'stepped_box: z range');
design = struct('mode', 'modular_precast', 'edges', e, 'vs', 0, 't', [NaN; NaN], 'z_ballast', e(1), 'solid_modules', [1 2]);
body = mwecmass.solid.build_body(geo, design, []);
st = body_checks(geo, body, struct('uhpc', 2500, 'air', 1.2), 'stepped_box precast, both modules solid');
B = body.brep;
V = [st.bp.modules.V];
% bilinear faces, Gauss order 8: exact up to the rounding of the face sums (body_checks bound)
check(abs(V(1) - 9) <= st.bound_module(1) && abs(V(2) - 1.5) <= st.bound_module(2), ...
    'stepped_box: module volumes %.17g, %.17g (9, 1.5)', V);
flat = flat_faces(B, zs);
deck = flat(strcmp({B.faces(flat).role}, 'outer'));
check(numel(deck) == 2 && all(arrayfun(@(k) isequal(B.faces(k).module, 1), deck)), ...
    'stepped_box: the low deck is not written once, as outer faces of module 1');
cap = flat(strcmp({B.faces(flat).role}, 'cap'));
check(numel(cap) == 1 && numel(flat) == 3, 'stepped_box: one cap at the step');
above = arrayfun(@(p) p.z_range(1) ~= p.z_range(2) && min(p.z_range) == zs, geo.outer);
La = mwecmass.solid.slice_bspline_surface(geo.outer(above), zs);
A = face_area(B, cap);
check(abs(A - 1.5) <= 64 * eps && abs(La.area - 1.5) <= 64 * eps, 'stepped_box: cap area %.17g, F4 from above %.17g', A, La.area);
% the cap's loop: the bottom rows of the faces above (riser, wall y = +-0.75 and wall x = 1 above
% the step), each one edge, which the face below takes as its top row (the walls) or the low deck
% as its boundary (the riser's rows)
E = abs(B.faces(cap).loops{1});
check(numel(B.faces(cap).loops) == 1 && numel(E) == 6, 'stepped_box: the cap loop has %d edges', numel(E));
for k = E
    users = find(arrayfun(@(f) any(abs([B.faces(f).loops{:}]) == k), 1:numel(B.faces)));
    zr = arrayfun(@(f) face_zrange(B, f), users, 'UniformOutput', false);
    check(numel(users) == 3 && any(cellfun(@(r) r(1) == zs && r(2) > zs, zr)) && sum(cellfun(@(r) r(2) == zs, zr)) == 2, ...
        'stepped_box: cap edge %d is not one row shared by a face above and a face below or in the plane', k);
end
s1 = mwecmass.solid.body_section(body, zs);
s2 = mwecmass.solid.body_section(body, zs, 'below');
check(s1.module == 2 && abs(s1.outer.area - 1.5) <= 64 * eps && s2.module == 1 && abs(s2.outer.area - 3) <= 64 * eps * 3, ...
    'stepped_box: F6b at the step: %.17g above, %.17g below', s1.outer.area, s2.outer.area);
fprintf('stepped_box: module volumes %.17g, %.17g m3; cap %.17g m2; F6b at z = 0.5: %.17g (above), %.17g (below) m2\n', ...
    V, A, s1.outer.area, s2.outer.area);
end

function idx = flat_faces(B, z)
% faces lying in the plane z (plane faces and constant-z B-spline faces)
idx = zeros(1, 0);
for k = 1:numel(B.faces)
    r = face_zrange(B, k);
    if r(1) == z && r(2) == z
        idx(end + 1) = k; %#ok<AGROW>
    end
end
end

function r = face_zrange(B, k)
s = B.surfaces{B.faces(k).surface};
if strcmp(s.type, 'plane')
    r = s.origin(3) * [1 1];
else
    zz = s.ctrl(:, :, 3);
    r = [min(zz(:)), max(zz(:))];
end
end

function A = face_area(B, k)
% area of a plane face z = const from its loops by Green's theorem (16-point Gauss per knot span);
% loops are counter-clockwise (outer) and clockwise (holes) seen from the face normal
s = B.surfaces{B.faces(k).surface};
A = 0;
for L = B.faces(k).loops
    for e = L{1}
        c = B.curves(B.edges(abs(e)).curve);
        A = A + sign(e) * green(c);
    end
end
A = A * sign(s.normal(3)) * (2 * B.faces(k).same_sense - 1);
end

function a = green(c)
n = 16;
b = (1:n - 1) ./ sqrt(4 * (1:n - 1).^2 - 1);
[V, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D));
w = 2 * V(1, i)'.^2;
ku = unique(c.knots);
a = 0;
for j = 1:numel(ku) - 1
    s = (ku(j) + ku(j + 1)) / 2 + (ku(j + 1) - ku(j)) / 2 * x;
    [C, Cs] = mwecmass.solid.eval_bspline_curve(c, s);
    a = a + (ku(j + 1) - ku(j)) / 2 * (w' * (C(:, 1) .* Cs(:, 2)));
end
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'step'));
addpath(fullfile(root, 'tests', 'solid'));
end

function check(cond, varargin)
if ~cond
    error('test_build_body_flat:fail', varargin{:});
end
end
