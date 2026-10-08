function test_build_body_flat_regions()
%TEST_BUILD_BODY_FLAT_REGIONS  F5, F6b, F7 with flat regions (S1b and S2 'flat', general path) in
%   place of constant-z patches. No deck before J1 has a flat region, so this test makes them from
%   all-exact fixtures: stepped_box with its low deck (z = 0.5, normal +z) and its bottom (z = -2.5,
%   normal -z) as flat regions, and the SK box's inner set (sti_inner_box) with its bottom and top
%   as inner flat regions; each row that bounded a removed patch names the region (seam [0 j]).
%   The same body with the constant-z patches is the reference: the lateral faces are the same, the
%   plane faces add 0 to F6, so the region volumes agree to rounding (asserted with the bound of
%   body_checks). Asserted besides body_checks (validate_brep, I2, I1, I7, STEP import): one plane
%   face per flat region with the expected loop count and area, F6b at the flat heights, and F7's
%   S_wet with flat regions equal to the patch version.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
evalc('model = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''tests'', ''solid'', ''fixtures'', ''stepped_box.ms2''));');
geo = mwecmass.solid.outer_nurbs(model);
names = {geo.outer.name};
deck = find(strncmp(names, 'deck_lo', 7));
bottom = find(strncmp(names, 'bottom_', 7));
geoF = to_flat(geo, {deck, bottom}, [0.5, -2.5], [1, -1]);
check(numel(geoF.outer) == 16 && numel(geoF.flat) == 2, 'stepped_box flat version: 16 patches, 2 flat regions');
rho = struct('uhpc', 2500, 'air', 1.2);
for e = {[-2.5; 0.5; 1.5], [-2.5; -1; 1.5]}
    design = struct('mode', 'modular_precast', 'edges', e{1}, 'vs', 0, 't', [NaN; NaN], 'z_ballast', -2.5, 'solid_modules', [1 2]);
    label = sprintf('stepped_box with flat regions, edge at z = %g', e{1}(2));
    body = mwecmass.solid.build_body(geoF, design, []);
    st = body_checks(geoF, body, rho, label);
    ref = mwecmass.solid.body_properties(mwecmass.solid.build_body(geo, design, []), rho, []);
    for i = 1:2
        check(abs(st.bp.modules(i).V - ref.modules(i).V) <= st.bound_module(i), '%s: module %d volume %.17g, patch version %.17g', ...
            label, i, st.bp.modules(i).V, ref.modules(i).V);
    end
    B = body.brep;
    pf = find(arrayfun(@(f) strcmp(B.surfaces{f.surface}.type, 'plane') && strcmp(f.role, 'outer'), B.faces));
    check(numel(pf) == 2, '%s: %d outer plane faces', label, numel(pf));
    for k = pf
        s = B.surfaces{B.faces(k).surface};
        A = face_area(B, k);
        if s.origin(3) == 0.5
            check(s.normal(3) == 1 && numel(B.faces(k).loops) == 1 && numel(B.faces(k).loops{1}) == 6 && abs(A - 1.5) <= 64 * eps, ...
                '%s: low deck face (area %.17g)', label, A);
        else
            check(s.origin(3) == -2.5 && s.normal(3) == -1 && numel(B.faces(k).loops{1}) == 8 && abs(A - 3) <= 64 * eps * 3, ...
                '%s: bottom face (area %.17g)', label, A);
        end
    end
    s1 = mwecmass.solid.body_section(body, -2.5);
    s2 = mwecmass.solid.body_section(body, 0.5, 'below');
    s3 = mwecmass.solid.body_section(body, 0.5);
    check(abs(s1.outer.area - 3) <= 64 * eps * 3 && abs(s2.outer.area - 3) <= 64 * eps * 3 && abs(s3.outer.area - 1.5) <= 64 * eps, ...
        '%s: F6b at the flat heights', label);
    fprintf('  plane faces of the flat regions: areas %s m2; F6b at z = -2.5: %.17g, at 0.5 below %.17g, above %.17g\n', ...
        sprintf('%.17g ', arrayfun(@(k) face_area(B, k), pf)), s1.outer.area, s2.outer.area, s3.outer.area);
end
design = struct('mode', 'thin_shell', 'edges', [-2.5; 0.5; 1.5], 'vs', 0, 't', [NaN; NaN], 'z_ballast', 0.5, 'solid_modules', []);
body = mwecmass.solid.build_body(geoF, design, []);
st = body_checks(geoF, body, struct('ballast', 7500, 'shell', 7850, 'air', 1.2), 'stepped_box with flat regions, thin shell, ballast up to the step');
check(abs(sum(st.bp.regions.ballast.V) - 9) <= sum(st.bound_module), 'thin shell: V_ballast %.17g (9)', sum(st.bp.regions.ballast.V));
for vs = [-1, 1, -2]
    h1 = mwecmass.solid.hydrostatics_at_draft(geoF, vs);
    h0 = mwecmass.solid.hydrostatics_at_draft(geo, vs);
    check(abs(h1.S_wet - h0.S_wet) <= 64 * eps * h0.S_wet && h1.V_sub == h0.V_sub, 'F7 vs = %g: S_wet %.17g with flat regions, %.17g with patches', ...
        vs, h1.S_wet, h0.S_wet);
    fprintf('F7 vs = %g: S_wet %.17g m2 with flat regions, %.17g m2 with constant-z patches\n', vs, h1.S_wet, h0.S_wet);
end

% inner flat regions: the SK box's inner set with its bottom and top as flat regions
evalc('mb = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''tests'', ''standins'', ''fixtures'', ''box.ms2''));');
gb = mwecmass.solid.outer_nurbs(mb);
ib = sti_inner_box(gb, 0.1, 0.1);
d = ib.d;
nb = {ib.patches.name};
ibF = ib;
pp = to_flat(struct('outer', ib.patches, 'flat', struct('z', {}, 'normal_z', {}, 'visible', {})), ...
    {find(strcmp(nb, 'box_bottom_inner')), find(strcmp(nb, 'box_top_inner'))}, [-2.5 + d, 0.5 - d], [1, -1]);
ibF.patches = pp.outer;
ibF.flat = pp.flat;
eb = [-2.5; -1.5; -0.5; 0.5];
cases = {'modular_precast', -2.5, 'box precast, void from the inner bottom to the inner top'
         'thin_shell', -2.45, 'box thin shell, ballast below the inner z_lo'
         'modular_precast', -2.5 + d, 'box precast, ballast at the inner z_lo (flat region replaced by ballast_top)'};
for c = 1:size(cases, 1)
    [mode, zb, label] = cases{c, :};
    design = struct('mode', mode, 'edges', eb, 'vs', 0, 't', 0.1 * ones(3, 1), 'z_ballast', zb, 'solid_modules', []);
    if strcmp(mode, 'modular_precast')
        r = rho;
        rn = {'uhpc', 'air'};
    else
        r = struct('ballast', 7500, 'shell', 7850, 'air', 1.2);
        rn = {'ballast', 'shell', 'air'};
    end
    body = mwecmass.solid.build_body(gb, design, ibF);
    st = body_checks(gb, body, r, label);
    ref = mwecmass.solid.body_properties(mwecmass.solid.build_body(gb, design, ib), r, []);
    for i = 1:3
        for q = 1:numel(rn)
            check(abs(st.bp.modules(i).(['V_' rn{q}]) - ref.modules(i).(['V_' rn{q}])) <= st.bound_module(i), ...
                '%s: module %d V_%s differs from the patch version', label, i, rn{q});
        end
    end
    B = body.brep;
    pf = find(arrayfun(@(f) strcmp(B.surfaces{f.surface}.type, 'plane') && strcmp(f.role, 'inner'), B.faces));
    n_expect = 2 - (zb == -2.5 + d);
    check(numel(pf) == n_expect, '%s: %d inner plane faces, expected %d', label, numel(pf), n_expect);
    for k = pf
        s = B.surfaces{B.faces(k).surface};
        check(any(s.origin(3) == [-2.5 + d, 0.5 - d]) && abs(face_area(B, k) - (2 - 2 * d) * (1.5 - 2 * d)) <= 64 * eps * 3, ...
            '%s: inner flat face at z = %.17g', label, s.origin(3));
    end
end
end

function g = to_flat(geo, groups, z, nz)
% geo with the patches of each group replaced by one flat region; rows that named them name it
keep = true(1, numel(geo.outer));
for j = 1:numel(groups)
    keep(groups{j}) = false;
end
newidx = cumsum(keep);
g = geo;
g.flat = struct('z', num2cell(z), 'normal_z', num2cell(nz), 'visible', groups);
P = geo.outer(keep);
f = {'seam_u0', 'seam_u1', 'seam_v0', 'seam_v1'};
for k = 1:numel(P)
    for q = 1:4
        sm = P(k).(f{q});
        if numel(sm) ~= 2 || sm(1) == 0
            continue
        end
        j = find(cellfun(@(G) any(G == sm(1)), groups), 1);
        if ~isempty(j)
            P(k).(f{q}) = [0 j];
        else
            P(k).(f{q}) = [newidx(sm(1)) sm(2)];
        end
    end
end
g.outer = P;
end

function A = face_area(B, k)
% area of a plane face from its loops by Green's theorem (16-point Gauss per knot span), positive
% for loops counter-clockwise seen from the face normal
s = B.surfaces{B.faces(k).surface};
A = 0;
for L = B.faces(k).loops
    for e = L{1}
        c = B.curves(B.edges(abs(e)).curve);
        A = A + sign(e) * green(c);
    end
end
A = A * sign(s.normal(3));
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
addpath(fullfile(root, 'tests', 'standins'), '-end');
addpath(fullfile(root, 'tests', 'standins', 'fixtures'), '-end');
end

function check(cond, varargin)
if ~cond
    error('test_build_body_flat_regions:fail', varargin{:});
end
end
