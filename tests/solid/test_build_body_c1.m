function test_build_body_c1()
%TEST_BUILD_BODY_C1  F5, F6, F6b on the C1 hull, modular precast, t = 76.2 mm in every hollow module.
%   Module edges as build_config makes them (four platform modules below a 1.8 m solid wall
%   module, the ends at the hull's own z range). Ballast level inside module 2, at the module edge
%   2|3, spilled into module 3, below the inner z_lo in module 1, and at the inner z_lo (a pole: no
%   cut). Asserted (body_checks): validate_brep, I2, I1, I7 and the STEP import; here also I8 (CG
%   x = y = 0 and products of inertia 0 to the rounding bound of the face sums) and the F6b
%   sections. Independent references, printed: the gmsh surface-mesh volume of every imported solid
%   (body_checks; element size 0.02 m in the first case), the C1 oracle of contract section 3
%   (stadium section of area 0.4 + 0.01 pi for z in [-0.5, 1.0]), and the change of every quantity
%   with the Gauss order of F6.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
geo = c1_kernel_cache();
zr = geo.z_range;
wall_height = 1.8;
e = [linspace(zr(1), zr(2) - wall_height, 5)'; zr(2)];
t = 0.0762;
[~, inner] = c1_kernel_cache(t, [zr(1) e(5)]);
fprintf('C1 inner set t = %.4f m (d = %.17g) over [%.4f, %.4f]: z_lo %.17g\n', t, inner.d, zr(1), e(5), inner.z_lo);
fprintf('C1 module edges (body) %s\n', sprintf('%.17g ', e));
rho = struct('uhpc', 2500, 'air', 1.2);
tt = [t; t; t; t; NaN];
cases = {-2.3, 'ballast inside module 2'
         e(3), 'ballast at the module edge 2|3'
         -1.6, 'ballast spilled into module 3'
         -3.2, 'ballast in module 1 below the inner z_lo'
         inner.z_lo, 'ballast at the inner z_lo (keel pole)'
         zr(1), 'no ballast'};
for c = 1:size(cases, 1)
    zb = cases{c, 1};
    label = sprintf('C1 precast, %s (z_ballast %.6f)', cases{c, 2}, zb);
    design = struct('mode', 'modular_precast', 'edges', e, 'vs', 0.9, 't', tt, 'z_ballast', zb, 'solid_modules', 5);
    tic;
    body = mwecmass.solid.build_body(geo, design, inner);
    tb = toc;
    bopts = struct();
    if c == 1
        bopts.mesh_size = 0.02;
    end
    st = body_checks(geo, body, rho, label, bopts);
    fprintf('  build_body %.2f s\n', tb);
    bp = st.bp;
    M = bp.total.mass;
    c3 = bp.total.CG_body;
    % I8: mirrored faces give exactly negated x (or y) integrals, so CG x, y and the products of
    % inertia are rounding residues of sums over n faces (recursive summation: (n + 10) eps times
    % the sum of absolute terms; 10 covers the region, module and total sums)
    nb = (st.n_faces + 10) * eps * rho.uhpc;
    bx = nb * st.abs_moments(2) / M;
    by = nb * st.abs_moments(3) / M;
    check(abs(c3(1)) <= bx && abs(c3(2)) <= by, '%s: CG x, y = %.3e, %.3e (bounds %.1e, %.1e)', label, c3(1:2), bx, by);
    I = bp.total.I_cg;
    bxy = nb * st.abs_moments(8) + M * bx * by;
    bxz = nb * st.abs_moments(9) + M * bx * abs(c3(3));
    byz = nb * st.abs_moments(10) + M * by * abs(c3(3));
    check(abs(I(1, 2)) <= bxy && abs(I(1, 3)) <= bxz && abs(I(2, 3)) <= byz, ...
        '%s: products of inertia %.3e %.3e %.3e', label, I(1, 2), I(1, 3), I(2, 3));
    fprintf('  I8: CG x %.2e (bound %.1e), y %.2e (bound %.1e); Ixy %.2e, Ixz %.2e, Iyz %.2e kg m2 (bounds %.1e %.1e %.1e)\n', ...
        c3(1), bx, c3(2), by, I(1, 2), I(1, 3), I(2, 3), bxy, bxz, byz);
    k_star = find(e(1:end - 1) <= zb & zb < e(2:end), 1);
    for i = 1:4
        if ~(zb >= e(i + 1)) && i >= k_star
            check(body.voids(i).z_hi > body.voids(i).z_lo, '%s: module %d has no void', label, i);
        end
    end
    if c == 1
        body_a = body;
    end
end

% F6b: sections of the first case at module mid-heights, against F4 of the written outer patches
fprintf('F6b sections (ballast inside module 2):\n');
for i = 1:5
    z = (e(i) + e(i + 1)) / 2;
    sec = mwecmass.solid.body_section(body_a, z);
    L = mwecmass.solid.slice_bspline_surface(geo.outer, z);
    check(sec.module == i && abs(sec.outer.area - L.area) <= 64 * eps * L.area, 'F6b outer section at z = %g', z);
    check(sec.solid == ~(i >= 2 && i <= 4 && z > -2.3), 'F6b solid flag at z = %g', z);
    ai = 0;
    if ~sec.solid
        ai = sec.inner.area;
        check(sec.inner.simple && ai < sec.outer.area, 'F6b void section at z = %g', z);
    end
    fprintf('  z = %8.4f (module %d): outer %.12f m2, void %.12f m2\n', z, sec.module, sec.outer.area, ai);
end
sec = mwecmass.solid.body_section(body_a, 0.2);
fprintf('  z = 0.2: outer area %.15f m2, C1 oracle 0.4 + 0.01 pi = %.15f m2 (difference %.2e)\n', ...
    sec.outer.area, 0.4 + 0.01 * pi, sec.outer.area - 0.4 - 0.01 * pi);
s1 = mwecmass.solid.body_section(body_a, e(4));
s2 = mwecmass.solid.body_section(body_a, e(4), 'below');
check(s1.module == 4 && s2.module == 3 && ~isempty(s1.inner) && ~isempty(s2.inner), 'F6b sides at a module edge');
s1 = mwecmass.solid.body_section(body_a, -2.3);
s2 = mwecmass.solid.body_section(body_a, -2.3, 'below');
check(~s1.solid && s2.solid && s1.module == 2 && s2.module == 2, 'F6b sides at z_ballast');

% convergence with the Gauss order of F6
fprintf('F6 Gauss order (ballast inside module 2): mass [kg], CG_z [m], Iyy [kg m2], V_air module 3 [m3]\n');
for ng = [2 4 6 8 10 12 16]
    b = mwecmass.solid.body_properties(body_a, rho, struct('n_gauss', ng));
    fprintf('  n = %2d: %.10f %.13f %.9f %.12f\n', ng, b.total.mass, b.total.CG_body(3), b.total.I_cg(2, 2), b.modules(3).V_air);
end
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests'));
addpath(fullfile(root, 'tests', 'step'));
addpath(fullfile(root, 'tests', 'solid'));
end

function check(cond, varargin)
if ~cond
    error('test_build_body_c1:fail', varargin{:});
end
end
