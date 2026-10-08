function test_build_body_c1_thin()
%TEST_BUILD_BODY_C1_THIN  F5, F6, F6b on the C1 hull, thin shell, t = t_min = 25.4 mm over the whole
%   hull (inner set up to the neck top). Module edges as Stage 2 of the thin-shell mode makes them
%   (density nodes at five equally spaced heights, edges at their midpoints, half-height end
%   modules). Ballast level inside module 1, at the module edge 1|2 and in module 3.
%   Asserted (body_checks): validate_brep, I2, I1, I7 and the STEP import (the ballast solid, the
%   shell sheet from z_ballast up, and the shell layer and void as closed solids); here also I8 and
%   the bodies of contract S4 (a ballast solid only below z_ballast, a sheet of outer faces above it
%   whose boundary is the outer section at z_ballast). Printed: the mesh volumes (body_checks; element
%   size 0.03 m in the first case).

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
evalc('model = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''Input'', ''C1.ms2''));');
geo = mwecmass.solid.outer_nurbs(model);
zr = geo.z_range;
nodes = linspace(zr(1), zr(2), 5)';
e = [zr(1); (nodes(1:end - 1) + nodes(2:end)) / 2; zr(2)];
t = 0.0254;
tic;
inner = mwecmass.solid.offset_surface(model, [], geo, t, zr, struct('t_min', t));
fprintf('C1 inner set t = %.4f m (d = %.17g) over the hull: %.1f s, z_range %s\n', t, inner.d, toc, mat2str(inner.z_range, 17));
fprintf('C1 thin-shell module edges (body) %s\n', sprintf('%.17g ', e));
rho = struct('ballast', 7500, 'shell', 7500, 'air', 1.2);
cases = {-2.9, 'ballast inside module 1'
         e(2), 'ballast at the module edge 1|2'
         -1.2, 'ballast up into module 3'};
for c = 1:size(cases, 1)
    zb = cases{c, 1};
    label = sprintf('C1 thin shell, %s (z_ballast %.6f)', cases{c, 2}, zb);
    design = struct('mode', 'thin_shell', 'edges', e, 'vs', 0.9, 't', t * ones(5, 1), 'z_ballast', zb, 'solid_modules', []);
    body = mwecmass.solid.build_body(geo, design, inner);
    bopts = struct();
    if c == 1
        bopts.mesh_size = 0.03;
    end
    st = body_checks(geo, body, rho, label, bopts);
    bp = st.bp;
    B = body.brep;
    check(numel(B.bodies) == 2 && strcmp(B.bodies(1).name, 'C1_STEEL_ballast') && strcmp(B.bodies(2).name, 'C1_STEEL_shell'), ...
        '%s: bodies', label);
    % the sheet: outer faces at or above z_ballast; its free boundary is the outer section at z_ballast
    sheet = B.bodies(2).shells{1};
    zlo = arrayfun(@(f) min(reshape(B.surfaces{B.faces(abs(f)).surface}.ctrl(:, :, 3), [], 1)), sheet);
    check(all(strcmp({B.faces(abs(sheet)).role}, 'outer')) && all(zlo >= zb), '%s: sheet faces', label);
    E = abs(cell2mat(arrayfun(@(f) [B.faces(abs(f)).loops{:}], sheet, 'UniformOutput', false)));
    free = unique(E(arrayfun(@(k) sum(E == k), E) == 1));
    check(all(B.vertices(unique([B.edges(free).vertices]), 3) == zb), '%s: the sheet boundary is not at z_ballast', label);
    ballast = B.bodies(1).shells{1};
    zhi = arrayfun(@(f) face_zmax(B, abs(f)), ballast);
    check(all(zhi <= zb), '%s: ballast solid above z_ballast', label);
    V_ballast = sum(bp.regions.ballast.V);
    fprintf('  V_ballast %.12f m3, V_shell %.12f m3, V_air %.12f m3\n', V_ballast, sum(bp.regions.shell.V), sum(bp.regions.air.V));
    M = bp.total.mass;
    c3 = bp.total.CG_body;
    nb = (st.n_faces + 10) * eps * max([rho.ballast, rho.shell]);
    bx = nb * st.abs_moments(2) / M;
    by = nb * st.abs_moments(3) / M;
    check(abs(c3(1)) <= bx && abs(c3(2)) <= by, '%s: CG x, y = %.3e, %.3e (bounds %.1e, %.1e)', label, c3(1:2), bx, by);
    I = bp.total.I_cg;
    check(abs(I(1, 2)) <= nb * st.abs_moments(8) + M * bx * by && abs(I(1, 3)) <= nb * st.abs_moments(9) + M * bx * abs(c3(3)) && ...
        abs(I(2, 3)) <= nb * st.abs_moments(10) + M * by * abs(c3(3)), '%s: products of inertia', label);
    fprintf('  I8: CG x %.2e (bound %.1e), y %.2e (bound %.1e); Ixy %.2e, Ixz %.2e, Iyz %.2e kg m2\n', ...
        c3(1), bx, c3(2), by, I(1, 2), I(1, 3), I(2, 3));
    s1 = mwecmass.solid.body_section(body, zb);
    s2 = mwecmass.solid.body_section(body, zb, 'below');
    check(~s1.solid && s2.solid, '%s: F6b at z_ballast', label);
end
end

function z = face_zmax(B, k)
s = B.surfaces{B.faces(k).surface};
if strcmp(s.type, 'plane')
    z = s.origin(3);
else
    z = max(reshape(s.ctrl(:, :, 3), [], 1));
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
    error('test_build_body_c1_thin:fail', varargin{:});
end
end
