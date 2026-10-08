function test_build_body_fixtures()
%TEST_BUILD_BODY_FIXTURES  F5, F6, F6b on the SK cylinder (real F1, F2) and the SK box (real F1,
%   inner sets from sti_inner_box) in both modes: ballast inside a module, at a module edge,
%   spilled, below the inner z_lo and at it, no ballast, all solid, void ends on module edges.
%   Asserted (body_checks): validate_brep, I2, I1, I7, the STEP import; here also the region
%   volumes, first and second moments of every module against the closed forms of
%   sti_closed_form (the independent oracle, contract section 3): the box's faces are bilinear, so
%   Gauss-Legendre of order 8 is exact and the difference is rounding (asserted); the cylinder's
%   faces are rational, so its differences are printed. F6b section areas likewise.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
fix = fullfile(root, 'tests', 'standins', 'fixtures');
evalc('mc = mwecmass.geometry.MS2Parser.parse(fullfile(fix, ''cylinder.ms2''));');
evalc('mb = mwecmass.geometry.MS2Parser.parse(fullfile(fix, ''box.ms2''));');
gc = mwecmass.solid.outer_nurbs(mc);
gb = mwecmass.solid.outer_nurbs(mb);
t_min = 0.1;
opts = struct('t_min', t_min);
ec = [-3; -2; -1; 0; 1];
eb = [-2.5; -1.5; -0.5; 0.5];
ic = [mwecmass.solid.offset_surface(mc, [], gc, 0.1, gc.z_range, opts), ...
      mwecmass.solid.offset_surface(mc, [], gc, 0.2, gc.z_range, opts)];
ib = [sti_inner_box(gb, 0.1, t_min), sti_inner_box(gb, 0.15, t_min)];
d1 = ic(1).d;
d2 = ic(2).d;
db = ib(1).d;

% {geo, inner sets, mode, edges, t, z_ballast, solid_modules, label}
cases = {
    gc, ic, 'modular_precast', ec, [0.1; 0.1; 0.2; NaN], -2.5, 4, 'cylinder precast, ballast inside module 1, joint_step, wall on top'
    gc, ic, 'modular_precast', ec, [0.1; 0.2; 0.2; NaN], -2, 4, 'cylinder precast, ballast at an edge'
    gc, ic, 'modular_precast', ec, [0.1; 0.1; 0.1; NaN], -1.5, 4, 'cylinder precast, ballast spilled into module 2'
    gc, ic, 'modular_precast', ec, [0.2; 0.1; 0.1; 0.1], -2.95, [], 'cylinder precast, ballast below the inner z_lo'
    gc, ic, 'modular_precast', ec, [0.1; 0.1; 0.1; 0.1], -3 + d1, [], 'cylinder precast, ballast at the inner z_lo'
    gc, ic, 'modular_precast', [-3; -3 + d1; -1; 0; 1], [0.1; 0.1; 0.1; NaN], -3, 4, 'cylinder precast, inner z_lo on a module edge'
    gc, ic, 'modular_precast', [-3; -2; 1 - d1; 1], [0.1; 0.1; 0.1], -2.5, [], 'cylinder precast, inner z_hi on a module edge'
    gc, ic, 'modular_precast', [-3; -2; 1 - d2; 1], [0.1; 0.2; 0.1], -2.5, [], 'cylinder precast, thicker lower z_hi on a joint'
    gc, [], 'modular_precast', ec, NaN(4, 1), -2, 1:4, 'cylinder precast, every module solid'
    gc, ic, 'thin_shell', ec, 0.1 * ones(4, 1), -2.5, [], 'cylinder thin shell, ballast inside module 1'
    gc, ic, 'thin_shell', ec, 0.1 * ones(4, 1), -1, [], 'cylinder thin shell, ballast at an edge'
    gc, ic, 'thin_shell', ec, 0.1 * ones(4, 1), -2.95, [], 'cylinder thin shell, ballast below the inner z_lo'
    gc, ic, 'thin_shell', ec, 0.1 * ones(4, 1), -3 + d1, [], 'cylinder thin shell, ballast at the inner z_lo'
    gc, ic, 'thin_shell', ec, 0.2 * ones(4, 1), -3, [], 'cylinder thin shell, no ballast'
    gb, ib, 'modular_precast', eb, [0.1; 0.15; NaN], -2, 3, 'box precast, ballast inside module 1'
    gb, ib, 'modular_precast', eb, [0.15; 0.1; 0.1], -1.5, [], 'box precast, ballast at an edge'
    gb, ib, 'modular_precast', eb, [0.1; 0.1; 0.1], -2.5 + db, [], 'box precast, ballast at the inner z_lo'
    gb, ib, 'modular_precast', eb, [0.1; 0.1; 0.1], -0.8, [], 'box precast, ballast spilled into module 2'
    gb, [], 'modular_precast', eb, NaN(3, 1), -2.5, 1:3, 'box precast, every module solid'
    gb, ib, 'thin_shell', eb, 0.1 * ones(3, 1), -1, [], 'box thin shell, ballast in module 2'
    gb, ib, 'thin_shell', eb, 0.15 * ones(3, 1), -2.45, [], 'box thin shell, ballast below the inner z_lo'
    gb, ib, 'thin_shell', eb, 0.1 * ones(3, 1), -2.5 + db, [], 'box thin shell, ballast at the inner z_lo'
    gb, ib, 'thin_shell', eb, 0.1 * ones(3, 1), -1.5, [], 'box thin shell, ballast at an edge'
    };
worst = struct('V', 0, 'S', 0, 'J', 0, 'sec', 0);
for c = 1:size(cases, 1)
    [geo, sets, mode, edges, t, zb, solid, label] = cases{c, :};
    design = struct('mode', mode, 'edges', edges, 'vs', 0.5, 't', t, 'z_ballast', zb, 'solid_modules', solid);
    body = mwecmass.solid.build_body(geo, design, sets);
    precast = strcmp(mode, 'modular_precast');
    if precast
        rho = struct('uhpc', 2500, 'air', 1.2);
        names = {'uhpc', 'air'};
    else
        rho = struct('ballast', 7500, 'shell', 7850, 'air', 1.2);
        names = {'ballast', 'shell', 'air'};
    end
    st = body_checks(geo, body, rho, label);
    fx = sti_closed_form('fixture', geo);
    N = numel(edges) - 1;
    d = NaN(N, 1);
    for i = 1:N
        if isfinite(t(i)) && ~any(solid == i)
            d(i) = sets(arrayfun(@(s) isequal(s.t, t(i)), sets)).d;
        end
    end
    reg = sti_closed_form('regions', fx, design, d);
    box = strcmp(fx.kind, 'box');
    L = max(abs([geo.z_range, 1.5]));
    % per node: products of at most five coordinate-sized factors from bilinear evaluations (a few
    % ulp each) and Gauss nodes and weights (2 ulp): 100 eps per term on the scale L^k of the
    % integrand times the face area; then the recursive summation over the nodes of a face
    unit = (100 + st.n_nodes) * eps * st.area;
    dev = [0 0 0];
    for i = 1:N
        for r = 1:numel(names)
            R = st.bp.regions.(names{r});
            Q = reg.(names{r});
            dV = abs(R.V(i) - Q.V(i));
            dS = max(abs(R.S(i, :) - Q.S(i, :)));
            dJ = max(max(abs(R.J(:, :, i) - Q.J(:, :, i))));
            dev = max(dev, [dV, dS, dJ]);
            if box
                check(dV <= unit * L && dS <= unit * L^2 && dJ <= unit * L^3, ...
                    '%s: module %d region %s differs from the closed form (V %.2e, S %.2e, J %.2e)', label, i, names{r}, dV, dS, dJ);
            end
        end
    end
    worst.V = max(worst.V, dev(1));
    worst.S = max(worst.S, dev(2));
    worst.J = max(worst.J, dev(3));
    fprintf('  F6 - closed form, largest over modules and regions: V %.2e m3, S %.2e m4, J %.2e m5%s\n', dev, ...
        repmat(' (asserted)', 1, box));

    % F6b inside every module, away from its planes
    so = sti_closed_form('section', fx, 0);
    for i = 1:N
        z = edges(i) + 0.6 * (edges(i + 1) - edges(i));
        sec = mwecmass.solid.body_section(body, z);
        rs = sti_closed_form('region_sections', fx, design, d, i, z);
        check(sec.module == i && sec.solid == (rs.air.A == 0), '%s: F6b module %d', label, i);
        e1 = abs(sec.outer.area - so.A) / so.A;
        e2 = 0;
        if ~sec.solid
            e2 = abs(sec.inner.area - rs.air.A) / rs.air.A;
        end
        if box
            % Green's theorem with 16-point Gauss on straight rows: exact up to rounding
            check(e1 <= 64 * eps && e2 <= 64 * eps, '%s: F6b section areas at z = %g', label, z);
        end
        worst.sec = max([worst.sec, e1, e2]);
    end
end
fprintf('largest |F6 - closed form| over all cases: V %.2e m3, S %.2e m4, J %.2e m5 (box asserted, cylinder printed)\n', ...
    worst.V, worst.S, worst.J);
fprintf('largest relative F6b section area difference from the closed form: %.2e\n', worst.sec);

% errors: edges not spanning the hull, a hollow module without its inner set, and a joint_step
% whose loops cross (the t = 0.15 set replaced by the t = 0.1 set moved 0.2 m along x)
design = struct('mode', 'modular_precast', 'edges', [-2.4; -1.5; -0.5; 0.5], 'vs', 0, 't', [0.1; 0.1; 0.1], ...
    'z_ballast', -2.5, 'solid_modules', []);
expect_error(@() mwecmass.solid.build_body(gb, design, ib), 'mwecmass:solid:BadEdges');
design.edges = eb;
design.t = [0.1; 0.2; 0.1];
expect_error(@() mwecmass.solid.build_body(gb, design, ib), 'mwecmass:solid:MissingInnerSet');
moved = ib(1);
moved.t = 0.15;
for k = 1:numel(moved.patches)
    moved.patches(k).surf.ctrl(:, :, 1) = moved.patches(k).surf.ctrl(:, :, 1) + 0.2;
end
design.t = [0.1; 0.15; 0.15];
expect_error(@() mwecmass.solid.build_body(gb, design, [ib(1), moved]), 'mwecmass:solid:JointNotNested');
fprintf('errors BadEdges, MissingInnerSet, JointNotNested raised\n');
end

function expect_error(fn, id)
try
    fn();
catch err
    check(strcmp(err.identifier, id), 'expected %s, got %s: %s', id, err.identifier, err.message);
    return
end
error('test_build_body_fixtures:fail', 'expected error %s', id);
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
    error('test_build_body_fixtures:fail', varargin{:});
end
end
