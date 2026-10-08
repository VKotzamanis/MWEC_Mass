function test_hydrostatics_at_draft()
%TEST_HYDROSTATICS_AT_DRAFT  F7 hydrostatics_at_draft (contract S7) on the SK box and cylinder, on
%   stepped_box and on C1.
%   Independent references: sti_closed_form('hydrostatics') for the SK fixtures (contract section
%   3; the box's faces are bilinear, so F6's Gauss rule and F4's Green integrals are exact and the
%   differences are rounding, asserted; the cylinder's are rational, printed); stepped_box's
%   waterplane above the step is the rectangle [0, 1] x [-0.75, 0.75] off x = 0, whose moments
%   about its centre of flotation are closed forms (I12, asserted); the C1 stadium section of
%   contract section 3 for a waterline in z in [-0.5, 1.0] (printed). Also asserted: I6 (CB = CB_body
%   + [0 0 vs], KM = CB(3) + BM_L, draft = -(z_min + vs), bitwise), the full and none cases of S7,
%   the waterline at a flat part's height (the faces below), C1 xF = yF = 0 (CB x, y printed; I8 of
%   the F6 sums is asserted in test_build_body_c1).

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
fix = fullfile(root, 'tests', 'standins', 'fixtures');
for name = {'box', 'cylinder'}
    evalc('model = mwecmass.geometry.MS2Parser.parse(fullfile(fix, [name{1} ''.ms2'']));');
    geo = mwecmass.solid.outer_nurbs(model);
    fx = sti_closed_form('fixture', geo);
    box = strcmp(name{1}, 'box');
    L = max(abs([geo.z_range, 1.5]));
    zr = geo.z_range;
    for vs = [-zr(1) + 0.5, -(zr(1) + zr(2)) / 2, -zr(2) + 0.25, -zr(2), -zr(2) - 0.5, -zr(1), -zr(1) + 0.5]
        hs = mwecmass.solid.hydrostatics_at_draft(geo, vs);
        ref = sti_closed_form('hydrostatics', fx, vs);
        frames(hs, geo, name{1});
        check(strcmp(hs.submersion, ref.submersion), '%s vs = %g: submersion %s, closed form %s', name{1}, vs, hs.submersion, ref.submersion);
        f = {'V_sub', 'S_wet', 'Aw', 'I_wp_xx', 'I_wp_yy', 'I_wp_yy_origin', 'BM_L', 'KM', 'xF', 'yF'};
        k = [1 0 0 2 2 2 -1 1 1 1];
        dev = zeros(1, numel(f));
        for j = 1:numel(f)
            a = hs.(f{j});
            b = ref.(f{j});
            check(isnan(a) == isnan(b), '%s vs = %g: %s is %g, closed form %g', name{1}, vs, f{j}, a, b);
            if ~isnan(a)
                dev(j) = abs(a - b);
            end
        end
        dcb = 0;
        if ~strcmp(hs.submersion, 'none')
            dcb = max(abs(hs.CB - ref.CB));
        end
        if box
            % per Gauss node a product of at most four coordinate-sized factors of bilinear
            % evaluations (100 eps per term on the scale L^k times the face area) and the summation
            % of 64 nodes per span; BM_L divides by V_sub, KM adds CB(3)
            unit = (100 + 64) * eps * max(hs.S_wet, 1);
            tol = unit * L.^max(k, 0);
            tol(7) = unit * L^5 / max(hs.V_sub, eps) + unit * L^2;
            tol(8) = unit * L^5;
            check(all(dev <= tol) && dcb <= unit * L^2, '%s vs = %g: differs from the closed form (%s)', name{1}, vs, ...
                sprintf('%.2e ', [dev, dcb]));
        end
        fprintf('%s vs = %6.3f (%s): V_sub %.12f, S_wet %.12f, Aw %.12f, KM %.12f; |F7 - closed form| max over S7 fields %.2e%s\n', ...
            name{1}, vs, hs.submersion, hs.V_sub, hs.S_wet, hs.Aw, hs.KM, max([dev, dcb]), repmat(' (asserted)', 1, box));
    end
end

% stepped_box: waterplane off x = 0, and the waterline at the step
evalc('model = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''tests'', ''solid'', ''fixtures'', ''stepped_box.ms2''));');
geo = mwecmass.solid.outer_nurbs(model);
vs = -1;
hs = mwecmass.solid.hydrostatics_at_draft(geo, vs);
frames(hs, geo, 'stepped_box');
V = 9 + 1.5 * 0.5;
ref = struct('V_sub', V, 'Aw', 1.5, 'xF', 0.5, 'yF', 0, 'I_wp_xx', 1 * 1.5^3 / 12, 'I_wp_yy', 1.5 * 1 / 12, ...
    'I_wp_yy_origin', 1.5 / 3, 'BM_L', (1.5 / 12) / V);
cbz = (9 * -1 + 0.75 * 0.75) / V;
cbx = 0.75 * 0.5 / V;
f = fieldnames(ref);
unit = (100 + 64) * eps * hs.S_wet * 2.5^3;
for j = 1:numel(f)
    check(abs(hs.(f{j}) - ref.(f{j})) <= unit, 'stepped_box vs = -1: %s %.17g, closed form %.17g', f{j}, hs.(f{j}), ref.(f{j}));
end
check(abs(hs.CB_body(1) - cbx) <= unit && abs(hs.CB_body(3) - cbz) <= unit && abs(hs.CB_body(2)) <= unit, 'stepped_box: CB');
fprintf(['stepped_box vs = -1: Aw %.17g, xF %.17g, I_wp_yy %.17g about the CF (closed form 0.125), I_wp_yy_origin %.17g ' ...
    '(0.5), BM_L %.17g, V_sub %.17g (9.75)\n'], hs.Aw, hs.xF, hs.I_wp_yy, hs.I_wp_yy_origin, hs.BM_L, hs.V_sub);
hs = mwecmass.solid.hydrostatics_at_draft(geo, -0.5);
check(strcmp(hs.submersion, 'partial') && abs(hs.Aw - 3) <= 64 * eps * 3 && abs(hs.V_sub - 9) <= unit && ...
    abs(hs.xF) <= 64 * eps, 'stepped_box: waterline at the step uses the faces below (Aw %.17g, V_sub %.17g)', hs.Aw, hs.V_sub);
fprintf('stepped_box waterline at the step z = 0.5: Aw %.17g (faces below: 3), V_sub %.17g (9)\n', hs.Aw, hs.V_sub);

% C1
evalc('model = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''Input'', ''C1.ms2''));');
geo = mwecmass.solid.outer_nurbs(model);
r = 0.1;
Aw = 0.4 + pi * r^2;
Iyy = 2 * 0.2^3 / 12 + pi * r^4 / 4;
Ixx = 0.2 * 2^3 / 12 + 2 * (pi * r^2 / 2 + 4 * r^3 / 3 + pi * r^4 / 8);
for vs = [0, -0.5, 0.9828, 0.8932, 2.0, -1.2, 3.3]
    hs = mwecmass.solid.hydrostatics_at_draft(geo, vs);
    frames(hs, geo, 'C1');
    if strcmp(hs.submersion, 'partial')
        % mirror pairs of rows give exactly negated Green integrands, so xF Aw = 1/2 loop int x^2 dy
        % (and yF Aw) is a summation residue: its node count (16 Gauss points per span, pts holds 8
        % points per span) times eps times the sum of the absolute terms, at most max|x|^2 times
        % the perimeter
        X = hs.waterline.pts(:, 1:2);
        per = sum(sqrt(sum(diff([X; X(1, :)]).^2, 2)));
        bxy = 2 * size(X, 1) * eps * max(abs(X(:)))^2 * per;
        check(abs(hs.xF) * hs.Aw <= bxy && abs(hs.yF) * hs.Aw <= bxy, 'C1 vs = %g: xF, yF = %.3e %.3e', vs, hs.xF, hs.yF);
        check(hs.waterline.simple, 'C1 vs = %g: waterline not simple', vs);
    end
    fprintf(['C1 vs = %7.4f (%s, draft %.4f): V_sub %.12f m3, CB_body [%.1e %.1e %.12f], S_wet %.9f m2, Aw %.12f m2, ' ...
        'I_wp_xx %.9f, I_wp_yy %.12f m4, BM_L %.9f, KM %.12f m\n'], vs, hs.submersion, hs.draft, hs.V_sub, hs.CB_body, ...
        hs.S_wet, hs.Aw, hs.I_wp_xx, hs.I_wp_yy, hs.BM_L, hs.KM);
    if any(vs == [0, -0.5])
        fprintf('  stadium oracle: Aw %.15f (diff %.2e), I_wp_yy %.15f (diff %.2e), I_wp_xx %.15f (diff %.2e)\n', ...
            Aw, hs.Aw - Aw, Iyy, hs.I_wp_yy - Iyy, Ixx, hs.I_wp_xx - Ixx);
    end
end
end

function frames(hs, geo, name)
% I6 and the S7 cases at the ends of the hull, bitwise
check(hs.draft == -(geo.z_range(1) + hs.vs), '%s: draft', name);
zw = -hs.vs;
if zw <= geo.z_range(1)
    check(strcmp(hs.submersion, 'none') && hs.V_sub == 0 && hs.S_wet == 0 && hs.Aw == 0 && isnan(hs.KM) && ...
        isnan(hs.BM_L) && all(isnan(hs.CB)) && isempty(hs.waterline), '%s: submersion none', name);
    return
end
check(isequal(hs.CB, hs.CB_body + [0 0 hs.vs]) && hs.KM == hs.CB(3) + hs.BM_L, '%s vs = %g: I6 frames', name, hs.vs);
if zw >= geo.z_range(2)
    check(strcmp(hs.submersion, 'full') && hs.Aw == 0 && hs.BM_L == 0 && hs.I_wp_yy == 0 && hs.I_wp_xx == 0 && ...
        hs.KM == hs.CB(3) && isempty(hs.waterline), '%s: submersion full', name);
else
    check(strcmp(hs.submersion, 'partial') && hs.waterline.z == zw && hs.Aw == hs.waterline.area && ...
        hs.BM_L == hs.I_wp_yy / hs.V_sub, '%s vs = %g: partial', name, hs.vs);
end
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'standins'), '-end');
addpath(fullfile(root, 'tests', 'standins', 'fixtures'), '-end');
end

function check(cond, varargin)
if ~cond
    error('test_hydrostatics_at_draft:fail', varargin{:});
end
end
