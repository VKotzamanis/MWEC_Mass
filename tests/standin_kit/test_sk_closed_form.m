function test_sk_closed_form()
%TEST_SK_CLOSED_FORM  sti_closed_form V, S, J against Gauss-Legendre integration in z of its own
%   sections (piecewise constant in z, so the 2-node rule is exact), and the F7 stand-in.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
fc = sti_closed_form('fixture', 'cylinder');
fb = sti_closed_form('fixture', 'box');
ec = [-3; -2; -1; 0; 1];
eb = [-2.5; -1.5; -0.5; 0.5];
eps_fit = 0.01 * 0.1;
D = @(t) t + eps_fit / 2;
cases = {
    fc, 'modular_precast', ec, D([0.1; 0.1; 0.2; NaN]), -2.5, 4
    fc, 'modular_precast', ec, D([0.1; 0.2; 0.2; NaN]), -2, 4
    fc, 'modular_precast', ec, D([0.1; 0.1; 0.1; NaN]), -1.5, 4
    fc, 'modular_precast', ec, D([0.2; 0.1; 0.1; 0.1]), -2.95, []
    fc, 'thin_shell', ec, D(0.1 * ones(4, 1)), -2.5, []
    fc, 'thin_shell', ec, D(0.1 * ones(4, 1)), -1, []
    fc, 'thin_shell', ec, D(0.2 * ones(4, 1)), -3, []
    fc, 'thin_shell', ec, D(0.1 * ones(4, 1)), 0.95, []
    fb, 'modular_precast', eb, D([0.1; 0.15; NaN]), -2, 3
    fb, 'modular_precast', eb, D([0.15; 0.1; 0.1]), -1.5, []
    fb, 'thin_shell', eb, D(0.1 * ones(3, 1)), -1, []
    fb, 'thin_shell', eb, D(0.15 * ones(3, 1)), -2.45, []
    };
xg = [-1 1] / sqrt(3);
worst = 0;
for c = 1:size(cases, 1)
    [fx, mode, edges, d, zb, solid] = cases{c, :};
    design = struct('mode', mode, 'edges', edges, 'vs', 0, 't', d - eps_fit / 2, 'z_ballast', zb, ...
        'solid_modules', solid);
    reg = sti_closed_form('regions', fx, design, d);
    lay = sti_closed_form('layout', fx, design, d);
    names = setdiff(fieldnames(reg), {'V_module'});
    for i = 1:numel(edges) - 1
        br = unique([edges(i); edges(i + 1); clip([zb; lay(i).a; lay(i).b], edges(i), edges(i + 1))]);
        for r = 1:numel(names)
            V = 0; S = zeros(1, 3); J = zeros(3); aV = 0; aS = zeros(1, 3); aJ = zeros(3);
            for k = 1:numel(br) - 1
                h = br(k + 1) - br(k);
                for z = br(k) + (xg + 1) / 2 * h
                    s = sti_closed_form('region_sections', fx, design, d, i, z).(names{r});
                    w = h / 2;
                    v = w * s.A;
                    sv = w * [s.Sx, s.Sy, s.A * z];
                    jv = w * [s.Iyy, s.Ixy, s.Sx * z; s.Ixy, s.Ixx, s.Sy * z; s.Sx * z, s.Sy * z, s.A * z^2];
                    V = V + v; S = S + sv; J = J + jv;
                    aV = aV + abs(v); aS = aS + abs(sv); aJ = aJ + abs(jv);
                end
            end
            R = reg.(names{r});
            err = [abs(R.V(i) - V) / max(aV, realmin), max(abs(R.S(i, :) - S) ./ max(aS, realmin)), ...
                max(max(abs(R.J(:, :, i) - J) ./ max(aJ, realmin)))];
            worst = max([worst, err]);
            % two exact rules on the same piecewise-constant sections: they differ by the rounding
            % of at most a few dozen products and sums per subinterval
            check(all(err <= 64 * eps), 'case %d module %d region %s: closed form differs from Gauss-Legendre', ...
                c, i, names{r});
        end
    end
end
fprintf('V, S, J closed form vs Gauss-Legendre in z: largest difference relative to the absolute sum %.3e\n', worst);

% F7 stand-in: partial, at the ends, full and none
evalc('mc = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''tests'', ''standins'', ''fixtures'', ''cylinder.ms2''));');
gc = mwecmass.solid.outer_nurbs(mc);
s7 = {'vs', 'draft', 'submersion', 'V_sub', 'CB_body', 'CB', 'S_wet', 'Aw', 'xF', 'yF', 'I_wp_xx', ...
    'I_wp_yy', 'I_wp_yy_origin', 'BM_L', 'KM', 'waterline'};
for vs = [0.5, -1, -1.5, 3, 3.5]
    hs = mwecmass.solid.hydrostatics_at_draft(gc, vs, struct());
    check(isequal(sort(fieldnames(hs))', sort(s7)), 'S7 fields');
    check(hs.draft == -(-3 + vs), 'draft');
    switch hs.submersion
        case 'partial'
            check(vs > -1 && vs < 3, 'partial range');
            check(abs(hs.waterline.area - hs.Aw) <= 64 * eps * hs.Aw && hs.waterline.z == -vs, 'waterline loop');
            check(hs.CB(3) == hs.CB_body(3) + vs && hs.KM == hs.CB(3) + hs.BM_L, 'frames (I6)');
            check(hs.V_sub == pi * 1.5^2 * (-vs + 3), 'V_sub of the prism');
        case 'full'
            check(vs <= -1 && hs.Aw == 0 && hs.BM_L == 0 && isempty(hs.waterline) && hs.KM == hs.CB(3), 'full');
        case 'none'
            check(vs >= 3 && hs.V_sub == 0 && hs.S_wet == 0 && isnan(hs.KM) && isempty(hs.waterline), 'none');
    end
    fprintf('F7 vs = %5.2f: %-7s V_sub %.6f m^3, S_wet %.6f m^2, KM %.6f m\n', vs, hs.submersion, hs.V_sub, hs.S_wet, hs.KM);
end
end

function z = clip(z, lo, hi)
z = z(isfinite(z) & z > lo & z < hi);
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
    error('test_sk_closed_form:fail', varargin{:});
end
end
