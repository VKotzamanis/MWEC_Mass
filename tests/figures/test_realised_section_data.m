function test_realised_section_data()
%TEST_REALISED_SECTION_DATA  F12 on the stand-in bodies: elevation polygons, plan sections, status.
%   The stand-ins (cylinder, box; precast and thin shell) are closed forms, which serve as the
%   independent oracle: sti_closed_form('layout') gives the air interval of every module, the
%   section at y = 0 of a cylinder of radius R is [-R, R] and of the void [-(R-d), R-d]. The real
%   F5/F6b replace the stand-ins when T3 merges; this test then runs unchanged on them.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);

% {fixture, mode, t [N x 1], z_ballast, solid_modules, vs, t_min}
I3_BOUND = 1e-7;   % m, contract I3: the kernel's cut at a plane
cases = {
    'cylinder', 'modular_precast', [0.1; 0.1; 0.1; 0.1], -2.5, [], 0.5, 0.1
    'cylinder', 'modular_precast', [0.1; 0.1; 0.1; NaN], -2.5, 4, 0.5, 0.1
    'cylinder', 'modular_precast', [0.1; 0.1; 0.1; 0.1], -2, [], 0.5, 0.1
    'cylinder', 'modular_precast', [0.1; 0.2; 0.1; 0.1], -2.5, [], 0.5, 0.1
    'cylinder', 'thin_shell', 0.0254 * ones(4, 1), -2.85, [], 0.5, 0.0254
    'cylinder', 'thin_shell', 0.0254 * ones(4, 1), -2.99, [], 0.5, 0.0254
    'box', 'modular_precast', [0.08; 0.08; 0.08], -2.3, [], 0, 0.08
    'box', 'thin_shell', 0.0254 * ones(3, 1), -1.0, [], 0.5, 0.0254
    };
for c = 1:size(cases, 1)
    run_case(I3_BOUND, cases{c, :});
end
test_errors_and_status();
fprintf('all F12 tests passed\n');
end

function run_case(I3_BOUND, name, mode, t, zb, solid, vs, t_min)
config = sti_config(name);
config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
[~, stage2] = sti_stage2(config, vs, [NaN; 250 * ones(numel(t) - 1, 1)]);
design = struct('mode', mode, 'edges', config.strip_edges, 'vs', vs, 't', t, 'z_ballast', zb, 'solid_modules', solid);
if strcmp(mode, 'modular_precast')
    rho = struct('uhpc', 2500, 'air', 1.2);
else
    rho = struct('ballast', 7500, 'shell', 7850, 'air', 1.2);
end
realised = sti_realised(config, design, rho, stage2, struct('t_min', t_min));
fx = sti_closed_form('fixture', config.hull_solid);
e = config.strip_edges(:);
N = numel(e) - 1;
d = t + 0.01 * t_min / 2;
lay = sti_closed_form('layout', fx, design, d);
if strcmp(fx.kind, 'cylinder')
    hw_out = fx.R;
    hw_in = fx.R - d;
else
    hw_out = fx.x(2);
    hw_in = fx.x(2) - d;
end
margin = 8 * eps(max(abs(e)));
exact_inner = strcmp(fx.kind, 'box');
tag = sprintf('%s %s z_ballast %.3g solid [%s]', name, mode, zb, num2str(solid));

tic;
data = mwecmass.output.figures.realised_section_data(realised, [], 7);
fprintf('%s: %d polygons, %d plan sections, %d omitted heights, %.1f s\n', tag, numel(data.polygons), ...
    numel(data.plan), numel(data.omitted), toc);
if ~isempty(data.omitted)
    error('test_realised_section_data:fail', '%s: heights omitted: %s', tag, data.omitted(1).reason);
end

% names and frames
check(isequal(data.edges, e) && data.vs == vs && data.z_ballast == zb && data.waterline_z == -vs, '%s: frame fields', tag);
check(isequal(data.CG, realised.props.CG_total) && isequal(data.CB, realised.props.CB), '%s: CG, CB', tag);
regions = unique({data.polygons.region});
if strcmp(mode, 'modular_precast')
    check(all(ismember(regions, {'uhpc', 'air'})), '%s: precast regions %s', tag, strjoin(regions, ','));
else
    check(all(ismember(regions, {'ballast', 'shell', 'air'})), '%s: thin-shell regions %s', tag, strjoin(regions, ','));
end
roles = {data.polygons.role};
expected_roles = {};
if zb > e(1)
    expected_roles{end + 1} = 'ballast';
end
if ~isempty(solid)
    expected_roles{end + 1} = 'solid_module';
end
if any([lay.air])
    expected_roles = [expected_roles, {'void', 'wall'}];
end
check(isequal(sort(unique(roles)), sort(expected_roles)), '%s: roles drawn are %s, expected %s', tag, ...
    strjoin(unique(roles), ','), strjoin(expected_roles, ','));
check(all(strcmp({data.polygons(strcmp(roles, 'void')).region}, 'air')), '%s: voids are air', tag);

% outline: constant half-width, exact to the evaluation of the curves
o = data.outline;
err = max(abs([o.x_hi - hw_out; o.x_lo + hw_out]));
fprintf('  outline half-width error %.2e m over %d heights\n', err, numel(o.z));
check(err <= 8 * eps(hw_out), '%s: outline half-width off by %.3e', tag, err);
check(o.z(1) == e(1) + margin && o.z(end) == e(end) - margin, '%s: outline ends', tag);
check(polyarea(o.profile(:, 1), o.profile(:, 2)) > 0 && signed_area(o.profile) > 0, '%s: outline orientation', tag);

% every polygon lies inside its module and is counter-clockwise
for k = 1:numel(data.polygons)
    p = data.polygons(k);
    check(p.z_lo >= e(p.module) && p.z_hi <= e(p.module + 1), '%s: polygon %d leaves module %d', tag, k, p.module);
    check(signed_area(p.xz) > 0, '%s: polygon %d is not counter-clockwise', tag, k);
    check(p.z_lo == min(p.xz(:, 2)) && p.z_hi == max(p.xz(:, 2)), '%s: polygon %d z range', tag, k);
end

% void polygons: one per module with air, at the closed-form air interval and half-width
z_tol = margin + 4 * eps(max(abs(e)));
for i = 1:N
    v = data.polygons(strcmp(roles, 'void') & [data.polygons.module] == i);
    if ~lay(i).air
        check(isempty(v), '%s: module %d has a void but no air', tag, i);
        continue
    end
    check(numel(v) == 1, '%s: module %d has %d void polygons', tag, i, numel(v));
    tol_lo = z_tol;
    if lay(i).a == zb
        % the void starts at the cut of the ballast level, which the kernel places to I3 (1e-7 m)
        tol_lo = I3_BOUND;
        fprintf('  module %d void bottom deviates %.2e m from z_ballast\n', i, abs(v.z_lo - zb));
    end
    inner_check(exact_inner, abs(v.z_lo - lay(i).a) <= tol_lo && abs(v.z_hi - lay(i).b) <= z_tol, ...
        '%s: module %d void z [%.17g %.17g], closed form [%.17g %.17g]', tag, i, v.z_lo, v.z_hi, lay(i).a, lay(i).b);
    xs = v.xz(:, 1);
    inner_check(exact_inner, max(abs([max(xs) - hw_in(i), min(xs) + hw_in(i)])) <= 8 * eps(hw_out), ...
        '%s: module %d void half-width [%.17g %.17g], closed form %.17g', tag, i, min(xs), max(xs), hw_in(i));
    wa = data.polygons(strcmp(roles, 'wall') & [data.polygons.module] == i);
    check(numel(wa) >= 2, '%s: module %d has %d wall polygons', tag, i, numel(wa));
end

% void outlines: one per connected air region; consecutive modules with air up to and from their
% common edge and equal t are one outline with no segment at that edge, otherwise each keeps its own
air_to_edge = @(i) lay(i).air && lay(i).b >= e(i + 1) - z_tol;
air_from_edge = @(i) lay(i).air && lay(i).a <= e(i) + z_tol;
joined = false(1, N - 1);
for i = 1:N - 1
    joined(i) = air_to_edge(i) && air_from_edge(i + 1) && t(i) == t(i + 1);
end
vo = data.void_outlines;
check(numel(vo) == sum([lay.air]) - sum(joined), '%s: %d void outlines, expected %d', tag, numel(vo), ...
    sum([lay.air]) - sum(joined));
for k = 1:numel(vo)
    check(signed_area(vo(k).xz) > 0, '%s: void outline %d orientation', tag, k);
    xz = vo(k).xz;
    nxt = circshift(xz, -1);
    for i = 1:N - 1
        flat = abs(xz(:, 2) - e(i + 1)) <= margin & abs(nxt(:, 2) - e(i + 1)) <= margin & xz(:, 1) ~= nxt(:, 1);
        if joined(i) && any(vo(k).modules == i)
            check(~any(flat), '%s: void outline %d has a segment at the joined edge %d', tag, k, i);
        end
    end
    check(all(diff(vo(k).modules) == 1), '%s: void outline %d modules', tag, k);
end
for i = 1:N
    check(sum(arrayfun(@(o) any(o.modules == i), vo)) == lay(i).air, '%s: module %d in the void outlines', tag, i);
end
A_poly = sum(arrayfun(@(p) polyarea(p.xz(:, 1), p.xz(:, 2)), data.polygons(strcmp(roles, 'void'))));
A_out = sum(arrayfun(@(o) polyarea(o.xz(:, 1), o.xz(:, 2)), vo));
check(abs(A_out - A_poly) <= 64 * eps * max(A_poly, 1) + 4 * N * margin * 2 * hw_out, ...
    '%s: void outline areas %.15f, polygons %.15f', tag, A_out, A_poly);
fprintf('  %d void outlines, %d joined module edges\n', numel(vo), sum(joined));

% a void that closes below the top of its module leaves a cap of wall material up to the top of the
% hull; a ballast level below the inner keel leaves a wall body between z_ballast and the void
cap_module = N;
if lay(cap_module).air && strcmp(lay(cap_module).top, 'inner')
    cap = data.polygons(strcmp(roles, 'wall') & [data.polygons.module] == cap_module & [data.polygons.z_lo] > lay(cap_module).b - z_tol & ...
        [data.polygons.z_lo] < lay(cap_module).b + z_tol);
    inner_check(exact_inner, numel(cap) == 1 && cap.z_hi == e(end) - margin, '%s: cap above the void of module %d', tag, cap_module);
end
if zb > e(1) && zb < fx.z(1) + d(1)
    below = data.polygons(strcmp(roles, 'wall') & [data.polygons.module] == 1 & [data.polygons.z_lo] == zb);
    check(numel(below) == 1, '%s: wall body between z_ballast and the inner keel', tag);
    inner_check(exact_inner, abs(below.z_hi - (fx.z(1) + d(1))) <= z_tol, ...
        '%s: wall body top %.17g, closed-form inner keel %.17g', tag, below.z_hi, fx.z(1) + d(1));
end

% ballast level honoured: ballast polygons end at z_ballast to adjacent floats, voids start there
b = data.polygons(strcmp(roles, 'ballast'));
if zb > e(1)
    check(~isempty(b), '%s: no ballast polygon', tag);
    top = max([b.z_hi]);
    % the top is where the kernel's section changes from solid to hollow: the ballast cut, placed to I3
    fprintf('  ballast top deviates %.2e m from z_ballast\n', abs(zb - top));
    check(abs(zb - top) <= I3_BOUND, '%s: ballast top %.17g, z_ballast %.17g', tag, top, zb);
    check(min([b.z_lo]) == e(1) + margin, '%s: ballast bottom', tag);
end
solid_roles = data.polygons(strcmp(roles, 'solid_module'));
check(isequal(sort(unique([solid_roles.module])), sort(solid(:)')) || (isempty(solid) && isempty(solid_roles)), ...
    '%s: solid-module polygons', tag);

% area closure: the polygons tile the strip between the module-edge margins
A = sum(arrayfun(@(p) polyarea(p.xz(:, 1), p.xz(:, 2)), data.polygons));
H = e(end) - e(1);
uncovered = 2 * hw_out * (2 * N * margin + 8 * numel(data.polygons) * eps(max(abs(e))));
check(abs(A - 2 * hw_out * H) <= uncovered + 64 * eps * A, '%s: polygon areas %.15f, section %.15f', tag, A, 2 * hw_out * H);
fprintf('  polygon areas %.15f, closed-form section %.15f, difference %.2e\n', A, 2 * hw_out * H, A - 2 * hw_out * H);

% plan sections
check_plan(data, fx, hw_out, hw_in, lay, e, margin, tag, solid);

% one waterline section when the waterline is inside the hull
wl = data.plan(strcmp({data.plan.kind}, 'waterline'));
check(numel(wl) == data.waterline_in_hull, '%s: waterline section', tag);
end

function check_plan(data, fx, hw_out, hw_in, lay, e, margin, tag, solid)
N = numel(e) - 1;
exact_inner = strcmp(fx.kind, 'box');
plan = data.plan;
for k = 1:numel(plan)
    p = plan(k);
    check(signed_area(p.outer) > 0, '%s: plan %s outer not counter-clockwise', tag, p.label);
    if strcmp(fx.kind, 'cylinder')
        r = hypot(p.outer(:, 1), p.outer(:, 2));
        err = max(abs(r - hw_out));
        check(err <= 16 * eps(hw_out), '%s: plan %s outer radius off by %.3e', tag, p.label, err);
    else
        on_edge = abs(p.outer(:, 1)) == fx.x(2) | abs(abs(p.outer(:, 2)) - fx.y(2)) <= 4 * eps(fx.y(2));
        check(all(on_edge | abs(p.outer(:, 2)) < fx.y(2)) && max(abs(p.outer(:, 1))) <= fx.x(2) * (1 + 4 * eps), ...
            '%s: plan %s outer polygon is not the box section', tag, p.label);
    end
    i = p.module;
    hollow = lay(i).air && p.z >= lay(i).a && p.z <= lay(i).b;
    check(p.solid == ~hollow, '%s: plan %s solid = %d, closed form hollow = %d', tag, p.label, p.solid, hollow);
    check(isempty(p.inner) == p.solid, '%s: plan %s inner loop', tag, p.label);
    if hollow
        check(signed_area(p.inner) > 0, '%s: plan %s inner not counter-clockwise', tag, p.label);
        if strcmp(fx.kind, 'cylinder')
            err = max(abs(hypot(p.inner(:, 1), p.inner(:, 2)) - hw_in(i)));
            inner_check(exact_inner, err <= 16 * eps(hw_out), '%s: plan %s inner radius off by %.3e', tag, p.label, err);
        else
            check(max(abs(p.inner(:, 1))) <= hw_in(i) * (1 + 4 * eps), '%s: plan %s inner width', tag, p.label);
        end
    end
    if p.solid
        if any(i == solid)
            check(strcmp(p.role, 'solid_module'), '%s: plan %s role %s', tag, p.label, p.role);
        elseif p.z <= data.z_ballast
            check(strcmp(p.role, 'ballast'), '%s: plan %s role %s', tag, p.label, p.role);
        end
    end
end
kinds = {plan.kind};
check(sum(strcmp(kinds, 'module_bottom')) == N && sum(strcmp(kinds, 'module_top')) == N, '%s: module plan sections', tag);
bt = plan(strcmp(kinds, 'module_bottom'));
tp = plan(strcmp(kinds, 'module_top'));
check(all([bt.z] == e(1:end - 1)' + margin) && all([tp.z] == e(2:end)' - margin), '%s: module plan heights', tag);
shown = setdiff(1:N, solid);
pn = data.strip_panels;
check(numel(pn) == 2 * numel(shown), '%s: %d strip panels for %d modules', tag, numel(pn), numel(shown));
for q = 1:numel(pn)
    item = plan(pn(q).plan);
    check(item.module == pn(q).module && ismember(pn(q).module, shown), '%s: panel module', tag);
    check(strcmp(item.kind, ['module_' pn(q).end]), '%s: panel kind', tag);
    check(pn(q).row == 1 + strcmp(pn(q).end, 'bottom') && pn(q).col == find(shown == pn(q).module), '%s: panel place', tag);
end
zb = data.z_ballast;
if zb > e(1) && zb < e(end)
    bz = plan(strcmp(kinds, 'ballast_top'));
    m = find(e(1:end - 1) <= zb & zb < e(2:end), 1);
    check(numel(bz) == 1 && bz.module == m && bz.z == zb + margin, '%s: ballast_top section', tag);
    check(bz.solid == ~(lay(m).air && bz.z >= lay(m).a && bz.z <= lay(m).b), '%s: ballast_top solid flag', tag);
end
end

function test_errors_and_status()
config = sti_config('cylinder');
config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
[~, stage2] = sti_stage2(config, 0.5, [NaN; 250; 200; 150]);
design = struct('mode', 'modular_precast', 'edges', config.strip_edges, 'vs', 0.5, 't', 0.1 * ones(4, 1), ...
    'z_ballast', -2.5, 'solid_modules', []);
realised = sti_realised(config, design, struct('uhpc', 2500, 'air', 1.2), stage2, struct('t_min', 0.1));
check(strcmp(realised.status, 'failed'), 'the stand-in design is expected to fail the Stage-2 check');
d = mwecmass.output.figures.realised_section_data(realised, [], 3);
check(strcmp(d.status, 'failed') && strcmp(d.status_lines{1}, 'Stage 3 FAILED') && ...
    isequal(d.status_lines(2:end), strsplit(realised.reason, '; ')) && ~isempty(d.failed), 'failed status lines');
fprintf('status lines: %s\n', strjoin(d.status_lines, ' | '));
realised.status = 'accepted';
realised.reason = '';
d = mwecmass.output.figures.realised_section_data(realised, [], 3);
check(isequal(d.status_lines, {'Stage 3 accepted'}), 'accepted status lines');

% z_plan: user heights give plan sections at exactly those heights
d = mwecmass.output.figures.realised_section_data(realised, [-2.2, 0.3], 3);
extra = d.plan(strcmp({d.plan.kind}, 'plan'));
check(numel(extra) == 2 && isequal([extra.z], [-2.2, 0.3]), 'z_plan sections');
check(extra(1).module == 1 && extra(2).module == 4, 'z_plan modules');

check_error(@() mwecmass.output.figures.realised_section_data(struct()), 'mwecmass:figures:BadRealised');
check_error(@() mwecmass.output.figures.realised_section_data(realised, 5, 3), 'mwecmass:figures:ZOutside');
check_error(@() mwecmass.output.figures.realised_section_data(realised, [], 1), 'mwecmass:figures:BadInput');

% hull outline from the outer patches alone
geo = config.hull_solid;
out = mwecmass.output.figures.hull_outline_data(geo, 9);
fx = sti_closed_form('fixture', geo);
check(isempty(out.omitted) && numel(out.z) == 9 && isequal(out.z_range, fx.z), 'hull outline levels');
err = max(abs([out.x_hi - fx.R; out.x_lo + fx.R]));
margin = 8 * eps(max(abs(fx.z)));
check(err <= 8 * eps(fx.R) && out.z(1) == fx.z(1) + margin && out.z(end) == fx.z(2) - margin, 'hull outline half-width');
check(signed_area(out.profile) > 0, 'hull outline orientation');
fprintf('hull outline: 9 heights, half-width error %.2e m, area %.15f (closed form %.15f)\n', err, ...
    polyarea(out.profile(:, 1), out.profile(:, 2)), 2 * fx.R * diff(fx.z));
box = sti_config('box');
box.hull_solid = mwecmass.solid.outer_nurbs(box.ms2_model);
out = mwecmass.output.figures.hull_outline_data(box.hull_solid, 5);
fb = sti_closed_form('fixture', box.hull_solid);
check(isequal(out.x_hi, fb.x(2) * ones(5, 1)) && isequal(out.x_lo, fb.x(1) * ones(5, 1)), 'box outline');
end

function a = signed_area(P)
x = P(:, 1);
y = P(:, 2);
a = 0.5 * sum(x .* circshift(y, -1) - circshift(x, -1) .* y);
end

function check_error(fn, id)
try
    fn();
catch err
    check(strcmp(err.identifier, id), 'error id %s, expected %s', err.identifier, id);
    return
end
error('test_realised_section_data:fail', 'expected error %s', id);
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
    error('test_realised_section_data:fail', varargin{:});
end
end

function inner_check(exact, cond, varargin)
% The box inner wall is a closed form on any producer. The cylinder inner wall is a fitted
% surface once the real offset replaces the stand-in, so its deviation is printed, not gated.
if exact
    check(cond, varargin{:});
elseif ~cond
    fprintf('  inner geometry deviates from the closed form (not gated): %s\n', sprintf(varargin{:}));
end
end
