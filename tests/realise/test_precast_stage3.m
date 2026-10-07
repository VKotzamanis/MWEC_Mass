function test_precast_stage3()
%TEST_PRECAST_STAGE3  Modular-precast Stage 3 (split, build, check, store) on the SK fixtures.
%   Kernel F1-F7 are the SK stand-ins until J1 (cylinder inner sets through F2, box inner sets
%   through sti_inner_box); Stage 2 is sti_stage2 (closed-form prisms), not an optimiser run.
%   After J1 the same cases run on the real kernel (join test); the hull-volume sum against the
%   closed form is then printed only (rational faces, Gauss quadrature).
%   Asserted: S8 layout, final_props from the realised body, volume closure (rule 11; z_ballast
%   inside k*, at the hull bottom e(1) below the inner z_lo in case E, at the interior joint
%   e(k*+1) in case E2), shell and ballast bounds, flotation held to the Stage-3 constraint
%   tolerance 1e-6 (AGENTS OD10), the status rule, the restart of the shell root search on a
%   changed knot structure. Root-finding residuals and Stage-2 deviations are printed.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
s8 = {'mode', 'hull_name', 'status', 'reason', 'escalation', 'vs', 'draft', 'stage2', 'rho', 'design', ...
    'k_star', 'V_uhpc_target', 'modules', 'props', 'check', 'solver', 'fit', 'body', 'step_files'};
mod_fields = {'z_lo', 'z_hi', 't', 'h_ballast', 'V', 'V_uhpc', 'V_air', 'mass', 'rho_eff', 'rho_stage2', ...
    'rho_floor', 'CG_world'};
tol_eq = 1e-6;

% A: cylinder through F14 (run.m), no wall module; F2 inner sets
config = fixture_config('cylinder', [], 2500);
vs = 0.5;
[f3, s2] = sti_stage2(config, vs, [NaN; 300; 270; 440]);
opt = struct('Final3D', f3, 'stage2_3d', struct('properties', f3), 'constructability', []);
[results, fp] = mwecmass.realise.modular_precast.run(config, [vs; s2.rho], opt);
r = results.stage3;
check(isequaln(results.Final3D, f3) && isequaln(results.stage2_3d.properties, fp), 'F14 results layout');
verify(r, fp, config, s2, s8, mod_fields, tol_eq, 'cylinder');
check(isequal(r.design.solid_modules, zeros(1, 0)) && r.k_star == 1, 'cylinder roles');
fx = sti_closed_form('fixture', 'cylinder');
eps_fit = 0.01 * config.constructability_t_min;
for i = 2:3
    % oracle: an interior hollow module of the cylinder holds air pi (R - d)^2 h (closed form)
    h = r.modules(i).z_hi - r.modules(i).z_lo;
    d = fx.R - sqrt((r.modules(i).V - r.V_uhpc_target(i)) / (pi * h));
    fprintf('cylinder module %d: t = %.12f m, closed-form t = %.12f m, difference %.3g m\n', i, ...
        r.modules(i).t, d - eps_fit / 2, r.modules(i).t - (d - eps_fit / 2));
end

% B: box with a solid wall module on top (test material rho_uhpc = 1200 so that the box floats)
config = fixture_config('box', 3, 1200);
[f3, s2] = sti_stage2(config, 0, [NaN; 450; 1200]);
r = stage3(config, [0; s2.rho], f3);
verify(r, [], config, s2, s8, mod_fields, tol_eq, 'box, wall');
check(isequal(r.design.solid_modules, 3) && r.modules(3).V_air == 0 && isnan(r.modules(3).t) && ...
    r.modules(3).mass == r.modules(3).V * 1200, 'wall module solid');

% C: box, no wall
config = fixture_config('box', [], 2500);
[f3, s2] = sti_stage2(config, 0, [NaN; 600; 600]);
x = [0; s2.rho];
r = stage3(config, x, f3);
verify(r, [], config, s2, s8, mod_fields, tol_eq, 'box');
r_c = r;

% D: a Stage 2 equal to the realised values is accepted with zero deviation (same computation)
f3e = f3;
f3e.CG_total = r.props.CG_total;
f3e.GM_L = r.props.GM_L;
f3e.periods = r.props.periods;
[re, fpe] = mwecmass.realise.modular_precast.solve_and_extract(config, x, f3e, inner_opts(config));
check(strcmp(re.status, 'accepted') && isempty(re.reason) && all([re.check.metrics.rel_dev] == 0) && ...
    strcmp(fpe.stage3_status, 'accepted') && isequal(re.design, r.design), 'identity accepted');

% E: flotation out of reach inside k*: closest level kept, failed, stored
[f3, s2] = sti_stage2(config, 2.4, [100; 100; 100]);
[r, fp] = mwecmass.realise.modular_precast.solve_and_extract(config, [2.4; s2.rho], f3, inner_opts(config));
check(strcmp(r.status, 'failed') && any(strcmp(r.check.failed, 'flotation')) && ...
    ~isempty(strfind(r.reason, 'flotation')) && r.design.z_ballast == r.design.edges(1) && ...
    r.solver.exitflag == -2, 'unreachable flotation reported');
check(fp.mass_total == r.props.mass_total && fp.mass_total ~= f3.mass_total && ...
    strcmp(fp.stage3_status, 'failed'), 'closest design stored, not Stage 2');
fprintf('box, flotation out of reach: mass %.3f kg vs displaced %.3f kg; reason: %s\n', ...
    r.props.mass_total, r.props.mass_buoyant_force, r.reason);
% z_ballast = e(1) is a module edge and lies below the inner z_lo (contract I1 positions)
inner_e = sti_inner_box(config.hull_solid, r.design.t(1), config.constructability_t_min);
check(r.design.z_ballast == r.design.edges(1) && r.design.z_ballast < inner_e.z_lo, 'case E positions');
closure(r, config, 'box, flotation out of reach');

% E2: too light even with k* full: the ballast is kept at the interior joint e(k*+1) (I1 module edge)
config_l = fixture_config('box', [], 1200);
[f3, s2] = sti_stage2(config_l, 0, [1000; 210; 285]);
[r, fp] = mwecmass.realise.modular_precast.solve_and_extract(config_l, [0; s2.rho], f3, ...
    inner_opts(config_l));
k = r.k_star;
m = r.modules(k);
check(r.design.z_ballast == r.design.edges(k + 1) && k + 1 < numel(r.design.edges), 'joint position');
check(m.V_air == 0 && isnan(m.t) && m.h_ballast == m.z_hi - m.z_lo, 'k* full of ballast');
check(strcmp(r.status, 'failed') && any(strcmp(r.check.failed, 'flotation')) && ...
    r.solver.exitflag == -2 && strcmp(fp.stage3_status, 'failed'), 'too light reported');
fprintf('box, too light: k* = %d, z_ballast %.6f m = e(%d); mass %.3f kg vs displaced %.3f kg\n', ...
    k, r.design.z_ballast, k + 1, r.props.mass_total, r.props.mass_buoyant_force);
closure(r, config_l, 'box, ballast at the joint e(k*+1)');

% F: the adaptive fit above t_thr carries an extra knot, so the root found on the t_min knots has
% another knot structure and the search restarts on it (contract section 7 item 4). t_thr lies
% between t_min and the module-2 root of case C (109.2 mm) and above the module-3 root (78.1 mm).
t_thr = 0.09;
geo = config.hull_solid;
t_min = config.constructability_t_min;
knotted_box('reset');
opts = struct('inner_fn', @(t, knots_from) knotted_box(geo, t, t_min, t_thr, knots_from));
[f3, s2] = sti_stage2(config, 0, [NaN; 600; 600]);
rk = mwecmass.realise.modular_precast.solve_and_extract(config, [0; s2.rho], f3, opts);
verify(rk, [], config, s2, s8, mod_fields, tol_eq, 'box, knot restart');
calls = knotted_box('log');
check(any(calls(:, 2) == 1), 'a root search ran on the knots of the adaptive fit at the root');
check(rk.design.t(2) > t_thr && rk.design.t(3) < t_thr, 'module roots on both sides of t_thr');
fprintf('box, knot restart: %d inner-set calls on the inserted knot; t %s mm (case C %s mm)\n', ...
    sum(calls(:, 2) == 1), mat2str(1000 * rk.design.t', 9), mat2str(1000 * r_c.design.t', 9));
dev_k = [rk.modules.V_uhpc]' - rk.V_uhpc_target;
dev_c = [r_c.modules.V_uhpc]' - r_c.V_uhpc_target;
fprintf('box, knot restart: V_uhpc - target %s m^3 (case C %s m^3)\n', mat2str(dev_k', 3), ...
    mat2str(dev_c', 3));
if ~isempty(geo.analytic)
    % stand-in kernel: the box volumes are closed forms of d, blind to the inserted knot, so the
    % restarted root equals the case-C root bitwise
    check(isequal(rk.design.t, r_c.design.t) && isequal(dev_k, dev_c), 'restart reproduces case C');
end
end

function out = knotted_box(geo, t, t_min, t_thr, knots_from)
% sti_inner_box with the knot u = 1/2 inserted in every side patch when t > t_thr (adaptive fit),
% or when knots_from carries it (refit on fixed knots). Calls are logged as [t, refit on it].
persistent calls
if ischar(geo)
    if strcmp(geo, 'reset')
        calls = zeros(0, 2);
    end
    out = calls;
    return
end
out = sti_inner_box(geo, t, t_min);
if isempty(knots_from)
    knotted = t > t_thr;
else
    knotted = any(arrayfun(@(q) numel(q.surf.knots{1}) == 5, knots_from.patches));
    out.refit = true;
end
calls(end + 1, :) = [t, ~isempty(knots_from) && knotted];
if ~knotted
    return
end
for k = 1:numel(out.patches)
    s = out.patches(k).surf;
    if s.ctrl(1, 1, 3) == s.ctrl(2, 1, 3)
        continue
    end
    % degree-1 insertion of u = 1/2 in [0 0 1 1]: the new row is the mean of the two rows
    s.ctrl = cat(1, s.ctrl(1, :, :), (s.ctrl(1, :, :) + s.ctrl(2, :, :)) / 2, s.ctrl(2, :, :));
    s.knots{1} = [0 0 0.5 1 1];
    out.patches(k).surf = s;
    out.report.patches(k).n_knots = [5 4];
end
end

function config = fixture_config(name, wall, rho_uhpc)
config = sti_config(name);
config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
config.wall_strip_index = wall;
config.constructability_rho_hull = rho_uhpc;
config.constructability_rho_air = 1.2;
config.constructability_t_min = 0.0762;
config.output.save.stage3 = struct('precast_midplane', false, 'precast_strips', false, 'step', false);
end

function opts = inner_opts(config)
opts = struct();
if strcmp(config.hull_solid.hull_name, 'box')
    t_min = config.constructability_t_min;
    geo = config.hull_solid;
    opts.inner_fn = @(t, knots_from) sti_inner_box(geo, t, t_min);
end
end

function r = stage3(config, x, f3)
r = mwecmass.realise.modular_precast.solve_and_extract(config, x, f3, inner_opts(config));
end

function verify(r, fp, config, s2, s8, mod_fields, tol_eq, label)
check(isequal(sort(fieldnames(r))', sort(s8)), [label ': S8 fields']);
check(isequal(fieldnames(r.modules)', mod_fields), [label ': S8 modules fields']);
check(strcmp(r.mode, 'modular_precast') && strcmp(r.escalation, 'split') && ...
    isequal(r.stage2, s2) && r.vs == s2.vs && r.design.vs == s2.vs, [label ': S8 header']);
bp = mwecmass.solid.body_properties(r.body, r.rho, struct());
check(isequaln(r.props.mass_total, bp.total.mass) && r.props.CG_total(3) == bp.total.CG_body(3) + r.vs, ...
    [label ': props from the stored body']);
if ~isempty(fp)
    check(isequal(rmfield(fp, {'stage3_status', 'stage3_check'}), r.props) && ...
        strcmp(fp.stage3_status, r.status) && isequal(fp.stage3_check, r.check), [label ': final_props']);
end
check(strcmp(r.status, 'accepted') == r.check.pass && isempty(r.reason) == r.check.pass, [label ': status rule']);
check(abs(r.check.equalities(1).residual) <= tol_eq, [label ': flotation']);
closure(r, config, label);
N = numel(r.modules);
e = r.design.edges;
k = r.k_star;
check(r.design.z_ballast >= e(k) && r.design.z_ballast <= e(k + 1), [label ': ballast inside k*']);
t_min = config.constructability_t_min;
hollow = setdiff(k + 1:N, r.design.solid_modules);
check(all(r.design.t(hollow) >= t_min) && r.design.t(k) == t_min, [label ': shells >= t_min, t_k* = t_min']);
check(all([r.modules(r.design.solid_modules).V_air] == 0), [label ': solid modules']);
dev = [r.modules.V_uhpc]' - r.V_uhpc_target;
fprintf('%s: status %s, z_ballast %.6f m, t %s mm, V_uhpc - target %s m^3, flotation residual %.3g\n', ...
    label, r.status, r.design.z_ballast, mat2str(1000 * r.design.t', 6), mat2str(dev', 3), ...
    r.check.equalities(1).residual);
fprintf('%s: rel. deviations Z_CG %.4g, GM %.4g, T_heave %.4g, T_pitch %.4g; GM equality residual %.4g\n', ...
    label, r.check.metrics.rel_dev, r.check.equalities(2).residual);
end

function closure(r, config, label)
N = numel(r.modules);
hull = sti_closed_form('hull', sti_closed_form('fixture', config.hull_solid));
Vsum = 0;
for i = 1:N
    m = r.modules(i);
    % V is the sum of the two region volumes: equal to one rounding of that sum
    check(abs(m.V_uhpc + m.V_air - m.V) <= 2 * eps * m.V, sprintf('%s: module %d closure', label, i));
    fprintf('%s: module %d V_uhpc + V_air - V = %.3g m^3\n', label, i, m.V_uhpc + m.V_air - m.V);
    Vsum = Vsum + m.V;
end
fprintf('%s: sum of module volumes - closed-form hull volume = %.3g m^3\n', label, Vsum - hull.V);
if ~isempty(config.hull_solid.analytic)
    % stand-in kernel: N closed-form prism volumes summed, within N roundings of the hull volume
    check(abs(Vsum - hull.V) <= N * eps * hull.V, [label ': module volumes sum to the hull']);
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

function check(cond, msg)
if ~cond
    error('test_precast_stage3:fail', '%s', msg);
end
end
