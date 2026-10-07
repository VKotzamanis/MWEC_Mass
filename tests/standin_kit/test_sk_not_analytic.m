function test_sk_not_analytic()
%TEST_SK_NOT_ANALYTIC  Every stand-in errors mwecmass:standin:NotAnalytic on C1 data.
%   The C1 inputs: the parsed Input/C1.ms2; a geo with analytic = [] (as a real F1 returns); an
%   S1-like patch built from C1's curve1 (quadratic B-spline control points from MS2Parser) ruled to
%   its projection on y = 0, u-degree 2; a body and props without the stand-in marker; a config
%   whose hull_solid is that geo. F3 (eval_bspline_curve, eval_bspline_surface) is exempt: its
%   signature carries neither a geo nor a patch and its evaluation is exact for every NURBS.
%   J1 deletes the rows of the F1-F7 and F6b stand-ins it removes, J2 the rows of F9 and F10; the
%   rows of the fixture helpers sti_closed_form and sti_inner_box stay.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'standins'), '-end');
addpath(fullfile(root, 'tests', 'standins', 'fixtures'), '-end');

evalc('model = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''Input'', ''C1.ms2''));');
cache = mwecmass.geometry.precompute_boundary_cache(model, 21);
e = model.entities('curve1');
P = zeros(numel(e.params.ctrl_pt_names), 3);
for k = 1:size(P, 1)
    P(k, :) = model.eval_point(e.params.ctrl_pt_names{k});
end
Q = P;
Q(:, 2) = 0;
n = size(P, 1);
surf = struct('type', 'bspline', 'degree', [2 1], 'ctrl', cat(2, reshape(P, n, 1, 3), reshape(Q, n, 1, 3)), ...
    'knots', {{[0 0 0 (1:n - 3) / (n - 2) 1 1 1], [0 0 1 1]}}, 'weights', []);
patch = struct('name', 'curve1_ruled', 'source', 'curve1', 'type', 'RuledSurf', 'flips', {{}}, 'surf', surf, ...
    'outward', true, 'exact', true, 'z_of_u', true, 'u_range', [0 1], 'z_range', [P(1, 3) P(end, 3)], ...
    'offset_kind', 'ruled_parallel', 'pole', [false false], 'c0_u', [], 'c0_v', [], 'seam_u0', [], ...
    'seam_u1', [], 'seam_v0', [], 'seam_v1', []);
geo = struct('hull_name', 'C1', 'outer', patch, 'z_range', [-3.25 1.1], 'analytic', []);
edges = [-3.25; -2.706; -1.619; -0.531; 0.556; 1.10];
design = struct('mode', 'modular_precast', 'edges', edges, 'vs', 0.9828, 't', [0.0762; 0.0762; 0.0762; 0.0762; NaN], ...
    'z_ballast', -2.3, 'solid_modules', 5);
body = struct('design', design, 'inner_t', 0.0762, 'planes', edges(2:end - 1)', 'brep', [], 'shells', [], 'analytic', []);
props = struct('CG_total', [0 0 -1.04], 'GM_L', 0.2, 'periods', struct('heave', 3, 'pitch', 5.03), ...
    'mass_total', 20387.7, 'V_sub', 20.16);
stage2 = struct('vs', 0.9828, 'rho', [2500; 1430; 369.4; 369.4; 2500], 'mass', 20387.7, 'Z_CG', -1.04, ...
    'GM', 0.2, 'T_heave', 3, 'T_pitch', 5.03);
config = struct('hull_solid', geo, 'RHO_WATER', 1025, 'G', 9.80665);
opts = struct('t_min', 0.0762);

calls = {
    'outer_nurbs', @() mwecmass.solid.outer_nurbs(model)
    'outer_nurbs (3 inputs)', @() mwecmass.solid.outer_nurbs(model, cache, opts)
    'offset_surface', @() mwecmass.solid.offset_surface(model, cache, geo, 0.0762, [-3.25 1.1], opts)
    'void_closing_distance', @() mwecmass.solid.void_closing_distance(model, cache, geo, [-2.8 0.65])
    'split_bspline_surface', @() mwecmass.solid.split_bspline_surface(patch, -1)
    'slice_bspline_surface', @() mwecmass.solid.slice_bspline_surface(patch, -1)
    'build_body', @() mwecmass.solid.build_body(geo, design, [])
    'body_properties', @() mwecmass.solid.body_properties(body, struct('uhpc', 2500, 'air', 1.2), struct())
    'body_section', @() mwecmass.solid.body_section(body, -1)
    'body_section (below)', @() mwecmass.solid.body_section(body, -1, 'below')
    'hydrostatics_at_draft', @() mwecmass.solid.hydrostatics_at_draft(geo, 0.9828, struct())
    'evaluate_realised', @() mwecmass.realise.evaluate_realised(struct(), struct(), design, config)
    'check_against_stage2', @() mwecmass.realise.check_against_stage2(props, stage2, 10, 1e-6, 1025)
    'sti_closed_form fixture', @() sti_closed_form('fixture', geo)
    'sti_inner_box', @() sti_inner_box(geo, 0.0762, 0.0762)
    };
for k = 1:size(calls, 1)
    id = '';
    try
        calls{k, 2}();
    catch err
        id = err.identifier;
    end
    fprintf('%-24s on C1: %s\n', calls{k, 1}, id);
    if ~strcmp(id, 'mwecmass:standin:NotAnalytic')
        error('test_sk_not_analytic:fail', '%s did not error NotAnalytic on C1 (got "%s")', calls{k, 1}, id);
    end
end
end
