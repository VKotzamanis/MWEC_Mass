function test_stage3_export_schema()
%TEST_STAGE3_EXPORT_SCHEMA  results.stage3 (S8) and final_props through export_results and
%   check_export_schema. The S8 record is sti_realised of a five-module cylinder design (the
%   schema's fixed-shape check expects five components, as C1); the kernel is the SK stand-ins.
%   Also prints the Report Stage-3 panel of that record.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
config = sti_config('cylinder');
config.strip_edges = [-3; -2; -1; 0; 0.5; 1];
config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
config.ms2_file = config.ms2_model.filename;
for f = {'Aw_table_z', 'Aw_table', 'I_wp_yy_table', 'V_sub_table', 'CB_z_table'}
    config.(f{1}) = [0; 1];
end
[f3, s2] = sti_stage2(config, 0.5, [NaN; 300; 270; 440; 440]);
design = struct('mode', 'modular_precast', 'edges', config.strip_edges, 'vs', 0.5, ...
    't', [0.0762; 0.09; 0.08; 0.0762; 0.0762], 'z_ballast', -2.5, 'solid_modules', zeros(1, 0));
[r, fp] = sti_realised(config, design, struct('uhpc', 2500, 'air', 1.2), s2, struct('t_min', 0.0762));
schema = mwecmass.output.export_schema();
check(isequal(sort(fieldnames(r))', sort({schema.stage3_fields.name})), 'schema lists the S8 fields');
check(~any(strcmp(schema.results_top_fields, 'constructability')) && ...
    ~any(strcmp(schema.results_top_fields, 'steel_data')) && any(strcmp(schema.results_top_fields, 'stage3')), ...
    'stage3 replaces constructability and steel_data');

in = struct();
in.materials.realisation_type = 'modular_precast';
in.context = struct('T_heave_goal', 3, 'T_pitch_goal', 5, 'T_surge_goal', 60, 'T_heave_range', [2 4], ...
    'T_pitch_range', [4 6], 'T_surge_range', [50 70], 'gm_target', 0.2, 'WIS_station', '', 'data_year', '');
in.bem.water_depth = -1;
results = struct('config', config, 'stage1_2d', struct(), 'stage2_3d', struct(), 'stage3', r, ...
    'Final3D', f3, 'optimization_time', 1);
file = [tempname() '.mat'];
cleanup = onCleanup(@() delete_if(file));
mwecmass.output.export_results(results, fp, in, file);
mwecmass.output.check_export_schema(file, 'modular_precast');
L = load(file);
check(isequal(L.results.stage3.design, r.design) && strcmp(L.final_props.stage3_status, r.status), ...
    'stage3 round trip');
fprintf('modular_precast export: stage3 status %s, %d S8 fields, %d modules\n', L.results.stage3.status, ...
    numel(fieldnames(L.results.stage3)), numel(L.results.stage3.modules));

txt = evalc('mwecmass.output.Report.stage3(struct(''stage3'', r), 1);');
check(~isempty(strfind(txt, 'STAGE 3: MODULAR PRECAST')) && ~isempty(strfind(txt, 'T_pitch')), 'Report stage3 panel');
fprintf('%s', txt);

% preliminary: stage3 empty, realisation-only final_props fields typed-empty
in.materials.realisation_type = 'preliminary';
results.stage3 = [];
mwecmass.output.export_results(rmfield(results, 'stage3'), f3, in, file);
mwecmass.output.check_export_schema(file, 'preliminary');
L = load(file);
check(isempty(L.results.stage3) && isempty(L.final_props.stage3_status), 'preliminary typed-empty');

% a realising type without stage3 is rejected
in.materials.realisation_type = 'modular_precast';
mwecmass.output.export_results(rmfield(results, 'stage3'), fp, in, file);
try
    mwecmass.output.check_export_schema(file, 'modular_precast');
    ok = false;
catch err
    ok = strcmp(err.identifier, 'mwec:massSchema:unexpectedlyEmpty');
end
check(ok, 'missing stage3 rejected');
end

function delete_if(file)
if exist(file, 'file')
    delete(file);
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
    error('test_stage3_export_schema:fail', '%s', msg);
end
end
