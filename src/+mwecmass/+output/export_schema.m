function schema = export_schema()
%EXPORT_SCHEMA Return shared results-export field/class/shape metadata.
%   This is the single source for export and validation inventories.
    schema = struct();
    schema.schema_version = '1.0';
    schema.types = {'preliminary', 'thin_shell', 'modular_precast'};
    schema.results_top_fields = {'config', 'stage1_2d', 'stage2_3d', 'stage3', 'Final3D', ...
        'optimization_time', 'schema_version', 'realisation_type', 'context'};
    schema.config_required_subfields = {'Aw_table_z', 'Aw_table', 'I_wp_yy_table', ...
        'V_sub_table', 'CB_z_table'};
    schema.config_required_char_subfields = {'ms2_deck_sha256'};

    cx = struct('name', {}, 'class', {}, 'size', {});
    cx = local_add3(cx, 'T_heave_goal',  'double', [1 1]);
    cx = local_add3(cx, 'T_pitch_goal',  'double', [1 1]);
    cx = local_add3(cx, 'T_surge_goal',  'double', [1 1]);
    cx = local_add3(cx, 'T_heave_range', 'double', [1 2]);
    cx = local_add3(cx, 'T_pitch_range', 'double', [1 2]);
    cx = local_add3(cx, 'T_surge_range', 'double', [1 2]);
    cx = local_add3(cx, 'gm_target',     'double', [1 1]);
    cx = local_add3(cx, 'WIS_station',   'char',   []);
    cx = local_add3(cx, 'data_year',     'char',   []);
    cx = local_add3(cx, 'water_depth',   'double', [1 1]);
    schema.context_fields = cx;
    schema.context_field_count = numel(cx);   % 10

    fp = struct('name', {}, 'class', {}, 'size', {}, 'realisation_only', {});
    fp = local_add(fp, 'vertical_shift',           'double', [1 1], false);
    fp = local_add(fp, 'draft',                    'double', [1 1], false);
    fp = local_add(fp, 'Aw',                       'double', [1 1], false);
    fp = local_add(fp, 'I_wp_yy',                  'double', [1 1], false);
    fp = local_add(fp, 'I_wp_xx',                  'double', [1 1], false);
    fp = local_add(fp, 'V_sub',                    'double', [1 1], false);
    fp = local_add(fp, 'CB',                       'double', [1 3], false);
    fp = local_add(fp, 'A_sub',                    'double', [1 1], false);
    fp = local_add(fp, 'mass_buoyant_force',       'double', [1 1], false);
    fp = local_add(fp, 'KM',                       'double', [1 1], false);
    fp = local_add(fp, 'mass_total',                'double', [1 1], false);
    fp = local_add(fp, 'CG_total',                  'double', [1 3], false);
    fp = local_add(fp, 'Ixx',                       'double', [1 1], false);
    fp = local_add(fp, 'Iyy',                       'double', [1 1], false);
    fp = local_add(fp, 'Izz',                       'double', [1 1], false);
    fp = local_add(fp, 'Inertia_Tensor',            'double', [3 3], false);
    fp = local_add(fp, 'GM_L',                      'double', [1 1], false);
    fp = local_add(fp, 'mass_discrepancy',          'double', [1 1], false);
    fp = local_add(fp, 'K_hydro',                   'double', [3 3], false);
    fp = local_add(fp, 'A11',                       'double', [1 1], false);
    fp = local_add(fp, 'A33',                       'double', [1 1], false);
    fp = local_add(fp, 'A55',                       'double', [1 1], false);
    fp = local_add(fp, 'A_full',                    'double', [3 3], false);
    fp = local_add(fp, 'B_full',                    'double', [3 3], false);
    fp = local_add(fp, 'K_pto',                     'double', [3 3], false);
    fp = local_add(fp, 'K_total',                   'double', [3 3], false);
    fp = local_add(fp, 'periods',                   'struct', [1 1], false);
    fp = local_add(fp, 'coupled_periods',           'double', [3 1], false);
    fp = local_add(fp, 'coupled_modes',              'double', [3 3], false);
    fp = local_add(fp, 'participation_factors',      'double', [3 3], false);
    fp = local_add(fp, 'MassMatrix_CG',              'double', [6 6], false);
    fp = local_add(fp, 'MassMatrix_Origin',          'double', [6 6], false);
    fp = local_add(fp, 'fill_method',                'char',   [], true);   % renamed from realisation_mode
    fp = local_add(fp, 'density_profile_source',     'char',   [], true);
    fp = local_add(fp, 'cross_section',              'double', [], false);   % [398x2] on the probed file; row count is geometry-dependent
    fp = local_add(fp, 'densities_at_nodes',         'double', [], false);   % [5x1] on the probed file; row count = num_ballast_sections
    fp = local_add(fp, 'realised_strip_density',     'double', [], true);
    fp = local_add(fp, 'realised_strip_edges',       'double', [], true);
    fp = local_add(fp, 'components',                 'struct', [], false);   % struct array; size/field-count is realisation-type-dependent, handled separately
    fp = local_add(fp, 'stage3_status',              'char',   [], true);    % 'accepted' | 'failed'
    fp = local_add(fp, 'stage3_check',               'struct', [1 1], true); % check_against_stage2 output
    schema.final_props_fields = fp;
    schema.final_props_field_count = numel(fp);   % 41

    % results.stage3 (contract S8), both realising types; an empty class is not checked (k_star,
    % V_uhpc_target and fit are [] where they do not apply).
    st = struct('name', {}, 'class', {}, 'size', {});
    st = local_add3(st, 'mode',          'char',   []);
    st = local_add3(st, 'hull_name',     'char',   []);
    st = local_add3(st, 'status',        'char',   []);
    st = local_add3(st, 'reason',        'char',   []);
    st = local_add3(st, 'escalation',    'char',   []);
    st = local_add3(st, 'vs',            'double', [1 1]);
    st = local_add3(st, 'draft',         'double', [1 1]);
    st = local_add3(st, 'stage2',        'struct', [1 1]);
    st = local_add3(st, 'rho',           'struct', [1 1]);
    st = local_add3(st, 'design',        'struct', [1 1]);
    st = local_add3(st, 'k_star',        '',       []);
    st = local_add3(st, 'V_uhpc_target', '',       []);
    st = local_add3(st, 'modules',       'struct', []);
    st = local_add3(st, 'props',         'struct', [1 1]);
    st = local_add3(st, 'check',         'struct', [1 1]);
    st = local_add3(st, 'solver',        'struct', []);
    st = local_add3(st, 'fit',           '',       []);
    st = local_add3(st, 'body',          'struct', [1 1]);
    st = local_add3(st, 'step_files',    'struct', []);
    schema.stage3_fields = st;
    schema.stage3_field_count = numel(st);   % 19

    fs = struct('name', {}, 'class', {}, 'size', {});
    fs = local_add3(fs, 'CB',                'double', [1 3]);
    fs = local_add3(fs, 'CG_total',          'double', [1 3]);
    fs = local_add3(fs, 'Inertia_Tensor',    'double', [3 3]);
    fs = local_add3(fs, 'K_hydro',           'double', [3 3]);
    fs = local_add3(fs, 'A_full',            'double', [3 3]);
    fs = local_add3(fs, 'MassMatrix_CG',     'double', [6 6]);
    fs = local_add3(fs, 'MassMatrix_Origin', 'double', [6 6]);
    fs = local_add3(fs, 'components',        'struct', [5 1]);
    schema.fixed_shape_checks = fs;
end

function s = local_add(s, name, class_, size_, realisation_only)
    k = numel(s) + 1;
    s(k).name = name; s(k).class = class_; s(k).size = size_; s(k).realisation_only = realisation_only;
end

function s = local_add3(s, name, class_, size_)
    k = numel(s) + 1;
    s(k).name = name; s(k).class = class_; s(k).size = size_;
end
