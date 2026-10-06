function schema = export_schema()
%EXPORT_SCHEMA Return shared results-export field/class/shape metadata.
%   This is the single source for export and validation inventories.
    schema = struct();
    schema.schema_version = '1.0';
    schema.types = {'preliminary', 'thin_shell', 'modular_precast'};
    schema.results_top_fields = {'config', 'stage1_2d', 'stage2_3d', 'constructability', ...
        'steel_data', 'Final3D', 'optimization_time', 'schema_version', 'realisation_type', 'context'};
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
    fp = local_add(fp, 'realised_strips',            'struct', [1 1], true);
    fp = local_add(fp, 'cross_section',              'double', [], false);   % [398x2] on the probed file; row count is geometry-dependent
    fp = local_add(fp, 'densities_at_nodes',         'double', [], false);   % [5x1] on the probed file; row count = num_ballast_sections
    fp = local_add(fp, 'realised_strip_density',     'double', [], true);
    fp = local_add(fp, 'realised_strip_edges',       'double', [], true);
    fp = local_add(fp, 'components',                 'struct', [], false);   % struct array; size/field-count is realisation-type-dependent, handled separately
    schema.final_props_fields = fp;
    schema.final_props_field_count = numel(fp);   % 40, per the spec's own re-derived count

    ct = struct('name', {}, 'class', {}, 'size', {});
    ct = local_add3(ct, 't_steel','double',[1 1]);            ct = local_add3(ct, 'z_fill','double',[1 1]);
    ct = local_add3(ct, 'draft','double',[1 1]);               ct = local_add3(ct, 'vertical_shift','double',[1 1]);
    ct = local_add3(ct, 'draft_optimiser','double',[1 1]);     ct = local_add3(ct, 'vs_optimiser','double',[1 1]);
    ct = local_add3(ct, 'rho_steel','double',[1 1]);           ct = local_add3(ct, 'rho_air','double',[1 1]);
    ct = local_add3(ct, 't_max','double',[1 1]);               ct = local_add3(ct, 't_min','double',[1 1]);
    ct = local_add3(ct, 't_min_active','logical',[1 1]);       ct = local_add3(ct, 'V_steel','double',[1 1]);
    ct = local_add3(ct, 'V_air','double',[1 1]);               ct = local_add3(ct, 'V_hull','double',[1 1]);
    ct = local_add3(ct, 'M_steel','double',[1 1]);             ct = local_add3(ct, 'M_air','double',[1 1]);
    ct = local_add3(ct, 'M_total','double',[1 1]);             ct = local_add3(ct, 'z_cg_steel','double',[1 1]);
    ct = local_add3(ct, 'z_cg_air','double',[1 1]);            ct = local_add3(ct, 'CG_z_body','double',[1 1]);
    ct = local_add3(ct, 'CG_z_world','double',[1 1]);          ct = local_add3(ct, 'Iyy_total_origin','double',[1 1]);
    ct = local_add3(ct, 'Iyy_about_cg','double',[1 1]);        ct = local_add3(ct, 'Ixx_total_origin','double',[1 1]);
    ct = local_add3(ct, 'Ixx_about_cg','double',[1 1]);        ct = local_add3(ct, 'Izz_total_origin','double',[1 1]);
    ct = local_add3(ct, 'Izz_about_cg','double',[1 1]);        ct = local_add3(ct, 'V_sub','double',[1 1]);
    ct = local_add3(ct, 'Aw','double',[1 1]);                  ct = local_add3(ct, 'I_wp_yy','double',[1 1]);
    ct = local_add3(ct, 'CB_z_world','double',[1 1]);          ct = local_add3(ct, 'KM_world','double',[1 1]);
    ct = local_add3(ct, 'GM_realised','double',[1 1]);         ct = local_add3(ct, 'T_heave_realised','double',[1 1]);
    ct = local_add3(ct, 'T_pitch_realised','double',[1 1]);    ct = local_add3(ct, 'K33_hydro','double',[1 1]);
    ct = local_add3(ct, 'K55_hydro','double',[1 1]);           ct = local_add3(ct, 'A11','double',[1 1]);
    ct = local_add3(ct, 'A33','double',[1 1]);                 ct = local_add3(ct, 'A55','double',[1 1]);
    ct = local_add3(ct, 'targets','struct',[1 1]);             ct = local_add3(ct, 'residuals','struct',[1 1]);
    ct = local_add3(ct, 'mass_balance_error_pct','double',[1 1]); ct = local_add3(ct, 'phi_star','double',[1 1]);
    ct = local_add3(ct, 'feasible','logical',[1 1]);           ct = local_add3(ct, 'exitflag','double',[1 1]);
    ct = local_add3(ct, 'solver','char',[]);                   ct = local_add3(ct, 'z_grid','double',[]);
    ct = local_add3(ct, 'A_outer_grid','double',[]);           ct = local_add3(ct, 'A_inner_grid','double',[]);
    ct = local_add3(ct, 'A_jacket_grid','double',[]);          ct = local_add3(ct, 'Iyy_outer_grid','double',[]);
    ct = local_add3(ct, 'Iyy_inner_grid','double',[]);         ct = local_add3(ct, 't_offset_strip','double',[]);
    ct = local_add3(ct, 'is_solid_strip','logical',[]);        ct = local_add3(ct, 'strip_edges','double',[]);
    ct = local_add3(ct, 'wall_strip_idx','double',[1 1]);      ct = local_add3(ct, 'iter_history','struct',[1 1]);
    ct = local_add3(ct, 'elapsed_seconds','double',[1 1]);     ct = local_add3(ct, 'fill_method','char',[]); % Task S7b: renamed from realisation_mode (see final_props_fields' 'fill_method' comment above): 'uhpc_fill', a fill-method flag, not a realisation type
    ct = local_add3(ct, 'mode','char',[]);                     ct = local_add3(ct, 'Z_max','double',[1 1]);   % this 'mode' = 'constructable_hull' (extract_strip_geometry.m:361), a fixed solver-identity tag, unrelated to the Stage-4 realisation-type selector and never read anywhere in this repository (grep verified) -- not touched by Task S7
    ct = local_add3(ct, 'Z_min','double',[1 1]);               ct = local_add3(ct, 'wall_z_bottom','double',[1 1]);
    ct = local_add3(ct, 'wall_z_top','double',[1 1]);          ct = local_add3(ct, 'rho_hull','double',[1 1]);
    ct = local_add3(ct, 'rho_fill','double',[1 1]);            ct = local_add3(ct, 'wall_height','double',[1 1]);
    ct = local_add3(ct, 't_UHPC','double',[1 1]);              ct = local_add3(ct, 'strip_z_lo','double',[]);
    ct = local_add3(ct, 'strip_z_hi','double',[]);             ct = local_add3(ct, 'strip_rho_eff','double',[]);
    ct = local_add3(ct, 'strip_scale_factor','double',[]);     ct = local_add3(ct, 'strip_V_total','double',[]);
    ct = local_add3(ct, 'strip_V_UHPC','double',[]);           ct = local_add3(ct, 'strip_V_void','double',[]);
    ct = local_add3(ct, 'strip_mass_UHPC','double',[]);        ct = local_add3(ct, 'strip_mass_void','double',[]);
    ct = local_add3(ct, 'strip_mass_total','double',[]);       ct = local_add3(ct, 'strip_t_min_actual','double',[]);
    ct = local_add3(ct, 'strip_r_min','double',[]);            ct = local_add3(ct, 'strip_z_cg','double',[]);
    ct = local_add3(ct, 'strip_is_wall','logical',[]);         ct = local_add3(ct, 'strip_is_solid','logical',[]);
    ct = local_add3(ct, 'strip_is_feasible','logical',[]);     ct = local_add3(ct, 'strip_Iyy_UHPC','double',[]);
    ct = local_add3(ct, 'strip_Iyy_void','double',[]);         ct = local_add3(ct, 'contours_outer','cell',[]);
    ct = local_add3(ct, 'contours_inner','cell',[]);           ct = local_add3(ct, 'total_mass','double',[1 1]);
    ct = local_add3(ct, 'total_V_UHPC','double',[1 1]);        ct = local_add3(ct, 'total_V_void','double',[1 1]);
    ct = local_add3(ct, 'total_V_hull','double',[1 1]);        ct = local_add3(ct, 'UHPC_volume_fraction','double',[1 1]);
    ct = local_add3(ct, 'feasibility','struct',[1 1]);
    schema.constructability_fields = ct;
    schema.constructability_field_count = numel(ct);   % 95, per the spec's own re-derived count

    sd = struct('name', {}, 'class', {}, 'size', {});
    sd = local_add3(sd, 't_steel',              'double',  [1 1]); sd = local_add3(sd, 'z_fill',               'double',  [1 1]);
    sd = local_add3(sd, 'draft',                'double',  [1 1]); sd = local_add3(sd, 'vertical_shift',       'double',  [1 1]);
    sd = local_add3(sd, 'draft_optimiser',      'double',  [1 1]); sd = local_add3(sd, 'vs_optimiser',         'double',  [1 1]);
    sd = local_add3(sd, 'rho_steel',            'double',  [1 1]); sd = local_add3(sd, 'rho_shell',            'double',  [1 1]);   % two-density model
    sd = local_add3(sd, 'rho_fill',             'double',  [1 1]); sd = local_add3(sd, 'rho_air',              'double',  [1 1]);   % two-density model
    sd = local_add3(sd, 't_max',                'double',  [1 1]); sd = local_add3(sd, 't_min',                'double',  [1 1]);
    sd = local_add3(sd, 't_min_active',         'logical', [1 1]); sd = local_add3(sd, 'V_steel',              'double',  [1 1]);
    sd = local_add3(sd, 'V_air',                'double',  [1 1]); sd = local_add3(sd, 'V_hull',               'double',  [1 1]);
    sd = local_add3(sd, 'V_shell',              'double',  [1 1]); sd = local_add3(sd, 'V_fill',               'double',  [1 1]);   % two-density model
    sd = local_add3(sd, 'M_steel',              'double',  [1 1]); sd = local_add3(sd, 'M_air',                'double',  [1 1]);
    sd = local_add3(sd, 'M_shell',              'double',  [1 1]); sd = local_add3(sd, 'M_fill',               'double',  [1 1]);   % two-density model
    sd = local_add3(sd, 'M_total',              'double',  [1 1]); sd = local_add3(sd, 'z_cg_steel',           'double',  [1 1]);
    sd = local_add3(sd, 'z_cg_fill',            'double',  [1 1]); sd = local_add3(sd, 'z_cg_shell',           'double',  [1 1]);   % two-density model
    sd = local_add3(sd, 'z_cg_air',             'double',  [1 1]); sd = local_add3(sd, 'CG_z_body',            'double',  [1 1]);
    sd = local_add3(sd, 'CG_z_world',           'double',  [1 1]); sd = local_add3(sd, 'Iyy_total_origin',     'double',  [1 1]);
    sd = local_add3(sd, 'Iyy_about_cg',         'double',  [1 1]); sd = local_add3(sd, 'Ixx_total_origin',     'double',  [1 1]);
    sd = local_add3(sd, 'Ixx_about_cg',         'double',  [1 1]); sd = local_add3(sd, 'Izz_total_origin',     'double',  [1 1]);
    sd = local_add3(sd, 'Izz_about_cg',         'double',  [1 1]); sd = local_add3(sd, 'V_sub',                'double',  [1 1]);
    sd = local_add3(sd, 'Aw',                   'double',  [1 1]); sd = local_add3(sd, 'I_wp_yy',              'double',  [1 1]);
    sd = local_add3(sd, 'CB_z_world',           'double',  [1 1]); sd = local_add3(sd, 'KM_world',             'double',  [1 1]);
    sd = local_add3(sd, 'GM_realised',          'double',  [1 1]); sd = local_add3(sd, 'T_heave_realised',     'double',  [1 1]);
    sd = local_add3(sd, 'T_pitch_realised',     'double',  [1 1]); sd = local_add3(sd, 'K33_hydro',            'double',  [1 1]);
    sd = local_add3(sd, 'K55_hydro',            'double',  [1 1]); sd = local_add3(sd, 'A11',                  'double',  [1 1]);
    sd = local_add3(sd, 'A33',                  'double',  [1 1]); sd = local_add3(sd, 'A55',                  'double',  [1 1]);
    sd = local_add3(sd, 'targets',              'struct',  [1 1]); sd = local_add3(sd, 'residuals',            'struct',  [1 1]);
    sd = local_add3(sd, 'mass_balance_error_pct','double', [1 1]); sd = local_add3(sd, 'phi_star',             'double',  [1 1]);
    sd = local_add3(sd, 'feasible',             'logical', [1 1]); sd = local_add3(sd, 'exitflag',             'double',  [1 1]);
    sd = local_add3(sd, 'solver',               'char',    []);   sd = local_add3(sd, 'z_grid',                'double',  []);
    sd = local_add3(sd, 'A_outer_grid',         'double',  []);   sd = local_add3(sd, 'A_inner_grid',          'double',  []);
    sd = local_add3(sd, 'A_jacket_grid',        'double',  []);   sd = local_add3(sd, 'Iyy_outer_grid',        'double',  []);
    sd = local_add3(sd, 'Iyy_inner_grid',       'double',  []);   sd = local_add3(sd, 'strip_rho_eff',         'double',  []);
    sd = local_add3(sd, 'strip_edges',          'double',  []);   sd = local_add3(sd, 'strip_V_env',           'double',  []);
    sd = local_add3(sd, 'strip_V_solid',        'double',  []);   sd = local_add3(sd, 'strip_V_void',          'double',  []);
    sd = local_add3(sd, 'strip_V_fill',         'double',  []);   sd = local_add3(sd, 'strip_V_shell',         'double',  []);      % two-density model
    sd = local_add3(sd, 'fill_method',          'char',    []);   sd = local_add3(sd, 'elapsed_seconds',       'double',  [1 1]); % Task S7b: renamed from realisation_mode (see final_props_fields' comment above): 'steel_fill', a fill-method flag, not a realisation type
    schema.steel_data_fields = sd;
    schema.steel_data_field_count = numel(sd);   % 70 (matches solve.m's own field count exactly, per the block comment above)

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
