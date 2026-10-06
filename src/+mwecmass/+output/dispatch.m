function selected = dispatch(results, final_props, config, out, context)
%DISPATCH Produce enabled post-run outputs.
%   Applies output flags and data gates, records selected artefacts, and exports results last.
    selected = {};
    validate_save_flags(out.save);

    dry_run = isfield(context, 'dry_run') && ~isempty(context.dry_run) && context.dry_run;
    x_opt           = context.x_opt;
    iteration       = context.iteration;
    ms2_file        = context.ms2_file;
    bem_cache_file  = context.bem_cache_file;
    hams_exe        = context.hams_exe;
    wp_target_edge  = context.wp_target_edge;
    in              = context.in;
    output_dir      = context.output_dir;
    type_output_dir = context.type_output_dir;

    prev_dir = pwd;
    if ~dry_run
        prev_dir = cd(type_output_dir);
    end
    try
        if offer('stage2.cross_section', true)
            mwecmass.output.figures.plot_optimised_cross_section(final_props, config, x_opt);
        end
        if offer('stage1.density_2d', true)
            mwecmass.output.figures.plot_equivalent_density_2d(final_props, config, x_opt);
        end
        if offer('stage3.final_results_log', true)
            fid_fr = mwecmass.output.open_log(out, in.materials.realisation_type, 'final_results');
            mwecmass.output.Report.final_results(final_props, config, report_fids(fid_fr));
            mwecmass.output.close_log(fid_fr);
        end

        if offer('stage2.density_3d', true)
            mwecmass.output.figures.plot_equivalent_density_3d(results.Final3D, config);
        end

        if offer('stage3.realised_density_3d', isfield(final_props, 'realised_strip_density') && ...
                ~isempty(final_props.realised_strip_density))
            mwecmass.output.figures.plot_equivalent_density_3d(final_props, config);
        end

        if offer('convergence', iteration > 1)
            try
                mwecmass.output.figures.plot_complete_convergence(results, config);
            catch
            end
        end

        if offer('stage2.summary_log', true)
            fid_s2 = mwecmass.output.open_log(out, in.materials.realisation_type, 'stage2_summary');
            report_stage2_summary(results.stage2_3d, report_fids(fid_s2));
            mwecmass.output.close_log(fid_s2);
        end

        if offer('stage3.diagnostic_panels_log', true)
            fid_dp = mwecmass.output.open_log(out, in.materials.realisation_type, 'diagnostic_panels');
            mwecmass.output.Report.all_panels(results, final_props, report_fids(fid_dp));
            mwecmass.output.close_log(fid_dp);
        end

        hams_available = exist(hams_exe, 'file') ~= 0;
        cache_file     = bem_cache_file;
        if offer({'hydrodynamics.coefficients', 'hydrodynamics.raos', ...
                  'hydrodynamics.added_mass_vs_draft'}, hams_available && exist(cache_file, 'file'))
            if isfield(results, 'Final3D')
                fp_for_mesh = results.Final3D;
            else
                fp_for_mesh = final_props;
            end
            loaded_cache = load(cache_file, 'hydro_table');
            config_final = results.config;
            mwecmass.output.figures.plot_hydrodynamics(loaded_cache.hydro_table, ...
                fp_for_mesh.vertical_shift, final_props, config_final);
        end

        if offer({'hydrodynamics.mesh_diagnostic', 'hydrodynamics.panel_normals'}, hams_available && ...
                ~isempty(ms2_file) && exist(ms2_file, 'file'))
            try
                parser_viz = mwecmass.geometry.MS2Parser.parse(ms2_file);
                if isfield(results, 'Final3D')
                    draft_viz = results.Final3D.vertical_shift;
                else
                    draft_viz = final_props.vertical_shift;
                end
                pan_opts_viz = struct('trim_wl', true, 'close_gaps', false, ...
                                      'cosine_spacing', false, 'verbose', false, ...
                                      'quarter_body', false);
                mesh_viz = mwecmass.mesh.generate(parser_viz, draft_viz, ...
                               config.mesh_Nu, config.mesh_Nv, pan_opts_viz);
                config_viz = struct('wp_target_edge', wp_target_edge);

                if out.save.hydrodynamics.mesh_diagnostic
                    mwecmass.output.figures.plot_mesh_diagnostic(mesh_viz, config_viz);
                end
                if out.save.hydrodynamics.panel_normals
                    mwecmass.output.figures.plot_panel_normals(mesh_viz, config);
                end
            catch ME_mesh
                warning('WEC:MeshPlotFailed', 'Mesh/normal plot failed: %s', ME_mesh.message);
            end
        end

        if offer('stage1.draft_landscape', strcmp(in.pid.stage1_mode, 'skip') && ...
                isfield(results.stage1_2d, 'convergence_data') && ...
                isfield(results.stage1_2d.convergence_data, 'sweep'))
            mwecmass.output.figures.plot_draft_landscape( ...
                results.stage1_2d.convergence_data.sweep, results.config);
        end

        if ~dry_run
            close all;
            cd(prev_dir);
        end
    catch ME_report
        if ~dry_run
            close all;
            cd(prev_dir);
        end
        rethrow(ME_report);
    end

    if offer('results_mat', true)
        [~, ms2_name] = fileparts(in.files.ms2_file);
        out_file = fullfile(output_dir, [ms2_name '_' in.materials.realisation_type '_results.mat']);
        mwecmass.output.export_results(results, final_props, in, out_file);
    end

    % --- Diagnostic figures (need the exported MAT file) ----------------------
    [~, ms2_name_diag] = fileparts(in.files.ms2_file);
    diag_mat = fullfile(output_dir, [ms2_name_diag '_' in.materials.realisation_type '_results.mat']);
    diag_dir = fullfile(output_dir, 'diagnostics');

    if exist(diag_mat, 'file') == 2
        repo_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
        diag_func_path = fullfile(repo_root, 'validation', 'diagnostics');
        addpath(diag_func_path);
        diag_cleanup = onCleanup(@() rmpath(diag_func_path));

        if offer('diagnostics.hull_at_draft', true)
            if ~exist(diag_dir, 'dir'), mkdir(diag_dir); end
            try
                hull_at_draft(diag_mat, diag_dir);
            catch ME_hull
                warning('mwecmass:output:HullAtDraftFailed', ...
                    'Hull-at-draft diagnostic failed: %s', ME_hull.message);
            end
        end

        if offer('diagnostics.uhpc_mass_balance', ...
                strcmp(in.materials.realisation_type, 'modular_precast'))
            if ~exist(diag_dir, 'dir'), mkdir(diag_dir); end
            try
                uhpc_mass_balance(diag_mat, diag_dir);
            catch ME_uhpc
                warning('mwecmass:output:UhpcBalanceFailed', ...
                    'UHPC mass-balance diagnostic failed: %s', ME_uhpc.message);
            end
        end

        if offer('diagnostics.stage_animations', true)
            if ~exist(diag_dir, 'dir'), mkdir(diag_dir); end
            try
                stage_animations(diag_mat, diag_dir);
            catch ME_anim
                warning('mwecmass:output:AnimationsFailed', ...
                    'Stage animations failed: %s', ME_anim.message);
            end
        end

        if ~dry_run
            close all;
        end
    end

    function tf = offer(flag_names, gate)
        if ischar(flag_names)
            flag_names = {flag_names};
        end
        tf = false;
        if ~gate
            return;
        end
        for k = 1:numel(flag_names)
            if save_flag(flag_names{k})
                selected{end+1} = flag_names{k}; %#ok<AGROW> one entry per artefact, at most three
                tf = true;
            end
        end
        tf = tf && ~dry_run;
    end

    function enabled = save_flag(flag_name)
        names = strsplit(flag_name, '.');
        enabled = out.save;
        for n = 1:numel(names)
            enabled = enabled.(names{n});
        end
    end

    function fids = report_fids(fid)
        if out.console_echo
            fids = [1 fid];
        else
            fids = fid;
        end
    end
end

function validate_save_flags(save)
    flags = [ ...
        save.stage1.density_2d; save.stage1.draft_landscape; ...
        save.stage2.density_3d; save.stage2.cross_section; save.stage2.summary_log; ...
        save.convergence; ...
        save.stage3.steel_solve; save.stage3.steel_solve_log; save.stage3.precast_midplane; ...
        save.stage3.precast_strips; save.stage3.realised_density_3d; save.stage3.final_results_log; ...
        save.stage3.diagnostic_panels_log; ...
        save.hydrodynamics.coefficients; save.hydrodynamics.raos; ...
        save.hydrodynamics.added_mass_vs_draft; save.hydrodynamics.mesh_diagnostic; ...
        save.hydrodynamics.panel_normals; save.hydrodynamics.cache_rewrite; save.results_mat; ...
        save.diagnostics.hull_at_draft; save.diagnostics.uhpc_mass_balance; ...
        save.diagnostics.stage_animations; ...
    ];
    if ~islogical(flags) || ~isvector(flags)
        error('mwecmass:output:InvalidSaveFlags', ...
              'Every out.save leaf must be a logical scalar.');
    end
end

function report_stage2_summary(stage2_data, fids)
    if nargin < 2, fids = 1; end
    qm = stage2_data.quality_metrics;
    mwecmass.output.emit(fids, '\n=== STAGE 2 SUMMARY ===\n\n');
    mwecmass.output.emit(fids, '  Convergence: %d\n', stage2_data.converged);
    mwecmass.output.emit(fids, '\n  Constraints:\n');
    mwecmass.output.emit(fids, '    Monotonic density: %d\n', qm.monotonic);
    mwecmass.output.emit(fids, '    Mass balance: %d (error: %.2e kg)\n', qm.mass_balance, qm.mass_balance_error_kg);
    mwecmass.output.emit(fids, '    GM constraint: %d (margin: %.3f m)\n', qm.GM_satisfied, qm.GM_margin);
    mwecmass.output.emit(fids, '\n  Solver:\n');
    mwecmass.output.emit(fids, '    Exit flag: %d\n', qm.exitflag);
    mwecmass.output.emit(fids, '    Optimal: %d\n', qm.fmincon_optimal);
    mwecmass.output.emit(fids, '    Constraint violation: %.2e\n', qm.constrviolation);
    mwecmass.output.emit(fids, '    First-order optimality: %.2e\n\n', qm.firstorderopt);
end
