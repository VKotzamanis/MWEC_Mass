function run(in, out)
%RUN Orchestrate input loading, BEM source selection, configuration, optimisation, realisation, and outputs.
% Syntax: run(in,out), with either struct optional; omitted values come from WEC_User_Input and
% WEC_Output_Options. Inputs are the author input and output-selection structs; no interactive input.
% The pipeline resolves paths relative to the repository root and dispatches realisation by type.

    repo_root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
    % Kept for a caller that invokes this function before src/ is on the path; harmless when src/
    % is already there.
    addpath(fullfile(repo_root, 'src'));

    % Read the two root files.
    if nargin < 1 || isempty(in), in = mwecmass.driver.load_inputs(); end
    if nargin < 2 || isempty(out), out = WEC_Output_Options(); end

    % Validate realisation_type.
    valid_types = {'preliminary', 'thin_shell', 'modular_precast'};
    if ~ismember(in.materials.realisation_type, valid_types)
        error('mwecmass:driver:badRealisationType', ...
            'in.materials.realisation_type = ''%s'' is not one of preliminary/thin_shell/modular_precast.', ...
            in.materials.realisation_type);
    end

    % Resolve input files.
    % Resolve the two input files under Input/ and the solver executable under the repository root.
    ms2_file       = fullfile(repo_root, 'Input', in.files.ms2_file);       % hull geometry deck
    bem_cache_file = fullfile(repo_root, 'Input', in.bem.bem_cache_file);   % BEM cache of record
    hams_exe       = fullfile(repo_root, in.bem.hams_exe);                  % HAMS-MREL binary

    % Resolve output directories.
    % Resolve output directories before the BEM source step so reporting and export share them.
    output_dir = fullfile(repo_root, 'Output');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    type_output_dir = fullfile(output_dir, in.materials.realisation_type);
    if ~exist(type_output_dir, 'dir')
        mkdir(type_output_dir);
    end

    % Resolve the BEM source.
    % mesh_sizing carries the BEM panel-grid counts of the run: empty until the cache (or, on the
    % HAMS-MREL branch, the sizing pass inside that route) supplies them. wp_target_edge is the
    % waterplane-lid panel edge target [m], author-set unless the cache records the one it was
    % built with. Neither is read by the optimiser -- only the .pnl export and the panel-mesh
    % diagnostic use them.
    mesh_sizing    = struct();
    wp_target_edge = in.geometry.wp_target_edge;   % [m]
    if in.bem.run_HAMS_MREL
        % Frozen behind the flag; +bem/+hams_mrel owns this branch and resolves its own
        % workspace paths from the repository root it is given.
        hydro_table = mwecmass.bem.hams_mrel.run(in, repo_root);
    else
        hydro_table = mwecmass.bem.load_hydro_cache(bem_cache_file);

        % The cache records the panel-grid counts and the waterplane edge length it was built
        % with; without restoring them, build_config's own mesh_sizing guard would default
        % config.mesh_Nu/config.mesh_Nv to [] instead of the cached values.
        [mesh_sizing.mesh_Nu, mesh_sizing.mesh_Nv, wp_edge_from_cache] = ...
            mwecmass.bem.wamit.restore_mesh_sizing_from_cache(hydro_table, in);
        if ~isempty(wp_edge_from_cache)
            wp_target_edge = wp_edge_from_cache;
        end
    end

    % Run the pipeline.
    % Wrapped in the pipeline's outer try/catch: mwecmass.output.Report.exception prints the
    % error and the likely missing input before it is rethrown. The BEM source's hydro and mesh-sizing
    % loading above stays outside this try.
    try
        % Build configuration.
        config = mwecmass.driver.build_config(in, hydro_table, mesh_sizing);

        % The output options travel with the configuration, so every stage that writes a figure,
        % a console panel or a file can reach them from the one struct it already receives.
        config.output = out;

        % Optimise.
        [opt_results, x_opt, iteration] = mwecmass.optim.run(config);

        % Realise.
        realise_fn = str2func(['mwecmass.realise.' in.materials.realisation_type '.run']);
        [results, final_props] = realise_fn(config, x_opt, opt_results);

        % Write outputs.
        % Run outputs include post-realisation figures, console reports and the
        % exported results file -- is produced by mwecmass.output.dispatch, which keeps the call
        % order and the gates of the two post-run blocks and offers each call only when its
        % own flag in out.save is true. It also manages the working-directory change into
        % Output/<realisation_type>/ that the bare relative image filenames need, and closes the
        % figures this run has open before the export. The context struct carries the values those
        % gates and the export read that are not already on results or config: optimisation's design
        % vector and iteration count, the three resolved paths the machine-dependent gates test
        % (the solver binary, the cached coefficients, the hull deck), the waterplane edge length
        % the diagnostic mesh is drawn at, the author-set input struct the export writes into the
        % file, and the two resolved output directories.
        dispatch_context = struct();
        dispatch_context.x_opt           = x_opt;            % [1x(1+N)] [m; kg/m^3 ...]
        dispatch_context.iteration       = iteration;        % 2-D density search's iteration count [-]
        dispatch_context.ms2_file        = ms2_file;         % hull geometry deck
        dispatch_context.bem_cache_file  = bem_cache_file;   % cached BEM coefficients
        dispatch_context.hams_exe        = hams_exe;         % HAMS-MREL binary (presence gate)
        dispatch_context.wp_target_edge  = wp_target_edge;   % [m] waterplane-lid panel edge
        dispatch_context.in              = in;               % author-set input struct
        dispatch_context.output_dir      = output_dir;       % <repository>/Output
        dispatch_context.type_output_dir = type_output_dir;  % <repository>/Output/<type>
        mwecmass.output.dispatch(results, final_props, config, out, dispatch_context);

    catch ME
        mwecmass.output.Report.exception(ME, ms2_file, bem_cache_file);
        rethrow(ME);
    end
end
