function hydro_table = run(in, repo_root)
%RUN  Build the HAMS-MREL hydro_table from geometry and draft sweeps.
% Inputs are the author-set in struct and absolute repo_root; the route resolves the MS2 deck,
% solver, workspace, and cache paths from them.  It sizes the panel mesh with
% mwecmass.mesh.panel_grid_counts, builds geometry-only config via build_config, then calls
% generate_draft_nodes for coarse plus adaptive refinement.  hydro_table is the 23-field BEM
% contract, with signed vertical_shift [m], mesh counts [-], and wp_target_edge [m], and is
% returned unsorted; mwecmass.bem.load_hydro_cache sorts loaded caches.

    ms2_file        = fullfile(repo_root, 'Input', in.files.ms2_file);
    hams_dir        = fullfile(repo_root, in.bem.hams_dir);
    hams_exe        = fullfile(repo_root, in.bem.hams_exe);
    % Write target of this route's own cache, under Output/hams_mrel/; never the WAMIT cache
    % under Input/, which this suite only reads.
    hams_cache_file = fullfile(repo_root, in.bem.hams_dir, in.bem.hams_cache_file);

    % Parse geometry and size the panel mesh.
    parser_pre = mwecmass.geometry.MS2Parser.parse(ms2_file);
    [Nu, Nv, sizing_report] = mwecmass.mesh.panel_grid_counts(parser_pre, in.geometry.panel_size);
    mesh_sizing = struct('mesh_Nu', Nu, ...                      % panel count, u direction [-]
                         'mesh_Nv', Nv, ...                      % panel count, v direction [-]
                         'wp_target_edge', sizing_report.wp_target_edge);   % [m]

    % Build the geometry-only configuration.
    geo_config = mwecmass.driver.build_config(in, [], mesh_sizing);   % empty 2nd arg = geometry-only

    % Run the coarse and adaptive HAMS sweep.
    hydro_table = mwecmass.bem.hams_mrel.generate_draft_nodes(in, mesh_sizing, geo_config, ...
        ms2_file, hams_dir, hams_exe, hams_cache_file);
end
