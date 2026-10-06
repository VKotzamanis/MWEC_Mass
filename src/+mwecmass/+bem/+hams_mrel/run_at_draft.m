function hams_data = run_at_draft(config, vertical_shift, hams_dir, hams_exe, options)
%RUN_AT_DRAFT  Run HAMS at one vertical shift and parse hydrodynamic output.

    if nargin < 5; options = struct(); end
    if ~isfield(options, 'Nu');             options.Nu = config.mesh_Nu; end
    if ~isfield(options, 'Nv');             options.Nv = config.mesh_Nv; end
    if ~isfield(options, 'wp_target_edge'); options.wp_target_edge = 0.4; end
    % Averaging band stays inside period grid HAMS solves; config.period_min/max is single source of truth.
    if ~isfield(options, 'T_min')
        if isfield(config, 'period_min') && ~isempty(config.period_min)
            options.T_min = config.period_min;
        else
            options.T_min = 2.0;
        end
    end
    if ~isfield(options, 'T_max')
        if isfield(config, 'period_max') && ~isempty(config.period_max)
            options.T_max = config.period_max;
        else
            options.T_max = 20.0;
        end
    end
    if ~isfield(options, 'verbose');        options.verbose = true; end
    % Full mesh is the default; symmetry reduction is opt-in.
    % Use full mesh (ISX=0, ISY=0) — no HAMS mirroring.
    if ~isfield(options, 'half_body');      options.half_body = false; end

    parser = config.ms2_model;
    draft  = vertical_shift;   % panelizer convention
    z_wl   = -vertical_shift;  % body-frame waterline

    if options.verbose
        fprintf('  HAMS run: vertical_shift=%+.4f m\n', vertical_shift);
    end

    % 1. Generate mesh (half-body for HAMS-MREL)
    pan_opts = struct('trim_wl', true, 'close_gaps', false, ...
                      'cosine_spacing', false, 'verbose', false, ...
                      'quarter_body', false, 'half_body', options.half_body);
    mesh = mwecmass.mesh.generate(parser, draft, ...
               options.Nu, options.Nv, pan_opts);

    % 2. Write mesh files
    setup_hams_directory(hams_dir);
    write_panelizer_hull_pnl(mesh, ...
        fullfile(hams_dir, 'Input', 'HullMesh.pnl'));

    % 2b. Build the waterplane lid from the written hull waterline; lock its
    % boundary nodes while using options.wp_target_edge for interior density.
    boundary_xy = mwecmass.mesh.hull_waterline_polygon(mesh);

    if ~isempty(boundary_xy) && size(boundary_xy, 1) >= 3
        [wp_nodes, wp_panels, wp_nverts] = mwecmass.mesh.waterplane_mesh_unstructured( ...
            boundary_xy, options.wp_target_edge, mesh.x_sym, mesh.y_sym, 0, true);
        if ~isempty(wp_nodes)
            wp_nodes(:, 3) = 0;   % force world z = 0 exactly
        end
        if options.verbose
            Aw_hull = abs(polyarea(boundary_xy(:,1), boundary_xy(:,2)));
            fprintf('  WP hull-match: %d boundary pts, %d nodes, %d panels, Aw=%.4f m²\n', ...
                size(boundary_xy,1), size(wp_nodes,1), numel(wp_nverts), Aw_hull);
        end
    else
        warning('mwecmass:hams_mrel:WPFailed', ...
            'hull_waterline_polygon empty at vs=%+.4f m. IRSP disabled.', vertical_shift);
        wp_nodes  = zeros(0,3);
        wp_panels = zeros(0,4);
        wp_nverts = zeros(0,1);
    end

    mwecmass.bem.hams_mrel.HamsWriter.pnl_file( ...
        fullfile(hams_dir, 'Input', 'WaterPlaneMesh.pnl'), ...
        wp_nodes, wp_panels, wp_nverts, mesh.x_sym, mesh.y_sym);

    % 3. Mass + restoring
    [CG_global, M_6x6, C_6x6, m_k, sub_k] = ...
        mwecmass.hydrostatics.hydrostatic_inputs(parser, z_wl, config);

    if options.verbose
        fprintf('    V_sub=%.4f m3, mass=%.1f kg, CG_z=%.4f m\n', ...
                sub_k.V_sub, m_k, CG_global(3));
    end

    % 4. Write HAMS input files
    hams_params = mwecmass.bem.hams_mrel.default_hams_params(config);
    % XR = [0,0,0]: HAMS outputs A, B, Fe at the global origin.
    % mwecmass.driver.build_config / rebuild_config_hydro apply the
    % single congruence transform origin → CG in post-processing.
    % retransform_at_actual_cg re-does origin → actual CG for the
    % converged draft in trained mode.
    hams_params.ref_body_center = [0, 0, 0];
    mwecmass.bem.hams_mrel.HamsWriter.control_file( ...
        fullfile(hams_dir, 'Input', 'ControlFile.in'), hams_params);
    mwecmass.bem.hams_mrel.HamsWriter.hydrostatic_file( ...
        fullfile(hams_dir, 'Input', 'Hydrostatic.in'), ...
        CG_global, M_6x6, zeros(6), zeros(6), C_6x6, zeros(6));

    % 5. Run HAMS
    [status, result] = mwecmass.bem.hams_mrel.run_solver(hams_exe, hams_dir);

    if status ~= 0
        warning('mwecmass:hams_mrel:SingleRunFailed', ...
                'HAMS failed at vs=%+.4f: %s', vertical_shift, result);
        hams_data.added_mass_inf = zeros(6);
        hams_data.radiation_damping_band_avg = zeros(6);
        hams_data.added_mass_omega = [];
        hams_data.radiation_damping_omega = [];
        hams_data.exciting_force_omega = [];
        hams_data.omega = [];
        hams_data.z_cg  = CG_global(3);
        hams_data.submerged_volume = sub_k.V_sub;
        hams_data.displaced_mass = m_k;
        hams_data.status = 'failed';
        return;
    end

    % 6. Parse output
    one_files = dir(fullfile(hams_dir, 'Output', 'Wamit_format', '*.1'));
    if isempty(one_files)
        warning('mwecmass:hams_mrel:NoOutput', ...
                'No .1 file at vs=%+.4f', vertical_shift);
        hams_data.added_mass_inf = zeros(6);
        hams_data.radiation_damping_band_avg = zeros(6);
        hams_data.added_mass_omega = [];
        hams_data.radiation_damping_omega = [];
        hams_data.exciting_force_omega = [];
        hams_data.omega = [];
        hams_data.z_cg  = CG_global(3);
        hams_data.submerged_volume = sub_k.V_sub;
        hams_data.displaced_mass = m_k;
        hams_data.status = 'no_output';
        return;
    end

    hd = mwecmass.bem.hams_mrel.parse_wamit_1_file( ...
        fullfile(hams_dir, 'Output', 'Wamit_format', one_files(1).name), ...
        1.0, 3);

    % 6b. Parse excitation force (.3 file)
    three_files = dir(fullfile(hams_dir, 'Output', 'Wamit_format', '*.3'));
    if ~isempty(three_files)
        fe_data = mwecmass.bem.hams_mrel.parse_wamit_3_file( ...
            fullfile(hams_dir, 'Output', 'Wamit_format', three_files(1).name), ...
            1.0, 3);
        Fe_complex = fe_data.Fe;   % [6×M] complex, dimensional [N or N·m]
    else
        Fe_complex = [];
        warning('mwecmass:hams_mrel:No3File', ...
                'No .3 file at vs=%+.4f — Fe unavailable', vertical_shift);
    end

    % 7. Assemble output
    hams_data.added_mass_inf = hd.A_inf;
    hams_data.added_mass_omega = hd.A;
    hams_data.radiation_damping_omega = hd.B;
    hams_data.omega  = hd.omega;
    hams_data.radiation_damping_band_avg = mwecmass.bem.hams_mrel.band_averaged_damping(hd, options.T_min, options.T_max);
    hams_data.exciting_force_omega = Fe_complex;
    hams_data.z_cg   = CG_global(3);
    hams_data.submerged_volume = sub_k.V_sub;
    hams_data.displaced_mass = m_k;
    hams_data.status = 'ok';

    % Guard against NaN values in A(∞).
    %  HAMS can write NaN in the A(∞) rows of the .1 file when the
    %  BEM solver diverges at infinite frequency.  Known triggers:
    %    • waterplane very close to the hull apex (converging-cone tip)
    %    • near-degenerate panels at the waterline trim boundary
    %
    %  The frequency-dependent A(ω) at high ω always converges to
    %  A(∞).  When HAMS returns NaN, estimate A(∞) from the mean of
    %  the top-5 frequency points (ω_max ≈ 2.9 rad/s, T ≈ 2.2 s).
    %  This is valid because all WEC natural periods of interest
    %  (T_heave, T_pitch) are far above 2.2 s, so A(ω_max) ≈ A(∞).
    %
    %  If A(ω) data is also unavailable (e.g. HAMS ran with zero
    %  frequencies), A_inf stays zero — which stage1_screen_draft treats as
    %  a degenerate (infeasible) draft, as intended.
    if any(isnan(hams_data.added_mass_inf(:)))
        n_hw = min(5, size(hams_data.added_mass_omega, 3));
        if n_hw >= 1
            warning('mwecmass:hams_mrel:AinfNaN', ...
                ['A(inf) = NaN at vs=%+.4f m — HAMS high-freq BEM failed.\n' ...
                 '  Estimating A(inf) from top-%d frequency points ' ...
                 '(omega=[%.2f..%.2f] rad/s).'], ...
                vertical_shift, n_hw, ...
                hams_data.omega(max(1, end-n_hw+1)), hams_data.omega(end));
            for ii = 1:6
                for jj = 1:6
                    hw = squeeze(hams_data.added_mass_omega(ii, jj, end-n_hw+1:end));
                    if all(isfinite(hw))
                        hams_data.added_mass_inf(ii, jj) = mean(hw);
                    else
                        hams_data.added_mass_inf(ii, jj) = 0;
                    end
                end
            end
            if options.verbose
                fprintf('    A(inf) fallback applied: A33=%.1f kg, A55=%.1f kg*m2\n', ...
                        hams_data.added_mass_inf(3,3), hams_data.added_mass_inf(5,5));
            end
        else
            warning('mwecmass:hams_mrel:AinfNaN', ...
                'A(inf) = NaN at vs=%+.4f m and no A(omega) data — setting A_inf=0.', ...
                vertical_shift);
            hams_data.added_mass_inf = zeros(6);
            hams_data.status = 'ainf_fallback_failed';
        end
    end

    if options.verbose
        fprintf('    A33(inf)=%.1f kg, A55(inf)=%.1f kg*m2\n', ...
                hams_data.added_mass_inf(3,3), hams_data.added_mass_inf(5,5));
    end
end

function setup_hams_directory(run_dir)
%SETUP_HAMS_DIRECTORY  Create the HAMS-compatible folder tree under run_dir.
    %
    %   Creates:  run_dir/Input/
    %             run_dir/Output/Hams_format/
    %             run_dir/Output/Hydrostar_format/
    %             run_dir/Output/Wamit_format/

    dirs = {fullfile(run_dir, 'Input'), ...
            fullfile(run_dir, 'Output', 'Hams_format'), ...
            fullfile(run_dir, 'Output', 'Hydrostar_format'), ...
            fullfile(run_dir, 'Output', 'Wamit_format')};

    for d = 1:length(dirs)
        if ~exist(dirs{d}, 'dir')
            mkdir(dirs{d});
        end
    end

    % Create empty ErrorCheck.txt (required by HAMS)
    err_file = fullfile(run_dir, 'Output', 'ErrorCheck.txt');
    if ~exist(err_file, 'file')
        fid = fopen(err_file, 'w');
        fclose(fid);
    end

    fprintf('  HAMS directory ready: %s\n', run_dir);
end

function write_panelizer_hull_pnl(mesh, filepath)
%WRITE_PANELIZER_HULL_PNL  Write a panel-mesh struct to a HAMS .pnl file.
% HAMS-only wrapper around HamsWriter.pnl_file. mesh comes from
% mwecmass.mesh.generate; triangles use v4=v3 and are emitted with nverts=3.

    nodes  = mesh.vertices;
    panels = mesh.panels;

    % Detect triangles (panel-mesher convention: v4 == v3)
    n_p = size(panels, 1);
    panel_nverts = 4 * ones(n_p, 1);
    for p = 1:n_p
        if panels(p, 3) == panels(p, 4)
            panel_nverts(p) = 3;
        end
    end

    % Symmetry from mesh (quarter_body → [1,1], full body → [0,0])
    if isfield(mesh, 'x_sym'), x_sym = mesh.x_sym; else, x_sym = 0; end
    if isfield(mesh, 'y_sym'), y_sym = mesh.y_sym; else, y_sym = 0; end

    mwecmass.bem.hams_mrel.HamsWriter.pnl_file(filepath, nodes, panels, ...
        panel_nverts, x_sym, y_sym);
end
