function cstr = solve_and_extract(config, x_opt, final_props)
%SOLVE_AND_EXTRACT Solve and extract a constructible modular-precast UHPC/void geometry.
% Uses a wall-aware constrained solve, then re-extracts per-strip contours and verifies mass/feasibility.
% x_opt is [vertical_shift;rho_1..rho_N]; final_props supplies target mass and hydrostatics.
% Returns cstr without modifying optimiser outputs.
% See docs/METHODS_ENGINE.md#realise-modular-precast
    fprintf('\n    mwecmass.realise.modular_precast.solve_and_extract (UHPC two-phase global solve):\n');

    %% RESOLVE UHPC MATERIAL PARAMETERS
    rho_UHPC     = config.constructability_rho_hull;
    rho_void_mat = config.constructability_rho_fill;
    t_min_uhpc   = config.constructability_t_min;

    % Build opts struct overriding Shell_Offset defaults with UHPC params.
    % opt_or_cfg in Shell_Offset reads opts field first, then config field,
    % so overriding here does NOT require changing config.
    uhpc_opts = struct();
    uhpc_opts.rho_steel = rho_UHPC;
    uhpc_opts.rho_air   = rho_void_mat;
    uhpc_opts.t_min     = t_min_uhpc;

    % Warm-start: UHPC needs ~rho_steel/rho_UHPC times more wall volume
    % than steel for the same mass, so scale t_init proportionally.
    if isfield(config, 'uhpc_t_init') && ~isempty(config.uhpc_t_init)
        uhpc_opts.t_init = config.uhpc_t_init;
    else
        uhpc_opts.t_init = config.steel_t_init * ...
                           (config.rho_steel / rho_UHPC);
    end
    if isfield(config, 'uhpc_max_slope_factor')
        uhpc_opts.max_slope_factor = config.uhpc_max_slope_factor;
    end
    if isfield(config, 'uhpc_n_z_grid')
        uhpc_opts.n_z_grid = config.uhpc_n_z_grid;
    end
    % The solve uses range-normalised heave/pitch penalties; mass and GM are constraints.

    fprintf('      rho_UHPC = %.0f kg/m³,  rho_void = %.2f kg/m³\n', ...
            rho_UHPC, rho_void_mat);
    fprintf('      t_min = %.4f m (%.1f mm),  t_init = %.4f m\n', ...
            t_min_uhpc, t_min_uhpc*1000, uhpc_opts.t_init);

    % Wall strip must be pinned solid (matching the solver's constraint) so pre-checks that
    % omit this constraint will understate M_min.
    N_strips = length(config.density_nodes_z);
    density_nodes_z = config.density_nodes_z(:);
    if isfield(config, 'strip_edges') && ~isempty(config.strip_edges)
        strip_edges_realize = config.strip_edges(:);
    else
        strip_edges_realize = zeros(N_strips + 1, 1);
        strip_edges_realize(1) = config.hull_z_min;
        for ii = 2:N_strips
            strip_edges_realize(ii) = 0.5 * ...
                (density_nodes_z(ii-1) + density_nodes_z(ii));
        end
        strip_edges_realize(N_strips + 1) = config.hull_z_max;
    end
    if isfield(config, 'wall_position') && strcmp(config.wall_position, 'bottom')
        wall_strip_idx = 1;
    else
        wall_strip_idx = N_strips;
    end

    %% FEASIBILITY PRE-CHECK
    %
    %  At t = t_min (wall strip pinned solid, matching solve.m's own seed state), M_min =
    %  jacket-only mass (void inside the non-wall strips) plus the solid wall strip's mass.
    %  At t = inf (fully solid everywhere), M_max = rho_UHPC * V_hull -- unaffected by the
    %  per-strip mask, since a fully-solid hull has no strip-dependent inner void either way.
    %  Target mass must lie within [M_min, M_max].
    %
    %  Uses mwecmass.realise.modular_precast.build_geometry_grid (wall strip solid, non-wall
    %  strips at t_min) rather than the uniform-thickness mwecmass.realise.thin_shell.build_geometry_grid:
    %  the latter applies t_min to the wall strip too, so it can UNDERSTATE M_min by the wall
    %  strip's own solid-vs-jacket mass difference and accept a target the modular solver's
    %  wall-pinned model cannot actually reach.

    n_z_pre   = 100;
    slope_f   = mwecmass.internal.option_or_config(uhpc_opts, 'max_slope_factor', ...
                    config.steel_max_slope_factor);
    is_solid_pre = false(N_strips, 1);
    is_solid_pre(wall_strip_idx) = true;
    t_offset_pre = t_min_uhpc * ones(N_strips, 1);
    t_offset_pre(wall_strip_idx) = inf;
    [grids_pre, ~] = mwecmass.realise.modular_precast.build_geometry_grid( ...
        config, strip_edges_realize, t_offset_pre, is_solid_pre, n_z_pre, slope_f);
    A_jacket_pre = grids_pre.A_outer - grids_pre.A_inner;
    M_min_uhpc   = rho_UHPC    * trapz(grids_pre.z, A_jacket_pre) + ...
                   rho_void_mat * trapz(grids_pre.z, grids_pre.A_inner);
    M_max_uhpc   = rho_UHPC    * trapz(grids_pre.z, grids_pre.A_outer);

    fprintf('      Achievable mass range: [%.0f, %.0f] kg  |  target: %.0f kg\n', ...
            M_min_uhpc, M_max_uhpc, final_props.mass_total);

    if final_props.mass_total < M_min_uhpc * 0.95
        error('mwecmass:modular_precast:MassTooLight', ...
              ['Target mass %.1f kg < minimum UHPC+void mass %.1f kg ', ...
               '(at t_min=%.4f m). Hull is too light for UHPC construction.'], ...
              final_props.mass_total, M_min_uhpc, t_min_uhpc);
    end
    if final_props.mass_total > M_max_uhpc * 1.05
        error('mwecmass:modular_precast:MassTooHeavy', ...
              ['Target mass %.1f kg > fully-solid UHPC mass %.1f kg. ', ...
               'Hull is too heavy for UHPC; use a denser material.'], ...
              final_props.mass_total, M_max_uhpc);
    end

    %% PHASE 1: WALL-AWARE GLOBAL SOLVE
    %
    % The constrained solve returns steel_data-compatible fields plus per-strip
    % thickness, solid flags, strip edges, and wall index.

    uhpc_opts.wall_strip_idx = wall_strip_idx;
    uhpc_opts.strip_edges    = strip_edges_realize;

    solve_data = mwecmass.realise.modular_precast.solve( ...
        config, x_opt, final_props, uhpc_opts);

    if ~solve_data.feasible
        warning('mwecmass:modular_precast:InfeasibleUHPCSolve', ...
                ['UHPC constructable solve returned infeasible result. ' ...
                 'Targets may not be achievable under UHPC+void model.']);
    end
    solve_data.fill_method = 'uhpc_fill';

    %% PHASE 2: PER-STRIP GEOMETRY EXTRACTION
    %
    %  Re-runs the strip loop with the per-strip t_offset vector
    %  AND is_solid flags from Phase 1b.  Produces contours_outer,
    %  contours_inner, per-strip volumes/masses for visualisation
    %  and the realised_strips substruct in build_realised_properties.

    cstr = mwecmass.realise.modular_precast.extract_strip_geometry( ...
               config, solve_data.t_offset_strip, solve_data.z_fill, ...
               rho_UHPC, rho_void_mat, solve_data);

    fprintf('    mwecmass.realise.modular_precast.solve_and_extract: done.\n');
end
