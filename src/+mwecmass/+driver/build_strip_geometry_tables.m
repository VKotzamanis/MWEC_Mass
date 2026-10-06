function [tables, strips] = build_strip_geometry_tables(parser, tables, strip_layout, hull_z_min, hull_z_max)
%BUILD_STRIP_GEOMETRY_TABLES Insert strip edges into hydrostatic tables and integrate strip geometry.
% Syntax: [tables,strips] = build_strip_geometry_tables(parser,tables,strip_layout,hull_z_min,hull_z_max).
% Inputs: parser; hydrostatic tables; strip layout with edges/nodes; hull z limits [m].
% Outputs: augmented tables and strips containing edges [m], volume [m^3], centroid z [m], and
% second volume moments Ixx/Iyy/Izz [m^5]. Geometry is independent of ballast density and draft.

    n_u_cache   = tables.n_u_cache;
    strip_edges = strip_layout.strip_edges;

    N = strip_layout.num_ballast_sections;

    % Determine strip edges (for the default realisation type, constructability already set)
    if ~strip_layout.enable_constructability || isempty(strip_edges)
        if N > 1
            dz_se = strip_layout.density_nodes_z(2) - strip_layout.density_nodes_z(1);
        else
            dz_se = hull_z_max - hull_z_min;
        end
        se = zeros(N + 1, 1);
        se(1) = strip_layout.density_nodes_z(1) - dz_se/2;
        for i = 2:N
            se(i) = 0.5 * (strip_layout.density_nodes_z(i-1) + strip_layout.density_nodes_z(i));
        end
        se(N+1) = strip_layout.density_nodes_z(N) + dz_se/2;
        se(1)   = max(se(1), hull_z_min);
        se(end) = min(se(end), hull_z_max);
        strip_edges = se;
    end

    % Augment Aw with exact strip-boundary evaluations.
    %
    %  The Aw_table was built on an adaptive grid of 45 z-levels.  Strip
    %  boundary z-values (e.g. z = −0.9 m, the platform-to-column shoulder)
    %  may fall BETWEEN table grid points.  The compute_strip fast path uses
    %  interp1 to get A(z_lo) and A(z_hi), but interpolation is inaccurate
    %  near rapid A(z) transitions (Aw jumps from 0.6 to 10.4 m² over 0.34 m
    %  at the shoulder).
    %
    %  Evaluate A(z) exactly at every strip edge using the parametric
    %  boundary definition (extract_isocurve_at_z), then insert those points
    %  into the Aw_table grid so interp1 finds an exact match.
    %
    %  V_sub_table is on the same z-grid (Aw_table_z); augmenting both
    %  together keeps them consistent for properties_3d, which uses
    %    interp1(config.Aw_table_z, config.V_sub_table, z_wl)
    %  V_sub at the new strip-edge z-values is computed with compute_submerged.

    n_edge_added  = 0;

    for ie = 1:length(strip_edges)
        z_edge = strip_edges(ie);

        % Skip if a table point already exists within numerical tolerance
        if any(abs(tables.Aw_table_z - z_edge) < 1e-8)
            continue;
        end

        % Exact isocurve at the strip edge
        wl_edge = mwecmass.geometry.extract_isocurve_at_z( ...
                      parser, z_edge, n_u_cache, tables.boundary_cache);
        if ~isempty(wl_edge) && size(wl_edge, 1) >= 3
            [Aw_edge, Ixx_edge, Iyy_edge] = ...
                mwecmass.hydrostatics.waterplane_properties(wl_edge);
        else
            Aw_edge  = 0;
            Ixx_edge = 0;
            Iyy_edge = 0;
        end

        % V_sub and CB_z at the edge from cumulative Aw integration.
        %   Same trapz method as the V_sub table builder — no compute_submerged.
        mask_edge = tables.Aw_table_z <= z_edge;
        z_int_e   = tables.Aw_table_z(mask_edge);
        A_int_e   = tables.Aw_table(mask_edge);

        % Sort the integration grid before trapz so cumulative volume remains
        % monotone after edge insertion.
        [z_int_e, sort_e] = sort(z_int_e);
        A_int_e           = A_int_e(sort_e);

        if length(z_int_e) >= 2
            Vsub_edge = trapz(z_int_e, A_int_e);
            if Vsub_edge > 1e-10
                CBz_edge = trapz(z_int_e, z_int_e .* A_int_e) / Vsub_edge;
            else
                Vsub_edge = 0;
                CBz_edge  = hull_z_min;
            end
        else
            Vsub_edge = 0;
            CBz_edge  = hull_z_min;
        end

        % Insert into tables
        tables.Aw_table_z    = [tables.Aw_table_z;    z_edge  ];
        tables.Aw_table      = [tables.Aw_table;      Aw_edge ];
        tables.I_wp_xx_table = [tables.I_wp_xx_table; Ixx_edge];
        tables.I_wp_yy_table = [tables.I_wp_yy_table; Iyy_edge];
        tables.V_sub_table   = [tables.V_sub_table;   Vsub_edge];
        tables.CB_z_table    = [tables.CB_z_table;    CBz_edge ];

        % Augment S_wet_table at this strip edge to keep it the same
        % length as Aw_table_z.  Without this, interp1 in
        % properties_3d crashes with "X and V must be of the
        % same length" because Aw_table_z grows but S_wet_table does not.
        if ~isempty(tables.S_wet_table)
            if z_edge <= hull_z_min
                Swet_edge = 0;
            else
                Swet_edge = mwecmass.hydrostatics.compute_wetted_surface_area( ...
                                parser, z_edge);
            end
            tables.S_wet_table = [tables.S_wet_table; Swet_edge];
        end

        n_edge_added = n_edge_added + 1;
    end

    % Re-sort all tables by z (insertion above was unsorted)
    [tables.Aw_table_z, sort_ie] = sort(tables.Aw_table_z);
    tables.Aw_table      = tables.Aw_table(sort_ie);
    tables.I_wp_xx_table = tables.I_wp_xx_table(sort_ie);
    tables.I_wp_yy_table = tables.I_wp_yy_table(sort_ie);
    tables.V_sub_table   = tables.V_sub_table(sort_ie);
    tables.CB_z_table    = tables.CB_z_table(sort_ie);
    if ~isempty(tables.S_wet_table)
        tables.S_wet_table = tables.S_wet_table(sort_ie);
    end

    % RE-COMPUTE V_sub_table and CB_z_table from scratch on the FINAL sorted
    % grid.  Required because the per-edge trapz above used
    % a DIFFERENT (possibly partial / earlier-state) grid for each entry,
    % producing values that are no longer mutually consistent after sorting.
    % Using the same cumulative-trapz semantics as the table builder guarantees a strictly
    % non-decreasing V_sub_table after augmentation.
    n_vtab2 = length(tables.Aw_table_z);
    for k_tab = 1:n_vtab2
        z_wl_k = tables.Aw_table_z(k_tab);
        mask   = tables.Aw_table_z <= z_wl_k;
        z_int  = tables.Aw_table_z(mask);
        A_int  = tables.Aw_table(mask);

        if length(z_int) < 2
            tables.V_sub_table(k_tab) = 0;
            tables.CB_z_table(k_tab)  = hull_z_min;
            continue;
        end

        Vk = trapz(z_int, A_int);
        if Vk > 1e-10
            tables.V_sub_table(k_tab) = Vk;
            tables.CB_z_table(k_tab)  = trapz(z_int, z_int .* A_int) / Vk;
        else
            tables.V_sub_table(k_tab) = 0;
            tables.CB_z_table(k_tab)  = hull_z_min;
        end
    end

    % Hard guard: assert monotonicity after the rebuild.  If this trips,
    % the Aw_table itself has issues (negative entries or duplicate z's).
    dVsub = diff(tables.V_sub_table);
    if any(dVsub < -1e-9)
        bad = find(dVsub < -1e-9, 1);
        error('mwecmass:driver:VsubNonMonotonic', ...
            'V_sub_table decreased by %.3em³ between z=%.4f and z=%.4f after rebuild — Aw_table likely has negative or duplicate-z entries.', ...
            -dVsub(bad), tables.Aw_table_z(bad), tables.Aw_table_z(bad+1));
    end

    fprintf('  Aw/V_sub tables augmented: %d strip edge(s) added → %d total points\n', ...
            n_edge_added, length(tables.Aw_table_z));
    fprintf('  V_sub_table rebuilt on final grid; monotonic check passed.\n');

    fprintf('  Precomputing strip geometry (%d strips, table-based trapz)...\n', N);
    strips.strip_V       = zeros(N, 1);
    strips.strip_CB_z    = zeros(N, 1);
    strips.strip_Iyy     = zeros(N, 1);
    strips.strip_Ixx     = zeros(N, 1);
    strips.strip_Izz     = zeros(N, 1);

    strip_quad_opts = struct('n_quad', 16, ...
                             'Aw_table_z', tables.Aw_table_z, ...
                             'Aw_table', tables.Aw_table, ...
                             'I_wp_xx_table', tables.I_wp_xx_table, ...
                             'I_wp_yy_table', tables.I_wp_yy_table);
    for i = 1:N
        z_lo_i = strip_edges(i);
        z_hi_i = strip_edges(i + 1);

        if z_hi_i <= z_lo_i + 1e-10
            continue;
        end

        strip_i = mwecmass.hydrostatics.compute_strip( ...
            parser, z_lo_i, z_hi_i, strip_quad_opts);

        strips.strip_V(i)    = strip_i.V;
        strips.strip_CB_z(i) = strip_i.CB_z;
        strips.strip_Iyy(i)  = strip_i.Iyy;
        strips.strip_Ixx(i)  = strip_i.Ixx;
        strips.strip_Izz(i)  = strip_i.Izz;

        fprintf('    Strip %2d [%+6.3f, %+6.3f]: V=%.4f m³, z̄=%.3f m\n', ...
                i, z_lo_i, z_hi_i, strip_i.V, strip_i.CB_z);
    end

    V_sum = sum(strips.strip_V);
    fprintf('  Strip volume sum: %.4f m³ (hull total: %.4f m³, diff: %.2f%%)\n', ...
            V_sum, tables.total_wec_volume, ...
            abs(V_sum - tables.total_wec_volume) / tables.total_wec_volume * 100);

    strips.strip_edges = strip_edges;
end
