function tables = build_hydrostatic_tables(parser, hull_z_min, hull_z_max, aw_table_dz, volume_dt)
%BUILD_HYDROSTATIC_TABLES Build shared hydrostatic tables on an adaptive z grid.
% Syntax: tables = mwecmass.driver.build_hydrostatic_tables(parser,hull_z_min,hull_z_max,aw_table_dz,volume_dt).
% Inputs: parser; z limits [m]; target z spacing [m]; divergence-theorem volume [m^3].
% Output tables share Aw_table_z [N x 1] m and include area [m^2], waterplane moments [m^4],
% perimeter [m], submerged volume [m^3], CB_z [m], wetted area [m^2], and a boundary cache.
% See docs/METHODS_ENGINE.md#hydrostatic-tables

    tables.total_wec_volume = volume_dt;   % [m^3]

    %% Boundary cache and waterplane area table
    %
    %  Boundary cache:
    %    Evaluates all source surface boundary curves ONCE at n_u
    %    u-samples.  Stores the results in a cache struct that
    %    extract_isocurve_at_z reuses for every z-level query.
    %    Cost: ~0.15 s.  Eliminates all redundant eval_curve/eval_snake
    %    calls from the Aw table construction, r_min computation,
    %    and any standalone compute_submerged/compute_strip calls.
    %
    %  Waterplane-area table:
    %    Same adaptive sampling strategy (coarse + refinement at high
    %    |dAw/dz|), but contour extraction now uses type-dispatched
    %    algebraic solvers instead of brute-force bisection sweeps.
    %    Validated: 0.00% Aw error vs reference, 1249× speedup on C0.
    %
    %  USAGE: properties_3d interpolates Aw from this table
    %    and passes it as Aw_override to compute_submerged.  No contour
    %    tracing in the optimizer loop.

    fprintf('  Precomputing boundary cache...\n');
    n_u_cache = 100;
    tables.boundary_cache = mwecmass.geometry.precompute_boundary_cache( ...
                                parser, n_u_cache);
    fprintf('    Sources: %d, mirrors: %d, u-samples: %d\n', ...
            length(tables.boundary_cache.sources), ...
            length(tables.boundary_cache.mirrors), n_u_cache);
    % The u-sample count travels with the cache: the strip-edge augmentation evaluates further
    % isocurves from the same cache and must use the same sampling.
    tables.n_u_cache = n_u_cache;

    fprintf('  Building waterplane area table (adaptive, fast iso-z)...\n');

    % Parametric resolution follows the target z spacing.
    if ~isempty(aw_table_dz) && aw_table_dz > 0
        n_coarse = max(20, ceil((hull_z_max - hull_z_min) / aw_table_dz));
    else
        n_coarse = 40;  % backward compatibility
    end
    z_coarse = linspace(hull_z_min + 1e-4, hull_z_max - 1e-4, n_coarse)';
    Aw_coarse  = zeros(n_coarse, 1);
    Ixx_coarse = zeros(n_coarse, 1);
    Iyy_coarse = zeros(n_coarse, 1);
    P_coarse   = zeros(n_coarse, 1);

    for k = 1:n_coarse
        wl_pts = mwecmass.geometry.extract_isocurve_at_z( ...
                     parser, z_coarse(k), n_u_cache, ...
                     tables.boundary_cache);
        if ~isempty(wl_pts) && size(wl_pts, 1) >= 3
            [Aw_coarse(k), Ixx_coarse(k), Iyy_coarse(k)] = ...
                mwecmass.hydrostatics.waterplane_properties(wl_pts);
            x_cl = [wl_pts(:,1); wl_pts(1,1)];
            y_cl = [wl_pts(:,2); wl_pts(1,2)];
            P_coarse(k) = sum(sqrt(diff(x_cl).^2 + diff(y_cl).^2));
        end
    end

    % Pass 2: adaptive refinement where Aw changes rapidly
    dAw = abs(diff(Aw_coarse));
    dz_c = diff(z_coarse);
    grad_Aw = dAw ./ dz_c;
    grad_thresh = 0.3 * max(grad_Aw);  % refine top 30% gradient

    z_refine   = [];
    Aw_refine  = [];
    Ixx_refine = [];
    Iyy_refine = [];
    P_refine   = [];

    for k = 1:length(grad_Aw)
        if grad_Aw(k) > grad_thresh
            z_mid = linspace(z_coarse(k), z_coarse(k+1), 5)';
            z_mid = z_mid(2:end-1);
            for j = 1:length(z_mid)
                wl_pts = mwecmass.geometry.extract_isocurve_at_z( ...
                             parser, z_mid(j), n_u_cache, ...
                             tables.boundary_cache);
                Aw_j = 0; Ixx_j = 0; Iyy_j = 0; P_j = 0;
                if ~isempty(wl_pts) && size(wl_pts, 1) >= 3
                    [Aw_j, Ixx_j, Iyy_j] = ...
                        mwecmass.hydrostatics.waterplane_properties(wl_pts);
                    x_cl = [wl_pts(:,1); wl_pts(1,1)];
                    y_cl = [wl_pts(:,2); wl_pts(1,2)];
                    P_j  = sum(sqrt(diff(x_cl).^2 + diff(y_cl).^2));
                end
                z_refine(end+1)   = z_mid(j); %#ok<AGROW>
                Aw_refine(end+1)  = Aw_j;     %#ok<AGROW>
                Ixx_refine(end+1) = Ixx_j;    %#ok<AGROW>
                Iyy_refine(end+1) = Iyy_j;    %#ok<AGROW>
                P_refine(end+1)   = P_j;       %#ok<AGROW>
            end
        end
    end

    % Merge coarse + refined, sort by z
    z_all_aw  = [z_coarse(:); z_refine(:)];
    Aw_all    = [Aw_coarse(:);  Aw_refine(:)];
    Ixx_all   = [Ixx_coarse(:); Ixx_refine(:)];
    Iyy_all   = [Iyy_coarse(:); Iyy_refine(:)];
    P_all     = [P_coarse(:);   P_refine(:)];

    [z_all_aw, sort_idx] = sort(z_all_aw);
    Aw_all  = Aw_all(sort_idx);
    Ixx_all = Ixx_all(sort_idx);
    Iyy_all = Iyy_all(sort_idx);
    P_all   = P_all(sort_idx);

    % Add boundary values: Aw=0 (and P=0) at exact hull limits
    z_all_aw  = [hull_z_min; z_all_aw; hull_z_max];
    Aw_all    = [0; Aw_all;  0];
    Ixx_all   = [0; Ixx_all; 0];
    Iyy_all   = [0; Iyy_all; 0];
    P_all     = [0; P_all;   0];

    tables.Aw_table_z     = z_all_aw;
    tables.Aw_table       = Aw_all;
    tables.I_wp_xx_table  = Ixx_all;
    tables.I_wp_yy_table  = Iyy_all;
    tables.P_table        = P_all;

    fprintf('    Aw table: %d points (%d coarse + %d refined)\n', ...
            length(z_all_aw), n_coarse, length(z_refine));
    fprintf('    Aw range: [%.4f, %.4f] m²\n', min(Aw_all), max(Aw_all));

    %% V_sub and CB_z tables
    %  V_sub and CB_z from cumulative Aw integration.
    %
    %  Cumulative integration avoids partial-surface quadrature:
    %    Called compute_submerged (divergence theorem on partial parametric
    %    surfaces) at each z-level.  This produces 13–33% errors for
    %    offset-axis hull families because GL quadrature cannot resolve
    %    the step discontinuity at the waterline crossing.
    %
    %    V_sub(z_wl) = ∫_{z_min}^{z_wl} A(z) dz         (trapz)
    %    CB_z(z_wl)  = ∫_{z_min}^{z_wl} z·A(z) dz / V_sub
    %
    %    Integrates the smooth, validated Aw(z) table.  No surface
    %    integrals, no orientation checks, no cancellation issues.
    %    Accuracy limited only by Aw table resolution (controlled by
    %    aw_table_dz).

    fprintf('  Precomputing V_sub / CB_z tables (%d z-levels, Aw-trapz)...\n', ...
            length(z_all_aw));

    n_vtab      = length(z_all_aw);
    V_sub_table = zeros(n_vtab, 1);
    CB_z_table  = zeros(n_vtab, 1);

    for k_tab = 1:n_vtab
        z_wl_k = z_all_aw(k_tab);

        % Cumulative integral of A(z) from hull bottom to z_wl_k
        mask = z_all_aw <= z_wl_k;
        z_int = z_all_aw(mask);
        A_int = Aw_all(mask);

        if length(z_int) < 2
            V_sub_table(k_tab) = 0;
            CB_z_table(k_tab)  = hull_z_min;
            continue;
        end

        V_sub_table(k_tab) = trapz(z_int, A_int);

        if V_sub_table(k_tab) > 1e-10
            CB_z_table(k_tab) = trapz(z_int, z_int .* A_int) / V_sub_table(k_tab);
        else
            V_sub_table(k_tab) = 0;
            CB_z_table(k_tab)  = hull_z_min;
        end
    end

    tables.V_sub_table = V_sub_table;
    tables.CB_z_table  = CB_z_table;

    fprintf('    V_sub range: [%.4f, %.4f] m^3\n', min(V_sub_table), max(V_sub_table));
    fprintf('    CB_z  range: [%.4f, %.4f] m\n',   min(CB_z_table),  max(CB_z_table));

    % Override: replace the divergence-theorem hull volume with the
    % Aw-trapz maximum.  The DT value (from hp_full.volume) overcounts
    % for hulls whose visible surfaces do NOT form a fully closed
    % boundary.  E1.ms2 Ellipsoid is an open bowl — the DT includes
    % a virtual cap volume (~4 m³) that does not physically exist.
    %
    % V_sub_max = integral Aw(z) dz is always correct: it counts only
    % the cross-sectional area present at each z-level.  For a fully
    % closed hull both values agree; for open hulls the Aw-trapz value
    % is the physically meaningful displaced volume.
    %
    % The hull centroid and the inertia integrals from the divergence
    % theorem are retained — they remain accurate for open hulls.
    V_sub_max_Aw = max(V_sub_table);
    if abs(V_sub_max_Aw - tables.total_wec_volume) / ...
            max(tables.total_wec_volume, 1e-6) > 0.05
        fprintf('    NOTE: DT volume (%.4f m3) differs from Aw-trapz (%.4f m3) by %.1f%%.\n', ...
                tables.total_wec_volume, V_sub_max_Aw, ...
                abs(V_sub_max_Aw - tables.total_wec_volume) / ...
                tables.total_wec_volume * 100);
        fprintf('    Using Aw-trapz value (open-hull / non-closed surface).\n');
    end
    tables.total_wec_volume = V_sub_max_Aw;

    % Wetted surface area table.
    %  S_wet(z_wl) = ∬_{z ≤ z_wl} ||Su×Sv|| du dv  (hull sides only)
    %
    %  Built on the same z-grid as Aw / V_sub so that
    %  properties_3d can interpolate with interp1 on Aw_table_z.
    %
    %  IMPORTANT: mwecmass.driver.build_strip_geometry_tables also
    %  extends S_wet_table, keeping it the same length as Aw_table_z.
    %
    %  Definition: wetted hull area EXCLUDING the waterplane.
    %  Standard use: frictional drag, viscous BEM corrections, Re-scaling.
    %
    %  Guard: only if parser is present (WAMIT-only PATH B/C has
    %  no geometry). Absence sets tables.S_wet_table = [].

    if ~isempty(parser)
        fprintf('  Precomputing S_wet table (%d z-levels, GL n=20)...\n', n_vtab);
        S_wet_table = zeros(n_vtab, 1);
        for k_tab = 1:n_vtab
            z_wl_k = z_all_aw(k_tab);
            if z_wl_k <= hull_z_min
                S_wet_table(k_tab) = 0;
            else
                S_wet_table(k_tab) = mwecmass.hydrostatics.compute_wetted_surface_area( ...
                    parser, z_wl_k);
            end
        end
        tables.S_wet_table = S_wet_table;
        fprintf('    S_wet range: [%.4f, %.4f] m^2\n', ...
                min(S_wet_table), max(S_wet_table));
    else
        tables.S_wet_table = [];   % WAMIT-only or geometry-only: no MS2 geometry
        fprintf('  S_wet table: skipped (no ms2_model)\n');
    end

    % ── Validation: monotonicity and bounds check ──
    dV = diff(V_sub_table);
    n_violations = sum(dV < -1e-6);
    if n_violations > 0
        warning('WEC:VsubTableNonMono', ...
                'V_sub table has %d non-monotonic intervals (max dV = %.4e m^3).', ...
                n_violations, min(dV));
    end
    if min(CB_z_table) < hull_z_min - 0.1 || max(CB_z_table) > hull_z_max + 0.1
        warning('WEC:CBzOutOfBounds', ...
                'CB_z table [%.3f, %.3f] exceeds hull bounds [%.3f, %.3f].', ...
                min(CB_z_table), max(CB_z_table), hull_z_min, hull_z_max);
    end
end
