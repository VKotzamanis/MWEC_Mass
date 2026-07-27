classdef HAMS_Pipeline
    % HAMS_PIPELINE  MATLAB interface for the HAMS-MREL BEM solver.
    %
    %   Static methods for generating input files, running HAMS, and
    %   parsing output for the WEC mass-distribution optimiser.
    %
    %   METHOD INVENTORY
    %   ────────────────────────────────────────────────────────────────
    %   MESH GENERATION
    %     generate_axisymmetric_mesh  — .pnl mesh from r(z) profile
    %     generate_ellipsoid_profile  — r(z) for a prolate ellipsoid
    %     generate_waterplane_mesh    — flat lid mesh at z=0
    %
    %   FILE WRITERS
    %     write_pnl_file       — HullMesh.pnl / WaterPlaneMesh.pnl
    %     write_control_file   — ControlFile.in
    %     write_hydrostatic_file — Hydrostatic.in
    %
    %   HYDROSTATICS (geometric, uniform density)
    %     compute_geometric_hydrostatics — V_sub, Aw, CB, GM, Iyy, radii
    %
    %   HAMS EXECUTION
    %     setup_hams_directory — create Input/, Output/ folder tree
    %     run_hams             — call HAMS executable
    %
    %   OUTPUT PARSING (WAMIT format — primary path)
    %     parse_wamit_1_file   — read .1 file → A(ω), B(ω), A(∞)
    %     parse_wamit_3_file   — read .3 file → Fe(ω)
    %     read_pnl_file        — read .pnl mesh back into MATLAB
    %
    %   CROSS-VALIDATION
    %     cross_validate_output — compare HAMS vs WAMIT format output
    %
    %   PANELIZER ADAPTERS (WEC_Panelizer → HAMS format)
    %     write_panelizer_hull_pnl — convert panelizer mesh to HAMS .pnl
    %     write_panelizer_wp_pnl   — generate + write waterplane from panelizer mesh
    %     generate_waterplane_mesh_structured — structured quad grid for arbitrary contour
    %
    %   HAMS INPUT BUILDER (from WEC geometry)
    %     compute_hams_inputs    — CG, M_6x6, C_6x6 from parser + draft
    %
    %   HYDRO CACHE (cache-aware single-draft HAMS runs)
    %     empty_hydro_cache      — create empty cache struct
    %     run_single_hams        — one HAMS run at one vertical_shift
    %     get_or_run_hams        — cache-aware wrapper (skip if cached)
    %     rebuild_config_hydro   — repopulate config.wamit_* from cache
    %     retransform_at_actual_cg — re-do origin→CG with optimizer's actual CG
    %
    %   PIPELINE DRIVERS
    %     run_single_draft     — one complete HAMS run at given draft (ellipsoid)
    %     run_5draft_sweep     — 5-draft sweep → A(∞) interpolation table
    %     interpolate_Ainf     — linear interpolation from table
    %     run_draft_sweep      — production draft sweep → hydro_table
    %   ────────────────────────────────────────────────────────────────
    %
    %   CONVENTIONS
    %     - SI units throughout (m, kg, s)
    %     - z = 0 is the still water level
    %     - Body shifted DOWN by draft d: body centre at z = -d
    %     - Wetted surface: z ≤ 0
    %     - Panel normals point INTO the fluid (outward from body)
    %     - Vertex ordering: counterclockwise when viewed from outside
    %
    %   See also: WEC_Configuration_Builder, WEC_File_IO
    %
    %   Author:  WEC Optimisation Team
    %   Version: 1.5 — Hydro cache, run_single_hams, get_or_run_hams (2026-03-23)
    %     v1.4: compute_hams_inputs, run_draft_sweep
    %     v1.3: General WP mesher, panelizer adapters
    %     v1.1: Parser fixes, wave_diffrac_soln, output_freq_type

    properties (Constant)
        % BUG FIX (B1): HAMS internal constants (WavDynMods.f90 line 126):
        %   DATA G,PI,RHO / 9.80665D0, pi, 1000.D0 /
        %   RHO=1000 (NOT 1025). HAMS normalises output by 1000.
        %   Using 1025 here inflated every A_ij and B_ij by 2.5%.
        %
        %   RHO_WATER = 1000 → use for denormalising HAMS output ONLY.
        %   Use your physical seawater density (1025) for your own
        %   hydrostatic calculations (mass, C33, etc.) passed to Hydrostatic.in.
        RHO_WATER = 1000.0;   % HAMS internal ρ [kg/m³] — DO NOT change to 1025
        G         = 9.80665;  % HAMS internal g [m/s²]  — matches WavDynMods.f90
    end

    methods (Static)

        %% ═══════════════════════════════════════════════════════════
        %%  MESH GENERATION
        %% ═══════════════════════════════════════════════════════════

        function [nodes, panels, panel_nverts] = generate_axisymmetric_mesh(...
                r_profile, z_profile, N_theta_half, use_x_symmetry)
            % GENERATE_AXISYMMETRIC_MESH  BEM panel mesh for body of revolution.
            %
            %   [nodes, panels, nverts] = generate_axisymmetric_mesh(r, z, Nt, x_sym)
            %
            %   The profile r(z) defines the half-section in the (r,z) plane.
            %   z must be monotonically increasing with z(end) = 0 (waterline).
            %   r(1) may be zero (closed bottom) or nonzero (open bottom — rare).
            %
            %   With Y-symmetry (always on), only the y ≥ 0 half is meshed.
            %   With X-symmetry (optional), only the x ≥ 0, y ≥ 0 quadrant.
            %
            %   INPUTS
            %     r_profile      [K×1]  radial coordinates, r ≥ 0
            %     z_profile      [K×1]  vertical coordinates, z(end) = 0
            %     N_theta_half   int    azimuthal panels in the y ≥ 0 half
            %                           (with x_sym: panels in the first quadrant)
            %     use_x_symmetry bool   true = mesh only x ≥ 0, y ≥ 0
            %
            %   OUTPUTS
            %     nodes          [N×3]  (x, y, z) coordinates
            %     panels         [M×4]  vertex indices (col 4 = 0 for triangles)
            %     panel_nverts   [M×1]  3 or 4
            %
            %   VALIDATION
            %     For a hemisphere of radius R:
            %       - Total surface area (full body) = 2πR²
            %       - Sum of panel areas should match to within mesh resolution

            K = length(r_profile);
            assert(length(z_profile) == K, 'r and z must be same length');
            assert(abs(z_profile(end)) < 1e-10, 'z(end) must be 0 (waterline)');
            assert(all(diff(z_profile) > 0), 'z must be monotonically increasing');

            % Azimuthal angles
            % Y-symmetry: θ ∈ [0, π].  X+Y symmetry: θ ∈ [0, π/2].
            if use_x_symmetry
                theta_max = pi/2;
            else
                theta_max = pi;
            end
            theta = linspace(0, theta_max, N_theta_half + 1);  % N_theta_half panels

            has_pole = (r_profile(1) < 1e-10);  % closed bottom (south pole)

            % --- Build node array ---
            % If pole exists: one pole node, then rings 2..K
            % If no pole: rings 1..K
            if has_pole
                n_rings = K - 1;          % rings 2..K
                ring_start_idx = 2;       % first ring in r_profile
                N_per_ring = N_theta_half + 1;
                n_nodes = 1 + n_rings * N_per_ring;
            else
                n_rings = K;
                ring_start_idx = 1;
                N_per_ring = N_theta_half + 1;
                n_nodes = n_rings * N_per_ring;
            end

            nodes = zeros(n_nodes, 3);
            node_idx = 0;

            % Pole node
            if has_pole
                node_idx = node_idx + 1;
                nodes(node_idx, :) = [0, 0, z_profile(1)];
            end

            % Ring nodes
            for k = ring_start_idx:K
                r_k = r_profile(k);
                z_k = z_profile(k);
                for j = 1:(N_theta_half + 1)
                    node_idx = node_idx + 1;
                    nodes(node_idx, :) = [r_k * cos(theta(j)), ...
                                          r_k * sin(theta(j)), ...
                                          z_k];
                end
            end
            assert(node_idx == n_nodes, 'Node count mismatch');

            % --- Build panel connectivity ---
            % Helper: node index for ring k (1-based in rings), azimuth j
            if has_pole
                ring_node = @(ring, j) 1 + (ring - 1) * N_per_ring + j;
                pole_id = 1;
            else
                ring_node = @(ring, j) (ring - 1) * N_per_ring + j;
            end

            % Count panels
            if has_pole
                n_panels = N_theta_half + (n_rings - 1) * N_theta_half;
            else
                n_panels = (n_rings - 1) * N_theta_half;
            end
            panels = zeros(n_panels, 4);
            panel_nverts = zeros(n_panels, 1);
            p_idx = 0;

            % Pole cap: triangles from pole to first ring
            if has_pole
                for j = 1:N_theta_half
                    p_idx = p_idx + 1;
                    % Vertices: pole, ring1(j), ring1(j+1)
                    % Order: CCW when viewed from outside (outward normal)
                    v1 = pole_id;
                    v2 = ring_node(1, j+1);   % higher theta first
                    v3 = ring_node(1, j);     % lower theta second
                    panels(p_idx, :) = [v1, v2, v3, 0];
                    panel_nverts(p_idx) = 3;
                end
            end

            % Body panels: quads between adjacent rings
            for ring = 1:(n_rings - 1)
                for j = 1:N_theta_half
                    p_idx = p_idx + 1;
                    % Bottom-left, bottom-right, top-right, top-left
                    % "Bottom" = lower ring (closer to pole), "top" = upper ring
                    v1 = ring_node(ring, j);
                    v2 = ring_node(ring, j+1);
                    v3 = ring_node(ring+1, j+1);
                    v4 = ring_node(ring+1, j);
                    panels(p_idx, :) = [v1, v2, v3, v4];
                    panel_nverts(p_idx) = 4;
                end
            end
            assert(p_idx == n_panels, 'Panel count mismatch');

            % --- Verify panel normals point outward ---
            % Spot-check a body panel: cross product should point away from z-axis
            if n_panels > N_theta_half + 1 && has_pole
                test_p = N_theta_half + 1;  % first quad panel
            elseif n_panels > 0
                test_p = 1;
            end
            v = panels(test_p, :);
            p1 = nodes(v(1), :);
            p2 = nodes(v(2), :);
            if v(4) > 0
                p4 = nodes(v(4), :);
            else
                p4 = nodes(v(3), :);  % triangle: use v3 as substitute
            end
            edge1 = p2 - p1;
            edge2 = p4 - p1;
            normal = cross(edge1, edge2);
            centroid = mean(nodes(v(v>0), :), 1);
            % Outward normal should have positive dot product with radial direction
            radial = [centroid(1), centroid(2), 0];
            if dot(normal, radial) < 0
                warning('HAMS_Pipeline:NormalCheck', ...
                    'Panel normals may point inward. Reversing vertex order.');
                panels(:, [1 2 3 4]) = panels(:, [1 4 3 2]);
                % For triangles with v4=0, swap v2 and v3
                tri_mask = panel_nverts == 3;
                panels(tri_mask, [2 3]) = panels(tri_mask, [3 2]);
            end

            fprintf('  Mesh: %d nodes, %d panels (%d tri + %d quad)\n', ...
                n_nodes, n_panels, sum(panel_nverts==3), sum(panel_nverts==4));
        end


        function [r_profile, z_profile] = generate_ellipsoid_profile(...
                a, c, draft, N_z)
            % GENERATE_ELLIPSOID_PROFILE  Meridional profile for a prolate ellipsoid.
            %
            %   [r, z] = generate_ellipsoid_profile(a, c, draft, N_z)
            %
            %   Prolate ellipsoid: x²/a² + y²/a² + Z²/c² = 1
            %   where Z is body-frame coordinate, c > a (elongated along z).
            %
            %   In world frame (waterline at z=0), body shifted down by draft:
            %     z_world = Z_body - draft
            %
            %   The profile is the wetted meridional curve (z_world ≤ 0).
            %
            %   INPUTS
            %     a      [m]  semi-minor axis (equatorial radius)
            %     c      [m]  semi-major axis (polar half-length)
            %     draft  [m]  downward shift of body centre (0 = symmetric)
            %     N_z    int  number of z-stations (including pole and waterline)
            %
            %   OUTPUTS
            %     r_profile  [N_z×1]  radial coordinates
            %     z_profile  [N_z×1]  vertical coordinates (world frame)
            %
            %   ANALYTICAL VALIDATION (draft = 0):
            %     V_sub = (2/3)π a² c
            %     Aw    = π a²
            %     CB_z  = -3c/8    (centroid of half-ellipsoid)

            assert(c > 0 && a > 0, 'Semi-axes must be positive');
            assert(abs(draft) < c, 'Draft magnitude must be less than semi-major axis');

            % Parametric angle φ: 0 = south pole, φ_max = waterline
            % z_world(φ) = -c cos(φ) - draft
            % At waterline: 0 = -c cos(φ_max) - draft → φ_max = acos(-draft/c)
            phi_max = acos(-draft / c);

            % Uniform spacing in φ gives near-uniform panel aspect ratios
            phi = linspace(0, phi_max, N_z)';

            r_profile = a * sin(phi);
            z_profile = -c * cos(phi) - draft;

            % Force exact values at boundaries
            r_profile(1) = 0;           % pole is exactly zero
            z_profile(end) = 0;         % waterline is exactly zero

            fprintf('  Ellipsoid profile: a=%.3f m, c=%.3f m, draft=%.3f m\n', a, c, draft);
            fprintf('  z range: [%.4f, 0] m, r_waterline = %.4f m\n', ...
                z_profile(1), r_profile(end));
        end


        %% ═══════════════════════════════════════════════════════════
        %%  FILE WRITERS
        %% ═══════════════════════════════════════════════════════════

        function write_pnl_file(filepath, nodes, panels, panel_nverts, ...
                x_sym, y_sym)
            % WRITE_PNL_FILE  Write HAMS-format .pnl mesh file.
            %
            %   write_pnl_file(path, nodes, panels, nverts, x_sym, y_sym)
            %
            %   FORMAT (from HAMS documentation):
            %     Line 1: header
            %     Line 3: # Panels  # Nodes  X-Symmetry  Y-Symmetry
            %     Nodes:  node_id  x  y  z
            %     Panels: panel_id  nverts  v1  v2  v3  [v4]

            n_nodes = size(nodes, 1);
            n_panels = size(panels, 1);

            fid = fopen(filepath, 'w');
            if fid == -1
                error('Cannot open file for writing: %s', filepath);
            end

            fprintf(fid, '    --------------Hull Mesh File---------------\n');
            fprintf(fid, ' \n');
            fprintf(fid, '    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry\n');
            fprintf(fid, '    %8d    %8d       %5d       %5d\n', ...
                n_panels, n_nodes, x_sym, y_sym);
            fprintf(fid, ' \n');

            % Node coordinates
            fprintf(fid, '    # Start Definition of Node Coordinates     ! node_number   x   y   z\n');
            for i = 1:n_nodes
                fprintf(fid, ' %4d    %14.6f    %14.6f    %14.6f\n', ...
                    i, nodes(i,1), nodes(i,2), nodes(i,3));
            end
            fprintf(fid, '    # End Definition of Node Coordinates\n');
            fprintf(fid, ' \n');

            % Panel connectivity
            fprintf(fid, '  # Start Definition of Node Relations   ! panel_number  number_of_vertices   Vertex1_ID   Vertex2_ID   Vertex3_ID   (Vertex4_ID\n');
            for i = 1:n_panels
                nv = panel_nverts(i);
                if nv == 4
                    fprintf(fid, ' %4d    %d    %6d    %6d    %6d    %6d\n', ...
                        i, 4, panels(i,1), panels(i,2), panels(i,3), panels(i,4));
                else
                    fprintf(fid, ' %4d    %d    %6d    %6d    %6d\n', ...
                        i, 3, panels(i,1), panels(i,2), panels(i,3));
                end
            end
            fprintf(fid, '    # End Definition of Node Relations\n');
            fprintf(fid, ' \n');
            fprintf(fid, '    --------------End Hull Mesh File---------------\n');

            fclose(fid);
            fprintf('  Wrote %s (%d nodes, %d panels)\n', filepath, n_nodes, n_panels);
        end


        function write_control_file(filepath, params)
            % WRITE_CONTROL_FILE  Write HAMS ControlFile.in.
            %
            %   write_control_file(path, params)
            %
            %   PARAMS struct fields:
            %     .depth              [m]   water depth (-1 for deep water)
            %     .zero_inf_limits    0|1   compute A(0) and A(∞)
            %     .input_freq_type    int   3=wave frequency (rad/s)
            %     .output_freq_type   int   3=wave frequency
            %     .n_frequencies      int   negative = auto-range
            %     .min_frequency      [depends on type]
            %     .freq_step          [depends on type]
            %     .n_headings         int   (negative = auto-range)
            %     .heading_values     [deg]  (if n_headings > 0)
            %     .min_heading        [deg]  (if n_headings < 0)
            %     .heading_step       [deg]  (if n_headings < 0)
            %     .ref_body_center    [3×1 m]
            %     .ref_body_length    [m]
            %     .wave_diffrac_soln  1|2   (1=diffraction+incident, 2=diffraction only)
            %     .remove_irr_freq    0|1
            %     .n_threads          int
            %     .n_field_points     int
            %     .field_points       [N×3 m]

            % CRITICAL: label widths match Fortran READ field widths in InputFiles.f90.
            % Verified against write_hams_controlfile.m (unit-tested).
            % Wrong widths cause HAMS to silently read 0 for IRSP/NTHREAD/ISOL.
            fid = fopen(filepath, 'w');
            if fid == -1
                error('Cannot open file for writing: %s', filepath);
            end

            % Lines 1-2: consumed by read(*) x2
            fprintf(fid, '   --------------HAMS Control file---------------\n');
            fprintf(fid, '\n');

            % Waterdepth  [14x,f30.15]
            fprintf(fid, '%-14s%30.15f\n', '   Waterdepth ', params.depth);
            fprintf(fid, '\n');
            fprintf(fid, '   #Start Definition of Wave Frequencies\n');

            % SYBO  [27x,i16]
            fprintf(fid, '%-27s%16d\n', '    0_inf_frequency_limits ', params.zero_inf_limits);
            % INFT  [25x,i16]
            fprintf(fid, '%-25s%16d\n', '    Input_frequency_type ', params.input_freq_type);
            % OUFT  [25x,i17]
            fprintf(fid, '%-25s%17d\n', '    Output_frequency_type', params.output_freq_type);
            % NPET  [26x,i16]
            fprintf(fid, '%-26s%16d\n', '    Number_of_frequencies ', params.n_frequencies);
            if params.n_frequencies < 0
                % WK1  [27x,f30.15]
                fprintf(fid, '%-27s%30.15f\n', '    Minimum_frequency_Wmin ', params.min_frequency);
                % DWK  [19x,f30.15]
                fprintf(fid, '%-19s%30.15f\n', '    Frequency_step ', params.freq_step);
            else
                freq_str = sprintf('%.6f ', params.frequency_list);
                fprintf(fid, '    %s\n', strtrim(freq_str));
            end
            fprintf(fid, '   #End Definition of Wave Frequencies\n');
            fprintf(fid, '\n');

            % Wave headings — always use auto-range (n_headings < 0)
            fprintf(fid, '   #Start Definition of Wave Headings\n');
            if params.n_headings < 0
                % NBETA  [23x,i16]
                fprintf(fid, '%-23s%16d\n', '    Number_of_headings ', params.n_headings);
                % BETA1  [20x,f30.15]
                fprintf(fid, '%-20s%30.15f\n', '    Minimum_heading ', params.min_heading);
                % DBETA  [17x,f30.15]
                fprintf(fid, '%-17s%30.15f\n', '    Heading_step ', params.heading_step);
            else
                % Positive n_headings: emit as auto-range with first value
                fprintf(fid, '%-23s%16d\n', '    Number_of_headings ', -1);
                h0 = params.heading_values(1);
                fprintf(fid, '%-20s%30.15f\n', '    Minimum_heading ', h0);
                fprintf(fid, '%-17s%30.15f\n', '    Heading_step ', 90.0);
            end
            fprintf(fid, '   #End Definition of Wave Headings\n');
            fprintf(fid, '\n');

            % XR   [28x,3f12.3]
            fprintf(fid, '%-28s%12.3f%12.3f%12.3f\n', ...
                '    Reference_body_center   ', ...
                params.ref_body_center(1), params.ref_body_center(2), params.ref_body_center(3));
            % REFL  [26x,f30.15]
            fprintf(fid, '%-26s%30.15f\n', '    Reference_body_length ', params.ref_body_length);
            % ISOL  [26x,i16]  — label EXACTLY 26 chars
            fprintf(fid, '%-26s%16d\n', '    Wave_diffrac_solution ', params.wave_diffrac_soln);
            % IRSP  [23x,i16]  — label EXACTLY 23 chars
            fprintf(fid, '%-23s%16d\n', '    If_remove_irr_freq ', params.remove_irr_freq);
            % NTHREAD [23x,i16]  — label EXACTLY 23 chars
            fprintf(fid, '%-23s%16d\n', '    Number of threads  ', params.n_threads);
            fprintf(fid, '\n');
            fprintf(fid, '   #Start Definition of Pressure and/or Elevation (PE)\n');

            % NFP  [27x,i16]
            fprintf(fid, '%-27s%16d\n', '    Number_of_field_points ', params.n_field_points);
            for i = 1:params.n_field_points
                fprintf(fid, '    %f    %f    %f    Global_coords_point_%d\n', ...
                    params.field_points(i,1), params.field_points(i,2), ...
                    params.field_points(i,3), i);
            end
            % Three lines required by InputFiles.f90 after NFP (lines 190,191,193)
            fprintf(fid, '   #End Definition of Pressure and/or Elevation\n');
            fprintf(fid, '   ----------End HAMS Control file----------\n');
            fprintf(fid, '\n');

            fclose(fid);
            fprintf('  Wrote %s\n', filepath);
        end


        function write_hydrostatic_file(filepath, CG, M_6x6, ...
                B_ext_lin, B_ext_quad, C_hydro, K_ext)
            % WRITE_HYDROSTATIC_FILE  Write HAMS Hydrostatic.in.
            %
            %   write_hydrostatic_file(path, CG, M, B_lin, B_quad, C, K)
            %
            %   INPUTS
            %     CG        [3×1 m]    centre of gravity
            %     M_6x6     [6×6]      body mass matrix
            %     B_ext_lin [6×6]      external linear damping
            %     B_ext_quad[6×6]      external quadratic damping
            %     C_hydro   [6×6]      hydrostatic restoring matrix
            %     K_ext     [6×6]      external restoring (PTO/mooring)

            fid = fopen(filepath, 'w');
            if fid == -1
                error('Cannot open file for writing: %s', filepath);
            end

            % CG — Fortran reads with format similar to matrices.
            %   Use E12.5 fields (12 chars, 5 decimal digits) with
            %   2-space prefix per value = 14 chars per field.
            %   %.15E produced 21-char numbers → Fortran parse failure.
            fprintf(fid, ' Center of Gravity:\n');
            fprintf(fid, '  %12.5E  %12.5E  %12.5E\n', CG(1), CG(2), CG(3));

            % Mass matrix
            fprintf(fid, ' Body Mass Matrix:\n');
            HAMS_Pipeline.write_6x6_matrix(fid, M_6x6);

            % External linear damping
            fprintf(fid, ' External Linear Damping Matrix:\n');
            HAMS_Pipeline.write_6x6_matrix(fid, B_ext_lin);

            % External quadratic damping
            fprintf(fid, ' External Quadratic Damping Matrix:\n');
            HAMS_Pipeline.write_6x6_matrix(fid, B_ext_quad);

            % Hydrostatic restoring
            fprintf(fid, ' Hydrostatic Restoring Matrix:\n');
            HAMS_Pipeline.write_6x6_matrix(fid, C_hydro);

            % External restoring
            fprintf(fid, ' External Restoring Matrix:\n');
            HAMS_Pipeline.write_6x6_matrix(fid, K_ext);

            fclose(fid);
            fprintf('  Wrote %s\n', filepath);
        end


        %% ═══════════════════════════════════════════════════════════
        %%  6x6 MATRIX BUILDERS (about origin, general body)
        %% ═══════════════════════════════════════════════════════════

        function M = build_mass_matrix_6x6(m, r_G, J_CG)
            % BUILD_MASS_MATRIX_6X6  General 6x6 mass matrix about ORIGIN.
            %
            %   M = build_mass_matrix_6x6(m, r_G, J_CG)
            %
            %   INPUTS
            %     m      scalar     body mass [kg]
            %     r_G    [3x1]      CG position in global frame [x_G, y_G, z_G] [m]
            %     J_CG   [3x3]      inertia tensor about CG (physics sign convention:
            %                        negative products of inertia on off-diagonals)
            %
            %   OUTPUT
            %     M      [6x6]      symmetric mass matrix about (0,0,0)
            %                        DOF order: [surge sway heave roll pitch yaw]
            %
            %   Uses the full Steiner (parallel axis) transformation:
            %     J^O = J_CG + m * [(r_G.r_G) I_3  -  r_G (x) r_G]
            %
            %   The translation-rotation coupling block is m * S(r_G)
            %   where S is the skew-symmetric matrix.
            %
            %   NO symmetry assumed. Works for arbitrary CG location
            %   and arbitrary inertia tensor.
            %
            %   See: HAMS_MREL_MASTER_REFERENCE_v2.md, Section 4.

            r_G = r_G(:);  % ensure column
            x = r_G(1); y = r_G(2); z = r_G(3);

            % Steiner correction: J^O = J_CG + m*((r.r)*I - r*r')
            d2 = x^2 + y^2 + z^2;
            J_O = J_CG + m * (d2 * eye(3) - r_G * r_G');

            % Skew-symmetric coupling: m * S(r_G)
            %   S(r) = [  0   z  -y ]
            %          [ -z   0   x ]
            %          [  y  -x   0 ]
            C_block = m * [  0,  z, -y;
                          -z,  0,  x;
                           y, -x,  0 ];

            % Assemble 6x6 (symmetric)
            M = [ m*eye(3),   C_block;
                  C_block',   J_O     ];
        end


        %% ═══════════════════════════════════════════════════════════
        %%  GEOMETRIC HYDROSTATICS (uniform density, shape only)
        %% ═══════════════════════════════════════════════════════════

        function hydro = compute_geometric_hydrostatics(...
                r_profile, z_profile, rho_water, g)
            % COMPUTE_GEOMETRIC_HYDROSTATICS  Hydrostatic properties from r(z).
            %
            %   hydro = compute_geometric_hydrostatics(r, z, rho_w, g)
            %
            %   Computes properties assuming UNIFORM density and equilibrium.
            %   Uses Pappus' theorem for bodies of revolution.
            %
            %   All results are geometric — independent of mass distribution.
            %   This is exactly what Hydrostatic.in needs: uniform-density
            %   mass matrix and hydrostatic restoring, which give correct A(∞).
            %
            %   OUTPUTS (struct)
            %     .V_sub       [m³]     submerged volume (z ≤ 0 portion)
            %     .Aw          [m²]     waterplane area = π r(0)²
            %     .CB_z        [m]      centre of buoyancy (z-coordinate)
            %     .mass        [kg]     equilibrium mass = rho_w * V_sub
            %     .CG_z_uniform [m]     CG of uniform-density submerged volume
            %     .Iyy_wp      [m⁴]     second moment of waterplane area about x-axis
            %     .KB          [m]      vertical distance from keel to CB (positive up)
            %     .BM          [m]      Iyy_wp / V_sub
            %     .KG          [m]      vertical distance from keel to CG (uniform density)
            %     .GM          [m]      metacentric height = KB + BM - KG
            %     .r_gyration_xx [m]    radius of gyration about x (geometric)
            %     .r_gyration_yy [m]    radius of gyration about y (geometric)
            %     .r_gyration_zz [m]    radius of gyration about z (geometric)
            %     .Ixx, Iyy, Izz [kg·m²] mass moments of inertia (uniform density)
            %     .C33         [N/m]    hydrostatic heave restoring
            %     .C44         [N·m]    hydrostatic roll restoring
            %     .C55         [N·m]    hydrostatic pitch restoring
            %     .M_6x6       [6×6]    mass matrix for Hydrostatic.in
            %     .C_6x6       [6×6]    hydrostatic restoring for Hydrostatic.in

            K = length(r_profile);

            % --- Submerged volume via Pappus (disk integration) ---
            % V = π ∫ r(z)² dz  (trapezoidal rule)
            r2 = r_profile.^2;
            hydro.V_sub = pi * trapz(z_profile, r2);

            % --- Waterplane area ---
            r_wp = r_profile(end);  % radius at z = 0
            hydro.Aw = pi * r_wp^2;

            % --- Centre of buoyancy ---
            % CB_z = ∫ z · π r² dz / V_sub
            z_r2 = z_profile .* r2;
            hydro.CB_z = pi * trapz(z_profile, z_r2) / hydro.V_sub;

            % --- Equilibrium mass ---
            hydro.mass = rho_water * hydro.V_sub;

            % --- CG for uniform density (same as CB for uniform) ---
            hydro.CG_z_uniform = hydro.CB_z;

            % --- Second moment of waterplane area ---
            % For circular waterplane: Iyy_wp = Ixx_wp = π r_wp⁴ / 4
            hydro.Iyy_wp = pi * r_wp^4 / 4;

            % --- Metacentric height (uniform density) ---
            z_keel = z_profile(1);
            hydro.KB = hydro.CB_z - z_keel;
            hydro.BM = hydro.Iyy_wp / hydro.V_sub;
            hydro.KG = hydro.CG_z_uniform - z_keel;
            hydro.GM = hydro.KB + hydro.BM - hydro.KG;

            % --- Mass moments of inertia (uniform density, body of revolution) ---
            % For a body of revolution about z-axis with density ρ:
            %
            %   Izz = ρ π/2 ∫ r⁴ dz          (polar moment about z)
            %   Ixx = Iyy = ρ π ∫ [r²/4 + z²] r² dz
            %                    = ρ [π/4 ∫ r⁴ dz + π ∫ z² r² dz]
            %
            % But these are about the ORIGIN. We need about CG.
            % Parallel axis: I_cg = I_origin - m * d²

            rho_uniform = hydro.mass / hydro.V_sub;  % = rho_water

            r4 = r_profile.^4;
            z2_r2 = (z_profile.^2) .* r2;

            Izz_origin = rho_uniform * (pi/2) * trapz(z_profile, r4);
            Ixx_origin = rho_uniform * (pi/4 * trapz(z_profile, r4) + ...
                         pi * trapz(z_profile, z2_r2));
            Iyy_origin = Ixx_origin;  % axisymmetric

            % Parallel axis to CG (CG is on z-axis at CG_z)
            cg_z = hydro.CG_z_uniform;
            hydro.Izz = Izz_origin;  % z-axis passes through CG for axisymmetric
            hydro.Ixx = Ixx_origin - hydro.mass * cg_z^2;
            hydro.Iyy = Iyy_origin - hydro.mass * cg_z^2;

            % Radii of gyration
            hydro.r_gyration_xx = sqrt(hydro.Ixx / hydro.mass);
            hydro.r_gyration_yy = sqrt(hydro.Iyy / hydro.mass);
            hydro.r_gyration_zz = sqrt(hydro.Izz / hydro.mass);

            % --- Hydrostatic restoring matrix ---
            % C(3,3) = ρ g Aw
            % C(4,4) = ρ g Iyy_wp + ρ g V_sub * CB_z - m g CG_z
            %        = ρ g [Iyy_wp + V_sub * (CB_z - CG_z)]
            %   (For uniform density, CB_z = CG_z, so C44 = ρ g Iyy_wp)
            % C(5,5) = same as C(4,4) for axisymmetric body
            % C(3,5) = ρ g ∫∫_wp x dA = 0 (symmetric waterplane)

            hydro.C33 = rho_water * g * hydro.Aw;
            hydro.C44 = rho_water * g * hydro.Iyy_wp + ...
                        rho_water * g * hydro.V_sub * hydro.CB_z - ...
                        hydro.mass * g * cg_z;
            hydro.C55 = hydro.C44;  % axisymmetric

            % --- Build 6x6 matrices for Hydrostatic.in (about ORIGIN) ---
            % The mass matrix must be about (0,0,0) per our locked convention.
            % compute_geometric_hydrostatics gives Ixx, Iyy, Izz about the CG.
            % We use build_mass_matrix_6x6 for the full Steiner transformation.
            J_CG = diag([hydro.Ixx, hydro.Iyy, hydro.Izz]);
            r_G = [0, 0, cg_z];
            hydro.M_6x6 = HAMS_Pipeline.build_mass_matrix_6x6(hydro.mass, r_G, J_CG);

            hydro.C_6x6 = zeros(6, 6);
            hydro.C_6x6(3,3) = hydro.C33;
            hydro.C_6x6(4,4) = hydro.C44;
            hydro.C_6x6(5,5) = hydro.C55;

            % CG vector (for Hydrostatic.in)
            hydro.CG = [0, 0, cg_z];
            hydro.R_wp = r_wp;  % waterplane radius

            fprintf('  Hydrostatics: V=%.4f m³, Aw=%.4f m², CB_z=%.4f m\n', ...
                hydro.V_sub, hydro.Aw, hydro.CB_z);
            fprintf('  GM=%.4f m, C33=%.1f N/m, C55=%.1f N·m\n', ...
                hydro.GM, hydro.C33, hydro.C55);
            fprintf('  Mass=%.1f kg, Iyy=%.1f kg·m², r_yy=%.4f m\n', ...
                hydro.mass, hydro.Iyy, hydro.r_gyration_yy);
        end


        %% ═══════════════════════════════════════════════════════════
        %%  HAMS EXECUTION
        %% ═══════════════════════════════════════════════════════════

        function setup_hams_directory(run_dir)
            % SETUP_HAMS_DIRECTORY  Create HAMS-compatible folder structure.
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


        function [status, result] = run_hams(hams_exe, run_dir)
            % RUN_HAMS  Execute HAMS-MREL solver.
            %
            %   [status, result] = run_hams(hams_exe, run_dir)
            %
            %   Uses explicit path mode: exe <InputDir> <OutputDir>
            %   where InputDir = run_dir/Input, OutputDir = run_dir/Output.
            %   Calling convention matches both the Windows and Linux
            %   binaries (see HAMS-MREL_Fedora/readme.md).
            %
            %   LINUX RUNTIME REQUIREMENT
            %   The Fedora binary is dynamically linked against Intel oneAPI
            %   (libmkl_*, libiomp5, libimf).  These are NOT installed in
            %   ldconfig — they're picked up via LD_LIBRARY_PATH set by
            %   /opt/intel/oneapi/setvars.sh.  MATLAB inherits the env from
            %   its launching shell, so MATLAB must be started from a
            %   terminal that has already sourced setvars.sh, e.g.:
            %     $ source /opt/intel/oneapi/setvars.sh && matlab
            %
            %   Returns status (0 = success) and console output.

            assert(exist(hams_exe, 'file') == 2, ...
                'HAMS executable not found: %s', hams_exe);

            % LINUX/MAC PRE-FLIGHT: verify Intel oneAPI is on LD_LIBRARY_PATH.
            % Bail out with an actionable message before invoking system();
            % otherwise the binary aborts with a cryptic "libmkl_intel_lp64.so.3:
            % cannot open shared object file" and the only signal upstream is
            % a non-zero status code.
            if ~ispc
                ld_path = getenv('LD_LIBRARY_PATH');
                if ~contains(ld_path, 'intel/oneapi') && ...
                        ~contains(ld_path, 'intel\oneapi')
                    error('HAMS_Pipeline:OneAPIEnvMissing', ...
                        ['Intel oneAPI runtime not on LD_LIBRARY_PATH.\n' ...
                         'HAMS-MREL is linked against Intel MKL + iomp5 and ' ...
                         'cannot load them without setvars.sh.\n\n' ...
                         'Fix: quit MATLAB, then from a terminal run\n' ...
                         '    source /opt/intel/oneapi/setvars.sh\n' ...
                         '    matlab\n' ...
                         '(adjust the setvars path if your oneAPI install is elsewhere).']);
                end
            end

            input_dir  = fullfile(run_dir, 'Input');
            output_dir = fullfile(run_dir, 'Output');

            cmd = sprintf('"%s" "%s" "%s"', hams_exe, input_dir, output_dir);

            fprintf('  Running HAMS: %s\n', cmd);
            tic;
            [status, result] = system(cmd);
            elapsed = toc;

            % HAMS-MREL returns exit code 0 even on input-file errors (it
            % prints the diagnostic, then `stop` falls through with status
            % 0).  Scan stdout for the known abort markers so callers
            % don't see "success" when no .1 file was produced.
            failure_markers = { ...
                'Terminating application', ...
                'Input files missing', ...
                'input file missing', ...
                'Error opening', ...
                'Error encountered reading'};
            stdout_failed = false;
            for k = 1:numel(failure_markers)
                if contains(result, failure_markers{k})
                    stdout_failed = true;
                    break;
                end
            end

            if status == 0 && ~stdout_failed
                fprintf('  HAMS completed in %.1f seconds\n', elapsed);
            else
                if status == 0 && stdout_failed
                    status = 1;   % surface stdout failure to caller
                end
                warning('HAMS_Pipeline:HAMSFailed', ...
                    'HAMS aborted (status=%d).\nOutput:\n%s', status, result);
            end
        end


        %% ═══════════════════════════════════════════════════════════
        %%  OUTPUT PARSING
        %% ═══════════════════════════════════════════════════════════

        function hams_data = parse_wamit_1_file(filepath, ref_body_length, output_freq_type)
            % PARSE_WAMIT_1_FILE  Read HAMS output in WAMIT .1 format.
            %
            %   hams_data = parse_wamit_1_file(filepath)
            %   hams_data = parse_wamit_1_file(filepath, L)
            %   hams_data = parse_wamit_1_file(filepath, L, freq_type)
            %
            %   CRITICAL: HAMS-MREL uses the output_frequency_type from
            %   ControlFile.in for the first column of the .1 file.
            %     freq_type = 3 -> first column is omega [rad/s]  (DEFAULT)
            %     freq_type = 4 -> first column is period [s]
            %   Standard WAMIT always uses period.  HAMS does NOT.
            %
            %   WAMIT non-dimensionalization:
            %     A_bar(i,j) = A(i,j) / (rho L^k)
            %     B_bar(i,j) = B(i,j) / (rho omega L^k)
            %   where k = 3 + (number of rotational DOFs in {i,j}).
            %
            %   OUTPUT (ALL DIMENSIONAL, SI units):
            %     .A_inf       [6x6]     added mass at omega -> inf
            %     .A_zero      [6x6]     added mass at omega -> 0
            %     .A           [6x6xN]   freq-dependent added mass
            %     .B           [6x6xN]   freq-dependent radiation damping
            %     .omega       [Nx1]     angular freq [rad/s], ASCENDING
            %     .periods     [Nx1]     wave periods [s]
            %     .L           scalar    reference body length used
            %     .output_freq_type  int frequency type used

            if nargin < 2 || isempty(ref_body_length)
                ref_body_length = 1.0;
            end
            if nargin < 3 || isempty(output_freq_type)
                output_freq_type = 3;  % DEFAULT: omega in rad/s
            end

            L = ref_body_length;
            rho = HAMS_Pipeline.RHO_WATER;

            assert(exist(filepath, 'file') == 2, ...
                'WAMIT .1 file not found: %s', filepath);
            assert(ismember(output_freq_type, [3, 4]), ...
                'output_freq_type must be 3 (rad/s) or 4 (period)');

            % -- Read line-by-line (handles mixed 4/5 column format) --
            fid = fopen(filepath, 'r');
            col1_raw = [];
            i_modes = [];
            j_modes = [];
            A_vals = [];
            B_vals = [];

            while ~feof(fid)
                line = fgetl(fid);
                if ~ischar(line) || isempty(strtrim(line))
                    continue;
                end
                vals = sscanf(line, '%f');
                if length(vals) >= 4
                    col1_raw(end+1,1) = vals(1); %#ok<AGROW>
                    i_modes(end+1,1) = round(vals(2)); %#ok<AGROW>
                    j_modes(end+1,1) = round(vals(3)); %#ok<AGROW>
                    A_vals(end+1,1) = vals(4); %#ok<AGROW>
                    if length(vals) >= 5
                        B_vals(end+1,1) = vals(5); %#ok<AGROW>
                    else
                        B_vals(end+1,1) = NaN; %#ok<AGROW>
                    end
                end
            end
            fclose(fid);

            % Store raw values
            hams_data.raw_nondim.col1 = col1_raw;
            hams_data.raw_nondim.i_modes = i_modes;
            hams_data.raw_nondim.j_modes = j_modes;
            hams_data.raw_nondim.A = A_vals;
            hams_data.raw_nondim.B = B_vals;
            hams_data.L = L;
            hams_data.output_freq_type = output_freq_type;

            % -- Denormalization exponent --
            k_exp = zeros(6, 6);
            for ii = 1:6
                for jj = 1:6
                    n_rot = (ii > 3) + (jj > 3);
                    k_exp(ii, jj) = 3 + n_rot;
                end
            end

            % -- Separate special and regular frequencies --
            inf_mask  = col1_raw < 0;    % omega -> inf
            zero_mask = col1_raw == 0;   % omega -> 0
            reg_mask  = col1_raw > 0;

            unique_col1 = unique(col1_raw(reg_mask));

            % -- Convert first column to omega --
            if output_freq_type == 3
                omega_all = unique_col1;          % already omega
            elseif output_freq_type == 4
                omega_all = 2 * pi ./ unique_col1; % period -> omega
            end

            % Sort omega ascending
            [omega_sorted, sort_idx] = sort(omega_all, 'ascend');
            col1_sorted = unique_col1(sort_idx);
            n_freq = length(omega_sorted);

            % Initialize dimensional matrices
            hams_data.A_inf  = zeros(6, 6);
            hams_data.A_zero = zeros(6, 6);
            hams_data.A = zeros(6, 6, n_freq);
            hams_data.B = zeros(6, 6, n_freq);

            % -- Populate A(inf) --
            for idx = find(inf_mask)'
                i = i_modes(idx);  j = j_modes(idx);
                if i >= 1 && i <= 6 && j >= 1 && j <= 6
                    hams_data.A_inf(i, j) = A_vals(idx) * rho * L^k_exp(i,j);
                end
            end

            % -- Populate A(0) --
            for idx = find(zero_mask)'
                i = i_modes(idx);  j = j_modes(idx);
                if i >= 1 && i <= 6 && j >= 1 && j <= 6
                    hams_data.A_zero(i, j) = A_vals(idx) * rho * L^k_exp(i,j);
                end
            end

            % -- Populate frequency-dependent A(omega), B(omega) --
            for idx = find(reg_mask)'
                c1 = col1_raw(idx);
                i = i_modes(idx);  j = j_modes(idx);
                if i < 1 || i > 6 || j < 1 || j > 6; continue; end

                f_idx = find(col1_sorted == c1, 1);
                if isempty(f_idx); continue; end

                omega_k = omega_sorted(f_idx);

                scale_A = rho * L^k_exp(i,j);
                scale_B = rho * omega_k * L^k_exp(i,j);

                hams_data.A(i, j, f_idx) = A_vals(idx) * scale_A;
                if ~isnan(B_vals(idx))
                    hams_data.B(i, j, f_idx) = B_vals(idx) * scale_B;
                end
            end

            hams_data.omega = omega_sorted;
            hams_data.periods = 2 * pi ./ omega_sorted;

            fprintf('  Parsed %s: %d frequencies, omega = [%.3f, %.3f] rad/s\n', ...
                filepath, n_freq, omega_sorted(1), omega_sorted(end));
            fprintf('    A33(inf) = %.1f kg,  A55(inf) = %.1f kg*m^2  [dimensional]\n', ...
                hams_data.A_inf(3,3), hams_data.A_inf(5,5));
            if output_freq_type == 3
                fprintf('    (first column = omega [rad/s])\n');
            else
                fprintf('    (first column = period [s])\n');
            end
        end


        %% ═══════════════════════════════════════════════════════════
        %%  WATERPLANE LID MESH (for irregular frequency removal)
        %% ═══════════════════════════════════════════════════════════

        function [nodes, panels, nverts] = generate_waterplane_mesh(...
                r_waterplane, N_theta_half, N_radial, use_x_symmetry)
            % GENERATE_WATERPLANE_MESH  Flat circular disk at z = 0.
            %
            %   [nodes, panels, nverts] = generate_waterplane_mesh(r_wp, N_th, N_r, x_sym)
            %
            %   HAMS requires a waterplane lid mesh when irregular frequency
            %   removal is enabled (If_remove_irr_freq = 1). The lid is a
            %   flat mesh at z = 0 covering the interior waterplane area.
            %
            %   Panel normals point DOWNWARD (into the fluid interior),
            %   consistent with the lid method convention.
            %
            %   Symmetry must match HullMesh.pnl: Y-symmetry always ON,
            %   X-symmetry optional.
            %
            %   INPUTS:
            %     r_waterplane  scalar    waterplane radius [m]
            %     N_theta_half  integer   azimuthal panels in y ≥ 0 half
            %     N_radial      integer   number of radial divisions
            %     use_x_symmetry  bool    if true, mesh first quadrant only
            %
            %   Reference: Liu (2019), JMSE 7(3):81
            %     "additional 10×20 panels at the interior waterplane"

            if nargin < 4; use_x_symmetry = false; end
            if nargin < 3; N_radial = 8; end

            % Azimuthal range
            if use_x_symmetry
                theta = linspace(0, pi/2, N_theta_half + 1);
            else
                theta = linspace(0, pi, N_theta_half + 1);
            end
            N_theta = length(theta) - 1;

            % Radial stations (uniform, pole to edge)
            r_stations = linspace(0, r_waterplane, N_radial + 1);

            % ── Build nodes ──
            N_per_ring = N_theta + 1;
            n_rings = N_radial;
            n_nodes = 1 + n_rings * N_per_ring;

            nodes = zeros(n_nodes, 3);
            ni = 0;

            % Pole node
            ni = ni + 1;
            nodes(ni, :) = [0, 0, 0];

            % Ring nodes
            for k = 1:n_rings
                rk = r_stations(k + 1);
                for j = 1:N_per_ring
                    ni = ni + 1;
                    nodes(ni, :) = [rk * cos(theta(j)), ...
                                    rk * sin(theta(j)), 0];
                end
            end

            % ── Helper: node index in ring k, position j ──
            ring_node = @(k, j) 1 + (k - 1) * N_per_ring + j;

            % ── Build panels ──
            n_panels = N_theta + (n_rings - 1) * N_theta;
            panels = zeros(n_panels, 4);
            nverts = zeros(n_panels, 1);
            pi_idx = 0;

            % Pole cap: triangles (pole → ring 1)
            % Winding order chosen so that cross(v2-v1, v3-v1) points DOWN
            for j = 1:N_theta
                pi_idx = pi_idx + 1;
                v1 = 1;                  % pole
                v2 = ring_node(1, j+1);  % reversed order for downward normal
                v3 = ring_node(1, j);
                panels(pi_idx, :) = [v1, v2, v3, 0];
                nverts(pi_idx) = 3;
            end

            % Body: quads (ring k → ring k+1)
            % Winding order for downward normal
            for k = 1:(n_rings - 1)
                for j = 1:N_theta
                    pi_idx = pi_idx + 1;
                    v1 = ring_node(k,   j);
                    v2 = ring_node(k,   j+1);
                    v3 = ring_node(k+1, j+1);
                    v4 = ring_node(k+1, j);
                    panels(pi_idx, :) = [v1, v2, v3, v4];
                    nverts(pi_idx) = 4;
                end
            end

            % ── Verify normal direction: must be UPWARD (nz > 0) ────────
            % Scan for first non-degenerate panel to avoid a silent no-op
            % on a degenerate triangle at index N_theta+1.
            % Formula: HAMS CalTransNormals convention —
            %   Quad: cross(V3-V1, V4-V2)   Tri: cross(V2-V1, V3-V2)
            flip_wp = false;
            for test_p = 1:n_panels
                v   = panels(test_p, 1:nverts(test_p));
                pts = nodes(v, :);
                if nverts(test_p) == 4
                    nrm = cross(pts(3,:)-pts(1,:), pts(4,:)-pts(2,:));
                else
                    nrm = cross(pts(2,:)-pts(1,:), pts(3,:)-pts(2,:));
                end
                if norm(nrm) > 1e-14
                    flip_wp = (nrm(3) < 0);   % flip if pointing DOWN (wrong)
                    break;
                end
            end
            if flip_wp
                warning('HAMS_Pipeline:WPNormals', ...
                    'Waterplane normals point down; reversing to UPWARD.');
                for p = 1:n_panels
                    if nverts(p) == 4
                        panels(p, 1:4) = panels(p, [1 4 3 2]);
                    else
                        panels(p, [2 3]) = panels(p, [3 2]);
                    end
                end
            end

            fprintf('  Waterplane mesh: %d nodes, %d panels (%d tri + %d quad)\n', ...
                n_nodes, n_panels, sum(nverts==3), sum(nverts==4));
        end


        %% ═══════════════════════════════════════════════════════════
        %%  WAMIT .3 FILE PARSER (exciting forces)
        %% ═══════════════════════════════════════════════════════════

        function fe_data = parse_wamit_3_file(filepath, ref_body_length, output_freq_type)
            % PARSE_WAMIT_3_FILE  Read HAMS output in WAMIT .3 format.
            %
            %   fe_data = parse_wamit_3_file(filepath)
            %   fe_data = parse_wamit_3_file(filepath, L, freq_type)
            %
            %   CRITICAL: first column interpretation depends on freq_type.
            %     freq_type = 3 -> omega [rad/s]  (DEFAULT)
            %     freq_type = 4 -> period [s]
            %
            %   .3 file format (7-col, auto-detected):
            %     freq  heading  i  |Fe_bar|  phase  Re(Fe_bar)  Im(Fe_bar)
            %
            %   Denormalization:
            %     Fe(i) = Fe_bar(i) * rho * g * L^m
            %     m = 2 for translation (i = 1,2,3)
            %     m = 3 for rotation    (i = 4,5,6)
            %
            %   OUTPUT (ALL DIMENSIONAL, SI units):
            %     .omega      [Nx1]    angular freq [rad/s], ASCENDING
            %     .periods    [Nx1]    wave periods [s]
            %     .Fe         [6xN]    complex exciting force [N or N*m]
            %     .Fe_mag     [6xN]    |Fe| magnitude
            %     .Fe_phase   [6xN]    phase [deg]
            %     .headings   [Hx1]    wave headings [deg]

            if nargin < 2 || isempty(ref_body_length)
                ref_body_length = 1.0;
            end
            if nargin < 3 || isempty(output_freq_type)
                output_freq_type = 3;
            end

            L = ref_body_length;
            rho = HAMS_Pipeline.RHO_WATER;
            grav = HAMS_Pipeline.G;

            assert(exist(filepath, 'file') == 2, ...
                'WAMIT .3 file not found: %s', filepath);

            % -- Read all numeric lines --
            fid = fopen(filepath, 'r');
            raw_data = [];
            while ~feof(fid)
                line = fgetl(fid);
                if ~ischar(line) || isempty(strtrim(line))
                    continue;
                end
                vals = sscanf(line, '%f');
                if length(vals) >= 6
                    raw_data(end+1, 1:length(vals)) = vals(:)'; %#ok<AGROW>
                end
            end
            fclose(fid);

            if isempty(raw_data)
                warning('HAMS_Pipeline:Empty3File', 'No data in .3 file: %s', filepath);
                fe_data = struct('omega', [], 'Fe', []);
                return;
            end

            n_cols = size(raw_data, 2);

            % -- Auto-detect 6-col vs 7-col format --
            col2_vals = unique(raw_data(:, 2));
            if n_cols >= 7 && any(col2_vals > 6 | col2_vals < 1 | ...
                    mod(col2_vals, 1) ~= 0)
                heading_col = 2; mode_col = 3;
                mag_col = 4; phase_col = 5; re_col = 6; im_col = 7;
                fmt = 7;
            else
                heading_col = 0; mode_col = 2;
                mag_col = 3; phase_col = 4; re_col = 5; im_col = 6;
                fmt = 6;
            end

            all_col1  = raw_data(:, 1);
            all_modes = round(raw_data(:, mode_col));
            if heading_col > 0
                all_headings = raw_data(:, heading_col);
            else
                all_headings = zeros(size(all_col1));
            end

            reg_mask = all_col1 > 0;
            unique_col1 = unique(all_col1(reg_mask));

            if output_freq_type == 3
                omega_all = unique_col1;
            elseif output_freq_type == 4
                omega_all = 2 * pi ./ unique_col1;
            end

            [omega_sorted, sort_idx] = sort(omega_all, 'ascend');
            col1_sorted = unique_col1(sort_idx);
            n_freq = length(omega_sorted);
            unique_headings = unique(all_headings(reg_mask));

            m_exp = [2, 2, 2, 3, 3, 3];

            Fe_complex = zeros(6, n_freq);
            Fe_mag     = zeros(6, n_freq);
            Fe_phase   = zeros(6, n_freq);

            target_heading = unique_headings(1);

            for idx = 1:size(raw_data, 1)
                c1 = all_col1(idx);
                if c1 <= 0; continue; end
                heading = all_headings(idx);
                if abs(heading - target_heading) > 1e-6; continue; end
                i = all_modes(idx);
                if i < 1 || i > 6; continue; end
                f_idx = find(col1_sorted == c1, 1);
                if isempty(f_idx); continue; end

                scale = rho * grav * L^m_exp(i);
                Fe_mag(i, f_idx)   = raw_data(idx, mag_col) * scale;
                Fe_phase(i, f_idx) = raw_data(idx, phase_col);
                re_val = raw_data(idx, re_col) * scale;
                im_val = raw_data(idx, im_col) * scale;
                Fe_complex(i, f_idx) = re_val + 1i * im_val;
            end

            fe_data.omega    = omega_sorted;
            fe_data.periods  = 2 * pi ./ omega_sorted;
            fe_data.Fe       = Fe_complex;
            fe_data.Fe_mag   = Fe_mag;
            fe_data.Fe_phase = Fe_phase;
            fe_data.headings = unique_headings;
            fe_data.L        = L;
            fe_data.output_freq_type = output_freq_type;
            fe_data.format_detected = fmt;

            fprintf('  Parsed %s: %d frequencies, %d headings (format: %d-col)\n', ...
                filepath, n_freq, length(unique_headings), fmt);
            nonzero3 = Fe_mag(3, Fe_mag(3,:) > 0);
            if ~isempty(nonzero3)
                fprintf('    |Fe3| range: [%.1f, %.1f] N  [dimensional]\n', ...
                    min(nonzero3), max(nonzero3));
            end
        end


        %% ═══════════════════════════════════════════════════════════
        %%  CROSS-VALIDATION (HAMS format vs WAMIT format)
        %% ═══════════════════════════════════════════════════════════

        function report = cross_validate_output(run_dir, ref_body_length, output_freq_type)
            % CROSS_VALIDATE_OUTPUT  Compare HAMS-format vs WAMIT-format.
            %
            %   report = cross_validate_output(run_dir)
            %   report = cross_validate_output(run_dir, L, freq_type)
            %
            % ── DEAD CODE ─────────────────────────────────────────────────────
            % Requires HAMS_Postprocess.read_all_dofs which is NOT present in
            % this codebase. Calling this function will throw an "Undefined
            % function" error at the HAMS_Postprocess call. Kept for reference.
            % ──────────────────────────────────────────────────────────────────
            %
            %   PASS criterion: max relative error < 1e-4.

            if nargin < 2 || isempty(ref_body_length); ref_body_length = 1.0; end
            if nargin < 3 || isempty(output_freq_type); output_freq_type = 3; end

            fprintf('\n=== Cross-validation: HAMS vs WAMIT format ===\n');

            hams_dir = fullfile(run_dir, 'Output', 'Hams_format');
            assert(exist(hams_dir, 'dir') == 7, 'Hams_format dir not found: %s', hams_dir);
            hams = HAMS_Postprocess.read_all_dofs(hams_dir);

            wamit_dir = fullfile(run_dir, 'Output', 'Wamit_format');
            one_files = dir(fullfile(wamit_dir, '*.1'));
            assert(~isempty(one_files), 'No .1 file in: %s', wamit_dir);
            wamit = HAMS_Pipeline.parse_wamit_1_file(...
                fullfile(wamit_dir, one_files(1).name), ref_body_length, output_freq_type);

            % Compare A(inf)
            err_Ainf = abs(hams.A_inf - wamit.A_inf);
            denom = max(abs(hams.A_inf), abs(wamit.A_inf));
            denom(denom < 1e-12) = 1;
            report.max_err_Ainf = max(err_Ainf(:) ./ denom(:));
            fprintf('  A(inf): max rel error = %.2e\n', report.max_err_Ainf);

            % Compare A(omega) and B(omega) at matching frequencies
            max_err_A = 0; max_err_B = 0; n_compared = 0;
            for kw = 1:length(wamit.omega)
                [md, kh] = min(abs(hams.omega - wamit.omega(kw)));
                if md > 0.01; continue; end
                n_compared = n_compared + 1;
                for ii = 1:6
                    for jj = 1:6
                        d = max(abs(hams.A(ii,jj,kh)), abs(wamit.A(ii,jj,kw)));
                        if d > 1e-10
                            max_err_A = max(max_err_A, abs(hams.A(ii,jj,kh)-wamit.A(ii,jj,kw))/d);
                        end
                        d = max(abs(hams.B(ii,jj,kh)), abs(wamit.B(ii,jj,kw)));
                        if d > 1e-10
                            max_err_B = max(max_err_B, abs(hams.B(ii,jj,kh)-wamit.B(ii,jj,kw))/d);
                        end
                    end
                end
            end
            report.max_err_A = max_err_A;
            report.max_err_B = max_err_B;
            report.n_freq_compared = n_compared;
            fprintf('  A(omega): max rel error = %.2e (%d freq)\n', max_err_A, n_compared);
            fprintf('  B(omega): max rel error = %.2e\n', max_err_B);

            tol = 1e-4;
            report.passed = report.max_err_Ainf < tol && max_err_A < tol && max_err_B < tol;
            if report.passed
                fprintf('  CROSS-VALIDATION PASSED (all < %.0e)\n', tol);
            else
                fprintf('  CROSS-VALIDATION FAILED\n');
            end
        end


        %% ═══════════════════════════════════════════════════════════
        %%  PNL FILE READER
        %% ═══════════════════════════════════════════════════════════

        function [n_panels, n_nodes, nodes, panels, nverts, symmetry] = ...
                read_pnl_file(filepath)
            % READ_PNL_FILE  Parse HAMS .pnl mesh file.
            %
            %   [np, nn, nodes, panels, nverts, sym] = read_pnl_file(path)

            fid = fopen(filepath, 'r');
            raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
            fclose(fid);
            all_lines = raw{1};
            n_lines = length(all_lines);

            header_line = 0; node_start = 0; node_end = 0;
            relation_start = 0; relation_end = 0;

            for il = 1:n_lines
                ln = strtrim(all_lines{il});
                ln_low = lower(ln);
                if header_line == 0 && ~isempty(ln) && ln(1) ~= '#' && ln(1) ~= '-'
                    vals = sscanf(ln, '%d');
                    if length(vals) == 4; header_line = il; end
                end
                if ~isempty(ln) && ln(1) == '#'
                    if contains(ln_low,'start') && contains(ln_low,'node coord')
                        node_start = il + 1;
                    end
                    if contains(ln_low,'end') && contains(ln_low,'node coord')
                        node_end = il - 1;
                    end
                    if contains(ln_low,'start') && contains(ln_low,'node relation')
                        relation_start = il + 1;
                    end
                    if contains(ln_low,'end') && contains(ln_low,'node relation')
                        relation_end = il - 1;
                    end
                end
            end

            assert(header_line > 0, 'No header in %s', filepath);
            vals = sscanf(all_lines{header_line}, '%d');
            n_panels = vals(1); n_nodes = vals(2);
            symmetry = [vals(3), vals(4)];

            nodes = zeros(n_nodes, 3);
            for il = node_start:node_end
                ln = strtrim(all_lines{il});
                if isempty(ln) || ln(1)=='#' || ln(1)=='-'; continue; end
                vals = sscanf(ln, '%f');
                if length(vals) >= 4
                    nid = round(vals(1));
                    if nid >= 1 && nid <= n_nodes
                        nodes(nid,:) = vals(2:4)';
                    end
                end
            end

            panels = zeros(n_panels, 4);
            nverts = zeros(n_panels, 1);
            pc = 0;
            for il = relation_start:relation_end
                ln = strtrim(all_lines{il});
                if isempty(ln) || ln(1)=='#' || ln(1)=='-'; continue; end
                vals = sscanf(ln, '%d');
                if length(vals) >= 5
                    pc = pc + 1;
                    nv_p = vals(2);
                    nverts(pc) = nv_p;
                    panels(pc, 1:nv_p) = vals(3:2+nv_p)';
                end
            end

            fprintf('  Read %s: %d nodes, %d panels, sym=[%d,%d]\n', ...
                filepath, n_nodes, n_panels, symmetry(1), symmetry(2));
        end


        %% ═══════════════════════════════════════════════════════════
        %%  GENERAL WATERPLANE MESH (non-circular contours)
        %% ═══════════════════════════════════════════════════════════

        function [wp_nodes, wp_panels, wp_nverts] = ...
                generate_waterplane_mesh_structured(contour_pts, target_edge, z_wl)
            % GENERATE_WATERPLANE_MESH_STRUCTURED  Structured quad lid mesh.
            %
            %   [nodes, panels, nverts] = generate_waterplane_mesh_structured(
            %       contour_pts, target_edge, z_wl)
            %
            %   INPUTS
            %     contour_pts  [N × 2]  ORDERED boundary polygon (chain walk)
            %     target_edge  scalar   max panel edge length [m] (default: 0.4)
            %     z_wl         scalar   waterplane elevation [m] (default: 0)
            %
            %   ALGORITHM: Two-curve transfinite interpolation (TFI)
            %     1. Find long axis (two farthest boundary vertices).
            %     2. Split boundary into left/right curves at those points.
            %     3. Arc-length parameterise both curves.
            %     4. Structured grid: node(k,j) = (1-s_j)*L(t_k) + s_j*R(t_k).
            %     5. Tips (t=0, t=1) → triangle fan.  Interior → quads.
            %     6. Enforce downward normals on all panels.
            %
            %   WHY structured quads, not CDT?
            %     The previous CDT mesher used ONLY boundary vertices — zero
            %     interior refinement.  MATLAB's delaunayTriangulation does
            %     not add Steiner points.  Result: 88 triangles spanning the
            %     entire domain (max edge 2.96 m, aspect ratio 44:1), unable
            %     to suppress the 17 interior Dirichlet eigenfrequencies
            %     within the HAMS frequency grid.  The structured TFI
            %     approach produces ~136 panels with max edge ≤ target_edge,
            %     adequate for irregular frequency removal at any aspect ratio.
            %
            %   REPLACES: generate_waterplane_mesh_general (v1.5)
            %
            %   VERIFIED: 2026-03-24, 8 diagnostic checks pass.
            %     Max edge 0.39 m (target 0.4 m), max aspect 4.0,
            %     boundary seal gap 1.1e-16 m, all normals DOWN.
            %
            %   See also: write_panelizer_wp_pnl, run_single_hams

            if nargin < 2 || isempty(target_edge); target_edge = 0.4; end
            if nargin < 3 || isempty(z_wl);        z_wl = 0; end

            xy = contour_pts(:, 1:2);
            N  = size(xy, 1);
            if N < 4
                wp_nodes = zeros(0,3); wp_panels = zeros(0,4);
                wp_nverts = zeros(0,1); return;
            end

            % ── Step 1: Long axis (two farthest boundary vertices) ─────
            D_max = 0; iA = 1; iB = 1;
            for i = 1:N; for j = i+1:N
                d = norm(xy(i,:) - xy(j,:));
                if d > D_max; D_max = d; iA = i; iB = j; end
            end; end
            if iA > iB; tmp = iA; iA = iB; iB = tmp; end

            % ── Step 2: Split boundary into two curves ─────────────────
            idx_fwd = (iA:iB)';
            idx_bwd = [iA:-1:1, N:-1:iB]';
            c_fwd = xy(idx_fwd, :);
            c_bwd = xy(idx_bwd, :);

            % Left/right by cross product with long axis
            ax_dir = xy(iB,:) - xy(iA,:);
            mid_fwd = c_fwd(round(size(c_fwd,1)/2), :) - xy(iA,:);
            cz = ax_dir(1)*mid_fwd(2) - ax_dir(2)*mid_fwd(1);
            if cz >= 0; cL = c_fwd; cR = c_bwd;
            else;       cL = c_bwd; cR = c_fwd; end

            % ── Step 3: Arc-length parameterisation ────────────────────
            sL = HAMS_Pipeline.wp_arc_len(cL);
            sR = HAMS_Pipeline.wp_arc_len(cR);

            % ── Step 4: Grid dimensions ────────────────────────────────
            half_perim = max(sL(end), sR(end));
            max_width = 0;
            for k = 1:50
                t_probe = (k-1)/49;
                pL = HAMS_Pipeline.wp_interp(cL, sL, t_probe);
                pR = HAMS_Pipeline.wp_interp(cR, sR, t_probe);
                max_width = max(max_width, norm(pR - pL));
            end

            Nt = max(6, ceil(half_perim / target_edge) + 1);
            Ns = max(4, ceil(max_width  / target_edge) + 1);
            if mod(Ns, 2) == 0; Ns = Ns + 1; end  % centre column

            t_grid = linspace(0, 1, Nt);
            s_grid = linspace(0, 1, Ns);

            % ── Step 5: Generate grid nodes ────────────────────────────
            gxy = zeros(Nt, Ns, 2);
            for k = 1:Nt
                pL = HAMS_Pipeline.wp_interp(cL, sL, t_grid(k));
                pR = HAMS_Pipeline.wp_interp(cR, sR, t_grid(k));
                for j = 1:Ns
                    gxy(k, j, :) = (1 - s_grid(j)) * pL + s_grid(j) * pR;
                end
            end

            % ── Step 6: Flatten nodes (tip A, interior rows, tip B) ────
            n_int_rows = Nt - 2;
            n_nodes = 2 + n_int_rows * Ns;
            nxy = zeros(n_nodes, 2);
            nxy(1, :) = squeeze(gxy(1, 1, :))';       % tip A
            for k = 2:Nt-1
                rs = 2 + (k-2) * Ns;
                for j = 1:Ns
                    nxy(rs + j - 1, :) = squeeze(gxy(k, j, :))';
                end
            end
            nxy(n_nodes, :) = squeeze(gxy(Nt, 1, :))'; % tip B

            ri = @(k, j) 2 + (k-2)*Ns + (j-1);  % row k (2-based), col j (1-based)

            % ── Step 7: Build panels ───────────────────────────────────
            n_tip = Ns - 1;
            n_quads = max(0, Nt - 3) * (Ns - 1);
            n_panels = 2 * n_tip + n_quads;
            panels = zeros(n_panels, 4);
            nverts = zeros(n_panels, 1);
            pi_idx = 0;

            % Tip A fan (triangles)
            for j = 1:n_tip
                pi_idx = pi_idx + 1;
                panels(pi_idx, :) = [1, ri(2, j+1), ri(2, j), ri(2, j)];
                nverts(pi_idx) = 3;
            end
            % Interior quads
            for k = 2:Nt-2
                for j = 1:Ns-1
                    pi_idx = pi_idx + 1;
                    panels(pi_idx, :) = [ri(k,j), ri(k,j+1), ri(k+1,j+1), ri(k+1,j)];
                    nverts(pi_idx) = 4;
                end
            end
            % Tip B fan (triangles)
            kl = Nt - 1;
            for j = 1:n_tip
                pi_idx = pi_idx + 1;
                panels(pi_idx, :) = [ri(kl, j), ri(kl, j+1), n_nodes, n_nodes];
                nverts(pi_idx) = 3;
            end

            % ── Step 8: 3D nodes, enforce normals UPWARD (+Z) ──────────
            % HAMS CalTransNormals convention:
            %   Quad: cross(V3-V1, V4-V2)   Tri: cross(V2-V1, V3-V2)
            % Scan for first non-degenerate panel.
            wp_nodes = [nxy, repmat(z_wl, n_nodes, 1)];

            flip_wp2 = false;
            for ti_scan = 1:n_panels
                vv2 = panels(ti_scan, :); nv2 = nverts(ti_scan);
                if nv2 == 4
                    nrm2 = cross(wp_nodes(vv2(3),:)-wp_nodes(vv2(1),:), ...
                                 wp_nodes(vv2(4),:)-wp_nodes(vv2(2),:));
                else
                    nrm2 = cross(wp_nodes(vv2(2),:)-wp_nodes(vv2(1),:), ...
                                 wp_nodes(vv2(3),:)-wp_nodes(vv2(2),:));
                end
                if norm(nrm2) > 1e-14
                    flip_wp2 = (nrm2(3) < 0);  % flip if DOWN (wrong)
                    break;
                end
            end
            if flip_wp2
                for p = 1:n_panels
                    if nverts(p) == 4
                        panels(p, 1:4) = panels(p, [1 4 3 2]);
                    else
                        panels(p, [2 3]) = panels(p, [3 2]);
                        panels(p, 4) = panels(p, 3);
                    end
                end
            end

            wp_panels = panels;
            wp_nverts = nverts;

            nq = sum(nverts == 4); nt = sum(nverts == 3);
            fprintf('  WP structured mesh: %d nodes, %d panels (%d quads + %d tri), z=%.3f\n', ...
                n_nodes, n_panels, nq, nt, z_wl);
        end


        function [wp_nodes, wp_panels, wp_nverts] = ...
                generate_wp_from_mesh_waterline(mesh, target_edge)
            % GENERATE_WP_FROM_MESH_WATERLINE  v2 - Topological walker + Blossom-Quad.
            %
            %   [wp_nodes, wp_panels, wp_nverts] = ...
            %       HAMS_Pipeline.generate_wp_from_mesh_waterline(mesh)
            %   [wp_nodes, wp_panels, wp_nverts] = ...
            %       HAMS_Pipeline.generate_wp_from_mesh_waterline(mesh, target_edge)
            %
            %   BUG FIX (B3): v1 used abs(z)<z_tol vertex filter, which for
            %   C0 at draft=1 m captured ALL 625 platform nodes (flat at z=0).
            %   atan2 sort of 625 mixed interior+boundary points produces a
            %   self-intersecting polygon -> corrupted WP mesh -> IRSP=1 fails
            %   -> irregular frequency suppression fails -> A15 != A51 spikes.
            %
            %   v2 uses topological open-edge walk to extract ONLY the ~21
            %   true boundary nodes -> correct polygon -> Blossom-Quad mesh.

            if nargin < 2 || isempty(target_edge), target_edge = 0.15; end

            x_sym = 0; y_sym = 0;
            if isfield(mesh,'x_sym'), x_sym = mesh.x_sym; end
            if isfield(mesh,'y_sym'), y_sym = mesh.y_sym; end

            % ── Phase 1: Topological boundary extraction ─────────────
            [boundary_xy, n_loops] = HAMS_Pipeline.extract_waterline_boundary_topo( ...
                mesh, 1e-3, target_edge);

            if isempty(boundary_xy) || size(boundary_xy,1) < 3
                warning('HAMS_Pipeline:WPNoWaterline', ...
                    'No z=0 open-edge boundary found. Check mesh trim_at_wl.');
                wp_nodes  = zeros(0,3);
                wp_panels = zeros(0,4);
                wp_nverts = zeros(0,1);
                return;
            end

            fprintf('  WP boundary: %d nodes (%d loop) via topological edge walk\n', ...
                size(boundary_xy,1), n_loops);

            % ── Phase 2: Blossom-Quad interior mesh ──────────────────
            [wp_nodes, wp_panels, wp_nverts] = HAMS_Pipeline.mesh_wp_blossomquad( ...
                boundary_xy, target_edge, x_sym, y_sym, 0);

            n_q = sum(wp_nverts == 4); n_t = sum(wp_nverts == 3);
            Aw  = abs(polyarea(boundary_xy(:,1), boundary_xy(:,2)));
            fprintf('  WP Blossom-Quad: %d nodes, %d panels (%d quad + %d tri)\n', ...
                size(wp_nodes,1), numel(wp_nverts), n_q, n_t);
            fprintf('    WP area %.4f m^2 | target edge %.3f m\n', Aw, target_edge);

            n_up = HAMS_Pipeline.count_wp_normals_up(wp_nodes, wp_panels, wp_nverts);
            if n_up ~= numel(wp_nverts)
                warning('HAMS_Pipeline:WPNormalBad', ...
                    '%d/%d WP panels have wrong (downward) normal — expected UPWARD.', ...
                    numel(wp_nverts)-n_up, numel(wp_nverts));
            end
        end


        function [boundary_xy, n_loops] = extract_waterline_boundary_topo( ...
                mesh, z_tol, target_edge_hint)
            % EXTRACT_WATERLINE_BOUNDARY_TOPO  Ordered waterline polygon via edge walk.
            %
            %   Finds hull mesh edges that are OPEN (belong to exactly 1 panel)
            %   AND have both endpoints at |z| < z_tol. These are the waterline
            %   boundary edges. Interior platform nodes are never open edges and
            %   are therefore never captured -- fixes the C0 platform node flood.
            %
            %   Chains boundary edges into an ordered polygon. For half-body
            %   (y_sym=1), adds a straight closure segment along y=0.

            if nargin < 2 || isempty(z_tol),           z_tol = 1e-3; end
            if nargin < 3 || isempty(target_edge_hint), target_edge_hint = 0.15; end

            x_sym = 0; y_sym = 0;
            if isfield(mesh,'x_sym'), x_sym = mesh.x_sym; end
            if isfield(mesh,'y_sym'), y_sym = mesh.y_sym; end

            verts = mesh.vertices;
            n_p   = mesh.n_panels;

            % ── Build edge -> count map ────────────────────────────────
            edge_count = containers.Map('KeyType','char','ValueType','int32');
            edge_ep    = containers.Map('KeyType','char','ValueType','any');

            for p = 1:n_p
                v  = mesh.panels(p,:);
                nv = 4; if v(3)==v(4), nv=3; end
                for e = 1:nv
                    va = v(e); vb = v(mod(e,nv)+1);
                    key = sprintf('%d_%d', min(va,vb), max(va,vb));
                    if edge_count.isKey(key)
                        edge_count(key) = edge_count(key) + int32(1);
                    else
                        edge_count(key) = int32(1);
                        edge_ep(key)    = [min(va,vb), max(va,vb)];
                    end
                end
            end

            % ── Collect open edges at z=0 ─────────────────────────────
            bnd_a = []; bnd_b = [];
            ks = edge_count.keys();
            for k = 1:numel(ks)
                if edge_count(ks{k}) ~= 1, continue; end
                ev = edge_ep(ks{k});
                if abs(verts(ev(1),3))<z_tol && abs(verts(ev(2),3))<z_tol
                    bnd_a(end+1) = ev(1); %#ok<AGROW>
                    bnd_b(end+1) = ev(2); %#ok<AGROW>
                end
            end

            if isempty(bnd_a)
                boundary_xy = []; n_loops = 0; return;
            end

            % ── Build adjacency list ───────────────────────────────────
            all_bv = unique([bnd_a, bnd_b]);
            n_bv   = numel(all_bv);
            mx     = max(all_bv)+1;
            g2l    = zeros(mx,1,'int32');
            for i = 1:n_bv, g2l(all_bv(i)) = int32(i); end
            adj = cell(n_bv,1);
            for e = 1:numel(bnd_a)
                la=g2l(bnd_a(e)); lb=g2l(bnd_b(e));
                adj{la}(end+1)=lb; adj{lb}(end+1)=la;
            end

            % ── Walk loops ────────────────────────────────────────────
            visited = false(n_bv,1);
            loops   = {};
            for sv = 1:n_bv
                if visited(sv), continue; end
                loop=[sv]; visited(sv)=true; prev=-1; cur=sv;
                while true
                    nbrs=adj{cur}; unv=nbrs(~visited(nbrs));
                    if isempty(unv), break; end
                    if prev>0, unv=unv(unv~=prev); end
                    if isempty(unv), break; end
                    prev=cur; cur=unv(1); visited(cur)=true; loop(end+1)=cur; %#ok<AGROW>
                end
                if numel(loop)>=3
                    loops{end+1} = verts(all_bv(loop),1:2); %#ok<AGROW>
                end
            end

            n_loops = numel(loops);
            if n_loops==0, boundary_xy=[]; return; end

            % Sort by area, take outer (largest)
            ar = cellfun(@(L) abs(polyarea(L(:,1),L(:,2))), loops);
            [~,si] = sort(ar,'descend');
            boundary_xy = loops{si(1)};

            % ── Half-body y=0 closure ─────────────────────────────────
            if y_sym==1 && x_sym==0
                boundary_xy(1,2)=0; boundary_xy(end,2)=0;
                p1=boundary_xy(1,:); p2=boundary_xy(end,:);
                span=abs(p1(1)-p2(1));
                nc=max(0, ceil(span/target_edge_hint)-1);
                if nc>0
                    xc=linspace(max(p1(1),p2(1)), min(p1(1),p2(1)), nc+2)';
                    xc=xc(2:end-1);
                    boundary_xy=[boundary_xy; xc, zeros(numel(xc),1)];
                    fprintf('    Half-body WP: y=0 closure %.2f m (%d pts)\n', span, nc);
                end
            elseif x_sym==1 && y_sym==0
                boundary_xy(1,1)=0; boundary_xy(end,1)=0;
            end

            % Remove duplicate endpoint if closed polygon supplied
            if norm(boundary_xy(1,:)-boundary_xy(end,:)) < 1e-10
                boundary_xy=boundary_xy(1:end-1,:);
            end
        end


        function boundary_xy = clean_isocurve_polygon(raw_xy, target_edge)
            % CLEAN_ISOCURVE_POLYGON  Fix non-simple polygon from extract_isocurve_at_z.
            %
            % PROBLEM
            %   extract_isocurve_at_z concatenates arcs from source + mirror surfaces
            %   in surface-table order, not angular order.  For a 4-surface body (C0),
            %   the four quarter-arcs are appended as Q1→Q4→Q2→Q3 (a figure-8, not a
            %   circle), producing a non-simple polygon where:
            %     · abs(polyarea) = 0  (signed areas of crossing loops cancel)
            %     · MATLAB CDT warns "Intersecting edge constraints have been split"
            %     · CDT inserts Steiner points → 823 nodes, 650 panels for a disk
            %
            % FIX — four steps
            %   1. Sort by azimuthal angle around centroid → simple polygon for any
            %      star-convex waterplane (all RevSurf WEC hulls)
            %   2. Remove near-duplicate seam points (consecutive pts ≤ 1e-5 m apart)
            %   3. Validate: if polyarea still ≈ 0, return raw input so caller can warn
            %   4. Downsample to target_edge/4 boundary resolution.
            %      extract_isocurve_at_z may return n_u≈200 pts/quadrant = 800 pts.
            %      resample_polygon (inside blossomquad) adds pts for long edges but
            %      never removes pts for short ones → 800 boundary points produce
            %      823 CDT nodes, 650 panels for a disk needing ~8 panels.
            %      Criterion: max(48, ceil(perimeter/(target_edge/4))).
            %        · n=48 for r=0.37m, target=0.30m → Aw error = 0.27%
            %        · n=84 for r=1.0m  target=0.30m → Aw error = 0.09%
            %      Using target_edge alone gives n=8 → octagon, Aw error = 10%.

            if isempty(raw_xy) || size(raw_xy,1) < 3
                boundary_xy = raw_xy;
                return;
            end

            xy = raw_xy(:, 1:2);

            % Step 1: sort by atan2 around centroid
            cx = mean(xy(:,1));
            cy = mean(xy(:,2));
            [~, ang_ord] = sort(atan2(xy(:,2)-cy, xy(:,1)-cx));
            xy = xy(ang_ord, :);

            % Step 2: drop near-duplicate seam points
            dists = sqrt(sum(diff([xy; xy(1,:)], 1, 1).^2, 2));
            keep  = dists > 1e-5;
            if sum(keep) >= 3
                xy = xy(keep, :);
            end

            % Step 3: validate
            Aw = abs(polyarea(xy(:,1), xy(:,2)));
            if Aw < 1e-8
                boundary_xy = raw_xy(:, 1:2);   % still bad — let caller warn
                return;
            end

            % Step 4: downsample to target_edge resolution
            % Criterion: target_edge/4 (not target_edge) because mesh_wp_blossomquad
            % calls resample_polygon which adds interior points but never removes
            % boundary points.  Using target_edge directly gives n≈8 for a small
            % WP (r=0.37m), producing a regular octagon with ~10% area error that
            % corrupts C33 and C55.  target_edge/4 gives n≈32-48, error <0.3%.
            r_est    = sqrt(Aw / pi);
            n_needed = max(48, ceil(2*pi*r_est / (target_edge/4)));
            n_pts    = size(xy, 1);
            if n_pts > n_needed * 2
                idx = round(linspace(1, n_pts, n_needed));
                xy  = xy(unique(idx, 'stable'), :);
            end

            boundary_xy = xy;
        end

        function [wp_nodes, wp_panels, wp_nverts] = mesh_wp_blossomquad( ...
                boundary_xy, target_edge, x_sym, y_sym, z_wl, lock_boundary)
            % MESH_WP_BLOSSOMQUAD  CDT + greedy quad matching (Blossom-Quad).
            %   quad coverage -- adequate for BEM irregular frequency removal.
            %
            %   HAMS normal convention (NormalProcess.f90, CalTransNormals):
            %     Quad: cross(v4-v2, v3-v1) -> nz<0 for CCW winding at z=0.
            %     Tri:  cross(v1-v2, v3-v2) -> nz<0 for CCW winding.
            %
            %   lock_boundary (optional, default false):
            %     When true, the outer boundary nodes are used verbatim — no
            %     resample_polygon call.  Use this when boundary_xy is sourced
            %     directly from hull mesh waterline vertices (hull_waterline_
            %     polygon) so that the WP outer boundary matches the hull mesh
            %     node positions exactly.  Interior Steiner points are still
            %     added from the meshgrid; only the boundary resampling is
            %     suppressed.  All existing 5-argument call sites are unaffected
            %     (lock_boundary defaults to false).

            if nargin<5 || isempty(z_wl),          z_wl = 0;          end
            if nargin<4 || isempty(y_sym),          y_sym = 0;         end  %#ok<NASGU>
            if nargin<3 || isempty(x_sym),          x_sym = 0;         end  %#ok<NASGU>
            if nargin<6 || isempty(lock_boundary),  lock_boundary = false; end

            xy = boundary_xy(:,1:2);
            if size(xy,1) < 3
                wp_nodes=zeros(0,3); wp_panels=zeros(0,4); wp_nverts=zeros(0,1); return;
            end

            % Resample boundary and add interior Steiner points.
            % When lock_boundary=true the caller has supplied hull mesh nodes
            % directly; resampling the boundary would move those nodes off the
            % hull mesh positions and recreate the boundary mismatch.  The
            % interior Steiner grid (meshgrid below) is always generated.
            if ~lock_boundary
                xy = HAMS_Pipeline.resample_polygon(xy, target_edge);
            end
            N    = size(xy,1);
            step = target_edge;
            [Xg,Yg] = meshgrid( ...
                (min(xy(:,1))+step/2):step:(max(xy(:,1))-step/2), ...
                (min(xy(:,2))+step/2):step:(max(xy(:,2))-step/2));
            pts_s = [Xg(:),Yg(:)];
            in_m  = inpolygon(pts_s(:,1),pts_s(:,2),xy(:,1),xy(:,2));
            pts_s = pts_s(in_m,:);
            % Remove points too close to boundary
            if ~isempty(pts_s)
                md = inf(size(pts_s,1),1);
                for i=1:N
                    j=mod(i,N)+1; ab=xy(j,:)-xy(i,:); lab=norm(ab);
                    if lab<1e-14, continue; end
                    t=max(0,min(1,((pts_s-xy(i,:))*ab')/(lab^2)));
                    d=sqrt(sum((pts_s-xy(i,:)-t.*ab).^2,2));
                    md=min(md,d);
                end
                pts_s=pts_s(md>=target_edge/3,:);
            end

            all_pts = [xy; pts_s];
            N_bnd   = N;
            bnd_c   = [(1:N_bnd-1)',(2:N_bnd)'; N_bnd,1];

            % CDT
            try
                dt = delaunayTriangulation(all_pts(:,1), all_pts(:,2), bnd_c);
            catch ME_dt
                warning('HAMS_Pipeline:CDTFailed','CDT failed: %s. TFI fallback.',ME_dt.message);
                [wp_nodes,wp_panels,wp_nverts] = ...
                    HAMS_Pipeline.generate_waterplane_mesh_structured(boundary_xy,target_edge,z_wl);
                return;
            end

            int_mask = isInterior(dt);
            tris     = dt.ConnectivityList(int_mask,:);
            pts2d    = dt.Points;
            nT       = size(tris,1);
            if nT==0
                wp_nodes=zeros(0,3); wp_panels=zeros(0,4); wp_nverts=zeros(0,1); return;
            end

            % Build dual adjacency
            bnd_set = containers.Map('KeyType','char','ValueType','logical');
            for e=1:size(bnd_c,1)
                bnd_set(sprintf('%d_%d',min(bnd_c(e,:)),max(bnd_c(e,:))))=true;
            end
            e2t = containers.Map('KeyType','char','ValueType','any');
            for ti=1:nT
                v=tris(ti,:);
                ees=[v(1),v(2);v(2),v(3);v(3),v(1)];
                for e=1:3
                    k=sprintf('%d_%d',min(ees(e,:)),max(ees(e,:)));
                    if e2t.isKey(k), e2t(k)=[e2t(k),ti]; else, e2t(k)=ti; end
                end
            end

            cp=[]; ce=[];
            ks3=e2t.keys();
            for k=1:numel(ks3)
                tl=e2t(ks3{k});
                if numel(tl)~=2, continue; end
                if bnd_set.isKey(ks3{k}), continue; end
                nm=sscanf(ks3{k},'%d_%d');
                cp(end+1,:)=[tl(1),tl(2)];  %#ok<AGROW>
                ce(end+1,:)=[nm(1),nm(2)];   %#ok<AGROW>
            end
            nC=size(cp,1);

            % Quality per pair
            cq=zeros(nC,1); ct=zeros(nC,2,'int32');
            for c=1:nC
                ti=cp(c,1); tj=cp(c,2);
                va=ce(c,1); vb=ce(c,2);
                vi=tris(ti,:); vj=tris(tj,:);
                ti_t=vi(vi~=va & vi~=vb); tj_t=vj(vj~=va & vj~=vb);
                if isempty(ti_t)||isempty(tj_t), continue; end
                ti_t=ti_t(1); tj_t=tj_t(1);
                ct(c,:)=int32([ti_t,tj_t]);
                cq(c)=HAMS_Pipeline.quad_scaled_jacobian( ...
                    pts2d(ti_t,:),pts2d(va,:),pts2d(tj_t,:),pts2d(vb,:));
            end

            % Greedy matching
            [~,si2]=sort(cq,'descend');
            matched=false(nT,1); nq=0; qv=zeros(nC,4,'int32');
            for kk=1:nC
                c=si2(kk);
                if cq(c)<0.01, continue; end
                ti=cp(c,1); tj=cp(c,2);
                if matched(ti)||matched(tj), continue; end
                matched(ti)=true; matched(tj)=true;
                va=ce(c,1); vb=ce(c,2);
                ti_t=ct(c,1); tj_t=ct(c,2);
                p1=pts2d(ti_t,:); p2=pts2d(va,:); p3=pts2d(tj_t,:); p4=pts2d(vb,:);
                pts4=[p1;p2;p3;p4];
                sa=0.5*sum(pts4(:,1).*pts4([2:4,1],2)-pts4([2:4,1],1).*pts4(:,2));
                nq=nq+1;
                if sa>=0, qv(nq,:)=int32([ti_t,va,tj_t,vb]);
                else,      qv(nq,:)=int32([ti_t,vb,tj_t,va]); end
            end

            % Assemble
            nt_rem=sum(~matched); pp_tot=nq+nt_rem;
            wp_panels=zeros(pp_tot,4); wp_nverts=zeros(pp_tot,1);
            pp=0;
            for q=1:nq, pp=pp+1; wp_panels(pp,:)=double(qv(q,:)); wp_nverts(pp)=4; end
            for ti=1:nT
                if matched(ti), continue; end
                pp=pp+1; v=tris(ti,:);
                wp_panels(pp,:)=[v(1),v(2),v(3),v(3)]; wp_nverts(pp)=3;
            end
            wp_panels=wp_panels(1:pp,:); wp_nverts=wp_nverts(1:pp);

            % 3D nodes
            wp_nodes=[pts2d, repmat(z_wl,size(pts2d,1),1)];

            % Enforce HAMS nz>0 (UPWARD — verified correct convention).
            % Formula: HAMS CalTransNormals — Quad: cross(V3-V1,V4-V2)
            %                                  Tri:  cross(V2-V1,V3-V2)
            flip_needed=false;
            for p=1:pp
                v=wp_panels(p,:); nv=wp_nverts(p);
                if nv==4, d1=wp_nodes(v(3),:)-wp_nodes(v(1),:); d2=wp_nodes(v(4),:)-wp_nodes(v(2),:);
                else,      d1=wp_nodes(v(2),:)-wp_nodes(v(1),:); d2=wp_nodes(v(3),:)-wp_nodes(v(2),:); end
                nrm=cross(d1,d2);
                if norm(nrm)>1e-14, flip_needed=(nrm(3)<0); break; end  % flip if DOWN (wrong)
            end
            if flip_needed
                for p=1:pp
                    if wp_nverts(p)==4, wp_panels(p,:)=wp_panels(p,[1 4 3 2]);
                    else, wp_panels(p,:)=[wp_panels(p,1),wp_panels(p,3),wp_panels(p,2),wp_panels(p,2)]; end
                end
            end

            fprintf('    CDT+Blossom: %d tri -> %d quad + %d tri (%.0f%% quad)\n', ...
                nT, nq, nt_rem, 100*nq*2/max(nT,1));
        end


        function Q = quad_scaled_jacobian(p1, p2, p3, p4)
            % QUAD_SCALED_JACOBIAN  Min scaled Jacobian over 4 corners [0,1].
            %   Measures worst corner angle: 0=degenerate, 1=perfect square.
            pts=[p1(1:2);p2(1:2);p3(1:2);p4(1:2)]; Q=1.0;
            for k=1:4
                pk=pts(k,:); pp=pts(mod(k-2,4)+1,:); pn=pts(mod(k,4)+1,:);
                e1=pn-pk; e2=pp-pk; l1=norm(e1); l2=norm(e2);
                if l1<1e-14||l2<1e-14, Q=0; return; end
                Q=min(Q, abs(e1(1)*e2(2)-e1(2)*e2(1))/(l1*l2));
            end
        end

        function xy_out = resample_polygon(xy, max_edge)
            % RESAMPLE_POLYGON  Insert midpoints so no edge > max_edge.
            N=size(xy,1); out=zeros(0,2);
            for i=1:N
                j=mod(i,N)+1; pa=xy(i,:); pb=xy(j,:); L=norm(pb-pa);
                out(end+1,:)=pa; %#ok<AGROW>
                if L>max_edge*1.5
                    ns=ceil(L/max_edge);
                    for s=1:ns-1, out(end+1,:)=pa+(s/ns)*(pb-pa); end %#ok<AGROW>
                end
            end
            xy_out=out;
        end

        function n_up = count_wp_normals_up(nodes, panels, nverts)
            % COUNT_WP_NORMALS_UP  Count WP panels with UPWARD normal (nz > 0).
            % Formula: HAMS CalTransNormals — Quad: cross(V3-V1,V4-V2)
            %                                  Tri:  cross(V2-V1,V3-V2)
            n_up=0;
            for p=1:size(panels,1)
                v=panels(p,:); nv=nverts(p);
                if nv==4, d1=nodes(v(3),:)-nodes(v(1),:); d2=nodes(v(4),:)-nodes(v(2),:);
                else,      d1=nodes(v(2),:)-nodes(v(1),:); d2=nodes(v(3),:)-nodes(v(2),:); end
                nrm_c=cross(d1,d2); if nrm_c(3) > 0, n_up=n_up+1; end
            end
        end



        function s = wp_arc_len(curve)
        % WP_ARC_LEN  Cumulative arc length along a 2D polyline.
            s = zeros(size(curve, 1), 1);
            for k = 2:size(curve, 1)
                s(k) = s(k-1) + norm(curve(k,:) - curve(k-1,:));
            end
        end


        function pt = wp_interp(curve, s_cum, t)
        % WP_INTERP  Interpolate position on arc-length parameterised curve.
            st = t * s_cum(end);
            if st <= 0;        pt = curve(1,:);   return; end
            if st >= s_cum(end); pt = curve(end,:); return; end
            idx = find(s_cum >= st, 1, 'first');
            if idx <= 1;       pt = curve(1,:);   return; end
            f = (st - s_cum(idx-1)) / (s_cum(idx) - s_cum(idx-1));
            pt = curve(idx-1,:) + f * (curve(idx,:) - curve(idx-1,:));
        end


        %% ═══════════════════════════════════════════════════════════
        %%  PANELIZER MESH ADAPTER (WEC_Panelizer → HAMS format)
        %% ═══════════════════════════════════════════════════════════

        function write_panelizer_hull_pnl(mesh, filepath)
            % WRITE_PANELIZER_HULL_PNL  Write a WEC_Panelizer mesh in HAMS format.
            %
            %   HAMS_Pipeline.write_panelizer_hull_pnl(mesh, filepath)
            %
            %   Takes a mesh struct from WEC_Panelizer.generate() and writes
            %   it using HAMS_Pipeline.write_pnl_file (verified HAMS format).
            %
            %   Handles the format differences:
            %     - WEC_Panelizer uses v4 = v3 for triangles
            %     - HAMS needs the nverts column (3 or 4)
            %     - HAMS header includes symmetry flags and section markers
            %     - WEC_Panelizer expands all mirrors → x_sym=0, y_sym=0
            %
            %   INPUTS
            %     mesh     — struct from WEC_Panelizer.generate()
            %     filepath — output .pnl file path

            nodes  = mesh.vertices;
            panels = mesh.panels;

            % Detect triangles (WEC_Panelizer convention: v4 == v3)
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

            HAMS_Pipeline.write_pnl_file(filepath, nodes, panels, ...
                panel_nverts, x_sym, y_sym);
        end


        function write_panelizer_wp_pnl(mesh, filepath, z_wl, target_edge)
            % WRITE_PANELIZER_WP_PNL  Structured waterplane lid from hull mesh.
            %
            %   Extracts the waterline boundary from hull open edges at z≈0,
            %   merges duplicate vertices (from split_panel_at_z), walks the
            %   edge chain for correct ordering, then calls the structured
            %   quad grid builder.
            %
            %   QUARTER-BODY MODE (mesh.x_sym=1, mesh.y_sym=1):
            %     The hull mesh covers Q1 only.  The waterline boundary is
            %     an open arc from the X-axis to the Y-axis.  This method
            %     closes the boundary by adding segments along the two
            %     symmetry planes (Y=0 and X=0) through the origin.  The
            %     resulting closed polygon is passed to the TFI builder.
            %     The .pnl file is written with [1,1] symmetry flags so
            %     HAMS mirrors the quarter WP mesh to the full waterplane.
            %
            %   INPUTS
            %     mesh         — struct from WEC_Panelizer.generate()
            %     filepath     — output .pnl file path
            %     z_wl         — waterline elevation [m] (default: 0)
            %     target_edge  — max panel edge length [m] (default: 0.4)

            if nargin < 3 || isempty(z_wl);         z_wl = 0; end
            if nargin < 4 || isempty(target_edge); target_edge = 0.4; end

            is_quarter = isfield(mesh, 'x_sym') && mesh.x_sym == 1 ...
                      && isfield(mesh, 'y_sym') && mesh.y_sym == 1;

            % ── Step 1: Find hull open edges at z ≈ 0 ────────────────
            z_tol = 0.01;
            edge_cnt = containers.Map('KeyType','char','ValueType','int32');
            edge_map = containers.Map('KeyType','char','ValueType','any');

            n_p = size(mesh.panels, 1);
            for p = 1:n_p
                v = mesh.panels(p,:);
                if v(3) == v(4)
                    ee = [v(1) v(2); v(2) v(3); v(3) v(1)];
                else
                    ee = [v(1) v(2); v(2) v(3); v(3) v(4); v(4) v(1)];
                end
                for e = 1:size(ee,1)
                    key = sprintf('%d_%d', min(ee(e,:)), max(ee(e,:)));
                    if edge_cnt.isKey(key)
                        edge_cnt(key) = edge_cnt(key) + 1;
                    else
                        edge_cnt(key) = 1;
                        edge_map(key) = ee(e,:);
                    end
                end
            end

            wl_edges = zeros(0, 2);
            keys_all = edge_cnt.keys();
            for i = 1:length(keys_all)
                if edge_cnt(keys_all{i}) == 1
                    ev = edge_map(keys_all{i});
                    if abs(mesh.vertices(ev(1),3)) < z_tol && ...
                       abs(mesh.vertices(ev(2),3)) < z_tol
                        wl_edges(end+1,:) = ev; %#ok<AGROW>
                    end
                end
            end

            if isempty(wl_edges)
                warning('HAMS_Pipeline:NoWLEdges', 'No waterline edges found.');
                HAMS_Pipeline.write_pnl_file(filepath, zeros(0,3), ...
                    zeros(0,4), zeros(0,1), 0, 0);
                return;
            end

            % ── Step 2: Extract vertices and merge duplicates ─────────
            wl_vert_idx = unique(wl_edges(:));
            wl_xyz = mesh.vertices(wl_vert_idx, :);
            n_wl = size(wl_xyz, 1);

            g2l = containers.Map('KeyType','int32','ValueType','int32');
            for i = 1:n_wl; g2l(wl_vert_idx(i)) = i; end

            local_edges = zeros(size(wl_edges));
            for e = 1:size(wl_edges,1)
                local_edges(e,1) = g2l(wl_edges(e,1));
                local_edges(e,2) = g2l(wl_edges(e,2));
            end

            merge_tol = 1e-6;
            canon = (1:n_wl)';
            for ii = 2:n_wl
                for jj = 1:ii-1
                    if canon(jj) ~= jj; continue; end
                    if norm(wl_xyz(ii,:) - wl_xyz(jj,:)) < merge_tol
                        canon(ii) = jj; break;
                    end
                end
            end

            for e = 1:size(local_edges,1)
                local_edges(e,1) = canon(local_edges(e,1));
                local_edges(e,2) = canon(local_edges(e,2));
            end
            local_edges(local_edges(:,1)==local_edges(:,2),:) = [];
            [~, ui] = unique(sort(local_edges,2), 'rows');
            local_edges = local_edges(ui,:);

            unique_ids = unique(canon);
            n_unique = length(unique_ids);
            new_idx = zeros(n_wl, 1);
            new_idx(unique_ids) = (1:n_unique)';

            compact_edges = zeros(size(local_edges));
            for e = 1:size(local_edges,1)
                compact_edges(e,1) = new_idx(canon(local_edges(e,1)));
                compact_edges(e,2) = new_idx(canon(local_edges(e,2)));
            end
            unique_xy = wl_xyz(unique_ids, 1:2);

            n_merged = n_wl - n_unique;
            fprintf('    WP: %d wl verts, %d merged → %d unique, %d edges\n', ...
                    n_wl, n_merged, n_unique, size(compact_edges,1));

            if is_quarter
                % ─── QUARTER BODY: structured mesh from hull WL arc ────
                %
                %  The hull waterline vertices in Q1 are exact spline
                %  evaluations (from wl_conform).  They form an arc from
                %  the Y-axis to the X-axis.
                %
                %  MESH CONSTRUCTION
                %    For each arc vertex (x_j, y_j), create a cross-wise
                %    row from the X=0 symmetry plane (0, y_j) to the hull
                %    boundary (x_j, y_j).  This gives a structured quad
                %    grid [N_cross+1 × N_arc] with:
                %      - Right edge = hull waterline (exact node match)
                %      - Left edge  = X=0 symmetry plane
                %      - Top edge   = near Y-axis (j=1)
                %      - Bottom edge = near X-axis (j=N_arc, ≈ Y=0)
                %    All quads, no fan triangles, no collapsed vertices.

                % Sort arc vertices: Y-axis (high θ) → X-axis (low θ)
                %   This ordering gives cross(P_i, P_j) in -Z (DOWN),
                %   which is the correct HAMS waterplane normal direction.
                angles = atan2(unique_xy(:,2), unique_xy(:,1));
                [~, order] = sort(angles, 'descend');
                arc_sorted = unique_xy(order, :);
                N_arc = size(arc_sorted, 1);

                % Determine N_cross from aspect ratio matching
                arc_len = 0;
                for k = 2:N_arc
                    arc_len = arc_len + norm(arc_sorted(k,:) - arc_sorted(k-1,:));
                end
                mean_ds  = arc_len / max(N_arc - 1, 1);
                max_width = max(arc_sorted(:, 1));
                N_cross = max(2, ceil(max_width / mean_ds));

                % Build structured grid [Ni × Nj × 3]
                Ni = N_cross + 1;
                Nj = N_arc;
                wp_grid = zeros(Ni, Nj, 3);
                for j = 1:Nj
                    xb = arc_sorted(j, 1);
                    yb = arc_sorted(j, 2);
                    for i = 1:Ni
                        frac = (i - 1) / N_cross;
                        wp_grid(i, j, 1) = frac * xb;
                        wp_grid(i, j, 2) = yb;
                        wp_grid(i, j, 3) = z_wl;
                    end
                end

                % Grid to quads (column-major indexing)
                wp_nodes = reshape(wp_grid, [], 3);
                n_wp = (Ni - 1) * (Nj - 1);
                wp_panels = zeros(n_wp, 4);
                pp = 0;
                for i = 1:Ni-1
                    for j = 1:Nj-1
                        pp = pp + 1;
                        wp_panels(pp, :) = [ ...
                            i     + (j-1)*Ni, ...
                            (i+1) + (j-1)*Ni, ...
                            (i+1) + j*Ni,     ...
                            i     + j*Ni];
                    end
                end
                wp_nverts = 4 * ones(n_wp, 1);

                % Verify normals point UPWARD (+z)
                mid_p = ceil(n_wp / 2);
                vi = wp_panels(mid_p, :);
                d1 = wp_nodes(vi(3),:) - wp_nodes(vi(1),:);
                d2 = wp_nodes(vi(4),:) - wp_nodes(vi(2),:);
                nz = cross(d1, d2);
                if nz(3) < 0
                    wp_panels = wp_panels(:, [1 4 3 2]);
                    fprintf('    WP normals: reversed winding to point UPWARD\n');
                end

                % Report element quality
                min_area = Inf;
                for pp_chk = 1:n_wp
                    vi_c = wp_panels(pp_chk, :);
                    d1_c = wp_nodes(vi_c(3),:) - wp_nodes(vi_c(1),:);
                    d2_c = wp_nodes(vi_c(4),:) - wp_nodes(vi_c(2),:);
                    a_c  = norm(cross(d1_c, d2_c)) / 2;
                    if a_c < min_area, min_area = a_c; end
                end

                fprintf('    WP quarter structured: %d×%d grid → %d quads\n', ...
                        Ni, Nj, n_wp);
                fprintf('    WP arc pts: %d, N_cross: %d, min panel area: %.2e m²\n', ...
                        N_arc, N_cross, min_area);

                HAMS_Pipeline.write_pnl_file(filepath, wp_nodes, wp_panels, ...
                    wp_nverts, 1, 1);
            else
                % ─── FULL BODY: existing edge-chain approach ──────────────

                % Walk edge chain → ordered boundary
                adj = cell(n_unique, 1);
                for e = 1:size(compact_edges,1)
                    v1 = compact_edges(e,1); v2 = compact_edges(e,2);
                    adj{v1}(end+1) = v2;
                    adj{v2}(end+1) = v1;
                end

                visited = false(n_unique, 1);
                order = zeros(n_unique, 1);
                order(1) = 1; visited(1) = true;
                for step = 2:n_unique
                    curr = order(step-1);
                    nbrs = adj{curr};
                    next = 0;
                    for ni = 1:length(nbrs)
                        if ~visited(nbrs(ni)); next = nbrs(ni); break; end
                    end
                    if next == 0; break; end
                    order(step) = next; visited(next) = true;
                end
                n_chain = find(order > 0, 1, 'last');
                order = order(1:n_chain);
                boundary_xy = unique_xy(order, :);

                fprintf('    WP chain: %d vertices in ordered loop\n', n_chain);

                % Build structured quad grid
                [wp_nodes, wp_panels, wp_nverts] = ...
                    HAMS_Pipeline.generate_waterplane_mesh_structured( ...
                        boundary_xy, target_edge, z_wl);

                HAMS_Pipeline.write_pnl_file(filepath, wp_nodes, wp_panels, ...
                    wp_nverts, 0, 0);
            end
        end


        %% ═══════════════════════════════════════════════════════════
        %%  HAMS INPUT BUILDER (from WEC geometry)
        %% ═══════════════════════════════════════════════════════════

        function [CG_global, M_6x6, C_6x6, mass, sub_props] = ...
                compute_hams_inputs(parser, z_wl, config)
            % COMPUTE_HAMS_INPUTS  Build Hydrostatic.in inputs from geometry.
            %
            %   [CG, M, C, m, sub] = HAMS_Pipeline.compute_hams_inputs(parser, z_wl, config)
            %
            %   Computes CG position, 6x6 mass matrix (about origin), and
            %   6x6 restoring matrix for a UNIFORM-DENSITY floating body at
            %   a given waterline elevation.
            %
            %   INPUTS
            %     parser  — WEC_MS2_Parser object
            %     z_wl    — waterline elevation in body frame [m]
            %     config  — struct with: .RHO_WATER, .G, .total_wec_volume,
            %               .hull_centroid, .hull_int_x2/y2/z2, .Aw_table_z/Aw_table,
            %               .I_wp_xx_table, .I_wp_yy_table
            %
            %   OUTPUTS
            %     CG_global  [1x3]  CG in global frame [m]
            %     M_6x6      [6x6]  mass matrix about origin
            %     C_6x6      [6x6]  restoring matrix about origin
            %     mass       scalar body mass [kg]
            %     sub_props  struct from compute_submerged (V_sub, CB, Aw, I_wp)

            rho_w = config.RHO_WATER;
            g     = config.G;
            V_total  = config.total_wec_volume;
            centroid = config.hull_centroid;
            draft = -z_wl;

            % 1. Submerged properties from precomputed tables (v9.0)
            %
            %  Uses the Aw-table trapz path validated in the config builder.
            %  Bypasses compute_submerged entirely — avoids the divergence-
            %  theorem partial-surface errors for offset-axis geometries.
            %
            %  Fallback: if V_sub_table is absent (geometry-only config or
            %  legacy struct), calls compute_submerged as before.
            if isfield(config, 'V_sub_table') && ~isempty(config.V_sub_table)
                z_wl_clamped = max(min(z_wl, config.Aw_table_z(end)), ...
                                       config.Aw_table_z(1));
                sub_props = struct();
                sub_props.V_sub   = max(0, interp1(config.Aw_table_z, ...
                                    config.V_sub_table,   z_wl_clamped, 'linear'));
                sub_props.CB      = [0, 0, interp1(config.Aw_table_z, ...
                                    config.CB_z_table,    z_wl_clamped, 'linear')];
                sub_props.Aw      = max(0, interp1(config.Aw_table_z, ...
                                    config.Aw_table,      z_wl_clamped, 'linear'));
                sub_props.I_wp_xx = max(0, interp1(config.Aw_table_z, ...
                                    config.I_wp_xx_table, z_wl_clamped, 'linear'));
                sub_props.I_wp_yy = max(0, interp1(config.Aw_table_z, ...
                                    config.I_wp_yy_table, z_wl_clamped, 'linear'));
            else
                % Legacy fallback: compute_submerged (may have errors for
                % offset-axis geometries — retained for backward compatibility)
                sub_opts = struct('n_quad', 16, 'verbose', false);
                if isfield(config, 'Aw_table_z') && ~isempty(config.Aw_table_z)
                    sub_opts.Aw_override = interp1( ...
                        config.Aw_table_z, config.Aw_table, z_wl, 'linear', 0);
                    sub_opts.I_wp_xx_override = interp1( ...
                        config.Aw_table_z, config.I_wp_xx_table, z_wl, 'linear', 0);
                    sub_opts.I_wp_yy_override = interp1( ...
                        config.Aw_table_z, config.I_wp_yy_table, z_wl, 'linear', 0);
                end
                sub_props = WEC_HydroProperties.compute_submerged(parser, z_wl, sub_opts);
            end

            % 2. Mass (uniform density, Archimedes)
            mass = rho_w * sub_props.V_sub;
            rho_body = mass / V_total;

            % 3. CG in global frame
            %   Bilateral symmetry: CG must lie on the z-axis.
            %   The parametric centroid has x,y ≈ 0 (O(1e-17) from
            %   numerical noise in the divergence theorem). Force exact.
            CG_global = [0, 0, centroid(3) + draft];

            % 4. Inertia tensor about CG (frame-invariant)
            %   Use cx = cy = 0 (bilateral symmetry) to avoid noise.
            cx = 0; cy = 0; cz = centroid(3);
            Ixx_cg = max(0, rho_body * (config.hull_int_y2 + config.hull_int_z2 ...
                     - V_total * (cy^2 + cz^2)));
            Iyy_cg = max(0, rho_body * (config.hull_int_x2 + config.hull_int_z2 ...
                     - V_total * (cx^2 + cz^2)));
            Izz_cg = max(0, rho_body * (config.hull_int_x2 + config.hull_int_y2 ...
                     - V_total * (cx^2 + cy^2)));
            J_CG = diag([Ixx_cg, Iyy_cg, Izz_cg]);

            % 5. 6x6 mass matrix about origin
            M_6x6 = HAMS_Pipeline.build_mass_matrix_6x6(mass, CG_global, J_CG);

            % 6. Restoring matrix about origin
            z_B = sub_props.CB(3) + draft;
            z_G = CG_global(3);
            C_6x6 = zeros(6);
            C_6x6(3,3) = rho_w * g * sub_props.Aw;
            C_6x6(4,4) = rho_w*g*sub_props.I_wp_xx + rho_w*g*sub_props.V_sub*z_B - mass*g*z_G;
            C_6x6(5,5) = rho_w*g*sub_props.I_wp_yy + rho_w*g*sub_props.V_sub*z_B - mass*g*z_G;
        end


        %% ═══════════════════════════════════════════════════════════
        %%  PIPELINE DRIVERS
        %% ═══════════════════════════════════════════════════════════
        %%  LEGACY ELLIPSOID PIPELINE
        %%  These functions (run_single_draft, run_5draft_sweep,
        %%  interpolate_Ainf) use generate_axisymmetric_mesh and the
        %%  simple circular waterplane mesher. They are superseded by
        %%  run_single_hams + WEC_Panelizer and are retained only for
        %%  reference. Do not use in new development.
        %% ═══════════════════════════════════════════════════════════

        function hams_data = run_single_draft(a, c, draft, ...
                hams_exe, base_dir, run_name, hams_params)
            % RUN_SINGLE_DRAFT  Full HAMS pipeline for one draft.
            %
            %   hams_data = run_single_draft(a, c, draft, exe, dir, name, params)
            %
            %   Steps:
            %     1. Generate ellipsoid profile at given draft
            %     2. Generate BEM mesh (hull + waterplane lid)
            %     3. Compute geometric hydrostatics
            %     4. Write all input files
            %     5. Run HAMS
            %     6. Parse output (with denormalization)
            %
            %   Returns hams_data struct with dimensional A(∞).

            rho_w = HAMS_Pipeline.RHO_WATER;
            g = HAMS_Pipeline.G;

            run_dir = fullfile(base_dir, run_name);
            fprintf('\n=== HAMS run: %s (draft = %.3f m) ===\n', run_name, draft);

            % 1. Profile
            N_z = 16;
            [r_prof, z_prof] = HAMS_Pipeline.generate_ellipsoid_profile(a, c, draft, N_z);

            % 2a. Hull mesh
            N_theta_half = 24;
            use_x_sym = false;
            [nodes, panels, nverts] = HAMS_Pipeline.generate_axisymmetric_mesh(...
                r_prof, z_prof, N_theta_half, use_x_sym);

            % 2b. Waterplane lid mesh (for irregular frequency removal)
            r_wp = r_prof(end);  % waterplane radius (last point is at z=0)
            N_radial_wp = 8;
            [wp_nodes, wp_panels, wp_nverts] = HAMS_Pipeline.generate_waterplane_mesh(...
                r_wp, N_theta_half, N_radial_wp, use_x_sym);

            % 3. Hydrostatics
            hydro = HAMS_Pipeline.compute_geometric_hydrostatics(...
                r_prof, z_prof, rho_w, g);

            % 4. Setup directory and write files
            HAMS_Pipeline.setup_hams_directory(run_dir);

            mesh_file = fullfile(run_dir, 'Input', 'HullMesh.pnl');
            HAMS_Pipeline.write_pnl_file(mesh_file, nodes, panels, nverts, ...
                0, 1);  % x_sym=0, y_sym=1

            wp_file = fullfile(run_dir, 'Input', 'WaterPlaneMesh.pnl');
            HAMS_Pipeline.write_pnl_file(wp_file, wp_nodes, wp_panels, wp_nverts, ...
                0, 1);  % same symmetry as hull

            ctrl_file = fullfile(run_dir, 'Input', 'ControlFile.in');
            ctrl_params = hams_params;
            % XR = [0,0,0]: HAMS outputs A, B, Fe at the global origin.
            % WEC_Configuration_Builder §5 / rebuild_config_hydro apply the
            % single congruence transform origin → CG in post-processing.
            ctrl_params.ref_body_center = [0, 0, 0];
            HAMS_Pipeline.write_control_file(ctrl_file, ctrl_params);

            hydro_file = fullfile(run_dir, 'Input', 'Hydrostatic.in');
            HAMS_Pipeline.write_hydrostatic_file(hydro_file, ...
                hydro.CG, hydro.M_6x6, ...
                zeros(6), zeros(6), ...
                hydro.C_6x6, zeros(6));

            % 5. Run HAMS
            [status, ~] = HAMS_Pipeline.run_hams(hams_exe, run_dir);

            % 6. Parse output (with WAMIT non-dimensional → SI conversion)
            ref_L = 1.0;  % reference_body_length from ControlFile.in
            if isfield(ctrl_params, 'ref_body_length')
                ref_L = ctrl_params.ref_body_length;
            end

            if status == 0
                wamit_dir = fullfile(run_dir, 'Output', 'Wamit_format');
                one_files = dir(fullfile(wamit_dir, '*.1'));
                if ~isempty(one_files)
                    hams_data = HAMS_Pipeline.parse_wamit_1_file(...
                        fullfile(wamit_dir, one_files(1).name), ref_L);
                else
                    warning('No .1 file found in %s', wamit_dir);
                    hams_data = struct('A_inf', zeros(6), 'status', 'no_output');
                end
            else
                hams_data = struct('A_inf', zeros(6), 'status', 'hams_failed');
            end

            % Attach hydrostatic data for reference
            hams_data.hydro = hydro;
            hams_data.draft = draft;
        end


        function sweep = run_5draft_sweep(a, c, draft_min, draft_max, ...
                hams_exe, base_dir, hams_params)
            % RUN_5DRAFT_SWEEP  Five-draft HAMS sweep for A(∞) table.
            %
            %   sweep = run_5draft_sweep(a, c, d_min, d_max, exe, dir, params)
            %
            %   Runs HAMS at 5 drafts: [d_min, d_min/2, 0, d_max/2, d_max]
            %   Stores A(∞) at each draft for linear interpolation.
            %
            %   OUTPUT struct:
            %     .drafts      [5×1]    draft values
            %     .A_inf       {5×1}    cell array of 6×6 A(∞) matrices
            %     .A33_inf     [5×1]    A_inf(3,3) at each draft (heave)
            %     .A55_inf     [5×1]    A_inf(5,5) at each draft (pitch)
            %     .hydro       {5×1}    cell array of hydrostatic structs

            drafts = [draft_min, draft_min/2, 0, draft_max/2, draft_max];
            n_drafts = length(drafts);

            sweep.drafts = drafts(:);
            sweep.A_inf = cell(n_drafts, 1);
            sweep.A33_inf = zeros(n_drafts, 1);
            sweep.A55_inf = zeros(n_drafts, 1);
            sweep.hydro = cell(n_drafts, 1);

            for k = 1:n_drafts
                run_name = sprintf('draft_%+.3f', drafts(k));
                hams_data = HAMS_Pipeline.run_single_draft(...
                    a, c, drafts(k), hams_exe, base_dir, run_name, hams_params);

                sweep.A_inf{k} = hams_data.A_inf;
                sweep.A33_inf(k) = hams_data.A_inf(3,3);
                sweep.A55_inf(k) = hams_data.A_inf(5,5);
                sweep.hydro{k} = hams_data.hydro;
            end

            fprintf('\n=== 5-Draft Sweep Complete ===\n');
            fprintf('  Draft [m]    A₃₃(∞) [kg]    A₅₅(∞) [kg·m²]    C₃₃ [N/m]     T_heave [s]    T_pitch [s]\n');
            g = HAMS_Pipeline.G;
            rho = HAMS_Pipeline.RHO_WATER;
            for k = 1:n_drafts
                C33_k = sweep.hydro{k}.C33;
                C55_k = sweep.hydro{k}.C55;
                m_k   = sweep.hydro{k}.mass;
                Iyy_k = sweep.hydro{k}.Iyy;
                T3 = 2*pi*sqrt((m_k + sweep.A33_inf(k)) / C33_k);
                T5 = 2*pi*sqrt((Iyy_k + sweep.A55_inf(k)) / C55_k);
                fprintf('  %+.3f       %8.1f        %8.1f         %8.1f       %6.3f         %6.3f\n', ...
                    drafts(k), sweep.A33_inf(k), sweep.A55_inf(k), C33_k, T3, T5);
            end
        end


        function A_inf_interp = interpolate_Ainf(draft, sweep)
            % INTERPOLATE_AINF  Linear interpolation of A(∞) from draft table.
            %
            %   A_inf = interpolate_Ainf(draft, sweep)
            %
            %   Interpolates each element of the 6×6 A(∞) matrix linearly
            %   between the 5 pre-computed drafts.

            A_inf_interp = zeros(6, 6);
            for i = 1:6
                for j = 1:6
                    vals = cellfun(@(A) A(i,j), sweep.A_inf);
                    A_inf_interp(i,j) = interp1(sweep.drafts, vals, draft, ...
                        'linear', 'extrap');
                end
            end
        end


        function hydro_table = run_draft_sweep(config, hams_exe, z_levels, ...
                output_dir, options)
            % RUN_DRAFT_SWEEP  Run HAMS at multiple drafts → hydro_table.
            %
            %   ht = HAMS_Pipeline.run_draft_sweep(config, exe, z_levels, dir)
            %   ht = HAMS_Pipeline.run_draft_sweep(config, exe, z_levels, dir, opts)
            %
            %   For each z_level (body-frame waterline):
            %     1. Generate hull mesh (WEC_Panelizer)
            %     2. Write hull + WP mesh (HAMS adapters)
            %     3. Compute mass + restoring (uniform density)
            %     4. Write ControlFile.in + Hydrostatic.in
            %     5. Run HAMS
            %     6. Parse output → A(inf), A(w), B(w)
            %     7. Average B over [T_min, T_max]
            %
            %   INPUTS
            %     config     — from WEC_Configuration_Builder
            %     hams_exe   — full path to HAMS executable
            %     z_levels   — [Nx1] waterline z in body frame [m]
            %     output_dir — base directory for run folders
            %     options    — (optional) struct:
            %       .Nu, .Nv         — panelizer grid (default: config values)
            %       .wp_target_edge  — WP mesh panel edge [m] (default: 0.4)
            %       .T_min, .T_max   — B averaging band [s] (default: [4, 16])
            %       .verbose         — print progress (default: true)
            %       .skip_existing   — skip completed drafts (default: false)
            %
            %   OUTPUT
            %     hydro_table — struct with:
            %       .z_levels, .drafts, .z_cg, .omega
            %       .A_inf {Nx1 of 6x6}, .B_avg {Nx1 of 6x6}
            %       .A {Nx1 of 6x6xM}, .B {Nx1 of 6x6xM}
            %       .V_sub [Nx1], .mass [Nx1]

            if nargin < 5; options = struct(); end
            if ~isfield(options, 'Nu');              options.Nu = config.mesh_Nu; end
            if ~isfield(options, 'Nv');              options.Nv = config.mesh_Nv; end
            if ~isfield(options, 'wp_target_edge');  options.wp_target_edge = 0.4; end
            if ~isfield(options, 'T_min');           options.T_min = 4.0; end
            if ~isfield(options, 'T_max');           options.T_max = 16.0; end
            if ~isfield(options, 'verbose');         options.verbose = true; end
            if ~isfield(options, 'skip_existing');   options.skip_existing = false; end

            parser = config.ms2_model;
            N = length(z_levels);
            z_levels = z_levels(:);
            hams_params = HAMS_Pipeline.default_hams_params(config);

            if options.verbose
                fprintf('\n  HAMS Draft Sweep: %d waterlines\n', N);
                fprintf('    Output: %s\n', output_dir);
                fprintf('    Grid: %dx%d, B band: [%.1f, %.1f] s\n', ...
                        options.Nu, options.Nv, options.T_min, options.T_max);
            end

            hydro_table.z_levels    = z_levels;
            hydro_table.drafts      = -z_levels;
            hydro_table.z_cg        = zeros(N, 1);
            hydro_table.V_sub       = zeros(N, 1);
            hydro_table.mass        = zeros(N, 1);
            hydro_table.A_inf       = cell(N, 1);
            hydro_table.B_avg       = cell(N, 1);
            hydro_table.A           = cell(N, 1);
            hydro_table.B           = cell(N, 1);
            hydro_table.omega       = [];
            hydro_table.hams_params = hams_params;
            hydro_table.T_band      = [options.T_min, options.T_max];

            if ~exist(output_dir, 'dir'); mkdir(output_dir); end

            t_sweep = tic;

            for k = 1:N
                z_wl  = z_levels(k);
                draft = -z_wl;
                run_name = sprintf('draft_%+.4f', draft);
                run_dir  = fullfile(output_dir, run_name);

                if options.verbose
                    fprintf('\n  [%d/%d] z_wl=%.4f m (vertical_shift=%+.4f)\n', ...
                            k, N, z_wl, draft);
                end

                % Skip if output already exists
                one_pat = fullfile(run_dir, 'Output', 'Wamit_format', '*.1');
                if options.skip_existing && ~isempty(dir(one_pat))
                    if options.verbose; fprintf('    Skipping (exists)\n'); end
                    one_f = dir(one_pat);
                    hd = HAMS_Pipeline.parse_wamit_1_file( ...
                        fullfile(run_dir,'Output','Wamit_format',one_f(1).name), 1.0, 3);
                    [CG_g, ~, ~, m_k, sub_k] = ...
                        HAMS_Pipeline.compute_hams_inputs(parser, z_wl, config);
                    hydro_table.z_cg(k) = CG_g(3);
                    hydro_table.V_sub(k) = sub_k.V_sub;
                    hydro_table.mass(k) = m_k;
                    hydro_table.A_inf{k} = hd.A_inf;
                    hydro_table.A{k} = hd.A;
                    hydro_table.B{k} = hd.B;
                    if isempty(hydro_table.omega); hydro_table.omega = hd.omega; end
                    hydro_table.B_avg{k} = HAMS_Pipeline.compute_B_avg( ...
                        hd, options.T_min, options.T_max);
                    continue;
                end

                % 1-2. Mesh (full body, ISX=0 ISY=0 — HAMS symmetry is buggy)
                pan_opts = struct('trim_wl', true, 'close_gaps', false, ...
                                  'cosine_spacing', false, 'verbose', false, ...
                                  'quarter_body', false, 'half_body', false);
                mesh = WEC_Panelizer.generate(parser, draft, ...
                           options.Nu, options.Nv, pan_opts);

                HAMS_Pipeline.setup_hams_directory(run_dir);
                HAMS_Pipeline.write_panelizer_hull_pnl(mesh, ...
                    fullfile(run_dir, 'Input', 'HullMesh.pnl'));

                % Waterplane lid — boundary inherited from hull mesh (WP-MATCH FIX)
                % Replaces extract_isocurve_at_z path (WP-3).  See run_single_hams
                % comment for full explanation.  WP-2 FIX in WEC_Panelizer ensures
                % seam vertices are merged, so the open-edge walk is topologically
                % connected.  lock_boundary=true preserves hull node positions on
                % the WP outer boundary; target_edge governs interior density only.
                boundary_xy_ds = HAMS_Pipeline.hull_waterline_polygon(mesh);

                if ~isempty(boundary_xy_ds) && size(boundary_xy_ds, 1) >= 3
                    [wpn, wpp, wpv] = HAMS_Pipeline.mesh_wp_blossomquad( ...
                        boundary_xy_ds, options.wp_target_edge, mesh.x_sym, mesh.y_sym, 0, true);
                    if ~isempty(wpn), wpn(:, 3) = 0; end
                else
                    warning('HAMS_Pipeline:WPFailed', ...
                        'hull_waterline_polygon empty at z_wl=%.4f m.', z_wl);
                    wpn = zeros(0,3); wpp = zeros(0,4); wpv = zeros(0,1);
                end

                HAMS_Pipeline.write_pnl_file( ...
                    fullfile(run_dir, 'Input', 'WaterPlaneMesh.pnl'), ...
                    wpn, wpp, wpv, mesh.x_sym, mesh.y_sym);

                % 3. Mass + restoring
                [CG_global, M_6x6, C_6x6, m_k, sub_k] = ...
                    HAMS_Pipeline.compute_hams_inputs(parser, z_wl, config);
                hydro_table.z_cg(k) = CG_global(3);
                hydro_table.V_sub(k) = sub_k.V_sub;
                hydro_table.mass(k) = m_k;

                if options.verbose
                    fprintf('    V_sub=%.4f m3, mass=%.1f kg, CG_z=%.4f m\n', ...
                            sub_k.V_sub, m_k, CG_global(3));
                end

                % 4. Write HAMS files
                ctrl = hams_params;
                % XR = [0,0,0]: HAMS outputs A, B, Fe at the global origin.
                % WEC_Configuration_Builder §5 applies the single
                % origin → CG congruence transform in post-processing.
                ctrl.ref_body_center = [0, 0, 0];
                HAMS_Pipeline.write_control_file( ...
                    fullfile(run_dir, 'Input', 'ControlFile.in'), ctrl);
                HAMS_Pipeline.write_hydrostatic_file( ...
                    fullfile(run_dir, 'Input', 'Hydrostatic.in'), ...
                    CG_global, M_6x6, zeros(6), zeros(6), C_6x6, zeros(6));

                % 5. Run HAMS
                [status, result] = HAMS_Pipeline.run_hams(hams_exe, run_dir);
                if status ~= 0
                    warning('HAMS_Pipeline:SweepRunFailed', ...
                            'HAMS failed at draft %.4f: %s', draft, result);
                    hydro_table.A_inf{k} = zeros(6);
                    hydro_table.B_avg{k} = zeros(6);
                    continue;
                end

                % 6. Parse output
                one_f = dir(fullfile(run_dir,'Output','Wamit_format','*.1'));
                if isempty(one_f)
                    warning('HAMS_Pipeline:NoOutput', 'No .1 file at draft %.4f', draft);
                    hydro_table.A_inf{k} = zeros(6);
                    hydro_table.B_avg{k} = zeros(6);
                    continue;
                end
                hd = HAMS_Pipeline.parse_wamit_1_file( ...
                    fullfile(run_dir,'Output','Wamit_format',one_f(1).name), 1.0, 3);

                hydro_table.A_inf{k} = hd.A_inf;
                hydro_table.A{k} = hd.A;
                hydro_table.B{k} = hd.B;
                if isempty(hydro_table.omega); hydro_table.omega = hd.omega; end

                % 7. B average
                hydro_table.B_avg{k} = HAMS_Pipeline.compute_B_avg( ...
                    hd, options.T_min, options.T_max);

                if options.verbose
                    fprintf('    A33(inf)=%.1f kg, A55(inf)=%.1f kg*m2\n', ...
                            hd.A_inf(3,3), hd.A_inf(5,5));
                end
            end

            elapsed = toc(t_sweep);
            if options.verbose
                fprintf('\n  Draft sweep complete: %.1f s (%.1f s/draft)\n', ...
                        elapsed, elapsed/N);
            end
        end


        %% ═══════════════════════════════════════════════════════════
        %%  HYDRO CACHE — CACHE-AWARE SINGLE-DRAFT HAMS RUNS
        %% ═══════════════════════════════════════════════════════════

        function cache = empty_hydro_cache()
            % EMPTY_HYDRO_CACHE  Create an empty hydro cache struct.
            %
            %   cache = HAMS_Pipeline.empty_hydro_cache()

            cache.drafts      = [];
            cache.z_cg        = [];
            cache.V_sub       = [];
            cache.mass        = [];
            cache.omega       = [];
            cache.A_inf       = {};
            cache.B_avg       = {};
            cache.A           = {};
            cache.B           = {};
            cache.Fe          = {};
            cache.T_band      = [4.0, 16.0];
            cache.hams_params = HAMS_Pipeline.default_hams_params();
            cache.timestamp   = datestr(now);
            cache.ms2_file    = '';
            cache.ms2_date    = '';
        end


        function hams_data = run_single_hams(config, vertical_shift, ...
                hams_dir, hams_exe, options)
            % RUN_SINGLE_HAMS  Run HAMS at one vertical_shift, return parsed data.
            %
            %   hams_data = HAMS_Pipeline.run_single_hams(config, vs, dir, exe)
            %   hams_data = HAMS_Pipeline.run_single_hams(config, vs, dir, exe, opts)
            %
            %   Steps:
            %     1. Generate panelizer mesh at this draft
            %     2. Write HullMesh.pnl, WaterPlaneMesh.pnl → dir/Input/
            %     3. Compute M_6x6, C_6x6 via compute_hams_inputs
            %     4. Write ControlFile.in, Hydrostatic.in → dir/Input/
            %     5. cd to dir, run HAMS executable
            %     6. Parse Output/Wamit_format/Buoy.1
            %     7. Compute B_avg over [T_min, T_max]
            %
            %   INPUTS
            %     config         — from WEC_Configuration_Builder
            %     vertical_shift — scalar [m]
            %     hams_dir       — path to HAMS directory (contains Input/, Output/)
            %     hams_exe       — full path to HAMS executable
            %     options        — (optional) struct:
            %       .Nu, .Nv         — panelizer grid (default: config values)
            %       .wp_target_edge  — WP panel edge [m] (default: 0.4)
            %       .T_min, .T_max   — B averaging band [s] (default: [4, 16])
            %       .verbose         — print progress (default: true)
            %
            %   OUTPUT
            %     hams_data — struct:
            %       .A_inf  [6×6], .B_avg [6×6], .A [6×6×M], .B [6×6×M]
            %       .omega [M×1], .z_cg scalar, .V_sub scalar, .mass scalar
            %       .status  char ('ok' or 'failed')

            if nargin < 5; options = struct(); end
            if ~isfield(options, 'Nu');             options.Nu = config.mesh_Nu; end
            if ~isfield(options, 'Nv');             options.Nv = config.mesh_Nv; end
            if ~isfield(options, 'wp_target_edge'); options.wp_target_edge = 0.4; end
            if ~isfield(options, 'T_min');          options.T_min = 4.0; end
            if ~isfield(options, 'T_max');          options.T_max = 16.0; end
            if ~isfield(options, 'verbose');        options.verbose = true; end
            % HAMS-MREL symmetry (ISX=1) is unvalidated and produces ~32%
            % error on added mass (hemisphere unit test, 2026-04-01).
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
            mesh = WEC_Panelizer.generate(parser, draft, ...
                       options.Nu, options.Nv, pan_opts);

            % 2. Write mesh files
            HAMS_Pipeline.setup_hams_directory(hams_dir);
            HAMS_Pipeline.write_panelizer_hull_pnl(mesh, ...
                fullfile(hams_dir, 'Input', 'HullMesh.pnl'));

            % 2b. Waterplane lid — boundary inherited from hull mesh (WP-MATCH FIX)
            %
            % HISTORY
            %   WP-3: replaced mesh-edge walk with extract_isocurve_at_z because
            %     split_panel_at_z created coincident-but-index-distinct seam
            %     vertices, breaking the topological walk.
            %   WP-MATCH: WP-2 FIX (post-trim merge in WEC_Panelizer, 1e-8 tol)
            %     already fuses those seam duplicates.  The isocurve approach
            %     sampled the parametric surface at 200 uniform parameter values
            %     and resampled to target_edge/4 spacing — completely decoupled
            %     from the hull mesh node positions (arc-length / Nu-controlled).
            %     This produced WP boundary nodes at different xy positions to
            %     the hull waterline, degrading IRFR quality.
            %
            % FIX: hull_waterline_polygon walks open edges of mesh.panels at
            %   |z| < 0.01 m — the same mesh that was just written as
            %   HullMesh.pnl.  Passing lock_boundary=true to mesh_wp_blossomquad
            %   suppresses resample_polygon on the outer ring so hull nodes are
            %   preserved verbatim.  Interior Steiner density is still governed
            %   by options.wp_target_edge.  Symmetry flags are forwarded from
            %   mesh.x_sym / mesh.y_sym, unchanged from the previous path.
            boundary_xy = HAMS_Pipeline.hull_waterline_polygon(mesh);

            if ~isempty(boundary_xy) && size(boundary_xy, 1) >= 3
                [wp_nodes, wp_panels, wp_nverts] = HAMS_Pipeline.mesh_wp_blossomquad( ...
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
                warning('HAMS_Pipeline:WPFailed', ...
                    'hull_waterline_polygon empty at vs=%+.4f m. IRSP disabled.', vertical_shift);
                wp_nodes  = zeros(0,3);
                wp_panels = zeros(0,4);
                wp_nverts = zeros(0,1);
            end

            HAMS_Pipeline.write_pnl_file( ...
                fullfile(hams_dir, 'Input', 'WaterPlaneMesh.pnl'), ...
                wp_nodes, wp_panels, wp_nverts, mesh.x_sym, mesh.y_sym);

            % 3. Mass + restoring
            [CG_global, M_6x6, C_6x6, m_k, sub_k] = ...
                HAMS_Pipeline.compute_hams_inputs(parser, z_wl, config);

            if options.verbose
                fprintf('    V_sub=%.4f m3, mass=%.1f kg, CG_z=%.4f m\n', ...
                        sub_k.V_sub, m_k, CG_global(3));
            end

            % 4. Write HAMS input files
            hams_params = HAMS_Pipeline.default_hams_params(config);
            % XR = [0,0,0]: HAMS outputs A, B, Fe at the global origin.
            % WEC_Configuration_Builder §5 / rebuild_config_hydro apply the
            % single congruence transform origin → CG in post-processing.
            % retransform_at_actual_cg re-does origin → actual CG for the
            % converged draft in trained mode.
            hams_params.ref_body_center = [0, 0, 0];
            HAMS_Pipeline.write_control_file( ...
                fullfile(hams_dir, 'Input', 'ControlFile.in'), hams_params);
            HAMS_Pipeline.write_hydrostatic_file( ...
                fullfile(hams_dir, 'Input', 'Hydrostatic.in'), ...
                CG_global, M_6x6, zeros(6), zeros(6), C_6x6, zeros(6));

            % 5. Run HAMS
            [status, result] = HAMS_Pipeline.run_hams(hams_exe, hams_dir);

            if status ~= 0
                warning('HAMS_Pipeline:SingleRunFailed', ...
                        'HAMS failed at vs=%+.4f: %s', vertical_shift, result);
                hams_data.A_inf = zeros(6);
                hams_data.B_avg = zeros(6);
                hams_data.A     = [];
                hams_data.B     = [];
                hams_data.Fe    = [];
                hams_data.omega = [];
                hams_data.z_cg  = CG_global(3);
                hams_data.V_sub = sub_k.V_sub;
                hams_data.mass  = m_k;
                hams_data.status = 'failed';
                return;
            end

            % 6. Parse output
            one_files = dir(fullfile(hams_dir, 'Output', 'Wamit_format', '*.1'));
            if isempty(one_files)
                warning('HAMS_Pipeline:NoOutput', ...
                        'No .1 file at vs=%+.4f', vertical_shift);
                hams_data.A_inf = zeros(6);
                hams_data.B_avg = zeros(6);
                hams_data.A     = [];
                hams_data.B     = [];
                hams_data.Fe    = [];
                hams_data.omega = [];
                hams_data.z_cg  = CG_global(3);
                hams_data.V_sub = sub_k.V_sub;
                hams_data.mass  = m_k;
                hams_data.status = 'no_output';
                return;
            end

            hd = HAMS_Pipeline.parse_wamit_1_file( ...
                fullfile(hams_dir, 'Output', 'Wamit_format', one_files(1).name), ...
                1.0, 3);

            % 6b. Parse excitation force (.3 file)
            three_files = dir(fullfile(hams_dir, 'Output', 'Wamit_format', '*.3'));
            if ~isempty(three_files)
                fe_data = HAMS_Pipeline.parse_wamit_3_file( ...
                    fullfile(hams_dir, 'Output', 'Wamit_format', three_files(1).name), ...
                    1.0, 3);
                Fe_complex = fe_data.Fe;   % [6×M] complex, dimensional [N or N·m]
            else
                Fe_complex = [];
                warning('HAMS_Pipeline:No3File', ...
                        'No .3 file at vs=%+.4f — Fe unavailable', vertical_shift);
            end

            % 7. Assemble output
            hams_data.A_inf  = hd.A_inf;
            hams_data.A      = hd.A;
            hams_data.B      = hd.B;
            hams_data.omega  = hd.omega;
            hams_data.B_avg  = HAMS_Pipeline.compute_B_avg(hd, options.T_min, options.T_max);
            hams_data.Fe     = Fe_complex;
            hams_data.z_cg   = CG_global(3);
            hams_data.V_sub  = sub_k.V_sub;
            hams_data.mass   = m_k;
            hams_data.status = 'ok';

            % 7b. A(∞) NaN guard ─────────────────────────────────────────────
            %
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
            %  frequencies), A_inf stays zero — which screen_draft treats as
            %  a degenerate (infeasible) draft, as intended.
            if any(isnan(hams_data.A_inf(:)))
                n_hw  = min(5, size(hams_data.A, 3));
                if n_hw >= 1
                    warning('HAMS_Pipeline:AinfNaN', ...
                        ['A(inf) = NaN at vs=%+.4f m — HAMS high-freq BEM failed.\n' ...
                         '  Estimating A(inf) from top-%d frequency points ' ...
                         '(omega=[%.2f..%.2f] rad/s).'], ...
                        vertical_shift, n_hw, ...
                        hams_data.omega(max(1, end-n_hw+1)), hams_data.omega(end));
                    for ii = 1:6
                        for jj = 1:6
                            hw = squeeze(hams_data.A(ii, jj, end-n_hw+1:end));
                            if all(isfinite(hw))
                                hams_data.A_inf(ii, jj) = mean(hw);
                            else
                                hams_data.A_inf(ii, jj) = 0;
                            end
                        end
                    end
                    if options.verbose
                        fprintf('    A(inf) fallback applied: A33=%.1f kg, A55=%.1f kg*m2\n', ...
                                hams_data.A_inf(3,3), hams_data.A_inf(5,5));
                    end
                else
                    warning('HAMS_Pipeline:AinfNaN', ...
                        'A(inf) = NaN at vs=%+.4f m and no A(omega) data — setting A_inf=0.', ...
                        vertical_shift);
                    hams_data.A_inf = zeros(6);
                    hams_data.status = 'ainf_fallback_failed';
                end
            end

            if options.verbose
                fprintf('    A33(inf)=%.1f kg, A55(inf)=%.1f kg*m2\n', ...
                        hams_data.A_inf(3,3), hams_data.A_inf(5,5));
            end
        end


        function [hams_data, cache] = get_or_run_hams(vertical_shift, ...
                cache, config, hams_dir, hams_exe, tol, options)
            % GET_OR_RUN_HAMS  Cache-aware HAMS run at one vertical_shift.
            %
            %   [data, cache] = HAMS_Pipeline.get_or_run_hams(vs, cache, config, dir, exe, tol)
            %
            %   If vertical_shift is within tol of an existing cache entry,
            %   returns cached data (no HAMS run).  Otherwise runs HAMS,
            %   appends result to cache, returns both.
            %
            %   INPUTS
            %     vertical_shift — scalar [m]
            %     cache          — hydro_cache struct (from empty_hydro_cache or loaded)
            %     config         — from WEC_Configuration_Builder
            %     hams_dir       — path to HAMS directory
            %     hams_exe       — full path to HAMS executable
            %     tol            — tolerance for cache lookup [m] (default: 0.01)
            %     options        — (optional) passed to run_single_hams
            %
            %   OUTPUTS
            %     hams_data — struct with A_inf, B_avg, etc.
            %     cache     — updated cache (may have one new entry)

            if nargin < 6 || isempty(tol); tol = 0.01; end
            if nargin < 7; options = struct(); end

            % Cache lookup
            if ~isempty(cache.drafts)
                [min_dist, idx] = min(abs(cache.drafts - vertical_shift));
                if min_dist < tol
                    % Cache hit
                    hams_data.A_inf  = cache.A_inf{idx};
                    hams_data.B_avg  = cache.B_avg{idx};
                    hams_data.A      = cache.A{idx};
                    hams_data.B      = cache.B{idx};
                    hams_data.Fe     = cache.Fe{idx};
                    hams_data.omega  = cache.omega;
                    hams_data.z_cg   = cache.z_cg(idx);
                    hams_data.V_sub  = cache.V_sub(idx);
                    hams_data.mass   = cache.mass(idx);
                    hams_data.status = 'cached';
                    if isfield(options, 'verbose') && options.verbose
                        fprintf('  Cache hit: vs=%+.4f (matched vs=%+.4f, dist=%.4f m)\n', ...
                                vertical_shift, cache.drafts(idx), min_dist);
                    end
                    return;
                end
            end

            % Cache miss — run HAMS
            hams_data = HAMS_Pipeline.run_single_hams( ...
                config, vertical_shift, hams_dir, hams_exe, options);

            % Append to cache
            n = length(cache.drafts) + 1;
            cache.drafts(n,1)  = vertical_shift;
            cache.z_cg(n,1)    = hams_data.z_cg;
            cache.V_sub(n,1)   = hams_data.V_sub;
            cache.mass(n,1)    = hams_data.mass;
            cache.A_inf{n,1}   = hams_data.A_inf;
            cache.B_avg{n,1}   = hams_data.B_avg;
            cache.A{n,1}       = hams_data.A;
            cache.B{n,1}       = hams_data.B;
            cache.Fe{n,1}      = hams_data.Fe;
            if isempty(cache.omega) && ~isempty(hams_data.omega)
                cache.omega = hams_data.omega;
            end
            cache.timestamp = datestr(now);
        end


        function config = rebuild_config_hydro(config, cache)
            % REBUILD_CONFIG_HYDRO  Repopulate config.wamit_* from hydro cache.
            %
            %   config = HAMS_Pipeline.rebuild_config_hydro(config, cache)
            %
            % ── ANNOTATION CORRECTION (was incorrectly marked dead code) ──────
            % WEC_File_IO.transform_hydrodynamic_matrices_3x3 EXISTS at line 333
            % of WEC_File_IO.m and is callable.  This function is live and called
            % from WEC_Main_Optimizer (line 751) in 'trained' stage-1 mode.
            % ──────────────────────────────────────────────────────────────────
            %   Sorts the cache by draft, extracts [1,3,5] submatrix from
            %   each 6×6, transforms from origin to CG, and populates:
            %     config.wamit_drafts, config.wamit_A, config.wamit_A_full,
            %     config.wamit_B, config.wamit_B_full
            %
            %   This makes the downstream interpolation code work unchanged.

            N = length(cache.drafts);
            if N == 0
                warning('HAMS_Pipeline:EmptyCache', 'Hydro cache is empty.');
                return;
            end

            % Sort by draft (ascending)
            [drafts_sorted, si] = sort(cache.drafts);
            z_cg_sorted  = cache.z_cg(si);
            A_inf_sorted = cache.A_inf(si);
            B_avg_sorted = cache.B_avg(si);

            config.wamit_drafts = drafts_sorted;
            config.wamit_z_cg   = z_cg_sorted;
            config.wamit_A      = zeros(N, 3);
            config.wamit_B      = zeros(N, 3);
            config.wamit_A_full = cell(N, 1);
            config.wamit_B_full = cell(N, 1);
            config.wamit_Fe_full = cell(N, 1);   % [3×M] Fe at CG per draft

            idx_3dof = [1, 3, 5];  % surge, heave, pitch

            for i = 1:N
                A_6x6 = A_inf_sorted{i};
                B_6x6 = B_avg_sorted{i};
                z_cg_i = z_cg_sorted(i);

                % Extract 3×3 submatrix
                A_3x3 = A_6x6(idx_3dof, idx_3dof);
                B_3x3 = B_6x6(idx_3dof, idx_3dof);

                % Transform from origin to CG
                [A_cg, B_cg] = WEC_File_IO.transform_hydrodynamic_matrices_3x3( ...
                    A_3x3, B_3x3, z_cg_i);

                config.wamit_A(i, :)   = [A_cg(1,1), A_cg(2,2), A_cg(3,3)];
                config.wamit_B(i, :)   = [B_cg(1,1), B_cg(2,2), B_cg(3,3)];
                config.wamit_A_full{i} = A_cg;
                config.wamit_B_full{i} = B_cg;

                % Fe vector transform: origin → CG
                %   Force vectors (surge, heave) are invariant under reference-point
                %   shift. Only the pitch moment changes:
                %     Fe_pitch_CG = Fe_pitch_origin − z_cg × Fe_surge_origin
                %   In matrix form: Fe_CG = T' * Fe_origin
                %   where T = [1 0 -z_cg; 0 1 0; 0 0 1]  (same T used for A/B).
                %   NOTE: this is the linear (vector) form, NOT the bilinear T'MT
                %   used for A and B.
                %
                %   The raw cache Fe is at the global origin (XR=[0,0,0]).
                %   This CG-shifted version is for downstream 2D/3D EOM consumers.
                %   WEC_Visualization.plot_hydrodynamics reads the raw cache directly
                %   (origin-frame EOM) and does NOT use this field.
                % Guard: config.hydro_cache is only populated on the HAMS path.
                % On the WAMIT .1 legacy path (PATH B) rebuild_config_hydro is
                % not called, so this branch never executes there.  The isfield
                % guard makes the function safe even if called defensively.
                if isfield(config, 'hydro_cache') && ...
                        ~isempty(config.hydro_cache) && ...
                        isfield(config.hydro_cache, 'drafts') && ...
                        ~isempty(config.hydro_cache.drafts)
                    cache_local = config.hydro_cache;
                    [~, cache_i] = min(abs(cache_local.drafts - drafts_sorted(i)));
                    if isfield(cache_local, 'Fe') && ...
                            ~isempty(cache_local.Fe) && ...
                            cache_i <= length(cache_local.Fe) && ...
                            ~isempty(cache_local.Fe{cache_i})
                        Fe_6xM = cache_local.Fe{cache_i};          % [6×M] at origin
                        Fe_3xM = Fe_6xM(idx_3dof, :);              % [3×M] surge/heave/pitch
                        T_fe   = [1, 0, -z_cg_i; 0, 1, 0; 0, 0, 1];
                        config.wamit_Fe_full{i} = T_fe' * Fe_3xM; % [3×M] at CG
                    else
                        config.wamit_Fe_full{i} = [];
                    end
                else
                    config.wamit_Fe_full{i} = [];
                end
            end

            config.hydro_ready = true;
        end


        function config = retransform_at_actual_cg(config, vertical_shift, actual_cg_z)
            % RETRANSFORM_AT_ACTUAL_CG  Re-do origin→CG transform with actual CG.
            %
            %   config = HAMS_Pipeline.retransform_at_actual_cg(config, vs, cg_z)
            %
            % ── ANNOTATION CORRECTION (was incorrectly marked dead code) ──────
            % WEC_File_IO.transform_hydrodynamic_matrices_3x3 EXISTS at line 333
            % of WEC_File_IO.m and is callable.  This function is live and called
            % from WEC_Main_Optimizer (line 776) in 'trained' stage-1 mode.
            % ──────────────────────────────────────────────────────────────────
            %   rebuild_config_hydro transforms A_origin → A_CG using the
            %   uniform-density CG stored in the hydro cache.  When the
            %   optimizer uses non-uniform density, the actual CG differs.
            %
            %   This method finds the cache entry closest to vertical_shift,
            %   retrieves A_origin (6×6), re-extracts the [1,3,5] submatrix,
            %   and re-transforms to actual_cg_z.  The result overwrites
            %   the matching entry in config.wamit_A_full and config.wamit_A.
            %
            %   WHY this matters:
            %     A33_CG = A33_origin                   → exact (no CG dependence)
            %     A11_CG = A11_origin                   → exact (no CG dependence)
            %     A55_CG = A55_origin − 2z·A15 + z²A11  → quadratic in z_CG
            %     A15_CG = A15_origin − z·A11            → linear in z_CG
            %
            %   For C0 with Δz_CG ≈ 0.2 m: ~8% error on A55, ~4% on T_pitch.
            %
            %   Cost: one 3×3 matrix multiply. No HAMS re-run.

            if isempty(config.wamit_drafts); return; end

            % Find the closest draft in the interpolation table
            [~, idx] = min(abs(config.wamit_drafts - vertical_shift));

            % Retrieve A_origin and B_origin from the cache (6×6, at origin)
            cache = config.hydro_cache;
            [~, cache_idx] = min(abs(cache.drafts - vertical_shift));

            A_6x6 = cache.A_inf{cache_idx};
            B_6x6 = cache.B_avg{cache_idx};

            % Extract 3×3 [surge, heave, pitch]
            idx_3dof = [1, 3, 5];
            A_3x3_origin = A_6x6(idx_3dof, idx_3dof);
            B_3x3_origin = B_6x6(idx_3dof, idx_3dof);

            % Re-transform with ACTUAL CG from non-uniform density
            [A_cg, B_cg] = WEC_File_IO.transform_hydrodynamic_matrices_3x3( ...
                A_3x3_origin, B_3x3_origin, actual_cg_z);

            % Overwrite the matching entry
            config.wamit_A(idx, :)   = [A_cg(1,1), A_cg(2,2), A_cg(3,3)];
            config.wamit_B(idx, :)   = [B_cg(1,1), B_cg(2,2), B_cg(3,3)];
            config.wamit_A_full{idx} = A_cg;
            config.wamit_B_full{idx} = B_cg;
            % Update wamit_z_cg so §9b in calculate_3d_properties sees
            % z_cg_hams == actual_cg_z → dz = 0 → skips the delta retransform.
            % Without this, §9b would apply a second congruence transform on a
            % matrix that is already at actual_cg_z, corrupting A55 and A11.
            % This entry is reset by rebuild_config_hydro at the start of the
            % next outer iteration, so there is no cross-iteration contamination.
            config.wamit_z_cg(idx) = actual_cg_z;
        end


        %% ═══════════════════════════════════════════════════════════
        %%  DEFAULT HAMS PARAMETERS
        %% ═══════════════════════════════════════════════════════════

        function params = default_hams_params(config)
            % DEFAULT_HAMS_PARAMS  Sensible defaults for WEC analysis.
            %
            %   params = HAMS_Pipeline.default_hams_params()
            %   params = HAMS_Pipeline.default_hams_params(config)
            %
            %   Deep water.  Single heading (0 deg, head seas).
            %   Irregular frequency removal ON.
            %   zero_inf_limits ON (needed for A(inf) and A(0)).
            %
            %   FREQUENCY GRID — PERIOD-UNIFORM (Input_frequency_type = 4)
            %     T = 2 s → 20 s in 0.5 s steps  (37 frequencies)
            %     Matches the 2D / 3D / climate-analysis bands so the
            %     hydrodynamics are sampled on the SAME grid as the rest
            %     of the pipeline.
            %     Output stays in omega [rad/s] (Output_frequency_type = 3)
            %     for downstream consumption — the .1 / .3 parser expects it.
            %
            %     Period-uniform spacing → resolution is much finer at
            %     LONG periods (low omega, design band) and coarser at
            %     SHORT periods (high omega, IRFR tail).  For the IRFR
            %     ringing seen at column-section drafts (T ≈ 2.2–2.6 s),
            %     refine panel_size — adding samples to the high-omega
            %     end of THIS grid does not push the irregular frequency
            %     up; only mesh density does.
            %
            %   CONFIG OVERRIDES
            %     If `config` is supplied with any of the fields
            %     hams_T_min / hams_T_max / hams_T_step, those values
            %     replace the defaults.  WEC_Driver sets these in §2.
            %
            %   wave_diffrac_soln = 1  ->  OEXFOR/.3 file contains TOTAL
            %   excitation force (Froude-Krylov + diffraction).  Required
            %   for correct RAO computation.
            %
            %   HISTORY
            %     2026-03-22 : wave_diffrac_soln 2 → 1
            %     2026-03-24 : ω-grid narrowed to ω = 0.2:0.1:2.9 rad/s
            %     2026-05-12 : switched to PERIOD-uniform (T = 2:0.5:20 s)
            %                  to match other pipelines; config overrides
            %                  added so the range is tunable from
            %                  WEC_Driver §2 via params.hams_T_min/max/step.

            if nargin < 1, config = struct(); end

            params.depth = 74;           % [m] NA site water depth (finite); -1 => deep water
            params.zero_inf_limits = 1;  % compute A(0) and A(inf)

            % Frequency grid: period-uniform.  HAMS reads (min, step)
            % as PERIOD in seconds when input_freq_type = 4.
            params.input_freq_type  = 4;   % period [s]
            params.output_freq_type = 3;   % omega [rad/s] (parser convention)
            T_min  = 2.0;    % [s]  shortest period
            T_max  = 20.0;   % [s]  longest period
            T_step = 0.5;    % [s]  period step

            % Apply config overrides if present
            if isfield(config, 'hams_T_min')  && ~isempty(config.hams_T_min)
                T_min  = config.hams_T_min;
            end
            if isfield(config, 'hams_T_max')  && ~isempty(config.hams_T_max)
                T_max  = config.hams_T_max;
            end
            if isfield(config, 'hams_T_step') && ~isempty(config.hams_T_step)
                T_step = config.hams_T_step;
            end

            n_T = round((T_max - T_min) / T_step) + 1;
            assert(T_step > 0,     'hams_T_step must be > 0 (got %g)', T_step);
            assert(T_max > T_min,  'hams_T_max (%g) must be > hams_T_min (%g)', T_max, T_min);
            assert(n_T >= 2,       'period grid must have >=2 entries (got %d)', n_T);

            params.min_frequency = T_min;    % T_min [s]  (label says "Wmin" but it's period when type=4)
            params.freq_step     = T_step;   % ΔT [s]
            params.n_frequencies = -n_T;     % negative → uniform stepping

            % Single heading (head seas) — auto-range mode (n_headings < 0).
            % write_control_file only supports the Minimum_heading+Heading_step
            % pair (verified Fortran format). n_headings > 0 writes a bare value
            % line that mis-parses in HAMS.
            params.n_headings  = -1;    % auto-range: 1 heading
            params.min_heading  = 0.0;  % heading start [deg]
            params.heading_step = 90.0; % heading step  [deg] (irrelevant for 1 heading)

            % Reference body center (rotation centre XR).
            % Always [0,0,0]: HAMS outputs A, B, Fe at the global origin.
            % The post-processor (rebuild_config_hydro / §5 PATH A) applies
            % the single congruence transform origin → CG.
            % DO NOT override this with CG — doing so causes a double-transform.
            params.ref_body_center = [0, 0, 0];
            params.ref_body_length = 1.0;

            % Solver settings
            %   wave_diffrac_soln = 1: total excitation (FK + diffraction)
            %   wave_diffrac_soln = 2: diffraction potential only (WRONG for RAO)
            params.wave_diffrac_soln = 1;   % FIXED (was 2)
            params.remove_irr_freq = 1;     % remove irregular frequencies
            params.n_threads = 4;

            % Minimal field points (required by format but not used)
            params.n_field_points = 1;
            params.field_points = [0, 0, 0];
        end


        %% ═══════════════════════════════════════════════════════════
        %%  VALIDATION UTILITIES
        %% ═══════════════════════════════════════════════════════════

        function report = validate_ellipsoid_hydrostatics(a, c, draft)
            % VALIDATE_ELLIPSOID_HYDROSTATICS  Compare numerical vs analytical.
            %
            %   report = validate_ellipsoid_hydrostatics(a, c, draft)
            %
            %   For a prolate ellipsoid at draft=0, closed-form solutions exist:
            %     V_sub = (2/3) π a² c
            %     Aw    = π a²
            %     CB_z  = -3c/8
            %
            %   For nonzero draft, computes numerically at two resolutions
            %   and checks convergence.

            rho_w = HAMS_Pipeline.RHO_WATER;
            g = HAMS_Pipeline.G;

            fprintf('\n=== Ellipsoid Hydrostatic Validation ===\n');
            fprintf('  a = %.3f m, c = %.3f m, draft = %.3f m\n', a, c, draft);

            % Numerical at two resolutions
            [r1, z1] = HAMS_Pipeline.generate_ellipsoid_profile(a, c, draft, 50);
            hydro1 = HAMS_Pipeline.compute_geometric_hydrostatics(r1, z1, rho_w, g);

            [r2, z2] = HAMS_Pipeline.generate_ellipsoid_profile(a, c, draft, 200);
            hydro2 = HAMS_Pipeline.compute_geometric_hydrostatics(r2, z2, rho_w, g);

            report.hydro_coarse = hydro1;
            report.hydro_fine = hydro2;

            if abs(draft) < 1e-10
                % Analytical solutions for half-ellipsoid
                V_exact = (2/3) * pi * a^2 * c;
                Aw_exact = pi * a^2;
                CB_exact = -3*c/8;

                % Izz for uniform half-ellipsoid about its centroid:
                % For full ellipsoid: Izz = (1/5)m_full(a² + a²) = (2/5)m_full a²
                % For half: Izz_origin = (1/5) m_half * 2a² (same formula by symmetry)
                % Actually let me derive properly.
                % For body of revolution about z, Izz = ρ π/2 ∫ r⁴ dz
                % r = a sin(φ), z = -c cos(φ), dz = c sin(φ) dφ
                % r⁴ = a⁴ sin⁴(φ)
                % Izz = ρ (π/2) a⁴ c ∫₀^{π/2} sin⁵(φ) dφ
                %      = ρ (π/2) a⁴ c × 8/15
                %      = ρ π a⁴ c × 4/15
                % With m = ρ (2/3) π a² c:
                % Izz = m × (2a²/5)
                m_exact = rho_w * V_exact;
                Izz_exact = m_exact * (2*a^2/5);

                % Ixx = Iyy about ORIGIN for half-ellipsoid:
                % Ixx_o = ρ [π/4 ∫ r⁴ dz + π ∫ z² r² dz]
                % Term 1 = ρ π/4 a⁴ c ∫₀^{π/2} sin⁵(φ) dφ = ρ π a⁴ c × 2/15
                % Term 2 = ρ π a² c³ ∫₀^{π/2} cos²(φ) sin³(φ) dφ
                %        = ρ π a² c³ × [2/15]  (∫₀^{π/2} cos²φ sin³φ dφ = 2/15)
                % Ixx_o = ρ π [a⁴ c × 2/15 + a² c³ × 2/15]
                %       = ρ π × 2/(15) × a² c (a² + c²)
                % Ixx_cg = Ixx_o - m × CB_z²

                Ixx_origin_exact = rho_w * pi * 2/15 * a^2 * c * (a^2 + c^2);
                Ixx_cg_exact = Ixx_origin_exact - m_exact * CB_exact^2;

                report.analytical.V_sub = V_exact;
                report.analytical.Aw = Aw_exact;
                report.analytical.CB_z = CB_exact;
                report.analytical.mass = m_exact;
                report.analytical.Izz = Izz_exact;
                report.analytical.Ixx = Ixx_cg_exact;
                report.analytical.Iyy = Ixx_cg_exact;

                % Print comparison
                fprintf('\n  Property         Analytical    Numerical(50)  Numerical(200)  Err(200)\n');
                fprintf('  V_sub [m³]       %10.6f    %10.6f     %10.6f     %.2e\n', ...
                    V_exact, hydro1.V_sub, hydro2.V_sub, abs(hydro2.V_sub - V_exact)/V_exact);
                fprintf('  Aw [m²]          %10.6f    %10.6f     %10.6f     %.2e\n', ...
                    Aw_exact, hydro1.Aw, hydro2.Aw, abs(hydro2.Aw - Aw_exact)/Aw_exact);
                fprintf('  CB_z [m]         %10.6f    %10.6f     %10.6f     %.2e\n', ...
                    CB_exact, hydro1.CB_z, hydro2.CB_z, abs(hydro2.CB_z - CB_exact)/abs(CB_exact));
                fprintf('  Izz [kg·m²]      %10.2f    %10.2f     %10.2f     %.2e\n', ...
                    Izz_exact, hydro1.Izz, hydro2.Izz, abs(hydro2.Izz - Izz_exact)/Izz_exact);
                fprintf('  Iyy [kg·m²]      %10.2f    %10.2f     %10.2f     %.2e\n', ...
                    Ixx_cg_exact, hydro1.Iyy, hydro2.Iyy, abs(hydro2.Iyy - Ixx_cg_exact)/Ixx_cg_exact);

                report.passed = (abs(hydro2.V_sub - V_exact)/V_exact < 1e-4) && ...
                                (abs(hydro2.CB_z - CB_exact)/abs(CB_exact) < 1e-3);

                if report.passed
                    fprintf('\n  ✓ VALIDATION PASSED (fine-grid errors < 0.01%%)\n');
                else
                    fprintf('\n  ✗ VALIDATION FAILED — check numerical integration\n');
                end
            else
                % No closed-form — check convergence between resolutions
                V_err = abs(hydro2.V_sub - hydro1.V_sub) / hydro2.V_sub;
                CB_err = abs(hydro2.CB_z - hydro1.CB_z) / abs(hydro2.CB_z);

                fprintf('\n  Convergence check (50 vs 200 points):\n');
                fprintf('  V_sub: %.6f vs %.6f  (rel diff: %.2e)\n', ...
                    hydro1.V_sub, hydro2.V_sub, V_err);
                fprintf('  CB_z:  %.6f vs %.6f  (rel diff: %.2e)\n', ...
                    hydro1.CB_z, hydro2.CB_z, CB_err);

                report.passed = V_err < 1e-3 && CB_err < 1e-3;
                if report.passed
                    fprintf('  ✓ CONVERGED\n');
                else
                    fprintf('  ✗ NOT CONVERGED — increase resolution\n');
                end
            end
        end


        function validate_mesh_areas(nodes, panels, panel_nverts, ...
                expected_area, symmetry_factor, tol)
            % VALIDATE_MESH_AREAS  Check total mesh surface area.
            %
            %   validate_mesh_areas(nodes, panels, nverts, A_expected, sym, tol)
            %
            %   Computes total panel area and compares against expected.
            %   symmetry_factor: multiply meshed area by this to get full body.
            %     Y-symmetry only: 2
            %     X+Y symmetry: 4

            total_area = 0;
            n_panels = size(panels, 1);

            for p = 1:n_panels
                v = panels(p, 1:panel_nverts(p));
                pts = nodes(v, :);

                if panel_nverts(p) == 3
                    % Triangle area
                    edge1 = pts(2,:) - pts(1,:);
                    edge2 = pts(3,:) - pts(1,:);
                    total_area = total_area + 0.5 * norm(cross(edge1, edge2));
                else
                    % Quad area (split into two triangles)
                    edge1 = pts(2,:) - pts(1,:);
                    edge2 = pts(3,:) - pts(1,:);
                    edge3 = pts(4,:) - pts(1,:);
                    total_area = total_area + 0.5 * norm(cross(edge1, edge2));
                    total_area = total_area + 0.5 * norm(cross(edge2, edge3));
                end
            end

            full_area = total_area * symmetry_factor;
            rel_err = abs(full_area - expected_area) / expected_area;

            fprintf('\n  Mesh area validation:\n');
            fprintf('    Meshed area:   %.6f m²\n', total_area);
            fprintf('    Full area:     %.6f m² (×%d symmetry)\n', full_area, symmetry_factor);
            fprintf('    Expected area: %.6f m²\n', expected_area);
            fprintf('    Relative error: %.4f%%\n', rel_err * 100);

            if rel_err < tol
                fprintf('    ✓ PASSED (error < %.2f%%)\n', tol * 100);
            else
                fprintf('    ✗ FAILED (error > %.2f%%)\n', tol * 100);
            end
        end

        %% ═══════════════════════════════════════════════════════════
        %%  HULL WATERLINE POLYGON (from mesh open-edge topology)
        %% ═══════════════════════════════════════════════════════════

        function boundary_xy = hull_waterline_polygon(mesh, z_tol)
        % HULL_WATERLINE_POLYGON  Ordered waterline polygon from hull mesh.
        %
        %   boundary_xy = HAMS_Pipeline.hull_waterline_polygon(mesh)
        %   boundary_xy = HAMS_Pipeline.hull_waterline_polygon(mesh, z_tol)
        %
        %   Extracts the waterplane boundary polygon whose vertices are exactly
        %   the hull mesh nodes that lie on the waterline.  This guarantees that
        %   the WP outer boundary is node-coincident with the hull mesh waterline
        %   for any Nu, Nv, or wp_target_edge setting.
        %
        %   ALGORITHM
        %     1. Build an edge-count map over all mesh.panels edges.
        %     2. Collect edges where count == 1 (open / boundary edges) AND
        %        both endpoints satisfy |z| < z_tol.
        %     3. Build adjacency list from those edges.
        %     4. Walk from each unvisited node, following unvisited neighbours
        %        (with prev-pointer to avoid immediate back-step).
        %     5. Keep any chain >= 3 nodes as a loop candidate.
        %     6. Return the loop with the largest abs(polyarea) — the outer
        %        boundary.  Remove duplicate closing vertex if present.
        %
        %   PREREQUISITES
        %     WP-2 FIX in WEC_Panelizer (post-trim merge at 1e-8 tolerance)
        %     must have run.  It fuses coincident-but-index-distinct seam
        %     vertices produced by split_panel_at_z, ensuring the waterline
        %     edge graph is topologically connected.
        %
        %   INPUTS
        %     mesh   — struct from WEC_Panelizer.generate() with trim_wl=true
        %     z_tol  — [m] half-band around z=0 for waterline detection
        %              (default 0.01 m; conservative for all expected geometries)
        %
        %   OUTPUT
        %     boundary_xy  [N×2]  ordered (x,y) polygon, no repeated endpoint.
        %                  Empty [0×2] with a warning if extraction fails.
        %
        %   FULL / HALF BODY
        %     Full body (x_sym=0, y_sym=0): open edges form a closed ring.
        %     Half body (y_sym=1, x_sym=0): open edges include the semicircular
        %       arc PLUS the y=0 symmetry-plane closure edges at z=0.  Together
        %       they form a closed half-disk polygon — correct for CDT meshing
        %       with HAMS y-symmetry mirroring.
        %
        %   See also: mesh_wp_blossomquad, run_single_hams, run_draft_sweep

            if nargin < 2 || isempty(z_tol), z_tol = 1e-6; end

            verts = mesh.vertices;
            n_p   = mesh.n_panels;

            % WL-FIX (2026-05-06): The waterline vertex set is the UNION of
            %   (a) vertices flagged by trim_at_wl/split_panel_at_z, AND
            %   (b) original mesh vertices that already sit at z = z_wl
            %       (within machine precision)
            % Why both: when one surface's parametric grid happens to land
            % a row at exactly z_wl, no trim runs for that surface, so its
            % z=0 vertices are NOT in waterline_verts — but they ARE on
            % the true waterline boundary. A second surface whose grid
            % straddles z_wl gets trimmed and its intersections ARE flagged.
            % Using only one source misses half the perimeter (Nu=23 / Nu=45
            % on C0: 22% / 12% area deficit). The previous |z|<0.01 heuristic
            % was too generous (swept in interior rows); 1e-6 is tight enough
            % that only true z=0 nodes pass.
            wl_set = abs(verts(:,3)) < z_tol;
            if isfield(mesh,'waterline_verts') && ~isempty(mesh.waterline_verts)
                wl_set(mesh.waterline_verts) = true;
            end
            on_wl = @(vi) wl_set(vi);

            % ── Step 1: Build edge → count map ────────────────────────────
            %  Key: 'minV_maxV' (unordered pair).  Count = number of panels
            %  sharing that edge.  Open (boundary) edges have count == 1.
            edge_cnt = containers.Map('KeyType','char','ValueType','int32');
            edge_ep  = containers.Map('KeyType','char','ValueType','any');

            for p = 1:n_p
                v  = mesh.panels(p,:);
                nv = 4;
                if v(3) == v(4), nv = 3; end
                for e = 1:nv
                    va  = v(e);
                    vb  = v(mod(e, nv) + 1);
                    key = sprintf('%d_%d', min(va,vb), max(va,vb));
                    if edge_cnt.isKey(key)
                        edge_cnt(key) = edge_cnt(key) + int32(1);
                    else
                        edge_cnt(key) = int32(1);
                        edge_ep(key)  = [va, vb];
                    end
                end
            end

            % ── Step 2: Collect open waterline edges ───────────────────────
            %  An edge qualifies when it is open (count==1) and both endpoints
            %  are members of the waterline vertex set (above).  This captures:
            %    - The waterline arc (hull trim boundary).
            %    - For half-body: the y=0 symmetry-plane closure at z=0.
            %  It excludes interior edges, submerged boundary edges, and
            %  near-z=0-but-not-on-trim edges that the old |z|<0.01 test let
            %  through.
            bnd_a = [];
            bnd_b = [];
            ks = edge_cnt.keys();
            for k = 1:numel(ks)
                if edge_cnt(ks{k}) ~= 1, continue; end
                ev = edge_ep(ks{k});
                if on_wl(ev(1)) && on_wl(ev(2))
                    bnd_a(end+1) = ev(1); %#ok<AGROW>
                    bnd_b(end+1) = ev(2); %#ok<AGROW>
                end
            end

            if isempty(bnd_a)
                warning('HAMS_Pipeline:HullWLEmpty', ...
                    'No hull waterline open edges at |z|<%.3f m. Check trim_wl=true and WP-2 merge.', z_tol);
                boundary_xy = zeros(0, 2);
                return;
            end

            % ── Step 3: Build adjacency list ───────────────────────────────
            all_bv = unique([bnd_a, bnd_b]);
            n_bv   = numel(all_bv);
            mx     = max(all_bv) + 1;
            g2l    = zeros(mx, 1, 'int32');
            for i = 1:n_bv
                g2l(all_bv(i)) = int32(i);
            end

            adj = cell(n_bv, 1);
            for e = 1:numel(bnd_a)
                la = g2l(bnd_a(e));
                lb = g2l(bnd_b(e));
                adj{la}(end+1) = lb;
                adj{lb}(end+1) = la;
            end

            % ── Step 4: Walk chains ────────────────────────────────────────
            %  Standard prev-pointer walk.  Terminates when all neighbours of
            %  the current node are already visited or only the previous node
            %  remains.  For a closed loop the walk naturally stops when it
            %  returns to a node whose only unvisited neighbour is the starting
            %  node — which is already marked visited — so the loop terminates.
            %  Maximum iterations = n_bv per call to this inner loop, so there
            %  is no infinite-loop risk.
            visited = false(n_bv, 1);
            loops   = {};

            for sv = 1:n_bv
                if visited(sv), continue; end

                loop = sv;
                visited(sv) = true;
                prev = -1;
                cur  = sv;

                while true
                    nbrs = adj{cur};
                    unv  = nbrs(~visited(nbrs));
                    if ~isempty(unv) && prev > 0
                        unv = unv(unv ~= prev);
                    end
                    if isempty(unv), break; end
                    prev = cur;
                    cur  = unv(1);
                    visited(cur) = true;
                    loop(end+1) = cur; %#ok<AGROW>
                end

                if numel(loop) >= 3
                    loops{end+1} = verts(all_bv(loop), 1:2); %#ok<AGROW>
                end
            end

            % ── Step 5: Select outer boundary ─────────────────────────────
            if isempty(loops)
                warning('HAMS_Pipeline:HullWLNoLoop', ...
                    'Edge walk produced no loop >= 3 nodes. Verify WP-2 FIX post-trim merge ran in WEC_Panelizer.');
                boundary_xy = zeros(0, 2);
                return;
            end

            ar = cellfun(@(L) abs(polyarea(L(:,1), L(:,2))), loops);
            [~, si] = max(ar);
            boundary_xy = loops{si};

            % Remove duplicate closing vertex (some walks append start again)
            if norm(boundary_xy(1,:) - boundary_xy(end,:)) < 1e-10
                boundary_xy = boundary_xy(1:end-1, :);
            end

            fprintf('    Hull WL polygon: %d nodes, Aw=%.4f m² (open-edge mesh walk)\n', ...
                size(boundary_xy,1), abs(polyarea(boundary_xy(:,1), boundary_xy(:,2))));
        end


    end  % methods (Static)


    methods (Static, Access = private)

        function write_6x6_matrix(fid, M)
            % WRITE_6X6_MATRIX  Write a 6×6 matrix in HAMS Hydrostatic.in format.
            %
            %   Fortran reads with format 6(2x,E12.5):
            %     2x  = skip 2 characters
            %     E12.5 = read 12 characters as scientific notation, 5 decimals
            %   Total per value: 14 characters.  Total per line: 84 characters.
            %
            %   MATLAB %12.5E produces exactly 12 characters:
            %     Positive: " 1.45760E+04" (1 leading space + 11 significant)
            %     Negative: "-6.70615E+03" (sign + 11 significant)
            %   Prefixed with '  ' (2 spaces): 14 chars per value.  Exact match.
            %
            %   PREVIOUS BUG: %.5E produced 11 chars (positive) or 12 chars
            %   (negative).  With 3-space gaps, the total width was correct
            %   for most cases but could misalign if two consecutive values
            %   were both negative.  Fixed 2026-04-01.
            for i = 1:6
                fprintf(fid, '  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E\n', ...
                    M(i,1), M(i,2), M(i,3), M(i,4), M(i,5), M(i,6));
            end
        end


        function B_avg = compute_B_avg(hams_data, T_min, T_max)
            % COMPUTE_B_AVG  Average B(omega) over period band [T_min, T_max].
            %
            %   Converts omega to period, masks, averages.  Returns 6x6.
            %   Falls back to closest frequency if no data in band.

            omega = hams_data.omega;
            T = 2 * pi ./ omega;
            mask = (T >= T_min) & (T <= T_max);

            if ~any(mask)
                T_center = 0.5 * (T_min + T_max);
                [~, closest] = min(abs(T - T_center));
                mask(closest) = true;
            end

            B_avg = mean(hams_data.B(:, :, mask), 3);
        end

    end  % methods (Static, Access = private)

end  % classdef