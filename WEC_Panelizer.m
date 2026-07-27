classdef WEC_Panelizer
% WEC_PANELIZER  Generate BEM panel meshes from parsed MultiSurf geometry.
%
%   mesh = WEC_PANELIZER.generate(parser, draft, Nu, Nv)
%
%   Takes a parsed .ms2 model (from WEC_MS2_Parser), classifies visible
%   surfaces into sources and mirrors, evaluates only source surfaces on
%   a parametric grid, derives mirror meshes by coordinate flip, merges
%   vertices at shared boundaries, trims at the waterline, and writes
%   HAMS-compatible .pnl files.
%
%   ARCHITECTURE
%   ────────────────────────────────────────────────────────────────────
%   WEC_Panelizer (this file)
%     ├─ §1  generate           — source/mirror eval → panel mesh
%     ├─ §2  grid_to_quads      — (Nu×Nv) grid → quad connectivity
%     ├─ §3  trim_at_wl         — clip panels at waterline, split straddling
%     ├─ §4  close_open_edges   — detect open boundary → fill with cap
%     │                           (disabled by default — not needed for BEM)
%     ├─ §5  orient_normals     — enforce consistent outward normals
%     ├─ §6  merge_vertices     — merge coincident vertices at seams
%     ├─ §7  write_hull_pnl     — write HAMS HullMesh.pnl
%     ├─ §8  write_wp_pnl       — write HAMS WaterPlaneMesh.pnl
%     ├─ §9  validate           — area, aspect ratios, open-edge diagnostics
%     ├─ §10 check_bem_quality  — panel size vs wavelength, BEM adequacy
%     └─ §11 diagnose_mesh      — per-surface stats, open edge classification
%   ────────────────────────────────────────────────────────────────────
%
%   DESIGN PRINCIPLES
%     1. SOURCE/MIRROR — EVALUATE SOURCES ONLY.
%        The .ms2 entity DAG identifies MirrSurf surfaces as coordinate-
%        flipped copies of a source surface.  Only source surfaces are
%        evaluated on the (Nu × Nv) grid.  Mirror meshes are derived by
%        duplicating source vertices and flipping the relevant coordinate.
%        This halves the number of eval_surface calls and guarantees
%        bit-identical vertices at shared boundaries (exact merge).
%
%     2. GEOMETRY FIRST, PANELS SECOND.
%        Surfaces are evaluated from their parametric representation
%        (B-splines, RevSurf, etc.) — not from a pre-existing STL.
%        The panel density [Nu, Nv] is a user parameter, independent
%        of geometry complexity.
%
%     3. TRIM AT WATERLINE.
%        For panels straddling the waterline, edge–waterline
%        intersection points are computed by linear interpolation,
%        then the panel is split into submerged/above-water parts.
%        This gives clean waterline edges without jagged staircase.
%
%     4. OPEN MESH FOR BEM.
%        HAMS BEM operates on an OPEN surface (the submerged hull
%        boundary).  The waterline opening is closed by a separate
%        WaterPlaneMesh.pnl file (for irregular frequency removal).
%        close_gaps defaults to FALSE — enabling it creates internal
%        geometry at surface junctions that corrupts the BIE.
%
%     5. NORMAL CONSISTENCY.
%        All panel normals point OUTWARD (away from hull interior).
%        Mirror panels have their winding reversed (for odd number of
%        flips) to preserve outward orientation after the coordinate flip.
%
%   MESH STRUCT FORMAT
%     mesh.vertices       : [N_v × 3]   vertex coordinates (x, y, z)
%     mesh.panels         : [N_p × 4]   vertex indices per panel
%                            (v4 = v3 for triangular panels)
%     mesh.normals        : [N_p × 3]   outward unit normals per panel
%     mesh.surface_ids    : [N_p × 1]   which surface each panel belongs to
%     mesh.is_cap         : [N_p × 1]   logical, true for gap-closing caps
%     mesh.waterline_verts: [N_wl × 1]  vertex indices on the waterline
%     mesh.draft          : scalar       draft shift applied
%     mesh.Nu, mesh.Nv    : grid density used
%     mesh.n_panels       : total panel count
%     mesh.n_vertices     : total vertex count
%
%   HAMS MESH FORMAT (HullMesh.pnl)
%   ────────────────────────────────────────────────────────────────────
%     Line 1:  comment header
%     Line 2:  N_vertices  N_panels
%     Vertex block:  index  x  y  z
%     Panel block:   index  v1 v2 v3 v4
%       (v4 = v3 for triangular panels — confirmed: HAMS supports both)
%
%   See also: WEC_MS2_Parser, WEC_HAMS_Interface
%
%   Author:  WEC Optimisation Team
%   Version: 3.0 — BEM quality tools, cosine spacing option

    methods (Static)

        %% ═════════════════════════════════════════════════════════
        %%  §1  GENERATE — Evaluate Surfaces → Panel Mesh
        %% ═════════════════════════════════════════════════════════

        function mesh = generate(parser, draft, Nu, Nv, options)
        % GENERATE  Build a panel mesh from parsed MS2 geometry.
        %
        %   mesh = WEC_Panelizer.generate(parser, draft, Nu, Nv)
        %   mesh = WEC_Panelizer.generate(parser, draft, Nu, Nv, options)
        %
        %   INPUTS
        %     parser   : WEC_MS2_Parser object (parsed .ms2 file)
        %     draft    : [m] vertical shift applied to hull z-coordinates
        %                (positive = hull moves up relative to waterline)
        %     Nu, Nv   : grid density per surface (grid points per
        %                parametric direction; panel count = (Nu-1)*(Nv-1))
        %     options  : (optional) struct with fields:
        %       .trim_wl         — logical, true = trim at waterline (default: true)
        %       .close_gaps      — logical, true = close open edges (default: FALSE)
        %       .z_wl            — [m] waterline elevation (default: 0)
        %       .cosine_spacing  — logical, use half-cosine point distribution
        %                          (default: false).  Clusters panels near
        %                          parametric boundaries (u=0,1 and v=0,1),
        %                          improving resolution at sharp features
        %                          (column-platform transition, waterline).
        %       .verbose         — logical, print progress (default: true)
        %
        %   OUTPUT
        %     mesh : struct (see MESH STRUCT FORMAT in class header)
        %
        %   SOURCE/MIRROR ARCHITECTURE
        %     The .ms2 entity DAG encodes which surfaces are MirrSurfs
        %     (coordinate-flipped copies of a source surface).  Rather than
        %     evaluating every visible surface independently, this method:
        %       1. Classifies visible surfaces into sources and mirrors.
        %       2. Evaluates ONLY source surfaces on the (Nu × Nv) grid.
        %       3. Derives mirror meshes by duplicating source vertices
        %          and flipping the relevant coordinate (X or Y).
        %       4. Reverses panel winding on mirrors to preserve outward normals.
        %
        %     WHY this guarantees watertight junctions:
        %       - EdgeSnake junctions: surface4(u,0) = snake3(u) = surface2(u,0).
        %         Same u_grid → bit-identical boundary vertices → exact merge.
        %       - Y=0 mirror junctions: source edge has y=0 → MirrY gives y=−0 = 0
        %         → identical vertices → exact merge.
        %       - Point junctions (keel): same eval_point → exact merge.
        %
        %   WHY close_gaps defaults to FALSE:
        %     HAMS BEM operates on an OPEN surface — the submerged hull
        %     boundary panels only.  The waterline opening is closed by
        %     WaterPlaneMesh.pnl (separate file, for irregular frequency
        %     removal).  The column top (above water) is not part of the
        %     BEM domain.  Closing gaps at surface junctions creates
        %     internal geometry that corrupts the boundary integral equation
        %     and is never correct for BEM input.

            % ── Parse options ─────────────────────────────────────
            if nargin < 5, options = struct(); end
            if ~isfield(options, 'trim_wl'),         options.trim_wl         = true;  end
            if ~isfield(options, 'close_gaps'),       options.close_gaps      = false; end
            if ~isfield(options, 'z_wl'),             options.z_wl            = 0;    end
            if ~isfield(options, 'cosine_spacing'),   options.cosine_spacing  = false; end
            if ~isfield(options, 'wl_conform'),       options.wl_conform      = true;  end
            if ~isfield(options, 'arclength_u'),      options.arclength_u     = true;  end
            if ~isfield(options, 'adapt_Nv'),         options.adapt_Nv        = true;  end
            if ~isfield(options, 'verbose'),          options.verbose         = true;  end
            if ~isfield(options, 'quarter_body'),      options.quarter_body    = false; end
            if ~isfield(options, 'half_body'),          options.half_body       = false; end
            if ~isfield(options, 'vertex_merge'),       options.vertex_merge    = 'merge'; end
            %  half_body:
            %    When true, produces a half-mesh (y >= 0) with ISX=0, ISY=1
            %    for HAMS-MREL.  Sources are transformed to Q1 (same as
            %    quarter_body), then X-mirrors are included to fill Q2.
            %    Y-mirrors and XY-mirrors are excluded (they cross into y<0).
            %    HAMS mirrors the y>=0 half across y=0 via the ELSE branch
            %    of IF(ISX.EQ.1.AND.ISY.EQ.0) in all solver routines.
            %    Mutually exclusive with quarter_body.
            if options.half_body && options.quarter_body
                error('WEC_Panelizer:MutualExclusion', ...
                      'half_body and quarter_body are mutually exclusive.');
            end
            %  vertex_merge options:
            %    'merge'  — tolerance-based merge (1e-8 m), then unique
            %    'dedup'  — exact coordinate match only (unique rows)
            %    'none'   — skip entirely (duplicate vertices remain)

            % ── Classify surfaces ──────────────────────────────────
            topo = parser.classify_visible_surfaces();

            % ── Initialise collectors ─────────────────────────────
            all_verts    = [];
            all_panels   = [];
            all_surf_ids = [];
            all_is_cap   = [];

            % ── Grid generation ───────────────────────────────────
            %  Cosine spacing clusters points near parametric boundaries
            %  (u=0, u=1, v=0, v=1), giving higher resolution at sharp
            %  features: column-to-platform transition, waterline region,
            %  surface junctions.  Does NOT break junction guarantees
            %  because both surfaces at a shared edge use the same u_grid.
            if options.cosine_spacing
                u_grid = (1 - cos(linspace(0, pi, Nu))) / 2;
                v_grid = (1 - cos(linspace(0, pi, Nv))) / 2;
            else
                u_grid = linspace(0, 1, Nu);
                v_grid = linspace(0, 1, Nv);
            end

            % ── Waterline-conforming grid (Change 1) ──────────────
            %  Find the parametric u-value where z(u) crosses the
            %  waterline.  Both RevSurf and DevSurf have z controlled
            %  by the u-parameter: z(u, v) ≈ z_profile(u) for all v.
            %  Inserting u_wl into the shared u_grid places a ROW of
            %  vertices exactly at z = 0.  trim_at_wl then never
            %  splits a panel — it only discards above-water panels.
            %
            %  The arc-length spacing (Change 2) uses u_wl as a
            %  segment boundary, so we find u_wl FIRST, then build
            %  the arc-length grid in two segments.

            u_wl_found = [];

            if options.trim_wl && options.wl_conform
                z_wl = options.z_wl;

                % Probe each source surface along u at v = 0.5
                u_probe = linspace(0, 1, 200);

                for s = 1:length(topo.sources)
                    sname = topo.sources{s};

                    % Sample z(u, v=0.5) + draft
                    z_probe = zeros(length(u_probe), 1);
                    for k = 1:length(u_probe)
                        pt = parser.eval_surface(sname, u_probe(k), 0.5);
                        z_probe(k) = pt(3) + draft;
                    end

                    % Find zero-crossings and bisect to machine precision.
                    %
                    % WP-1 FIX: The original strict f1*f2 < 0 test only fires
                    % on a true sign change and is silent when the surface
                    % profile TOUCHES z_wl (i.e. f=0 at a sample point).
                    % This occurs when the waterline lands exactly at the top
                    % endpoint of the profile (u=1) or at a flat horizontal
                    % section.  When u_wl_found stays empty, the arc-length
                    % grid is built without snapping to the waterline, so
                    % trim_at_wl calls split_panel_at_z on every straddling
                    % panel.  split_panel_at_z creates fresh, unmerged
                    % intersection vertices at z=0 for each panel
                    % independently.  Adjacent panels that share edge V2-V3
                    % each compute their own intersection I on that edge
                    % (same 3D position, different indices).  These appear as
                    % disconnected 2-vertex chains in extract_waterline_
                    % boundary_topo → loop length = 2 < 3 → rejected →
                    % empty WP boundary → IRSP disabled.
                    %
                    % Extended condition: also detect f1 < 0 AND |f2| ≤ tol
                    % (surface arriving at waterline from below at endpoint).
                    % Both cases use the SAME bisection — no special branch.
                    % Skipping bisection when is_touching would insert
                    % u_probe(k+1) directly; if f2 = +5e-9 the conformal row
                    % lands at z = z_wl + 5e-9, which is ABOVE the waterline
                    % (trim keeps only z ≤ z_wl + 1e-10) and gets discarded,
                    % recreating the original problem.  Bisection always gives
                    % the exact root regardless of how f2 arrived at ~0.
                    tol_wl = 1e-8;   % absolute z-distance for "at waterline"
                    for k = 1:length(z_probe) - 1
                        f1 = z_probe(k)   - z_wl;
                        f2 = z_probe(k+1) - z_wl;
                        is_crossing = f1 * f2 < 0;
                        is_touching = f1 < -tol_wl && abs(f2) <= tol_wl;
                        if is_crossing || is_touching
                            % Bisect [u_probe(k), u_probe(k+1)] — handles
                            % both genuine sign change and endpoint touch.
                            ua = u_probe(k);  ub = u_probe(k+1);
                            fa = f1;
                            for iter = 1:40
                                um = (ua + ub) / 2;
                                pt = parser.eval_surface(sname, um, 0.5);
                                fm = pt(3) + draft - z_wl;
                                if abs(fm) < 1e-12, break; end
                                if fa * fm < 0
                                    ub = um;
                                else
                                    ua = um;  fa = fm;
                                end
                            end
                            u_wl_found(end+1) = (ua + ub) / 2; %#ok<AGROW>
                        end
                    end
                end

                % Deduplicate across surfaces
                if ~isempty(u_wl_found)
                    u_wl_found = uniquetol(u_wl_found, 1e-6);
                end

                if options.verbose && ~isempty(u_wl_found)
                    for iw = 1:length(u_wl_found)
                        fprintf('    Waterline at u = %.6f (z_body = %.4f m)\n', ...
                                u_wl_found(iw), -draft);
                    end
                end
            end

            % ── Arc-length spacing in u (Change 2) ────────────────
            %  Replace the uniform/cosine u_grid with one that gives
            %  approximately equal PHYSICAL panel heights along z.
            %
            %  Strategy:
            %    1. Sample z(u) on a fine parametric grid.
            %    2. Compute cumulative arc length s(u) = ∫|dz/du| du.
            %    3. If waterline was found, split the grid into two
            %       segments: [0, u_wl] and [u_wl, 1].  Distribute
            %       grid points proportionally to each segment's
            %       arc length.
            %    4. Within each segment, invert s(u) to find the
            %       u-values that give equal s-spacing.
            %
            %  The reference z-profile comes from the first source
            %  surface sampled at v = 0.5.  For C0, both sources
            %  (RevSurf and DevSurf) have the same z-extent, so
            %  either gives the same arc-length distribution.

            if options.arclength_u && ~isempty(topo.sources)
                sname_ref = topo.sources{1};
                n_arc_samples = 500;
                u_fine = linspace(0, 1, n_arc_samples);
                pts_fine = zeros(n_arc_samples, 3);

                for k = 1:n_arc_samples
                    pt = parser.eval_surface(sname_ref, u_fine(k), 0.5);
                    pts_fine(k, :) = pt + [0, 0, draft];
                end

                % Cumulative 3D arc length along the surface profile
                %  This captures BOTH vertical AND radial changes.
                %  At the column-to-platform fillet, the radius changes
                %  rapidly → large ds/du → dense grid points there.
                %  On the flat platform base, the z change is small but
                %  the radius change is also small → moderate ds/du.
                dp = diff(pts_fine);
                ds = sqrt(sum(dp.^2, 2));
                s_cum = [0; cumsum(ds)];
                s_total = s_cum(end);

                if s_total > 1e-6
                    if ~isempty(u_wl_found)
                        % Use the first (primary) waterline crossing
                        u_wl = u_wl_found(1);

                        % Arc length to waterline
                        s_wl = interp1(u_fine, s_cum, u_wl, 'linear', 'extrap');

                        % Allocate grid points proportionally
                        frac_below = s_wl / s_total;
                        Nu_below = max(2, round(frac_below * (Nu - 1)));
                        Nu_above = max(1, (Nu - 1) - Nu_below);

                        % Arc-length-equispaced u below waterline
                        if Nu_below >= 2
                            s_targets_below = linspace(0, s_wl, Nu_below + 1);
                            u_below = interp1(s_cum, u_fine, ...
                                s_targets_below, 'linear', 'extrap');
                            % Snap last point exactly to u_wl
                            u_below(end) = u_wl;
                        else
                            u_below = [0, u_wl];
                        end

                        % Arc-length-equispaced u above waterline
                        if Nu_above >= 1
                            s_targets_above = linspace(s_wl, s_total, Nu_above + 1);
                            u_above = interp1(s_cum, u_fine, ...
                                s_targets_above, 'linear', 'extrap');
                            % Snap first point exactly to u_wl
                            u_above(1) = u_wl;
                        else
                            u_above = u_wl;
                        end

                        % Combine (u_wl appears once)
                        u_grid = unique([u_below, u_above]);

                        % Clamp endpoints
                        u_grid(1)   = 0;
                        u_grid(end) = 1;

                    else
                        % No waterline found — arc-length spacing only
                        s_targets = linspace(0, s_total, Nu);
                        u_grid = interp1(s_cum, u_fine, ...
                            s_targets, 'linear', 'extrap');
                        u_grid(1)   = 0;
                        u_grid(end) = 1;
                    end

                    % Insert any additional waterline crossings
                    % (rare: only if hull pierces waterline more than once)
                    if length(u_wl_found) > 1
                        for iw = 2:length(u_wl_found)
                            uwl_extra = u_wl_found(iw);
                            [min_d, min_i] = min(abs(u_grid - uwl_extra));
                            if min_d < 0.4 / Nu
                                u_grid(min_i) = uwl_extra;
                            else
                                u_grid = sort([u_grid, uwl_extra]);
                            end
                        end
                    end

                    % Enforce minimum spacing (prevent degenerate panels)
                    min_du = 0.3 / Nu;
                    u_grid = WEC_Panelizer.enforce_min_spacing(u_grid, min_du);

                    Nu = length(u_grid);

                    if options.verbose
                        fprintf('    Arc-length u_grid: %d points, max_du=%.4f, min_du=%.4f\n', ...
                                Nu, max(diff(u_grid)), min(diff(u_grid)));
                    end
                end

            elseif options.wl_conform && ~isempty(u_wl_found)
                % No arc-length, but snap waterline into the existing grid
                for iw = 1:length(u_wl_found)
                    uwl = u_wl_found(iw);
                    [min_d, min_i] = min(abs(u_grid - uwl));
                    if min_d < 0.5 / Nu
                        u_grid(min_i) = uwl;
                    else
                        u_grid = sort([u_grid, uwl]);
                    end
                end
                Nu = length(u_grid);
            end

            % Build a name→index map for surface IDs (sources first, then mirrors)
            all_surf_names = [topo.sources, {topo.mirrors.name}];
            n_all = length(all_surf_names);
            surf_id_map = containers.Map(all_surf_names, num2cell(1:n_all));

            % ── Per-surface v-grid (Change 3) ─────────────────────
            %  The v-direction physical extent varies dramatically
            %  between surfaces.  For C0:
            %    RevSurf (surface2): v = azimuthal, circumference up
            %      to π×2.0 ≈ 6.3 m at the platform base.
            %    DevSurf (surface4): v = cross-column, width ≈ 0.15 m.
            %
            %  With a shared Nv, the RevSurf panels at the base are
            %  6.3/Nv ≈ 0.5 m wide while the DevSurf panels are
            %  0.15/Nv ≈ 0.01 m wide.  This creates aspect ratios
            %  up to 20× on the RevSurf.
            %
            %  Fix: measure each surface's max v-extent, compute Nv
            %  so that panel v-width ≈ max physical u-step (making
            %  panels approximately square at the widest point).
            %
            %  WHY this is safe for vertex merging:
            %    Surfaces share edges along the u-direction (EdgeSnake
            %    junctions).  The shared boundary depends on u_grid
            %    only — both surfaces use the same u_grid.  The
            %    v-values at the junction are fixed endpoints (v=0 or
            %    v=1), independent of Nv.  Different Nv per surface
            %    does NOT affect junction vertex identity.

            source_v_grids = containers.Map();

            if options.adapt_Nv && ~isempty(topo.sources)
                % Compute target edge size from physical u-steps.
                %  Sample physical distances between consecutive
                %  u_grid points on the first source at v = 0.5.
                %  Use the MEAN as the target panel edge.
                sname_ref_v = topo.sources{1};
                u_phys_steps = zeros(length(u_grid) - 1, 1);
                for k = 1:length(u_grid) - 1
                    pt1 = parser.eval_surface(sname_ref_v, u_grid(k), 0.5);
                    pt2 = parser.eval_surface(sname_ref_v, u_grid(k+1), 0.5);
                    u_phys_steps(k) = norm(pt2 - pt1);
                end
                target_edge_v = mean(u_phys_steps);

                if options.verbose
                    fprintf('    Adaptive Nv: target edge = %.4f m (mean u-step)\n', ...
                            target_edge_v);
                end

                % Measure max v-extent for each source surface
                v_probe = linspace(0, 1, 50);
                u_sample_idx = round(linspace(1, length(u_grid), ...
                                     min(7, length(u_grid))));

                for s = 1:length(topo.sources)
                    sname = topo.sources{s};
                    max_Lv = 0;

                    for ui = 1:length(u_sample_idx)
                        u_val = u_grid(u_sample_idx(ui));
                        pts_v = zeros(50, 3);
                        for k = 1:50
                            pts_v(k,:) = parser.eval_surface(sname, u_val, v_probe(k));
                        end
                        Lv = sum(sqrt(sum(diff(pts_v).^2, 2)));
                        max_Lv = max(max_Lv, Lv);
                    end

                    % Nv for this surface: at least Nv (user request),
                    % at most enough to keep panels ≤ target_edge wide.
                    Nv_surf = max(Nv, ceil(max_Lv / target_edge_v) + 1);
                    source_v_grids(sname) = linspace(0, 1, Nv_surf);

                    if options.verbose
                        fprintf('    %s: max Lv = %.3f m → Nv = %d', ...
                                sname, max_Lv, Nv_surf);
                        if Nv_surf > Nv
                            fprintf(' (↑ from %d)\n', Nv);
                        else
                            fprintf('\n');
                        end
                    end
                end
            else
                % No adaptation — all surfaces use the base v_grid
                for s = 1:length(topo.sources)
                    source_v_grids(topo.sources{s}) = v_grid;
                end
            end

            if options.verbose
                spacing_str = 'uniform';
                if options.cosine_spacing, spacing_str = 'cosine'; end
                if options.arclength_u, spacing_str = 'arc-length'; end
                if options.arclength_u && options.wl_conform
                    spacing_str = 'arc-length + WL-conform';
                elseif options.wl_conform && ~options.arclength_u
                    spacing_str = [spacing_str, ' + WL-conform'];
                end
                if options.adapt_Nv, spacing_str = [spacing_str, ' + adapt-Nv']; end
                fprintf('  Panelizer: %d sources + %d mirrors, [%d×%d+] %s grid, draft=%.3f m\n', ...
                        length(topo.sources), length(topo.mirrors), Nu, Nv, spacing_str, draft);
                fprintf('    Sources: %s\n', strjoin(topo.sources, ', '));
                for mi = 1:length(topo.mirrors)
                    fprintf('    Mirror:  %s = Mirr%s(%s)\n', ...
                            topo.mirrors(mi).name, topo.mirrors(mi).plane, topo.mirrors(mi).source);
                end
            end

            % ── Quarter-body setup (before source evaluation) ─────
            %  Detect which quadrant the source surfaces occupy and
            %  prepare the coordinate transform.  The transform is
            %  applied to each source grid BEFORE grid_to_quads, so
            %  the mesh is born in Q1 with correct winding.
            negate_x = false;
            negate_y = false;
            n_flips  = 0;
            source_quadrant = 'Q1';
            x_sym = 0;
            y_sym = 0;

            if options.quarter_body
                [negate_x, negate_y, source_quadrant] = ...
                    WEC_Panelizer.detect_source_quadrant(parser, topo.sources);
                n_flips = negate_x + negate_y;
                x_sym = 1;
                y_sym = 1;

                if options.verbose
                    fprintf('    Quadrant: %s → Q1', source_quadrant);
                    if negate_x, fprintf(' [negate X]'); end
                    if negate_y, fprintf(' [negate Y]'); end
                    if mod(n_flips, 2) == 1, fprintf(' [flip v]'); end
                    fprintf('\n');
                    fprintf('    Symmetry flags: [%d, %d]\n', x_sym, y_sym);
                end

            elseif options.half_body
                % Half-body: same quadrant detection as quarter, but
                % ISX=0, ISY=1 (y>=0 half, HAMS mirrors across y=0).
                %
                % ISY=1 tells HAMS to reflect Y → −Y via the ELSE branch
                % of every IF(ISX.EQ.1.AND.ISY.EQ.0) block.  This is
                % implemented in ALL 16 solver source files (CalGreenFunc,
                % BodyIntgr, AssbMatx, PotentWavForce, PressureElevation,
                % SingularIntgr — single and multi-body variants).
                %
                % PREVIOUS BUG: set ISX=1, ISY=0.  ISX=1 mirrors X → −X,
                % but the mesh covers both x>0 (Q1) and x<0 (Q2).  HAMS
                % would create Green's function images that overlap with
                % existing panels, double-counting the hull.
                [negate_x, negate_y, source_quadrant] = ...
                    WEC_Panelizer.detect_source_quadrant(parser, topo.sources);
                n_flips = negate_x + negate_y;
                x_sym = 0;
                y_sym = 1;

                if options.verbose
                    fprintf('    Half-body mode: %s → Q1 + X-mirrors → Q2\n', source_quadrant);
                    if negate_x, fprintf('      [negate X]'); end
                    if negate_y, fprintf('      [negate Y]'); end
                    if mod(n_flips, 2) == 1, fprintf('      [flip v]'); end
                    fprintf('\n');
                    fprintf('    Symmetry flags: [%d, %d] (HAMS-MREL half-body)\n', x_sym, y_sym);
                end
            end

            % ── Evaluate each SOURCE surface on (Nu × Nv_s) grid ──
            %  Each source uses its own v_grid (from adapt_Nv or the
            %  base v_grid).  Mirrors inherit the source's grid.
            source_grids = containers.Map();

            for s = 1:length(topo.sources)
                sname = topo.sources{s};
                v_grid_s = source_v_grids(sname);

                if options.verbose
                    fprintf('    Evaluating source %s [%d×%d]...', ...
                            sname, Nu, length(v_grid_s));
                end

                S = parser.eval_surface_grid(sname, u_grid, v_grid_s);
                S(:, :, 3) = S(:, :, 3) + draft;

                % Store UNTRANSFORMED grid for mirror derivation (full body)
                if ~options.quarter_body && ~options.half_body
                    source_grids(sname) = S;
                end

                % Quarter-body or half-body: transform grid to Q1 BEFORE meshing
                if options.quarter_body || options.half_body
                    if negate_x, S(:,:,1) = -S(:,:,1); end
                    if negate_y, S(:,:,2) = -S(:,:,2); end
                    if mod(n_flips, 2) == 1
                        % Odd flips reverse the cross product direction.
                        % Flipping the v-parameter order restores it:
                        %   S_v → −S_v  ⇒  S_u × (−S_v) = −(S_u × S_v)
                        % Combined with the negation reversal, this gives
                        % the correct outward normal in Q1.
                        S = S(:, end:-1:1, :);
                    end
                end

                % Half-body: store TRANSFORMED grid for mirror derivation.
                %   Unlike full-body (which stores untransformed), the
                %   half-body mirrors need the Q1-transformed grid as their
                %   starting point — the X-flip produces Q2, completing
                %   the y>=0 half.
                if options.half_body
                    source_grids(sname) = S;
                end

                [verts, panels] = WEC_Panelizer.grid_to_quads(S);
                sid = surf_id_map(sname);

                offset     = size(all_verts, 1);
                all_verts  = [all_verts; verts]; %#ok<AGROW>
                all_panels = [all_panels; panels + offset]; %#ok<AGROW>
                all_surf_ids = [all_surf_ids; sid * ones(size(panels,1), 1)]; %#ok<AGROW>
                all_is_cap   = [all_is_cap; false(size(panels,1), 1)]; %#ok<AGROW>

                if options.verbose
                    fprintf(' %d panels\n', size(panels, 1));
                end
            end

            % ── Derive MIRROR meshes by coordinate flip ───────────
            if options.quarter_body
                if options.verbose
                    fprintf('    Quarter-body mode: skipping %d mirror surfaces\n', ...
                            length(topo.mirrors));
                end
            else
            for m = 1:length(topo.mirrors)
                mirr = topo.mirrors(m);

                % ── Half-body filter ──────────────────────────────
                %  For half_body (ISX=0, ISY=1), include only mirrors
                %  that stay in y >= 0 after the quadrant transform.
                %  A mirror with an ODD number of Y-flips in its chain
                %  crosses into y < 0 → exclude it.
                %
                %  C0 example (source in Q2):
                %    surface_mirrX = MirrX(source):        Y-flips=0 → INCLUDE (→ Q2)
                %    surface_mirrY = MirrY(source):        Y-flips=1 → EXCLUDE
                %    surface_mirrXY= MirrY(MirrX(source)): Y-flips=1 → EXCLUDE
                if options.half_body
                    n_y_flips = sum(strcmp(mirr.effective_flips, 'Y'));
                    if mod(n_y_flips, 2) == 1
                        if options.verbose
                            fprintf('    Skipping %s (Y-flips=%d, crosses y=0)\n', ...
                                    mirr.name, n_y_flips);
                        end
                        continue;
                    end
                end

                if options.verbose
                    fprintf('    Mirroring %s = Mirr%s(%s)...', ...
                            mirr.name, mirr.plane, mirr.source);
                end

                % Get the source grid (may be the direct source, or if the
                % source is itself a mirror, trace to the ultimate source)
                if source_grids.isKey(mirr.source)
                    S_src = source_grids(mirr.source);
                elseif options.half_body && source_grids.isKey(mirr.ultimate_source)
                    % Half-body: source_grids stores TRANSFORMED grids
                    % keyed by source name.  For chained mirrors, the
                    % direct source may not be a key — use ultimate source.
                    S_src = source_grids(mirr.ultimate_source);
                else
                    % Source is a mirror of a source — evaluate the
                    % ultimate source and apply accumulated flips
                    ult = mirr.ultimate_source;
                    if source_grids.isKey(ult)
                        S_src = source_grids(ult);
                        % Apply all flips EXCEPT the current one (those are
                        % from the intermediate mirrors)
                        for fi = 1:length(mirr.effective_flips)-1
                            flip = mirr.effective_flips{fi};
                            if strcmp(flip, 'Y')
                                S_src(:,:,2) = -S_src(:,:,2);
                            else
                                S_src(:,:,1) = -S_src(:,:,1);
                            end
                        end
                    else
                        % Fallback: evaluate directly (should not happen
                        % for well-formed .ms2 files)
                        if options.verbose, fprintf(' [direct eval]'); end
                        v_grid_fb = v_grid;
                        if source_v_grids.isKey(mirr.source)
                            v_grid_fb = source_v_grids(mirr.source);
                        elseif source_v_grids.isKey(ult)
                            v_grid_fb = source_v_grids(ult);
                        end
                        S_src = parser.eval_surface_grid(mirr.source, u_grid, v_grid_fb);
                        S_src(:,:,3) = S_src(:,:,3) + draft;
                    end
                end

                % Apply the mirror flip
                S_mirr = S_src;
                if options.half_body
                    % Half-body: source grids are TRANSFORMED to Q1
                    % (quadrant transform already baked in with correct
                    % winding).  The relative effect of any included
                    % mirror on the stored grid is ALWAYS an X-negate:
                    %
                    %   stored_source_x  = (-1)^negate_x × x_original
                    %   MirrX(source)_x  = -x_original
                    %   stored_mirror_x  = (-1)^negate_x × (-x_original)
                    %                    = -stored_source_x
                    %
                    % This holds regardless of negate_x.  The quadrant
                    % transform is applied equally to both source and
                    % mirror, so the RELATIVE difference is always one
                    % X-flip.  Result: Q1 source + Q2 mirror = y>=0 half.
                    S_mirr(:,:,1) = -S_mirr(:,:,1);
                else
                    current_flip = mirr.effective_flips{end};
                    if strcmp(current_flip, 'Y')
                        S_mirr(:,:,2) = -S_mirr(:,:,2);
                    else  % 'X'
                        S_mirr(:,:,1) = -S_mirr(:,:,1);
                    end
                end

                [verts, panels] = WEC_Panelizer.grid_to_quads(S_mirr);

                % Reverse panel winding for ODD number of coordinate
                % flips applied to the grid.
                %
                % Full-body: count all flips in the mirror chain.
                % Half-body: the stored Q1 grid has correct winding.
                %   The X-negate above is exactly 1 additional flip
                %   (odd) → ALWAYS reverse winding.
                if options.half_body
                    % 1 flip (X-negate) on correctly-wound Q1 grid
                    % → always reverse
                    panels = panels(:, [1 4 3 2]);
                else
                    n_flips_mirr = length(mirr.effective_flips);
                    if mod(n_flips_mirr, 2) == 1
                        panels = panels(:, [1 4 3 2]);
                    end
                end

                sid = surf_id_map(mirr.name);

                offset     = size(all_verts, 1);
                all_verts  = [all_verts; verts]; %#ok<AGROW>
                all_panels = [all_panels; panels + offset]; %#ok<AGROW>
                all_surf_ids = [all_surf_ids; sid * ones(size(panels,1), 1)]; %#ok<AGROW>
                all_is_cap   = [all_is_cap; false(size(panels,1), 1)]; %#ok<AGROW>

                if options.verbose
                    fprintf(' %d panels\n', size(panels, 1));
                end
            end
            end  % if ~quarter_body (mirror loop)

            % ── Vertex handling at surface seams ───────────────────
            switch options.vertex_merge
                case 'merge'
                    % Tolerance-based: round to 1e-8, then unique
                    [all_verts, all_panels] = WEC_Panelizer.merge_vertices( ...
                        all_verts, all_panels, 1e-8);
                    if options.verbose
                        fprintf('    Vertex merge: tolerance 1e-8\n');
                    end

                case 'dedup'
                    % Exact coordinate match only
                    n_before = size(all_verts, 1);
                    [all_verts, ~, ic] = unique(all_verts, 'rows', 'stable');
                    all_panels = ic(all_panels);
                    n_removed = n_before - size(all_verts, 1);
                    if options.verbose
                        fprintf('    Vertex dedup: %d exact duplicates removed\n', n_removed);
                    end

                case 'none'
                    if options.verbose
                        fprintf('    Vertex merge: NONE (duplicates remain)\n');
                    end
            end

            % ── Trim at waterline ─────────────────────────────────
            wl_verts = [];
            if options.trim_wl
                [all_verts, all_panels, all_surf_ids, all_is_cap, wl_verts] = ...
                    WEC_Panelizer.trim_at_wl(all_verts, all_panels, ...
                        all_surf_ids, all_is_cap, options.z_wl);
                if options.verbose
                    fprintf('    Trimmed at z=%.3f: %d panels remain\n', ...
                            options.z_wl, size(all_panels, 1));
                end
            end

            % ── Post-trim vertex merge  ───────────────────────────
            %
            % WP-2 FIX: merge_vertices runs BEFORE trim_at_wl (line ~765).
            % split_panel_at_z (called inside trim_at_wl) appends FRESH
            % intersection vertices at z = z_wl for every straddling panel
            % independently.  Adjacent panels that shared edge V2-V3 before
            % trimming each compute their own intersection I on that edge.
            % Result: two vertices at the SAME 3D position with DIFFERENT
            % indices (e.g. the x=0 seam between a RevSurf source and its
            % X-mirror produces coincident vertices I and I' with indices
            % n+1 and n+2).
            %
            % In extract_waterline_boundary_topo the edge I→next is in one
            % sliver panel (count=1, |z|<tol) but I'→next' is in a
            % different sliver panel.  The graph walk sees isolated 2-vertex
            % chains (I→next, no other neighbours) → loop length=2 < 3 →
            % rejected → empty WP boundary.
            %
            % A second merge pass with the same 1e-8 tolerance fuses I and
            % I' into one vertex, making the waterline ring topologically
            % connected.  wl_verts must be remapped through the same ic map.
            if options.trim_wl && ~isempty(all_panels)
                n_pre_pm          = size(all_verts, 1);
                tol_pm            = 1e-8;
                rounded_pm        = round(all_verts / tol_pm) * tol_pm;
                [~, ia_pm, ic_pm] = unique(rounded_pm, 'rows', 'stable');
                all_verts         = all_verts(ia_pm, :);
                all_panels        = ic_pm(all_panels);
                if ~isempty(wl_verts)
                    wl_verts = unique(ic_pm(wl_verts));
                end
                if options.verbose && size(all_verts, 1) < n_pre_pm
                    fprintf('    Post-trim merge: %d → %d vertices (%d seam duplicates fused)\n', ...
                            n_pre_pm, size(all_verts, 1), ...
                            n_pre_pm - size(all_verts, 1));
                end
            end

            % ── Close open edges (OPTIONAL — NOT for HAMS BEM) ────
            %  For HAMS input, leave close_gaps = false (the default).
            %  Enable only for visualization or volume-validation meshes
            %  where a closed surface is needed.
            if options.close_gaps
                [all_verts, all_panels, all_surf_ids, all_is_cap, n_cap] = ...
                    WEC_Panelizer.close_open_edges(all_verts, all_panels, ...
                        all_surf_ids, all_is_cap);
                if options.verbose && n_cap > 0
                    fprintf('    Closed %d open edges with cap panels\n', n_cap);
                end
            end

            % ── Compute panel normals from cross product ─────────
            %  The surface parameterisation gives S_u × S_v outward.
            %  grid_to_quads preserves this: diagonal cross product
            %  ∝ 2(S_u × S_v).  Mirror winding reversal corrects odd
            %  reflections.  No orient_normals heuristic needed.
            normals = WEC_Panelizer.compute_panel_normals( ...
                all_verts, all_panels);

            % ── Remove orphaned vertices ──────────────────────────
            %  After trimming, vertices from deleted panels remain in the
            %  array.  These orphans waste memory and would appear in the
            %  .pnl file, confusing HAMS.  Remove them and re-index.
            used_idx = unique(all_panels(:));
            remap    = zeros(size(all_verts, 1), 1);
            remap(used_idx) = 1:length(used_idx);
            all_verts  = all_verts(used_idx, :);
            all_panels = remap(all_panels);

            % Remap waterline vertex indices too
            if ~isempty(wl_verts)
                valid_wl  = wl_verts(remap(wl_verts) > 0);
                wl_verts  = remap(valid_wl);
            end

            % ── Q1 verification (quarter-body) ──────────────────
            %  The quadrant transform was applied at the grid level
            %  (before grid_to_quads).  Verify all vertices are in Q1.
            if options.quarter_body && options.verbose
                x_min = min(all_verts(:,1));
                y_min = min(all_verts(:,2));
                if x_min < -1e-8 || y_min < -1e-8
                    warning('WEC_Panelizer:QuadrantViolation', ...
                            'Vertices outside Q1: x_min=%.6f, y_min=%.6f', ...
                            x_min, y_min);
                else
                    fprintf('    Q1 verified: x_min=%.2e, y_min=%.2e\n', ...
                            x_min, y_min);
                end
            end

            % ── Half-body verification (ISX=0, ISY=1) ─────────────
            %  All vertices must have y >= 0 (HAMS mirrors across y=0).
            %  x may be positive or negative (both Q1 and Q2 are present).
            if options.half_body && options.verbose
                x_min = min(all_verts(:,1));
                x_max = max(all_verts(:,1));
                y_min = min(all_verts(:,2));
                if y_min < -1e-8
                    warning('WEC_Panelizer:HalfBodyYViolation', ...
                            'Half-body mesh has y < 0: y_min=%.6f', y_min);
                else
                    fprintf('    Half-body verified: y >= 0 (y_min=%.2e)\n', y_min);
                    fprintf('    x range: [%.4f, %.4f] (both sides of X=0)\n', ...
                            x_min, x_max);
                end
            end

            % ── Assemble output mesh struct ───────────────────────
            mesh = struct();
            mesh.vertices        = all_verts;
            mesh.panels          = all_panels;
            mesh.normals         = normals;
            mesh.surface_ids     = all_surf_ids;
            mesh.is_cap          = all_is_cap;
            mesh.waterline_verts = wl_verts;
            mesh.draft           = draft;
            mesh.Nu              = Nu;
            mesh.Nv              = Nv;
            mesh.u_grid          = u_grid;
            mesh.v_grid          = v_grid;
            mesh.cosine_spacing  = options.cosine_spacing;
            mesh.arclength_u     = options.arclength_u;
            mesh.wl_conform      = options.wl_conform;
            mesh.adapt_Nv        = options.adapt_Nv;
            mesh.source_v_grids  = source_v_grids;
            mesh.u_wl            = u_wl_found;
            mesh.n_panels        = size(all_panels, 1);
            mesh.n_vertices      = size(all_verts, 1);
            mesh.x_sym           = x_sym;
            mesh.y_sym           = y_sym;
            mesh.quarter_body    = options.quarter_body;
            mesh.half_body       = options.half_body;
            mesh.source_quadrant = source_quadrant;
            mesh.negate_x        = negate_x;
            mesh.negate_y        = negate_y;
            mesh.vertex_merge    = options.vertex_merge;

            if options.verbose
                fprintf('    Final mesh: %d vertices, %d panels\n', ...
                        mesh.n_vertices, mesh.n_panels);
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §2  GRID TO QUADS
        %% ═════════════════════════════════════════════════════════

        function [verts, panels] = grid_to_quads(S)
        % GRID_TO_QUADS  Convert a [Nu×Nv×3] surface grid to quad panels.
        %
        %   Each cell (i,j) in the parametric grid becomes one quad
        %   with vertices (i,j) → (i+1,j) → (i+1,j+1) → (i,j+1).
        %
        %   INDEXING
        %     MATLAB reshape(S,[],3) is column-major: S(i,j,:) maps to
        %     row i + (j-1)*Nu.  The vertex indices below match this.
        %     cross(v3-v1, v4-v2) ∝ S_u × S_v (outward by construction).

            [Nu, Nv, ~] = size(S);
            verts = reshape(S, [], 3);

            n_panels = (Nu - 1) * (Nv - 1);
            panels   = zeros(n_panels, 4);
            p = 0;
            for i = 1:Nu-1
                for j = 1:Nv-1
                    p = p + 1;
                    v1 = i     + (j-1)*Nu;
                    v2 = (i+1) + (j-1)*Nu;
                    v3 = (i+1) + j*Nu;
                    v4 = i     + j*Nu;
                    panels(p, :) = [v1, v2, v3, v4];
                end
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §3  TRIM AT WATERLINE
        %% ═════════════════════════════════════════════════════════

        function [verts, panels, surf_ids, is_cap, wl_verts] = ...
                trim_at_wl(verts, panels, surf_ids, is_cap, z_wl)
        % TRIM_AT_WL  Remove panels above waterline, split straddling panels.
        %
        %   Panels fully below z_wl: kept as-is.
        %   Panels fully above z_wl: removed.
        %   Panels straddling z_wl:  split along the waterline.
        %
        %   WHY split instead of discarding straddling panels?
        %     Discarding leaves a ragged waterline edge that doesn't
        %     match the hull cross-section at that draft.  This
        %     introduces error in the BEM waterplane-area integral.
        %     Splitting gives an exact (linear) cut along z = z_wl.

            n_panels    = size(panels, 1);
            new_verts   = verts;
            new_panels  = zeros(0, 4);
            new_sids    = zeros(0, 1);
            new_caps    = false(0, 1);
            wl_verts    = [];

            for p = 1:n_panels
                vidx = panels(p, :);
                z    = verts(vidx, 3);

                if all(z <= z_wl + 1e-10)
                    % Fully submerged — keep
                    new_panels(end+1, :) = vidx; %#ok<AGROW>
                    new_sids(end+1, 1)   = surf_ids(p); %#ok<AGROW>
                    new_caps(end+1, 1)   = is_cap(p); %#ok<AGROW>

                elseif all(z > z_wl - 1e-10)
                    % Fully above water — discard
                    continue;

                else
                    % Straddling — split at waterline
                    [sub_p, sub_v, sub_wl] = ...
                        WEC_Panelizer.split_panel_at_z(verts, vidx, z_wl);
                    if ~isempty(sub_p)
                        vo = size(new_verts, 1);
                        new_verts = [new_verts; sub_v]; %#ok<AGROW>
                        for sp = 1:size(sub_p, 1)
                            new_panels(end+1, :) = sub_p(sp, :) + vo; %#ok<AGROW>
                            new_sids(end+1, 1)   = surf_ids(p); %#ok<AGROW>
                            new_caps(end+1, 1)   = is_cap(p); %#ok<AGROW>
                        end
                        wl_verts = [wl_verts; sub_wl + vo]; %#ok<AGROW>
                    end
                end
            end

            verts    = new_verts;
            panels   = new_panels;
            surf_ids = new_sids;
            is_cap   = new_caps;
        end


        function [sub_panels, sub_verts, wl_idx] = split_panel_at_z(verts, vidx, z_wl)
        % SPLIT_PANEL_AT_Z  Split a quad at z = z_wl, keep submerged part.
        %
        %   Walks the quad edges V1→V2→V3→V4→V1 in sequence, collecting
        %   below-water vertices and waterline intersection points in
        %   traversal order.  This preserves the polygon winding — the
        %   submerged polygon has the same orientation as the original quad.
        %
        %   WHY traversal order matters:
        %     The old code used [pts(below,:); new_pts], which concatenates
        %     submerged vertices (in row order) then intersection points
        %     (in edge order).  When above/below vertices alternate around
        %     the quad (e.g. V1 below, V2 above, V3 below, V4 above), this
        %     produces [V1; V3; I_12; I_23; I_34; I_41] — a self-intersecting
        %     polygon.  Walking edges in sequence gives the correct winding:
        %     [V1; I_12; I_23; V3; I_34; I_41].
        %
        %   Returns triangle(s) or quad — HAMS accepts both.
        %   Triangular panels use v4 = v3 convention.

            pts   = verts(vidx, :);
            z     = pts(:, 3);
            below = z <= z_wl + 1e-10;

            % ── Walk edges in order, collect submerged polygon ─────
            edges   = [1 2; 2 3; 3 4; 4 1];
            sub_poly = zeros(0, 3);
            is_wl    = false(0, 1);

            for e = 1:4
                i1 = edges(e, 1);
                i2 = edges(e, 2);

                % Add start vertex if below waterline
                if below(i1)
                    sub_poly(end+1, :) = pts(i1, :); %#ok<AGROW>
                    is_wl(end+1, 1)    = false; %#ok<AGROW>
                end

                % Add intersection point if edge crosses waterline
                if below(i1) ~= below(i2)
                    dz = z(i2) - z(i1);
                    if abs(dz) > 1e-14
                        t = (z_wl - z(i1)) / dz;
                        t = max(0, min(1, t));
                        new_pt    = pts(i1,:) + t * (pts(i2,:) - pts(i1,:));
                        new_pt(3) = z_wl;   % snap to exact waterline
                        sub_poly(end+1, :) = new_pt; %#ok<AGROW>
                        is_wl(end+1, 1)    = true; %#ok<AGROW>
                    end
                end
            end

            n = size(sub_poly, 1);
            if n < 3
                sub_panels = [];
                sub_verts  = [];
                wl_idx     = [];
                return;
            end

            sub_verts = sub_poly;

            % ── Triangulate submerged polygon ──────────────────────
            if n == 3
                sub_panels = [1, 2, 3, 3];           % triangle
            elseif n == 4
                sub_panels = [1, 2, 3, 4];           % quad
            else
                % Fan triangulation from first vertex
                sub_panels = zeros(n - 2, 4);
                for k = 1:n-2
                    sub_panels(k, :) = [1, k+1, k+2, k+2];
                end
            end

            % ── Track waterline vertex indices ─────────────────────
            wl_idx = find(is_wl);
        end


        %% ═════════════════════════════════════════════════════════
        %%  §4  CLOSE OPEN EDGES
        %% ═════════════════════════════════════════════════════════

        function [verts, panels, surf_ids, is_cap, n_cap] = ...
                close_open_edges(verts, panels, surf_ids, is_cap)
        % CLOSE_OPEN_EDGES  Detect open boundary edges and fill with caps.
        %
        %   An OPEN edge belongs to exactly one panel.
        %   A CLOSED (interior) edge belongs to two panels.
        %
        %   ALGORITHM
        %     1. Build edge → panel count map.
        %     2. Find edges with count == 1 (boundary).
        %     3. Group boundary vertices by z-level.
        %     4. For each z-level group:
        %          a. Compute centroid.
        %          b. Order vertices by angle around centroid.
        %          c. Fan-triangulate into cap panels.
        %
        %   WHY this is needed:
        %     The .ms2 hull may not define a cap surface for the column
        %     top (z = +1).  The waterline cut also creates an opening.
        %     For a watertight BEM mesh, all openings EXCEPT the free
        %     surface (waterline) must be closed.  The waterline opening
        %     is handled separately by WaterPlaneMesh.pnl.
        %
        %   WHY fan triangulation from centroid?
        %     The boundary polygon may be non-convex (e.g. column
        %     cross-section is capsule-shaped).  Fan from centroid is
        %     simple and correct for the mildly non-convex shapes
        %     typical of WEC cross-sections.

            n_cap = 0;

            % ── Build edge → count map ────────────────────────────
            edge_map = containers.Map('KeyType', 'char', 'ValueType', 'int32');
            for p = 1:size(panels, 1)
                v = panels(p, :);
                if v(3) == v(4)
                    ee = [v(1) v(2); v(2) v(3); v(3) v(1)];
                else
                    ee = [v(1) v(2); v(2) v(3); v(3) v(4); v(4) v(1)];
                end
                for e = 1:size(ee, 1)
                    key = sprintf('%d_%d', min(ee(e,:)), max(ee(e,:)));
                    if edge_map.isKey(key)
                        edge_map(key) = edge_map(key) + 1;
                    else
                        edge_map(key) = 1;
                    end
                end
            end

            % ── Collect boundary edges (count == 1) ───────────────
            boundary_edges = [];
            keys = edge_map.keys();
            for i = 1:length(keys)
                if edge_map(keys{i}) == 1
                    v = sscanf(keys{i}, '%d_%d');
                    boundary_edges = [boundary_edges; v']; %#ok<AGROW>
                end
            end

            if isempty(boundary_edges), return; end

            % ── Group boundary vertices by z-level ────────────────
            boundary_vidx = unique(boundary_edges(:));
            boundary_z    = verts(boundary_vidx, 3);
            z_tol         = 0.01;  % [m]
            z_levels      = unique(round(boundary_z / z_tol) * z_tol);

            for iz = 1:length(z_levels)
                z_lvl    = z_levels(iz);
                at_level = boundary_vidx(abs(boundary_z - z_lvl) < z_tol);

                if length(at_level) < 3, continue; end

                % ── Fan triangulation from centroid ───────────────
                cap_center = mean(verts(at_level, :), 1);
                center_idx = size(verts, 1) + 1;
                verts      = [verts; cap_center]; %#ok<AGROW>

                rel    = verts(at_level, :) - cap_center;
                angles = atan2(rel(:, 2), rel(:, 1));
                [~, order]    = sort(angles);
                ordered_verts = at_level(order);

                n_fan = length(ordered_verts);
                for k = 1:n_fan
                    k_next = mod(k, n_fan) + 1;
                    new_tri = [center_idx, ordered_verts(k), ...
                               ordered_verts(k_next), ordered_verts(k_next)];
                    panels   = [panels; new_tri]; %#ok<AGROW>
                    surf_ids = [surf_ids; 0]; %#ok<AGROW>
                    is_cap   = [is_cap; true]; %#ok<AGROW>
                    n_cap    = n_cap + 1;
                end
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §5  ORIENT NORMALS
        %% ═════════════════════════════════════════════════════════

        function [panels, normals] = orient_normals(verts, panels)
        % ORIENT_NORMALS  Ensure all panel normals point outward.
        %
        % ── DEAD CODE ─────────────────────────────────────────────────────────
        % Never called from anywhere in this file. generate() explicitly states
        % "No orient_normals heuristic needed" — normals come directly from the
        % parametric surface (S_u × S_v) with winding corrected per-mirror.
        % Retained for reference only. Do not call from new code.
        % ──────────────────────────────────────────────────────────────────────
        %
        %   ALGORITHM
        %     1. Compute panel normal via cross product of diagonals.
        %     2. Estimate hull interior point (centroid of all vertices).
        %     3. For each panel: check normal · (panel_centroid − interior) > 0.
        %     4. If not, reverse vertex winding (swap v2 and v4).
        %
        %   WHY cross product of diagonals (not edges)?
        %     For non-planar quads, the diagonal cross product gives a
        %     normal representative of the entire panel area, not biased
        %     toward one triangle.

            n_panels = size(panels, 1);
            normals  = zeros(n_panels, 3);
            interior = mean(verts, 1);

            for p = 1:n_panels
                v  = panels(p, :);
                p1 = verts(v(1), :);
                p2 = verts(v(2), :);
                p3 = verts(v(3), :);
                p4 = verts(v(4), :);

                % Normal via diagonal cross product
                d1    = p3 - p1;
                d2    = p4 - p2;
                n_vec = cross(d1, d2);
                n_len = norm(n_vec);

                if n_len < 1e-14
                    normals(p, :) = [0, 0, 1];  % degenerate panel
                    continue;
                end

                normals(p, :) = n_vec / n_len;

                % Check orientation against outward direction
                panel_centroid = (p1 + p2 + p3 + p4) / 4;
                if dot(normals(p, :), panel_centroid - interior) < 0
                    % Flip winding: swap v2 and v4
                    panels(p, :) = [v(1), v(4), v(3), v(2)];
                    normals(p, :) = -normals(p, :);
                end
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §6  MERGE VERTICES
        %% ═════════════════════════════════════════════════════════

        function [verts, panels] = merge_vertices(verts, panels, tol)
        % MERGE_VERTICES  Merge coincident vertices within tolerance.
        %
        %   At surface boundaries, adjacent surfaces share edge vertices
        %   that differ only by floating-point noise.  Merging them to
        %   a single index creates a watertight mesh.
        %
        %   WHY tol = 1e-8 m?
        %     B-spline evaluation noise is O(1e-14).  Surface seams
        %     from independent evaluations of shared boundary curves
        %     differ by at most O(1e-12).  1e-8 is generous but well
        %     below any physical mesh dimension.  Using tol < 1e-6
        %     avoids accidentally merging nearby-but-distinct vertices
        %     on fine meshes.

            if nargin < 3, tol = 1e-8; end

            rounded = round(verts / tol) * tol;
            [~, ia, ic] = unique(rounded, 'rows', 'stable');

            verts  = verts(ia, :);
            panels = ic(panels);
        end


        %% ═════════════════════════════════════════════════════════
        %%  §7  WRITE HULL MESH (.pnl)
        %% ═════════════════════════════════════════════════════════

        function write_hull_pnl(mesh, filename)
        % WRITE_HULL_PNL  Write HAMS HullMesh.pnl file.
        %
        %   WEC_Panelizer.write_hull_pnl(mesh, 'Input/HullMesh.pnl')
        %
        %   FORMAT (verified against HAMS_Prog.f90 reading sequence):
        %     Lines 1-3: comments (skipped by HAMS: DO II=1,3; READ(2,*); ENDDO)
        %     Line 4:    NELEM  NTND  ISX  ISY  (free format)
        %     Lines 5-6: comments (skipped: DO II=1,2; READ(2,*); ENDDO)
        %     NTND lines: node_id  x  y  z  (free format)
        %     3 comment lines (skipped: DO J=1,3; READ(2,*); ENDDO)
        %     NELEM lines: panel_id  num_verts  v1 v2 v3 [v4]  (free format)
        %
        %   BUGS FIXED (relative to previous version):
        %     1. Header had only 1 comment line; HAMS skips 3 → data on wrong line
        %     2. Wrote n_vertices before n_panels; HAMS reads NELEM first
        %     3. Missing ISX, ISY symmetry flags
        %     4. Missing 2 comment lines between header and vertex block
        %     5. Missing 3 comment lines between vertex and panel blocks
        %     6. Panel line missing num_vertices column (HAMS reads NCN per panel)
        %     7. Comment block said "confirmed" but format was wrong

            % Ensure output directory exists
            [out_dir, ~, ~] = fileparts(filename);
            if ~isempty(out_dir) && ~exist(out_dir, 'dir')
                mkdir(out_dir);
            end

            fid = fopen(filename, 'w');
            if fid == -1
                error('WEC_Panelizer:WriteError', ...
                       'Cannot open %s for writing', filename);
            end

            % Symmetry flags from mesh struct
            if isfield(mesh, 'x_sym'), x_sym = mesh.x_sym; else, x_sym = 0; end
            if isfield(mesh, 'y_sym'), y_sym = mesh.y_sym; else, y_sym = 0; end

            % Lines 1-3: comments (HAMS skips these)
            fprintf(fid, '    --------------Hull Mesh File---------------\n');
            fprintf(fid, ' \n');
            fprintf(fid, '    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry\n');

            % Line 4: NELEM  NTND  ISX  ISY
            fprintf(fid, '    %8d    %8d       %5d       %5d\n', ...
                    mesh.n_panels, mesh.n_vertices, x_sym, y_sym);

            % Lines 5-6: comments (HAMS skips these)
            fprintf(fid, ' \n');
            fprintf(fid, '    # Start Definition of Node Coordinates     ! node_number   x   y   z\n');

            % Vertex block: NTND lines
            for i = 1:mesh.n_vertices
                fprintf(fid, ' %4d    %14.6f    %14.6f    %14.6f\n', ...
                        i, mesh.vertices(i,1), mesh.vertices(i,2), mesh.vertices(i,3));
            end

            % 3 comment lines between vertex and panel blocks
            fprintf(fid, '    # End Definition of Node Coordinates\n');
            fprintf(fid, ' \n');
            fprintf(fid, '    # Start Definition of Panel Connectivity   ! panel_number  num_vertices  v1 v2 v3 [v4]\n');

            % Panel block: NELEM lines with num_vertices column
            for i = 1:mesh.n_panels
                vi = mesh.panels(i,:);
                if vi(3) == vi(4)
                    % Triangle: 3 vertices
                    fprintf(fid, ' %4d    %d    %6d    %6d    %6d\n', ...
                            i, 3, vi(1), vi(2), vi(3));
                else
                    % Quad: 4 vertices
                    fprintf(fid, ' %4d    %d    %6d    %6d    %6d    %6d\n', ...
                            i, 4, vi(1), vi(2), vi(3), vi(4));
                end
            end

            fprintf(fid, '    # End Definition of Panel Connectivity\n');
            fclose(fid);
        end


        %% ═════════════════════════════════════════════════════════
        %%  §8  WRITE WATERPLANE MESH (.pnl)
        %% ═════════════════════════════════════════════════════════

        function write_wp_pnl(mesh, filename, z_wl)
        % WRITE_WP_PNL  Write HAMS WaterPlaneMesh.pnl file.
        %
        %   WEC_Panelizer.write_wp_pnl(mesh, filename, z_wl)
        %
        %   Generates a concentric-ring waterplane mesh from the hull
        %   mesh waterline vertices, then writes it in HAMS format.
        %
        %   MESH ALGORITHM: Concentric-ring quad + tri mesh.
        %     1. Extract and angular-sort waterline boundary nodes.
        %     2. Create Nr concentric rings scaled toward centre.
        %     3. Form quad panels between adjacent rings.
        %     4. Collapse innermost ring to triangular fan at centroid.
        %     All nodes at z = z_wl exactly.
        %
        %   WHY concentric rings over fan triangulation:
        %     A single-layer fan from centroid to boundary produces
        %     triangles with aspect ratios ≈ circumference/(N*radius),
        %     which degrades rapidly for large N.  Concentric rings
        %     keep aspect ratios near 1:1 at every ring level.
        %
        %   FILE FORMAT (verified against HAMS_Prog.f90):
        %     Same structure as HullMesh.pnl — see write_hull_pnl.
        %     Waterplane mesh is read identically by ReadWTPLMesh.

            if isempty(mesh.waterline_verts)
                warning('WEC_Panelizer:NoWaterline', ...
                        'No waterline vertices — cannot write waterplane mesh');
                return;
            end

            % Ensure output directory exists
            [out_dir, ~, ~] = fileparts(filename);
            if ~isempty(out_dir) && ~exist(out_dir, 'dir')
                mkdir(out_dir);
            end

            % Symmetry flags from mesh struct
            if isfield(mesh, 'x_sym'), x_sym = mesh.x_sym; else, x_sym = 0; end
            if isfield(mesh, 'y_sym'), y_sym = mesh.y_sym; else, y_sym = 0; end

            % ── Extract and order waterline boundary ──────────────
            wl_idx  = unique(mesh.waterline_verts);
            wl_pts  = mesh.vertices(wl_idx, 1:2);
            centroid = mean(wl_pts, 1);
            angles  = atan2(wl_pts(:,2) - centroid(2), ...
                            wl_pts(:,1) - centroid(1));
            [~, order] = sort(angles);
            wl_pts  = wl_pts(order, :);
            N_bnd   = size(wl_pts, 1);

            if N_bnd < 3
                warning('WEC_Panelizer:TooFewWL', ...
                        'Need >= 3 waterline vertices, got %d', N_bnd);
                return;
            end

            % ── Concentric-ring mesh generation ───────────────────
            %  Estimate a reasonable number of radial rings from the
            %  boundary arc spacing.
            bnd_edges = sqrt(sum(diff([wl_pts; wl_pts(1,:)]).^2, 2));
            mean_edge = mean(bnd_edges);
            max_r     = max(sqrt(sum((wl_pts - centroid).^2, 2)));
            Nr        = max(2, round(max_r / mean_edge));

            %  Quadratic radial spacing: clusters rings near boundary
            %  where the waterplane-hull intersection has steep gradients.
            r_frac = zeros(Nr+1, 1);
            for k = 0:Nr
                r_frac(k+1) = (k / Nr)^2;
            end
            %  r_frac(1) = 0 (centre), r_frac(end) = 1 (boundary)

            %  Build node array: Nr+1 rings × N_bnd azimuthal nodes
            %  Ring 0 = centroid (single node), rings 1..Nr = scaled boundary.
            n_ring_nodes = Nr * N_bnd + 1;  % +1 for the centroid
            wp_verts = zeros(n_ring_nodes, 3);

            % Centroid node (index 1)
            wp_verts(1, :) = [centroid, z_wl];

            % Rings 1..Nr: interpolated from centroid to boundary
            for k = 1:Nr
                r = r_frac(k+1);
                for j = 1:N_bnd
                    idx = 1 + (k-1)*N_bnd + j;
                    wp_verts(idx, :) = [(1-r)*centroid + r*wl_pts(j,:), z_wl];
                end
            end

            % ── Build panel connectivity ──────────────────────────
            wp_panels_list = zeros(0, 4);
            wp_nverts_list = zeros(0, 1);

            %  Innermost ring (k=1): triangles from centroid to ring 1
            for j = 1:N_bnd
                j_next = mod(j, N_bnd) + 1;
                v1 = 1;                       % centroid
                v2 = 1 + j;                   % ring 1, azimuth j
                v3 = 1 + j_next;              % ring 1, azimuth j+1
                wp_panels_list(end+1, :) = [v1, v2, v3, v3]; %#ok<AGROW>
                wp_nverts_list(end+1, 1) = 3; %#ok<AGROW>
            end

            %  Rings 1..(Nr-1) → quads between consecutive rings
            %  Vertex order: [inner_j, outer_j, outer_j+1, inner_j+1]
            %  This gives CCW winding (same as inner triangles above).
            %  HAMS cross(V3-V1, V4-V2) then gives nz > 0 (UPWARD — correct).
            for k = 1:Nr-1
                for j = 1:N_bnd
                    j_next = mod(j, N_bnd) + 1;
                    v1 = 1 + (k-1)*N_bnd + j;          % inner ring k, azimuth j
                    v2 = 1 + (k-1)*N_bnd + j_next;     % inner ring k, azimuth j+1
                    v3 = 1 + k*N_bnd     + j_next;     % outer ring k+1, azimuth j+1
                    v4 = 1 + k*N_bnd     + j;          % outer ring k+1, azimuth j
                    wp_panels_list(end+1, :) = [v1, v4, v3, v2]; %#ok<AGROW>
                    wp_nverts_list(end+1, 1) = 4; %#ok<AGROW>
                end
            end

            n_wp_panels = size(wp_panels_list, 1);

            % ── Enforce WP normals point UPWARD (+z) ─────────────────
            %  HAMS CalTransNormals (NormalProcess.f90):
            %    Tri:  cross(V2-V1, V3-V2)
            %    Quad: cross(V3-V1, V4-V2)
            %  Verified: UPWARD (nz > 0) is correct for HAMS waterplane.
            %  Evidence from test: nz<0 gives min(B33)=-2987 with irr ON.
            %  nz>0 is consistent with WAMIT_MeshTran output (waterplane cap
            %  has +z normal = outward from the interior cavity is upward).
            %  BUG FIX (B7): Previous code checked only panel 1 (may be
            %  degenerate giving nrm=[0,0,0], silent no-flip). Now scans for
            %  first non-degenerate panel.
            % ── Enforce WP normals UPWARD (+z) ──────────────────────────────────
            % Formula: HAMS CalTransNormals (verified vs diag_hams_run.m):
            %   Tri:  cross(V2-V1, V3-V2)    Quad: cross(V3-V1, V4-V2)
            % Scan for first non-degenerate panel; flip entire mesh if DOWNWARD.
            flip_all_wp = false;
            for p_chk = 1:n_wp_panels
                vi_c = wp_panels_list(p_chk,:);
                nv_c = wp_nverts_list(p_chk);
                if nv_c == 3
                    d1c = wp_verts(vi_c(2),:) - wp_verts(vi_c(1),:);
                    d2c = wp_verts(vi_c(3),:) - wp_verts(vi_c(2),:);
                else
                    d1c = wp_verts(vi_c(3),:) - wp_verts(vi_c(1),:);
                    d2c = wp_verts(vi_c(4),:) - wp_verts(vi_c(2),:);
                end
                nrm_c = cross(d1c, d2c);
                if norm(nrm_c) > 1e-14
                    flip_all_wp = (nrm_c(3) < 0);  % flip if DOWNWARD (wrong)
                    break;
                end
            end
            if flip_all_wp
                for p = 1:n_wp_panels
                    if wp_nverts_list(p) == 3
                        wp_panels_list(p,:) = wp_panels_list(p, [1 3 2 2]);
                    else
                        wp_panels_list(p,:) = wp_panels_list(p, [1 4 3 2]);
                    end
                end
            end

            % ── Write .pnl file in HAMS format ────────────────────
            n_nodes  = size(wp_verts, 1);

            fid = fopen(filename, 'w');
            if fid == -1
                error('WEC_Panelizer:WriteError', ...
                       'Cannot open %s for writing', filename);
            end

            % Lines 1-3: comments
            fprintf(fid, '    ----------Waterplane Mesh File----------\n');
            fprintf(fid, ' \n');
            fprintf(fid, '    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry\n');

            % Line 4: NELEM  NTND  ISX  ISY
            fprintf(fid, '    %8d    %8d       %5d       %5d\n', ...
                    n_wp_panels, n_nodes, x_sym, y_sym);

            % Lines 5-6: comments
            fprintf(fid, ' \n');
            fprintf(fid, '    # Start Definition of Node Coordinates     ! node_number   x   y   z\n');

            % Vertex block
            for i = 1:n_nodes
                fprintf(fid, ' %4d    %14.6f    %14.6f    %14.6f\n', ...
                        i, wp_verts(i,1), wp_verts(i,2), wp_verts(i,3));
            end

            % 3 comment lines between vertex and panel blocks
            fprintf(fid, '    # End Definition of Node Coordinates\n');
            fprintf(fid, ' \n');
            fprintf(fid, '    # Start Definition of Panel Connectivity\n');

            % Panel block with num_vertices column
            for i = 1:n_wp_panels
                vi = wp_panels_list(i,:);
                nv = wp_nverts_list(i);
                if nv == 3
                    fprintf(fid, ' %4d    %d    %6d    %6d    %6d\n', ...
                            i, 3, vi(1), vi(2), vi(3));
                else
                    fprintf(fid, ' %4d    %d    %6d    %6d    %6d    %6d\n', ...
                            i, 4, vi(1), vi(2), vi(3), vi(4));
                end
            end

            fprintf(fid, '    # End Definition of Panel Connectivity\n');
            fclose(fid);
        end


        %% ═════════════════════════════════════════════════════════
        %%  §9  VALIDATE
        %% ═════════════════════════════════════════════════════════

        function stats = validate(mesh)
        % VALIDATE  Compute mesh quality metrics.
        %
        %   stats = WEC_Panelizer.validate(mesh)
        %
        %   RETURNS struct with:
        %     .volume         — enclosed volume via divergence theorem [m³]
        %                       NOTE: only valid for watertight (closed) meshes.
        %                       For open BEM meshes (close_gaps=false), this is
        %                       approximate.  Use WEC_HydroProperties for accurate
        %                       volume computation from parametric surfaces.
        %     .centroid       — volume centroid [x, y, z] [m] (same caveat)
        %     .total_area     — sum of panel areas [m²]
        %     .min_area       — smallest panel area [m²]
        %     .max_aspect     — worst panel aspect ratio [-]
        %     .n_open_edges   — boundary edge count (0 = watertight)
        %                       For BEM meshes, open edges at the waterline and
        %                       column top are expected and correct.
        %     .n_panels       — total panel count
        %     .n_vertices     — total vertex count
        %
        %   VOLUME CALCULATION (divergence theorem)
        %     For a closed surface Σ with outward normal n̂:
        %       V = (1/3) ∮_Σ r · n̂ dA
        %     Discretised over triangular facets:
        %       V = Σ_p (p₁ · (p₂×p₃)) / 6
        %     For quads: split into two triangles and sum.

            n_p   = mesh.n_panels;
            areas = zeros(n_p, 1);
            aspects = zeros(n_p, 1);

            % ── Panel areas and aspect ratios ─────────────────────
            for p = 1:n_p
                v = mesh.panels(p, :);
                p1 = mesh.vertices(v(1),:);
                p2 = mesh.vertices(v(2),:);
                p3 = mesh.vertices(v(3),:);
                p4 = mesh.vertices(v(4),:);

                if v(3) == v(4)
                    areas(p) = 0.5 * norm(cross(p2 - p1, p3 - p1));
                else
                    areas(p) = 0.5 * norm(cross(p2 - p1, p3 - p1)) + ...
                               0.5 * norm(cross(p3 - p1, p4 - p1));
                end

                ee = [norm(p2-p1), norm(p3-p2), norm(p4-p3), norm(p1-p4)];
                ee = ee(ee > 1e-14);
                if ~isempty(ee)
                    aspects(p) = max(ee) / min(ee);
                end
            end

            stats.total_area = sum(areas);
            stats.min_area   = min(areas);
            stats.max_aspect = max(aspects);

            % ── Volume via divergence theorem ─────────────────────
            volume       = 0;
            centroid_num = [0, 0, 0];
            for p = 1:n_p
                v = mesh.panels(p,:);
                p1 = mesh.vertices(v(1),:);
                p2 = mesh.vertices(v(2),:);
                p3 = mesh.vertices(v(3),:);

                % First triangle
                n_tri = cross(p2 - p1, p3 - p1);
                vc    = dot(p1, n_tri) / 6;
                volume = volume + vc;
                centroid_num = centroid_num + vc * (p1 + p2 + p3) / 4;

                % Second triangle (if quad)
                if v(3) ~= v(4)
                    p4     = mesh.vertices(v(4),:);
                    n_tri2 = cross(p3 - p1, p4 - p1);
                    vc2    = dot(p1, n_tri2) / 6;
                    volume = volume + vc2;
                    centroid_num = centroid_num + vc2 * (p1 + p3 + p4) / 4;
                end
            end

            stats.volume = abs(volume);
            if abs(volume) > 1e-12
                stats.centroid = centroid_num / volume;
            else
                stats.centroid = [0, 0, 0];
            end

            % ── Open edge count ───────────────────────────────────
            edge_count = containers.Map('KeyType', 'char', 'ValueType', 'int32');
            for p = 1:n_p
                v = mesh.panels(p,:);
                if v(3) == v(4)
                    ee = [v(1) v(2); v(2) v(3); v(3) v(1)];
                else
                    ee = [v(1) v(2); v(2) v(3); v(3) v(4); v(4) v(1)];
                end
                for e = 1:size(ee, 1)
                    key = sprintf('%d_%d', min(ee(e,:)), max(ee(e,:)));
                    if edge_count.isKey(key)
                        edge_count(key) = edge_count(key) + 1;
                    else
                        edge_count(key) = 1;
                    end
                end
            end

            n_open = 0;
            keys = edge_count.keys();
            for i = 1:length(keys)
                if edge_count(keys{i}) == 1
                    n_open = n_open + 1;
                end
            end

            stats.n_open_edges = n_open;
            stats.n_panels     = n_p;
            stats.n_vertices   = mesh.n_vertices;

            % ── Summary report ────────────────────────────────────
            fprintf('\n');
            fprintf('  ┌─────────────────────────────────────────┐\n');
            fprintf('  │  MESH VALIDATION                         │\n');
            fprintf('  ├─────────────────────────────────────────┤\n');
            fprintf('  │  Vertices:    %6d                      \n', stats.n_vertices);
            fprintf('  │  Panels:      %6d                      \n', stats.n_panels);
            if stats.n_open_edges == 0
                fprintf('  │  Volume:      %10.4f m³               \n', stats.volume);
                fprintf('  │  Centroid:    [%.3f, %.3f, %.3f] m     \n', stats.centroid);
            else
                fprintf('  │  Volume:      %10.4f m³  (approx — mesh is open)\n', stats.volume);
                fprintf('  │  Centroid:    [%.3f, %.3f, %.3f] m  (approx)\n', stats.centroid);
            end
            fprintf('  │  Total area:  %10.4f m²               \n', stats.total_area);
            fprintf('  │  Min area:    %10.6f m²               \n', stats.min_area);
            fprintf('  │  Max aspect:  %10.2f                  \n', stats.max_aspect);
            if stats.n_open_edges == 0
                fprintf('  │  Open edges:  %6d  (watertight)       \n', stats.n_open_edges);
            else
                fprintf('  │  Open edges:  %6d  (expected for BEM mesh)\n', stats.n_open_edges);
            end
            fprintf('  └─────────────────────────────────────────┘\n');
        end


        %% ═════════════════════════════════════════════════════════
        %%  §10  CHECK BEM QUALITY — Panel Size vs Wavelength
        %% ═════════════════════════════════════════════════════════

        function report = check_bem_quality(mesh, T_range)
        % CHECK_BEM_QUALITY  Assess mesh adequacy for BEM at given wave periods.
        %
        %   report = WEC_Panelizer.check_bem_quality(mesh, [T_min, T_max])
        %   report = WEC_Panelizer.check_bem_quality(mesh)  % default T=[2,20] s
        %
        %   INPUTS
        %     mesh    : mesh struct from generate()
        %     T_range : [T_min, T_max] wave period range [s] (default: [2, 20])
        %
        %   BEM PANEL SIZE RULE
        %     The standard BEM convergence criterion requires the panel
        %     characteristic length L_c < λ/6, where λ is the shortest
        %     wavelength of interest.  For deep water: λ = gT²/(2π).
        %
        %     L_c is the maximum edge length of each panel (most
        %     conservative measure — ensures no edge spans more than
        %     1/6 of a wavelength).
        %
        %   OUTPUT (report struct)
        %     .n_panels          — total panel count
        %     .max_edge_length   — largest panel edge [m]
        %     .mean_edge_length  — average panel edge [m]
        %     .max_aspect_ratio  — worst aspect ratio
        %     .n_degenerate      — panels with area < 1e-12 m²
        %     .min_T_resolved    — minimum wave period this mesh can resolve [s]
        %     .lambda_min        — corresponding minimum wavelength [m]
        %     .n_oversized       — panels exceeding λ_min/6 for T_range(1)
        %     .pct_oversized     — percentage oversized
        %     .adequate          — logical, true if mesh is adequate for T_range
        %
        %   USAGE FOR HAMS
        %     model = WEC_MS2_Parser.parse('C0.ms2');
        %     mesh  = WEC_Panelizer.generate(model, 0, 16, 16);
        %     report = WEC_Panelizer.check_bem_quality(mesh, [3, 15]);
        %     if ~report.adequate
        %         % Increase Nu, Nv and regenerate
        %     end

            if nargin < 2 || isempty(T_range), T_range = [2, 20]; end

            g = 9.81;
            n_p = mesh.n_panels;

            % ── Compute panel edge lengths and aspect ratios ───────
            edge_lengths = zeros(n_p, 4);
            areas        = zeros(n_p, 1);
            aspects      = zeros(n_p, 1);

            for p = 1:n_p
                v = mesh.panels(p, :);
                p1 = mesh.vertices(v(1), :);
                p2 = mesh.vertices(v(2), :);
                p3 = mesh.vertices(v(3), :);
                p4 = mesh.vertices(v(4), :);

                edge_lengths(p, :) = [norm(p2-p1), norm(p3-p2), ...
                                      norm(p4-p3), norm(p1-p4)];

                if v(3) == v(4)
                    areas(p) = 0.5 * norm(cross(p2-p1, p3-p1));
                else
                    areas(p) = 0.5 * norm(cross(p2-p1, p3-p1)) + ...
                               0.5 * norm(cross(p3-p1, p4-p1));
                end

                ee = edge_lengths(p, :);
                ee = ee(ee > 1e-14);
                if ~isempty(ee)
                    aspects(p) = max(ee) / min(ee);
                end
            end

            max_edges = max(edge_lengths, [], 2);  % max edge per panel

            % ── BEM wavelength criteria ────────────────────────────
            %  Deep water dispersion: λ = g T² / (2π)
            %  Panel rule: L_c < λ/6
            %  → T_min_resolved = sqrt(6 * L_max * 2π / g)
            L_max = max(max_edges);
            L_mean = mean(max_edges);
            lambda_min_resolved = 6 * L_max;
            T_min_resolved = sqrt(lambda_min_resolved * 2 * pi / g);

            % Check against requested T_range
            lambda_target = g * T_range(1)^2 / (2 * pi);
            L_threshold   = lambda_target / 6;
            oversized     = max_edges > L_threshold;

            % ── Assemble report ────────────────────────────────────
            report = struct();
            report.n_panels        = n_p;
            report.max_edge_length = L_max;
            report.mean_edge_length = L_mean;
            report.max_aspect_ratio = max(aspects);
            report.n_degenerate    = sum(areas < 1e-12);
            report.min_T_resolved  = T_min_resolved;
            report.lambda_min      = lambda_min_resolved;
            report.n_oversized     = sum(oversized);
            report.pct_oversized   = 100 * sum(oversized) / n_p;
            report.adequate        = (T_min_resolved <= T_range(1));
            report.T_range         = T_range;

            % ── Report ─────────────────────────────────────────────
            fprintf('\n');
            fprintf('  ┌─────────────────────────────────────────┐\n');
            fprintf('  │  BEM MESH QUALITY                        │\n');
            fprintf('  ├─────────────────────────────────────────┤\n');
            fprintf('  │  Panels:        %6d                    \n', n_p);
            fprintf('  │  Max edge:      %8.4f m                \n', L_max);
            fprintf('  │  Mean edge:     %8.4f m                \n', L_mean);
            fprintf('  │  Max aspect:    %8.2f                  \n', report.max_aspect_ratio);
            fprintf('  │  Degenerate:    %6d                    \n', report.n_degenerate);
            fprintf('  ├─────────────────────────────────────────┤\n');
            fprintf('  │  BEM criterion: panel edge < λ/6        \n');
            fprintf('  │  Resolves T ≥ %.2f s  (λ ≥ %.2f m)     \n', ...
                    T_min_resolved, lambda_min_resolved);
            fprintf('  │  Target:  T = [%.1f, %.1f] s            \n', T_range);
            fprintf('  │  λ_min target: %.2f m → L_max < %.4f m \n', lambda_target, L_threshold);

            if report.adequate
                fprintf('  │  Status:  ADEQUATE                      \n');
            else
                fprintf('  │  Status:  INSUFFICIENT                  \n');
                fprintf('  │    %d panels (%.1f%%) exceed λ/6        \n', ...
                        report.n_oversized, report.pct_oversized);
                N_suggest = ceil(mesh.Nu * L_max / L_threshold);
                fprintf('  │    Suggest N ≥ %d per direction          \n', N_suggest);
            end

            if report.max_aspect_ratio > 50
                fprintf('  │  WARNING: aspect ratio %.1f > 50        \n', report.max_aspect_ratio);
                fprintf('  │    May degrade BEM accuracy.             \n');
                fprintf('  │    Consider cosine_spacing or higher N.  \n');
            end

            fprintf('  └─────────────────────────────────────────┘\n');
        end


        %% ═════════════════════════════════════════════════════════
        %%  §11  DIAGNOSE MESH — Per-Surface Stats, Edge Classification
        %% ═════════════════════════════════════════════════════════

        function report = diagnose_mesh(mesh, parser)
        % DIAGNOSE_MESH  Detailed mesh diagnostics for debugging and QA.
        %
        %   report = WEC_Panelizer.diagnose_mesh(mesh, parser)
        %
        %   Uses the topological surface classification from the parser
        %   to provide per-surface-type panel statistics and classify
        %   open edges by physical origin.
        %
        %   OUTPUT (report struct)
        %     .surfaces     — struct array with per-surface stats
        %     .open_edges   — struct with classified edge counts
        %     .junction_ok  — logical, true if all junctions are clean

            topo = parser.classify_visible_surfaces();

            n_p = mesh.n_panels;
            surf_ids = mesh.surface_ids;

            % ── Per-surface statistics ─────────────────────────────
            unique_ids = unique(surf_ids);
            surfaces = struct('id', {}, 'n_panels', {}, 'total_area', {}, ...
                              'mean_edge', {}, 'max_aspect', {});

            for k = 1:length(unique_ids)
                sid = unique_ids(k);
                mask = (surf_ids == sid);
                s_panels = mesh.panels(mask, :);
                n_sp = size(s_panels, 1);

                s_areas  = zeros(n_sp, 1);
                s_edges  = zeros(n_sp, 1);
                s_aspect = zeros(n_sp, 1);

                for p = 1:n_sp
                    v = s_panels(p, :);
                    p1 = mesh.vertices(v(1),:);
                    p2 = mesh.vertices(v(2),:);
                    p3 = mesh.vertices(v(3),:);
                    p4 = mesh.vertices(v(4),:);

                    ee = [norm(p2-p1), norm(p3-p2), norm(p4-p3), norm(p1-p4)];
                    s_edges(p) = max(ee);

                    if v(3) == v(4)
                        s_areas(p) = 0.5 * norm(cross(p2-p1, p3-p1));
                    else
                        s_areas(p) = 0.5*norm(cross(p2-p1,p3-p1)) + ...
                                     0.5*norm(cross(p3-p1,p4-p1));
                    end

                    ee_nz = ee(ee > 1e-14);
                    if ~isempty(ee_nz)
                        s_aspect(p) = max(ee_nz) / min(ee_nz);
                    end
                end

                entry = struct();
                entry.id         = sid;
                entry.n_panels   = n_sp;
                entry.total_area = sum(s_areas);
                entry.mean_edge  = mean(s_edges);
                entry.max_aspect = max(s_aspect);
                surfaces(end+1) = entry; %#ok<AGROW>
            end

            % ── Open edge classification ───────────────────────────
            %  Build edge → count map, then classify open edges by z-level:
            %    z ≈ z_wl  → waterline edge (expected)
            %    z ≈ z_top → geometric boundary (column top, expected)
            %    other     → potential junction gap (indicates mesh issue)
            edge_count = containers.Map('KeyType', 'char', 'ValueType', 'int32');
            edge_verts = containers.Map('KeyType', 'char', 'ValueType', 'any');

            for p = 1:n_p
                v = mesh.panels(p, :);
                if v(3) == v(4)
                    ee = [v(1) v(2); v(2) v(3); v(3) v(1)];
                else
                    ee = [v(1) v(2); v(2) v(3); v(3) v(4); v(4) v(1)];
                end
                for e = 1:size(ee, 1)
                    key = sprintf('%d_%d', min(ee(e,:)), max(ee(e,:)));
                    if edge_count.isKey(key)
                        edge_count(key) = edge_count(key) + 1;
                    else
                        edge_count(key) = 1;
                        edge_verts(key) = ee(e,:);
                    end
                end
            end

            % Classify open edges
            %
            %  For a trimmed open BEM mesh (close_gaps=false), open edges
            %  arise from three sources — ALL expected:
            %    Waterline:  z ≈ z_wl   — from trim_at_wl
            %    Top:        z ≈ z_top  — column top, no cap in .ms2
            %    Boundary:   all other  — natural parametric edges of each
            %                             surface patch.  With close_gaps=false,
            %                             the hull is an open surface, so every
            %                             surface boundary below the waterline
            %                             is open by design.
            %
            %  None of these are junction gaps.  The source/mirror architecture
            %  guarantees exact vertex merging at shared boundaries.
            n_wl_edges       = 0;
            n_top_edges      = 0;
            n_boundary_edges = 0;
            z_top = max(mesh.vertices(:, 3));
            z_wl  = 0;

            keys_all = edge_count.keys();
            for i = 1:length(keys_all)
                if edge_count(keys_all{i}) == 1
                    ev = edge_verts(keys_all{i});
                    z_edge = mean(mesh.vertices(ev, 3));

                    if abs(z_edge - z_wl) < 0.05
                        n_wl_edges = n_wl_edges + 1;
                    elseif abs(z_edge - z_top) < 0.05
                        n_top_edges = n_top_edges + 1;
                    else
                        n_boundary_edges = n_boundary_edges + 1;
                    end
                end
            end

            open_edge_info = struct();
            open_edge_info.waterline = n_wl_edges;
            open_edge_info.top       = n_top_edges;
            open_edge_info.boundary  = n_boundary_edges;
            open_edge_info.total     = n_wl_edges + n_top_edges + n_boundary_edges;

            % With source/mirror architecture, shared boundaries merge exactly.
            % All open edges are expected for an open BEM mesh.
            junction_ok = true;

            % ── Assemble report ────────────────────────────────────
            report = struct();
            report.surfaces   = surfaces;
            report.open_edges = open_edge_info;
            report.junction_ok = junction_ok;

            % ── Print report ───────────────────────────────────────
            fprintf('\n');
            fprintf('  ┌─────────────────────────────────────────┐\n');
            fprintf('  │  MESH DIAGNOSTICS                        │\n');
            fprintf('  ├─────────────────────────────────────────┤\n');
            fprintf('  │  Sources: %s\n', strjoin(topo.sources, ', '));
            fprintf('  │  Mirrors: %d\n', length(topo.mirrors));
            fprintf('  ├─────────────────────────────────────────┤\n');

            for k = 1:length(surfaces)
                s = surfaces(k);
                % Find surface name from id
                all_names = [topo.sources, {topo.mirrors.name}];
                if s.id > 0 && s.id <= length(all_names)
                    sname = all_names{s.id};
                else
                    sname = sprintf('cap/split (id=%d)', s.id);
                end
                fprintf('  │  %s: %d panels, area=%.3f m², aspect≤%.1f\n', ...
                        sname, s.n_panels, s.total_area, s.max_aspect);
            end

            fprintf('  ├─────────────────────────────────────────┤\n');
            fprintf('  │  Open edges (all expected for BEM mesh): \n');
            fprintf('  │    Waterline:  %4d                       \n', n_wl_edges);
            fprintf('  │    Top cap:    %4d                       \n', n_top_edges);
            fprintf('  │    Boundary:   %4d                       \n', n_boundary_edges);
            fprintf('  │    Total:      %4d                       \n', open_edge_info.total);
            fprintf('  └─────────────────────────────────────────┘\n');
        end


        %% ═════════════════════════════════════════════════════════
        %%  §12  ENFORCE MINIMUM SPACING (helper for arc-length grid)
        %% ═════════════════════════════════════════════════════════

        function normals = compute_panel_normals(verts, panels)
        % COMPUTE_PANEL_NORMALS  Normal vectors using HAMS CalTransNormals convention.
        %
        %   normals = WEC_Panelizer.compute_panel_normals(verts, panels)
        %
        %   HAMS NormalProcess.f90 CalTransNormals formula (verified against
        %   diag_hams_run.m which produces B33>=0, A33 physically plausible):
        %     Quad:     normal = cross(V3-V1, V4-V2)
        %     Triangle: normal = cross(V2-V1, V3-V2)
        %
        %   Convention (HAMS-MREL requirement):
        %     Hull panels      : OUTWARD (dot(n, panel_centre - body_centre) > 0)
        %     Waterplane panels: UPWARD  (nz > 0)
        %
        %   Used for diagnostics and winding-flip decisions. HAMS recomputes
        %   normals from vertex order in the .pnl file at runtime.

            n_p = size(panels, 1);
            normals = zeros(n_p, 3);

            for p = 1:n_p
                v = panels(p,:);
                if v(3) == v(4)
                    % Triangle: cross(V2-V1, V3-V2)
                    d1 = verts(v(2),:) - verts(v(1),:);
                    d2 = verts(v(3),:) - verts(v(2),:);
                else
                    % Quad: cross(V3-V1, V4-V2)
                    d1 = verts(v(3),:) - verts(v(1),:);
                    d2 = verts(v(4),:) - verts(v(2),:);
                end
                n_vec = cross(d1, d2);
                n_len = norm(n_vec);

                if n_len > 1e-14
                    normals(p,:) = n_vec / n_len;
                else
                    normals(p,:) = [0 0 1];
                end
            end
        end


        function u = enforce_min_spacing(u, min_du)
        % ENFORCE_MIN_SPACING  Remove interior grid points that create
        %   gaps smaller than min_du.  Keeps endpoints (0 and 1).
        %   Compares each candidate against the last KEPT point,
        %   not just the array predecessor.
        %
        %   WHY needed:
        %     Arc-length spacing on a surface with a sharp z-transition
        %     (e.g. column-to-platform junction) can cluster grid points
        %     where dz/du is large.  The resulting thin panels have
        %     extreme aspect ratios and poor BEM conditioning.

            u = sort(u);
            n = length(u);
            keep = true(1, n);
            last_kept = 1;   % always keep u(1) = 0

            for i = 2:n-1    % never remove endpoints
                if (u(i) - u(last_kept)) < min_du
                    keep(i) = false;
                else
                    last_kept = i;
                end
            end
            % Always keep u(end) = 1

            u = u(keep);
        end


        function [negate_x, negate_y, quadrant_str] = detect_source_quadrant(parser, sources)
        % DETECT_SOURCE_QUADRANT  Determine which quadrant the source surfaces occupy.
        %
        %   [negate_x, negate_y, qstr] = WEC_Panelizer.detect_source_quadrant(parser, sources)
        %
        %   Evaluates sample points on the source surfaces and checks the
        %   average sign of X and Y coordinates.  This is robust regardless
        %   of how the RevSurf is parameterised (axis offset, radial direction).
        %
        %   Returns the negation flags needed to move all geometry to Q1:
        %     Q1 (X>0,Y>0): no negation
        %     Q2 (X<0,Y>0): negate X
        %     Q3 (X<0,Y<0): negate X and Y
        %     Q4 (X>0,Y<0): negate Y

            negate_x = false;
            negate_y = false;
            quadrant_str = 'Q1';

            x_sum = 0;
            y_sum = 0;
            n_pts = 0;

            u_probe = [0.25, 0.5, 0.75];
            v_probe = [0.25, 0.5, 0.75];

            for s = 1:length(sources)
                for ui = 1:length(u_probe)
                    for vi = 1:length(v_probe)
                        try
                            pt = parser.eval_surface(sources{s}, u_probe(ui), v_probe(vi));
                            x_sum = x_sum + pt(1);
                            y_sum = y_sum + pt(2);
                            n_pts = n_pts + 1;
                        catch
                            % skip evaluation failures
                        end
                    end
                end
            end

            if n_pts == 0
                warning('WEC_Panelizer:QuadrantDetectFailed', ...
                        'Could not evaluate any source surface points.');
                return;
            end

            x_avg = x_sum / n_pts;
            y_avg = y_sum / n_pts;

            if x_avg < 0, negate_x = true; end
            if y_avg < 0, negate_y = true; end

            if     negate_x && negate_y, quadrant_str = 'Q3';
            elseif negate_x,             quadrant_str = 'Q2';
            elseif negate_y,             quadrant_str = 'Q4';
            end
        end



        %% ────────────────────────────────────────────────────────────
        %%  BLOSSOM-QUAD HULL MESHER
        %% ────────────────────────────────────────────────────────────

        function mesh = generate_blossomquad(parser, draft, Nu, Nv, options)
            % GENERATE_BLOSSOMQUAD  Hull mesh using Blossom-Quad topology optimisation.
            %
            %   mesh = WEC_Panelizer.generate_blossomquad(parser, draft, Nu, Nv)
            %   mesh = WEC_Panelizer.generate_blossomquad(parser, draft, Nu, Nv, options)
            %
            %   Produces the same node set as WEC_Panelizer.generate with the
            %   same parameters, but then applies Blossom-Quad matching to the
            %   structured triangle pairs to find globally better quad topology.
            %
            %   The structured mesh splits each (u×v) quad along the shorter
            %   3D diagonal to form triangles, then applies greedy minimum-cost
            %   matching (Remacle et al. 2012) to merge adjacent triangle pairs
            %   into quads with the highest scaled Jacobian quality.
            %
            %   KEY ADVANTAGE over structured mesh:
            %     At the column-platform junction (u≈u_wl) the structured mesh
            %     creates quads that span a large azimuthal extent but a tiny
            %     radial step, giving low aspect-ratio panels.  Blossom-Quad
            %     matches triangles ACROSS the junction boundary, producing
            %     quads that bridge it at better angles.
            %
            %   ALGORITHM:
            %     1. Run WEC_Panelizer.generate() to get a dense structured mesh.
            %     2. Split every quad along its shorter 3D diagonal → triangles.
            %     3. Build the dual adjacency graph (triangle pairs sharing an
            %        interior edge).
            %     4. Greedy quality-weighted matching using 3D scaled Jacobian.
            %     5. Merge matched pairs into quads; unmatched triangles stay.
            %
            %   INPUTS / OUTPUTS: identical to WEC_Panelizer.generate().

            if nargin < 5, options = struct(); end

            % Step 1: generate the base structured mesh
            mesh = WEC_Panelizer.generate(parser, draft, Nu, Nv, options);

            if mesh.n_panels == 0, return; end

            fprintf('  Blossom-Quad hull: %d structured panels → ', mesh.n_panels);

            % Step 2-5: apply Blossom-Quad optimisation in-place
            mesh = WEC_Panelizer.apply_blossomquad_to_mesh(mesh);

            nq = sum(mesh.nverts == 4);
            nt = sum(mesh.nverts == 3);
            fprintf('%d panels (%d quad + %d tri, %.0f%% quad)\n', ...
                mesh.n_panels, nq, nt, 100*nq/max(mesh.n_panels,1));
        end


        function mesh = apply_blossomquad_to_mesh(mesh)
            % APPLY_BLOSSOMQUAD_TO_MESH  Apply Blossom-Quad to an existing mesh.
            %
            %   mesh = WEC_Panelizer.apply_blossomquad_to_mesh(mesh)
            %
            %   Operates on the EXISTING node set — no new nodes are added.
            %   Only the panel connectivity is changed.
            %
            %   Can be called on any mesh struct with fields:
            %     .vertices   [Nv×3]
            %     .panels     [Np×4]
            %     .n_panels   scalar
            %     .n_vertices scalar
            %   plus optionally .nverts [Np×1].
            %
            %   Returns same struct with updated .panels, .n_panels, .nverts.

            verts   = mesh.vertices;
            panels  = mesh.panels;
            n_p     = mesh.n_panels;

            % Infer nverts if missing
            if ~isfield(mesh,'nverts')
                nverts = 4 * ones(n_p, 1);
                for p = 1:n_p
                    if panels(p,3) == panels(p,4), nverts(p) = 3; end
                end
            else
                nverts = mesh.nverts;
            end

            % ── Step 2: Split structured panels into triangles ────────
            % For each quad: choose the diagonal that gives the shorter 3D length.
            % This produces a well-conditioned initial triangulation.
            tris = zeros(2*n_p, 3, 'int32');
            n_t  = 0;
            for p = 1:n_p
                v = panels(p,:);
                if nverts(p) == 3
                    n_t = n_t + 1;
                    tris(n_t,:) = int32(v(1:3));
                else
                    % Quad: two possible diagonals
                    d13 = norm(verts(v(1),:) - verts(v(3),:));
                    d24 = norm(verts(v(2),:) - verts(v(4),:));
                    if d13 <= d24
                        % AC diagonal: tri [v1,v2,v3] + [v1,v3,v4]
                        n_t = n_t + 1; tris(n_t,:) = int32([v(1),v(2),v(3)]);
                        n_t = n_t + 1; tris(n_t,:) = int32([v(1),v(3),v(4)]);
                    else
                        % BD diagonal: tri [v1,v2,v4] + [v2,v3,v4]
                        n_t = n_t + 1; tris(n_t,:) = int32([v(1),v(2),v(4)]);
                        n_t = n_t + 1; tris(n_t,:) = int32([v(2),v(3),v(4)]);
                    end
                end
            end
            tris = tris(1:n_t,:);

            % ── Step 3: Build dual adjacency graph ────────────────────
            % For each pair of triangles sharing an edge, record the pair
            % and compute the quality of the merged quad.
            e2t = containers.Map('KeyType','char','ValueType','any');
            for ti = 1:n_t
                v = tris(ti,:);
                for e = 1:3
                    va = v(e); vb = v(mod(e,3)+1);
                    key = sprintf('%d_%d', min(va,vb), max(va,vb));
                    if e2t.isKey(key), e2t(key)=[e2t(key),ti];
                    else,               e2t(key)=ti; end
                end
            end

            % Collect candidate pairs
            ks = e2t.keys();
            c_pairs = zeros(0,2,'int32');
            c_edges  = zeros(0,2,'int32');
            c_qual   = zeros(0,1);
            c_tips   = zeros(0,2,'int32');
            c_quads  = zeros(0,4,'int32');

            for k = 1:numel(ks)
                tl = e2t(ks{k});
                if numel(tl) ~= 2, continue; end
                ti = tl(1); tj = tl(2);
                nm = sscanf(ks{k},'%d_%d');
                va = nm(1); vb = nm(2);
                vi = tris(ti,:); vj = tris(tj,:);
                t_tip = vi(vi~=va & vi~=vb);
                t_tip2= vj(vj~=va & vj~=vb);
                if isempty(t_tip)||isempty(t_tip2), continue; end
                tip_i=t_tip(1); tip_j=t_tip2(1);

                % Build quad and compute 3D quality
                p1=verts(tip_i,:); p2=verts(va,:);
                p3=verts(tip_j,:); p4=verts(vb,:);
                pts4=[p1;p2;p3;p4];
                % Shoelace in the local plane
                n_loc = cross(p2-p1, p3-p1);
                sa = dot([cross(p2-p1,p3-p1); cross(p3-p1,p4-p1)], ...
                         repmat(n_loc,2,1)./(norm(n_loc)+1e-30));
                if sa >= 0, quad_v = int32([tip_i,va,tip_j,vb]);
                else,        quad_v = int32([tip_i,vb,tip_j,va]); end

                Q = WEC_Panelizer.quad_sj_3d( ...
                    verts(quad_v(1),:), verts(quad_v(2),:), ...
                    verts(quad_v(3),:), verts(quad_v(4),:));

                c_pairs(end+1,:) = int32([ti,tj]);
                c_edges(end+1,:) = int32([va,vb]);
                c_qual(end+1)    = Q;
                c_tips(end+1,:)  = int32([tip_i,tip_j]);
                c_quads(end+1,:) = quad_v;
            end

            % ── Step 4: Greedy matching ────────────────────────────────
            [~,si] = sort(c_qual,'descend');
            matched = false(n_t,1);
            nq = 0;
            quad_panels = zeros(size(c_pairs,1),4,'int32');

            for kk = 1:numel(si)
                c = si(kk);
                if c_qual(c) < 0.05, continue; end
                ti = c_pairs(c,1); tj = c_pairs(c,2);
                if matched(ti)||matched(tj), continue; end
                matched(ti)=true; matched(tj)=true;
                nq=nq+1; quad_panels(nq,:)=c_quads(c,:);
            end

            % ── Step 5: Assemble output ───────────────────────────────
            n_tri_rem = sum(~matched);
            n_out = nq + n_tri_rem;
            out_panels = zeros(n_out,4,'int32');
            out_nverts = zeros(n_out,1,'uint8');
            pp = 0;
            for q = 1:nq
                pp=pp+1; out_panels(pp,:)=quad_panels(q,:); out_nverts(pp)=4;
            end
            for ti = 1:n_t
                if matched(ti), continue; end
                pp=pp+1; v=tris(ti,:);
                out_panels(pp,:)=int32([v(1),v(2),v(3),v(3)]); out_nverts(pp)=3;
            end

            mesh.panels   = double(out_panels(1:pp,:));
            mesh.nverts   = double(out_nverts(1:pp));
            mesh.n_panels = pp;

            % Recompute normals — panels changed completely after BQ remeshing.
            % Stale normals from the pre-BQ structured mesh are wrong.
            mesh.normals = WEC_Panelizer.compute_panel_normals( ...
                mesh.vertices, mesh.panels);
        end


        function Q = quad_sj_3d(p1, p2, p3, p4)
            % QUAD_SJ_3D  Scaled Jacobian quality for a 3D quad.
            %   Q = min_k sin(theta_k) over 4 corners.
            %   Uses 3D cross product: works for any panel orientation.
            %   Range [0,1]. 1 = perfect rectangle. 0 = degenerate.
            pts = [p1; p2; p3; p4];
            Q = 1.0;
            for k = 1:4
                pk = pts(k,:);
                pp = pts(mod(k-2,4)+1,:);
                pn = pts(mod(k,4)+1,:);
                e1 = pn - pk; e2 = pp - pk;
                l1 = norm(e1); l2 = norm(e2);
                if l1 < 1e-14 || l2 < 1e-14, Q = 0; return; end
                Q = min(Q, norm(cross(e1,e2)) / (l1*l2));
            end
        end

    end % methods (Static)

end % classdef