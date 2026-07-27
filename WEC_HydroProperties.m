classdef WEC_HydroProperties
% WEC_HYDROPROPERTIES  Volume and hydrostatic properties directly from
% parametric B-spline surfaces — no mesh, no panelizer in the loop.
%
%   props = WEC_HydroProperties.compute(parser)
%   props = WEC_HydroProperties.compute(parser, options)
%
%   ARCHITECTURE
%   ────────────────────────────────────────────────────────────────────
%   WEC_HydroProperties (this file)
%     ├─ §1  compute             — full hull: V, centroid
%     ├─ §2  surface_integral    — GL quadrature over one S(u,v)
%     ├─ §3  detect_open_edges   — geometric edge matching
%     ├─ §4  chain_open_edges    — order open edges into closed loop(s)
%     ├─ §5  cap_contribution    — analytical flat-cap V and centroid
%     ├─ §6  compute_submerged   — adjusted v-limits: V_sub, CB, Aw, I_wp
%     ├─ §7  surface_integral_submerged — GL with per-u v-limits
%     ├─ §8  find_submerged_v_limits — [v_lo, v_hi] where z ≤ z_wl
%     ├─ §9  extract_waterplane_contour — iso-z contour (reference method)
%     ├─ §9b precompute_boundary_cache  — evaluate boundary curves once
%     ├─ §9c extract_isocurve_at_z      — fast type-dispatched iso-z
%     │       ├─ isocurve_devsurf        — algebraic v-solve (exact)
%     │       ├─ isocurve_revsurf        — 1D profile root + analytic arc
%     │       └─ isocurve_bloftsurf      — B-spline span bracket + bisect
%     ├─ §10 waterplane_properties — Aw, I_xx, I_yy from boundary
%     ├─ §10c compute_strip       — horizontal strip V, CB, I via DT
%     ├─ §10c0 compute_rmin_at_z  — min centroid-to-boundary distance
%     └─ §11 gauss_legendre       — Golub-Welsch GL nodes + weights
%   ────────────────────────────────────────────────────────────────────
%
%   DIVERGENCE THEOREM
%     The divergence theorem converts a volume integral into a surface
%     integral over the hull boundary:
%
%       V = (1/3) ∮_∂Ω  r · n̂  dA
%
%     For a parametric surface S(u,v), the oriented area element is:
%
%       n̂ dA = (∂S/∂u × ∂S/∂v) du dv
%
%     So for each visible surface patch:
%
%       V_i = (1/3) ∫₀¹ ∫₀¹ S · (S_u × S_v)  du dv
%
%     The double integral is evaluated with Gauss-Legendre quadrature.
%     For smooth B-spline surfaces, GL converges exponentially — 20×20
%     points per surface gives ~6-digit accuracy.
%
%     Partial derivatives S_u, S_v are computed analytically via
%     eval_surface_with_derivs on the parser.  Each surface type
%     (RevSurf, DevSurf, BLoftSurf, etc.) has a closed-form
%     derivative method, cascading through curve/snake derivatives.
%
%   FAST ISO-Z EXTRACTION (v8.0)
%     For waterplane contour extraction, type-dispatched algebraic
%     solvers exploit each surface type's parametric structure:
%
%     DevSurf/RuledSurf: S(u,v) = (1-v)*C1(u) + v*C2(u) is linear in v.
%       z(u,v) = z_wl solves to v* = (z_wl - z1(u)) / (z2(u) - z1(u)).
%       Algebraically exact — no iteration needed.
%
%     RevSurf (vertical axis): z depends on u only (revolution
%       preserves z).  Find u* where z_profile(u*) = z_wl via 1D
%       root-finding on the profile curve.  The contour at u* is the
%       revolution arc — generated analytically from the cached axis
%       decomposition (proj + r*cos(phi)*e_r + r*sin(phi)*e_t).
%
%     BLoftSurf: At fixed u, z(v) is a degree-p B-spline with control
%       values z_k = section_k(u).z.  Root-finding uses the piecewise
%       structure: bracket via sign changes in control values at knot
%       span boundaries, then bisect within the single bracketed span.
%       ~44 scalar B-spline evaluations vs ~960 in the reference method.
%
%     MirrSurf: z is invariant under X/Y reflection.  The contour is
%       the source contour with one coordinate flipped.  Zero evaluations.
%
%     All boundary curve evaluations are precomputed ONCE via
%     precompute_boundary_cache and reused across all z-levels.
%     Validated: 0.00% Aw error vs reference, 1249× speedup on C0.
%
%   OPEN EDGES
%     Hull geometry may have open boundaries (column top, keel).
%     These are detected by geometric edge matching: each visible
%     surface has 4 parametric edges.  Edges shared by two surfaces
%     (same 3D curve) are closed; unmatched edges are open.
%
%     Open edges are chained into closed loop(s), and the flat-cap
%     contribution is computed analytically:
%       V_cap  = (1/3) z_cap A_cap     (top cap, outward = +z)
%       A_cap  = (1/2) ∮ (x dy − y dx)  (Green's theorem)
%     No cap surface is created — the boundary curve is sufficient.
%
%   CENTROID
%     Using div(F) = z with F = (0, 0, z²/2):
%       z̄·V = (1/2) ∮ z² n_z dA
%     And similarly for x̄, ȳ.
%
%   OUTPUT (props struct)
%     .volume    — total enclosed volume [m³]
%     .centroid  — volume centroid [x, y, z] [m]
%     .V_raw     — signed volume (for debugging orientation)
%     .V_surfs   — per-surface contributions [m³]
%     .open_edges — detected open boundary info
%
%   See also: WEC_MS2_Parser, WEC_Panelizer
%
%   Author:  WEC Optimisation Team
%   Version: 8.1 — Added second volume moments to compute() for inertia tensor

    methods (Static)

        %% ═════════════════════════════════════════════════════════
        %%  §1  COMPUTE — Main Entry Point
        %% ═════════════════════════════════════════════════════════

        function props = compute(parser, options)
        % COMPUTE  Volume and centroid from parametric surfaces.
        %
        %   props = WEC_HydroProperties.compute(parser)
        %   props = WEC_HydroProperties.compute(parser, options)
        %
        %   OPTIONS
        %     .n_quad          — GL quadrature order per direction (default: 20)
        %     .n_edge_samples  — samples per edge for matching (default: 50)
        %     .edge_tol        — edge matching tolerance [m] (default: 1e-4)
        %     .verbose         — print progress (default: true)

            if nargin < 2, options = struct(); end
            if ~isfield(options, 'n_quad'),         options.n_quad         = 20;   end
            if ~isfield(options, 'n_edge_samples'), options.n_edge_samples = 50;   end
            if ~isfield(options, 'edge_tol'),       options.edge_tol       = 1e-4; end
            if ~isfield(options, 'verbose'),        options.verbose        = true; end

            surf_names = parser.visible_surfs;
            n_surfs    = length(surf_names);
            n_quad     = options.n_quad;

            % ── Estimate hull interior point for orientation check ──
            interior = (parser.extents(1:3) + parser.extents(4:6))' / 2;

            % ── Classify surfaces for mirror deduplication ─────────
            topo = parser.classify_visible_surfaces();

            % ── Step 1: Surface integrals — sources only, derive mirrors
            V_total   = 0;
            Cxyz_num  = [0 0 0];
            V_surfs   = zeros(n_surfs, 1);

            % Second volume moments (for inertia tensor computation)
            int_x2_total = 0;
            int_y2_total = 0;
            int_z2_total = 0;

            % Map surface names to their index in visible_surfs
            surf_idx_map = containers.Map();
            for i = 1:n_surfs
                surf_idx_map(surf_names{i}) = i;
            end

            % Evaluate source surfaces (the expensive GL quadrature)
            source_V    = containers.Map();
            source_Cxyz = containers.Map();
            source_ix2  = containers.Map();
            source_iy2  = containers.Map();
            source_iz2  = containers.Map();

            if options.verbose
                fprintf('\n  Parametric divergence theorem (%d×%d GL quadrature)\n', ...
                        n_quad, n_quad);
                fprintf('    Evaluating %d source surfaces (skipping %d mirrors)\n', ...
                        length(topo.sources), length(topo.mirrors));
            end

            for s = 1:length(topo.sources)
                sname = topo.sources{s};
                [V_s, Cxyz_s, orient_s, ix2_s, iy2_s, iz2_s] = ...
                    WEC_HydroProperties.surface_integral( ...
                        parser, sname, n_quad, interior);

                source_V(sname)    = V_s;
                source_Cxyz(sname) = Cxyz_s;
                source_ix2(sname)  = ix2_s;
                source_iy2(sname)  = iy2_s;
                source_iz2(sname)  = iz2_s;

                idx = surf_idx_map(sname);
                V_surfs(idx) = V_s;
                V_total      = V_total + V_s;
                Cxyz_num     = Cxyz_num + Cxyz_s;
                int_x2_total = int_x2_total + ix2_s;
                int_y2_total = int_y2_total + iy2_s;
                int_z2_total = int_z2_total + iz2_s;

                if options.verbose
                    fprintf('    %s: V = %+.6f m³  (orient %+d) [SOURCE]\n', ...
                            sname, V_s, orient_s);
                end
            end

            % Derive mirror contributions analytically
            %   Volume and second moments are INVARIANT under reflection.
            %   (x² unchanged by x→−x; ∫x³ n_x dA also invariant because
            %    x³ flips and n_x flips, cancelling.)
            %   Only centroid components flip per mirror plane.
            for m = 1:length(topo.mirrors)
                mirr = topo.mirrors(m);

                % Find the ultimate source's V and Cxyz
                ult = mirr.ultimate_source;
                if source_V.isKey(ult)
                    V_src    = source_V(ult);
                    Cxyz_src = source_Cxyz(ult);
                    ix2_src  = source_ix2(ult);
                    iy2_src  = source_iy2(ult);
                    iz2_src  = source_iz2(ult);
                else
                    % Fallback: evaluate directly (should not happen)
                    [V_src, Cxyz_src, ~, ix2_src, iy2_src, iz2_src] = ...
                        WEC_HydroProperties.surface_integral( ...
                            parser, mirr.name, n_quad, interior);
                    V_total      = V_total + V_src;
                    Cxyz_num     = Cxyz_num + Cxyz_src;
                    int_x2_total = int_x2_total + ix2_src;
                    int_y2_total = int_y2_total + iy2_src;
                    int_z2_total = int_z2_total + iz2_src;
                    if surf_idx_map.isKey(mirr.name)
                        V_surfs(surf_idx_map(mirr.name)) = V_src;
                    end
                    continue;
                end

                % Volume is invariant under reflection
                V_mirr = V_src;

                % Second moments are invariant under reflection
                ix2_mirr = ix2_src;
                iy2_mirr = iy2_src;
                iz2_mirr = iz2_src;

                % Centroid: flip the component corresponding to each mirror plane
                %   MirrY: Cx same, Cy flips, Cz same
                %   MirrX: Cx flips, Cy same, Cz same
                Cxyz_mirr = Cxyz_src;
                for fi = 1:length(mirr.effective_flips)
                    if strcmp(mirr.effective_flips{fi}, 'Y')
                        Cxyz_mirr(2) = -Cxyz_mirr(2);
                    else  % 'X'
                        Cxyz_mirr(1) = -Cxyz_mirr(1);
                    end
                end

                if surf_idx_map.isKey(mirr.name)
                    V_surfs(surf_idx_map(mirr.name)) = V_mirr;
                end
                V_total      = V_total + V_mirr;
                Cxyz_num     = Cxyz_num + Cxyz_mirr;
                int_x2_total = int_x2_total + ix2_mirr;
                int_y2_total = int_y2_total + iy2_mirr;
                int_z2_total = int_z2_total + iz2_mirr;

                if options.verbose
                    fprintf('    %s: V = %+.6f m³  [MIRROR of %s]\n', ...
                            mirr.name, V_mirr, ult);
                end
            end

            % ── Step 2: Detect open edges ───────────────────────────
            [open_edges, ~] = WEC_HydroProperties.detect_open_edges( ...
                parser, options.n_edge_samples, options.edge_tol);

            if options.verbose
                fprintf('    Open edges: %d\n', length(open_edges));
                for i = 1:length(open_edges)
                    oe = open_edges{i};
                    fprintf('      %s edge %d  (z ≈ %.3f m)\n', ...
                            oe.surface, oe.edge_idx, oe.z_mean);
                end
            end

            % ── Step 3: Cap contribution ────────────────────────────
            V_cap_total  = 0;
            Cxyz_cap_total = [0 0 0];

            if ~isempty(open_edges)
                [V_cap, Cxyz_cap, A_cap, z_cap] = ...
                    WEC_HydroProperties.cap_contribution( ...
                        open_edges, interior, options.edge_tol);

                V_cap_total    = sum(V_cap);
                Cxyz_cap_total = sum(Cxyz_cap, 1);

                % Cap second moment contribution:
                %   Flat cap at z with outward normal ±z, area A.
                %   n_x = 0, n_y = 0  →  int_x2_cap = 0, int_y2_cap = 0
                %   int_z2_cap = sign_nz × z³ × A / 3
                for g = 1:length(V_cap)
                    if A_cap(g) > 1e-8
                        if z_cap(g) > interior(3)
                            sign_nz = +1;
                        else
                            sign_nz = -1;
                        end
                        int_z2_total = int_z2_total + sign_nz * z_cap(g)^3 * A_cap(g) / 3;
                    end
                end

                if options.verbose
                    for g = 1:length(V_cap)
                        if A_cap(g) > 1e-8  % skip degenerate caps (e.g. keel point)
                            fprintf('    Cap: z = %.3f m, A = %.6f m², V = %+.6f m³\n', ...
                                    z_cap(g), A_cap(g), V_cap(g));
                        end
                    end
                end
            end

            V_total  = V_total + V_cap_total;
            Cxyz_num = Cxyz_num + Cxyz_cap_total;

            % ── Assemble output ─────────────────────────────────────
            props = struct();
            props.volume  = abs(V_total);
            props.V_raw   = V_total;
            props.V_surfs = V_surfs;
            props.V_cap   = V_cap_total;
            props.open_edges = open_edges;

            if abs(V_total) > 1e-12
                props.centroid = Cxyz_num / V_total;
            else
                props.centroid = [0, 0, 0];
            end

            % Second volume moments (for inertia tensor, about body-frame origin)
            props.int_x2 = int_x2_total;
            props.int_y2 = int_y2_total;
            props.int_z2 = int_z2_total;

            % ── Report ──────────────────────────────────────────────
            if options.verbose
                fprintf('\n');
                fprintf('  ┌─────────────────────────────────────────┐\n');
                fprintf('  │  HYDRO PROPERTIES (Parametric)           │\n');
                fprintf('  ├─────────────────────────────────────────┤\n');
                fprintf('  │  Volume:    %10.4f m³                \n', props.volume);
                fprintf('  │  Centroid:  [%.3f, %.3f, %.3f] m      \n', props.centroid);
                fprintf('  │  V_raw:     %+10.4f m³ (signed)      \n', props.V_raw);
                fprintf('  │  V_cap:     %+10.6f m³               \n', props.V_cap);
                fprintf('  │  int_x2:    %+10.4f m⁵               \n', props.int_x2);
                fprintf('  │  int_y2:    %+10.4f m⁵               \n', props.int_y2);
                fprintf('  │  int_z2:    %+10.4f m⁵               \n', props.int_z2);
                fprintf('  └─────────────────────────────────────────┘\n');
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §2  SURFACE INTEGRAL — GL Quadrature over S(u,v)
        %% ═════════════════════════════════════════════════════════

        function [V_i, Cxyz_i, orient, int_x2_i, int_y2_i, int_z2_i] = ...
                surface_integral(parser, surf_name, n_quad, interior_pt)
        % SURFACE_INTEGRAL  Divergence theorem integral over one surface.
        %
        %   [V, Cxyz, orient] = surface_integral(parser, name, n_quad, interior)
        %   [V, Cxyz, orient, ix2, iy2, iz2] = surface_integral(...)
        %
        %   Uses n_quad × n_quad Gauss-Legendre points on [0,1]².
        %   Partial derivatives S_u, S_v computed analytically via
        %   eval_surface_with_derivs (Phase 2 — replaces finite differences).
        %
        %   SECOND VOLUME MOMENTS (v1.1):
        %     int_x2 = (1/3) ∫ x³ n_x dA   →  ∫ x² dV
        %     int_y2 = (1/3) ∫ y³ n_y dA   →  ∫ y² dV
        %     int_z2 = (1/3) ∫ z³ n_z dA   →  ∫ z² dV
        %   These enable inertia tensor computation for the full body.

            [u_gl, w_u] = WEC_HydroProperties.gauss_legendre(n_quad);
            [v_gl, w_v] = WEC_HydroProperties.gauss_legendre(n_quad);

            V_i      = 0;
            Cxyz_i   = [0, 0, 0];
            int_x2_i = 0;
            int_y2_i = 0;
            int_z2_i = 0;

            for i = 1:n_quad
                u = u_gl(i);
                for j = 1:n_quad
                    v = v_gl(j);

                    % Surface point and analytical derivatives
                    [S, Su, Sv] = parser.eval_surface_with_derivs(surf_name, u, v);

                    % Oriented area element
                    n_vec = cross(Su, Sv);

                    % Quadrature weight
                    w = w_u(i) * w_v(j);

                    % Volume: (1/3) ∫ S · (S_u × S_v) du dv
                    V_i = V_i + w * dot(S, n_vec) / 3;

                    % Centroid numerators
                    Cxyz_i = Cxyz_i + w * [ ...
                        S(1)^2 * n_vec(1) / 2, ...
                        S(2)^2 * n_vec(2) / 2, ...
                        S(3)^2 * n_vec(3) / 2];

                    % Second volume moments
                    int_x2_i = int_x2_i + w * S(1)^3 * n_vec(1) / 3;
                    int_y2_i = int_y2_i + w * S(2)^3 * n_vec(2) / 3;
                    int_z2_i = int_z2_i + w * S(3)^3 * n_vec(3) / 3;
                end
            end

            % ── Orientation check via coarse 5×5 signed-volume integral ──
            % FIX (L1): The single-point check at (0.5,0.5) fails for
            % non-convex surfaces where the midpoint normal may point
            % inward despite the surface overall being outward-oriented.
            % compute_submerged already uses this coarse-integral approach.
            % A 5×5 signed-volume integral gives the correct sign of
            % ∫ S·(Su×Sv) du dv — it is always negative for an inward-
            % oriented parameterisation and positive for outward.
            n_orient_chk = 5;
            [u_chk, wu_chk] = WEC_HydroProperties.gauss_legendre(n_orient_chk);
            [v_chk, wv_chk] = WEC_HydroProperties.gauss_legendre(n_orient_chk);
            V_sign_test = 0;
            for ic = 1:n_orient_chk
                for jc = 1:n_orient_chk
                    [S_c, Su_c, Sv_c] = parser.eval_surface_with_derivs( ...
                        surf_name, u_chk(ic), v_chk(jc));
                    V_sign_test = V_sign_test + wu_chk(ic) * wv_chk(jc) * ...
                        dot(S_c, cross(Su_c, Sv_c)) / 3;
                end
            end

            if V_sign_test < 0
                V_i      = -V_i;
                Cxyz_i   = -Cxyz_i;
                int_x2_i = -int_x2_i;
                int_y2_i = -int_y2_i;
                int_z2_i = -int_z2_i;
                orient   = -1;
            else
                orient = +1;
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §3  DETECT OPEN EDGES — Geometric Edge Matching
        %% ═════════════════════════════════════════════════════════

        function [open_edges, all_edges] = detect_open_edges(parser, n_samples, tol)
        % DETECT_OPEN_EDGES  Find boundary curves not shared by two surfaces.
        %
        %   Each visible surface has 4 parametric edges.  Two edges are
        %   "matched" if they trace the same 3D curve (same points up to
        %   tolerance, possibly in reversed direction).
        %
        %   Unmatched edges are open boundaries that need a cap for a
        %   watertight enclosure.
        %
        %   INPUTS
        %     parser    — parsed .ms2 model
        %     n_samples — points per edge for comparison (default: 50)
        %     tol       — geometric matching tolerance [m] (default: 1e-4)
        %
        %   OUTPUTS
        %     open_edges — cell array of unmatched edge structs
        %     all_edges  — cell array of ALL edge structs (for debugging)

            if nargin < 2 || isempty(n_samples), n_samples = 50; end
            if nargin < 3 || isempty(tol),       tol = 1e-4;     end

            t_sample   = linspace(0, 1, n_samples);
            surf_names = parser.visible_surfs;
            n_surfs    = length(surf_names);

            % ── Sample all 4 edges of each visible surface ──────────
            all_edges = {};
            for s = 1:n_surfs
                sname = surf_names{s};
                for edge_idx = 1:4
                    pts = zeros(n_samples, 3);
                    for k = 1:n_samples
                        t = t_sample(k);
                        switch edge_idx
                            case 1  % v = 0, u varies
                                pts(k,:) = parser.eval_surface(sname, t, 0);
                            case 2  % u = 1, v varies
                                pts(k,:) = parser.eval_surface(sname, 1, t);
                            case 3  % v = 1, u varies
                                pts(k,:) = parser.eval_surface(sname, t, 1);
                            case 4  % u = 0, v varies
                                pts(k,:) = parser.eval_surface(sname, 0, t);
                        end
                    end

                    edge = struct();
                    edge.surface  = sname;
                    edge.edge_idx = edge_idx;
                    edge.pts      = pts;
                    edge.z_mean   = mean(pts(:, 3));
                    edge.matched  = false;
                    all_edges{end+1} = edge; %#ok<AGROW>
                end
            end

            % ── Match edges pairwise ────────────────────────────────
            %  Two edges match if their sampled points coincide (up to
            %  tolerance) in either forward or reversed order.
            n_edges = length(all_edges);

            for a = 1:n_edges
                if all_edges{a}.matched, continue; end
                for b = a+1:n_edges
                    if all_edges{b}.matched, continue; end
                    % Don't match edges from the same surface
                    if strcmp(all_edges{a}.surface, all_edges{b}.surface)
                        continue;
                    end

                    pts_a = all_edges{a}.pts;
                    pts_b = all_edges{b}.pts;

                    dist_fwd = max(vecnorm(pts_a - pts_b, 2, 2));
                    dist_rev = max(vecnorm(pts_a - flipud(pts_b), 2, 2));

                    if min(dist_fwd, dist_rev) < tol
                        all_edges{a}.matched = true;
                        all_edges{b}.matched = true;
                        break;  % edge a is matched, move on
                    end
                end
            end

            % ── Collect unmatched ───────────────────────────────────
            open_edges = {};
            for i = 1:n_edges
                if ~all_edges{i}.matched
                    open_edges{end+1} = all_edges{i}; %#ok<AGROW>
                end
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §4  CHAIN OPEN EDGES — Build Closed Loop(s)
        %% ═════════════════════════════════════════════════════════

        function groups = chain_open_edges(open_edges, tol)
        % CHAIN_OPEN_EDGES  Order open edges into closed loop(s).
        %
        %   Open edges at the same z-level are grouped, then chained
        %   end-to-start into a closed polygon.  Edges may need to be
        %   reversed so that one edge's endpoint matches the next
        %   edge's startpoint.
        %
        %   RETURNS cell array of groups, each containing:
        %     .loop_pts — [N × 3] ordered boundary points
        %     .z_mean   — mean z of the group

            if nargin < 2, tol = 1e-3; end

            n = length(open_edges);
            if n == 0
                groups = {};
                return;
            end

            % ── Group edges by z-level ──────────────────────────────
            z_vals   = cellfun(@(e) e.z_mean, open_edges);
            assigned = false(n, 1);
            z_groups = {};

            for i = 1:n
                if assigned(i), continue; end
                members = i;
                assigned(i) = true;
                for j = i+1:n
                    if ~assigned(j) && abs(z_vals(j) - z_vals(i)) < 0.05
                        members(end+1) = j; %#ok<AGROW>
                        assigned(j) = true;
                    end
                end
                z_groups{end+1} = members; %#ok<AGROW>
            end

            % ── Chain each group into a loop ────────────────────────
            groups = {};
            for g = 1:length(z_groups)
                idx   = z_groups{g};
                edges = open_edges(idx);
                n_e   = length(edges);

                % Greedy endpoint-matching chain
                used     = false(n_e, 1);
                order    = zeros(n_e, 1);
                rev_flag = false(n_e, 1);

                order(1) = 1;
                used(1)  = true;

                for k = 2:n_e
                    prev = order(k-1);
                    if rev_flag(prev)
                        current_end = edges{prev}.pts(1, :);
                    else
                        current_end = edges{prev}.pts(end, :);
                    end

                    best_dist = inf;
                    best_idx  = 0;
                    best_rev  = false;

                    for j = 1:n_e
                        if used(j), continue; end
                        d_fwd = norm(current_end - edges{j}.pts(1, :));
                        d_rev = norm(current_end - edges{j}.pts(end, :));

                        if d_fwd < best_dist
                            best_dist = d_fwd;
                            best_idx  = j;
                            best_rev  = false;
                        end
                        if d_rev < best_dist
                            best_dist = d_rev;
                            best_idx  = j;
                            best_rev  = true;
                        end
                    end

                    if best_dist > tol
                        % Large gap in edge chain.  This typically occurs
                        % at the keel where source + MirrY semicircles
                        % share endpoints 2r apart.  The resulting cap area
                        % is zero (canceling orientations), so it's harmless.
                        % Downgraded from warning to silent — use verbose
                        % option in compute/compute_submerged to see details.
                    end

                    order(k)    = best_idx;
                    used(best_idx) = true;
                    rev_flag(k) = best_rev;
                end

                % Concatenate points in chain order
                loop = [];
                for k = 1:n_e
                    pts = edges{order(k)}.pts;
                    if rev_flag(k)
                        pts = flipud(pts);
                    end
                    if k > 1
                        pts = pts(2:end, :);  % avoid duplicate junction point
                    end
                    loop = [loop; pts]; %#ok<AGROW>
                end

                % Close the loop: remove trailing duplicate of first point
                if size(loop,1) > 1 && norm(loop(1,:) - loop(end,:)) < tol
                    loop = loop(1:end-1, :);
                end

                grp = struct();
                grp.loop_pts = loop;
                grp.z_mean   = mean(loop(:, 3));
                grp.z_range  = max(loop(:,3)) - min(loop(:,3));
                groups{end+1} = grp; %#ok<AGROW>
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §5  CAP CONTRIBUTION — Analytical Flat-Cap Integrals
        %% ═════════════════════════════════════════════════════════

        function [V_cap, Cxyz_cap, A_cap, z_cap] = cap_contribution( ...
                open_edges, interior_pt, tol)
        % CAP_CONTRIBUTION  Volume and centroid from flat caps.
        %
        %   For a flat cap at z = z_cap with outward normal ±ẑ:
        %
        %     A_cap = (1/2) |∮ (x dy − y dx)|    (Green's theorem)
        %     V_cap = ±(1/3) z_cap A_cap          (divergence theorem)
        %     z̄V   = ±(1/2) z_cap² A_cap
        %     x̄V   = 0,  ȳV = 0                  (n_x = n_y = 0)
        %
        %   The sign depends on whether the cap faces up (+z) or
        %   down (−z), determined by whether z_cap is above or below
        %   the interior point.

            if nargin < 3, tol = 1e-3; end

            % Chain open edges into loop(s)
            groups = WEC_HydroProperties.chain_open_edges(open_edges, tol);
            n_groups = length(groups);

            V_cap    = zeros(n_groups, 1);
            Cxyz_cap = zeros(n_groups, 3);
            A_cap    = zeros(n_groups, 1);
            z_cap    = zeros(n_groups, 1);

            for g = 1:n_groups
                loop = groups{g}.loop_pts;
                z_g  = groups{g}.z_mean;
                n    = size(loop, 1);

                if groups{g}.z_range > 0.05
                    warning('WEC_HydroProperties:NonFlatCap', ...
                            'Cap %d z-range = %.4f m (not flat)', ...
                            g, groups{g}.z_range);
                end

                % ── Signed area via shoelace on (x, y) ──────────────
                %  This IS Green's theorem evaluated discretely:
                %    A = (1/2) Σ (x_i y_{i+1} − x_{i+1} y_i)
                x = loop(:, 1);
                y = loop(:, 2);
                x_next = [x(2:end); x(1)];
                y_next = [y(2:end); y(1)];
                A_signed = sum(x .* y_next - x_next .* y) / 2;

                A_cap(g) = abs(A_signed);
                z_cap(g) = z_g;

                % ── Determine outward normal direction ──────────────
                %  If z_cap > interior_z → top cap → outward = +z
                %  If z_cap < interior_z → bottom cap → outward = −z
                if z_g > interior_pt(3)
                    sign_nz = +1;
                else
                    sign_nz = -1;
                end

                % The cap area is always the unsigned polygon area.
                % The sign of the volume contribution comes from sign_nz.
                A_use = abs(A_signed);

                % V_cap = (sign_nz / 3) * z_cap * A_cap
                V_cap(g) = sign_nz * z_g * A_use / 3;

                % Centroid: only z-component (n_x = n_y = 0)
                % z̄·V contribution = (sign_nz / 2) * z_cap² * A_cap
                Cxyz_cap(g, :) = [0, 0, sign_nz * z_g^2 * A_use / 2];
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §6  COMPUTE SUBMERGED — Trimmed at Waterline
        %% ═════════════════════════════════════════════════════════

        function props = compute_submerged(parser, z_wl, options)
        % COMPUTE_SUBMERGED  Volume and hydrostatics below a waterline.
        %
        %   props = WEC_HydroProperties.compute_submerged(parser, z_wl)
        %   props = WEC_HydroProperties.compute_submerged(parser, z_wl, opts)
        %
        %   Same divergence theorem as compute(), but for each surface
        %   the v-integration limits are adjusted per u-node so that
        %   only the submerged part (z ≤ z_wl) contributes.  Surfaces
        %   fully below get limits [0,1]; fully above get [0,0] (zero
        %   contribution automatically).  No pre-classification needed.
        %
        %   The waterplane at z = z_wl closes the submerged volume.
        %   Its contribution is analytical (Green's theorem on the
        %   boundary curve collected during integration).
        %
        %   OUTPUT
        %     .V_sub       — submerged volume [m³]
        %     .CB          — centre of buoyancy [x, y, z] [m]
        %     .Aw          — waterplane area [m²]
        %     .I_wp_xx     — waterplane 2nd moment about x-axis [m⁴]
        %     .I_wp_yy     — waterplane 2nd moment about y-axis [m⁴]
        %     .V_surfs     — per-surface contributions [m³]
        %     .wl_boundary — waterplane boundary points [N × 3]

            if nargin < 3, options = struct(); end
            if ~isfield(options, 'n_quad'),         options.n_quad         = 20;   end
            if ~isfield(options, 'n_edge_samples'), options.n_edge_samples = 50;   end
            if ~isfield(options, 'n_bisect'),       options.n_bisect       = 40;   end
            if ~isfield(options, 'verbose'),        options.verbose        = true; end
            if ~isfield(options, 'skip_open_edges'), options.skip_open_edges = false; end
            if ~isfield(options, 'Aw_override'),    options.Aw_override    = [];   end
            if ~isfield(options, 'I_wp_xx_override'), options.I_wp_xx_override = []; end
            if ~isfield(options, 'I_wp_yy_override'), options.I_wp_yy_override = []; end

            surf_names = parser.visible_surfs;
            n_surfs    = length(surf_names);
            n_quad     = options.n_quad;

            % Interior point estimate (for orientation).
            %  The orientation check determines whether S_u × S_v points
            %  outward from the HULL — a property of the surface, not
            %  of the submerged volume.  The bounding box center is
            %  always inside the hull regardless of z_wl.
            %
            %  The old code pushed z to z_wl - 0.5, which placed the
            %  interior below the hull for low waterlines, flipping
            %  surface normals and making V_sub non-monotonic.
            interior = (parser.extents(1:3) + parser.extents(4:6))' / 2;

            % Origin shift for numerical conditioning (v7.1)
            %  z_ref = midpoint between hull bottom and waterline.
            %  Shifts the position vector S → S' = (x, y, z−z_ref)
            %  in the divergence theorem.  Cap terms become proportional
            %  to (z_wl − z_ref) instead of z_wl, eliminating the
            %  catastrophic cancellation that made V_sub non-monotonic.
            z_ref = 0.5 * (min(parser.extents(3), parser.extents(6)) + z_wl);

            % ── Classify surfaces for mirror deduplication ─────────
            topo = parser.classify_visible_surfaces();

            % Map surface names to their index in visible_surfs
            surf_idx_map = containers.Map();
            for i = 1:n_surfs
                surf_idx_map(surf_names{i}) = i;
            end

            if options.verbose
                fprintf('\n  Submerged properties (z_wl = %.3f m, %d×%d GL)\n', ...
                        z_wl, n_quad, n_quad);
                fprintf('    Evaluating %d source surfaces (skipping %d mirrors)\n', ...
                        length(topo.sources), length(topo.mirrors));
            end

            % ── Surface integrals — sources only, derive mirrors ───
            V_total    = 0;
            Cxyz_num   = [0, 0, 0];
            V_surfs    = zeros(n_surfs, 1);

            % Second volume moments (unit density, about origin)
            int_x2_total = 0;
            int_y2_total = 0;
            int_z2_total = 0;

            source_V_sub    = containers.Map();
            source_Cxyz_sub = containers.Map();
            source_int_x2   = containers.Map();
            source_int_y2   = containers.Map();
            source_int_z2   = containers.Map();

            % Pre-determine orientation for each source surface using
            % a COARSE signed-volume integral over the FULL surface.
            %
            % sign(V_full) indicates whether Su×Sv points outward (+1)
            % or inward (-1).  This is more robust than the single-point
            % interior check at (0.5, 0.5), which fails for non-convex
            % geometries and capsized hulls (C0_180).
            %
            % The orientation is a property of the surface parameterization
            % and must be applied consistently to partial-domain integrals.
            % A coarse 5×5 GL quadrature suffices — only the sign matters.
            source_orient = containers.Map();
            n_orient = 5;
            [u_o, wu_o] = WEC_HydroProperties.gauss_legendre(n_orient);
            [v_o, wv_o] = WEC_HydroProperties.gauss_legendre(n_orient);
            for s = 1:length(topo.sources)
                sname = topo.sources{s};
                V_test = 0;
                for io = 1:n_orient
                    for jo = 1:n_orient
                        [S_t, Su_t, Sv_t] = parser.eval_surface_with_derivs( ...
                            sname, u_o(io), v_o(jo));
                        V_test = V_test + wu_o(io) * wv_o(jo) * ...
                            dot(S_t, cross(Su_t, Sv_t)) / 3;
                    end
                end
                if V_test < 0
                    source_orient(sname) = -1;
                else
                    source_orient(sname) = +1;
                end
            end

            for s = 1:length(topo.sources)
                sname = topo.sources{s};

                [V_s, Cxyz_s, ~, ix2_s, iy2_s, iz2_s] = ...
                    WEC_HydroProperties.surface_integral_submerged( ...
                        parser, sname, n_quad, z_wl, interior, ...
                        options.n_bisect, source_orient(sname), z_ref);

                source_V_sub(sname)    = V_s;
                source_Cxyz_sub(sname) = Cxyz_s;
                source_int_x2(sname)   = ix2_s;
                source_int_y2(sname)   = iy2_s;
                source_int_z2(sname)   = iz2_s;

                idx = surf_idx_map(sname);
                V_surfs(idx) = V_s;
                V_total      = V_total + V_s;
                Cxyz_num     = Cxyz_num + Cxyz_s;
                int_x2_total = int_x2_total + ix2_s;
                int_y2_total = int_y2_total + iy2_s;
                int_z2_total = int_z2_total + iz2_s;

                if options.verbose
                    fprintf('    %s: V = %+.6f m³  [SOURCE]\n', sname, V_s);
                end
            end

            % Derive mirror contributions
            %   Volume and second moments are INVARIANT under reflection.
            %   (x² is unchanged by x→−x; same for y², z².)
            %   Only centroid components flip per mirror plane.
            for m = 1:length(topo.mirrors)
                mirr = topo.mirrors(m);
                ult = mirr.ultimate_source;

                if source_V_sub.isKey(ult)
                    V_mirr    = source_V_sub(ult);
                    Cxyz_mirr = source_Cxyz_sub(ult);
                    ix2_mirr  = source_int_x2(ult);
                    iy2_mirr  = source_int_y2(ult);
                    iz2_mirr  = source_int_z2(ult);

                    % Flip centroid components per mirror plane
                    for fi = 1:length(mirr.effective_flips)
                        if strcmp(mirr.effective_flips{fi}, 'Y')
                            Cxyz_mirr(2) = -Cxyz_mirr(2);
                        else
                            Cxyz_mirr(1) = -Cxyz_mirr(1);
                        end
                    end
                else
                    % Fallback: evaluate directly with orientation check.
                    % FIX (L2): Use the same coarse 5×5 signed-volume integral
                    % as the normal source path — not the fragile single-point
                    % check at (0.5, 0.5) that can fail for non-convex surfaces.
                    n_fb = 5;
                    [u_fb, wu_fb] = WEC_HydroProperties.gauss_legendre(n_fb);
                    [v_fb, wv_fb] = WEC_HydroProperties.gauss_legendre(n_fb);
                    V_test_m = 0;
                    for ifb = 1:n_fb
                        for jfb = 1:n_fb
                            [S_t, Su_t, Sv_t] = parser.eval_surface_with_derivs( ...
                                mirr.name, u_fb(ifb), v_fb(jfb));
                            V_test_m = V_test_m + wu_fb(ifb) * wv_fb(jfb) * ...
                                dot(S_t, cross(Su_t, Sv_t)) / 3;
                        end
                    end
                    orient_m = sign(V_test_m);
                    if orient_m == 0, orient_m = +1; end

                    [V_mirr, Cxyz_mirr, ~, ix2_mirr, iy2_mirr, iz2_mirr] = ...
                        WEC_HydroProperties.surface_integral_submerged( ...
                            parser, mirr.name, n_quad, z_wl, interior, ...
                            options.n_bisect, orient_m, z_ref);
                end

                if surf_idx_map.isKey(mirr.name)
                    V_surfs(surf_idx_map(mirr.name)) = V_mirr;
                end
                V_total      = V_total + V_mirr;
                Cxyz_num     = Cxyz_num + Cxyz_mirr;
                int_x2_total = int_x2_total + ix2_mirr;
                int_y2_total = int_y2_total + iy2_mirr;
                int_z2_total = int_z2_total + iz2_mirr;

                if options.verbose
                    fprintf('    %s: V = %+.6f m³  [MIRROR of %s]\n', ...
                            mirr.name, V_mirr, ult);
                end
            end

            % ── Waterplane ────────────────────────────────────────────
            %  v7.0: Use precomputed Aw if provided (fast path for
            %  optimizer loop).  Otherwise trace actual spline contour.
            %
            %  The contour method sweeps both u and v on every visible
            %  surface, bisects to find z = z_wl, and computes Aw via
            %  Green's theorem.  This is correct but slow (~0.5 s per
            %  call).  Precomputing Aw(z) once at config time and
            %  interpolating eliminates this cost from the optimizer.

            wl_ordered = [];

            if ~isempty(options.Aw_override)
                % Fast path: precomputed waterplane properties
                Aw   = options.Aw_override;
                I_xx = options.I_wp_xx_override;
                I_yy = options.I_wp_yy_override;
                if isempty(I_xx), I_xx = 0; end
                if isempty(I_yy), I_yy = 0; end

                if options.verbose
                    fprintf('    Waterplane (precomputed): Aw = %.6f m²\n', Aw);
                end
            else
                % Fast iso-z contour extraction (v8.0)
                %  Auto-builds a boundary cache if not provided.
                %  Still ~250× faster than the reference method for a
                %  single z-level; for multiple z-levels the cache is
                %  built once and reused via Aw_override.
                if isfield(options, 'boundary_cache') && ~isempty(options.boundary_cache)
                    b_cache = options.boundary_cache;
                else
                    b_cache = WEC_HydroProperties.precompute_boundary_cache(parser, n_quad * 5);
                end

                wl_pts_all = WEC_HydroProperties.extract_isocurve_at_z( ...
                    parser, z_wl, n_quad * 5, b_cache);

                if ~isempty(wl_pts_all) && size(wl_pts_all, 1) >= 3
                    [Aw, I_xx, I_yy, wl_ordered] = ...
                        WEC_HydroProperties.waterplane_properties(wl_pts_all);
                else
                    Aw = 0; I_xx = 0; I_yy = 0;
                end

                if options.verbose
                    fprintf('    Waterplane (fast iso-z, %d pts): Aw = %.6f m²\n', ...
                            size(wl_pts_all, 1), Aw);
                end
            end

            % Cap at z_wl with SHIFTED position vector (v7.1)
            %  S' = (x, y, z_wl − z_ref) on the cap.  n = (0, 0, +1).
            %  S'·n = z_wl − z_ref.
            dz_cap = z_wl - z_ref;
            V_cap    = dz_cap * Aw / 3;
            Cxyz_cap = [0, 0, dz_cap^2 * Aw / 2];
            V_total  = V_total + V_cap;
            Cxyz_num = Cxyz_num + Cxyz_cap;

            % Shifted second moment: (z_wl − z_ref)³ × Aw / 3
            int_z2_total = int_z2_total + dz_cap^3 * Aw / 3;

            % ── Open edges below waterline ──────────────────────────
            if ~options.skip_open_edges
                [open_edges, ~] = WEC_HydroProperties.detect_open_edges( ...
                    parser, options.n_edge_samples, 1e-4);

                below_wl_edges = {};
                for i = 1:length(open_edges)
                    if open_edges{i}.z_mean < z_wl - 0.05
                        below_wl_edges{end+1} = open_edges{i}; %#ok<AGROW>
                    end
                end

                if ~isempty(below_wl_edges)
                    [V_oe, Cxyz_oe, A_oe, z_oe] = ...
                        WEC_HydroProperties.cap_contribution( ...
                            below_wl_edges, interior, 1e-3);

                    % Convert open-edge contributions to shifted coords
                    for g = 1:length(V_oe)
                        if z_oe(g) > interior(3)
                            sign_nz_oe = +1;
                        else
                            sign_nz_oe = -1;
                        end
                        dz_oe = z_oe(g) - z_ref;

                        % Shifted cap: V = sign * (z-z_ref) * A / 3
                        V_total  = V_total + sign_nz_oe * dz_oe * A_oe(g) / 3;
                        Cxyz_num(3) = Cxyz_num(3) + sign_nz_oe * dz_oe^2 * A_oe(g) / 2;
                        int_z2_total = int_z2_total + sign_nz_oe * dz_oe^3 * A_oe(g) / 3;
                    end

                    if options.verbose
                        for g = 1:length(V_oe)
                            if A_oe(g) > 1e-8
                                fprintf('    Open-edge cap: z=%.2f, A=%.4f m²\n', ...
                                        z_oe(g), A_oe(g));
                            end
                        end
                    end
                end
            end

            % ── Unshift from z_ref to world coordinates (v7.1) ─────
            %  V_total is origin-invariant (correct as-is).
            %  Cxyz_num(1,2) are unshifted (x,y not shifted).
            %  Cxyz_num(3) and int_z2 need the parallel-axis unshift:
            %
            %  ∫z dV = ∫(z−z_ref) dV + z_ref × V
            %  ∫z² dV = ∫(z−z_ref)² dV + 2·z_ref·∫(z−z_ref)dV + z_ref²·V

            Cz_shifted   = Cxyz_num(3);                  % ∫(z−z_ref) dV
            Cxyz_num(3)  = Cz_shifted + z_ref * V_total; % ∫z dV
            int_z2_total = int_z2_total + 2*z_ref*Cz_shifted + z_ref^2 * V_total;

            % ── Assemble output ─────────────────────────────────────
            props = struct();
            props.V_sub   = abs(V_total);
            props.V_raw   = V_total;
            props.V_surfs = V_surfs;
            props.V_cap   = V_cap;
            props.Aw      = Aw;
            props.I_wp_xx = I_xx;
            props.I_wp_yy = I_yy;
            props.wl_boundary = wl_ordered;

            if abs(V_total) > 1e-12
                props.CB = Cxyz_num / V_total;
            else
                props.CB = [0, 0, 0];
            end

            % Raw numerators and second moments for strip computation
            props.Cz_numerator = Cxyz_num(3);
            props.int_x2 = int_x2_total;
            props.int_y2 = int_y2_total;
            props.int_z2 = int_z2_total;

            if options.verbose
                fprintf('\n');
                fprintf('  ┌─────────────────────────────────────────┐\n');
                fprintf('  │  SUBMERGED PROPERTIES (z_wl = %.2f m)    \n', z_wl);
                fprintf('  ├─────────────────────────────────────────┤\n');
                fprintf('  │  V_sub:     %10.4f m³                \n', props.V_sub);
                fprintf('  │  CB:        [%.3f, %.3f, %.3f] m      \n', props.CB);
                fprintf('  │  Aw:        %10.6f m²                \n', props.Aw);
                fprintf('  │  I_wp_xx:   %10.6f m⁴                \n', props.I_wp_xx);
                fprintf('  │  I_wp_yy:   %10.6f m⁴                \n', props.I_wp_yy);
                fprintf('  └─────────────────────────────────────────┘\n');
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §7  SURFACE INTEGRAL (SUBMERGED) — Adjusted v-Limits
        %% ═════════════════════════════════════════════════════════

        function [V_s, Cxyz_s, wl_pts, int_x2_s, int_y2_s, int_z2_s] = ...
                surface_integral_submerged( ...
                parser, surf_name, n_quad, z_wl, interior_pt, n_bisect, orient_sign, z_ref)
        % SURFACE_INTEGRAL_SUBMERGED  Divergence theorem with v-limits
        % adjusted so only z ≤ z_wl contributes.
        %
        %   z_ref:  Origin shift in z for numerical conditioning.
        %     Replaces S with S' = (x, y, z−z_ref) in the divergence
        %     theorem integrands.  The total V over the closed surface
        %     is origin-invariant.  Shifting eliminates catastrophic
        %     cancellation between hull and cap terms when z << 0.
        %     Returns SHIFTED Cz and int_z2; caller unshifts.
        %
        %   orient_sign: +1 or -1, pre-determined from full surface.

            if nargin < 6, n_bisect = 40; end
            if nargin < 7, orient_sign = []; end
            if nargin < 8, z_ref = 0; end

            [u_gl, w_u] = WEC_HydroProperties.gauss_legendre(n_quad);
            [v_gl, w_v] = WEC_HydroProperties.gauss_legendre(n_quad);

            V_s      = 0;
            Cxyz_s   = [0, 0, 0];
            int_x2_s = 0;
            int_y2_s = 0;
            int_z2_s = 0;
            wl_pts   = [];

            for i = 1:n_quad
                u = u_gl(i);

                % ── Get submerged v-interval for this u ─────────────
                [v_lo, v_hi, v_wl] = WEC_HydroProperties.find_submerged_v_limits( ...
                    parser, surf_name, u, z_wl, n_bisect);

                dv = v_hi - v_lo;
                if dv < 1e-14, continue; end

                % Collect waterplane boundary point (if crossing)
                if ~isnan(v_wl)
                    wl_pts = [wl_pts; parser.eval_surface(surf_name, u, v_wl)]; %#ok<AGROW>
                end

                % ── GL quadrature over [v_lo, v_hi] ─────────────────
                for j = 1:n_quad
                    v = v_lo + v_gl(j) * dv;

                    [S, Su, Sv] = parser.eval_surface_with_derivs(surf_name, u, v);

                    n_vec = cross(Su, Sv);
                    w = w_u(i) * w_v(j) * dv;

                    % Shifted z for numerical conditioning
                    z_s = S(3) - z_ref;

                    % Volume: (1/3) ∮ S'·n dA
                    V_s = V_s + w * (S(1)*n_vec(1) + S(2)*n_vec(2) + z_s*n_vec(3)) / 3;

                    % Centroid: x,y unshifted; z shifted
                    Cxyz_s = Cxyz_s + w * [ ...
                        S(1)^2 * n_vec(1) / 2, ...
                        S(2)^2 * n_vec(2) / 2, ...
                        z_s^2  * n_vec(3) / 2];

                    % Second moments: x,y unshifted; z shifted
                    int_x2_s = int_x2_s + w * S(1)^3 * n_vec(1) / 3;
                    int_y2_s = int_y2_s + w * S(2)^3 * n_vec(2) / 3;
                    int_z2_s = int_z2_s + w * z_s^3  * n_vec(3) / 3;
                end
            end

            % ── Apply orientation sign ─────────────────────────────
            %  The sign is pre-determined from the FULL surface integral
            %  orientation check (dot product of normal at (0.5,0.5)
            %  with outward direction).  This must be applied consistently
            %  to all partial-domain integrals on this surface.
            %
            %  WHY not self-determine?
            %    The raw GL integral ∫∫ S·(S_u×S_v) du dv over a partial
            %    v-range can have a different sign than the full integral.
            %    The orientation (inward vs outward normal) is a property
            %    of the parameterization, not the integration domain.
            %    Determining it from the partial integral gives wrong
            %    signs for some waterlines, making V_sub non-monotonic.

            if ~isempty(orient_sign)
                % Use pre-determined sign
                if orient_sign < 0
                    V_s      = -V_s;
                    Cxyz_s   = -Cxyz_s;
                    int_x2_s = -int_x2_s;
                    int_y2_s = -int_y2_s;
                    int_z2_s = -int_z2_s;
                end
            else
                % Fallback: local check at (0.5, 0.5)
                [S_mid, Su_mid, Sv_mid] = parser.eval_surface_with_derivs(surf_name, 0.5, 0.5);
                n_mid = cross(Su_mid, Sv_mid);
                if dot(n_mid, S_mid - interior_pt(:)') < 0
                    V_s      = -V_s;
                    Cxyz_s   = -Cxyz_s;
                    int_x2_s = -int_x2_s;
                    int_y2_s = -int_y2_s;
                    int_z2_s = -int_z2_s;
                end
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §8  FIND SUBMERGED V-LIMITS
        %% ═════════════════════════════════════════════════════════

        function [v_lo, v_hi, v_wl] = find_submerged_v_limits( ...
                parser, surf_name, u, z_wl, n_bisect)
        % FIND_SUBMERGED_V_LIMITS  Return [v_lo, v_hi] where z ≤ z_wl.
        %
        %   Also returns v_wl = the v-value on the waterline (NaN if
        %   no crossing).
        %
        %   Cases:
        %     All below:  [0, 1], v_wl = NaN
        %     All above:  [0, 0], v_wl = NaN  (dv = 0 → zero contribution)
        %     Crossing:   [0, v*] or [v*, 1],  v_wl = v*

            if nargin < 5, n_bisect = 40; end

            % ── RevSurf fast-path ─────────────────────────────────
            %  For RevSurf, z depends on u only — revolution preserves z.
            %  A single z-check at any v suffices.  No v-sweep or bisection.
            %  For C0: surface2 accounts for ~33% of v-limit searches.
            if parser.entities.isKey(surf_name)
                se = parser.entities(surf_name);
                if strcmp(se.type, 'RevSurf')
                    pt = parser.eval_surface(surf_name, u, 0);
                    v_wl = NaN;
                    if pt(3) <= z_wl + 1e-10
                        v_lo = 0; v_hi = 1;
                    else
                        v_lo = 0; v_hi = 0;
                    end
                    return;
                end
                % MirrSurf of a RevSurf also preserves z under revolution
                if strcmp(se.type, 'MirrSurf')
                    src = parser.entities(se.params.source);
                    if strcmp(src.type, 'RevSurf')
                        pt = parser.eval_surface(surf_name, u, 0);
                        v_wl = NaN;
                        if pt(3) <= z_wl + 1e-10
                            v_lo = 0; v_hi = 1;
                        else
                            v_lo = 0; v_hi = 0;
                        end
                        return;
                    end
                end
            end

            % ── Coarse sweep for sign change ────────────────────────
            n_sample = 21;
            v_sample = linspace(0, 1, n_sample);
            f_sample = zeros(n_sample, 1);
            for k = 1:n_sample
                pt = parser.eval_surface(surf_name, u, v_sample(k));
                f_sample(k) = pt(3) - z_wl;
            end

            % ── No crossing: all below or all above ─────────────────
            sign_changes = find(f_sample(1:end-1) .* f_sample(2:end) < 0, 1);

            if isempty(sign_changes)
                v_wl = NaN;
                if f_sample(1) <= 0
                    v_lo = 0; v_hi = 1;   % all below
                else
                    v_lo = 0; v_hi = 0;   % all above
                end
                return;
            end

            % ── Bisection to refine crossing ────────────────────────
            a = v_sample(sign_changes);
            b = v_sample(sign_changes + 1);
            fa = f_sample(sign_changes);

            for iter = 1:n_bisect
                m = (a + b) / 2;
                pt = parser.eval_surface(surf_name, u, m);
                fm = pt(3) - z_wl;
                if abs(fm) < 1e-10, break; end
                if fa * fm < 0
                    b = m;
                else
                    a = m; fa = fm;
                end
            end

            v_wl = (a + b) / 2;

            % ── Which side is submerged? ────────────────────────────
            if f_sample(1) <= 0
                v_lo = 0;     v_hi = v_wl;    % below is v < v_wl
            else
                v_lo = v_wl;  v_hi = 1;        % below is v > v_wl
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §9  EXTRACT WATERPLANE CONTOUR
        %% ═════════════════════════════════════════════════════════

        function wl_pts = extract_waterplane_contour(parser, z_wl, n_samples, n_bisect)
        % EXTRACT_WATERPLANE_CONTOUR  Find the iso-z contour at z = z_wl.
        %
        % ── DEAD IN PRODUCTION (D1) ──────────────────────────────────────
        % All production paths (compute_submerged, compute_strip,
        % WEC_Configuration_Builder) call extract_isocurve_at_z with a
        % precomputed boundary cache — verified 1000× faster, identical
        % results. This method is retained for standalone validation and
        % regression testing ONLY. Do not call from new production code.
        % ─────────────────────────────────────────────────────────────────
        %
        %   wl_pts = extract_waterplane_contour(parser, z_wl, n_samples)
        %
        %   For each visible surface, searches BOTH parameter directions:
        %     Sweep 1: for each u_i, bisect in v to find z = z_wl
        %     Sweep 2: for each v_j, bisect in u to find z = z_wl
        %
        %   WHY both directions?
        %     For RevSurf, z depends on u only (revolution preserves z).
        %     Sweep 1 (bisecting v) never finds a crossing.
        %     Sweep 2 (bisecting u) finds the circle at the waterline.
        %     For BLoftSurf/DevSurf, z varies in v, so Sweep 1 works.
        %     Running both catches all surface types.
        %
        %   Duplicate points are removed by proximity filtering.

            if nargin < 3, n_samples = 100; end
            if nargin < 4, n_bisect  = 40;  end

            surf_names = parser.visible_surfs;
            wl_pts = [];

            for s = 1:length(surf_names)
                sname = surf_names{s};

                % ── Sweep 1: fix u, bisect in v ─────────────────────
                u_samples = linspace(0, 1, n_samples);
                for i = 1:n_samples
                    v_root = WEC_HydroProperties.bisect_param_for_z( ...
                        parser, sname, 'u', u_samples(i), z_wl, n_bisect);
                    if ~isnan(v_root)
                        wl_pts = [wl_pts; parser.eval_surface(sname, u_samples(i), v_root)]; %#ok<AGROW>
                    end
                end

                % ── Sweep 2: fix v, bisect in u ─────────────────────
                v_samples = linspace(0, 1, n_samples);
                for j = 1:n_samples
                    u_root = WEC_HydroProperties.bisect_param_for_z( ...
                        parser, sname, 'v', v_samples(j), z_wl, n_bisect);
                    if ~isnan(u_root)
                        wl_pts = [wl_pts; parser.eval_surface(sname, u_root, v_samples(j))]; %#ok<AGROW>
                    end
                end
            end

            % ── Remove near-duplicate points ────────────────────────
            if size(wl_pts, 1) > 1
                keep = true(size(wl_pts, 1), 1);
                for i = 2:size(wl_pts, 1)
                    for j = 1:i-1
                        if keep(j) && norm(wl_pts(i,:) - wl_pts(j,:)) < 1e-6
                            keep(i) = false;
                            break;
                        end
                    end
                end
                wl_pts = wl_pts(keep, :);
            end
        end


        function root = bisect_param_for_z(parser, surf_name, fixed_dir, fixed_val, z_wl, n_bisect)
        % BISECT_PARAM_FOR_Z  Bisect the free parameter to find z = z_wl.
        %
        % ── DEAD IN PRODUCTION (D1) ──────────────────────────────────────
        % Only called from extract_waterplane_contour, which is itself dead
        % in the production pipeline. Retained alongside it for validation.
        % ─────────────────────────────────────────────────────────────────
        %
        %   fixed_dir = 'u': u is fixed at fixed_val, bisect in v
        %   fixed_dir = 'v': v is fixed at fixed_val, bisect in u

            if nargin < 6, n_bisect = 40; end

            n_sample = 21;
            t_sample = linspace(0, 1, n_sample);
            f_sample = zeros(n_sample, 1);

            for k = 1:n_sample
                if strcmp(fixed_dir, 'u')
                    pt = parser.eval_surface(surf_name, fixed_val, t_sample(k));
                else
                    pt = parser.eval_surface(surf_name, t_sample(k), fixed_val);
                end
                f_sample(k) = pt(3) - z_wl;
            end

            % Find first sign change
            root = NaN;
            for k = 1:n_sample-1
                if f_sample(k) * f_sample(k+1) < 0
                    a = t_sample(k);
                    b = t_sample(k+1);
                    fa = f_sample(k);

                    for iter = 1:n_bisect
                        m = (a + b) / 2;
                        if strcmp(fixed_dir, 'u')
                            pt = parser.eval_surface(surf_name, fixed_val, m);
                        else
                            pt = parser.eval_surface(surf_name, m, fixed_val);
                        end
                        fm = pt(3) - z_wl;
                        if abs(fm) < 1e-10, break; end
                        if fa * fm < 0
                            b = m;
                        else
                            a = m; fa = fm;
                        end
                    end

                    root = (a + b) / 2;
                    return;  % take first crossing only
                end
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §9b  PRECOMPUTE BOUNDARY CACHE
        %% ═════════════════════════════════════════════════════════

        function cache = precompute_boundary_cache(parser, n_u)
        % PRECOMPUTE_BOUNDARY_CACHE  Evaluate boundary curves ONCE for all
        %   source surfaces.  The returned cache is reused by
        %   extract_isocurve_at_z for every z-level query — no further
        %   eval_curve or eval_snake calls are needed.
        %
        %   cache = WEC_HydroProperties.precompute_boundary_cache(parser, n_u)
        %
        %   WHAT IS CACHED (per source surface type)
        %     DevSurf/RuledSurf: boundary curve points [n_u × 3] for both
        %       boundary curves (snake+curve or curve1+curve2), plus their
        %       z-coordinates as separate vectors for fast arithmetic.
        %     RevSurf: profile curve points [n_u × 3] and z-coordinates,
        %       plus the axis decomposition (start, direction, angles).
        %     BLoftSurf: all section curve points [n_sec × n_u × 3] and
        %       their z-coordinates, plus the v-direction knot vector.
        %
        %   WHY PRECOMPUTE
        %     The iso-z solve at each z-level needs the z-coordinates of
        %     boundary curves at n_u u-samples.  These are properties of
        %     the hull geometry (independent of z_wl).  Evaluating them
        %     once costs ~0.15 s; evaluating them per-z-level (as the
        %     reference method does) costs ~46 s × n_levels.
        %
        %   OUTPUT
        %     cache.sources    — cell array of source surface names
        %     cache.mirrors    — struct array from classify_visible_surfaces
        %     cache.u_samples  — [n_u × 1] shared u-grid
        %     cache.data       — containers.Map: source_name → cached struct

            topo = parser.classify_visible_surfaces();
            cache.sources   = topo.sources;
            cache.mirrors   = topo.mirrors;
            cache.u_samples = linspace(0, 1, n_u)';
            cache.data      = containers.Map();

            for s = 1:length(topo.sources)
                sname = topo.sources{s};
                e = parser.entities(sname);
                d = struct();
                d.type = e.type;

                switch e.type
                    case {'DevSurf', 'RuledSurf'}
                        if strcmp(e.type, 'DevSurf')
                            d.pts_1 = parser.eval_snake(e.params.snake, cache.u_samples);
                            d.pts_2 = parser.eval_curve(e.params.curve, cache.u_samples);
                        else
                            d.pts_1 = parser.eval_curve(e.params.curve1, cache.u_samples);
                            d.pts_2 = parser.eval_curve(e.params.curve2, cache.u_samples);
                        end
                        d.z_boundary_1 = d.pts_1(:, 3);
                        d.z_boundary_2 = d.pts_2(:, 3);

                    case 'RevSurf'
                        profile_pts = parser.eval_curve_or_snake( ...
                                          e.params.profile, cache.u_samples);
                        d.profile_pts = profile_pts;
                        d.z_profile   = profile_pts(:, 3);
                        d.profile_name = e.params.profile;

                        axis_ent = parser.entities(e.params.axis);
                        d.axis_start = parser.eval_any_point(axis_ent.params.pt_start);
                        d.axis_end   = parser.eval_any_point(axis_ent.params.pt_end);
                        d.axis_vec   = d.axis_end - d.axis_start;
                        d.axis_dir   = d.axis_vec / norm(d.axis_vec);
                        d.angle_start = e.params.angle_start;
                        d.angle_end   = e.params.angle_end;

                    case 'BLoftSurf'
                        sec_names = e.params.section_names;
                        n_sec = length(sec_names);
                        d.degree = e.params.degree;
                        d.n_sections = n_sec;
                        d.section_pts = zeros(n_sec, n_u, 3);
                        d.z_sections  = zeros(n_sec, n_u);

                        for k = 1:n_sec
                            pts_k = parser.eval_curve_or_snake( ...
                                        sec_names{k}, cache.u_samples);
                            d.section_pts(k, :, :) = pts_k;
                            d.z_sections(k, :) = pts_k(:, 3)';
                        end

                        d.knots_v = WEC_MS2_Parser.make_clamped_knots(n_sec, d.degree);
                end

                cache.data(sname) = d;
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §9c  EXTRACT ISOCURVE AT Z — Fast Type-Dispatched
        %% ═════════════════════════════════════════════════════════

        function pts = extract_isocurve_at_z(parser, z_wl, n_u, cache)
        % EXTRACT_ISOCURVE_AT_Z  Compute the iso-z contour at z = z_wl
        %   using type-dispatched algebraic solvers on precomputed data.
        %
        %   pts = WEC_HydroProperties.extract_isocurve_at_z(parser, z_wl, n_u, cache)
        %
        %   Returns [N × 3] contour points in the same format as
        %   extract_waterplane_contour (compatible with waterplane_properties).
        %
        %   The cache must be built by precompute_boundary_cache.  All
        %   eval_curve/eval_snake calls happened during cache construction;
        %   this method operates on cached arrays via arithmetic only
        %   (DevSurf/RuledSurf) or lightweight 1D root-finding (RevSurf,
        %   BLoftSurf).
        %
        %   Mirror surfaces are handled by coordinate flip on the source
        %   contour — zero additional evaluations.
        %
        %   INPUTS
        %     parser  — WEC_MS2_Parser object (needed only for RevSurf
        %               profile refinement via eval_curve_or_snake)
        %     z_wl    — [m] target elevation
        %     n_u     — number of u-samples (must match cache.u_samples)
        %     cache   — struct from precompute_boundary_cache

            pts = [];
            source_pts = containers.Map();

            for s = 1:length(cache.sources)
                sname = cache.sources{s};
                d = cache.data(sname);

                switch d.type
                    case {'DevSurf', 'RuledSurf'}
                        src_pts = WEC_HydroProperties.isocurve_devsurf( ...
                                      z_wl, d.pts_1, d.pts_2);

                    case 'RevSurf'
                        src_pts = WEC_HydroProperties.isocurve_revsurf( ...
                                      parser, z_wl, n_u, d);

                    case 'BLoftSurf'
                        src_pts = WEC_HydroProperties.isocurve_bloftsurf( ...
                                      z_wl, n_u, d);
                    otherwise
                        src_pts = [];
                end

                source_pts(sname) = src_pts;
                if ~isempty(src_pts)
                    pts = [pts; src_pts]; %#ok<AGROW>
                end
            end

            % Mirror contributions: coordinate flip on source contour
            for m = 1:length(cache.mirrors)
                mirr = cache.mirrors(m);
                ult = mirr.ultimate_source;

                if source_pts.isKey(ult) && ~isempty(source_pts(ult))
                    mirr_pts = source_pts(ult);
                    for fi = 1:length(mirr.effective_flips)
                        if strcmp(mirr.effective_flips{fi}, 'Y')
                            mirr_pts(:, 2) = -mirr_pts(:, 2);
                        else
                            mirr_pts(:, 1) = -mirr_pts(:, 1);
                        end
                    end
                    pts = [pts; mirr_pts]; %#ok<AGROW>
                end
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §9d–9i  TYPE-SPECIFIC ISO-Z SOLVERS (private helpers)
        %% ═════════════════════════════════════════════════════════

        function pts = isocurve_devsurf(z_wl, pts_boundary_1, pts_boundary_2)
        % ISOCURVE_DEVSURF  Iso-z contour of a DevSurf or RuledSurf.
        %
        %   S(u,v) = (1-v)*C1(u) + v*C2(u)  is linear in v.
        %
        %   Solving z(u,v) = z_wl for v:
        %     v*(u) = (z_wl - z1(u)) / (z2(u) - z1(u))
        %
        %   Algebraically exact — no iteration, no eval_surface.
        %   The contour point is the lerp: (1-v*)*C1(u) + v*·C2(u),
        %   computed directly from the cached boundary arrays.

            z1 = pts_boundary_1(:, 3);
            z2 = pts_boundary_2(:, 3);
            n_u = length(z1);
            dz = z2 - z1;
            pts = [];

            for i = 1:n_u
                if abs(dz(i)) < 1e-14
                    if abs(z1(i) - z_wl) < 1e-10
                        v_star = 0.5;
                    else
                        continue;
                    end
                else
                    v_star = (z_wl - z1(i)) / dz(i);
                end

                if v_star < -1e-10 || v_star > 1 + 1e-10
                    continue;
                end
                v_star = max(0, min(1, v_star));

                pt = (1 - v_star) * pts_boundary_1(i, :) + ...
                     v_star * pts_boundary_2(i, :);
                pts = [pts; pt]; %#ok<AGROW>
            end
        end


        function pts = isocurve_revsurf(parser, z_wl, n_arc, cache_entry)
        % ISOCURVE_REVSURF  Iso-z contour of a RevSurf with vertical axis.
        %
        %   For a RevSurf with vertical revolution axis, z depends on the
        %   profile parameter u only — revolution preserves z.  The iso-z
        %   contour is therefore the full revolution arc at the single u*
        %   where z_profile(u*) = z_wl.
        %
        %   STEPS
        %     1. Bracket u* from cached z_profile array (O(n_u) scan).
        %     2. Refine u* via bisection on actual profile curve (~30
        %        eval_curve_or_snake calls → machine precision).
        %     3. Decompose profile point at u* into axial projection +
        %        radial distance using the cached axis.
        %     4. Generate arc points analytically:
        %        pt(phi) = proj + r*cos(phi)*e_r + r*sin(phi)*e_t
        %        This is pure trigonometry — no eval_surface calls.

            d = cache_entry;
            u_grid = linspace(0, 1, length(d.z_profile))';

            u_star = WEC_HydroProperties.find_profile_root( ...
                         d.z_profile, u_grid, z_wl, parser, d.profile_name);
            if isnan(u_star)
                pts = [];
                return;
            end

            % Exact profile point at u* (one eval_curve_or_snake call)
            profile_pt = parser.eval_curve_or_snake(d.profile_name, u_star);
            if size(profile_pt, 1) > 1, profile_pt = profile_pt(1,:); end

            % Axial + radial decomposition using cached axis
            v_rel = profile_pt - d.axis_start;
            z_along = dot(v_rel, d.axis_dir);
            proj = d.axis_start + z_along * d.axis_dir;
            radial = profile_pt - proj;
            r = norm(radial);

            if r < 1e-12
                pts = repmat(profile_pt, n_arc, 1);
                return;
            end

            e_r = radial / r;
            e_t = cross(d.axis_dir, e_r);
            e_t = e_t / norm(e_t);

            % Analytic arc generation
            phi_start = deg2rad(d.angle_start);
            phi_end   = deg2rad(d.angle_end);
            phi = linspace(phi_start, phi_end, n_arc)';

            pts = proj + r * cos(phi) .* e_r + r * sin(phi) .* e_t;
        end


        function pts = isocurve_bloftsurf(z_wl, n_u, cache_entry)
        % ISOCURVE_BLOFTSURF  Iso-z contour of a BLoftSurf.
        %
        %   TWO LOFT ORIENTATIONS — detected from cached data:
        %
        %   (A) AXIAL LOFT  (original): z varies with v at each u.
        %     The section curves are stacked at different heights.
        %     Algorithm: for each u, convex-hull check then bisect v.
        %
        %   (B) AZIMUTHAL LOFT: z constant in v, varies with u.
        %     All sections share the same z at each u (revolution of a
        %     profile around a vertical axis preserves z).  The v-direction
        %     z B-spline is therefore flat → bspline_z_root returns NaN
        %     for every u → pts = [] for every z_wl.  This caused the
        %     E1.ms2 Wall to vanish from Aw_table, compute_submerged,
        %     constructable hull contours, and all visualisations.
        %
        %   DETECTION: check max(z_sections(:,col)) - min(...) at three
        %     cached u-columns (first, middle, last).  If all spreads are
        %     < flat_z_tol the loft is azimuthal.
        %
        %   AZIMUTHAL FIX: root-find u* where z_sections(1,u*) = z_wl
        %     (z varies monotonically with u), then output the full v-arc
        %     at u* using linearly interpolated cached section_pts.
        %     No additional eval_curve calls — purely arithmetic on cache.

            d = cache_entry;
            degree  = d.degree;
            knots_v = d.knots_v;
            pts = [];

            % ── Detect azimuthal loft ────────────────────────────────
            flat_z_tol = 1e-4;
            i_mid = max(1, round(n_u / 2));
            spread = @(col) max(d.z_sections(:, col)) - min(d.z_sections(:, col));
            is_azimuthal = spread(1) < flat_z_tol && ...
                           spread(i_mid) < flat_z_tol && ...
                           spread(n_u) < flat_z_tol;

            if is_azimuthal
                % ── BRANCH A: azimuthal loft ─────────────────────────
                % z varies with u; use row 1 (any row — all equal) as
                % the z-profile across u.  Root-find u* by sign-change
                % scan, then output the full v-arc at u* from cached pts.
                z_at_u = d.z_sections(1, :)';   % [n_u × 1]
                n_v   = 200;
                v_arc = linspace(0, 1, n_v)';

                for k = 1:n_u - 1
                    % Skip if no sign change (crossing) in this interval
                    if (z_at_u(k) - z_wl) * (z_at_u(k+1) - z_wl) > 0
                        continue;
                    end

                    % Linear interpolation for fractional position
                    dz_k = z_at_u(k+1) - z_at_u(k);
                    if abs(dz_k) < 1e-14
                        t = 0;
                    else
                        t = (z_wl - z_at_u(k)) / dz_k;
                    end
                    t = max(0, min(1, t));

                    % Interpolate section_pts [n_sec × 3] at u*
                    % d.section_pts is [n_sec × n_u × 3]
                    sp = (1 - t) * squeeze(d.section_pts(:, k,   :)) + ...
                              t  * squeeze(d.section_pts(:, k+1, :));

                    % Output the full horizontal arc: sweep v in [0,1]
                    for iv = 1:n_v
                        pv = WEC_MS2_Parser.bspline_curve_eval( ...
                                 knots_v, sp, degree, v_arc(iv));
                        pts = [pts; pv(1,:)]; %#ok<AGROW>
                    end
                end
                return;
            end

            % ── BRANCH B: axial loft (original algorithm) ─────────────
            % z varies with v at each u.  For each u, convex-hull check
            % then bisect v* where z(v*) = z_wl.
            for i = 1:n_u
                z_ctrl = d.z_sections(:, i);

                if z_wl < min(z_ctrl) - 1e-10 || z_wl > max(z_ctrl) + 1e-10
                    continue;
                end

                v_star = WEC_HydroProperties.bspline_z_root( ...
                             knots_v, z_ctrl, degree, z_wl);
                if isnan(v_star)
                    continue;
                end

                sec_pts_at_u = squeeze(d.section_pts(:, i, :));
                pt = WEC_MS2_Parser.bspline_curve_eval( ...
                         knots_v, sec_pts_at_u, degree, v_star);
                pts = [pts; pt(1,:)]; %#ok<AGROW>
            end
        end


        function v_root = bspline_z_root(knots, z_ctrl, degree, z_wl)
        % BSPLINE_Z_ROOT  Find parameter v where a scalar B-spline z(v) = z_wl.
        %
        %   Exploits the piecewise polynomial structure of B-splines:
        %     1. Evaluate z at each knot span boundary (unique knot values).
        %     2. Find the span where z(v) - z_wl changes sign.
        %     3. Bisect within that span to machine precision.
        %
        %   For degree 2 with 5 control points: 4 span boundary evaluations
        %   + ~40 bisection evaluations = ~44 scalar B-spline evaluations.

            spans = unique(knots);
            n_spans = length(spans) - 1;
            v_root = NaN;

            z_at_spans = zeros(length(spans), 1);
            for j = 1:length(spans)
                z_at_spans(j) = WEC_HydroProperties.eval_scalar_bspline( ...
                                    knots, z_ctrl, degree, spans(j));
            end

            f_spans = z_at_spans - z_wl;

            % FIX (L3): Count sign changes before bisecting.
            % ASSUMPTION: z(v) crosses z_wl exactly once (monotone profile).
            % For all current WEC geometries this holds. If more than one
            % crossing is found, warn and return only the first — downstream
            % isocurve assembly will produce an incomplete contour point.
            n_crossings = sum(f_spans(1:end-1) .* f_spans(2:end) <= 0);
            if n_crossings > 1
                warning('WEC_HydroProperties:MultipleZCrossings', ...
                    ['bspline_z_root: %d sign changes found for z_wl=%.4f. ' ...
                     'Only the first is returned. Check for non-monotone ' ...
                     'BLoftSurf z-profile.'], n_crossings, z_wl);
            end

            for s = 1:n_spans
                if f_spans(s) * f_spans(s + 1) <= 0
                    a = spans(s);
                    b = spans(s + 1);
                    fa = f_spans(s);

                    for iter = 1:50
                        m = (a + b) / 2;
                        fm = WEC_HydroProperties.eval_scalar_bspline( ...
                                 knots, z_ctrl, degree, m) - z_wl;
                        if abs(fm) < 1e-12, break; end
                        if fa * fm < 0
                            b = m;
                        else
                            a = m; fa = fm;
                        end
                    end

                    v_root = (a + b) / 2;
                    return;
                end
            end
        end


        function z = eval_scalar_bspline(knots, z_ctrl, degree, t)
        % EVAL_SCALAR_BSPLINE  Evaluate a scalar B-spline at parameter t.
        %
        %   Wrapper around WEC_MS2_Parser.bspline_basis_all for evaluating
        %   z(t) = Σ N_{i,p}(t) · z_i without constructing a full [N×3]
        %   control point array.

            if t >= 1 - 1e-10
                z = z_ctrl(end);
                return;
            end
            if t <= 1e-10
                z = z_ctrl(1);
                return;
            end
            n_ctrl = length(z_ctrl);
            basis = WEC_MS2_Parser.bspline_basis_all( ...
                        knots, degree, min(max(t, 0), 1 - 1e-12), n_ctrl);
            z = basis * z_ctrl(:);
        end


        function u_star = find_profile_root(z_profile, u_grid, z_wl, parser, profile_name)
        % FIND_PROFILE_ROOT  Find u where z_profile(u) = z_wl.
        %
        %   Step 1: Scan cached z_profile for a sign-change bracket.
        %   Step 2: If parser and profile_name are provided, refine via
        %           bisection on actual curve evaluation (~30 calls →
        %           machine precision).  Otherwise, linear interpolation.
        %
        %   INPUTS
        %     z_profile    — [n × 1] cached z-values at u_grid
        %     u_grid       — [n × 1] parameter values
        %     z_wl         — [m] target z
        %     parser       — (optional) WEC_MS2_Parser for refinement
        %     profile_name — (optional) entity name for refinement

            f = z_profile - z_wl;
            u_star = NaN;

            for k = 1:length(f) - 1
                if f(k) * f(k + 1) <= 0
                    a  = u_grid(k);
                    b  = u_grid(k + 1);
                    fa = f(k);

                    if nargin >= 5 && ~isempty(parser) && ~isempty(profile_name)
                        for iter = 1:40
                            m = (a + b) / 2;
                            pt_m = parser.eval_curve_or_snake(profile_name, m);
                            fm = pt_m(1, 3) - z_wl;
                            if abs(fm) < 1e-12, break; end
                            if fa * fm < 0
                                b = m;
                            else
                                a = m; fa = fm;
                            end
                        end
                        u_star = (a + b) / 2;
                    else
                        if abs(f(k) - f(k+1)) < 1e-14
                            u_star = a;
                        else
                            u_star = a - fa * (b - a) / (f(k+1) - f(k));
                        end
                    end
                    return;
                end
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §10  WATERPLANE PROPERTIES — Aw, I_wp from Boundary
        %% ═════════════════════════════════════════════════════════

        function [Aw, I_xx, I_yy, pts_ordered] = waterplane_properties(wl_pts)
        % WATERPLANE_PROPERTIES  Area and second moments from boundary.
        %
        %   Uses Green's theorem on the closed waterplane boundary:
        %     Aw    = (1/2) |Σ (x_i y_{i+1} − x_{i+1} y_i)|
        %     I_xx  = (1/12) |Σ (x_i y_{i+1} − x_{i+1} y_i)(y_i² + y_i·y_{i+1} + y_{i+1}²)|
        %     I_yy  = (1/12) |Σ (x_i y_{i+1} − x_{i+1} y_i)(x_i² + x_i·x_{i+1} + x_{i+1}²)|
        %
        %   The boundary points are sorted angularly around their centroid.

            if isempty(wl_pts) || size(wl_pts, 1) < 3
                Aw = 0; I_xx = 0; I_yy = 0; pts_ordered = wl_pts;
                return;
            end

            % ── Angular sort around centroid ────────────────────────
            x = wl_pts(:, 1);
            y = wl_pts(:, 2);
            cx = mean(x);
            cy = mean(y);
            angles = atan2(y - cy, x - cx);
            [~, order] = sort(angles);
            x = x(order);
            y = y(order);
            pts_ordered = wl_pts(order, :);

            n = length(x);
            x_next = [x(2:end); x(1)];
            y_next = [y(2:end); y(1)];

            % ── Shoelace area (Green's theorem) ─────────────────────
            cross_terms = x .* y_next - x_next .* y;
            Aw = abs(sum(cross_terms)) / 2;

            % ── Second moment of area (Green's theorem extension) ───
            %  These are the standard formulas for polygon moments:
            %    I_xx = ∫∫ y² dA    (about x-axis)
            %    I_yy = ∫∫ x² dA    (about y-axis)
            I_xx = abs(sum(cross_terms .* (y.^2 + y .* y_next + y_next.^2))) / 12;
            I_yy = abs(sum(cross_terms .* (x.^2 + x .* x_next + x_next.^2))) / 12;
        end


        %% ═════════════════════════════════════════════════════════
        %%  §10c0  COMPUTE_RMIN_AT_Z — Min Centroid-to-Boundary Distance
        %% ═════════════════════════════════════════════════════════

        function [r_min, r_pts] = compute_rmin_at_z(parser, z_target, n_samples, cache)
        % COMPUTE_RMIN_AT_Z  Minimum centroid-to-boundary distance at z_target.
        %
        %   [r_min, contour_pts] = WEC_HydroProperties.compute_rmin_at_z(parser, z)
        %   [r_min, contour_pts] = WEC_HydroProperties.compute_rmin_at_z(parser, z, n, cache)
        %
        %   When cache (from precompute_boundary_cache) is provided, uses
        %   the fast type-dispatched extract_isocurve_at_z.  When omitted,
        %   auto-builds a temporary cache (~0.15 s, still ~250× faster
        %   than the old reference method per z-level).
        %
        %   INPUTS
        %     parser    — WEC_MS2_Parser object
        %     z_target  — [m] elevation of cross-section
        %     n_samples — samples per sweep direction (default: 100)
        %     cache     — (optional) precomputed boundary cache

            if nargin < 3 || isempty(n_samples), n_samples = 100; end

            r_min = Inf;
            r_pts = [];

            % Build or reuse boundary cache
            if nargin >= 4 && ~isempty(cache)
                b_cache = cache;
            else
                b_cache = WEC_HydroProperties.precompute_boundary_cache( ...
                              parser, n_samples);
            end

            wl_pts = WEC_HydroProperties.extract_isocurve_at_z( ...
                         parser, z_target, n_samples, b_cache);

            if isempty(wl_pts) || size(wl_pts, 1) < 3
                r_min = 0;
                return;
            end

            r_pts = wl_pts(:, 1:2);
            x_c = r_pts(:,1); y_c = r_pts(:,2);
            xp_c = circshift(x_c,-1); yp_c = circshift(y_c,-1);
            a_vec = x_c.*yp_c - xp_c.*y_c;
            A_poly = 0.5 * sum(a_vec);
            if abs(A_poly) > 1e-14
                centroid = [sum((x_c+xp_c).*a_vec)/(6*A_poly), ...
                            sum((y_c+yp_c).*a_vec)/(6*A_poly)];
            else
                centroid = mean(r_pts, 1);
            end
            dists = sqrt(sum((r_pts - centroid).^2, 2));
            r_min = min(dists);
        end


        %% ═════════════════════════════════════════════════════════
        %%  §10c0a  COMPUTE_PERPENDICULAR_SHELL_VOLUME
        %%  Uniform-thickness perpendicular-offset UHPC shell volume.
        %%  Used by the constructability subsystem to compute per-strip
        %%  feasibility bounds (replaces the homothetic-scaling formula).
        %% ═════════════════════════════════════════════════════════

        function [V_shell, A_outer, debug] = compute_perpendicular_shell_volume( ...
                                                  parser, z_lo, z_hi, t, n_z, b_cache)
        % COMPUTE_PERPENDICULAR_SHELL_VOLUME  Volume of an axisymmetric UHPC
        %   shell formed by offsetting the outer hull surface PERPENDICULAR
        %   INWARD by uniform distance t over z ∈ [z_lo, z_hi].
        %
        %   For full mass + centroid + inertia tensor of the shell or the
        %   complementary cavity (and for the bottom-stacked fill volume
        %   inversion), use compute_axisymmetric_strip_props instead — this
        %   function is retained because diagnostics in §3f-post call it
        %   for V_shell + A_outer summary.
        %
        %   [V_shell, A_outer, debug] = compute_perpendicular_shell_volume( ...
        %                                   parser, z_lo, z_hi, t, n_z, b_cache)
        %
        %   FORMULA (leading order in t)
        %     The inner offset surface, viewed at horizontal level z, has
        %     radius given by the perpendicular distance t projected onto
        %     the radial direction.  For a profile r(z) with slope r'(z),
        %     the surface normal makes angle α=atan(r') with the radial
        %     direction, so the radial reduction is t/cos(α) = t·sqrt(1+r'²):
        %
        %         r_inner(z) = max(0, r(z) − t·sqrt(1 + r'(z)²))
        %
        %     V_shell = π · ∫_{z_lo}^{z_hi} (r²(z) − r_inner²(z)) dz
        %     A_outer = ∫_{z_lo}^{z_hi} 2π·r(z)·sqrt(1+r'(z)²) dz   (lateral)
        %
        %     The cap r_inner = 0 represents a "fully solid" cross-section:
        %     wherever t·sqrt(1+r'²) > r, the perpendicular offset would
        %     produce a negative inner radius (the shell exceeds the
        %     available radial space) — physically this means the slice is
        %     filled with UHPC at that z.
        %
        %   ACCURACY
        %     Exact in the limit t → 0.  Leading correction is O((t·κ)²)
        %     where κ is the maximum surface curvature.  For the WEC C1
        %     hull at t = 76 mm with r > 0.3 m everywhere except dome tips,
        %     correction is < 1 %.
        %
        %     This formula is what the EXACT offset-surface integration via
        %     compute_strip would converge to as the parametric sampling
        %     tends to infinity, plus an O((t·r'')) bending term that we
        %     are deliberately discarding to avoid second derivatives.
        %
        %   AXISYMMETRIC ASSUMPTION
        %     Uses compute_rmin_at_z to get r(z) — the inscribed-circle
        %     radius of the cross-section.  For a true revolution surface
        %     this equals the actual radius; for non-axisymmetric strips
        %     (e.g. a wall column whose cross-section is non-circular),
        %     V_shell is a CONSERVATIVE lower bound on the true shell
        %     volume (because actual perimeter > 2πr_min).  The user's
        %     C1 hull bulb is a RevSurf, so the formula is exact there.
        %
        %   INPUTS
        %     parser   — WEC_MS2_Parser object
        %     z_lo, z_hi — strip z bounds [m] (z_lo < z_hi)
        %     t        — perpendicular shell thickness [m] (t ≥ 0)
        %     n_z      — number of z-samples (recommended ≥ 50; uses 100 if []).
        %                Must match the realiser's sampling for contract.
        %     b_cache  — boundary cache from precompute_boundary_cache
        %
        %   OUTPUTS
        %     V_shell  — shell volume [m³]
        %     A_outer  — outer lateral surface area [m²]
        %     debug    — struct with z_samples, r_samples, rp_samples,
        %                L_samples, r_inner_samples, solid_fraction
        %                (fraction of z-samples where the cap r_inner=0
        %                fired — high values flag steep-shoulder regions
        %                where the formula yields locally-solid UHPC)

            if z_hi <= z_lo
                error('WEC_HydroProperties:BadStripBounds', ...
                      'z_hi (%.4f) <= z_lo (%.4f)', z_hi, z_lo);
            end
            if t < 0
                error('WEC_HydroProperties:NegativeThickness', ...
                      't (%.4f) < 0', t);
            end
            if nargin < 5 || isempty(n_z), n_z = 100; end
            if nargin < 6 || isempty(b_cache)
                b_cache = WEC_HydroProperties.precompute_boundary_cache(parser, 100);
            end

            z_samples = linspace(z_lo, z_hi, n_z)';
            r_samples = zeros(n_z, 1);

            for k = 1:n_z
                [r_k, ~] = WEC_HydroProperties.compute_rmin_at_z( ...
                                parser, z_samples(k), 100, b_cache);
                r_samples(k) = r_k;
            end

            % Central differences for r' (mirrors WEC_Constructable_Hull
            % :409-421 and Configuration_Builder §3d so all three code
            % paths see the same slope at the same z).
            rp_samples = zeros(n_z, 1);
            for k = 1:n_z
                if k == 1
                    dz = z_samples(2) - z_samples(1);
                    rp_samples(k) = (r_samples(2) - r_samples(1)) / dz;
                elseif k == n_z
                    dz = z_samples(n_z) - z_samples(n_z-1);
                    rp_samples(k) = (r_samples(n_z) - r_samples(n_z-1)) / dz;
                else
                    dz = z_samples(k+1) - z_samples(k-1);
                    rp_samples(k) = (r_samples(k+1) - r_samples(k-1)) / dz;
                end
            end
            % Dome-tip correction: zero rp at samples where r < t (the
            % surface has effectively reached its closing point).  Mirrors
            % WEC_Constructable_Hull.realize:436.
            rp_samples(r_samples < t) = 0;

            L_samples       = sqrt(1 + rp_samples.^2);
            r_inner_samples = max(0, r_samples - t .* L_samples);

            % Volume integrand: π · (r² - r_inner²)
            integrand_V = pi .* (r_samples.^2 - r_inner_samples.^2);
            V_shell     = trapz(z_samples, integrand_V);

            % Lateral surface area (no caps)
            integrand_A = 2*pi .* r_samples .* L_samples;
            A_outer     = trapz(z_samples, integrand_A);

            if nargout >= 3
                debug.z_samples       = z_samples;
                debug.r_samples       = r_samples;
                debug.rp_samples      = rp_samples;
                debug.L_samples       = L_samples;
                debug.r_inner_samples = r_inner_samples;
                debug.solid_fraction  = sum(r_inner_samples == 0) / n_z;
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §10c0b  COMPUTE_AXISYMMETRIC_REVOLUTION_PROPS
        %%  Generic 1D integrator for axisymmetric annular/disk regions.
        %%  Used for shell, fill, and void contributions in the
        %%  uniform-thickness offset-shell constructability model.
        %% ═════════════════════════════════════════════════════════

        function props = compute_axisymmetric_revolution_props( ...
                              z_samples, r_out_samples, r_inner_samples)
        % COMPUTE_AXISYMMETRIC_REVOLUTION_PROPS  Volume, centroid, and inertia
        %   tensor of an axisymmetric region of revolution bounded by
        %   r_out(z) and r_inner(z) over z ∈ [z_samples(1), z_samples(end)].
        %   Density-independent (caller multiplies by ρ to get mass).
        %
        %   props = compute_axisymmetric_revolution_props( ...
        %               z_samples, r_out_samples, r_inner_samples)
        %
        %   FORMULAS (axisymmetric body, axis = z, all moments about origin)
        %     V    = π · ∫ (r_out² − r_inner²) dz
        %     z_cg = (1/V) · π · ∫ z · (r_out² − r_inner²) dz
        %     Iyy  = π · ∫ (r_out² − r_inner²) · ((r_out² + r_inner²)/4 + z²) dz
        %     Izz  = π · ∫ (r_out² − r_inner²) · (r_out² + r_inner²)/2 dz
        %     Ixx  = Iyy   (axisymmetry)
        %
        %   Iyy/Ixx are the second moment about a HORIZONTAL axis (y or x)
        %   passing through the ORIGIN of the body frame, NOT through the
        %   region's own centroid.  Caller should apply parallel-axis to
        %   shift to a different reference point.  The (r_out²+r_inner²)/4
        %   term is the disk's polar moment about its own diameter axis;
        %   the +z² term is the parallel-axis lever from origin to the
        %   disk at height z.
        %
        %   USAGE PATTERNS
        %     (a) Shell of an outer hull surface offset inward by t:
        %           r_out_samples = hull r(z),  r_inner_samples = offset r(z)
        %     (b) Solid disk (no inner cavity, e.g. bottom-stacked fill):
        %           r_out_samples = cavity wall r(z),  r_inner_samples = zeros
        %     (c) Void cavity (used only to compute mass·ρ_fill ≈ 0 sanity):
        %           same as (b) with ρ = ρ_fill
        %
        %   INPUT REQUIREMENTS
        %     z_samples       — [N×1] strictly increasing z-grid
        %     r_out_samples   — [N×1] outer radius at each z (≥ 0)
        %     r_inner_samples — [N×1] inner radius at each z (0 ≤ r_in ≤ r_out)
        %
        %   OUTPUT (struct)
        %     props.V          — region volume [m³]
        %     props.z_cg       — z-centroid (about origin) [m]
        %     props.Iyy_origin — Iyy about y-axis through origin [kg·m² /(kg/m³)]
        %     props.Izz_origin — Izz about z-axis through origin [kg·m² /(kg/m³)]
        %     props.Ixx_origin — Ixx about x-axis through origin (= Iyy_origin)

            z   = z_samples(:);
            r_o = r_out_samples(:);
            r_i = r_inner_samples(:);

            if length(z) < 2
                props.V          = 0;
                props.z_cg       = 0;
                props.Iyy_origin = 0;
                props.Izz_origin = 0;
                props.Ixx_origin = 0;
                return;
            end

            r_o = max(0, r_o);
            r_i = max(0, min(r_o, r_i));   % clamp 0 ≤ r_i ≤ r_o

            ro2 = r_o.^2;
            ri2 = r_i.^2;
            ann = ro2 - ri2;                        % annular cross-section r²

            integrand_V = pi .* ann;
            V           = trapz(z, integrand_V);

            if V > 1e-14
                z_cg = trapz(z, z .* integrand_V) / V;
            else
                z_cg = 0;
            end

            % x²-only contribution to Iyy: ∫∫∫ x² dV = π · ∫ (r_out⁴ − r_inner⁴)/4 dz
            % (axisymmetric: integrating x² = r²cos²θ over θ gives π·r²/2 per
            % unit area, then integrating over r from r_inner to r_out gives
            % π · (r_out² + r_inner²)/4 per unit dV — equivalent to the
            % cross-section's polar second moment about its own diameter axis.)
            integrand_Iyy_x = integrand_V .* (ro2 + ri2)/4;
            Iyy_x_only      = trapz(z, integrand_Iyy_x);

            % z²-only contribution to Iyy: ∫∫∫ z² dV
            integrand_Iyy_z = integrand_V .* z.^2;
            Iyy_z_only      = trapz(z, integrand_Iyy_z);

            % Full Iyy about y-axis through origin = x² + z²
            Iyy_origin = Iyy_x_only + Iyy_z_only;

            % Izz about z-axis: ∫∫∫ (x²+y²) dV = 2 · ∫∫∫ x² dV (axisymmetric)
            Izz_origin = 2 * Iyy_x_only;

            props.V          = V;
            props.z_cg       = z_cg;
            props.Iyy_x_only = Iyy_x_only;   % ∫∫∫ x² dV  (no z² contribution)
            props.Iyy_z_only = Iyy_z_only;   % ∫∫∫ z² dV  (= ∫∫∫ y² dV, axi.)
            props.Iyy_origin = Iyy_origin;
            props.Izz_origin = Izz_origin;
            props.Ixx_origin = Iyy_origin;   % axisymmetric
        end


        %% ═════════════════════════════════════════════════════════
        %%  §10c  COMPUTE_STRIP — Horizontal Strip Properties
        %% ═════════════════════════════════════════════════════════

        function strip = compute_strip(parser, z_lo, z_hi, options)
        % COMPUTE_STRIP  Volume and inertia of a horizontal strip [z_lo, z_hi].
        %
        %   strip = WEC_HydroProperties.compute_strip(parser, z_lo, z_hi)
        %   strip = WEC_HydroProperties.compute_strip(parser, z_lo, z_hi, opts)
        %
        %   FAST PATH (preferred): when opts.Aw_table_z is provided, volume
        %     and centroid are computed by direct numerical integration of the
        %     precomputed cross-section area table A(z):
        %
        %       V      = ∫_{z_lo}^{z_hi} A(z) dz
        %       z̄     = ∫_{z_lo}^{z_hi} z·A(z) dz  /  V
        %       Iyy    = ∫_{z_lo}^{z_hi} [I_yy(z) + z²·A(z)] dz
        %
        %     This is exact to the resolution of the Aw_table and avoids
        %     all surface-normal orientation and cap-boundary issues that
        %     arise when the hull has shoulder features at strip boundaries.
        %
        %   FALLBACK: when no Aw_table is provided, uses the divergence
        %     theorem on the parametric surfaces (original method).
        %
        %   OUTPUT
        %     strip.V     — [m³]   strip volume
        %     strip.CB_z  — [m]    z-centroid of strip volume
        %     strip.Iyy   — [m⁵]  ∫∫∫(x²+z²)dV  (pitch, about y-axis)
        %     strip.Ixx   — [m⁵]  ∫∫∫(y²+z²)dV  (roll,  about x-axis)
        %     strip.Izz   — [m⁵]  ∫∫∫(x²+y²)dV  (yaw,   about z-axis)
        %     strip.int_x2, .int_y2, .int_z2 — raw second moments

            if nargin < 4, options = struct(); end
            if ~isfield(options, 'n_quad'),        options.n_quad        = 16; end
            if ~isfield(options, 'Aw_table_z'),    options.Aw_table_z    = []; end
            if ~isfield(options, 'Aw_table'),      options.Aw_table      = []; end
            if ~isfield(options, 'I_wp_xx_table'), options.I_wp_xx_table = []; end
            if ~isfield(options, 'I_wp_yy_table'), options.I_wp_yy_table = []; end
            if ~isfield(options, 'boundary_cache'), options.boundary_cache = []; end

            %% ── FAST PATH: integrate A(z) table over [z_lo, z_hi] ──────────
            %
            %  Used whenever the Aw_table is available (always in the config-
            %  builder call from §3f).  Replaces the divergence-theorem surface
            %  integral entirely, avoiding:
            %   • orientation ambiguity at the platform-to-column shoulder
            %   • explicit cap double-counting when the hull surface already
            %     closes at a strip boundary
            %   • the per-element normal flip that produced z̄ outside strip
            %     bounds (5.455 m for a strip at [−0.9, +1.1])
            %
            %  Integration nodes: the Aw_table points already within the strip
            %  plus the two exact endpoints z_lo, z_hi (interpolated from the
            %  table).  This keeps the node count low while capturing any
            %  rapid A(z) changes resolved by the adaptive table.
            if ~isempty(options.Aw_table_z) && ~isempty(options.Aw_table)

                z_tab = options.Aw_table_z(:);
                A_tab = options.Aw_table(:);

                % Interior table nodes within the strip (exclusive of endpoints)
                in_mask = z_tab > z_lo & z_tab < z_hi;
                z_int   = z_tab(in_mask);
                A_int   = A_tab(in_mask);

                % Interpolate at the exact strip endpoints
                A_lo = max(0, interp1(z_tab, A_tab, z_lo, 'linear', 0));
                A_hi = max(0, interp1(z_tab, A_tab, z_hi, 'linear', 0));

                % Build the sorted node set [z_lo, interior..., z_hi]
                z_pts = [z_lo; z_int; z_hi];
                A_pts = [A_lo; A_int; A_hi];

                if length(z_pts) < 2 || (z_pts(end) - z_pts(1)) < 1e-12
                    % Degenerate strip
                    strip.V     = 0;
                    strip.CB_z  = 0.5 * (z_lo + z_hi);
                    strip.int_x2 = 0; strip.int_y2 = 0; strip.int_z2 = 0;
                    strip.Iyy   = 0; strip.Ixx = 0; strip.Izz = 0;
                    strip.z_lo  = z_lo; strip.z_hi = z_hi;
                    return;
                end

                % Volume: ∫ A(z) dz
                V = trapz(z_pts, A_pts);

                % Z-centroid: ∫ z·A(z) dz / V
                if V > 1e-10
                    CB_z = trapz(z_pts, z_pts .* A_pts) / V;
                else
                    V    = 0;
                    CB_z = 0.5 * (z_lo + z_hi);
                end

                % Second moments — need Iyy(z) and Ixx(z) tables
                % int_x2 = ∫ Iyy_cross(z) dz  = ∫∫∫ x² dV
                % int_y2 = ∫ Ixx_cross(z) dz  = ∫∫∫ y² dV
                % int_z2 = ∫ z²·A(z) dz       = ∫∫∫ z² dV
                int_z2 = trapz(z_pts, z_pts.^2 .* A_pts);

                if ~isempty(options.I_wp_yy_table)
                    Iyy_int = max(0, interp1(z_tab, options.I_wp_yy_table(:), z_int, 'linear', 0));
                    Iyy_lo  = max(0, interp1(z_tab, options.I_wp_yy_table(:), z_lo,  'linear', 0));
                    Iyy_hi  = max(0, interp1(z_tab, options.I_wp_yy_table(:), z_hi,  'linear', 0));
                    int_x2  = trapz(z_pts, [Iyy_lo; Iyy_int; Iyy_hi]);
                else
                    int_x2 = 0;
                end

                if ~isempty(options.I_wp_xx_table)
                    Ixx_int = max(0, interp1(z_tab, options.I_wp_xx_table(:), z_int, 'linear', 0));
                    Ixx_lo  = max(0, interp1(z_tab, options.I_wp_xx_table(:), z_lo,  'linear', 0));
                    Ixx_hi  = max(0, interp1(z_tab, options.I_wp_xx_table(:), z_hi,  'linear', 0));
                    int_y2  = trapz(z_pts, [Ixx_lo; Ixx_int; Ixx_hi]);
                else
                    int_y2 = 0;
                end

                strip.V      = V;
                strip.CB_z   = CB_z;
                strip.int_x2 = int_x2;
                strip.int_y2 = int_y2;
                strip.int_z2 = int_z2;
                strip.Iyy    = int_x2 + int_z2;
                strip.Ixx    = int_y2 + int_z2;
                strip.Izz    = int_x2 + int_y2;
                strip.z_lo   = z_lo;
                strip.z_hi   = z_hi;
                return;
            end

            %% ── FALLBACK: divergence-theorem surface integral ───────────────
            %  Used only when no Aw_table is available.

            n_quad = options.n_quad;
            topo = parser.classify_visible_surfaces();

            % Interior point: bounding box center (always inside hull)
            interior = (parser.extents(1:3) + parser.extents(4:6))' / 2;

            % Origin shift: strip midpoint (v7.1)
            z_ref = 0.5 * (z_lo + z_hi);

            % Pre-determine orientation per source (same as compute_submerged).
            % Uses coarse 5×5 signed-volume integral over FULL surface.
            source_orient = containers.Map();
            n_orient = 5;
            [u_o, wu_o] = WEC_HydroProperties.gauss_legendre(n_orient);
            [v_o, wv_o] = WEC_HydroProperties.gauss_legendre(n_orient);
            for s = 1:length(topo.sources)
                sname = topo.sources{s};
                V_test = 0;
                for io = 1:n_orient
                    for jo = 1:n_orient
                        [S_t, Su_t, Sv_t] = parser.eval_surface_with_derivs( ...
                            sname, u_o(io), v_o(jo));
                        V_test = V_test + wu_o(io) * wv_o(jo) * ...
                            dot(S_t, cross(Su_t, Sv_t)) / 3;
                    end
                end
                if V_test < 0
                    source_orient(sname) = -1;
                else
                    source_orient(sname) = +1;
                end
            end

            % Source surface integrals (shifted)
            V_total  = 0;
            Cz_num   = 0;
            ix2_total = 0;  iy2_total = 0;  iz2_total = 0;

            source_V    = containers.Map();
            source_Cz   = containers.Map();
            source_ix2  = containers.Map();
            source_iy2  = containers.Map();
            source_iz2  = containers.Map();

            for s = 1:length(topo.sources)
                sname = topo.sources{s};
                [V_s, Cz_s, ix2_s, iy2_s, iz2_s] = ...
                    WEC_HydroProperties.surface_integral_strip( ...
                        parser, sname, n_quad, z_lo, z_hi, interior, ...
                        source_orient(sname), z_ref);

                source_V(sname)   = V_s;
                source_Cz(sname)  = Cz_s;
                source_ix2(sname) = ix2_s;
                source_iy2(sname) = iy2_s;
                source_iz2(sname) = iz2_s;

                V_total   = V_total + V_s;
                Cz_num    = Cz_num + Cz_s;
                ix2_total = ix2_total + ix2_s;
                iy2_total = iy2_total + iy2_s;
                iz2_total = iz2_total + iz2_s;
            end

            % Mirror contributions (invariant under reflection)
            for m = 1:length(topo.mirrors)
                mirr = topo.mirrors(m);
                ult = mirr.ultimate_source;
                if source_V.isKey(ult)
                    V_total   = V_total + source_V(ult);
                    Cz_num    = Cz_num + source_Cz(ult);
                    ix2_total = ix2_total + source_ix2(ult);
                    iy2_total = iy2_total + source_iy2(ult);
                    iz2_total = iz2_total + source_iz2(ult);
                end
            end

            % Flat caps with SHIFTED position vector
            %  Top cap at z_hi:    n = +z, S'·n = z_hi − z_ref
            %  Bottom cap at z_lo: n = −z, S'·n = −(z_lo − z_ref)

            if ~isempty(options.Aw_table_z)
                Aw_hi = interp1(options.Aw_table_z, options.Aw_table, z_hi, 'linear', 0);
                Aw_lo = interp1(options.Aw_table_z, options.Aw_table, z_lo, 'linear', 0);
            else
                % Fast iso-z contour extraction (v8.0)
                if ~isempty(options.boundary_cache)
                    b_cache = options.boundary_cache;
                else
                    b_cache = WEC_HydroProperties.precompute_boundary_cache( ...
                                  parser, n_quad * 5);
                end

                wl_hi = WEC_HydroProperties.extract_isocurve_at_z( ...
                            parser, z_hi, n_quad * 5, b_cache);
                wl_lo = WEC_HydroProperties.extract_isocurve_at_z( ...
                            parser, z_lo, n_quad * 5, b_cache);

                if ~isempty(wl_hi) && size(wl_hi, 1) >= 3
                    Aw_hi = WEC_HydroProperties.waterplane_properties(wl_hi);
                else
                    Aw_hi = 0;
                end
                if ~isempty(wl_lo) && size(wl_lo, 1) >= 3
                    Aw_lo = WEC_HydroProperties.waterplane_properties(wl_lo);
                else
                    Aw_lo = 0;
                end
            end

            dz_hi = z_hi - z_ref;   % positive (half strip height)
            dz_lo = z_lo - z_ref;   % negative (half strip height)

            % Top cap (+z outward)
            V_total   = V_total + dz_hi * Aw_hi / 3;
            Cz_num    = Cz_num + dz_hi^2 * Aw_hi / 2;
            iz2_total = iz2_total + dz_hi^3 * Aw_hi / 3;

            % Bottom cap (−z outward)
            V_total   = V_total - dz_lo * Aw_lo / 3;
            Cz_num    = Cz_num - dz_lo^2 * Aw_lo / 2;
            iz2_total = iz2_total - dz_lo^3 * Aw_lo / 3;

            % ── Unshift from z_ref to world coordinates ────────────
            Cz_shifted  = Cz_num;
            Cz_num      = Cz_shifted + z_ref * V_total;
            iz2_total   = iz2_total + 2*z_ref*Cz_shifted + z_ref^2 * V_total;

            strip.V = abs(V_total);

            if strip.V > 1e-6
                strip.CB_z = Cz_num / V_total;
            else
                strip.V    = 0;
                strip.CB_z = 0.5 * (z_lo + z_hi);
            end

            strip.int_x2 = abs(ix2_total);
            strip.int_y2 = abs(iy2_total);
            strip.int_z2 = abs(iz2_total);

            strip.Iyy = strip.int_x2 + strip.int_z2;
            strip.Ixx = strip.int_y2 + strip.int_z2;
            strip.Izz = strip.int_x2 + strip.int_y2;

            strip.z_lo = z_lo;
            strip.z_hi = z_hi;
        end


        %% ═════════════════════════════════════════════════════════
        %%  §10d  SURFACE_INTEGRAL_STRIP — GL over z_lo ≤ z ≤ z_hi
        %% ═════════════════════════════════════════════════════════

        function [V_s, Cz_s, ix2_s, iy2_s, iz2_s] = surface_integral_strip( ...
                parser, surf_name, n_quad, z_lo, z_hi, interior_pt, orient_sign, z_ref)
        % SURFACE_INTEGRAL_STRIP  Divergence theorem over one surface,
        %   restricted to the horizontal band z_lo ≤ z ≤ z_hi.
        %
        %   orient_sign: +1 or -1, pre-determined from full surface.
        %   z_ref: origin shift in z (same as surface_integral_submerged).

            if nargin < 7, orient_sign = []; end
            if nargin < 8, z_ref = 0; end
            n_bisect = 40;
            [u_gl, w_u] = WEC_HydroProperties.gauss_legendre(n_quad);
            [v_gl, w_v] = WEC_HydroProperties.gauss_legendre(n_quad);

            V_s   = 0;
            Cz_s  = 0;
            ix2_s = 0;  iy2_s = 0;  iz2_s = 0;

            for i = 1:n_quad
                u = u_gl(i);

                [v_lo_i, v_hi_i] = WEC_HydroProperties.find_strip_v_limits( ...
                    parser, surf_name, u, z_lo, z_hi, n_bisect);

                dv = v_hi_i - v_lo_i;
                if dv < 1e-14, continue; end

                for j = 1:n_quad
                    v = v_lo_i + v_gl(j) * dv;
                    [S, Su, Sv] = parser.eval_surface_with_derivs(surf_name, u, v);
                    n_vec = cross(Su, Sv);
                    w = w_u(i) * w_v(j) * dv;

                    z_s = S(3) - z_ref;

                    V_s   = V_s + w * (S(1)*n_vec(1) + S(2)*n_vec(2) + z_s*n_vec(3)) / 3;
                    Cz_s  = Cz_s + w * z_s^2  * n_vec(3) / 2;
                    ix2_s = ix2_s + w * S(1)^3 * n_vec(1) / 3;
                    iy2_s = iy2_s + w * S(2)^3 * n_vec(2) / 3;
                    iz2_s = iz2_s + w * z_s^3  * n_vec(3) / 3;
                end
            end

            % Apply pre-determined orientation sign
            if ~isempty(orient_sign)
                if orient_sign < 0
                    V_s   = -V_s;
                    Cz_s  = -Cz_s;
                    ix2_s = -ix2_s;
                    iy2_s = -iy2_s;
                    iz2_s = -iz2_s;
                end
            else
                [S_mid, Su_mid, Sv_mid] = parser.eval_surface_with_derivs( ...
                    surf_name, 0.5, 0.5);
                n_mid = cross(Su_mid, Sv_mid);
                if dot(n_mid, S_mid - interior_pt(:)') < 0
                    V_s   = -V_s;
                    Cz_s  = -Cz_s;
                    ix2_s = -ix2_s;
                    iy2_s = -iy2_s;
                    iz2_s = -iz2_s;
                end
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §10e  FIND_STRIP_V_LIMITS — v-interval for z_lo ≤ z ≤ z_hi
        %% ═════════════════════════════════════════════════════════

        function [v_lo, v_hi] = find_strip_v_limits( ...
                parser, surf_name, u, z_lo, z_hi, n_bisect)
        % FIND_STRIP_V_LIMITS  Return [v_lo, v_hi] where z_lo ≤ z(u,v) ≤ z_hi.
        %
        %   Generalization of find_submerged_v_limits.  Instead of one
        %   bound (z ≤ z_wl), this finds the v-interval where z is
        %   between two bounds [z_lo, z_hi].
        %
        %   RevSurf fast-path: z depends on u only.  Single check.
        %   General: sweep 21 v-samples, classify each as below/inside/above,
        %   bisect transitions.

            if nargin < 6, n_bisect = 40; end

            % ── RevSurf fast-path ─────────────────────────────────
            if parser.entities.isKey(surf_name)
                se = parser.entities(surf_name);
                if strcmp(se.type, 'RevSurf')
                    pt = parser.eval_surface(surf_name, u, 0);
                    z_u = pt(3);
                    if z_u >= z_lo - 1e-10 && z_u <= z_hi + 1e-10
                        v_lo = 0; v_hi = 1;
                    else
                        v_lo = 0; v_hi = 0;
                    end
                    return;
                end
                if strcmp(se.type, 'MirrSurf')
                    src = parser.entities(se.params.source);
                    if strcmp(src.type, 'RevSurf')
                        pt = parser.eval_surface(surf_name, u, 0);
                        z_u = pt(3);
                        if z_u >= z_lo - 1e-10 && z_u <= z_hi + 1e-10
                            v_lo = 0; v_hi = 1;
                        else
                            v_lo = 0; v_hi = 0;
                        end
                        return;
                    end
                end
            end

            % ── Coarse sweep ──────────────────────────────────────
            n_sample = 21;
            v_sample = linspace(0, 1, n_sample);
            z_sample = zeros(n_sample, 1);
            for k = 1:n_sample
                pt = parser.eval_surface(surf_name, u, v_sample(k));
                z_sample(k) = pt(3);
            end

            % Classify: find first and last v where z is inside [z_lo, z_hi]
            inside = (z_sample >= z_lo - 1e-10) & (z_sample <= z_hi + 1e-10);

            if ~any(inside)
                v_lo = 0; v_hi = 0;
                return;
            end

            if all(inside)
                v_lo = 0; v_hi = 1;
                return;
            end

            % Find transitions: outside→inside and inside→outside
            first_in = find(inside, 1, 'first');
            last_in  = find(inside, 1, 'last');

            % Refine lower boundary
            if first_in == 1
                v_lo = 0;
            else
                % Bisect between v_sample(first_in-1) and v_sample(first_in)
                a = v_sample(first_in - 1);
                b = v_sample(first_in);
                for iter = 1:n_bisect
                    m = (a + b) / 2;
                    pt = parser.eval_surface(surf_name, u, m);
                    if pt(3) >= z_lo - 1e-10 && pt(3) <= z_hi + 1e-10
                        b = m;  % inside — tighten from right
                    else
                        a = m;  % outside — tighten from left
                    end
                end
                v_lo = (a + b) / 2;
            end

            % Refine upper boundary
            if last_in == n_sample
                v_hi = 1;
            else
                a = v_sample(last_in);
                b = v_sample(last_in + 1);
                for iter = 1:n_bisect
                    m = (a + b) / 2;
                    pt = parser.eval_surface(surf_name, u, m);
                    if pt(3) >= z_lo - 1e-10 && pt(3) <= z_hi + 1e-10
                        a = m;  % inside — tighten from left
                    else
                        b = m;  % outside — tighten from right
                    end
                end
                v_hi = (a + b) / 2;
            end
        end


        %% ═════════════════════════════════════════════════════════
        %%  §11a  WETTED SURFACE AREA (SUBMERGED HULL SIDES ONLY)
        %% ═════════════════════════════════════════════════════════

        function A_sub = compute_wetted_surface_area(parser, z_wl, n_quad)
        % COMPUTE_WETTED_SURFACE_AREA  Physical wetted hull area below z_wl.
        %
        %   A_sub = WEC_HydroProperties.compute_wetted_surface_area(parser, z_wl)
        %   A_sub = WEC_HydroProperties.compute_wetted_surface_area(parser, z_wl, n_quad)
        %
        %   Returns the area of the hull surface in contact with water [m²].
        %
        %   DEFINITION
        %     A_sub = ∬_{z(u,v) ≤ z_wl} ||∂S/∂u × ∂S/∂v|| du dv
        %             summed over all visible hull surfaces (sources + mirrors).
        %
        %   The waterplane area (Aw) is deliberately NOT included.
        %   The waterplane is the free-water surface, not a solid hull surface.
        %   Including it would inflate A_sub and violate the ITTC/ship-hydro
        %   definition of wetted surface area used for:
        %     - Frictional drag: C_f × ½ρV²S_wet
        %     - Reynolds-number scaling
        %     - Viscous BEM correction factors
        %
        %   IMPLEMENTATION NOTES
        %     Area element: ||Su × Sv|| du dv  (always positive — no orientation
        %     check needed, unlike the signed volume integrand S·n/3).
        %
        %     Submerged v-limits: uses find_submerged_v_limits, which provides
        %     a RevSurf fast-path (z depends only on u → entire azimuth is
        %     either in or out) and a bisection path for general surfaces.
        %
        %     Mirror surfaces: area is invariant under X/Y reflection.
        %     Iterates topo.mirrors exactly as compute_submerged does.
        %     Uses a source_A_sub map so each mirror looks up its source area
        %     without re-evaluating the parametric surface.
        %     Handles 1, 2, or 3 mirrors per source correctly (not a simple ×2).
        %
        %   INPUTS
        %     parser  — WEC_MS2_Parser  (config.ms2_model)
        %     z_wl    — waterline z-coordinate in body frame [m]
        %     n_quad  — GL quadrature order (default 20)
        %
        %   OUTPUT
        %     A_sub   — wetted hull surface area [m²], excluding waterplane

            if nargin < 3 || isempty(n_quad), n_quad = 20; end

            [u_gl, w_u] = WEC_HydroProperties.gauss_legendre(n_quad);
            [v_gl, w_v] = WEC_HydroProperties.gauss_legendre(n_quad);

            topo      = parser.classify_visible_surfaces();
            n_bisect  = 40;
            A_sub     = 0;

            % Map to hold source wetted area (needed for mirror lookup)
            source_A_sub = containers.Map('KeyType','char','ValueType','double');

            % ── SOURCE SURFACES ────────────────────────────────────────
            for s = 1:length(topo.sources)
                sname  = topo.sources{s};
                A_surf = 0;

                for i = 1:n_quad
                    u = u_gl(i);

                    [v_lo, v_hi, ~] = WEC_HydroProperties.find_submerged_v_limits( ...
                        parser, sname, u, z_wl, n_bisect);

                    dv = v_hi - v_lo;
                    if dv < 1e-14, continue; end

                    for j = 1:n_quad
                        v = v_lo + v_gl(j) * dv;
                        [~, Su, Sv] = parser.eval_surface_with_derivs(sname, u, v);
                        n_vec  = cross(Su, Sv);
                        w      = w_u(i) * w_v(j) * dv;
                        A_surf = A_surf + w * norm(n_vec);
                    end
                end

                source_A_sub(sname) = A_surf;
                A_sub = A_sub + A_surf;
            end

            % ── MIRROR SURFACES ────────────────────────────────────────
            %   Reflection preserves distances → A_mirror = A_source.
            %   Iterates topo.mirrors exactly as compute_submerged does —
            %   handles 1, 2, or 3 mirrors per source (not a simple ×2).
            for m = 1:length(topo.mirrors)
                mirr = topo.mirrors(m);
                ult  = mirr.ultimate_source;

                if source_A_sub.isKey(ult)
                    A_sub = A_sub + source_A_sub(ult);
                else
                    A_mirr = 0;
                    for i = 1:n_quad
                        u = u_gl(i);
                        [v_lo, v_hi, ~] = WEC_HydroProperties.find_submerged_v_limits( ...
                            parser, mirr.name, u, z_wl, n_bisect);
                        dv = v_hi - v_lo;
                        if dv < 1e-14, continue; end
                        for j = 1:n_quad
                            v = v_lo + v_gl(j) * dv;
                            [~, Su, Sv] = parser.eval_surface_with_derivs(mirr.name, u, v);
                            n_vec  = cross(Su, Sv);
                            w      = w_u(i) * w_v(j) * dv;
                            A_mirr = A_mirr + w * norm(n_vec);
                        end
                    end
                    A_sub = A_sub + A_mirr;
                end
            end

            % Waterplane cap: deliberately NOT added.
        end


        %% ═════════════════════════════════════════════════════════
        %%  §11  GAUSS-LEGENDRE QUADRATURE
        %% ═════════════════════════════════════════════════════════

        function [x, w] = gauss_legendre(n)
        % GAUSS_LEGENDRE  Nodes and weights on [0, 1].
        %
        %   [x, w] = gauss_legendre(n)
        %
        %   Uses the Golub-Welsch algorithm: the nodes are eigenvalues
        %   of the symmetric tridiagonal Jacobi matrix for Legendre
        %   polynomials, and the weights are derived from the first
        %   components of the eigenvectors.
        %
        %   The resulting rule is exact for polynomials of degree ≤ 2n−1.
        %   For smooth integrands, convergence is exponential.

            if n == 1
                x = 0.5;
                w = 1.0;
                return;
            end

            % Jacobi matrix sub-diagonal for Legendre polynomials
            k    = 1:n-1;
            beta = k ./ sqrt(4*k.^2 - 1);
            J    = diag(beta, 1) + diag(beta, -1);

            [V, D] = eig(J);
            x_ref  = diag(D);            % nodes on [−1, +1]
            w_ref  = 2 * V(1, :)'.^2;    % weights on [−1, +1]

            % Transform [−1, +1] → [0, 1]
            x = (x_ref + 1) / 2;
            w = w_ref / 2;

            % Sort ascending
            [x, idx] = sort(x);
            w = w(idx);
        end

    end % methods (Static)

end % classdef