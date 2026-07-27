classdef WEC_Core_Functions
    % WEC_CORE_FUNCTIONS  Static utility library for WEC geometry and hydrodynamics.
    %
    %   This class collects all pure-geometry and pure-hydrodynamic helper
    %   functions used by the optimisation pipeline.  Every method is static
    %   (no instance state).
    %
    %   REFACTORED (v4.0): All mesh-based methods replaced by parametric
    %   B-spline evaluation via the WEC_MS2_Parser.  Property calculations
    %   (volume, CG, inertia) now use strip integration over cross-section
    %   areas derived directly from the parent curves of each parametric
    %   surface — no triangulated mesh in the loop.
    %
    %   METHOD INVENTORY
    %   ────────────────────────────────────────────────────────────────────
    %   SPLINE CROSS-SECTION (NEW)
    %     evaluateCrossSectionMS2  — parser + z_target → ordered polygon(s)
    %     extractProfileMS2        — parser → midplane [x,z] profile
    %
    %   2D POLYGON OPERATIONS (UNCHANGED)
    %     clipPolygon              — Sutherland–Hodgman clip against z-plane
    %     polygeom                 — Area, centroid, second moments (shoelace)
    %     find_waterline_intersections — Profile × waterline crossing points
    %     calculatePolygonProperties   — Aggregate area + Ixx + Iyy over cells
    %
    %   STRIP INTEGRATION (ADAPTED)
    %     compute_strip_bspline    — Volume, CB, Iyy of one z-strip via
    %                                cubic spline integration.  Now accepts
    %                                parser (not mesh) for cross-sections.
    %
    %   HYDRODYNAMIC INTERPOLATION (UNCHANGED)
    %     interpolate_wamit_added_mass — Draft-interpolated A, B (diag + 3×3)
    %     interpolate_3x3_matrix       — Element-wise interp of 3×3 cell array
    %
    %   RIGID-BODY DYNAMICS (UNCHANGED)
    %     calculate6x6MassMatrix   — 6×6 mass/inertia matrix from CG + I_tensor
    %
    %   PIECEWISE-POLYNOMIAL INTEGRATION (UNCHANGED)
    %     integrate_pp             — ∫ f(z) dz
    %     integrate_pp_times_zp    — ∫ z^p f(z) dz
    %   ────────────────────────────────────────────────────────────────────
    %
    %   REMOVED (v4.0)
    %     loadAndCacheMesh             — STL loading (replaced by MS2 parser)
    %     clipMeshByPlane              — mesh half-space clip (not needed)
    %     calculateMeshProperties      — divergence theorem on triangles
    %     calculateWaterplanePolygon   — mesh slicing (replaced by spline eval)
    %     chainSegments                — mesh segment chaining
    %
    %   DEPENDENCIES
    %     WEC_MS2_Parser — required by evaluateCrossSectionMS2
    %
    %   See also: WEC_MS2_Parser, WEC_Configuration_Builder,
    %             calculate_2d_properties, calculate_3d_properties
    %
    %   Author:  WEC Optimisation Team
    %   Version: 4.0 — Spline-based refactor (all mesh operations removed)

    methods (Static)

        %% ═════════════════════════════════════════════════════════════
        %%  HULL SURFACE EXTENTS
        %% ═════════════════════════════════════════════════════════════

        function [z_min, z_max] = compute_surface_z_range(parser, n_sample)
        % COMPUTE_SURFACE_Z_RANGE  True z-extents from surface evaluation.
        %
        %   [z_min, z_max] = WEC_Core_Functions.compute_surface_z_range(parser)
        %
        %   The .ms2 file header extents include B-spline CONTROL POINTS,
        %   which lie outside the actual surface (B-splines approximate,
        %   not interpolate interior control points).  This method samples
        %   all visible surfaces on a coarse grid to find the true z-range.
        %
        %   WHY this matters:
        %     For C0, a control point at z = -3.0 makes hull_z_min = -3.0,
        %     but the actual keel surface only reaches z ≈ -2.5.  Using
        %     the file-header extent wastes a density strip on empty space.

            if nargin < 2, n_sample = 50; end

            z_all = zeros(n_sample * n_sample * length(parser.visible_surfs), 1);
            idx = 0;
            t_s = linspace(0, 1, n_sample);

            for s = 1:length(parser.visible_surfs)
                sname = parser.visible_surfs{s};
                for ui = 1:n_sample
                    for vi = 1:n_sample
                        pt = parser.eval_surface(sname, t_s(ui), t_s(vi));
                        idx = idx + 1;
                        z_all(idx) = pt(3);
                    end
                end
            end

            z_all = z_all(1:idx);
            z_min = min(z_all);
            z_max = max(z_all);
        end


        %% ═════════════════════════════════════════════════════════════
        %%  SPLINE-BASED CROSS-SECTION EVALUATION
        %% ═════════════════════════════════════════════════════════════

        function polygons = evaluateCrossSectionMS2(parser, z_target, n_pts)
        % EVALUATECROSSSECTIONMS2  Extract cross-section polygon(s) at z_target
        %   by evaluating the parent B-spline curves of each visible surface.
        %
        %   PERFORMANCE OPTIMIZATIONS (v4.1):
        %     1. REVOLUTION PAIRS: When a RevSurf + its MirrSurf(Y=0) together
        %        cover 360°, and no other surface contributes at this z, the
        %        cross-section is a perfect circle.  Use A=πr², Iyy=πr⁴/4.
        %        No polygon generated — the return value is a circle polygon
        %        for downstream consumers that need vertices.
        %     2. SOURCE CACHING: Each unique source surface is evaluated once
        %        per z-level.  MirrSurf flips the cached result (zero cost).
        %     3. PROFILE CACHING: RevSurf profile curves are sampled once
        %        (persistent cache) and reused for all z-levels.
        %
        %   INPUTS
        %     parser   : WEC_MS2_Parser object
        %     z_target : [m]  z-coordinate of the cross-section plane
        %     n_pts    : [-]  number of points per curve segment (default 80)
        %
        %   OUTPUT
        %     polygons : {M×1 cell}  each cell is [K×2] polygon [x, y]

            if nargin < 3, n_pts = 80; end

            try
                % --- Build surface analysis (cached in persistent) ---
                persistent surf_analysis;
                persistent analysis_parser_file;

                % --- Z-level result cache (v4.2) ---
                %  Cross-section polygons depend only on the geometry
                %  (parser) and z_target.  For the same parser + z, the
                %  result is always identical.  Caching eliminates redundant
                %  evaluations during strip integration and PID iterations.
                %
                %  Cache is invalidated when the parser file changes.
                %  Key: z_target rounded to 1 µm (6 decimal places).
                %  Call WEC_Core_Functions.clear_cross_section_cache() to
                %  reset manually (e.g., if control points are modified).
                persistent z_cache;
                persistent z_cache_parser_file;

                if isempty(z_cache) || isempty(z_cache_parser_file) || ...
                        ~strcmp(z_cache_parser_file, parser.filename)
                    z_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
                    z_cache_parser_file = parser.filename;
                end

                z_key = sprintf('%.6f', z_target);
                if z_cache.isKey(z_key)
                    polygons = z_cache(z_key);
                    return;
                end

                if isempty(surf_analysis) || ~strcmp(analysis_parser_file, parser.filename)
                    surf_analysis = WEC_Core_Functions.build_surface_analysis(parser);
                    analysis_parser_file = parser.filename;
                end

                % --- Revolution fast-path DISABLED (v4.2) ---
                %  The fast-path analytical circle was incorrectly
                %  classifying z-levels as "only revolution" when DevSurf
                %  and MirrSurf(DevSurf) also contribute.  For C0, this
                %  caused the platform to be completely missed, reducing
                %  volume from 18.12 to 3.36 m³ (5.4× error).
                %
                %  The standard path below handles all surface types
                %  correctly.  With the z-level cache, performance is
                %  acceptable (each z-level is computed once).
                %
                %  TODO: Fix the z-range detection in build_surface_analysis
                %  to account for MirrSurf chains, then re-enable.

                % --- STANDARD PATH: evaluate with source caching ---
                source_cache = containers.Map();
                all_segments = {};

                for s = 1:length(parser.visible_surfs)
                    sname = parser.visible_surfs{s};
                    e = parser.entities(sname);

                    seg = [];
                    switch e.type
                        case 'RevSurf'
                            seg = WEC_Core_Functions.cross_section_revsurf( ...
                                      parser, e, z_target, n_pts, surf_analysis);
                            source_cache(sname) = seg;

                        case 'BLoftSurf'
                            seg = WEC_Core_Functions.cross_section_bloftsurf( ...
                                      parser, e, z_target, n_pts);
                            source_cache(sname) = seg;

                        case 'RuledSurf'
                            seg = WEC_Core_Functions.cross_section_ruledsurf( ...
                                      parser, e, z_target, n_pts);
                            source_cache(sname) = seg;

                        case 'DevSurf'
                            seg = WEC_Core_Functions.cross_section_devsurf( ...
                                      parser, e, z_target, n_pts);
                            source_cache(sname) = seg;

                        case 'MirrSurf'
                            % Check if source is already cached
                            src_name = e.params.source;
                            if source_cache.isKey(src_name)
                                src_seg = source_cache(src_name);
                            else
                                % Evaluate source
                                src_ent = parser.entities(src_name);
                                switch src_ent.type
                                    case 'RevSurf'
                                        src_seg = WEC_Core_Functions.cross_section_revsurf( ...
                                                      parser, src_ent, z_target, n_pts, surf_analysis);
                                    case 'BLoftSurf'
                                        src_seg = WEC_Core_Functions.cross_section_bloftsurf( ...
                                                      parser, src_ent, z_target, n_pts);
                                    case 'RuledSurf'
                                        src_seg = WEC_Core_Functions.cross_section_ruledsurf( ...
                                                      parser, src_ent, z_target, n_pts);
                                    case 'DevSurf'
                                        src_seg = WEC_Core_Functions.cross_section_devsurf( ...
                                                      parser, src_ent, z_target, n_pts);
                                    case 'MirrSurf'
                                        src_seg = WEC_Core_Functions.cross_section_mirrsurf_recursive( ...
                                                      parser, src_ent, z_target, n_pts, surf_analysis);
                                    otherwise
                                        src_seg = [];
                                end
                                source_cache(src_name) = src_seg;
                            end

                            % Flip mirror coordinate
                            if ~isempty(src_seg)
                                seg = src_seg;
                                switch e.params.mirror_plane
                                    case 'Y', seg(:, 2) = -seg(:, 2);
                                    case 'X', seg(:, 1) = -seg(:, 1);
                                end
                            end
                    end

                    if ~isempty(seg) && size(seg, 1) >= 2
                        all_segments{end+1} = seg; %#ok<AGROW>
                    end
                end

                if isempty(all_segments)
                    polygons = {};
                    z_cache(z_key) = polygons;
                    return;
                end

                polygon = WEC_Core_Functions.chain_curve_segments(all_segments);

                if isempty(polygon) || size(polygon, 1) < 3
                    polygons = {};
                else
                    polygons = {polygon(:, 1:2)};
                end

                z_cache(z_key) = polygons;

            catch ME
                warning('WEC_Core_Functions:CrossSectionFailed', ...
                        'Cross-section at z=%.3f failed: %s', z_target, ME.message);
                polygons = {};
            end
        end


        %% ─── CROSS-SECTION CACHE MANAGEMENT ──────────────────────

        function clear_cross_section_cache()
        % CLEAR_CROSS_SECTION_CACHE  Reset the z-level result cache.
        %
        %   WEC_Core_Functions.clear_cross_section_cache()
        %
        %   Call this when MS2 control points have been modified
        %   (geometry changed at the same filename).  The cache
        %   auto-invalidates when the parser filename changes, so
        %   this is only needed for same-file geometry modifications.
        %
        %   NOTE: This clears ALL persistent variables in the class,
        %   including the surface analysis cache.  Both are rebuilt
        %   automatically on the next evaluateCrossSectionMS2 call.

            clear WEC_Core_Functions;
            fprintf('  Cross-section and surface analysis caches cleared.\n');
        end


        %% ─── SURFACE ANALYSIS (computed once, cached) ────────────

        function sa = build_surface_analysis(parser)
        % BUILD_SURFACE_ANALYSIS  Pre-analyze visible surfaces for fast
        %   cross-section evaluation.
        %
        %   Identifies:
        %     - RevSurf + MirrSurf(Y=0) pairs that together cover 360°
        %     - Non-revolution surfaces and their z-ranges
        %     - Pre-sampled RevSurf profile curves

            sa = struct();
            sa.has_rev_pairs      = false;
            sa.rev_pairs          = {};        % {revsurf_name, profile_name, axis_name}
            sa.non_rev_surfaces   = {};
            sa.non_rev_z_ranges   = [];
            sa.profile_cache      = containers.Map();

            % Identify RevSurf → MirrSurf pairs
            rev_names    = {};
            mirror_names = {};
            for s = 1:length(parser.visible_surfs)
                sname = parser.visible_surfs{s};
                e = parser.entities(sname);
                if strcmp(e.type, 'RevSurf')
                    rev_names{end+1} = sname; %#ok<AGROW>
                elseif strcmp(e.type, 'MirrSurf')
                    mirror_names{end+1} = sname; %#ok<AGROW>
                end
            end

            % Check each RevSurf for a matching mirror
            rev_paired = false(length(rev_names), 1);
            mir_paired = false(length(mirror_names), 1);
            for ir = 1:length(rev_names)
                re = parser.entities(rev_names{ir});
                sweep = abs(re.params.angle_end - re.params.angle_start);
                for im = 1:length(mirror_names)
                    if mir_paired(im), continue; end
                    me = parser.entities(mirror_names{im});
                    if strcmp(me.params.source, rev_names{ir}) && ...
                       strcmp(me.params.mirror_plane, 'Y') && ...
                       abs(sweep - 180) < 1
                        % RevSurf 0→180° + MirrSurf(Y=0) = full circle
                        sa.rev_pairs{end+1} = struct( ...
                            'rev_name', rev_names{ir}, ...
                            'profile', re.params.profile, ...
                            'axis', re.params.axis);
                        rev_paired(ir) = true;
                        mir_paired(im) = true;
                        break;
                    end
                end
            end
            sa.has_rev_pairs = ~isempty(sa.rev_pairs);

            % Pre-sample profile curves for all RevSurf pairs
            for ip = 1:length(sa.rev_pairs)
                rp = sa.rev_pairs{ip};
                n_sample = 500;
                t_s = linspace(0, 1, n_sample)';
                pts = parser.eval_curve(rp.profile, t_s);
                sa.profile_cache(rp.rev_name) = struct( ...
                    't', t_s, 'pts', pts, 'z', pts(:,3));
            end

            % Identify non-revolution visible surfaces (not in any pair, not a mirror of a paired rev)
            paired_rev_names = cellfun(@(x) x.rev_name, sa.rev_pairs, 'UniformOutput', false);
            paired_mir_sources = paired_rev_names;  % mirrors of paired revs

            for s = 1:length(parser.visible_surfs)
                sname = parser.visible_surfs{s};
                e = parser.entities(sname);

                % Skip paired RevSurfs and their mirrors
                if any(strcmp(sname, paired_rev_names)), continue; end
                if strcmp(e.type, 'MirrSurf') && any(strcmp(e.params.source, paired_mir_sources))
                    continue;
                end

                sa.non_rev_surfaces{end+1} = sname;

                % Estimate z-range from boundary curves
                z_range = WEC_Core_Functions.estimate_surface_z_range(parser, sname);
                sa.non_rev_z_ranges = [sa.non_rev_z_ranges; z_range];
            end
        end


        function z_range = estimate_surface_z_range(parser, sname)
        % ESTIMATE_SURFACE_Z_RANGE  Get [z_min, z_max] of a visible surface
        %   by evaluating its boundary curves.

            e = parser.entities(sname);
            z_vals = [];
            t_s = linspace(0, 1, 50)';

            try
                switch e.type
                    case 'RevSurf'
                        pts = parser.eval_curve(e.params.profile, t_s);
                        z_vals = pts(:, 3);
                    case 'BLoftSurf'
                        for k = 1:length(e.params.section_names)
                            pts = parser.eval_curve_or_snake(e.params.section_names{k}, t_s);
                            z_vals = [z_vals; pts(:, 3)]; %#ok<AGROW>
                        end
                    case 'RuledSurf'
                        pts1 = parser.eval_curve(e.params.curve1, t_s);
                        pts2 = parser.eval_curve(e.params.curve2, t_s);
                        z_vals = [pts1(:,3); pts2(:,3)];
                    case 'DevSurf'
                        pts1 = parser.eval_snake(e.params.snake, t_s);
                        pts2 = parser.eval_curve(e.params.curve, t_s);
                        z_vals = [pts1(:,3); pts2(:,3)];
                    case 'MirrSurf'
                        z_range = WEC_Core_Functions.estimate_surface_z_range( ...
                                      parser, e.params.source);
                        return;
                end
            catch
                z_vals = [];
            end

            if isempty(z_vals)
                z_range = [-inf, inf];
            else
                z_range = [min(z_vals), max(z_vals)];
            end
        end


        function polygons = revolution_circle_at_z(parser, sa, z_target, n_pts)
        % REVOLUTION_CIRCLE_AT_Z  Analytical circle for RevSurf+mirror pair.
        %
        %   When only revolution pairs contribute at z_target, the cross-section
        %   is a perfect circle.  Returns a polygon of n_pts vertices for
        %   downstream consumers that need vertices.

            polygons = {};

            for ip = 1:length(sa.rev_pairs)
                rp = sa.rev_pairs{ip};
                pc = sa.profile_cache(rp.rev_name);

                % Root-find on cached profile
                t_roots = WEC_Core_Functions.find_all_roots_sampled( ...
                              pc.t, pc.z, z_target);

                if isempty(t_roots), continue; end

                % Evaluate at root to get radius
                profile_pt = parser.eval_curve(rp.profile, t_roots(1));
                profile_pt = profile_pt(1, :);

                % Get axis
                axis_ent   = parser.entities(rp.axis);
                axis_start = parser.eval_any_point(axis_ent.params.pt_start);
                axis_end   = parser.eval_any_point(axis_ent.params.pt_end);
                axis_vec   = axis_end - axis_start;
                axis_dir   = axis_vec / norm(axis_vec);

                v_rel   = profile_pt - axis_start;
                proj    = axis_start + dot(v_rel, axis_dir) * axis_dir;
                r       = norm(profile_pt - proj);

                if r < 1e-10, continue; end

                % Generate full circle polygon
                theta = linspace(0, 2*pi, n_pts + 1)';
                theta = theta(1:end-1);  % remove duplicate endpoint
                circle_xy = [proj(1) + r * cos(theta), ...
                             proj(2) + r * sin(theta)];
                polygons = {circle_xy};
                return;  % one circle is the whole cross-section
            end
        end


        %% ─── PER-SURFACE-TYPE CROSS-SECTION HELPERS ──────────────

        function seg = cross_section_revsurf(parser, e, z_target, n_pts, surf_analysis)
        % CROSS_SECTION_REVSURF  Cross-section of a RevSurf at z = z_target.
        %
        %   OPTIMIZED (v4.1): Uses pre-cached profile curve sampling from
        %   surf_analysis when available.  Falls back to fresh sampling.

            p = e.params;

            % Use cached profile if available
            if nargin >= 5 && ~isempty(surf_analysis) && ...
                    surf_analysis.profile_cache.isKey(parser.visible_surfs{1})
                % Try to find this RevSurf in cache
                rev_name = '';
                for ir = 1:length(surf_analysis.rev_pairs)
                    if strcmp(surf_analysis.rev_pairs{ir}.profile, p.profile)
                        rev_name = surf_analysis.rev_pairs{ir}.rev_name;
                        break;
                    end
                end
                if ~isempty(rev_name) && surf_analysis.profile_cache.isKey(rev_name)
                    pc = surf_analysis.profile_cache(rev_name);
                    t_samples   = pc.t;
                    profile_pts = pc.pts;
                    z_profile   = pc.z;
                else
                    n_sample    = 500;
                    t_samples   = linspace(0, 1, n_sample)';
                    profile_pts = parser.eval_curve(p.profile, t_samples);
                    z_profile   = profile_pts(:, 3);
                end
            else
                n_sample    = 500;
                t_samples   = linspace(0, 1, n_sample)';
                profile_pts = parser.eval_curve(p.profile, t_samples);
                z_profile   = profile_pts(:, 3);
            end

            % Find all crossings where z_profile = z_target
            t_roots = WEC_Core_Functions.find_all_roots_sampled( ...
                          t_samples, z_profile, z_target);

            if isempty(t_roots)
                seg = [];
                return;
            end

            % Evaluate axis
            axis_ent   = parser.entities(p.axis);
            axis_start = parser.eval_any_point(axis_ent.params.pt_start);
            axis_end   = parser.eval_any_point(axis_ent.params.pt_end);
            axis_vec   = axis_end - axis_start;
            axis_len   = norm(axis_vec);
            if axis_len < 1e-12, seg = []; return; end
            axis_dir   = axis_vec / axis_len;

            % Use the first valid crossing (typical for single-valued profiles)
            % For multi-valued profiles (e.g. ellipsoid equator), multiple
            % crossings produce multiple arcs — handle if needed.
            seg = zeros(0, 3);
            for ic = 1:length(t_roots)
                t_star = t_roots(ic);
                profile_pt = parser.eval_curve(p.profile, t_star);
                profile_pt = profile_pt(1, :);

                % Project onto axis → get radial distance
                v_rel   = profile_pt - axis_start;
                z_along = dot(v_rel, axis_dir);
                proj    = axis_start + z_along * axis_dir;
                radial  = profile_pt - proj;
                r       = norm(radial);

                if r < 1e-10
                    % Point on axis — degenerate (tip of hull)
                    continue;
                end

                e_r = radial / r;
                e_t = cross(axis_dir, e_r);
                e_t = e_t / norm(e_t);

                % Generate arc from angle_start to angle_end
                angles = linspace(deg2rad(p.angle_start), ...
                                  deg2rad(p.angle_end), n_pts)';
                arc = proj + r * cos(angles) .* e_r + ...
                             r * sin(angles) .* e_t;
                seg = [seg; arc]; %#ok<AGROW>
            end
        end


        function seg = cross_section_bloftsurf(parser, e, z_target, n_pts)
        % CROSS_SECTION_BLOFTSURF  Cross-section of a BLoftSurf at z = z_target.
        %
        %   TWO LOFT ORIENTATIONS — detected automatically:
        %
        %   (A) AXIAL LOFT  (z varies with v, u moves along height):
        %     Sections are stacked at different heights in z.
        %     The original algorithm applies: for each u, root-find v*
        %     where z(v*) = z_target, then collect 3D points.
        %     Example: old-style BLoftSurf with CopyCurve sections at
        %     different z elevations.
        %
        %   (B) AZIMUTHAL LOFT  (z constant in v, varies with u):
        %     Sections wrap azimuthally at the same height.  Height is
        %     controlled by u, not v.  The v root-finder sees a constant
        %     signal and returns no roots — producing an empty segment
        %     for every u, even when z(u) = z_target.  This is the bug
        %     that caused the new E1.ms2 Wall column to vanish.
        %
        %     FIX: detect flat-z sections, then root-find in u instead.
        %     At u* where z(u*) = z_target, sweep all v values to produce
        %     the full horizontal arc (the correct cross-section).
        %     Example: new E1.ms2 Wall — snake3/snake1/snake2 all at the
        %     same z for any given u; revolution preserves z exactly.
        %
        %   DETECTION:
        %     Probe 5 u-values.  If max(z_ctrl)-min(z_ctrl) < flat_z_tol
        %     at every probe, classify as azimuthal.  The two branches are
        %     then applied globally for that BLoftSurf entity.
        %
        %   PROOF OF BUG (azimuthal case):
        %     snake1(u) = curve3(u) = [x(u), 1.0, z(u)]  — y-shifted copy
        %     snake2(u) = curve2(u) = [x(u), 0.0, z(u)]  — same z
        %     snake3(u) = RevSurf(snake1,line2,270°)      — revolution
        %                 preserves z → [0, 1±r(u), z(u)]
        %     ∴ z_ctrl = [z(u), z(u), z(u)] for ALL u → z_bspline = const
        %     → find_all_roots_sampled sees no sign changes → v_roots = []
        %     → seg = [] for every z_target → Wall absent from ALL plots.

            p = e.params;
            section_names = p.section_names;
            n_sec  = length(section_names);
            degree = p.degree;
            knots  = WEC_MS2_Parser.make_clamped_knots(n_sec, degree);

            u_samples = linspace(0, 1, n_pts)';
            seg = zeros(0, 3);

            % ── Classify loft orientation (5 probe points) ────────────
            %  flat_z_tol: sections are co-planar in z if z-spread < this.
            %  1e-4 m matches MultiSurf coordinate precision and is well
            %  below any physical z-variation in an axial loft.
            flat_z_tol = 1e-4;
            is_azimuthal = true;
            for ip = 1:5
                u_p = (ip - 1) / 4;
                sec_z_probe = zeros(n_sec, 1);
                for k = 1:n_sec
                    pk = parser.eval_curve_or_snake(section_names{k}, u_p);
                    sec_z_probe(k) = pk(1, 3);
                end
                if (max(sec_z_probe) - min(sec_z_probe)) > flat_z_tol
                    is_azimuthal = false;
                    break;
                end
            end

            % ── BRANCH A: azimuthal loft ──────────────────────────────
            %  z varies with u, not v.
            %  Step 1: collect z_mean at every u (one section eval per u).
            %  Step 2: root-find u* where z(u*) = z_target.
            %  Step 3: at each u*, sweep all v to produce the horizontal arc.
            if is_azimuthal
                z_at_u = zeros(n_pts, 1);
                for iu = 1:n_pts
                    pk = parser.eval_curve_or_snake( ...
                             section_names{1}, u_samples(iu));
                    z_at_u(iu) = pk(1, 3);
                end

                u_roots = WEC_Core_Functions.find_all_roots_sampled( ...
                              u_samples, z_at_u, z_target);

                if isempty(u_roots), return; end

                n_v_arc = 40;
                v_arc   = linspace(0, 1, n_v_arc)';

                for ir = 1:length(u_roots)
                    u_star = u_roots(ir);

                    % Evaluate all sections at u* once
                    sec_pts = zeros(n_sec, 3);
                    for k = 1:n_sec
                        pk = parser.eval_curve_or_snake( ...
                                 section_names{k}, u_star);
                        sec_pts(k, :) = pk(1, :);
                    end

                    % Output the full horizontal arc (sweep v)
                    for iv = 1:n_v_arc
                        pt_3d = WEC_MS2_Parser.bspline_curve_eval( ...
                                    knots, sec_pts, degree, v_arc(iv));
                        seg(end+1, :) = pt_3d; %#ok<AGROW>
                    end
                end
                return;
            end

            % ── BRANCH B: axial loft (original algorithm) ─────────────
            %  z varies with v.  For each u, root-find v* where z(v*)=z_target.
            for iu = 1:length(u_samples)
                u = u_samples(iu);

                % Evaluate all section curves at this u
                sec_pts = zeros(n_sec, 3);
                for k = 1:n_sec
                    pts_k = parser.eval_curve_or_snake(section_names{k}, u);
                    sec_pts(k, :) = pts_k(1, :);
                end

                % z-component control values for the v-direction B-spline
                z_ctrl = sec_pts(:, 3);

                % Check if z_target is reachable (within the B-spline hull)
                if z_target < min(z_ctrl) - 0.5 || z_target > max(z_ctrl) + 0.5
                    continue;
                end

                % Root-find: z(v*) = z_target on the 1D B-spline
                n_v_sample = 200;
                v_samples  = linspace(0, 1, n_v_sample)';
                z_bspline  = zeros(n_v_sample, 1);
                for iv = 1:n_v_sample
                    pt_v = WEC_MS2_Parser.bspline_curve_eval( ...
                               knots, sec_pts, degree, v_samples(iv));
                    z_bspline(iv) = pt_v(3);
                end

                v_roots = WEC_Core_Functions.find_all_roots_sampled( ...
                              v_samples, z_bspline, z_target);

                if isempty(v_roots), continue; end

                % Evaluate full 3D point at (u, v*)
                for ir = 1:length(v_roots)
                    pt_3d = WEC_MS2_Parser.bspline_curve_eval( ...
                                knots, sec_pts, degree, v_roots(ir));
                    seg(end+1, :) = pt_3d; %#ok<AGROW>
                end
            end
        end


        function seg = cross_section_ruledsurf(parser, e, z_target, n_pts)
        % CROSS_SECTION_RULEDSURF  Cross-section of a RuledSurf at z = z_target.
        %
        %   S(u,v) = (1−v)·C₁(u) + v·C₂(u)
        %   z(u,v) = (1−v)·z₁(u) + v·z₂(u)
        %   → v* = (z_target − z₁(u)) / (z₂(u) − z₁(u))

            p = e.params;
            u_samples = linspace(0, 1, n_pts)';

            pts1 = parser.eval_curve(p.curve1, u_samples);
            pts2 = parser.eval_curve(p.curve2, u_samples);
            z1   = pts1(:, 3);
            z2   = pts2(:, 3);

            seg = zeros(0, 3);
            for iu = 1:n_pts
                dz = z2(iu) - z1(iu);
                if abs(dz) < 1e-12
                    % Both curves at same z — check if it matches target
                    if abs(z1(iu) - z_target) < 1e-6
                        v_star = 0.5;
                    else
                        continue;
                    end
                else
                    v_star = (z_target - z1(iu)) / dz;
                end

                if v_star < -1e-6 || v_star > 1 + 1e-6
                    continue;
                end
                v_star = max(0, min(1, v_star));

                pt = (1 - v_star) * pts1(iu, :) + v_star * pts2(iu, :);
                seg(end+1, :) = pt; %#ok<AGROW>
            end
        end


        function seg = cross_section_devsurf(parser, e, z_target, n_pts)
        % CROSS_SECTION_DEVSURF  Cross-section of a DevSurf at z = z_target.
        %
        %   DevSurf is a ruled surface between a snake and a curve:
        %     S(u,v) = (1−v)·snake(u) + v·curve(u)
        %   Same algorithm as RuledSurf.

            p = e.params;
            u_samples = linspace(0, 1, n_pts)';

            pts_snake = parser.eval_snake(p.snake, u_samples);
            pts_curve = parser.eval_curve(p.curve, u_samples);
            z1 = pts_snake(:, 3);
            z2 = pts_curve(:, 3);

            seg = zeros(0, 3);
            for iu = 1:n_pts
                dz = z2(iu) - z1(iu);
                if abs(dz) < 1e-12
                    if abs(z1(iu) - z_target) < 1e-6
                        v_star = 0.5;
                    else
                        continue;
                    end
                else
                    v_star = (z_target - z1(iu)) / dz;
                end

                if v_star < -1e-6 || v_star > 1 + 1e-6
                    continue;
                end
                v_star = max(0, min(1, v_star));

                pt = (1 - v_star) * pts_snake(iu, :) + v_star * pts_curve(iu, :);
                seg(end+1, :) = pt; %#ok<AGROW>
            end
        end


        function seg = cross_section_mirrsurf(parser, e, z_target, n_pts)
        % CROSS_SECTION_MIRRSURF  Cross-section of a MirrSurf at z = z_target.
        %
        %   Evaluate source surface cross-section, flip the mirror coordinate.
        %   Handles chained MirrSurf (e.g. surface7 → surface5 → surface4).

            p = e.params;
            source_ent = parser.entities(p.source);

            % Dispatch to source surface type
            switch source_ent.type
                case 'RevSurf'
                    seg = WEC_Core_Functions.cross_section_revsurf( ...
                              parser, source_ent, z_target, n_pts);
                case 'BLoftSurf'
                    seg = WEC_Core_Functions.cross_section_bloftsurf( ...
                              parser, source_ent, z_target, n_pts);
                case 'RuledSurf'
                    seg = WEC_Core_Functions.cross_section_ruledsurf( ...
                              parser, source_ent, z_target, n_pts);
                case 'DevSurf'
                    seg = WEC_Core_Functions.cross_section_devsurf( ...
                              parser, source_ent, z_target, n_pts);
                case 'MirrSurf'
                    seg = WEC_Core_Functions.cross_section_mirrsurf_recursive( ...
                              parser, source_ent, z_target, n_pts);
                otherwise
                    seg = [];
                    return;
            end

            if isempty(seg), return; end

            % Apply mirror flip
            switch p.mirror_plane
                case 'Y'
                    seg(:, 2) = -seg(:, 2);
                case 'X'
                    seg(:, 1) = -seg(:, 1);
            end
        end


        function seg = cross_section_mirrsurf_recursive(parser, e, z_target, n_pts, surf_analysis)
        % CROSS_SECTION_MIRRSURF_RECURSIVE  Recursively resolve chained MirrSurf.
        %
        %   C0 topology: surface7 → MirrSurf(surface5) → MirrSurf(surface4, DevSurf)
        %   Follows the chain until reaching a primary surface type,
        %   evaluates it, then applies mirror flips in reverse order.

            if nargin < 5, surf_analysis = []; end

            p = e.params;
            source_ent = parser.entities(p.source);

            switch source_ent.type
                case 'RevSurf'
                    if ~isempty(surf_analysis)
                        seg = WEC_Core_Functions.cross_section_revsurf( ...
                                  parser, source_ent, z_target, n_pts, surf_analysis);
                    else
                        seg = WEC_Core_Functions.cross_section_revsurf( ...
                                  parser, source_ent, z_target, n_pts);
                    end
                case 'BLoftSurf'
                    seg = WEC_Core_Functions.cross_section_bloftsurf( ...
                              parser, source_ent, z_target, n_pts);
                case 'RuledSurf'
                    seg = WEC_Core_Functions.cross_section_ruledsurf( ...
                              parser, source_ent, z_target, n_pts);
                case 'DevSurf'
                    seg = WEC_Core_Functions.cross_section_devsurf( ...
                              parser, source_ent, z_target, n_pts);
                case 'MirrSurf'
                    seg = WEC_Core_Functions.cross_section_mirrsurf_recursive( ...
                              parser, source_ent, z_target, n_pts, surf_analysis);
                otherwise
                    seg = [];
                    return;
            end

            if isempty(seg), return; end

            switch p.mirror_plane
                case 'Y', seg(:, 2) = -seg(:, 2);
                case 'X', seg(:, 1) = -seg(:, 1);
            end
        end


        %% ─── CURVE SEGMENT CHAINING ──────────────────────────────

        function polygon = chain_curve_segments(segments)
        % CHAIN_CURVE_SEGMENTS  Chain multiple open curve segments into a
        %   closed polygon by nearest-endpoint matching.
        %
        %   INPUTS
        %     segments : cell array, each cell is [N_k × 3] curve points
        %
        %   OUTPUT
        %     polygon  : [M × 3] ordered closed polygon
        %
        %   WHY 3D output?
        %     The chaining works in 3D (matching [x,y,z] endpoints).
        %     The caller extracts [x,y] columns for polygon operations.

            if isempty(segments)
                polygon = [];
                return;
            end

            n_segs = length(segments);
            used   = false(n_segs, 1);

            % Start with the first segment
            polygon = segments{1};
            used(1) = true;

            tol = 0.5;  % [m] matching tolerance (conservative)

            for iter = 1:n_segs - 1
                cur_end = polygon(end, :);
                best_d   = inf;
                best_idx = 0;
                flip     = false;

                for s = 1:n_segs
                    if used(s), continue; end

                    seg = segments{s};
                    d_start = norm(cur_end - seg(1, :));
                    d_end   = norm(cur_end - seg(end, :));

                    if d_start < best_d
                        best_d   = d_start;
                        best_idx = s;
                        flip     = false;
                    end
                    if d_end < best_d
                        best_d   = d_end;
                        best_idx = s;
                        flip     = true;
                    end
                end

                if best_d > tol || best_idx == 0
                    break;
                end

                used(best_idx) = true;
                seg = segments{best_idx};
                if flip
                    seg = seg(end:-1:1, :);
                end
                polygon = [polygon; seg]; %#ok<AGROW>
            end

            % Close if endpoints are near
            if size(polygon, 1) > 2 && ...
               norm(polygon(1,:) - polygon(end,:)) < tol
                polygon(end+1, :) = polygon(1, :);
            end
        end


        %% ─── ROOT-FINDING ON SAMPLED CURVES ─────────────────────

        function roots = find_all_roots_sampled(t_samples, values, target)
        % FIND_ALL_ROOTS_SAMPLED  Find all t-values where a sampled function
        %   crosses a target value.  Uses linear interpolation between samples
        %   followed by bisection refinement.
        %
        %   INPUTS
        %     t_samples : [N×1] parameter values (monotonically increasing)
        %     values    : [N×1] function values f(t)
        %     target    : scalar target value
        %
        %   OUTPUT
        %     roots     : [M×1] refined t-values where f(t) ≈ target

            residual = values - target;
            n = length(residual);
            roots = [];

            for i = 1:n-1
                if residual(i) * residual(i+1) < 0
                    % Sign change — bracket found
                    t_lo = t_samples(i);
                    t_hi = t_samples(i+1);

                    % Bisection refinement (20 iterations → ~1e-6 precision)
                    for iter = 1:20
                        t_mid = 0.5 * (t_lo + t_hi);
                        f_mid = interp1(t_samples, values, t_mid, 'linear');
                        if (f_mid - target) * (interp1(t_samples, values, t_lo, 'linear') - target) < 0
                            t_hi = t_mid;
                        else
                            t_lo = t_mid;
                        end
                    end
                    roots(end+1) = 0.5 * (t_lo + t_hi); %#ok<AGROW>

                elseif abs(residual(i)) < 1e-8
                    roots(end+1) = t_samples(i); %#ok<AGROW>
                end
            end
        end


        %% ═════════════════════════════════════════════════════════════
        %%  MIDPLANE PROFILE EXTRACTION
        %% ═════════════════════════════════════════════════════════════

        function profile = extractProfileMS2(parser)
        % EXTRACTPROFILEMS2  Extract the xz-midplane profile from the parser.
        %
        %   profile = WEC_Core_Functions.extractProfileMS2(parser)
        %
        %   The profile is the hull's outer boundary projected onto the
        %   xz-plane (y = 0).  For symmetric hulls, this is the silhouette
        %   at the symmetry plane.
        %
        %   ALGORITHM
        %     For each visible surface, evaluate boundary edges and
        %     collect points near y = 0.  The profile is the union of
        %     these boundary points, sorted by angle from centroid.
        %
        %   WHY boundary edges?
        %     The profile lies on the hull boundary in the symmetry plane.
        %     For RevSurf, this is the profile curve at the start/end angle.
        %     For BLoftSurf, this is an edge where the loft meets the
        %     symmetry plane (typically u = 0 or v = 0/1).
        %     For MirrSurf (Y=0), the junction between source and mirror
        %     lies at y = 0 — this is an edge of the source surface.
        %
        %   OUTPUT
        %     profile : [M × 2]  polygon vertices [x, z], angle-sorted

            try
                n_edge_pts = 200;
                t_edge     = linspace(0, 1, n_edge_pts)';
                all_xz     = zeros(0, 2);

                for s = 1:length(parser.visible_surfs)
                    sname = parser.visible_surfs{s};
                    e = parser.entities(sname);

                    switch e.type
                        case 'RevSurf'
                            % Profile curve IS the y=0 silhouette (at start angle).
                            pts = parser.eval_curve(e.params.profile, t_edge);
                            all_xz = [all_xz; pts(:, [1, 3])]; %#ok<AGROW>

                        case 'BLoftSurf'
                            % Evaluate edges at u=0 and u=1 (boundary curves)
                            for ie = [1, 3]  % edges 1 and 3: v=0,1 with u varying
                                for it = 1:n_edge_pts
                                    if ie == 1
                                        pt = parser.eval_surface(sname, t_edge(it), 0);
                                    else
                                        pt = parser.eval_surface(sname, t_edge(it), 1);
                                    end
                                    if abs(pt(2)) < 0.05  % near y=0
                                        all_xz(end+1, :) = [pt(1), pt(3)]; %#ok<AGROW>
                                    end
                                end
                            end
                            % Also check edges at v=0 and v=1
                            for ie = [2, 4]
                                for it = 1:n_edge_pts
                                    if ie == 2
                                        pt = parser.eval_surface(sname, 1, t_edge(it));
                                    else
                                        pt = parser.eval_surface(sname, 0, t_edge(it));
                                    end
                                    if abs(pt(2)) < 0.05
                                        all_xz(end+1, :) = [pt(1), pt(3)]; %#ok<AGROW>
                                    end
                                end
                            end

                        case 'RuledSurf'
                            pts1 = parser.eval_curve(e.params.curve1, t_edge);
                            pts2 = parser.eval_curve(e.params.curve2, t_edge);
                            near1 = abs(pts1(:,2)) < 0.05;
                            near2 = abs(pts2(:,2)) < 0.05;
                            all_xz = [all_xz; pts1(near1, [1,3]); pts2(near2, [1,3])]; %#ok<AGROW>

                        case 'DevSurf'
                            pts_s = parser.eval_snake(e.params.snake, t_edge);
                            pts_c = parser.eval_curve(e.params.curve, t_edge);
                            near_s = abs(pts_s(:,2)) < 0.05;
                            near_c = abs(pts_c(:,2)) < 0.05;
                            all_xz = [all_xz; pts_s(near_s, [1,3]); pts_c(near_c, [1,3])]; %#ok<AGROW>

                        case 'MirrSurf'
                            % MirrSurf about Y=0: the source boundary AT y=0
                            % is the profile.  Already captured by the source
                            % surface's edge evaluation.  Skip to avoid
                            % duplicating points.
                            continue;
                    end
                end

                if size(all_xz, 1) < 3
                    warning('WEC_Core_Functions:ProfileEmpty', ...
                            'Profile extraction found < 3 points');
                    profile = all_xz;
                    return;
                end

                % Mirror across x = 0 to get the full symmetric profile.
                % The source surfaces (RevSurf, DevSurf, etc.) produce
                % x >= 0 points only.  MirrSurf(Y=0) creates the x < 0
                % half of the 3D hull, but doesn't appear in the y=0
                % silhouette because MirrSurf flips y, not x.  For the
                % XZ profile view, bilateral symmetry means the left half
                % is a mirror of the right half about x = 0.
                mirrored = all_xz;
                mirrored(:, 1) = -mirrored(:, 1);
                % Remove mirrored points that are at x ≈ 0 (would duplicate)
                not_on_axis = abs(mirrored(:, 1)) > 1e-6;
                all_xz = [all_xz; mirrored(not_on_axis, :)];

                % Remove duplicates
                all_xz = unique(round(all_xz * 1e6) / 1e6, 'rows', 'stable');

                % Sort by angle from centroid to form proper polygon
                cx = mean(all_xz(:, 1));
                cz = mean(all_xz(:, 2));
                angles = atan2(all_xz(:, 2) - cz, all_xz(:, 1) - cx);
                [~, order] = sort(angles);
                profile = all_xz(order, :);

            catch ME
                warning('WEC_Core_Functions:ProfileFailed', ...
                        'Profile extraction failed: %s', ME.message);
                profile = [0 0; 1 0; 0.5 1];
            end
        end


        %% ═════════════════════════════════════════════════════════════
        %%  2D POLYGON OPERATIONS (UNCHANGED)
        %% ═════════════════════════════════════════════════════════════

        function poly_out = clipPolygon(poly, z_plane, keep_side)
        % CLIPPOLYGON  Sutherland–Hodgman clip of a 2D polygon against a
        %   horizontal z-plane.
        %
        %   poly_out = CLIPPOLYGON(poly, z_plane, keep_side)
        %
        %   INPUTS
        %     poly      : [N×2]  polygon vertices [x, z]
        %     z_plane   : [m]    z-coordinate of the clipping plane
        %     keep_side : 'below' | 'above'  which half-space to keep
        %
        %   OUTPUT
        %     poly_out  : [M×2]  clipped polygon (empty if fully clipped)
            try
                if isempty(poly) || size(poly, 1) < 2
                    poly_out = [];
                    return;
                end

                n = size(poly, 1);
                output = zeros(2*n, 2);
                out_count = 0;
                epsilon = 1e-9;

                for i = 1:n
                    p1 = poly(i, :);
                    next_idx = mod(i, n) + 1;
                    p2 = poly(next_idx, :);

                    if strcmp(keep_side, 'below')
                        p1_in = p1(2) <= z_plane + epsilon;
                        p2_in = p2(2) <= z_plane + epsilon;
                    else
                        p1_in = p1(2) >= z_plane - epsilon;
                        p2_in = p2(2) >= z_plane - epsilon;
                    end

                    if p1_in
                        out_count = out_count + 1;
                        output(out_count, :) = p1;
                    end

                    if xor(p1_in, p2_in)
                        dz = p2(2) - p1(2);
                        if abs(dz) > epsilon
                            t = (z_plane - p1(2)) / dz;
                            t = max(0, min(1, t));
                            intersection_point = p1 + t * (p2 - p1);
                            out_count = out_count + 1;
                            output(out_count, :) = intersection_point;
                        end
                    end
                end

                if out_count > 2
                    poly_out_temp = output(1:out_count, :);
                    [~, idx] = unique(round(poly_out_temp*1e6)/1e6, 'rows', 'stable');
                    poly_out = poly_out_temp(sort(idx), :);
                else
                    poly_out = [];
                end

            catch
                poly_out = [];
            end
        end

        function [geom, iner, cpmo] = polygeom(x, y)
            % POLYGEOM - Calculate 2D polygon geometric properties
            %
            % INPUTS:
            %   x, y: [Nx1] polygon vertex coordinates
            %
            % OUTPUTS:
            %   geom: [Area, xc, yc, 0] - area and centroid
            %   iner: [Ixx, Iyy, Ixy, 0, 0, 0] - moments about origin
            %   cpmo: [Iuu, Ivv, Iuv, 0, 0, 0] - moments about centroid

            try
                if nargin < 2
                    error('polygeom requires x and y coordinates');
                end

                x = x(:);
                y = y(:);

                % Close polygon if not already closed
                if (x(1) ~= x(end)) || (y(1) ~= y(end))
                    x(end+1) = x(1);
                    y(end+1) = y(1);
                end

                n = length(x);

                % Require at least 3 unique vertices
                if n < 4
                    geom = [0, 0, 0, 0];
                    iner = [0, 0, 0, 0, 0, 0];
                    cpmo = [0, 0, 0, 0, 0, 0];
                    return;
                end

                xi = x(1:n-1);
                yi = y(1:n-1);
                xip1 = x(2:n);
                yip1 = y(2:n);

                % Shoelace formula
                a = xi .* yip1 - xip1 .* yi;
                Area = 0.5 * sum(a);

                if abs(Area) < 1e-12
                    xc = mean(xi);
                    yc = mean(yi);
                    Ixx = 0; Iyy = 0; Ixy = 0;
                    Iuu = 0; Ivv = 0; Iuv = 0;
                else
                    xc = sum((xi + xip1) .* a) / (6 * Area);
                    yc = sum((yi + yip1) .* a) / (6 * Area);

                    Ixx = sum((yi.^2 + yi.*yip1 + yip1.^2) .* a) / 12;
                    Iyy = sum((xi.^2 + xi.*xip1 + xip1.^2) .* a) / 12;
                    Ixy = sum((xi.*yip1 + 2*xi.*yi + 2*xip1.*yip1 + xip1.*yi) .* a) / 24;

                    Iuu = Ixx - Area * yc^2;
                    Ivv = Iyy - Area * xc^2;
                    Iuv = Ixy - Area * xc * yc;
                end

                geom = [abs(Area), xc, yc, 0];
                iner = [Ixx, Iyy, Ixy, 0, 0, 0];
                cpmo = [abs(Iuu), abs(Ivv), Iuv, 0, 0, 0];

            catch ME
                warning('WEC_Core_Functions:PolygeomFailed', ...
                        'Polygon geometry calculation failed: %s', ME.message);
                geom = [0, 0, 0, 0];
                iner = [0, 0, 0, 0, 0, 0];
                cpmo = [0, 0, 0, 0, 0, 0];
            end
        end

        function intersections = find_waterline_intersections(profile, z_level)
            % FIND_WATERLINE_INTERSECTIONS - Find points where polygon crosses horizontal line
            %
            % INPUTS:
            %   profile: [Nx2] polygon vertices [x, z]
            %   z_level: Scalar z-coordinate of waterline
            %
            % OUTPUT:
            %   intersections: [Mx2] intersection points [x, z_level]

            try
                intersections = [];

                if isempty(profile) || size(profile, 1) < 2
                    return;
                end

                if ~isscalar(z_level) || ~isnumeric(z_level)
                    return;
                end

                n_vertices = size(profile, 1);
                epsilon = 1e-10;

                for i = 1:n_vertices
                    p1 = profile(i, :);
                    p2_idx = mod(i, n_vertices) + 1;
                    p2 = profile(p2_idx, :);

                    z1 = p1(2);
                    z2 = p2(2);

                    if (z1 - z_level) * (z2 - z_level) < -epsilon
                        dz = z2 - z1;
                        if abs(dz) > epsilon
                            t = (z_level - z1) / dz;
                            t = max(0, min(1, t));
                            x_intersect = p1(1) + t * (p2(1) - p1(1));
                            intersections = [intersections; x_intersect, z_level]; %#ok<AGROW>
                        end
                    end
                end

            catch
                intersections = [];
            end
        end

        function [total_area, total_Ixx, total_Iyy] = calculatePolygonProperties(polygons)
            % CALCULATEPOLYGONPROPERTIES - Calculate aggregate properties of multiple polygons
            %
            % INPUT:
            %   polygons: Cell array of [Nx2] polygon vertices
            %
            % OUTPUTS:
            %   total_area, total_Ixx, total_Iyy: Summed properties

            try
                total_area = 0;
                total_Ixx = 0;
                total_Iyy = 0;

                for i = 1:length(polygons)
                    poly = polygons{i};

                    if size(poly,1) < 3
                        continue;
                    end

                    % Close polygon if not already
                    if any(poly(1,:) ~= poly(end,:))
                        poly(end+1,:) = poly(1,:); %#ok<AGROW>
                    end

                    x = poly(:,1);
                    y = poly(:,2);

                    % Shoelace formula
                    a = x(1:end-1).*y(2:end) - x(2:end).*y(1:end-1);
                    area = 0.5 * sum(a);

                    % Second moments
                    Ixx = sum((y(1:end-1).^2 + y(1:end-1).*y(2:end) + y(2:end).^2).*a)/12;
                    Iyy = sum((x(1:end-1).^2 + x(1:end-1).*x(2:end) + x(2:end).^2).*a)/12;

                    total_area = total_area + abs(area);
                    total_Ixx = total_Ixx + abs(Ixx);
                    total_Iyy = total_Iyy + abs(Iyy);
                end

            catch ME
                warning('WEC_Core_Functions:PolygonPropertiesFailed', ...
                        'Polygon properties calculation failed: %s', ME.message);
                total_area = 0;
                total_Ixx = 0;
                total_Iyy = 0;
            end
        end


        %% ═════════════════════════════════════════════════════════════
        %%  HYDRODYNAMIC INTERPOLATION (UNCHANGED)
        %% ═════════════════════════════════════════════════════════════

        function [A11, A33, A55, K_pto, A_full, B_full] = ...
                interpolate_wamit_added_mass(vertical_shift, config, z_cg_target)
            % INTERPOLATE_WAMIT_ADDED_MASS  Vertical-shift-interpolated added mass & damping.
            %
            %   [A11, A33, A55, K_pto, A_full, B_full] = ...
            %       interpolate_wamit_added_mass(vs, config)
            %   [A11, A33, A55, K_pto, A_full, B_full] = ...
            %       interpolate_wamit_added_mass(vs, config, z_cg_target)
            %
            %   Linearly interpolates between pre-computed WAMIT cases
            %   stored in config (sorted by vertical_shift in
            %   WEC_Configuration_Builder).  Clamps to data range when
            %   vertical_shift is outside the WAMIT grid.
            %
            %   REFERENCE-FRAME CONVENTION
            %     config.wamit_A_full{i} is stored at the CG that was used
            %     to compute the HAMS hydrostatics for cache entry i —
            %     i.e., config.wamit_z_cg(i), the uniform-density CG at
            %     that draft.  This is referred to as "HAMS-CG".
            %
            %     If z_cg_target is provided (and non-empty), the function
            %     applies a single delta congruence transform from HAMS-CG
            %     to z_cg_target:
            %       T  = [1 0 −dz; 0 1 0; 0 0 1]     dz = z_cg_target − z_cg_hams
            %       M_target = T' · M_hams_cg · T   (A55, A15, A51 corrected)
            %     A11, A33 are translation-invariant; only A55 (and the
            %     surge-pitch cross terms in A_full / B_full) change.
            %
            %     If z_cg_target is empty or omitted, values are returned
            %     at HAMS-CG (legacy behaviour, suitable for diagnostic /
            %     visualisation consumers that don't tie A to a specific
            %     mass distribution).
            %
            %   This is the SINGLE point where the delta transform is
            %   applied.  Callers used to inline this logic in three
            %   places (calculate_3d_properties §9b, WEC_Shell_Offset.
            %   build_realised_props, and a missing site in evaluate_at)
            %   — those duplicates have been removed in favour of passing
            %   z_cg_target here.

            if nargin < 3, z_cg_target = []; end

            try
                if ~isfield(config, 'wamit_drafts') || ~isfield(config, 'wamit_A')
                    error('WAMIT data not found in config');
                end

                drafts = config.wamit_drafts;
                A_data = config.wamit_A;

                if isfield(config, 'wamit_K_pto')
                    K_data = config.wamit_K_pto;
                else
                    K_data = ones(size(drafts)) * 1e6;
                end

                % Clamp to data range
                draft_clamped = max(min(vertical_shift, max(drafts)), min(drafts));

                % Handle single-draft case
                if length(drafts) == 1
                    A11 = A_data(1, 1);
                    A33 = A_data(1, 2);
                    A55 = A_data(1, 3);
                    K_pto = K_data(1);

                    if isfield(config, 'wamit_A_full') && ~isempty(config.wamit_A_full)
                        A_full = config.wamit_A_full{1};
                    else
                        A_full = diag([A11, A33, A55]);
                    end

                    if isfield(config, 'wamit_B_full') && ~isempty(config.wamit_B_full)
                        B_full = config.wamit_B_full{1};
                    else
                        B_full = zeros(3, 3);
                    end
                else
                    A11 = interp1(drafts, A_data(:,1), draft_clamped, 'linear', 'extrap');
                    A33 = interp1(drafts, A_data(:,2), draft_clamped, 'linear', 'extrap');
                    A55 = interp1(drafts, A_data(:,3), draft_clamped, 'linear', 'extrap');
                    K_pto = interp1(drafts, K_data, draft_clamped, 'linear', 'extrap');

                    if isfield(config, 'wamit_A_full') && ~isempty(config.wamit_A_full)
                        A_full = WEC_Core_Functions.interpolate_3x3_matrix(...
                            drafts, config.wamit_A_full, draft_clamped);
                    else
                        A_full = diag([A11, A33, A55]);
                    end

                    if isfield(config, 'wamit_B_full') && ~isempty(config.wamit_B_full)
                        B_full = WEC_Core_Functions.interpolate_3x3_matrix(...
                            drafts, config.wamit_B_full, draft_clamped);
                    else
                        B_full = zeros(3, 3);
                    end
                end

                % Ensure non-negative diagonal terms
                A11 = max(A11, 0);
                A33 = max(A33, 0);
                A55 = max(A55, 0);
                K_pto = max(K_pto, 0);

                for ii = 1:3
                    A_full(ii, ii) = max(A_full(ii, ii), 0);
                    B_full(ii, ii) = max(B_full(ii, ii), 0);
                end

                %% --- Optional delta transform: HAMS-CG → z_cg_target ---
                %  Centralized congruence transform.  Mathematically
                %  equivalent to (and replaces) the inline blocks that
                %  used to live in calculate_3d_properties §9b and in
                %  WEC_Shell_Offset.build_realised_props.
                if ~isempty(z_cg_target) && ...
                        isfield(config, 'wamit_z_cg') && ...
                        ~isempty(config.wamit_z_cg)
                    if length(drafts) == 1
                        z_cg_hams = config.wamit_z_cg(1);
                    else
                        z_cg_hams = interp1(drafts, config.wamit_z_cg, ...
                                            draft_clamped, 'linear', 'extrap');
                    end
                    dz = z_cg_target - z_cg_hams;
                    if abs(dz) > 1e-4
                        % WEC_File_IO.transform_hydrodynamic_matrices_3x3
                        % applies M_target = T'·M·T where T encodes the
                        % CG shift specified by its third argument.
                        % Composes correctly: applying T(z_cg_hams) (during
                        % rebuild_config_hydro) followed by T(dz) here is
                        % algebraically identical to T(z_cg_target) applied
                        % to the origin-frame matrices, so the net result
                        % is M at z_cg_target irrespective of the HAMS-CG
                        % intermediate.
                        [A_full, B_full] = ...
                            WEC_File_IO.transform_hydrodynamic_matrices_3x3( ...
                                A_full, B_full, dz);
                        % Re-extract scalars from the transformed A_full.
                        % A11 / A33 are translation-invariant (recompute for
                        % safety against numerical noise from the congruence
                        % round-trip), A55 picks up the real correction.
                        A11 = max(0, A_full(1, 1));
                        A33 = max(0, A_full(2, 2));
                        A55 = max(0, A_full(3, 3));
                    end
                end

            catch ME
                warning('WEC_Core_Functions:WamitInterpolationFailed', ...
                        'WAMIT interpolation failed: %s. Returning zeros.', ME.message);
                A11 = 0; A33 = 0; A55 = 0; K_pto = 1e6;
                A_full = zeros(3, 3);
                B_full = zeros(3, 3);
            end
        end

        function M_interp = interpolate_3x3_matrix(drafts, M_cell, draft_target)
            % INTERPOLATE_3X3_MATRIX  Element-wise linear interpolation of
            %   3×3 matrices stored in a cell array.
            %
            %   Vectorised: pre-extracts all elements into an
            %   [n_cases × 9] matrix so interp1 is called ONCE across
            %   all nine columns simultaneously, instead of 9 × n_cases
            %   separate scalar calls.
            %
            %   Element ordering (ii outer, jj inner):
            %     col 1 = M(1,1), col 2 = M(1,2), col 3 = M(1,3),
            %     col 4 = M(2,1), ..., col 9 = M(3,3)
            %   After interp1: reshape to [3×3] row-major, then transpose
            %   to recover the correct [3×3] matrix.
            %
            %   Verified: reshape([a1..a9],3,3) fills column-major in
            %   MATLAB, so transposing recovers the row-major ordering
            %   used when filling elem_mat.

            try
                n_cases = length(drafts);

                % Stack all matrix elements: rows = drafts, cols = (i,j)
                elem_mat = zeros(n_cases, 9);
                col = 1;
                for ii = 1:3
                    for jj = 1:3
                        for k = 1:n_cases
                            elem_mat(k, col) = M_cell{k}(ii, jj);
                        end
                        col = col + 1;
                    end
                end

                % Single vectorised interp1 call across all 9 elements
                row_interp = interp1(drafts, elem_mat, draft_target, 'linear', 'extrap');

                % Reshape and symmetrize
                M_interp = reshape(row_interp, 3, 3)';
                M_interp = 0.5 * (M_interp + M_interp');

            catch
                M_interp = zeros(3, 3);
            end
        end


        %% ═════════════════════════════════════════════════════════════
        %%  RIGID-BODY DYNAMICS (UNCHANGED)
        %% ═════════════════════════════════════════════════════════════

        function M = calculate6x6MassMatrix(mass, cg_position, inertia_tensor)
            % CALCULATE6X6MASSMATRIX  Construct the 6×6 rigid-body mass matrix.
            %
            %   M = [ m·I₃        −m·[r_cg×] ]
            %       [ m·[r_cg×]'   I_cg       ]

            try
                if mass <= 0
                    M = eye(6);
                    return;
                end

                rx = cg_position(1);
                ry = cg_position(2);
                rz = cg_position(3);

                r_skew = [  0,  -rz,   ry;
                           rz,    0,  -rx;
                          -ry,   rx,    0];

                M = zeros(6, 6);
                M(1:3, 1:3) = mass * eye(3);
                M(1:3, 4:6) = -mass * r_skew;
                M(4:6, 1:3) = mass * r_skew';
                M(4:6, 4:6) = inertia_tensor;

            catch ME
                warning('WEC_Core_Functions:MassMatrixFailed', ...
                        'Mass matrix calculation failed: %s', ME.message);
                M = eye(6);
            end
        end


        %% ═════════════════════════════════════════════════════════════
        %%  STRIP INTEGRATION (ADAPTED — parser instead of mesh)
        %% ═════════════════════════════════════════════════════════════

        function [V, CB_z, Iyy_origin, contours, A_samples, z_samples, ...
                 Ixx_origin, Izz_origin] = ...
                compute_strip_bspline(parser, z_lo, z_hi, n_sub)
        % COMPUTE_STRIP_BSPLINE  [DEPRECATED] Use WEC_HydroProperties.compute_strip instead.
        %
        %   This method uses cross-section polygon assembly which produces
        %   incorrect volumes due to incomplete mirror handling.  It is
        %   retained for backward compatibility only.
        %
        %   The replacement WEC_HydroProperties.compute_strip uses the
        %   parametric divergence theorem (same as compute_submerged) with
        %   exact source/mirror dedup — no polygon assembly, no root-finding,
        %   no segment chaining.
        %
        %   See also: WEC_HydroProperties.compute_strip

            warning('WEC_Core_Functions:Deprecated', ...
                    ['compute_strip_bspline is deprecated. ' ...
                     'Use WEC_HydroProperties.compute_strip instead.']);
        %
        %   INPUTS
        %     parser : WEC_MS2_Parser object (parsed geometry)
        %     z_lo   : [m]  strip bottom boundary
        %     z_hi   : [m]  strip top boundary
        %     n_sub  : [-]  number of z-samples (min 5; ≥ 15 recommended)
        %
        %   OUTPUTS
        %     V          : [m³]        strip volume = ∫ A(z) dz
        %     CB_z       : [m]         z-centroid of strip volume
        %     Iyy_origin : [m⁵]        ∫(Iyy_section + z² A) dz  (pitch)
        %     contours   : {n_sub×1}   xy-polygon at each z-sample ([] if failed)
        %     A_samples  : [n_sub×1]   cross-section area at each z
        %     z_samples  : [n_sub×1]   z-coordinates sampled
        %     Ixx_origin : [m⁵]        ∫(Ixx_section + z² A) dz  (roll)
        %     Izz_origin : [m⁵]        ∫(Iyy_section + Ixx_section) dz  (yaw, no z² term)

            if n_sub < 5, n_sub = 5; end
            z_samples    = linspace(z_lo, z_hi, n_sub)';
            A_samples    = NaN(n_sub, 1);
            Iyy_sec_samp = NaN(n_sub, 1);
            Ixx_sec_samp = NaN(n_sub, 1);
            contours     = cell(n_sub, 1);

            for k = 1:n_sub
                wp = WEC_Core_Functions.evaluateCrossSectionMS2( ...
                        parser, z_samples(k));
                if isempty(wp) || isempty(wp{1}) || size(wp{1}, 1) < 3
                    continue;
                end
                poly_xy = wp{1}(:, 1:2);
                contours{k}      = poly_xy;
                A_samples(k)     = polyarea(poly_xy(:,1), poly_xy(:,2));
                [~, iner, ~]     = WEC_Core_Functions.polygeom( ...
                                       poly_xy(:,1), poly_xy(:,2));
                Ixx_sec_samp(k)  = abs(iner(1));   % ∫∫ y² dA  (roll)
                Iyy_sec_samp(k)  = abs(iner(2));   % ∫∫ x² dA  (pitch)
            end

            %% Gap-fill: replace interior NaN by linear interpolation.
            valid = ~isnan(A_samples);
            if sum(valid) < 2
                V = 0;  CB_z = 0.5*(z_lo + z_hi);
                Iyy_origin = 0;  Ixx_origin = 0;  Izz_origin = 0;
                A_samples(isnan(A_samples)) = 0;
                return;
            end
            if any(~valid)
                A_samples(~valid) = max(0, interp1( ...
                    z_samples(valid), A_samples(valid), ...
                    z_samples(~valid), 'linear', 0));
                Iyy_sec_samp(~valid) = max(0, interp1( ...
                    z_samples(valid), Iyy_sec_samp(valid), ...
                    z_samples(~valid), 'linear', 0));
                Ixx_sec_samp(~valid) = max(0, interp1( ...
                    z_samples(valid), Ixx_sec_samp(valid), ...
                    z_samples(~valid), 'linear', 0));
            end

            %% Fit cubic splines and integrate analytically.
            pp_A     = spline(z_samples, A_samples);
            pp_Iyy_s = spline(z_samples, Iyy_sec_samp);
            pp_Ixx_s = spline(z_samples, Ixx_sec_samp);

            V           = WEC_Core_Functions.integrate_pp(pp_A, z_lo, z_hi);
            moment_z    = WEC_Core_Functions.integrate_pp_times_zp( ...
                              pp_A, z_lo, z_hi, 1);
            moment_z2   = WEC_Core_Functions.integrate_pp_times_zp( ...
                              pp_A, z_lo, z_hi, 2);
            Iyy_sec_int = WEC_Core_Functions.integrate_pp(pp_Iyy_s, z_lo, z_hi);
            Ixx_sec_int = WEC_Core_Functions.integrate_pp(pp_Ixx_s, z_lo, z_hi);

            if V > 1e-12
                CB_z = moment_z / V;
            else
                CB_z = 0.5*(z_lo + z_hi);
            end

            %% Moments of inertia about the BODY ORIGIN (unit density).
            %
            %  Iyy_origin = ∫ [∫∫ x² dA + z² A(z)] dz   (pitch, about y-axis)
            %  Ixx_origin = ∫ [∫∫ y² dA + z² A(z)] dz   (roll,  about x-axis)
            %  Izz_origin = ∫ [∫∫ x² dA + ∫∫ y² dA] dz  (yaw,   about z-axis)
            %
            %  Note: Izz has NO z² term because yaw rotation is about the
            %  z-axis.  The perpendicular-axis theorem gives
            %  Izz = Ixx_planar + Iyy_planar at each cross-section.
            Iyy_origin = Iyy_sec_int + moment_z2;
            Ixx_origin = Ixx_sec_int + moment_z2;
            Izz_origin = Iyy_sec_int + Ixx_sec_int;
        end


        %% ─────────────────────────────────────────────────────────────
        %%  PIECEWISE-POLYNOMIAL INTEGRATION HELPERS (UNCHANGED)
        %% ─────────────────────────────────────────────────────────────

        function I = integrate_pp(pp, a, b)
        % INTEGRATE_PP  Definite integral of a piecewise polynomial from a to b.

            I = 0;
            breaks = pp.breaks(:)';
            coefs  = pp.coefs;
            for k = 1:length(breaks)-1
                t_lo = max(a, breaks(k))   - breaks(k);
                t_hi = min(b, breaks(k+1)) - breaks(k);
                if t_hi <= t_lo + 1e-14, continue; end
                c = coefs(k, :);
                I = I + c(1)*(t_hi^4 - t_lo^4)/4 ...
                      + c(2)*(t_hi^3 - t_lo^3)/3 ...
                      + c(3)*(t_hi^2 - t_lo^2)/2 ...
                      + c(4)*(t_hi   - t_lo);
            end
        end


        function I = integrate_pp_times_zp(pp, a, b, p)
        % INTEGRATE_PP_TIMES_ZP  Definite integral of z^p * f(z) dz from a to b,
        %   where f is a piecewise polynomial.

            I = 0;
            breaks = pp.breaks(:)';
            coefs  = pp.coefs;
            for k = 1:length(breaks)-1
                z_k  = breaks(k);
                t_lo = max(a, breaks(k))   - z_k;
                t_hi = min(b, breaks(k+1)) - z_k;
                if t_hi <= t_lo + 1e-14, continue; end
                c = coefs(k, :);
                if p == 1
                    I = I + c(1)*(t_hi^5 - t_lo^5)/5 ...
                          + (c(2) + z_k*c(1))*(t_hi^4 - t_lo^4)/4 ...
                          + (c(3) + z_k*c(2))*(t_hi^3 - t_lo^3)/3 ...
                          + (c(4) + z_k*c(3))*(t_hi^2 - t_lo^2)/2 ...
                          + z_k*c(4)*(t_hi - t_lo);
                elseif p == 2
                    zk2 = z_k^2;
                    I = I + c(1)*(t_hi^6 - t_lo^6)/6 ...
                          + (c(2) + 2*z_k*c(1))*(t_hi^5 - t_lo^5)/5 ...
                          + (c(3) + 2*z_k*c(2) + zk2*c(1))*(t_hi^4 - t_lo^4)/4 ...
                          + (c(4) + 2*z_k*c(3) + zk2*c(2))*(t_hi^3 - t_lo^3)/3 ...
                          + (2*z_k*c(4) + zk2*c(3))*(t_hi^2 - t_lo^2)/2 ...
                          + zk2*c(4)*(t_hi - t_lo);
                else
                    error('WEC_Core_Functions:integrate_pp_times_zp', ...
                          'Only p = 1 or p = 2 are supported; got p = %d', p);
                end
            end
        end

    end
end