function mesh = generate(parser, draft, Nu, Nv, options)
%GENERATE Build a panel mesh from parsed MS2 surfaces at a draft and grid density.
% Syntax: mesh = mwecmass.mesh.generate(parser,draft,Nu,Nv,options), where draft [m]
% is vertical shift and Nu,Nv are grid points per parametric direction. options
% controls trimming, mirroring, spacing, merging, and verbosity; output contains
% vertices [N×3 m], panels [P×4], normals, IDs, caps, and waterline data.

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
        error('mwecmass:mesh:MutualExclusion', ...
              'half_body and quarter_body are mutually exclusive.');
    end
    %  vertex_merge options:
    %    'merge'  — tolerance-based merge (1e-8 m), then unique
    %    'dedup'  — exact coordinate match only (unique rows)
    %    'none'   — skip entirely (duplicate vertices remain)

    topo = parser.classify_visible_surfaces();

    all_verts    = [];
    all_panels   = [];
    all_surf_ids = [];
    all_is_cap   = [];

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

    % Waterline-conforming grid: insert the shared crossing into u_grid
    % before arc-length spacing so trimming does not split panels.

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

            % Find zero-crossings and endpoint contacts, then bisect to
            % machine precision so the shared waterline row is exact.
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

    % Arc-length spacing uses the first source profile and preserves the
    % waterline as a segment boundary.

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
            u_grid = enforce_min_spacing(u_grid, min_du);

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

    % Per-surface v-grid.
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
    % Measure each surface's max v-extent and compute Nv
    %  so that panel v-width ≈ max physical u-step (making
    %  panels approximately square at the widest point).
    %
    %  Rationale for independent Nv values:
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

    % Quarter-body setup before source evaluation.
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
            mwecmass.mesh.detect_source_quadrant(parser, topo.sources);
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
        % ISX=1 mirrors X → −X,
        % but the mesh covers both x>0 (Q1) and x<0 (Q2).  HAMS
        % would create Green's function images that overlap with
        % existing panels, double-counting the hull.
        [negate_x, negate_y, source_quadrant] = ...
            mwecmass.mesh.detect_source_quadrant(parser, topo.sources);
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

        % Store the untransformed grid for mirror construction (full body)
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

        % Half-body: store the transformed grid for mirror construction.
        %   Unlike full-body (which stores untransformed), the
        %   half-body mirrors need the Q1-transformed grid as their
        %   starting point — the X-flip produces Q2, completing
        %   the y>=0 half.
        if options.half_body
            source_grids(sname) = S;
        end

        [verts, panels] = grid_to_quads(S);
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

    % Build mirror meshes by coordinate flip.
    if options.quarter_body
        if options.verbose
            fprintf('    Quarter-body mode: skipping %d mirror surfaces\n', ...
                    length(topo.mirrors));
        end
    else
    for m = 1:length(topo.mirrors)
        mirr = topo.mirrors(m);

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
            % This applies regardless of negate_x.  The quadrant
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

        [verts, panels] = grid_to_quads(S_mirr);

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

    % Vertex handling at surface seams.
    switch options.vertex_merge
        case 'merge'
            % Tolerance-based: round to 1e-8, then unique
            [all_verts, all_panels] = merge_vertices( ...
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

    % Trim at waterline.
    wl_verts = [];
    if options.trim_wl
        [all_verts, all_panels, all_surf_ids, all_is_cap, wl_verts] = ...
            mwecmass.mesh.trim_at_waterline(all_verts, all_panels, ...
                all_surf_ids, all_is_cap, options.z_wl);
        if options.verbose
            fprintf('    Trimmed at z=%.3f: %d panels remain\n', ...
                    options.z_wl, size(all_panels, 1));
        end
    end

    % Merge again after trimming: independently generated intersections
    % must share indices for a connected waterline ring.
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

    % Close open edges only for visualization or volume validation.
    %  For HAMS input, leave close_gaps = false (the default).
    %  Enable only for visualization or volume-validation meshes
    %  where a closed surface is needed.
    if options.close_gaps
        [all_verts, all_panels, all_surf_ids, all_is_cap, n_cap] = ...
            mwecmass.mesh.close_open_edges(all_verts, all_panels, ...
                all_surf_ids, all_is_cap);
        if options.verbose && n_cap > 0
            fprintf('    Closed %d open edges with cap panels\n', n_cap);
        end
    end

    % Compute panel normals from the cross product.
    %  The surface parameterisation gives S_u × S_v outward.
    %  grid_to_quads preserves this: diagonal cross product
    %  ∝ 2(S_u × S_v).  Mirror winding reversal corrects odd
    %  reflections.  No orient_normals heuristic needed.
    normals = mwecmass.mesh.compute_panel_normals( ...
        all_verts, all_panels);

    % Remove orphaned vertices.
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

    % Verify Q1 coordinates for quarter-body meshes.
    %  The quadrant transform was applied at the grid level
    %  (before grid_to_quads).  Verify all vertices are in Q1.
    if options.quarter_body && options.verbose
        x_min = min(all_verts(:,1));
        y_min = min(all_verts(:,2));
        if x_min < -1e-8 || y_min < -1e-8
            warning('mwecmass:mesh:QuadrantViolation', ...
                    'Vertices outside Q1: x_min=%.6f, y_min=%.6f', ...
                    x_min, y_min);
        else
            fprintf('    Q1 verified: x_min=%.2e, y_min=%.2e\n', ...
                    x_min, y_min);
        end
    end

    % Verify half-body coordinates.
    %  All vertices must have y >= 0 (HAMS mirrors across y=0).
    %  x may be positive or negative (both Q1 and Q2 are present).
    if options.half_body && options.verbose
        x_min = min(all_verts(:,1));
        x_max = max(all_verts(:,1));
        y_min = min(all_verts(:,2));
        if y_min < -1e-8
            warning('mwecmass:mesh:HalfBodyYViolation', ...
                    'Half-body mesh has y < 0: y_min=%.6f', y_min);
        else
            fprintf('    Half-body verified: y >= 0 (y_min=%.2e)\n', y_min);
            fprintf('    x range: [%.4f, %.4f] (both sides of X=0)\n', ...
                    x_min, x_max);
        end
    end

    % Assemble output mesh struct.
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

function u = enforce_min_spacing(u, min_du)
%ENFORCE_MIN_SPACING Remove interior points closer than min_du to the last kept
% point; endpoints are always retained. Inputs/outputs are vectors in [0,1],
% and min_du is the minimum parametric spacing.

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

function [verts, panels] = grid_to_quads(S)
%GRID_TO_QUADS Convert an [Nu×Nv×3] surface grid to quad connectivity.
% reshape(S,[],3) is column-major, so vertex i+(j-1)*Nu corresponds to S(i,j,:);
% each cell is ordered (i,j),(i+1,j),(i+1,j+1),(i,j+1).

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

function [verts, panels] = merge_vertices(verts, panels, tol)
%MERGE_VERTICES Merge coincident vertices within tol [m] (default 1e-8) and
% reindex panels to close evaluation seams without merging distinct geometry.

    if nargin < 3, tol = 1e-8; end

    rounded = round(verts / tol) * tol;
    [~, ia, ic] = unique(rounded, 'rows', 'stable');

    verts  = verts(ia, :);
    panels = ic(panels);
end
