function props = compute_submerged(parser, z_wl, options)
%COMPUTE_SUBMERGED Compute volume and hydrostatics below z_wl [m].
% Surface integration limits are clipped to z <= z_wl; the waterplane and
% detected open edges close the divergence-theorem domain. options controls
% quadrature, edge matching, bisection, verbosity, and optional waterplane
% moment overrides. Outputs include V_sub [m^3], CB [1x3 m], Aw [m^2],
% I_wp_xx/I_wp_yy [m^4], per-surface volumes, and waterline boundary points.

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
    %  Keep the test point at the hull bounding-box center so it remains
    %  an interior point even when the waterline is below the hull center.
    interior = (parser.extents(1:3) + parser.extents(4:6))' / 2;

    % Origin shift for numerical conditioning.
    %  z_ref = midpoint between hull bottom and waterline.
    %  Shifts the position vector S → S' = (x, y, z−z_ref)
    %  in the divergence theorem.  Cap terms become proportional
    %  to (z_wl − z_ref) instead of z_wl, eliminating the
    %  catastrophic cancellation that made V_sub non-monotonic.
    z_ref = 0.5 * (min(parser.extents(3), parser.extents(6)) + z_wl);

    % Classify surfaces for mirror deduplication.
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

    % Integrate source surfaces and derive mirrored contributions.
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
    [u_o, wu_o] = mwecmass.internal.gauss_legendre(n_orient);
    [v_o, wv_o] = mwecmass.internal.gauss_legendre(n_orient);
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
            mwecmass.hydrostatics.surface_integral_submerged( ...
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
            % Use the same coarse signed-volume orientation test as the
            % source path; a single midpoint normal is unreliable on
            % non-convex parameterized surfaces.
            n_fb = 5;
            [u_fb, wu_fb] = mwecmass.internal.gauss_legendre(n_fb);
            [v_fb, wv_fb] = mwecmass.internal.gauss_legendre(n_fb);
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
                mwecmass.hydrostatics.surface_integral_submerged( ...
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

    % Waterplane closure.
    % Use precomputed waterplane moments when supplied; otherwise trace the spline contour.

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
        % Trace the iso-z contour, using a boundary cache when supplied or building one as needed.
        if isfield(options, 'boundary_cache') && ~isempty(options.boundary_cache)
            b_cache = options.boundary_cache;
        else
            b_cache = mwecmass.geometry.precompute_boundary_cache(parser, n_quad * 5);
        end

        wl_pts_all = mwecmass.geometry.extract_isocurve_at_z( ...
            parser, z_wl, n_quad * 5, b_cache);

        if ~isempty(wl_pts_all) && size(wl_pts_all, 1) >= 3
            [Aw, I_xx, I_yy, wl_ordered] = ...
                mwecmass.hydrostatics.waterplane_properties(wl_pts_all);
        else
            Aw = 0; I_xx = 0; I_yy = 0;
        end

        if options.verbose
            fprintf('    Waterplane (fast iso-z, %d pts): Aw = %.6f m²\n', ...
                    size(wl_pts_all, 1), Aw);
        end
    end

    % Cap at z_wl with the shifted position vector.
    %  S' = (x, y, z_wl − z_ref) on the cap.  n = (0, 0, +1).
    %  S'·n = z_wl − z_ref.
    dz_cap = z_wl - z_ref;
    V_cap    = dz_cap * Aw / 3;
    Cxyz_cap = [0, 0, dz_cap^2 * Aw / 2];
    V_total  = V_total + V_cap;
    Cxyz_num = Cxyz_num + Cxyz_cap;

    % Shifted second moment: (z_wl − z_ref)³ × Aw / 3
    int_z2_total = int_z2_total + dz_cap^3 * Aw / 3;

    % Close detected open edges below the waterline.
    if ~options.skip_open_edges
        [open_edges, ~] = mwecmass.geometry.detect_open_edges( ...
            parser, options.n_edge_samples, 1e-4);

        below_wl_edges = {};
        for i = 1:length(open_edges)
            if open_edges{i}.z_mean < z_wl - 0.05
                below_wl_edges{end+1} = open_edges{i}; %#ok<AGROW>
            end
        end

        if ~isempty(below_wl_edges)
            [V_oe, ~, A_oe, z_oe] = ...
                mwecmass.hydrostatics.cap_contribution( ...
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

    % Transform z moments from shifted to world coordinates.
    %  V_total is origin-invariant (correct as-is).
    %  Cxyz_num(1,2) are unshifted (x,y not shifted).
    %  Cxyz_num(3) and int_z2 need the parallel-axis unshift:
    %
    %  ∫z dV = ∫(z−z_ref) dV + z_ref × V
    %  ∫z² dV = ∫(z−z_ref)² dV + 2·z_ref·∫(z−z_ref)dV + z_ref²·V

    Cz_shifted   = Cxyz_num(3);                  % ∫(z−z_ref) dV
    Cxyz_num(3)  = Cz_shifted + z_ref * V_total; % ∫z dV
    int_z2_total = int_z2_total + 2*z_ref*Cz_shifted + z_ref^2 * V_total;

    % Assemble output.
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
