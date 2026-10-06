function strip = compute_strip(parser, z_lo, z_hi, options)
%COMPUTE_STRIP Compute volume and inertia of the horizontal band [z_lo,z_hi].
% With Aw_table_z/Aw_table, integrates cross-section area (and optional
% waterplane moments) directly; otherwise falls back to surface integration.
% Output strip fields V [m^3], CB_z [m], Ixx/Iyy/Izz [m^5], and raw second
% moments int_x2/int_y2/int_z2. options.n_quad defaults to 16.

    if nargin < 4, options = struct(); end
    if ~isfield(options, 'n_quad'),        options.n_quad        = 16; end
    if ~isfield(options, 'Aw_table_z'),    options.Aw_table_z    = []; end
    if ~isfield(options, 'Aw_table'),      options.Aw_table      = []; end
    if ~isfield(options, 'I_wp_xx_table'), options.I_wp_xx_table = []; end
    if ~isfield(options, 'I_wp_yy_table'), options.I_wp_yy_table = []; end
    if ~isfield(options, 'boundary_cache'), options.boundary_cache = []; end

    %% Integrate the cross-section table over [z_lo, z_hi].
    %
    %  Used whenever the Aw_table is available (always in the config-
    %  builder call).  Replaces the divergence-theorem surface
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

    %% Surface-integral fallback when no Aw_table is available.

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
            mwecmass.hydrostatics.surface_integral_strip( ...
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
            b_cache = mwecmass.geometry.precompute_boundary_cache( ...
                          parser, n_quad * 5);
        end

        wl_hi = mwecmass.geometry.extract_isocurve_at_z( ...
                    parser, z_hi, n_quad * 5, b_cache);
        wl_lo = mwecmass.geometry.extract_isocurve_at_z( ...
                    parser, z_lo, n_quad * 5, b_cache);

        if ~isempty(wl_hi) && size(wl_hi, 1) >= 3
            Aw_hi = mwecmass.hydrostatics.waterplane_properties(wl_hi);
        else
            Aw_hi = 0;
        end
        if ~isempty(wl_lo) && size(wl_lo, 1) >= 3
            Aw_lo = mwecmass.hydrostatics.waterplane_properties(wl_lo);
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

    % Transform z moments from shifted to world coordinates.
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
