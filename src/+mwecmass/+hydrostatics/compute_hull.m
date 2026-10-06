function props = compute_hull(parser, options)
%COMPUTE_HULL Integrate a closed parametric hull and return volume properties.
% Uses divergence-theorem quadrature on source surfaces, derives reflected
% mirror contributions, and closes detected open boundaries with flat caps.
% options.n_quad (default 20), n_edge_samples (50), edge_tol [m] (1e-4),
% and verbose (true) control quadrature, edge matching, and progress output.

    if nargin < 2, options = struct(); end
    if ~isfield(options, 'n_quad'),         options.n_quad         = 20;   end
    if ~isfield(options, 'n_edge_samples'), options.n_edge_samples = 50;   end
    if ~isfield(options, 'edge_tol'),       options.edge_tol       = 1e-4; end
    if ~isfield(options, 'verbose'),        options.verbose        = true; end

    surf_names = parser.visible_surfs;
    n_surfs    = length(surf_names);
    n_quad     = options.n_quad;

    % Estimate a hull interior point for orientation checks.
    interior = (parser.extents(1:3) + parser.extents(4:6))' / 2;

    % Classify surfaces for mirror deduplication.
    topo = parser.classify_visible_surfaces();

    % Integrate source surfaces and derive mirrored contributions.
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
            mwecmass.hydrostatics.surface_integral( ...
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
                mwecmass.hydrostatics.surface_integral( ...
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

    % Detect open edges.
    [open_edges, ~] = mwecmass.geometry.detect_open_edges( ...
        parser, options.n_edge_samples, options.edge_tol);

    if options.verbose
        fprintf('    Open edges: %d\n', length(open_edges));
        for i = 1:length(open_edges)
            oe = open_edges{i};
            fprintf('      %s edge %d  (z ≈ %.3f m)\n', ...
                    oe.surface, oe.edge_idx, oe.z_mean);
        end
    end

    % Add cap contributions.
    V_cap_total  = 0;
    Cxyz_cap_total = [0 0 0];

    if ~isempty(open_edges)
        [V_cap, Cxyz_cap, A_cap, z_cap] = ...
            mwecmass.hydrostatics.cap_contribution( ...
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

    % Assemble output.
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

    % Report.
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
