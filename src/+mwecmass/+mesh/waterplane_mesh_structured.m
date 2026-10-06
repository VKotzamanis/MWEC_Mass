function [wp_nodes, wp_panels, wp_nverts] = waterplane_mesh_structured(contour_pts, target_edge, z_wl)
%WATERPLANE_MESH_STRUCTURED Create a structured waterplane lid from an ordered contour.
% Syntax: [wp_nodes,wp_panels,wp_nverts] = waterplane_mesh_structured(contour_pts,target_edge,z_wl).
% Inputs: contour_pts [N x 2] ordered boundary in m; target_edge [m] (default 0.4); z_wl [m].
% Outputs: nodes [Nv x 3] m, panels [Np x 4] (triangles use repeated v3/v4), and counts.
% A transfinite grid follows the two boundary curves; panel winding is set for downward normals.

    if nargin < 2 || isempty(target_edge); target_edge = 0.4; end
    if nargin < 3 || isempty(z_wl);        z_wl = 0; end

    xy = contour_pts(:, 1:2);
    N  = size(xy, 1);
    if N < 4
        wp_nodes = zeros(0,3); wp_panels = zeros(0,4);
        wp_nverts = zeros(0,1); return;
    end

    D_max = 0; iA = 1; iB = 1;
    for i = 1:N; for j = i+1:N %#ok<ALIGN>
        d = norm(xy(i,:) - xy(j,:));
        if d > D_max; D_max = d; iA = i; iB = j; end
    end; end
    if iA > iB; tmp = iA; iA = iB; iB = tmp; end

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

    sL = wp_arc_len(cL);
    sR = wp_arc_len(cR);

    half_perim = max(sL(end), sR(end));
    max_width = 0;
    for k = 1:50
        t_probe = (k-1)/49;
        pL = wp_interp(cL, sL, t_probe);
        pR = wp_interp(cR, sR, t_probe);
        max_width = max(max_width, norm(pR - pL));
    end

    Nt = max(6, ceil(half_perim / target_edge) + 1);
    Ns = max(4, ceil(max_width  / target_edge) + 1);
    if mod(Ns, 2) == 0; Ns = Ns + 1; end  % centre column

    t_grid = linspace(0, 1, Nt);
    s_grid = linspace(0, 1, Ns);

    gxy = zeros(Nt, Ns, 2);
    for k = 1:Nt
        pL = wp_interp(cL, sL, t_grid(k));
        pR = wp_interp(cR, sR, t_grid(k));
        for j = 1:Ns
            gxy(k, j, :) = (1 - s_grid(j)) * pL + s_grid(j) * pR;
        end
    end

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

function s = wp_arc_len(curve)
%WP_ARC_LEN  Cumulative arc length along a 2D polyline.
% Inputs:  curve  double [Nx2], ordered polyline vertices
% Outputs: s      double [Nx1], cumulative arc length, s(1) = 0
    s = zeros(size(curve, 1), 1);
    for k = 2:size(curve, 1)
        s(k) = s(k-1) + norm(curve(k,:) - curve(k-1,:));
    end
end

function pt = wp_interp(curve, s_cum, t)
%WP_INTERP  Interpolates a position on an arc-length-parameterised curve at fraction t in [0,1].
% Inputs:  curve  double [Nx2], ordered polyline vertices
%          s_cum  double [Nx1], cumulative arc length along curve (from wp_arc_len)
%          t      double scalar in [0,1], fractional arc-length position
% Outputs: pt     double [1x2], interpolated (x,y) position at fractional arc length t
    st = t * s_cum(end);
    if st <= 0;        pt = curve(1,:);   return; end
    if st >= s_cum(end); pt = curve(end,:); return; end
    idx = find(s_cum >= st, 1, 'first');
    if idx <= 1;       pt = curve(1,:);   return; end
    f = (st - s_cum(idx-1)) / (s_cum(idx) - s_cum(idx-1));
    pt = curve(idx-1,:) + f * (curve(idx,:) - curve(idx-1,:));
end
