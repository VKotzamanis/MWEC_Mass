function boundary_xy = hull_waterline_polygon(mesh, z_tol)
%HULL_WATERLINE_POLYGON Extract the largest ordered open-edge loop on the hull waterline.
% Syntax: boundary_xy = hull_waterline_polygon(mesh,z_tol).
% Input: mesh from mesh.generate; optional z_tol [m] (default 1e-6). Output boundary_xy [N x 2] m,
% without a repeated endpoint; returns [0 x 2] with a warning if no loop has at least three nodes.
% Edges must be connected after trimming/merging; largest area is selected as the outer loop.
% See docs/METHODS_ENGINE.md#waterline-boundary

    if nargin < 2 || isempty(z_tol), z_tol = 1e-6; end

    verts = mesh.vertices;
    n_p   = mesh.n_panels;

    % Waterline vertices are the union of
    %   (a) vertices flagged by trim_at_wl/split_panel_at_z, AND
    %   (b) original mesh vertices that already sit at z = z_wl
    %       (within machine precision)
    % Both sets are needed because a parametric grid can land
    % a row at exactly z_wl, no trim runs for that surface, so its
    % z=0 vertices are NOT in waterline_verts — but they ARE on
    % the true waterline boundary. A second surface whose grid
    % straddles z_wl gets trimmed and its intersections ARE flagged.
    % Using only one source misses half the perimeter (Nu=23 / Nu=45
    % on C0: 22% / 12% area deficit). A 1e-6 tolerance excludes
    % that only true z=0 nodes pass.
    wl_set = abs(verts(:,3)) < z_tol;
    if isfield(mesh,'waterline_verts') && ~isempty(mesh.waterline_verts)
        wl_set(mesh.waterline_verts) = true;
    end
    on_wl = @(vi) wl_set(vi);

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

    %  An edge qualifies when it is open (count==1) and both endpoints
    %  are members of the waterline vertex set (above).  This captures:
    %    - The waterline arc (hull trim boundary).
    %    - For half-body: the y=0 symmetry-plane closure at z=0.
    %  It excludes interior edges, submerged boundary edges, and
    %  near-z=0-but-not-on-trim edges that a broad tolerance would let
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
        warning('mwecmass:mesh:HullWLEmpty', ...
            'No hull waterline open edges at |z|<%.3f m. Check trim_wl=true and WP-2 merge.', z_tol);
        boundary_xy = zeros(0, 2);
        return;
    end

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

    if isempty(loops)
        warning('mwecmass:mesh:HullWLNoLoop', ...
            'Edge walk produced no loop >= 3 nodes. Verify the post-trim merge ran in the panel mesher.');
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
