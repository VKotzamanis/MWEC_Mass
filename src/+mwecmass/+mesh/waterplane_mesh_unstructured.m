function [wp_nodes, wp_panels, wp_nverts] = waterplane_mesh_unstructured(boundary_xy, target_edge, x_sym, y_sym, z_wl, lock_boundary)
%WATERPLANE_MESH_UNSTRUCTURED Mesh a waterplane polygon with constrained Delaunay points and quad matching.
% Syntax: [wp_nodes,wp_panels,wp_nverts] = waterplane_mesh_unstructured(boundary_xy,target_edge,x_sym,y_sym,z_wl,lock_boundary).
% Inputs: boundary_xy [N x 2] m; target_edge [m]; symmetry flags; z_wl [m]; optional lock_boundary.
% Outputs: nodes [Nv x 3] m, panels [Np x 4], and per-panel vertex counts. Panels are wound downward.
% lock_boundary=true preserves supplied hull waterline nodes while still refining the interior.

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
        xy = resample_polygon(xy, target_edge);
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
        % Structured fallback cannot preserve lock_boundary, breaking run_at_draft's WP-hull match.
        if lock_boundary
            error('mwecmass:mesh:CDTFailedLockedBoundary', ...
                ['CDT failed (%s) with lock_boundary=true. The structured-mesh fallback cannot ' ...
                 'honour a locked boundary, so it is not used here.'], ME_dt.message);
        end
        warning('mwecmass:mesh:CDTFailed','CDT failed: %s. TFI fallback.',ME_dt.message);
        [wp_nodes,wp_panels,wp_nverts] = ...
            mwecmass.mesh.waterplane_mesh_structured(boundary_xy,target_edge,z_wl);
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
        cq(c)=quad_scaled_jacobian( ...
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
%QUAD_SCALED_JACOBIAN  Minimum scaled-Jacobian quality metric over the 4 corners of a quad: 0 = degenerate, 1 = perfect square.
% Inputs:  p1,p2,p3,p4  double [1x2] each, the quad's 4 corner (x,y) coordinates in order
% Outputs: Q            double scalar in [0,1], the worst-corner scaled-Jacobian quality metric
    pts=[p1(1:2);p2(1:2);p3(1:2);p4(1:2)]; Q=1.0;
    for k=1:4
        pk=pts(k,:); pp=pts(mod(k-2,4)+1,:); pn=pts(mod(k,4)+1,:);
        e1=pn-pk; e2=pp-pk; l1=norm(e1); l2=norm(e2);
        if l1<1e-14||l2<1e-14, Q=0; return; end
        Q=min(Q, abs(e1(1)*e2(2)-e1(2)*e2(1))/(l1*l2));
    end
end

function xy_out = resample_polygon(xy, max_edge)
%RESAMPLE_POLYGON  Inserts midpoints into a closed polygon so that no edge exceeds max_edge.
% Inputs:  xy       double [Nx2], ordered closed-polygon vertices
%          max_edge double scalar, maximum allowed edge length before subdivision
% Outputs: xy_out   double [Mx2], resampled polygon vertices (M >= N)
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
