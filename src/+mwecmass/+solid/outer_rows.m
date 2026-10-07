function [rows, grids] = outer_rows(model, cache, z, grids)
%OUTER_ROWS Closed horizontal sections of the hull surface at the heights z.
%   rows = OUTER_ROWS(model, cache, z) with model = MS2Parser.parse output,
%   cache = precompute_boundary_cache(model, n) and z [m, body frame].
%   [rows, grids] = OUTER_ROWS(model, cache, z, grids) reuses the grids from an
%   earlier call (OUTER_ROWS(model, cache, []) builds them without sections).
%
%   Every source patch is traced on the parameter grid cache.u_samples (n x n)
%   by marching squares on z(u,v) - z_k. Each crossing of a grid edge is
%   refined to rounding on the exact surface, so every point is exact and its
%   (u,v) is known. Mirror patches are the source pieces with the coordinate
%   flip of the deck. Pieces are joined end to end by mutually nearest
%   endpoints; the joined endpoints are one point stored once.
%   The section is not assumed to be of any shape; only a section made of one
%   closed loop is supported (error mwecmass:solid:SectionNotClosed otherwise).
%   RuledSurf and RevSurf evaluate through their cached u-samples; any other
%   surface type runs through MS2Parser.eval_surface (exact, slower).
%
%   rows(k) fields:
%     z, n          height and number of points
%     pts           [n x 3] counter-clockwise, implicitly closed (no repeated first point)
%     patch, u, v   [n x 1] source: index into model.visible_surfs and parameters
%     seam          [n x 1] logical, point lies on a joint between two patches
%     patch2,u2,v2  [n x 1] the same point as seen from the other patch (0, NaN at non-seams)
%     area          enclosed area [m^2] (> 0)
%     seam_gap      largest distance between two endpoints identified at a joint [m]
%     degenerate    true if z is not inside the open height range of the hull
%                   (then n = 0 and pts is empty)

    if nargin < 4 || isempty(grids)
        grids = build_grids(model, cache);
    end
    z = z(:);
    rows = repmat(empty_row(0), 0, 1);
    for k = 1:numel(z)
        rows(k, 1) = section_at(model, cache, grids, z(k));
    end
end

function row = empty_row(zk)
    row = struct('z', zk, 'n', 0, 'pts', zeros(0, 3), 'patch', zeros(0, 1), ...
                 'u', zeros(0, 1), 'v', zeros(0, 1), 'seam', false(0, 1), ...
                 'patch2', zeros(0, 1), 'u2', zeros(0, 1), 'v2', zeros(0, 1), ...
                 'area', 0, 'seam_gap', 0, 'degenerate', true);
end

% ---------------------------------------------------------------- grids

function grids = build_grids(model, cache)
    U = cache.u_samples(:);
    n = numel(U);
    grids.U = U;
    grids.names = cache.sources(:)';
    grids.patch = struct('name', {}, 'kind', {}, 'Z', {}, 'c1', {}, 'c2', {}, ...
                         'pp', {}, 'axis_start', {}, 'axis_dir', {}, ...
                         'phi0', {}, 'phi1', {}, 'profile', {}, 'curve1', {}, ...
                         'curve2', {});
    zr = [Inf, -Inf];
    for s = 1:numel(cache.sources)
        sname = cache.sources{s};
        d = cache.data(sname);
        g = struct('name', sname, 'kind', 'generic', 'Z', [], 'c1', [], 'c2', [], ...
                   'pp', [], 'axis_start', [], 'axis_dir', [], 'phi0', 0, 'phi1', 0, ...
                   'profile', '', 'curve1', '', 'curve2', '');
        e = model.entities(sname);
        switch d.type
            case 'RuledSurf'
                g.kind = 'ruled';
                g.c1 = d.pts_1;
                g.c2 = d.pts_2;
                g.curve1 = e.params.curve1;
                g.curve2 = e.params.curve2;
            case 'RevSurf'
                g.kind = 'rev';
                g.pp = d.profile_pts;
                g.axis_start = d.axis_start;
                g.axis_dir = d.axis_dir;
                g.phi0 = d.angle_start * pi / 180;
                g.phi1 = d.angle_end * pi / 180;
                g.profile = e.params.profile;
        end
        Z = zeros(n, n);
        for i = 1:n
            S = patch_F(g, grid_P(g, i, U), U, model);
            Z(i, :) = S(:, 3)';
        end
        g.Z = Z;
        zr = [min(zr(1), min(Z(:))), max(zr(2), max(Z(:)))];
        grids.patch(s) = g;
    end
    grids.z_range = zr;
end

function P = grid_P(g, i, U)
    switch g.kind
        case 'ruled'
            P = [g.c1(i, :); g.c2(i, :)];
        case 'rev'
            P = g.pp(i, :);
        otherwise
            P = U(i);
    end
end

function P = patch_P(model, g, u)
% Parameter-dependent part of the surface point (the expensive curve evaluations).
    switch g.kind
        case 'ruled'
            P = [model.eval_curve(g.curve1, u); model.eval_curve(g.curve2, u)];
        case 'rev'
            P = model.eval_curve(g.profile, u);
        otherwise
            P = u;
    end
end

function S = patch_F(g, P, v, model)
% Surface points [numel(v) x 3] for the parameter-dependent part P.
    v = v(:);
    switch g.kind
        case 'ruled'
            S = (1 - v) * P(1, :) + v * P(2, :);
        case 'rev'
            rel = P - g.axis_start;
            proj = g.axis_start + dot(rel, g.axis_dir) * g.axis_dir;
            radial = P - proj;
            r = norm(radial);
            if r < 1e-12
                S = repmat(P, numel(v), 1);
                return;
            end
            e_r = radial / r;
            a = g.axis_dir;
            e_t = [a(2) * e_r(3) - a(3) * e_r(2), a(3) * e_r(1) - a(1) * e_r(3), ...
                   a(1) * e_r(2) - a(2) * e_r(1)];
            e_t = e_t / norm(e_t);
            phi = g.phi0 + v * (g.phi1 - g.phi0);
            S = proj + r * cos(phi) * e_r + r * sin(phi) * e_t;
        otherwise
            S = zeros(numel(v), 3);
            for k = 1:numel(v)
                S(k, :) = model.eval_surface(g.name, P, v(k));
            end
    end
end

% ------------------------------------------------------------- sections

function row = section_at(model, cache, grids, zk)
    row = empty_row(zk);
    if zk <= grids.z_range(1) || zk >= grids.z_range(2)
        return;
    end
    pieces = cell(1, 0);
    src_pieces = cell(numel(grids.patch), 1);
    for s = 1:numel(grids.patch)
        src_pieces{s} = march_patch(model, grids.patch(s), grids.U, zk);
    end
    names = model.visible_surfs(:)';
    for s = 1:numel(grids.patch)
        id = find(strcmp(names, grids.names{s}), 1);
        for q = 1:numel(src_pieces{s})
            pc = src_pieces{s}{q};
            pc.patch = id;
            pieces{end + 1} = pc; %#ok<AGROW>
        end
    end
    for m = 1:numel(cache.mirrors)
        s = find(strcmp(grids.names, cache.mirrors(m).ultimate_source), 1);
        id = find(strcmp(names, cache.mirrors(m).name), 1);
        for q = 1:numel(src_pieces{s})
            pc = src_pieces{s}{q};
            pc.patch = id;
            for f = 1:numel(cache.mirrors(m).effective_flips)
                pc.pts(:, mirror_column(cache.mirrors(m).effective_flips{f})) = ...
                    -pc.pts(:, mirror_column(cache.mirrors(m).effective_flips{f}));
            end
            pieces{end + 1} = pc; %#ok<AGROW>
        end
    end
    row = join_pieces(pieces, zk);
end

function c = mirror_column(plane)
    switch plane
        case 'X'
            c = 1;
        case 'Y'
            c = 2;
        otherwise
            error('mwecmass:solid:UnsupportedMirror', ...
                  'Mirror plane %s cannot map a horizontal section onto itself.', plane);
    end
end

function pieces = march_patch(model, g, U, zk)
% Iso-z curve of one patch as ordered pieces with parameters.
    n = numel(U);
    f = g.Z - zk;
    pos = f >= 0;
    % z of a patch point carries about 5 ulp of rounding (projected edge curves of C1); 16 ulp stops the root search there
    tolz = 16 * eps(max(1, abs(zk)));

    cu = pos(1:n-1, :) ~= pos(2:n, :);
    cv = pos(:, 1:n-1) ~= pos(:, 2:n);
    [iu, ju] = find(cu);
    [iv, jv] = find(cv);
    ord = sortrows([iu, ju]);
    iu = ord(:, 1); ju = ord(:, 2);
    nu = numel(iu);
    nv = numel(iv);
    cnt = nu + nv;
    pieces = cell(1, 0);
    if cnt == 0
        return;
    end

    xu = zeros(cnt, 1);
    xv = zeros(cnt, 1);
    gi = zeros(cnt, 1);
    Pu = cell(cnt, 1);
    id_u = zeros(n - 1, n);
    id_v = zeros(n, n - 1);

    x_prev = NaN;
    i_prev = 0;
    P_prev = [];
    for k = 1:nu
        i = iu(k); j = ju(k);
        id_u(i, j) = k;
        if i ~= i_prev
            x_prev = NaN;
        end
        x = NaN;
        if ~isnan(x_prev)
            if abs(patch_F_z(g, P_prev, U(j), model) - zk) <= tolz
                x = x_prev;
            end
        end
        if isnan(x)
            fun = @(x) edge_z_u(model, g, x, U(j)) - zk;
            x = solve_edge(fun, U(i), U(i + 1), f(i, j), f(i + 1, j), tolz, x_prev);
            P_prev = patch_P(model, g, x);
        end
        Pu{k} = P_prev;
        x_prev = x;
        i_prev = i;
        xu(k) = x;
        xv(k) = U(j);
    end
    for k = 1:nv
        i = iv(k); j = jv(k);
        id_v(i, j) = nu + k;
        P = grid_P(g, i, U);
        fun = @(x) patch_F_z(g, P, x, model) - zk;
        x = solve_edge(fun, U(j), U(j + 1), f(i, j), f(i, j + 1), tolz, NaN);
        xu(nu + k) = U(i);
        xv(nu + k) = x;
        gi(nu + k) = i;
    end

    seg = zeros(0, 2);
    [ci, cj] = find((cu(:, 1:n-1) + cu(:, 2:n) + cv(1:n-1, :) + cv(2:n, :)) >= 2);
    for q = 1:numel(ci)
        i = ci(q); j = cj(q);
        e = [id_u(i, j), id_u(i, j + 1), id_v(i, j), id_v(i + 1, j)];
        % edge order: bottom, top, left, right
        have = e > 0;
        if nnz(have) == 2
            seg(end + 1, :) = e(have); %#ok<AGROW>
        else
            b = e(1); t = e(2); l = e(3); r = e(4);
            um = 0.5 * (U(i) + U(i + 1));
            vm = 0.5 * (U(j) + U(j + 1));
            zc = patch_F_z(g, patch_P(model, g, um), vm, model) - zk;
            if (zc >= 0) == pos(i, j)
                seg(end + 1, :) = [b, r]; %#ok<AGROW>
                seg(end + 1, :) = [t, l]; %#ok<AGROW>
            else
                seg(end + 1, :) = [l, b]; %#ok<AGROW>
                seg(end + 1, :) = [r, t]; %#ok<AGROW>
            end
        end
    end

    nb = zeros(cnt, 2);
    deg = zeros(cnt, 1);
    for q = 1:size(seg, 1)
        for c = 1:2
            a = seg(q, c); b = seg(q, 3 - c);
            deg(a) = deg(a) + 1;
            nb(a, deg(a)) = b;
        end
    end
    seen = false(cnt, 1);
    chains = cell(1, 0);
    starts = find(deg == 1)';
    for s0 = [starts, find(deg == 2)']
        if seen(s0), continue; end
        chain = s0; seen(s0) = true;
        prev = 0; cur = s0; closed = false;
        while true
            nxt = nb(cur, nb(cur, :) ~= prev & nb(cur, :) > 0);
            if isempty(nxt), break; end
            nxt = nxt(1);
            if nxt == s0, closed = true; break; end
            if seen(nxt), break; end
            seen(nxt) = true;
            chain(end + 1) = nxt; %#ok<AGROW>
            prev = cur; cur = nxt;
        end
        chains{end + 1} = struct('idx', chain, 'closed', closed); %#ok<AGROW>
    end

    for c = 1:numel(chains)
        idx = chains{c}.idx(:);
        pts = zeros(numel(idx), 3);
        for q = 1:numel(idx)
            if idx(q) <= nu
                P = Pu{idx(q)};
            else
                P = grid_P(g, gi(idx(q)), U);
            end
            pts(q, :) = patch_F(g, P, xv(idx(q)), model);
        end
        keep = [true; any(diff(pts, 1, 1) ~= 0, 2)];
        pieces{end + 1} = struct('pts', pts(keep, :), 'u', xu(idx(keep)), ...
                                 'v', xv(idx(keep)), 'closed', chains{c}.closed, ...
                                 'patch', 0); %#ok<AGROW>
    end
end

function z = edge_z_u(model, g, u, v)
    S = patch_F(g, patch_P(model, g, u), v, model);
    z = S(3);
end

function z = patch_F_z(g, P, v, model)
    S = patch_F(g, P, v, model);
    z = S(3);
end

function x = solve_edge(fun, a, b, fa, fb, tolz, x0)
% Root of fun on [a,b] (fun(a) = fa, fun(b) = fb of opposite sign, zero counts as >= 0).
% Illinois regula falsi; an optional start x0 (a root of a neighbouring edge) is tried first.
% It stops when successive iterates agree to 2 ulp: the z evaluation has a rounding floor above tolz.
    if fa == 0, x = a; return; end
    if fb == 0, x = b; return; end
    if ~isnan(x0) && x0 > a && x0 < b
        f0 = fun(x0);
        if abs(f0) <= tolz, x = x0; return; end
        if (f0 > 0) == (fa > 0)
            a = x0; fa = f0;
        else
            b = x0; fb = f0;
        end
    end
    side = 0;
    x_old = NaN;
    for it = 1:80
        x = (fa * b - fb * a) / (fa - fb);
        if abs(x - x_old) <= 2 * eps(abs(x))
            return;
        end
        x_old = x;
        fx = fun(x);
        if abs(fx) <= tolz || abs(b - a) <= eps(max(abs(a), abs(b)))
            return;
        end
        if (fx > 0) == (fb > 0)
            b = x; fb = fx;
            if side == -1, fa = fa / 2; end
            side = -1;
        else
            a = x; fa = fx;
            if side == 1, fb = fb / 2; end
            side = 1;
        end
    end
end

function row = join_pieces(pieces, zk)
    row = empty_row(zk);
    m = numel(pieces);
    if m == 0
        return;
    end
    if m == 1 && pieces{1}.closed
        order = 1; entry = 0; gaps = 0;
    else
        if any(cellfun(@(p) p.closed, pieces))
            error('mwecmass:solid:SectionNotClosed', ...
                  'Section at z = %g has more than one loop.', zk);
        end
        E = zeros(2 * m, 3);
        for p = 1:m
            E(2 * p - 1, :) = pieces{p}.pts(1, :);
            E(2 * p, :) = pieces{p}.pts(end, :);
        end
        D = sqrt((E(:, 1) - E(:, 1)').^2 + (E(:, 2) - E(:, 2)').^2 + (E(:, 3) - E(:, 3)').^2);
        D(1:2*m+1:end) = Inf;
        [dmin, partner] = min(D, [], 2);
        if any(partner(partner) ~= (1:2*m)')
            error('mwecmass:solid:SectionNotClosed', ...
                  'Section at z = %g: patch boundary points do not pair up.', zk);
        end
        order = zeros(1, 0); entry = zeros(1, 0); gaps = zeros(1, 0);
        e = 1;
        for it = 1:m
            p = ceil(e / 2);
            order(end + 1) = p; %#ok<AGROW>
            entry(end + 1) = mod(e, 2) == 0; %#ok<AGROW>
            ex = 2 * p - 1 + mod(e, 2);
            gaps(end + 1) = dmin(ex); %#ok<AGROW>
            e = partner(ex);
            if e == 1, break; end
        end
        if e ~= 1 || numel(unique(order)) ~= m
            error('mwecmass:solid:SectionNotClosed', ...
                  'Section at z = %g has more than one loop.', zk);
        end
    end
    row = assemble(pieces, order, entry, gaps, zk);
end

function row = assemble(pieces, order, entry, gaps, zk)
    row = empty_row(zk);
    P = zeros(0, 3); patch = zeros(0, 1); u = zeros(0, 1); v = zeros(0, 1);
    seam = false(0, 1); patch2 = zeros(0, 1); u2 = zeros(0, 1); v2 = zeros(0, 1);
    for t = 1:numel(order)
        pc = pieces{order(t)};
        pts = pc.pts; pu = pc.u; pv = pc.v;
        if entry(t)
            pts = flipud(pts); pu = flipud(pu); pv = flipud(pv);
        end
        k = size(pts, 1);
        P = [P; pts]; %#ok<AGROW>
        patch = [patch; repmat(pc.patch, k, 1)]; %#ok<AGROW>
        u = [u; pu]; v = [v; pv]; %#ok<AGROW>
        seam = [seam; false(k, 1)]; %#ok<AGROW>
        patch2 = [patch2; zeros(k, 1)]; %#ok<AGROW>
        u2 = [u2; NaN(k, 1)]; v2 = [v2; NaN(k, 1)]; %#ok<AGROW>
    end
    if numel(order) > 1 || ~pieces{order(1)}.closed
        starts = cumsum([1, cellfun(@(p) size(p.pts, 1), pieces(order(1:end-1)))]);
        ends = [starts(2:end) - 1, size(P, 1)];
        drop = false(size(P, 1), 1);
        for t = 1:numel(order)
            a = ends(t);
            nxt = mod(t, numel(order)) + 1;
            b = starts(nxt);
            seam(a) = true;
            patch2(a) = patch(b); u2(a) = u(b); v2(a) = v(b);
            drop(b) = true;
        end
        keep = ~drop;
        P = P(keep, :); patch = patch(keep); u = u(keep); v = v(keep);
        seam = seam(keep); patch2 = patch2(keep); u2 = u2(keep); v2 = v2(keep);
    end
    nx = [2:size(P, 1), 1];
    dup = all(P == P(nx, :), 2);
    if any(dup)
        P = P(~dup, :); patch = patch(~dup); u = u(~dup); v = v(~dup);
        seam = seam(~dup); patch2 = patch2(~dup); u2 = u2(~dup); v2 = v2(~dup);
        nx = [2:size(P, 1), 1];
    end
    A = 0.5 * sum(P(:, 1) .* P(nx, 2) - P(nx, 1) .* P(:, 2));
    if A < 0
        P = flipud(P); patch = flipud(patch); u = flipud(u); v = flipud(v);
        seam = flipud(seam); patch2 = flipud(patch2); u2 = flipud(u2); v2 = flipud(v2);
        A = -A;
    end
    row.n = size(P, 1);
    row.pts = P; row.patch = patch; row.u = u; row.v = v;
    row.seam = seam; row.patch2 = patch2; row.u2 = u2; row.v2 = v2;
    row.area = A;
    row.seam_gap = max(gaps);
    row.degenerate = false;
end
