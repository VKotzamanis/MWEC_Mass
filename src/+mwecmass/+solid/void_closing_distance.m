function d_close = void_closing_distance(model, cache, geo, z_range)
%VOID_CLOSING_DISTANCE  Offset distance at which offset layers from opposite sides meet (contract F2b).
%
%   d_close = mwecmass.solid.void_closing_distance(model, cache, geo, z_range)
%
%   Two inward offset layers of the outer surface (geo.outer, S1) meet when the offset distance
%   reaches half the length of a double normal: a chord P-X inside the hull that is normal to the
%   surface at both ends, with opposite inward normals (the ball of radius |X - P|/2 about its
%   midpoint touches the surface at two antipodal points). A fold of one layer (a curvature centre,
%   where P and X coincide) is not a double normal. d_close is the smallest half-length of a
%   double normal whose midpoint height lies in z_range [m, body]; Inf if there is none.
%   Double normals are critical points of |X - P|^2 over pairs of surface points: seeds are pairs of
%   sample points (4 per knot span and direction on every face) whose chord is closest to both
%   normals, refined by damped Newton steps on (u, v) of both points to convergence of the
%   gradient (relative 1e-10, the solver's stopping rule); a chord is kept when its midpoint and
%   quarter points lie inside the hull's sections. Where an end of z_range lies inside the hull,
%   layers that converge without a double normal (a funnel) meet first on that end plane: the
%   smallest medial ball centred on the plane is found as well (plane_meeting). The design bound is
%   t_max = d_close - eps_fit/2.
%   model and cache are not used by the exact path (the faces are geo.outer).

persistent store
if isempty(store)
    store = containers.Map();
end
zr = sort(z_range(:)');
key = fingerprint(geo, zr);
if store.isKey(key)
    d_close = store(key);
    return
end
E = geo.outer;
E = E(arrayfun(@(e) ~isempty(e.surf), E));
flats = cell(1, numel(E));
% samples: points, inward unit normals, face and parameters; flat: index of the sample's
% constant-z interval in flats{face} (0 off them)
S = zeros(0, 3);
N = zeros(0, 3);
id = zeros(0, 3);
flat = zeros(0, 1);
for k = 1:numel(E)
    s = E(k).surf;
    flats{k} = flat_intervals(s);
    us = samples(s.knots{1});
    vs = samples(s.knots{2});
    [UU, VV] = ndgrid(us, vs);
    [Pk, Su, Sv] = mwecmass.solid.eval_bspline_surface(s, UU(:), VV(:));
    n = cross(Su, Sv, 2);
    ln = sqrt(sum(n.^2, 2));
    good = ln > 0;
    n = n(good, :) ./ ln(good);
    if E(k).outward
        n = -n;
    end
    fk = zeros(nnz(good), 1);
    ug = UU(good);
    for j = 1:size(flats{k}, 1)
        fk(ug >= flats{k}(j, 1) & ug <= flats{k}(j, 2)) = j;
    end
    S = [S; Pk(good, :)]; %#ok<AGROW>
    N = [N; n]; %#ok<AGROW>
    id = [id; [repmat(k, nnz(good), 1), ug, VV(good)]]; %#ok<AGROW>
    flat = [flat; fk]; %#ok<AGROW>
end
m = size(S, 1);
best = zeros(m, 1);
for i = 1:m
    c = S - S(i, :);
    L = sqrt(sum(c.^2, 2));
    ch = c ./ max(L, realmin);
    ok = L > 0 & ch * N(i, :)' > 0 & sum(ch .* N, 2) < 0 & N * N(i, :)' < 0;
    if ~any(ok)
        continue
    end
    r = sqrt(sum(cross(ch, repmat(N(i, :), m, 1), 2).^2, 2)) + sqrt(sum(cross(ch, N, 2).^2, 2));
    r(~ok) = Inf;
    [~, best(i)] = min(r);
end
pairs = unique(sort([(1:m)', best], 2), 'rows');
pairs = pairs(pairs(:, 1) > 0, :);
[P, X, nP, nX, conv] = refine_all(E, id(pairs(:, 1), :), id(pairs(:, 2), :));
cand = zeros(0, 1);
for q = 1:size(pairs, 1)
    c = X(q, :) - P(q, :);
    L = norm(c);
    mid = (P(q, :) + X(q, :)) / 2;
    if conv(q) && L > 0 && c * nP(q, :)' > 0 && c * nX(q, :)' < 0 && mid(3) >= zr(1) && mid(3) <= zr(2)
        cand(end + 1, 1) = q; %#ok<AGROW>
    end
end
[~, o] = sort(sqrt(sum((X(cand, :) - P(cand, :)).^2, 2)));
cand = cand(o);
d_close = Inf;
for q = cand'
    c = X(q, :) - P(q, :);
    if inside(E, P(q, :) + c / 4) && inside(E, P(q, :) + c / 2) && inside(E, P(q, :) + 3 * c / 4)
        d_close = norm(c) / 2;
        break
    end
end
% where z_range ends inside the hull, layers that converge without a double normal (a funnel)
% meet first on that end plane
zz = reshape([E.z_range], 2, [])';
for h = zr(zr > min(zz(:)) & zr < max(zz(:)))
    d_close = min(d_close, plane_meeting(E, S, N, id, h, d_close));
end
d_close = min(d_close, contact_events(E, flats, S, N, id, flat, zr, d_close));
store(key) = d_close;
end

function key = fingerprint(geo, zr)
f = 0;
g = 0;
for k = 1:numel(geo.outer)
    s = geo.outer(k).surf;
    if isempty(s)
        continue
    end
    x = s.ctrl(:);
    f = f + sum(x .* (1:numel(x))') + k * sum(s.knots{1});
    g = g + sum(x.^2) + numel(x) + geo.outer(k).outward;
end
key = [geo.hull_name, '|', num2hex(f), num2hex(g), '|', num2hex(zr(1)), num2hex(zr(2))];
end

function s = samples(k)
ku = unique(k);
s = zeros(0, 1);
for j = 1:numel(ku) - 1
    s = [s; ku(j) + ((1:4)' - 0.5) / 4 * (ku(j + 1) - ku(j))]; %#ok<AGROW>
end
end

function [P, X, nP, nX, conv] = refine_all(E, A, B)
% critical points of |X(u2, v2) - P(u1, v1)|^2 for every seed pair (rows of A, B: face, u, v):
% damped Newton (Levenberg-Marquardt) on the gradient; nP, nX: inward unit normals at the ends
m = size(A, 1);
x = [A(:, 2:3), B(:, 2:3)];
lo = zeros(m, 4);
hi = zeros(m, 4);
for q = 1:m
    sa = E(A(q, 1)).surf;
    sb = E(B(q, 1)).surf;
    lo(q, :) = [sa.knots{1}(1), sa.knots{2}(1), sb.knots{1}(1), sb.knots{2}(1)];
    hi(q, :) = [sa.knots{1}(end), sa.knots{2}(end), sb.knots{1}(end), sb.knots{2}(end)];
end
conv = false(m, 1);
active = true(m, 1);
for it = 1:60
    [P, Pu, Pv, Puu, Puv, Pvv] = eval_faces(E, A(:, 1), x(:, 1), x(:, 2));
    [X, Xu, Xv, Xuu, Xuv, Xvv] = eval_faces(E, B(:, 1), x(:, 3), x(:, 4));
    r = X - P;
    for q = find(active)'
        f = [r(q, :) * Pu(q, :)'; r(q, :) * Pv(q, :)'; r(q, :) * Xu(q, :)'; r(q, :) * Xv(q, :)'];
        scale = norm(r(q, :)) * max([norm(Pu(q, :)), norm(Pv(q, :)), norm(Xu(q, :)), norm(Xv(q, :)), realmin]);
        if norm(f) <= 1e-10 * scale
            conv(q) = true;
            active(q) = false;
            continue
        end
        a1 = Pu(q, :); a2 = Pv(q, :); b1 = Xu(q, :); b2 = Xv(q, :); rq = r(q, :);
        J = [-a1 * a1' + rq * Puu(q, :)', -a2 * a1' + rq * Puv(q, :)', b1 * a1', b2 * a1'; ...
            -a1 * a2' + rq * Puv(q, :)', -a2 * a2' + rq * Pvv(q, :)', b1 * a2', b2 * a2'; ...
            -a1 * b1', -a2 * b1', b1 * b1' + rq * Xuu(q, :)', b2 * b1' + rq * Xuv(q, :)'; ...
            -a1 * b2', -a2 * b2', b1 * b2' + rq * Xuv(q, :)', b2 * b2' + rq * Xvv(q, :)'];
        lam = 1e-12 * max(trace(J' * J), realmin);
        step = -(J' * J + lam * eye(4)) \ (J' * f);
        xn = min(max(x(q, :) + step', lo(q, :)), hi(q, :));
        if isequal(xn, x(q, :))
            active(q) = false;
        end
        x(q, :) = xn;
    end
    if ~any(active)
        break
    end
end
[P, Pu, Pv] = eval_faces(E, A(:, 1), x(:, 1), x(:, 2));
[X, Xu, Xv] = eval_faces(E, B(:, 1), x(:, 3), x(:, 4));
nP = inward(E, A(:, 1), Pu, Pv);
nX = inward(E, B(:, 1), Xu, Xv);
end

function [S, Su, Sv, Suu, Suv, Svv] = eval_faces(E, f, u, v)
m = numel(f);
S = zeros(m, 3); Su = S; Sv = S; Suu = S; Suv = S; Svv = S;
for k = unique(f)'
    i = f == k;
    [S(i, :), Su(i, :), Sv(i, :), Suu(i, :), Suv(i, :), Svv(i, :)] = ...
        mwecmass.solid.eval_bspline_surface(E(k).surf, u(i), v(i));
end
end

function rho = plane_meeting(E, S, N, id, h, rho_max)
% smallest radius of a ball centred on the plane z = h that touches the surface at two points
% (a medial point of the hull on that plane): unknowns (u1, v1, u2, v2, rho), equations
% P + rho nP = X + rho nX, P_z + rho nP_z = h, and stationarity of rho along the medial curve in
% the plane, ((nP - nX) x e_z) . (nP + nX) = 0; damped Newton with a forward-difference Jacobian
% from seeds (a sample point P, its ball centre on the plane, the sample X whose distance to that
% centre is closest to rho with a normal toward it). A solution counts when the two normals differ
% (not the trivial P = X) and no surface point is closer to the centre than rho, both to the
% resolution of the iteration (sqrt(eps)).
rho = Inf;
m = size(S, 1);
cand = find(N(:, 3) ~= 0);
r0 = (h - S(cand, 3)) ./ N(cand, 3);
keep = r0 > 0 & r0 < rho_max;
cand = cand(keep);
r0 = r0(keep);
if isempty(cand)
    return
end
A = zeros(numel(cand), 3);
B = zeros(numel(cand), 3);
for q = 1:numel(cand)
    i = cand(q);
    c = S(i, :) + r0(q) * N(i, :);
    toward = sum((c - S) .* N, 2) > 0 & sqrt(sum((N - N(i, :)).^2, 2)) > 0.5;
    if ~any(toward)
        A(q, 1) = 0;
        continue
    end
    dev = abs(sqrt(sum((S - c).^2, 2)) - r0(q));
    dev(~toward) = Inf;
    [~, j] = min(dev);
    A(q, :) = id(i, :);
    B(q, :) = id(j, :);
end
ok = A(:, 1) > 0;
A = A(ok, :);
B = B(ok, :);
r0 = r0(ok);
if isempty(r0)
    return
end
x = [A(:, 2:3), B(:, 2:3), r0];
n = size(x, 1);
lo = zeros(n, 5);
hi = zeros(n, 5);
for q = 1:n
    sa = E(A(q, 1)).surf;
    sb = E(B(q, 1)).surf;
    lo(q, :) = [sa.knots{1}(1), sa.knots{2}(1), sb.knots{1}(1), sb.knots{2}(1), 0];
    hi(q, :) = [sa.knots{1}(end), sa.knots{2}(end), sb.knots{1}(end), sb.knots{2}(end), Inf];
end
for it = 1:40
    F0 = medial_eq(E, A, B, x, h);
    J = zeros(n, 5, 5);
    for k = 1:5
        dx = 1e-7 * max(1, abs(x(:, k)));
        xp = x;
        xp(:, k) = xp(:, k) + dx;
        J(:, :, k) = (medial_eq(E, A, B, xp, h) - F0) ./ dx;
    end
    for q = 1:n
        Jq = reshape(J(q, :, :), 5, 5);
        lam = 1e-12 * max(trace(Jq' * Jq), realmin);
        x(q, :) = min(max(x(q, :) - ((Jq' * Jq + lam * eye(5)) \ (Jq' * F0(q, :)'))', lo(q, :)), hi(q, :));
    end
end
[F0, P, X, nP, nX] = medial_eq(E, A, B, x, h);
res = sqrt(sum(F0.^2, 2));
good = res <= sqrt(eps) * (1 + x(:, 5)) & sum(nP .* nX, 2) < 0 & x(:, 5) > 0;
for q = find(good)'
    c = P(q, :) + x(q, 5) * nP(q, :);
    if x(q, 5) < rho && surface_distance(E, c) >= x(q, 5) * (1 - sqrt(eps))
        rho = x(q, 5);
    end
end
end

function [F, P, X, nP, nX] = medial_eq(E, A, B, x, h)
[P, Pu, Pv] = eval_faces(E, A(:, 1), x(:, 1), x(:, 2));
[X, Xu, Xv] = eval_faces(E, B(:, 1), x(:, 3), x(:, 4));
nP = inward(E, A(:, 1), Pu, Pv);
nX = inward(E, B(:, 1), Xu, Xv);
r = x(:, 5);
a = nP - nX;
b = nP + nX;
F = [P + r .* nP - X - r .* nX, P(:, 3) + r .* nP(:, 3) - h, a(:, 2) .* b(:, 1) - a(:, 1) .* b(:, 2)];
end

function dmin = surface_distance(E, X)
% distance from every point of X to the faces E: on every knot-span cell of every face, Newton on
% (u, v) clamped to the cell, with descent safeguard and an active bound, from the three nearest of
% its 9 x 9 grid points (the distance may have several local minima in a cell, e.g. behind a fold;
% a cell is polynomial or rational without interior knots, so the iteration does not cross a
% crease). The grid points of all cells give an upper bound first; a cell is then searched only
% for the points whose distance to its grid's bounding box, grown by the box's size, is below it.
m = size(X, 1);
cells = struct('k', {}, 'lo', {}, 'hi', {}, 'U', {}, 'V', {}, 'G', {});
dmin = Inf(m, 1);
for k = 1:numel(E)
    s = E(k).surf;
    if isempty(s)
        continue
    end
    ku = unique(s.knots{1});
    kv = unique(s.knots{2});
    for a = 1:numel(ku) - 1
        for b = 1:numel(kv) - 1
            % at an interior knot the evaluator takes the next span: stay just below it
            hi = [ku(a + 1) - (a + 1 < numel(ku)) * eps(ku(a + 1)), kv(b + 1) - (b + 1 < numel(kv)) * eps(kv(b + 1))];
            us = ku(a) + (0:8)' / 8 * (hi(1) - ku(a));
            vs = kv(b) + (0:8)' / 8 * (hi(2) - kv(b));
            [UU, VV] = ndgrid(us, vs);
            G = mwecmass.solid.eval_bspline_surface(s, UU(:), VV(:));
            cells(end + 1) = struct('k', k, 'lo', [ku(a) kv(b)], 'hi', hi, 'U', UU(:), 'V', VV(:), 'G', G); %#ok<AGROW>
            D = (X(:, 1) - G(:, 1)').^2 + (X(:, 2) - G(:, 2)').^2 + (X(:, 3) - G(:, 3)').^2;
            dmin = min(dmin, sqrt(min(D, [], 2)));
        end
    end
end
for c = cells
    s = E(c.k).surf;
    G = c.G;
    ext = max(max(G, [], 1) - min(G, [], 1));
    lb = sqrt(sum(max(0, max(min(G, [], 1) - X, X - max(G, [], 1))).^2, 2));
    idx = find(lb <= dmin + ext);
    if isempty(idx)
        continue
    end
    Xi = X(idx, :);
    D = (Xi(:, 1) - G(:, 1)').^2 + (Xi(:, 2) - G(:, 2)').^2 + (Xi(:, 3) - G(:, 3)').^2;
    [~, order] = sort(D, 2);
    for start = 1:3
        j = order(:, start);
        u = c.U(j);
        v = c.V(j);
        dcur = sqrt(sum((mwecmass.solid.eval_bspline_surface(s, u, v) - Xi).^2, 2));
        for it = 1:25
            [S, Su, Sv, Suu, Suv, Svv] = mwecmass.solid.eval_bspline_surface(s, u, v);
            r = S - Xi;
            g1 = sum(r .* Su, 2);
            g2 = sum(r .* Sv, 2);
            h11 = sum(Su .* Su, 2) + sum(r .* Suu, 2);
            h12 = sum(Su .* Sv, 2) + sum(r .* Suv, 2);
            h22 = sum(Sv .* Sv, 2) + sum(r .* Svv, 2);
            bad = ~(h11 .* h22 - h12.^2 > 0 & h11 > 0);
            h11(bad) = sum(Su(bad, :).^2, 2) + realmin;
            h22(bad) = sum(Sv(bad, :).^2, 2) + realmin;
            h12(bad) = 0;
            dt = h11 .* h22 - h12.^2;
            du = -(h22 .* g1 - h12 .* g2) ./ dt;
            dv = -(h11 .* g2 - h12 .* g1) ./ dt;
            % a parameter held at a bound of the cell: Newton in the other one alone
            hu = (u <= c.lo(1) & du < 0) | (u >= c.hi(1) & du > 0);
            hv = (v <= c.lo(2) & dv < 0) | (v >= c.hi(2) & dv > 0);
            du(hv & ~hu) = -g1(hv & ~hu) ./ h11(hv & ~hu);
            dv(hv) = 0;
            dv(hu & ~hv) = -g2(hu & ~hv) ./ h22(hu & ~hv);
            du(hu) = 0;
            du(~isfinite(du)) = 0;
            dv(~isfinite(dv)) = 0;
            % descent: a step toward a stationary point that is not a minimum is halved until the
            % distance decreases, or dropped
            lam = ones(size(u));
            for h = 1:12
                un = min(max(u + lam .* du, c.lo(1)), c.hi(1));
                vn = min(max(v + lam .* dv, c.lo(2)), c.hi(2));
                dn = sqrt(sum((mwecmass.solid.eval_bspline_surface(s, un, vn) - Xi).^2, 2));
                worse = dn > dcur;
                if ~any(worse)
                    break
                end
                lam(worse) = lam(worse) / 2;
            end
            un(worse) = u(worse);
            vn(worse) = v(worse);
            dn(worse) = dcur(worse);
            if isequal(un, u) && isequal(vn, v)
                break
            end
            u = un;
            v = vn;
            dcur = dn;
        end
        dmin(idx) = min(dmin(idx), dcur);
    end
end
end

function n = inward(E, f, Su, Sv)
n = cross(Su, Sv, 2);
n = n ./ max(sqrt(sum(n.^2, 2)), realmin);
out = arrayfun(@(k) E(k).outward, f);
n(out, :) = -n(out, :);
end

function tf = inside(E, Q)
zz = reshape([E.z_range], 2, [])';
if Q(3) <= min(zz(:)) || Q(3) >= max(zz(:))
    tf = false;
    return
end
lat = zz(:, 1) ~= zz(:, 2);
try
    L = mwecmass.solid.slice_bspline_surface(E(lat), Q(3));
catch
    tf = false;
    return
end
tf = inpolygon(Q(1), Q(2), L.pts(:, 1), L.pts(:, 2));
end

function F = flat_intervals(s)
% u-intervals (whole knot spans) on which every control point acting there has one z, bitwise:
% the face is the plane z = F(:, 3) there and its normal is +-e_z exactly. Rows [ua ub z].
p = s.degree(1);
U = s.knots{1};
Z = s.ctrl(:, :, 3);
nu = size(Z, 1);
ku = unique(U);
F = zeros(0, 3);
for a = 1:numel(ku) - 1
    rows = find(U(1:nu) < ku(a + 1) & U(p + 2:nu + p + 1) > ku(a));
    z = Z(rows, :);
    if all(z(:) == z(1))
        if ~isempty(F) && F(end, 2) == ku(a) && F(end, 3) == z(1)
            F(end, 2) = ku(a + 1);
        else
            F(end + 1, :) = [ku(a) ku(a + 1) z(1)]; %#ok<AGROW>
        end
    end
end
end

function rho = contact_events(E, flats, S, N, id, flat, zr, rho_cap)
% smallest radius of a ball inscribed in the hull, centred in z_range, that touches the surface at
% three or four points c = q_i + rho n_i (n_i inward unit normals) of one of these kinds:
%   critical: 0 lies in the convex hull of the n_i (three coplanar normals, det = 0, or four): the
%     void pinches or vanishes there without a double normal;
%   flat face: one contact on a constant-z interval (n = +-e_z) and 0 in the convex hull of the
%     horizontal parts of the others (two antiparallel ones with det = 0, or three): that face's
%     offset shrinks to nothing (or splits) because the layers of the other contacts meet over it.
% A fold of one layer is neither: its contacts lie on one connected arc of normals, within an open
% hemisphere, and a flat face does not fold. A sign of a convex weight that rounding decides
% belongs to a configuration that reduces to fewer contacts, a double normal at the same radius.
% Seeds: the ball of each sample point that touches the surface there and passes through another
% sample point, smallest radius (at most twice rho_cap, the smallest event found so far), with the
% sample points nearest to its sphere, one per normal direction; refined by Levenberg-Marquardt
% steps on the contact parameters, the centre and the radius to the solver's stopping rule
% (residual relative 1e-10 of the coordinates involved). The ball counts when no surface point is
% closer to its centre than rho, to the same stopping rule.
rho = Inf;
m = size(S, 1);
cap = 2 * rho_cap;
seeds = zeros(0, 6);
c0 = zeros(0, 3);
r0 = zeros(0, 1);
for i = 1:m
    dS = S - S(i, :);
    w = dS * N(i, :)';
    ok = w > 0;
    if ~any(ok)
        continue
    end
    ri = min(sum(dS(ok, :).^2, 2) ./ (2 * w(ok)));
    if ri > cap
        continue
    end
    c = S(i, :) + ri * N(i, :);
    dev = sqrt(sum((S - c).^2, 2)) - ri;
    dev(sum((c - S) .* N, 2) <= 0) = Inf;
    dev(i) = Inf;
    [ds, o] = sort(dev);
    o = o(isfinite(ds));
    pick = i;
    for j = o'
        if numel(pick) == 6
            break
        end
        if all(N(pick, :) * N(j, :)' < cos(pi / 6))
            pick(end + 1) = j; %#ok<AGROW>
        end
    end
    others = pick(2:end);
    for k = 3:4
        if numel(others) < k - 1
            continue
        end
        sub = nchoosek(others, k - 1);
        for r = 1:size(sub, 1)
            q = [i, sub(r, :)];
            if event_kind(N(q, :), flat(q) > 0)
                seeds(end + 1, :) = [sort(q), zeros(1, 4 - k), k, 0]; %#ok<AGROW>
                c0(end + 1, :) = c; %#ok<AGROW>
                r0(end + 1, 1) = ri; %#ok<AGROW>
            end
        end
    end
end
if isempty(seeds)
    return
end
[~, u] = unique(seeds(:, 1:5), 'rows');
seeds = seeds(u, :);
c0 = c0(u, :);
r0 = r0(u);
for k = 3:4
    sel = seeds(:, 5) == k;
    if ~any(sel)
        continue
    end
    q = seeds(sel, 1:k);
    [x, ok, nn, L] = refine_contacts(E, flats, id, flat, N, q, c0(sel, :), r0(sel), k);
    c = x(:, 2 * k + 1:2 * k + 3);
    r = x(:, end);
    good = ok & r > 0 & r < min(rho, rho_cap) & c(:, 3) >= zr(1) & c(:, 3) <= zr(2);
    for g = find(good)'
        good(g) = event_kind(reshape(nn(g, :, :), k, 3), flat(q(g, :)) > 0);
    end
    if ~any(good)
        continue
    end
    g = find(good);
    dmin = surface_distance(E, c(g, :));
    in = dmin >= r(g) - 1e-10 * L(g);
    if any(in)
        rho = min(rho, min(r(g(in))));
    end
end
end

function tf = event_kind(n, isflat)
% n: [k x 3] inward unit normals of the contacts; isflat: contacts on constant-z intervals
k = size(n, 1);
tf = same_sign(null_vector(n'));
if tf
    return
end
for f = find(isflat(:)')
    h = n(setdiff(1:k, f), 1:2);
    if k == 3
        tf = h(1, :) * h(2, :)' < 0;
    else
        tf = same_sign(null_vector(h'));
    end
    if tf
        return
    end
end
end

function x = null_vector(A)
[~, ~, V] = svd(A);
x = V(:, end);
end

function tf = same_sign(x)
tf = all(x > 0) || all(x < 0);
end

function [x, conv, nn, L] = refine_contacts(E, flats, id, flat, N, q, c0, r0, k)
% Levenberg-Marquardt on x = [u_1 v_1 ... u_k v_k, c, rho] for q_i + rho n_i = c (and, for three
% contacts, det(n_1, n_2, n_3) = 0: rho stationary along the curve of balls touching all three).
% A contact seeded on a constant-z interval stays on it (its u clamped to the interval) and takes
% the plane's exact height and normal.
n = size(q, 1);
nx = 2 * k + 4;
x = zeros(n, nx);
F = id(q, 1);
F = reshape(F, n, k);
FL = false(n, k);
ZF = zeros(n, k);
SG = zeros(n, k);
lo = -Inf(n, nx);
hi = Inf(n, nx);
for j = 1:k
    x(:, 2 * j - 1) = id(q(:, j), 2);
    x(:, 2 * j) = id(q(:, j), 3);
    for r = 1:n
        s = E(F(r, j)).surf;
        lo(r, 2 * j - 1:2 * j) = [s.knots{1}(1), s.knots{2}(1)];
        hi(r, 2 * j - 1:2 * j) = [s.knots{1}(end), s.knots{2}(end)];
        fi = flat(q(r, j));
        if fi > 0
            iv = flats{F(r, j)}(fi, :);
            lo(r, 2 * j - 1) = iv(1);
            hi(r, 2 * j - 1) = iv(2);
            FL(r, j) = true;
            ZF(r, j) = iv(3);
            SG(r, j) = sign(N(q(r, j), 3));
        end
    end
end
x(:, 2 * k + 1:2 * k + 3) = c0;
x(:, end) = r0;
lo(:, end) = 0;
[R, J] = k_system(E, x, F, FL, ZF, SG, k);
% damping relative to the largest squared singular value of J; steps by the SVD, so that a
% singular J (contacts on a continuum, as on a surface of revolution) takes the minimum-norm step
lam = 1e-12 * ones(n, 1);
active = true(n, 1);
for it = 1:200
    % iterate until no step lowers the residual, then judge the stopping rule
    active = active & lam < 1;
    if ~any(active)
        break
    end
    xn = x;
    for r = find(active)'
        [Uj, Sj, Vj] = svd(reshape(J(r, :, :), [], nx), 'econ');
        sj = diag(Sj);
        step = -Vj * ((sj ./ (sj.^2 + lam(r) * sj(1)^2)) .* (Uj' * R(r, :)'));
        xn(r, :) = min(max(x(r, :) + step', lo(r, :)), hi(r, :));
    end
    a = find(active);
    Rn = k_system(E, xn(a, :), F(a, :), FL(a, :), ZF(a, :), SG(a, :), k);
    better = sum(Rn.^2, 2) < sum(R(a, :).^2, 2);
    b = a(better);
    x(b, :) = xn(b, :);
    lam(b) = max(lam(b) / 10, eps);
    lam(a(~better)) = lam(a(~better)) * 100;
    if any(better)
        [R(b, :), J(b, :, :)] = k_system(E, x(b, :), F(b, :), FL(b, :), ZF(b, :), SG(b, :), k);
    end
end
[R, ~, nn] = k_system(E, x, F, FL, ZF, SG, k);
[conv, L] = converged(R, x, k);
end

function [conv, L] = converged(R, x, k)
c = x(:, 2 * k + 1:2 * k + 3);
P = zeros(size(x, 1), 3 * k);
for j = 1:k
    P(:, 3 * j - 2:3 * j) = c - R(:, 3 * j - 2:3 * j);
end
L = max(abs([c, P, x(:, end)]), [], 2);
conv = sqrt(sum(R(:, 1:3 * k).^2, 2)) <= 1e-10 * L;
if size(R, 2) > 3 * k
    conv = conv & abs(R(:, end)) <= 1e-10;
end
end

function [R, J, nn] = k_system(E, x, F, FL, ZF, SG, k)
n = size(x, 1);
nx = 2 * k + 4;
ne = 3 * k + (k == 3);
ow = double([E.outward]);
c = x(:, 2 * k + 1:2 * k + 3);
rho = x(:, end);
R = zeros(n, ne);
J = zeros(n, ne, nx);
nn = zeros(n, k, 3);
nU = zeros(n, k, 3);
nV = zeros(n, k, 3);
for j = 1:k
    [P, Pu, Pv, Puu, Puv, Pvv] = eval_faces(E, F(:, j), x(:, 2 * j - 1), x(:, 2 * j));
    m = cross(Pu, Pv, 2);
    Lm = sqrt(sum(m.^2, 2));
    nh = m ./ Lm;
    mu = cross(Puu, Pv, 2) + cross(Pu, Puv, 2);
    mv = cross(Puv, Pv, 2) + cross(Pu, Pvv, 2);
    nu = (mu - nh .* sum(nh .* mu, 2)) ./ Lm;
    nv = (mv - nh .* sum(nh .* mv, 2)) ./ Lm;
    sg = 1 - 2 * reshape(ow(F(:, j)), [], 1);
    nh = nh .* sg;
    nu = nu .* sg;
    nv = nv .* sg;
    fl = FL(:, j);
    if any(fl)
        P(fl, 3) = ZF(fl, j);
        Pu(fl, 3) = 0;
        Pv(fl, 3) = 0;
        nh(fl, :) = [zeros(nnz(fl), 2), SG(fl, j)];
        nu(fl, :) = 0;
        nv(fl, :) = 0;
    end
    rows = 3 * j - 2:3 * j;
    R(:, rows) = P + rho .* nh - c;
    J(:, rows, 2 * j - 1) = Pu + rho .* nu;
    J(:, rows, 2 * j) = Pv + rho .* nv;
    for a = 1:3
        J(:, rows(a), 2 * k + a) = -1;
    end
    J(:, rows, nx) = nh;
    nn(:, j, :) = reshape(nh, n, 1, 3);
    nU(:, j, :) = reshape(nu, n, 1, 3);
    nV(:, j, :) = reshape(nv, n, 1, 3);
end
if k == 3
    g = @(a, b, c3) sum(a .* cross(b, c3, 2), 2);
    n1 = reshape(nn(:, 1, :), n, 3);
    n2 = reshape(nn(:, 2, :), n, 3);
    n3 = reshape(nn(:, 3, :), n, 3);
    R(:, ne) = g(n1, n2, n3);
    J(:, ne, 1) = g(reshape(nU(:, 1, :), n, 3), n2, n3);
    J(:, ne, 2) = g(reshape(nV(:, 1, :), n, 3), n2, n3);
    J(:, ne, 3) = g(n1, reshape(nU(:, 2, :), n, 3), n3);
    J(:, ne, 4) = g(n1, reshape(nV(:, 2, :), n, 3), n3);
    J(:, ne, 5) = g(n1, n2, reshape(nU(:, 3, :), n, 3));
    J(:, ne, 6) = g(n1, n2, reshape(nV(:, 3, :), n, 3));
end
end
