function [n, orient, limit] = surface_normals(model, cache, patch, u, v, orient)
%SURFACE_NORMALS Outward unit normals S_u x S_v / |S_u x S_v| of hull patches.
%   [n, orient, limit] = SURFACE_NORMALS(model, cache, patch, u, v) for the
%   points (patch(k), u(k), v(k)); patch is a name, a cellstr of names or
%   indices into model.visible_surfs (as in outer_rows rows(k).patch). Scalars
%   broadcast. n is [N x 3].
%
%   S_u and S_v are the analytic derivatives of MS2Parser.eval_surface_with_derivs
%   (central differences of eval_surface, parameter step 1e-6, only for surface
%   types the parser has no derivatives for; C1 uses none).
%
%   Orientation: for each patch the sign of S_u x S_v is fixed once by the exact
%   closed sections of outer_rows. Their counter-clockwise order gives the
%   outward horizontal direction (t_y, -t_x), and the sign is chosen so that
%   the horizontal part of the normal points that way (this covers mirror
%   patches, whose parameter orientation is reversed). orient (fields names,
%   sign) can be passed back in to skip that step.
%
%   limit(k) is true where |S_u x S_v| <= sqrt(eps)|S_u||S_v| (keel point, top of
%   the neck): n is then the Richardson limit 2 n(h) - n(2h) along the
%   parameter line leading into the patch, h = 1e-4.

    if nargin < 6 || isempty(orient)
        orient = patch_orientation(model, cache);
    end
    names = model.visible_surfs(:)';
    if isnumeric(patch)
        patch = names(patch);
    elseif ischar(patch)
        patch = {patch};
    end
    N = max([numel(patch), numel(u), numel(v)]);
    if numel(patch) == 1, patch = repmat(patch, 1, N); end
    if numel(u) == 1, u = repmat(u, N, 1); end
    if numel(v) == 1, v = repmat(v, N, 1); end

    n = zeros(N, 3);
    limit = false(N, 1);
    for k = 1:N
        pk = find(strcmp(names, patch{k}), 1);
        [c, ok] = raw_normal(model, patch{k}, u(k), v(k));
        if ~ok
            limit(k) = true;
            c = limit_normal(model, patch{k}, u(k), v(k));
        end
        n(k, :) = orient.sign(pk) * c;
    end
end

function [c, ok] = raw_normal(model, name, u, v)
    [~, Su, Sv] = model.eval_surface_with_derivs(name, u, v);
    c = cross(Su(1, :), Sv(1, :));
    nc = norm(c);
    ok = nc > sqrt(eps) * norm(Su) * norm(Sv);
    if nc > 0
        c = c / nc;
    end
end

function c = limit_normal(model, name, u, v)
    h = 1e-4;
    dirs = [sign(0.5 - u) + (u == 0.5), 0; 0, sign(0.5 - v) + (v == 0.5)];
    for d = 1:2
        [c1, ok1] = raw_normal(model, name, u + dirs(d, 1) * h, v + dirs(d, 2) * h);
        [c2, ok2] = raw_normal(model, name, u + dirs(d, 1) * 2 * h, v + dirs(d, 2) * 2 * h);
        if ok1 && ok2
            c = 2 * c1 - c2;
            c = c / norm(c);
            return;
        end
    end
    error('mwecmass:solid:DegenerateNormal', ...
          'Normal of %s at (u,v) = (%g,%g) has no limit along either parameter line.', ...
          name, u, v);
end

function orient = patch_orientation(model, cache)
    [~, grids] = mwecmass.solid.outer_rows(model, cache, []);
    zr = grids.z_range;
    rows = mwecmass.solid.outer_rows(model, cache, ...
               zr(1) + (1:9)' / 10 * (zr(2) - zr(1)), grids);
    names = model.visible_surfs(:)';
    best = zeros(numel(names), 1);
    for k = 1:numel(rows)
        r = rows(k);
        if r.degenerate, continue; end
        prev = [r.n, 1:r.n - 1];
        next = [2:r.n, 1];
        for p = unique(r.patch)'
            idx = find(r.patch == p & ~r.seam);
            if isempty(idx), continue; end
            i = idx(ceil(numel(idx) / 2));
            t = r.pts(next(i), 1:2) - r.pts(prev(i), 1:2);
            out = [t(2), -t(1)] / norm(t);
            c = raw_normal(model, names{p}, r.u(i), r.v(i));
            cs = dot(c(1:2), out);
            if abs(cs) > abs(best(p))
                best(p) = cs;
            end
        end
    end
    if any(abs(best) < 0.5)
        error('mwecmass:solid:OrientationUndecided', ...
              'Patch %s: no sampled section point has a horizontal normal component.', ...
              names{find(abs(best) < 0.5, 1)});
    end
    orient = struct('names', {names}, 'sign', sign(best));
end
