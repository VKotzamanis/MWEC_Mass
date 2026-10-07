function [n, orient, limit] = surface_normals(model, cache, patch, u, v, orient)
%SURFACE_NORMALS Outward unit normals S_u x S_v / |S_u x S_v| of hull patches.
%   [n, orient, limit] = SURFACE_NORMALS(model, cache, patch, u, v) for the
%   points (patch(k), u(k), v(k)); patch is a name, a cellstr of names or
%   indices into model.visible_surfs (as in outer_rows rows(k).patch). Scalars
%   broadcast. n is [N x 3].
%
%   S_u and S_v come from MS2Parser.eval_surface_with_derivs (analytic for every
%   surface type the parser evaluates; its central-difference fallback uses a
%   parameter step of 1e-6).
%
%   Orientation: the sign of S_u x S_v is fixed once per patch (mirror patches
%   have reversed parameter orientation) by a 3-D step test on the exact closed
%   sections of outer_rows. At sample points p of the patch with raw normal c, a
%   step d along c must leave the hull (outside the section at the stepped
%   height, or beyond the hull's z range) and a step -d must stay inside; d runs
%   through 1e-2, 1e-3, 1e-4 of the hull height. This works for horizontal and
%   near-horizontal patches as well as steep ones. Up to 3 decisive samples per
%   patch must agree, otherwise error mwecmass:solid:OrientationUndecided.
%   orient (fields names, sign) can be passed back in to skip that step.
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
% Sign of S_u x S_v per patch from a 3-D step test on the exact closed sections:
% a decisive sample has p + d*c outside the section at its height (or above/below
% the hull) and p - d*c inside; then c points outward, and the reverse for -c.
    [~, grids] = mwecmass.solid.outer_rows(model, cache, []);
    ctx = struct('model', model, 'cache', cache, 'grids', grids, ...
                 'rows', containers.Map('KeyType', 'double', 'ValueType', 'any'));
    H = grids.z_range(2) - grids.z_range(1);
    steps = H * [1e-2, 1e-3, 1e-4];
    samples = [0.5 0.5; 0.3 0.7; 0.7 0.3; 0.2 0.2; 0.8 0.8; 0.1 0.6; 0.6 0.1];
    names = model.visible_surfs(:)';
    sgn = zeros(numel(names), 1);
    for p = 1:numel(names)
        votes = zeros(1, 0);
        for q = 1:size(samples, 1)
            [c, ok] = raw_normal(model, names{p}, samples(q, 1), samples(q, 2));
            if ~ok, continue; end
            P = model.eval_surface(names{p}, samples(q, 1), samples(q, 2));
            for d = steps
                [ctx, out_p] = in_hull(ctx, P + d * c);
                [ctx, in_p] = in_hull(ctx, P - d * c);
                if out_p == 0 && in_p == 1
                    votes(end + 1) = 1; %#ok<AGROW>
                    break;
                elseif out_p == 1 && in_p == 0
                    votes(end + 1) = -1; %#ok<AGROW>
                    break;
                end
            end
            if numel(votes) >= 3, break; end
        end
        if isempty(votes) || any(votes ~= votes(1))
            error('mwecmass:solid:OrientationUndecided', ...
                  'Patch %s: the step test along S_u x S_v is %s.', names{p}, ...
                  'inconsistent between samples or never decisive');
        end
        sgn(p) = votes(1);
    end
    orient = struct('names', {names}, 'sign', sgn);
end

function [ctx, r] = in_hull(ctx, q)
% 1 inside, 0 outside the exact section at the height of q; NaN if the section cannot be built.
    z = q(3);
    if z <= ctx.grids.z_range(1) || z >= ctx.grids.z_range(2)
        r = 0;
        return;
    end
    if isKey(ctx.rows, z)
        row = ctx.rows(z);
    else
        try
            row = mwecmass.solid.outer_rows(ctx.model, ctx.cache, z, ctx.grids);
        catch err
            if ~strcmp(err.identifier, 'mwecmass:solid:SectionNotClosed'), rethrow(err); end
            row = [];
        end
        ctx.rows(z) = row;
    end
    if isempty(row)
        r = NaN;
    else
        r = double(inpolygon(q(1), q(2), row.pts(:, 1), row.pts(:, 2)));
    end
end
