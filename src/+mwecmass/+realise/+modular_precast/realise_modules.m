function [ev, ctx] = realise_modules(ctx, design)
%REALISE_MODULES  Build the UHPC modules of a design on the exact kernel and integrate them.
%
%   [ev, ctx] = mwecmass.realise.modular_precast.realise_modules(ctx, design)
%
%   ctx: kernel context of solve_and_extract (geo, rho, z_hollow, inner_fn, knots_from, sets).
%   design: S3 (mode 'modular_precast'). One S2 inner set per distinct finite t of the modules
%   that are not in design.solid_modules (contract section 7 item 2): taken from ctx.sets when a
%   set of that t (bitwise) is there, else made by ctx.inner_fn(t, ctx.knots_from) over the whole
%   hollow range ctx.z_hollow and kept in ctx.sets. With ctx.knots_from empty only adaptive sets
%   (refit false) are used, so the stored design is always built on its own adaptive fit; with
%   ctx.knots_from set only sets of its pieces and knot vectors.
%   ev.body: S4 (mwecmass.solid.build_body), ev.bp: S6 (mwecmass.solid.body_properties with
%   ctx.rho), ev.inner: the S2 sets the body uses.

N = numel(design.edges) - 1;
open = true(N, 1);
open(design.solid_modules) = false;
ts = unique(design.t(open & isfinite(design.t(:))));
inner = [];
for k = 1:numel(ts)
    [set, ctx] = inner_set(ctx, ts(k));
    inner = [inner, set]; %#ok<AGROW>
end
ev.body = mwecmass.solid.build_body(ctx.geo, design, inner);
ev.bp = mwecmass.solid.body_properties(ev.body, ctx.rho, struct());
ev.inner = inner;
end

function [set, ctx] = inner_set(ctx, t)
for k = 1:numel(ctx.sets)
    s = ctx.sets(k);
    if isequal(s.t, t)
        if isempty(ctx.knots_from)
            usable = ~s.refit;
        else
            usable = same_knots(s, ctx.knots_from);
        end
        if usable
            set = s;
            return
        end
    end
end
set = ctx.inner_fn(t, ctx.knots_from);
ctx.sets = [ctx.sets, set];
end

function same = same_knots(a, b)
same = numel(a.patches) == numel(b.patches);
for p = 1:numel(a.patches)
    if ~same
        return
    end
    same = isequal(a.patches(p).surf.degree, b.patches(p).surf.degree) && ...
        isequal(a.patches(p).surf.knots, b.patches(p).surf.knots);
end
end
