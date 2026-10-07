function [ev, ctx] = realise_modules(ctx, design)
%REALISE_MODULES  Build the UHPC modules of a design on the exact kernel and integrate them.
%
%   [ev, ctx] = mwecmass.realise.modular_precast.realise_modules(ctx, design)
%
%   ctx: kernel context of solve_and_extract (geo, rho, z_hollow, t_max_hollow, inner_fn,
%   knots_from, sets, set_ranges). design: S3 (mode 'modular_precast'). One S2 inner set per
%   distinct finite t of the modules that are not in design.solid_modules (contract section 7
%   item 2), made by ctx.inner_fn(t, ctx.knots_from, z_range) and kept in ctx.sets with its z_range
%   in ctx.set_ranges. z_range is the whole hollow range ctx.z_hollow while t < ctx.t_max_hollow
%   (one set serves every module of that t); a thicker t closes the void somewhere in the hollow
%   range, so its set covers only the span of the modules that use it, whose own bound t_max,i
%   allows it. A kept set of that t (bitwise) whose z_range covers the span is reused. With
%   ctx.knots_from empty only adaptive sets (refit false) are used, so the stored design is always
%   built on its own adaptive fit; with ctx.knots_from set only sets of its pieces and knot vectors.
%   ev.body: S4 (mwecmass.solid.build_body), ev.bp: S6 (mwecmass.solid.body_properties with
%   ctx.rho), ev.inner: the S2 sets the body uses.

e = design.edges(:);
N = numel(e) - 1;
has_void = true(N, 1);
has_void(design.solid_modules) = false;
has_void = has_void & isfinite(design.t(:));
ts = unique(design.t(has_void));
inner = [];
for k = 1:numel(ts)
    users = find(has_void & design.t(:) == ts(k));
    span = [e(users(1)), e(users(end) + 1)];
    [set, ctx] = inner_set(ctx, ts(k), span);
    inner = [inner, set]; %#ok<AGROW>
end
ev.body = mwecmass.solid.build_body(ctx.geo, design, inner);
ev.bp = mwecmass.solid.body_properties(ev.body, ctx.rho, struct());
ev.inner = inner;
end

function [set, ctx] = inner_set(ctx, t, span)
for k = 1:numel(ctx.sets)
    s = ctx.sets(k);
    zr = ctx.set_ranges(k, :);
    if isequal(s.t, t) && zr(1) <= span(1) && zr(2) >= span(2)
        if isempty(ctx.knots_from)
            usable = ~s.refit;
        else
            usable = mwecmass.realise.modular_precast.same_knots(s, ctx.knots_from);
        end
        if usable
            set = s;
            return
        end
    end
end
z_range = ctx.z_hollow;
if t >= ctx.t_max_hollow
    z_range = span;
end
set = ctx.inner_fn(t, ctx.knots_from, z_range);
ctx.sets = [ctx.sets, set];
ctx.set_ranges = [ctx.set_ranges; z_range];
end
