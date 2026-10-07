function ev = evaluate_design_point(ctx, vs, t, z_ballast, adaptive)
%EVALUATE_DESIGN_POINT  Thin-shell design (vs, t, z_ballast) evaluated on the exact body.
%
%   ev = mwecmass.realise.thin_shell.evaluate_design_point(ctx, vs, t, z_ballast, adaptive)
%
%   ctx is built by mwecmass.realise.thin_shell.solve: config, geo (S1b), edges, rho (ballast,
%   shell, air), stage2 (S8.stage2), t_min, tol_eq and the handle caches sets, adaptive_sets,
%   hydro and state. vs [m], t [m] (one shell thickness for every module), z_ballast [m, body].
%   adaptive = true takes the inner set fitted adaptively at t (F2 without knots_from), the set a
%   reported design is built from; otherwise the set is refitted on the knot vectors of
%   ctx.state('ref'), so the properties stay differentiable in t within one solve (contract
%   section 7 item 4). Inner sets are cached by the bitwise value of t, hydrostatics by the
%   bitwise value of vs.
%
%   ev: design (S3), inner (S2), body (S4), bp (S6), hs (S7), props (F9), check (F10),
%   objective = sum over X in {Z_CG = CG_total(3), coupled T_heave, coupled T_pitch} of
%   ((X3 - X2)/X2)^2 (AGENTS section 3 item 27), ceq = [flotation; GM] residuals of F10.

N = numel(ctx.edges) - 1;
ev = struct();
ev.design = struct('mode', 'thin_shell', 'edges', ctx.edges, 'vs', vs, 't', t * ones(N, 1), ...
    'z_ballast', z_ballast, 'solid_modules', []);
ev.inner = inner_set(ctx, t, adaptive);
ev.body = mwecmass.solid.build_body(ctx.geo, ev.design, ev.inner);
ev.bp = mwecmass.solid.body_properties(ev.body, ctx.rho, struct());
key = num2hex(vs);
if isKey(ctx.hydro, key)
    ev.hs = ctx.hydro(key);
else
    ev.hs = mwecmass.solid.hydrostatics_at_draft(ctx.geo, vs, struct());
    ctx.hydro(key) = ev.hs;
end
ev.props = mwecmass.realise.evaluate_realised(ev.bp, ev.hs, ev.design, ctx.config);
ev.check = mwecmass.realise.check_against_stage2(ev.props, ctx.stage2, ...
    ctx.config.mass_acceptable_pct, ctx.tol_eq, ctx.config.RHO_WATER);
s2 = ctx.stage2;
x3 = [ev.props.CG_total(3), ev.props.periods.heave, ev.props.periods.pitch];
x2 = [s2.Z_CG, s2.T_heave, s2.T_pitch];
ev.objective = sum(((x3 - x2) ./ x2).^2);
ev.ceq = [ev.check.equalities.residual]';
end

function fitted = inner_set(ctx, t, adaptive)
key = num2hex(t);
if adaptive
    store = ctx.adaptive_sets;
else
    store = ctx.sets;
end
if isKey(store, key)
    fitted = store(key);
    return
end
opts = struct('t_min', ctx.t_min);
if ~adaptive
    opts.knots_from = ctx.state('ref');
end
fitted = mwecmass.solid.offset_surface(ctx.config.ms2_model, ctx.config.boundary_cache, ctx.geo, t, ...
    ctx.geo.z_range, opts);
store(key) = fitted; %#ok<NASGU> containers.Map is a handle: the assignment fills the cache
end
