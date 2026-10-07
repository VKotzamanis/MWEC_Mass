function [inner, rep] = offset_surface(model, cache, geo, t, z_range, opts) %#ok<INUSL>
%OFFSET_SURFACE  Stand-in of contract F2: the exact S2 inner set of a stand-in fixture.
%
%   [inner, rep] = mwecmass.solid.offset_surface(model, cache, geo, t, z_range, opts)
%
%   geo must carry geo.analytic (stand-in F1 of cylinder.ms2), else mwecmass:standin:NotAnalytic.
%   opts.t_min sets eps_fit = 0.01 t_min and d = t + eps_fit/2 (contract section 0). The inner
%   surface is the cylinder of radius R - d between z0 + d and z1 - d: source, the normal offset of
%   the side by d and of the end disks by d, which meet at the sharp rims (convex creases) without
%   a fold (contract section 3, Stand-in kit SK). The rims are u-creases, so every quarter patch
%   gives three crease pieces (bottom disk, side, top disk), each its own patch in the parent's
%   parameter, normals into the void. The set covers the whole void [z0 + d, z1 - d]; z_range only
%   has to lie in the hull. The fit report is exact: t_local = d everywhere, no nodes, no passes.
%   d >= d_close errors mwecmass:solid:VoidClosed. For the box the set is sti_inner_box (its six
%   planar faces shifted inward by d), as the real F2 must build it. Each piece carries the `visible`
%   of the outer entry it is offset from; flat is empty. opts.knots_from keeps the same pieces and
%   knot vectors (refit = true).

if ~isstruct(geo) || ~isfield(geo, 'analytic') || isempty(geo.analytic)
    error('mwecmass:standin:NotAnalytic', 'offset_surface stand-in: geo.analytic is empty (not a stand-in fixture)');
end
standin_path();
fx = geo.analytic;
zr = sort(z_range(:)');
if numel(zr) ~= 2 || zr(1) < fx.z(1) || zr(2) > fx.z(2)
    error('mwecmass:solid:ZOutside', 'offset_surface stand-in: z_range [%g %g] outside the hull', z_range);
end
refit = nargin >= 6 && isfield(opts, 'knots_from') && ~isempty(opts.knots_from);
if strcmp(fx.kind, 'box')
    inner = sti_inner_box(geo, t, opts.t_min);
    inner.refit = refit;
    rep = inner.report;
    return
end
eps_fit = 0.01 * opts.t_min;
d = t + eps_fit / 2;
d_close = sti_closed_form('d_close', fx);
if d >= d_close
    error('mwecmass:solid:VoidClosed', 'offset_surface stand-in: d = %g m >= d_close = %g m', d, d_close);
end
full = sti_closed_form('patches', fx, d);
rows = {[1 2], [2 3], [3 4]};
u = {[0 1/3], [1/3 2/3], [2/3 1]};
n = numel(full);
patches = repmat(full(1), 1, 3 * n);
for k = 1:n
    f = full(k);
    nb0 = f.seam_v0(1);
    nb1 = f.seam_v1(1);
    for p = 1:3
        q = f;
        idx = 3 * (k - 1) + p;
        q.name = sprintf('%s_inner%d', f.name, p);
        q.visible = geo.outer(k).visible;
        q.surf.ctrl = f.surf.ctrl(rows{p}, :, :);
        q.surf.weights = f.surf.weights(rows{p}, :);
        q.surf.knots{1} = u{p}([1 1 2 2]);
        q.outward = ~f.outward;
        q.exact = false;
        q.u_range = u{p};
        zz = [f.surf.ctrl(rows{p}(1), 1, 3), f.surf.ctrl(rows{p}(2), 1, 3)];
        q.z_range = zz;
        q.pole = [p == 1, p == 3];
        q.c0_u = [];
        q.c0_v = [];
        q.seam_u0 = [];
        q.seam_u1 = [];
        if p > 1
            q.seam_u0 = [idx - 1 2];
        end
        if p < 3
            q.seam_u1 = [idx + 1 4];
        end
        q.seam_v0 = [3 * (nb0 - 1) + p, 1];
        q.seam_v1 = [3 * (nb1 - 1) + p, 3];
        patches(idx) = q;
    end
end
pr = struct('n_nodes', 0, 'n_knots', [6 6], 'n_passes', 0, 'n_removed', 0, 'n_check', 0, ...
    't_local_min', d, 't_local_max', d, 'M1', d >= opts.t_min, 'M2', t <= d && d <= t + eps_fit, ...
    'M3', true, 'M3_reason', '');
for k = 1:numel(patches)
    pr(k) = pr(1);
    pr(k).n_knots = [numel(patches(k).surf.knots{1}), numel(patches(k).surf.knots{2})];
end
rep = struct('patches', pr, 'ok', all([pr.M1]) && all([pr.M2]) && all([pr.M3]), 'cap_reached', false);
inner = struct('t', t, 'd', d, 'eps_fit', eps_fit, 'z_range', [fx.z(1) + d, fx.z(2) - d], ...
    'z_lo', fx.z(1) + d, 'refit', refit, 'patches', patches, 'flat', [], 'report', rep);
end

function standin_path()
if exist('sti_closed_form', 'file') ~= 2
    addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'fixtures'), '-end');
end
end
