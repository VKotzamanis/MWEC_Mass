function inner = sti_inner_box(geo, t, t_min)
%STI_INNER_BOX  Exact S2 inner set of the box fixture at design thickness t.
%
%   inner = sti_inner_box(geo, t, t_min)
%
%   geo is the S1b of box.ms2 (stand-in or real F1). The six planar faces are moved inward by
%   d = t + eps_fit/2, eps_fit = 0.01 t_min (contract section 0), and trimmed to each other: the
%   inner box [x0+d, x1-d] x [y0+d, y1-d] x [z0+d, z1-d] (source: the normal offset of a plane is
%   the parallel plane; contract section 3, Stand-in kit SK). Same patch layout, parameters and
%   seams as the outer faces, normals into the void (closed void); each patch carries the
%   `visible` of the outer face it is offset from. d >= d_close errors mwecmass:solid:VoidClosed.

fx = sti_closed_form('fixture', geo);
if ~strcmp(fx.kind, 'box')
    error('sti_inner_box: %s is not the box fixture', fx.name);
end
eps_fit = 0.01 * t_min;
d = t + eps_fit / 2;
if d >= sti_closed_form('d_close', fx)
    error('mwecmass:solid:VoidClosed', 'sti_inner_box: d = %g m closes the void', d);
end
patches = sti_closed_form('patches', fx, d);
for k = 1:numel(patches)
    patches(k).name = [patches(k).name '_inner'];
    patches(k).outward = ~patches(k).outward;
end
pr = struct('n_nodes', 0, 'n_knots', [4 4], 'n_passes', 0, 'n_removed', 0, 'n_check', 0, ...
    't_local_min', d, 't_local_max', d, 'M1', d >= t_min, 'M2', t <= d && d <= t + eps_fit, ...
    'M3', true, 'M3_reason', '');
pr = repmat(pr, 1, numel(patches));
rep = struct('patches', pr, 'ok', pr(1).M1 && pr(1).M2, 'cap_reached', false);
inner = struct('t', t, 'd', d, 'eps_fit', eps_fit, 'z_range', [fx.z(1) + d, fx.z(2) - d], ...
    'z_lo', fx.z(1) + d, 'refit', false, 'patches', patches, 'flat', [], 'report', rep);
end
