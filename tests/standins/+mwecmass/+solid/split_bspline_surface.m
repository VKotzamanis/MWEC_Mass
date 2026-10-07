function [lo, hi, u_star] = split_bspline_surface(patch, z)
%SPLIT_BSPLINE_SURFACE  Stand-in of contract F3b: cut a z_of_u patch of u-degree 1 at height z.
%
%   [lo, hi, u_star] = mwecmass.solid.split_bspline_surface(patch, z)
%
%   Closed form for u-degree 1 (the fixtures): z(u) on the span between control rows i and i+1 with
%   row weights a_i, a_(i+1) is (1-s) a_i z_i + s a_(i+1) z_(i+1) over (1-s) a_i + s a_(i+1), so
%   s = a_i (z - z_i) / (a_i (z - z_i) + a_(i+1) (z_(i+1) - z)); knot insertion of u* to
%   multiplicity 1 adds the row interpolated with alpha = s in homogeneous coordinates (Piegl &
%   Tiller, The NURBS Book, 2nd ed., eq. (5.15)). lo and hi keep the parent parameter and share the
%   cut row bitwise; the cut boundaries get empty seams. A patch of u-degree other than 1 errors
%   mwecmass:standin:NotAnalytic. Errors as F3b: ZNotOneParameter, ZOutside, ZNotMonotonic.

surf = patch.surf;
if surf.degree(1) ~= 1
    error('mwecmass:standin:NotAnalytic', 'split_bspline_surface stand-in: u-degree %d (closed form for degree 1 only)', surf.degree(1));
end
[nu, nv, ~] = size(surf.ctrl);
Z = surf.ctrl(:, :, 3);
W = surf.weights;
if isempty(W)
    W = ones(nu, nv);
end
if ~patch.z_of_u || any(any(Z ~= Z(:, 1))) || any(any(W .* W(1, 1) ~= W(:, 1) * W(1, :)))
    error('mwecmass:solid:ZNotOneParameter', 'split_bspline_surface: %s is not z_of_u', patch.name);
end
Z = Z(:, 1);
if ~((z > Z(1) && z < Z(end)) || (z < Z(1) && z > Z(end)))
    error('mwecmass:solid:ZOutside', 'split_bspline_surface: z = %g not strictly inside (%g, %g)', z, Z(1), Z(end));
end
dZ = diff(Z);
if any(dZ > 0) && any(dZ < 0)
    error('mwecmass:solid:ZNotMonotonic', 'split_bspline_surface: z(u) of %s is not monotonic', patch.name);
end
t = surf.knots{1}(:)';
hit = find(Z == z);
if ~isempty(hit)
    if numel(hit) > 1
        error('mwecmass:solid:ZNotMonotonic', 'split_bspline_surface: z(u) = %g on an interval', z);
    end
    i = hit;
    u_star = t(i + 1);
    lo_ctrl = surf.ctrl(1:i, :, :);
    lo_w = W(1:i, :);
    lo_k = [t(1:i + 1) u_star];
    hi_ctrl = surf.ctrl(i:end, :, :);
    hi_w = W(i:end, :);
    hi_k = [u_star t(i + 1:end)];
else
    i = find((Z(1:end - 1) < z & Z(2:end) > z) | (Z(1:end - 1) > z & Z(2:end) < z), 1);
    a0 = W(i, 1);
    a1 = W(i + 1, 1);
    s = a0 * (z - Z(i)) / (a0 * (z - Z(i)) + a1 * (Z(i + 1) - z));
    u_star = t(i + 1) + s * (t(i + 2) - t(i + 1));
    [row, wrow] = interp_row(surf.ctrl(i, :, :), surf.ctrl(i + 1, :, :), W(i, :), W(i + 1, :), s);
    lo_ctrl = cat(1, surf.ctrl(1:i, :, :), row);
    lo_w = [W(1:i, :); wrow];
    lo_k = [t(1:i + 1) u_star u_star];
    hi_ctrl = cat(1, row, surf.ctrl(i + 1:end, :, :));
    hi_w = [wrow; W(i + 1:end, :)];
    hi_k = [u_star u_star t(i + 2:end)];
end
if isempty(surf.weights)
    lo_w = [];
    hi_w = [];
end
lo = patch;
hi = patch;
lo.surf.ctrl = lo_ctrl;
lo.surf.weights = lo_w;
lo.surf.knots{1} = lo_k;
hi.surf.ctrl = hi_ctrl;
hi.surf.weights = hi_w;
hi.surf.knots{1} = hi_k;
lo.u_range = [patch.u_range(1) u_star];
hi.u_range = [u_star patch.u_range(2)];
lo.z_range = [patch.z_range(1) z];
hi.z_range = [z patch.z_range(2)];
lo.pole = [patch.pole(1) false];
hi.pole = [false patch.pole(2)];
lo.c0_u = patch.c0_u(patch.c0_u > lo.u_range(1) & patch.c0_u < u_star);
hi.c0_u = patch.c0_u(patch.c0_u > u_star & patch.c0_u < hi.u_range(2));
lo.seam_u1 = [];
hi.seam_u0 = [];
end

function [row, wrow] = interp_row(Pa, Pb, wa, wb, s)
if isequal(wa, wb)
    row = (1 - s) * Pa + s * Pb;
    wrow = wa;
else
    wrow = (1 - s) * wa + s * wb;
    row = ((1 - s) * Pa .* wa + s * Pb .* wb) ./ wrow;
end
end
