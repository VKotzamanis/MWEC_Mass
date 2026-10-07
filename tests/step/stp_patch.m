function [b, f] = stp_patch(b, degree, ctrl, knots, weights)
%STP_PATCH  Add a B-spline surface face bounded by the four boundary curves of its control net.
%
%   The face normal is S_u x S_v. Corners are the net corners; each boundary curve is the
%   boundary row of the net, so no trimming is needed. Edges that exist already are reused.

[nu, nv, ~] = size(ctrl);
if isempty(weights)
    w = @(i, j) [];
else
    w = @(i, j) weights(i, j);
end
corner = @(i, j) squeeze(ctrl(i, j, :))';
v = zeros(1, 4);
[b, v(1)] = stp_vertex(b, corner(1, 1));
[b, v(2)] = stp_vertex(b, corner(nu, 1));
[b, v(3)] = stp_vertex(b, corner(nu, nv));
[b, v(4)] = stp_vertex(b, corner(1, nv));

cu = @(j) struct('degree', degree(1), 'ctrl', squeeze(ctrl(:, j, :)), 'knots', knots{1}, 'weights', col(w(1:nu, j)));
cv = @(i) struct('degree', degree(2), 'ctrl', squeeze(ctrl(i, :, :)), 'knots', knots{2}, 'weights', col(w(i, 1:nv)));
loop = zeros(1, 4);
[b, loop(1)] = stp_edge(b, v(1), v(2), cu(1));
[b, loop(2)] = stp_edge(b, v(2), v(3), cv(nu));
[b, e] = stp_edge(b, v(4), v(3), cu(nv));
loop(3) = -e;
[b, e] = stp_edge(b, v(1), v(4), cv(1));
loop(4) = -e;

b.surfaces{end + 1} = struct('type', 'bspline', 'degree', degree, 'ctrl', ctrl, 'knots', {knots}, 'weights', weights);
b.faces(end + 1) = struct('surface', numel(b.surfaces), 'same_sense', true, 'loops', {{loop}});
f = numel(b.faces);
end

function c = col(w)
c = w(:);
end
