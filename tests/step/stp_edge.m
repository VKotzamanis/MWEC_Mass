function [b, e] = stp_edge(b, v1, v2, curve)
%STP_EDGE  Signed edge index for traversing from vertex v1 to v2.
%
%   The edge between the two vertices is added with the given curve (running v1 to v2; a straight
%   degree-1 line when empty) unless it exists; the existing curve is then kept. The sign is
%   negative when the stored edge runs from v2 to v1.

for k = 1:numel(b.edges)
    s = b.edges(k).vertices;
    if s(1) == v1 && s(2) == v2
        e = k;
        return
    elseif s(1) == v2 && s(2) == v1
        e = -k;
        return
    end
end
if isempty(curve)
    curve = struct('degree', 1, 'ctrl', b.vertices([v1 v2], :), 'knots', [0 0 1 1], 'weights', []);
end
b.curves(end + 1) = curve;
b.edges(end + 1) = struct('vertices', [v1 v2], 'curve', numel(b.curves));
e = numel(b.edges);
end
