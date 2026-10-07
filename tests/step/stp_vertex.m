function [b, i] = stp_vertex(b, xyz)
%STP_VERTEX  Index of the vertex at xyz, added if no vertex has exactly these coordinates.

xyz = xyz(:)';
i = find(all(b.vertices == xyz, 2), 1);
if isempty(i)
    b.vertices(end + 1, :) = xyz;
    i = size(b.vertices, 1);
end
end
