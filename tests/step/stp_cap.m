function [b, f] = stp_cap(b, corners, normal)
%STP_CAP  Signed index of the planar face through the vertices in corners.
%
%   corners are vertex indices in counter-clockwise order seen from normal. Edges between
%   consecutive corners are reused or added as straight lines. If a face on the same vertices
%   exists, it is reused and the sign is negative when its normal opposes the given one.

key = sort(corners);
for k = 1:numel(b.faces)
    s = b.surfaces{b.faces(k).surface};
    if ~strcmp(s.type, 'plane')
        continue
    end
    L = b.faces(k).loops{1};
    first = arrayfun(@(e) b.edges(abs(e)).vertices(1 + (e < 0)), L);
    if isequal(sort(first), key)
        f = k * sign(s.normal * normal(:));
        return
    end
end
n = numel(corners);
loop = zeros(1, n);
for q = 1:n
    [b, loop(q)] = stp_edge(b, corners(q), corners(mod(q, n) + 1), []);
end
b.surfaces{end + 1} = struct('type', 'plane', 'origin', b.vertices(corners(1), :), 'normal', normal(:)');
b.faces(end + 1) = struct('surface', numel(b.surfaces), 'same_sense', true, 'loops', {{loop}});
f = numel(b.faces);
end
