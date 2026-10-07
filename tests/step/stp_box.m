function [b, shell] = stp_box(b, lo, hi)
%STP_BOX  Add an axis-aligned box from planar faces; shell is its signed face index vector.
%
%   Vertices, edges and faces that already exist in b are reused (shared faces get a negative
%   sign in the second box).

c = zeros(2, 2, 2);
for i = 1:2
    for j = 1:2
        for k = 1:2
            xyz = [lo(1) + (i - 1) * (hi(1) - lo(1)), lo(2) + (j - 1) * (hi(2) - lo(2)), lo(3) + (k - 1) * (hi(3) - lo(3))];
            [b, c(i, j, k)] = stp_vertex(b, xyz);
        end
    end
end
shell = zeros(1, 6);
n = 0;
for a = 1:3
    for s = [-1 1]
        n = n + 1;
        ib = mod(a, 3) + 1;
        ic = mod(a + 1, 3) + 1;
        corners = zeros(1, 4);
        pat = [0 0; 1 0; 1 1; 0 1];
        for q = 1:4
            idx = ones(1, 3);
            idx(a) = 1 + (s > 0);
            idx(ib) = 1 + pat(q, 1);
            idx(ic) = 1 + pat(q, 2);
            corners(q) = c(idx(1), idx(2), idx(3));
        end
        nrm = zeros(1, 3);
        nrm(a) = s;
        if s < 0
            corners = corners([1 4 3 2]);
        end
        [b, shell(n)] = stp_cap(b, corners, nrm);
    end
end
end
