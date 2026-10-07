function [b, lateral, caps, volume] = stp_bulged_box(b)
%STP_BULGED_BOX  Unit cube whose four side faces are bicubic non-rational B-spline patches.
%
%   Each side is the cube face with its interior control points moved outward along the face
%   normal; the patch boundary rows stay on the cube edges. The control points of the in-plane
%   coordinates sit at the Greville abscissae, so the in-plane map is the identity and the
%   patch is a height field of the displacement. Its volume is then 1 plus the integral of the
%   displacements, from the basis integrals int N_i = (t(i+p+1) - t(i)) / (p + 1), which uses
%   none of the writer or of OpenCASCADE.

p = 3;
t = [0 0 0 0 0.5 1 1 1 1];
n = numel(t) - p - 1;
greville = arrayfun(@(i) mean(t(i + 1:i + p)), 1:n);
basis_integral = arrayfun(@(i) (t(i + p + 1) - t(i)) / (p + 1), 1:n);
pattern = [1 2 1; 3 4 2; 1 2 3] * 0.01;
volume = 1;
lateral = zeros(1, 4);
for f = 1:4
    d = zeros(n, n);
    d(2:4, 2:4) = pattern * f;
    volume = volume + basis_integral * d * basis_integral';
    ctrl = zeros(n, n, 3);
    for i = 1:n
        for j = 1:n
            a = greville(i);
            c = greville(j);
            switch f
                case 1
                    xyz = [1 + d(i, j), a, c];
                case 2
                    xyz = [c, 1 + d(i, j), a];
                case 3
                    xyz = [-d(i, j), c, a];
                case 4
                    xyz = [a, -d(i, j), c];
            end
            ctrl(i, j, :) = xyz;
        end
    end
    [b, lateral(f)] = stp_patch(b, [p p], ctrl, {t, t}, []);
end
caps = zeros(1, 2);
v = zeros(1, 4);
for q = 1:4
    [b, v(q)] = stp_vertex(b, [mod(q - 1, 3) > 0, q > 2, 0]);
end
[b, caps(1)] = stp_cap(b, v([1 4 3 2]), [0 0 -1]);
for q = 1:4
    [b, v(q)] = stp_vertex(b, [mod(q - 1, 3) > 0, q > 2, 1]);
end
[b, caps(2)] = stp_cap(b, v, [0 0 1]);
end
