function sec = body_section(body, z)
%BODY_SECTION  Test mock for test_realised_section_gap: a one-module body whose section along y = 0
%   changes type at the heights body.mock.z_from.
%
%   Outer loop: the rectangle |x| <= 1, |y| <= 1. From z_from(k) up to z_from(k+1) the section is of
%   type body.mock.type(k): 0 solid, 1 hollow with the inner rectangle |x| <= 0.5, |y| <= 0.3,
%   2 hollow with one U-shaped inner loop, |x| <= 0.8, |y| <= 0.3, less the notch |x| < 0.4, y > -0.1,
%   whose y = 0 section is the two intervals 0.4 <= |x| <= 0.8. Every void slice is one simple closed
%   loop (M3).
if z < body.design.edges(1) || z > body.design.edges(end)
    error('mwecmass:solid:ZOutside', 'body_section mock: z = %g outside the hull', z);
end
outer = rectangle_loop(-1, 1, -1, 1);
type = body.mock.type(find(z >= body.mock.z_from, 1, 'last'));
inner = [];
switch type
    case 1
        inner = rectangle_loop(-0.5, 0.5, -0.3, 0.3);
    case 2
        inner = polygon_loop([-0.8 -0.3; 0.8 -0.3; 0.8 0.3; 0.4 0.3; 0.4 -0.1; -0.4 -0.1; -0.4 0.3; -0.8 0.3]);
end
sec = struct('z', z, 'module', 1, 'outer', outer, 'inner', inner, 'solid', type == 0);
end

function loop = rectangle_loop(x0, x1, y0, y1)
loop = polygon_loop([x0 y0; x1 y0; x1 y1; x0 y1]);
end

function loop = polygon_loop(V)
% Counter-clockwise loop of straight pieces through the vertices V [n x 2].
P = [V; V(1, :)];
pieces = struct('curve', {}, 'dir', {});
for k = 1:size(V, 1)
    c = struct('degree', 1, 'knots', [0 0 1 1], 'ctrl', [P(k, :), 0; P(k + 1, :), 0], 'weights', [1; 1]);
    pieces(k) = struct('curve', c, 'dir', 1);
end
loop = struct('pieces', pieces);
end
