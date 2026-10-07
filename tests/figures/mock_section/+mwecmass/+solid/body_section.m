function sec = body_section(body, z)
%BODY_SECTION  Test mock for test_realised_section_gap: a one-module body whose section along y = 0
%   changes type at the heights body.mock.z_from.
%
%   Outer loop: the rectangle |x| <= 1, |y| <= 1. From z_from(k) up to z_from(k+1) the section is of
%   type body.mock.type(k): 0 solid, 1 hollow with one inner rectangle |x| <= 0.5, |y| <= 0.3,
%   2 hollow with two, 0.4 <= |x| <= 0.8, |y| <= 0.3.
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
        a = rectangle_loop(-0.8, -0.4, -0.3, 0.3);
        b = rectangle_loop(0.4, 0.8, -0.3, 0.3);
        inner = struct('pieces', [a.pieces, b.pieces]);
end
sec = struct('z', z, 'module', 1, 'outer', outer, 'inner', inner, 'solid', type == 0);
end

function loop = rectangle_loop(x0, x1, y0, y1)
P = [x0 y0; x1 y0; x1 y1; x0 y1; x0 y0];
pieces = struct('curve', {}, 'dir', {});
for k = 1:4
    c = struct('degree', 1, 'knots', [0 0 1 1], 'ctrl', [P(k, :), 0; P(k + 1, :), 0], 'weights', [1; 1]);
    pieces(k) = struct('curve', c, 'dir', 1);
end
loop = struct('pieces', pieces);
end
