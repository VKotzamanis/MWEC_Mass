function sec = body_section(body, z, side)
%BODY_SECTION  Outer and void loops of a B-rep body at one height (contract F6b, S5).
%
%   sec = mwecmass.solid.body_section(body, z)
%   sec = mwecmass.solid.body_section(body, z, side)
%
%   body: S4 (build_body); z [m, body]. side: 'above' (default) or 'below'. With 'above' a height
%   belongs to the lateral faces and the module above it (half-open [z_lo, z_hi)), z_max to those
%   below; 'below' takes the faces and the module below it ((z_lo, z_hi]), z_min those above. So
%   at a module edge, at z_ballast and at the height of a flat part the side chooses the module,
%   its void and the faces whose section is drawn. Loops are slice_bspline_surface (F4) of the
%   chosen outer faces and of the chosen inner faces of that module; pieces(k).patch is the index
%   of the B-rep face.
%   sec: z, module, outer (S5; empty at z_min or z_max where the hull ends in a pole or a line),
%   inner (S5; empty where the module has no void at z, or where its void ends in a pole or a line
%   at z), solid (true when inner is empty).

if nargin < 3 || isempty(side)
    side = 'above';
end
above = strcmp(side, 'above');
if ~above && ~strcmp(side, 'below')
    error('mwecmass:solid:BadSide', 'body_section: side must be ''above'' or ''below''');
end
e = body.design.edges(:);
N = numel(e) - 1;
zmin = e(1);
zmax = e(end);
if z < zmin || z > zmax
    error('mwecmass:solid:ZOutside', 'body_section: z = %.17g outside the hull [%.17g, %.17g]', z, zmin, zmax);
end
if (above && z < zmax) || z == zmin
    m = find(e(1:end - 1) <= z & z < e(2:end), 1);
    pick = @(lo, hi) lo <= z && z < hi;
else
    m = find(e(1:end - 1) < z & z <= e(2:end), 1);
    pick = @(lo, hi) lo < z && z <= hi;
end
if isempty(m)
    m = N;
end
B = body.brep;
outer_idx = [];
inner_idx = [];
flat_here = false;
for k = 1:numel(B.faces)
    f = B.faces(k);
    s = B.surfaces{f.surface};
    if ~strcmp(s.type, 'bspline')
        flat_here = flat_here || (strcmp(f.role, 'outer') && s.origin(3) == z);
        continue
    end
    zr = [s.ctrl(1, 1, 3), s.ctrl(end, 1, 3)];
    lo = min(zr);
    hi = max(zr);
    if lo == hi
        flat_here = flat_here || (strcmp(f.role, 'outer') && lo == z);
        continue
    end
    if ~pick(lo, hi)
        continue
    end
    if strcmp(f.role, 'outer')
        outer_idx(end + 1) = k; %#ok<AGROW>
    elseif strcmp(f.role, 'inner') && f.module(1) == m
        inner_idx(end + 1) = k; %#ok<AGROW>
    end
end
outer = [];
if ~((z == zmin || z == zmax) && ~flat_here)
    outer = slice(outer_idx);
end
inner = [];
v = body.voids(m);
open_end = (z == v.z_lo && v.open_lo) || (z == v.z_hi && v.open_hi);
if ~isempty(inner_idx) && ~open_end
    inner = slice(inner_idx);
end
sec = struct('z', z, 'module', m, 'outer', outer, 'inner', inner, 'solid', isempty(inner));

    function loop = slice(idx)
        pats = struct('name', {}, 'surf', {}, 'z_of_u', {}, 'u_range', {}, 'z_range', {}, ...
            'pole', {}, 'c0_u', {}, 'seam_u0', {}, 'seam_u1', {});
        for q = 1:numel(idx)
            s = B.surfaces{B.faces(idx(q)).surface};
            ku = s.knots{1};
            pats(q) = struct('name', sprintf('face %d', idx(q)), 'surf', s, 'z_of_u', true, ...
                'u_range', [ku(1) ku(end)], 'z_range', [s.ctrl(1, 1, 3), s.ctrl(end, 1, 3)], ...
                'pole', [collapsed(s.ctrl(1, :, :)), collapsed(s.ctrl(end, :, :))], 'c0_u', [], ...
                'seam_u0', [], 'seam_u1', []);
        end
        loop = mwecmass.solid.slice_bspline_surface(pats, z);
        for p = 1:numel(loop.pieces)
            loop.pieces(p).patch = idx(loop.pieces(p).patch);
        end
    end
end

function tf = collapsed(r)
r = reshape(r, [], 3);
tf = all(all(r == r(1, :)));
end
