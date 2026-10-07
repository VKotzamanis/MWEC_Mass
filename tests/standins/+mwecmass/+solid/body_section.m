function sec = body_section(body, z, side)
%BODY_SECTION  Stand-in of contract F6b: outer and void loops of a stand-in body at height z.
%
%   sec = mwecmass.solid.body_section(body, z)
%   sec = mwecmass.solid.body_section(body, z, side)
%
%   body: S4 from the build_body stand-in (body.analytic set, else mwecmass:standin:NotAnalytic).
%   The loops are slices (slice_bspline_surface stand-in) of the body's lateral B-rep faces;
%   pieces(k).patch is the index of the brep face. side 'above' (the default) takes the faces with
%   z_lo <= z < z_hi, so a plane face height belongs to the faces above it (z_max to the faces
%   below it), and the module with edges(i) <= z < edges(i+1) (the last module at z_max); 'below'
%   takes the faces with z_lo < z <= z_hi (z_min: those above it) and the module with
%   edges(i) < z <= edges(i+1) (the first module at z_min). The two sides can differ only at a
%   module edge (the module, and the inner loop where the voids of the two modules differ), at
%   z_ballast (the inner loop: the void on the side where it is, none on the solid side) and at an
%   end of a void (the height of a constant-z piece of an inner set: the void loop on the void's
%   side only); elsewhere they agree. Neither fixture has a flat part of the outer surface inside
%   its z range. inner is empty and solid is true where the module has no air at z on that side.

if ~isstruct(body) || ~isfield(body, 'analytic') || isempty(body.analytic)
    error('mwecmass:standin:NotAnalytic', 'body_section stand-in: body.analytic is empty (not a stand-in body)');
end
if nargin < 3
    side = 'above';
end
if ~any(strcmp(side, {'above', 'below'}))
    error('mwecmass:solid:BadSide', 'body_section: side must be ''above'' or ''below''');
end
above = strcmp(side, 'above');
e = body.design.edges(:);
zmin = e(1);
zmax = e(end);
if z < zmin || z > zmax
    error('mwecmass:solid:ZOutside', 'body_section: z = %g outside the hull', z);
end
if above
    m = find(e(1:end - 1) <= z & z < e(2:end), 1);
    if isempty(m)
        m = numel(e) - 1;
    end
else
    m = find(e(1:end - 1) < z & z <= e(2:end), 1);
    if isempty(m)
        m = 1;
    end
end
faces = body.brep.faces;
pick_outer = [];
pick_inner = [];
for k = 1:numel(faces)
    f = faces(k);
    s = body.brep.surfaces{f.surface};
    if ~strcmp(s.type, 'bspline')
        continue
    end
    Z = s.ctrl(:, 1, 3);
    lo = min(Z);
    hi = max(Z);
    if lo == hi
        continue
    end
    if above
        inside = (lo <= z && z < hi) || (z == zmax && hi == zmax);
    else
        inside = (lo < z && z <= hi) || (z == zmin && lo == zmin);
    end
    if ~inside
        continue
    end
    if strcmp(f.role, 'outer')
        pick_outer(end + 1) = k; %#ok<AGROW>
    elseif strcmp(f.role, 'inner') && f.module(1) == m
        pick_inner(end + 1) = k; %#ok<AGROW>
    end
end
sec = struct('z', z, 'module', m, 'outer', slice(pick_outer), 'inner', [], 'solid', isempty(pick_inner));
if ~isempty(pick_inner)
    sec.inner = slice(pick_inner);
end

    function loop = slice(idx)
        pats = struct('surf', cell(1, numel(idx)));
        for q = 1:numel(idx)
            pats(q).surf = body.brep.surfaces{faces(idx(q)).surface};
        end
        loop = mwecmass.solid.slice_bspline_surface(pats, z);
        for p = 1:numel(loop.pieces)
            loop.pieces(p).patch = idx(loop.pieces(p).patch);
        end
    end
end
