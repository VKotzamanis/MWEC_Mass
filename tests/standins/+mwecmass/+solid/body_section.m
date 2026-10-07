function sec = body_section(body, z)
%BODY_SECTION  Stand-in of contract F6b: outer and void loops of a stand-in body at height z.
%
%   sec = mwecmass.solid.body_section(body, z)
%
%   body: S4 from the build_body stand-in (body.analytic set, else mwecmass:standin:NotAnalytic).
%   The loops are slices (slice_bspline_surface stand-in) of the body's lateral B-rep faces that
%   contain z, taken half-open [z_lo, z_hi) so a plane face height belongs to the face above it
%   (z_max to the faces below it); pieces(k).patch is the index of the brep face. module is the
%   module with edges(i) <= z < edges(i+1) (the last module at z_max); inner is empty and solid is
%   true where the module has no air at z.

if ~isstruct(body) || ~isfield(body, 'analytic') || isempty(body.analytic)
    error('mwecmass:standin:NotAnalytic', 'body_section stand-in: body.analytic is empty (not a stand-in body)');
end
e = body.design.edges(:);
zmax = e(end);
if z < e(1) || z > zmax
    error('mwecmass:solid:ZOutside', 'body_section: z = %g outside the hull', z);
end
m = find(e(1:end - 1) <= z & z < e(2:end), 1);
if isempty(m)
    m = numel(e) - 1;
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
    inside = (lo <= z && z < hi) || (z == zmax && hi == zmax);
    if ~inside
        continue
    end
    if strcmp(f.role, 'outer')
        pick_outer(end + 1) = k; %#ok<AGROW>
    elseif strcmp(f.role, 'inner') && f.module(1) == m && z < hi
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
