function out = hull_outline_data(geo, n_z)
%HULL_OUTLINE_DATA Exact y = 0 outline of the hull, from its outer surface patches.
%
%   out = mwecmass.output.figures.hull_outline_data(geo, n_z)
%
%   geo: S1b hull geometry (config.hull_solid). The outline is the set of points where the section
%   loops of the outer patches (mwecmass.solid.slice_bspline_surface) cross y = 0, found on the
%   exact curves (section_y0_crossings), at n_z heights (default 121) spread evenly from 8 ulp of the
%   largest |z| above the keel to the same distance below the top (a pole is not sectioned exactly
%   at its point). A height at which the kernel reports no section is listed in out.omitted; a
%   section with more than two crossings errors mwecmass:figures:SectionTopology (the outline would
%   need several intervals).
%   Fields, body frame [m]: z [n x 1], x_lo, x_hi [n x 1], profile [2n x 2] = [x z] counter-clockwise
%   (right side up, left side down), z_range, omitted (struct array z, reason).

if nargin < 2 || isempty(n_z)
    n_z = 121;
end
z_min = geo.z_range(1);
z_max = geo.z_range(2);
margin = 8 * eps(max(abs([z_min z_max])));
levels = linspace(z_min + margin, z_max - margin, n_z)';
z = zeros(0, 1);
x_lo = zeros(0, 1);
x_hi = zeros(0, 1);
omitted = struct('z', {}, 'reason', {});
for k = 1:n_z
    zk = levels(k);
    [loop, reason] = try_slice(geo.outer, zk);
    if isempty(loop)
        omitted(end + 1) = struct('z', zk, 'reason', reason); %#ok<AGROW>
        continue
    end
    xs = mwecmass.output.figures.section_y0_crossings(loop);
    if isempty(xs)
        omitted(end + 1) = struct('z', zk, 'reason', 'the section does not reach y = 0'); %#ok<AGROW>
        continue
    end
    if numel(xs) > 2
        error('mwecmass:figures:SectionTopology', ...
            'hull_outline_data: the section at z = %g crosses y = 0 at %d points', zk, numel(xs));
    end
    z(end + 1, 1) = zk; %#ok<AGROW>
    x_lo(end + 1, 1) = xs(1); %#ok<AGROW>
    x_hi(end + 1, 1) = xs(end); %#ok<AGROW>
end
out = struct('z', z, 'x_lo', x_lo, 'x_hi', x_hi, 'profile', [x_hi, z; flipud(x_lo), flipud(z)], ...
    'z_range', [z_min z_max], 'omitted', omitted);
end

function [loop, reason] = try_slice(patches, z)
loop = [];
reason = '';
try
    loop = mwecmass.solid.slice_bspline_surface(patches, z);
catch err
    if ~any(strcmp(err.identifier, {'mwecmass:solid:SectionNotClosed', 'mwecmass:solid:ZOutside'}))
        rethrow(err);
    end
    reason = err.message;
end
end
