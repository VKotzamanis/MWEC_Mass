function data = realised_section_data(realised, z_plan, n_z)
%REALISED_SECTION_DATA Exact sections of the realised (or closest-fail) solid for the Stage-3 figures.
%
%   data = mwecmass.output.figures.realised_section_data(realised, z_plan, n_z)
%
%   realised: results.stage3 (S8; uses hull_name, design, body, props, status, reason, check, vs,
%   mode). z_plan [m, body]: extra heights for plan sections (default none). n_z: heights per module
%   in the elevation (default 41, at least 2). Every section is the mwecmass.solid.body_section of
%   the realised body (the solid that is also written to STEP); nothing is offset, scaled or clipped
%   here. Sections at a module edge are taken 8 ulp of the largest |z| of the module edges: the faces
%   on the two sides of a module edge end at rounding-level offsets from the edge plane, so the
%   section exactly at the plane is not reliable.
%
%   Elevation (y = 0 plane). Each module is sampled at n_z heights between its two edges and at the
%   ballast level. Where two consecutive sampled heights differ in type (solid or hollow, and the
%   number of intervals the void has along y = 0), the height of the change is bisected to adjacent
%   floating-point numbers, and every later change up to the next sampled height is found the same
%   way; the cells next to a change get more heights. A feature that starts and ends inside one cell
%   between two sampled heights of the same type is not resolved; n_z sets that resolution. The
%   ballast top, the closing of the void at its poles and the changes of the void's interval count
%   that are found are exact, and no polygon leaves a gap. The x of every crossing of y = 0 is found
%   on the exact curves (section_y0_crossings).
%   data.polygons(k): module, role ('solid_module' | 'ballast' | 'shell' | 'void'), region
%   (precast 'uhpc' | 'air'; thin shell 'ballast' | 'shell' | 'air'), z_lo, z_hi, xz [n x 2] = [x z]
%   counter-clockwise, body frame. A role 'shell' is the shell around the air: the material around a
%   void, its top and bottom layers included.
%   data.void_outlines(k): modules (the modules the void passes through), xz [n x 2] the outer loop,
%   counter-clockwise, holes (cell of [n x 2] loops, clockwise: material islands inside the air),
%   body frame: one outer loop per connected air region, drawn only where the realised solid has a
%   boundary. The void polygons tile the air exactly: the area of xz minus the areas of the holes is
%   the area of the region's void polygons. Void polygons that abut with overlapping x ranges, at a
%   module edge or at the height inside a module where the number of void intervals changes, are one
%   region: the segment of their common edge is kept only over the x that is air on one side and
%   material on the other. Where t is the same on both sides of a module edge nothing lies there;
%   where it differs the outline keeps only the jog between the two inner sections (the joint_step
%   annulus). Where the void splits in two at a height inside a module, the horizontal boundary is the
%   part of each interval that the other side does not cover; where the two parts join again, the
%   material between them is a hole. Separate regions remain where one side of the edge is solid (the
%   joint cap or the disk above a void) or the void closes inside a module.
%   data.outline: z, x_lo, x_hi, profile (closed [x z] polygon). data.omitted: heights where the
%   kernel gave no usable section (z, module, reason); a section whose outer loop crosses y = 0 at
%   more than two points errors mwecmass:figures:SectionTopology.
%
%   Plan sections. data.plan(k): kind ('module_bottom' | 'module_top' | 'ballast_top' | 'waterline' |
%   'plan'), label, module, z (the height evaluated; 'ballast_top' is just above the ballast level),
%   outer [n x 2] counter-clockwise polyline on the exact curves, inner (same, [] where the section is
%   solid), solid, role (the material of the section, as the polygon roles: 'solid_module', 'ballast'
%   or 'shell'). data.strip_panels(k): plan (index into data.plan), module, end, row, col, the
%   panels of the precast per-module figure: the modules that are not solid modules by design, top
%   views in row 1, bottom views in row 2, one column per module; a view the kernel could not
%   section is absent.
%
%   Frame: body frame [m] throughout; world z = body z + data.vs. Also data.edges, data.z_ballast
%   (body), data.waterline_z (body, = -vs), data.waterline_in_hull, data.solid_modules, data.CG,
%   data.CB (world, from props), data.status, data.reason, data.failed, data.status_lines (cellstr: 'Stage 3 accepted', or
%   'Stage 3 FAILED' followed by the reason split at '; ', one line per failed metric).

if nargin < 2
    z_plan = [];
end
if nargin < 3 || isempty(n_z)
    n_z = 41;
end
if n_z < 2
    error('mwecmass:figures:BadInput', 'realised_section_data: n_z = %g, at least 2 heights are required', n_z);
end
required = {'hull_name', 'design', 'body', 'props', 'status', 'reason', 'check', 'vs', 'mode'};
if ~isstruct(realised) || ~all(isfield(realised, required))
    error('mwecmass:figures:BadRealised', ...
        'realised_section_data: realised must be a struct holding %s (results.stage3)', strjoin(required, ', '));
end
design = realised.design;
body = realised.body;
e = design.edges(:);
N = numel(e) - 1;
z_ballast = design.z_ballast;
solid_modules = design.solid_modules(:)';
margin = 8 * eps(max(abs(e)));

omitted = struct('z', {}, 'module', {}, 'reason', {});
by_module = cell(N, 1);
seams = zeros(0, 3);
polygons = struct('module', {}, 'role', {}, 'region', {}, 'z_lo', {}, 'z_hi', {}, 'xz', {});
for m = 1:N
    [by_module{m}, om] = module_levels(body, e, m, n_z, z_ballast, margin);
    omitted = [omitted, om]; %#ok<AGROW>
    [poly_m, seam_m] = module_polygons(by_module{m}, m, solid_modules, z_ballast, realised.mode);
    polygons = [polygons, poly_m]; %#ok<AGROW>
    seams = [seams; seam_m]; %#ok<AGROW>
end
all_levels = [by_module{:}];
[~, order] = sort([all_levels.z]);
all_levels = all_levels(order);
z_o = [all_levels.z]';
xo = vertcat(all_levels.x_outer);
outline = struct('z', z_o, 'x_lo', xo(:, 1), 'x_hi', xo(:, 2), ...
    'profile', [xo(:, 2), z_o; flipud(xo(:, 1)), flipud(z_o)]);

[plan, plan_omitted] = plan_sections(body, e, z_ballast, realised.vs, z_plan, margin, solid_modules);
omitted = [omitted, plan_omitted];
panels = strip_panel_list(plan, N, solid_modules);
void_outlines = merge_voids(polygons, e, margin, seams);

failed = {};
if isfield(realised.check, 'failed')
    failed = realised.check.failed;
end
if strcmp(realised.status, 'accepted')
    status_lines = {'Stage 3 accepted'};
else
    status_lines = [{'Stage 3 FAILED'}, strsplit(realised.reason, '; ')];
end
data = struct('hull_name', realised.hull_name, 'mode', realised.mode, 'vs', realised.vs, ...
    'edges', e, 'z_range', [e(1) e(end)], 'z_ballast', z_ballast, ...
    'waterline_z', -realised.vs, 'waterline_in_hull', -realised.vs > e(1) && -realised.vs < e(end), ...
    'solid_modules', solid_modules, 'polygons', polygons, 'void_outlines', void_outlines, ...
    'outline', outline, 'plan', plan, ...
    'strip_panels', panels, 'omitted', omitted, 'CG', realised.props.CG_total, 'CB', realised.props.CB, ...
    'status', realised.status, 'reason', realised.reason, 'failed', {failed}, 'status_lines', {status_lines});
end

%% Elevation levels

function [recs, omitted] = module_levels(body, e, m, n_z, z_ballast, margin)
% Ascending level records of module m: uniform heights, the ballast level, the exact heights at which
% the section changes between solid and hollow or in its number of void intervals, and more heights in
% the cells next to those.
z_lo = e(m) + margin;
z_hi = e(m + 1) - margin;
dz = (z_hi - z_lo) / (n_z - 1);
zs = linspace(z_lo, z_hi, n_z)';
if z_ballast > z_lo && z_ballast < z_hi
    zs = sort([zs; z_ballast]);
end
omitted = struct('z', {}, 'module', {}, 'reason', {});
recs = [];
for z = zs'
    [rec, reason] = level_record(body, z, m);
    if isempty(rec)
        omitted(end + 1) = struct('z', z, 'module', m, 'reason', reason); %#ok<AGROW>
    else
        recs = [recs, rec]; %#ok<AGROW>
    end
end
if isempty(recs)
    error('mwecmass:figures:NoSection', 'realised_section_data: module %d has no usable section', m);
end
k = 1;
while k < numel(recs)
    if ~isequal(level_key(recs(k)), level_key(recs(k + 1)))
        [lo_rec, hi_rec] = locate_change(body, m, recs(k), recs(k + 1));
        recs = [recs(1:k), lo_rec, hi_rec, recs(k + 1:end)];
        k = k + 2;
    else
        k = k + 1;
    end
end
recs = recs([true, diff([recs.z]) > 0]);
dense = recs(1);
for k = 1:numel(recs) - 1
    p = recs(k);
    q = recs(k + 1);
    if isequal(level_key(p), level_key(q)) && (p.refined || q.refined)
        n_int = max(3, ceil((q.z - p.z) / dz) - 1);
        zi = linspace(p.z, q.z, n_int + 2);
        for z = zi(2:end - 1)
            [rec, reason] = level_record(body, z, m);
            if isempty(rec)
                omitted(end + 1) = struct('z', z, 'module', m, 'reason', reason); %#ok<AGROW>
            elseif isequal(level_key(rec), level_key(p))
                dense = [dense, rec]; %#ok<AGROW>
            end
        end
    end
    dense = [dense, q]; %#ok<AGROW>
end
recs = dense([true, diff([dense.z]) > 0]);
end

function [lo_rec, hi_rec] = locate_change(body, m, a, b)
% Bisect the height at which the level key changes between levels a and b of module m, to adjacent
% floating-point heights, and return the level records on both sides.
lo = a.z;
hi = b.z;
for it = 1:2000
    mid = lo + (hi - lo) / 2;
    if mid == lo || mid == hi
        break
    end
    key = section_key(body, mid, m);
    if isempty(key) || isequal(key, level_key(a))
        lo = mid;
    else
        hi = mid;
    end
end
lo_rec = a;
if lo ~= a.z
    lo_rec = level_record(body, lo, m);
    if ~isempty(lo_rec)
        lo_rec.refined = true;
    end
end
hi_rec = b;
if hi ~= b.z
    hi_rec = level_record(body, hi, m);
    if ~isempty(hi_rec)
        hi_rec.refined = true;
    end
end
if isempty(lo_rec) || isempty(hi_rec)
    error('mwecmass:figures:NoSection', ...
        'realised_section_data: no section at the change of the section type near z = %.17g (module %d)', lo, m);
end
end

function key = level_key(rec)
% What the polygons of a run share: the solid flag and the number of crossings of the void at y = 0.
key = [rec.solid, numel(rec.x_inner)];
end

function key = section_key(body, z, m)
key = [];
rec = level_record(body, z, m);
if ~isempty(rec)
    key = level_key(rec);
end
end

function [rec, reason] = level_record(body, z, m)
% Level record at height z, or [] when the kernel has no section there or it belongs to another module.
rec = [];
[sec, reason] = try_section(body, z);
if isempty(sec)
    return
end
if sec.module ~= m
    reason = sprintf('the section at z = %.17g belongs to module %d', z, sec.module);
    return
end
xs = mwecmass.output.figures.section_y0_crossings(sec.outer);
if isempty(xs)
    reason = 'the outer section does not reach y = 0';
    return
end
if numel(xs) > 2
    error('mwecmass:figures:SectionTopology', ...
        'realised_section_data: the outer section at z = %g crosses y = 0 at %d points', z, numel(xs));
end
solid = logical(sec.solid) || isempty(sec.inner);
xi = zeros(1, 0);
if ~solid
    xi = mwecmass.output.figures.section_y0_crossings(sec.inner)';
end
rec = struct('z', z, 'solid', solid, 'x_outer', [xs(1), xs(end)], 'x_inner', xi, 'refined', false);
end

function [sec, reason] = try_section(body, z)
sec = [];
reason = '';
try
    sec = mwecmass.solid.body_section(body, z);
catch err
    if ~any(strcmp(err.identifier, {'mwecmass:solid:SectionNotClosed', 'mwecmass:solid:ZOutside'}))
        rethrow(err);
    end
    reason = err.message;
end
end

%% Polygons

function [polygons, seams] = module_polygons(recs, m, solid_modules, z_ballast, mode)
% seams: rows [m, z of the last level of a group, z of the first level of the next group], the heights
% inside module m where the level key changes.
seams = zeros(0, 3);
polygons = struct('module', {}, 'role', {}, 'region', {}, 'z_lo', {}, 'z_hi', {}, 'xz', {});
n = numel(recs);
first = 1;
while first < n
    last = first;
    key = level_key(recs(first));
    while last < n && isequal(level_key(recs(last + 1)), key)
        last = last + 1;
    end
    group = recs(first:last);
    if group(1).solid
        pieces = split_at_ballast(group, z_ballast);
        for q = 1:numel(pieces)
            piece = pieces{q};
            if numel(piece) < 2
                continue
            end
            role = material_role(m, piece(end).z, solid_modules, z_ballast);
            xo = vertcat(piece.x_outer);
            polygons(end + 1) = make_polygon(m, role, mode, [piece.z]', xo(:, 1), xo(:, 2)); %#ok<AGROW>
        end
    elseif numel(group) >= 2
        zz = [group.z]';
        xo = vertcat(group.x_outer);
        xi = vertcat(group.x_inner);
        lefts = [xo(:, 1), xi(:, 2:2:end)];
        rights = [xi(:, 1:2:end), xo(:, 2)];
        for c = 1:size(lefts, 2)
            polygons(end + 1) = make_polygon(m, 'shell', mode, zz, lefts(:, c), rights(:, c)); %#ok<AGROW>
        end
        for c = 1:size(xi, 2) / 2
            polygons(end + 1) = make_polygon(m, 'void', mode, zz, xi(:, 2 * c - 1), xi(:, 2 * c)); %#ok<AGROW>
        end
    end
    if last < n
        seams(end + 1, :) = [m, recs(last).z, recs(last + 1).z]; %#ok<AGROW>
    end
    first = last + 1;
end
end

function role = material_role(m, z_top, solid_modules, z_ballast)
% Role of material that fills a whole section whose top is at z_top in module m.
if any(m == solid_modules)
    role = 'solid_module';
elseif z_top <= z_ballast
    role = 'ballast';
else
    role = 'shell';
end
end

function pieces = split_at_ballast(group, z_ballast)
% A solid run that contains the ballast level strictly inside is two bodies of material.
k = find([group.z] == z_ballast, 1);
if ~isempty(k) && k > 1 && k < numel(group)
    pieces = {group(1:k), group(k:end)};
else
    pieces = {group};
end
end

function p = make_polygon(m, role, mode, z, x_left, x_right)
p = struct('module', m, 'role', role, 'region', region_of(role, mode), 'z_lo', z(1), 'z_hi', z(end), ...
    'xz', [x_right, z; flipud(x_left), flipud(z)]);
end

function region = region_of(role, mode)
if strcmp(role, 'void')
    region = 'air';
elseif strcmp(mode, 'modular_precast')
    region = 'uhpc';
elseif strcmp(role, 'shell')
    region = 'shell';
else
    region = 'ballast';
end
end

function outlines = merge_voids(polygons, e, margin, seams)
% One outer loop, with its holes, per connected air region. Void polygon q lies above void polygon p
% when q is in the next module and the two meet at the module edge (p ends at e - margin, q starts at
% e + margin), or both are in one module and meet at a seam of its level groups; they are one region
% when their x ranges overlap on the common edge. The boundary of a region is the closed chains of the
% polygon sides and of the parts of the common edges that the other side does not cover. A polygon is
% [right side, bottom to top; left side, top to bottom].
outlines = struct('modules', {}, 'xz', {}, 'holes', {});
voids = polygons(strcmp({polygons.role}, 'void'));
nv = numel(voids);
if nv == 0
    return
end
links = zeros(0, 4);
for p = 1:nv
    for q = 1:nv
        if p == q || voids(q).module < voids(p).module
            continue
        end
        m = voids(p).module;
        if voids(q).module == m
            joined = any(seams(:, 1) == m & seams(:, 2) == voids(p).z_hi & seams(:, 3) == voids(q).z_lo);
        else
            joined = voids(q).module == m + 1 && m < numel(e) - 1 && ...
                voids(p).z_hi == e(m + 1) - margin && voids(q).z_lo == e(m + 1) + margin;
        end
        if ~joined
            continue
        end
        top = top_edge(voids(p));
        bottom = bottom_edge(voids(q));
        lo = max(top(1), bottom(1));
        hi = min(top(2), bottom(2));
        if hi > lo
            links(end + 1, :) = [p, q, lo, hi]; %#ok<AGROW>
        end
    end
end
region = 1:nv;
changed = true;
while changed
    changed = false;
    for k = 1:size(links, 1)
        r = min(region(links(k, 1)), region(links(k, 2)));
        if region(links(k, 1)) ~= r || region(links(k, 2)) ~= r
            region(region == region(links(k, 1)) | region == region(links(k, 2))) = r;
            changed = true;
        end
    end
end
found = struct('modules', {}, 'xz', {}, 'holes', {}, 'key', {});
for r = unique(region)
    idx = find(region == r);
    loops = region_loops(voids, idx, links);
    area = cellfun(@signed_area, loops);
    outer = find(area > 0);
    if numel(outer) ~= 1
        error('mwecmass:figures:OutlineTopology', ...
            'realised_section_data: a connected air region has %d counter-clockwise boundary loops', numel(outer));
    end
    xz = loops{outer};
    found(end + 1) = struct('modules', unique([voids(idx).module]), 'xz', xz, ...
        'holes', {loops(area < 0)}, 'key', [min(xz(:, 2)), min(xz(:, 1))]); %#ok<AGROW>
end
[~, order] = sortrows(vertcat(found.key));
outlines = rmfield(found(order), 'key');
end

function a = signed_area(xz)
a = 0.5 * sum(xz(:, 1) .* xz([2:end, 1], 2) - xz([2:end, 1], 1) .* xz(:, 2));
end

function x = top_edge(v)
n = size(v.xz, 1) / 2;
x = [v.xz(n + 1, 1), v.xz(n, 1)];
end

function x = bottom_edge(v)
n = size(v.xz, 1) / 2;
x = [v.xz(2 * n, 1), v.xz(1, 1)];
end

function loops = region_loops(voids, idx, links)
% Closed boundary chains of the void polygons voids(idx). Edges are matched by (x, z), where the z of
% the two sides of a shared edge is taken as the z of the lower one.
S = zeros(0, 2);    % key of the start of each edge
T = zeros(0, 2);    % key of the end
P = zeros(0, 4);    % actual [x z] of the start and of the end
for i = idx
    v = voids(i);
    n = size(v.xz, 1) / 2;
    kz = v.xz(:, 2);
    below_idx = links(links(:, 2) == i, 1);
    if ~isempty(below_idx)
        kz([1, 2 * n]) = voids(below_idx(1)).z_hi;
    end
    for j = [1:n - 1, n + 1:2 * n - 1]
        S(end + 1, :) = [v.xz(j, 1), kz(j)]; %#ok<AGROW>
        T(end + 1, :) = [v.xz(j + 1, 1), kz(j + 1)]; %#ok<AGROW>
        P(end + 1, :) = [v.xz(j, :), v.xz(j + 1, :)]; %#ok<AGROW>
    end
    above = links(links(:, 1) == i, :);
    for q = uncovered(top_edge(v), above(:, 3:4))'
        S(end + 1, :) = [q(2), kz(n)]; %#ok<AGROW>
        T(end + 1, :) = [q(1), kz(n)]; %#ok<AGROW>
        P(end + 1, :) = [q(2), v.xz(n, 2), q(1), v.xz(n, 2)]; %#ok<AGROW>
    end
    below = links(links(:, 2) == i, :);
    for q = uncovered(bottom_edge(v), below(:, 3:4))'
        S(end + 1, :) = [q(1), kz(1)]; %#ok<AGROW>
        T(end + 1, :) = [q(2), kz(1)]; %#ok<AGROW>
        P(end + 1, :) = [q(1), v.xz(1, 2), q(2), v.xz(1, 2)]; %#ok<AGROW>
    end
end
used = false(size(S, 1), 1);
loops = {};
while any(~used)
    first = find(~used, 1);
    cur = first;
    pts = zeros(0, 2);
    while true
        used(cur) = true;
        pts(end + 1, :) = P(cur, 1:2); %#ok<AGROW>
        nxt = find(~used & S(:, 1) == T(cur, 1) & S(:, 2) == T(cur, 2), 1);
        if isempty(nxt)
            if ~isequal(T(cur, :), S(first, :))
                error('mwecmass:figures:OutlineOpen', ...
                    'realised_section_data: the void outline does not close at (x, z) = (%.17g, %.17g)', T(cur, 1), T(cur, 2));
            end
            nxt_start = P(first, 1:2);
        else
            nxt_start = P(nxt, 1:2);
        end
        if ~isequal(P(cur, 3:4), nxt_start)
            pts(end + 1, :) = P(cur, 3:4); %#ok<AGROW>
        end
        if isempty(nxt)
            break
        end
        cur = nxt;
    end
    loops{end + 1} = pts; %#ok<AGROW>
end
end

function free = uncovered(span, covers)
% Parts [lo hi] of the interval span = [a b] that no row of covers = [lo hi] covers, as rows.
free = zeros(0, 2);
pos = span(1);
if ~isempty(covers)
    covers = sortrows(covers);
end
for k = 1:size(covers, 1)
    if covers(k, 1) > pos
        free(end + 1, :) = [pos, covers(k, 1)]; %#ok<AGROW>
    end
    pos = max(pos, covers(k, 2));
end
if span(2) > pos
    free(end + 1, :) = [pos, span(2)];
end
end

%% Plan sections

function [plan, omitted] = plan_sections(body, e, z_ballast, vs, z_plan, margin, solid_modules)
N = numel(e) - 1;
plan = struct('kind', {}, 'label', {}, 'module', {}, 'z', {}, 'outer', {}, 'inner', {}, 'solid', {}, 'role', {});
omitted = struct('z', {}, 'module', {}, 'reason', {});
ctx = {body, e, solid_modules, z_ballast};
for m = 1:N
    [plan, omitted] = add_plan(plan, omitted, ctx, 'module_bottom', sprintf('module %d, bottom', m), e(m) + margin);
    [plan, omitted] = add_plan(plan, omitted, ctx, 'module_top', sprintf('module %d, top', m), e(m + 1) - margin);
end
if z_ballast > e(1) && z_ballast < e(end)
    [plan, omitted] = add_plan(plan, omitted, ctx, 'ballast_top', 'ballast level', z_ballast + margin);
end
z_wl = -vs;
if z_wl > e(1) && z_wl < e(end)
    [plan, omitted] = add_plan(plan, omitted, ctx, 'waterline', 'waterline', z_wl);
end
for z = z_plan(:)'
    if z < e(1) || z > e(end)
        error('mwecmass:figures:ZOutside', 'realised_section_data: z_plan = %g is outside the hull [%g, %g]', ...
            z, e(1), e(end));
    end
    [plan, omitted] = add_plan(plan, omitted, ctx, 'plan', sprintf('z = %.3g m (body)', z), z);
end
end

function [plan, omitted] = add_plan(plan, omitted, ctx, kind, label, z)
[body, e, solid_modules, z_ballast] = ctx{:};
m = find(e(1:end - 1) <= z & z < e(2:end), 1);
if isempty(m)
    m = numel(e) - 1;
end
[sec, reason] = try_section(body, z);
if ~isempty(sec) && sec.module ~= m
    reason = sprintf('the section at z = %.17g belongs to module %d', z, sec.module);
    sec = [];
end
if isempty(sec)
    omitted(end + 1) = struct('z', z, 'module', m, 'reason', reason);
    return
end
inner = [];
solid = logical(sec.solid) || isempty(sec.inner);
if ~solid
    inner = loop_polyline(sec.inner, 25);
end
plan(end + 1) = struct('kind', kind, 'label', label, 'module', m, 'z', z, ...
    'outer', loop_polyline(sec.outer, 25), 'inner', inner, 'solid', solid, ...
    'role', material_role(m, z, solid_modules, z_ballast));
end

function pts = loop_polyline(loop, n)
% Counter-clockwise points on the exact curves of the loop pieces, n per piece, no repeated point.
pts = zeros(0, 2);
for k = 1:numel(loop.pieces)
    c = loop.pieces(k).curve;
    s = linspace(c.knots(1), c.knots(end), n + 1)';
    if loop.pieces(k).dir < 0
        s = flipud(s);
    end
    p = mwecmass.solid.eval_bspline_curve(c, s);
    pts = [pts; p(1:end - 1, 1:2)]; %#ok<AGROW>
end
end

function panels = strip_panel_list(plan, N, solid_modules)
panels = struct('plan', {}, 'module', {}, 'end', {}, 'row', {}, 'col', {});
shown = setdiff(1:N, solid_modules);
for c = 1:numel(shown)
    m = shown(c);
    for q = {'top', 1; 'bottom', 2}'
        k = find(strcmp({plan.kind}, ['module_' q{1}]) & [plan.module] == m, 1);
        if ~isempty(k)
            panels(end + 1) = struct('plan', k, 'module', m, 'end', q{1}, 'row', q{2}, 'col', c); %#ok<AGROW>
        end
    end
end
end
