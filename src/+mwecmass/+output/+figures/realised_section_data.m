function data = realised_section_data(realised, z_plan, n_z)
%REALISED_SECTION_DATA Exact sections of the realised (or closest-fail) solid for the Stage-3 figures.
%
%   data = mwecmass.output.figures.realised_section_data(realised, z_plan, n_z)
%
%   realised: results.stage3 (S8; uses hull_name, design, body, props, status, reason, check, vs,
%   mode). z_plan [m, body]: extra heights for plan sections (default none). n_z: heights per module
%   in the elevation (default 41, at least 2). Every section is the mwecmass.solid.body_section of
%   the realised body (the solid that is also written to STEP); nothing is offset, scaled or clipped
%   here. Sections at a module edge are taken 8 ulp of the largest |z| inside the module: the faces
%   on the two sides of a module edge end at rounding-level offsets from the edge plane, so the
%   section exactly at the plane is not reliable.
%
%   Elevation (y = 0 plane). Each module is sampled between its two edges and at the ballast level;
%   wherever the section changes between solid and hollow between two heights, the height of the
%   change is bisected to adjacent floating-point numbers, so the ballast top and the closing of the
%   void at its poles are exact, and the cell above such a change gets more heights. The x of every
%   crossing of y = 0 is found on the exact curves (section_y0_crossings).
%   data.polygons(k): module, role ('solid_module' | 'ballast' | 'wall' | 'void'), region
%   (precast 'uhpc' | 'air'; thin shell 'ballast' | 'shell' | 'air'), z_lo, z_hi, xz [n x 2] = [x z]
%   counter-clockwise, body frame. A role 'wall' is the material around a void, or the cap above it.
%   data.void_outlines(k): modules (the modules the void passes through), xz [n x 2] counter-clockwise,
%   body frame: the boundary of the air region as the realised solid has it. The void polygons of
%   consecutive modules are one outline where both hold air at the module edge and the two modules
%   have the same t (no face between them: no joint cap, no joint_step), so no segment lies at that
%   edge; where t differs, or one side is solid, the outlines stay separate and each keeps its face
%   at the edge. data.outline: z, x_lo, x_hi, profile (closed [x z] polygon). data.omitted: heights where the
%   kernel gave no usable section (z, module, reason); a section whose outer loop crosses y = 0 at
%   more than two points errors mwecmass:figures:SectionTopology.
%
%   Plan sections. data.plan(k): kind ('module_bottom' | 'module_top' | 'ballast_top' | 'waterline' |
%   'plan'), label, module, z (the height evaluated; 'ballast_top' is just above the ballast level),
%   outer [n x 2] counter-clockwise polyline on the exact curves, inner (same, [] where the section is
%   solid), solid, role (the material of the section, as the polygon roles: 'solid_module', 'ballast'
%   or 'wall'). data.strip_panels(k): plan (index into data.plan), module, end, row, col, the
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
polygons = struct('module', {}, 'role', {}, 'region', {}, 'z_lo', {}, 'z_hi', {}, 'xz', {});
for m = 1:N
    [by_module{m}, om] = module_levels(body, e, m, n_z, z_ballast, margin);
    omitted = [omitted, om]; %#ok<AGROW>
    polygons = [polygons, module_polygons(by_module{m}, m, solid_modules, z_ballast, realised.mode)]; %#ok<AGROW>
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
void_outlines = merge_voids(polygons, e, margin, design.t);

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
% the section changes between solid and hollow, and more heights in the cells those leave short.
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
    if recs(k).solid ~= recs(k + 1).solid
        [lo_rec, hi_rec] = locate_change(body, m, recs(k), recs(k + 1));
        recs = [recs(1:k), lo_rec, hi_rec, recs(k + 1:end)];
        k = k + 3;
    else
        k = k + 1;
    end
end
recs = recs([true, diff([recs.z]) > 0]);
dense = recs(1);
for k = 1:numel(recs) - 1
    p = recs(k);
    q = recs(k + 1);
    if p.solid == q.solid && (p.refined || q.refined)
        n_int = max(3, ceil((q.z - p.z) / dz) - 1);
        zi = linspace(p.z, q.z, n_int + 2);
        for z = zi(2:end - 1)
            [rec, reason] = level_record(body, z, m);
            if isempty(rec)
                omitted(end + 1) = struct('z', z, 'module', m, 'reason', reason); %#ok<AGROW>
            elseif rec.solid == p.solid
                dense = [dense, rec]; %#ok<AGROW>
            end
        end
    end
    dense = [dense, q]; %#ok<AGROW>
end
recs = dense([true, diff([dense.z]) > 0]);
end

function [lo_rec, hi_rec] = locate_change(body, m, a, b)
% Bisect the height at which the solid flag changes between levels a and b of module m, to adjacent
% floating-point heights, and return the level records on both sides.
lo = a.z;
hi = b.z;
for it = 1:2000
    mid = lo + (hi - lo) / 2;
    if mid == lo || mid == hi
        break
    end
    flag = section_is_solid(body, mid, m);
    if isempty(flag) || flag == a.solid
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
        'realised_section_data: no section at the change of the solid flag near z = %.17g (module %d)', lo, m);
end
end

function flag = section_is_solid(body, z, m)
flag = [];
sec = try_section(body, z);
if ~isempty(sec) && sec.module == m
    flag = logical(sec.solid) || isempty(sec.inner);
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

function polygons = module_polygons(recs, m, solid_modules, z_ballast, mode)
polygons = struct('module', {}, 'role', {}, 'region', {}, 'z_lo', {}, 'z_hi', {}, 'xz', {});
n = numel(recs);
first = 1;
while first < n
    last = first;
    key = [recs(first).solid, numel(recs(first).x_inner)];
    while last < n && isequal([recs(last + 1).solid, numel(recs(last + 1).x_inner)], key)
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
            polygons(end + 1) = make_polygon(m, 'wall', mode, zz, lefts(:, c), rights(:, c)); %#ok<AGROW>
        end
        for c = 1:size(xi, 2) / 2
            polygons(end + 1) = make_polygon(m, 'void', mode, zz, xi(:, 2 * c - 1), xi(:, 2 * c)); %#ok<AGROW>
        end
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
    role = 'wall';
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
elseif strcmp(role, 'wall')
    region = 'shell';
else
    region = 'ballast';
end
end

function outlines = merge_voids(polygons, e, margin, t)
% One outline per connected air region: the void polygon of module m that reaches the top of the
% module joins the void polygon of module m + 1 that starts at the bottom of m + 1 when the two
% x ranges overlap at the edge and t(m) == t(m + 1): the solid then has the same inner surface on both
% sides of the edge. A polygon is [right side, bottom to top; left side, top to bottom].
outlines = struct('modules', {}, 'xz', {});
voids = polygons(strcmp({polygons.role}, 'void'));
used = false(1, numel(voids));
[~, order] = sort([voids.z_lo]);
for a = order
    if used(a)
        continue
    end
    used(a) = true;
    chain = voids(a);
    right = chain.xz(1:size(chain.xz, 1) / 2, :);
    left = chain.xz(size(chain.xz, 1) / 2 + 1:end, :);
    modules = chain.module;
    while true
        m = modules(end);
        top = max(right(:, 2));
        nxt = 0;
        if m < numel(e) - 1 && top == e(m + 1) - margin && t(m) == t(m + 1)
            for b = find(~used)
                q = voids(b);
                if q.module == m + 1 && q.z_lo == e(m + 1) + margin && ...
                        max(left(1, 1), q.xz(end, 1)) < min(right(end, 1), q.xz(1, 1))
                    nxt = b;
                    break
                end
            end
        end
        if nxt == 0
            break
        end
        used(nxt) = true;
        n = size(voids(nxt).xz, 1) / 2;
        right = [right; voids(nxt).xz(1:n, :)]; %#ok<AGROW>
        left = [voids(nxt).xz(n + 1:end, :); left]; %#ok<AGROW>
        modules(end + 1) = voids(nxt).module; %#ok<AGROW>
    end
    outlines(end + 1) = struct('modules', modules, 'xz', [right; left]); %#ok<AGROW>
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
