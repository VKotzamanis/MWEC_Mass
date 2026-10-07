function validate_brep(brep)
%VALIDATE_BREP  Check the structural consistency of a B-rep struct (see write_step).
%
%   Errors with identifier mwecmass:step:invalid on the first inconsistency: index ranges,
%   spline data, loop closure, curve end points against vertices, planar faces against their
%   plane, and for every shell the pairing of edge uses (closed shell: each edge once in each
%   direction; open shell: at most once in each direction).

unc = mwecmass.output.step.declared_uncertainty(brep);

need(brep, {'vertices', 'curves', 'edges', 'surfaces', 'faces', 'bodies'}, 'brep');
V = brep.vertices;
if ~isnumeric(V) || (size(V, 2) ~= 3 && ~isempty(V)) || ~all(isfinite(V(:)))
    fail('vertices must be a finite [nv x 3] array');
end
if ~iscell(brep.surfaces)
    fail('surfaces must be a cell array');
end

for k = 1:numel(brep.curves)
    check_curve(brep.curves(k), k);
end
for k = 1:numel(brep.surfaces)
    check_surface(brep.surfaces{k}, k);
end

ne = numel(brep.edges);
for k = 1:ne
    ed = brep.edges(k);
    if numel(ed.vertices) ~= 2 || any(~isindex(ed.vertices, size(V, 1)))
        fail(sprintf('edge %d: vertices must be two valid vertex indices', k));
    end
    if ~isindex(ed.curve, numel(brep.curves))
        fail(sprintf('edge %d: invalid curve index', k));
    end
    c = brep.curves(ed.curve);
    if norm(c.ctrl(1, :) - V(ed.vertices(1), :)) > unc || ...
            norm(c.ctrl(end, :) - V(ed.vertices(2), :)) > unc
        fail(sprintf('edge %d: curve end points differ from its vertices by more than the stated uncertainty', k));
    end
end

for k = 1:numel(brep.faces)
    f = brep.faces(k);
    if ~isindex(f.surface, numel(brep.surfaces))
        fail(sprintf('face %d: invalid surface index', k));
    end
    if ~iscell(f.loops) || isempty(f.loops)
        fail(sprintf('face %d: loops must be a non-empty cell array', k));
    end
    onplane = strcmp(brep.surfaces{f.surface}.type, 'plane');
    if onplane
        s = brep.surfaces{f.surface};
        n = s.normal(:)' / norm(s.normal);
    end
    for j = 1:numel(f.loops)
        L = f.loops{j}(:)';
        if isempty(L) || any(L == 0) || any(L ~= fix(L)) || any(abs(L) > ne)
            fail(sprintf('face %d loop %d: invalid signed edge indices', k, j));
        end
        [a, b] = loop_ends(brep.edges, L);
        if any(b ~= a([2:end 1]))
            fail(sprintf('face %d loop %d: edges do not form a closed chain', k, j));
        end
        if onplane
            d = (V(a, :) - s.origin(:)') * n';
            if any(abs(d) > unc)
                fail(sprintf('face %d loop %d: vertices lie off the plane by more than the stated uncertainty', k, j));
            end
        end
    end
end

names = cell(1, numel(brep.bodies));
for k = 1:numel(brep.bodies)
    bd = brep.bodies(k);
    if ~ischar(bd.name) || isempty(bd.name) || any(bd.name < 32 | bd.name > 126)
        fail(sprintf('body %d: name must be a non-empty printable ASCII string', k));
    end
    names{k} = bd.name;
    if ~any(strcmp(bd.kind, {'solid', 'sheet'}))
        fail(sprintf('body %s: kind must be solid or sheet', bd.name));
    end
    if ~iscell(bd.shells) || isempty(bd.shells)
        fail(sprintf('body %s: shells must be a non-empty cell array', bd.name));
    end
    for j = 1:numel(bd.shells)
        S = bd.shells{j}(:)';
        if isempty(S) || any(S == 0) || any(S ~= fix(S)) || any(abs(S) > numel(brep.faces))
            fail(sprintf('body %s shell %d: invalid signed face indices', bd.name, j));
        end
        [pos, neg] = edge_uses(brep, S);
        if strcmp(bd.kind, 'solid')
            used = pos + neg > 0;
            if any(pos(used) ~= 1) || any(neg(used) ~= 1)
                fail(sprintf('body %s shell %d: not closed, every edge needs one use in each direction', bd.name, j));
            end
        elseif any(pos > 1) || any(neg > 1)
            fail(sprintf('body %s shell %d: an edge is used twice in the same direction', bd.name, j));
        end
    end
end
if numel(unique(names)) ~= numel(names)
    fail('body names must be unique (they name the layers)');
end

end

function fail(msg)
    error('mwecmass:step:invalid', 'write_step: %s', msg);
end

function check_curve(c, k)
    need(c, {'degree', 'ctrl', 'knots', 'weights'}, sprintf('curve %d', k));
    n = size(c.ctrl, 1);
    if size(c.ctrl, 2) ~= 3 || ~all(isfinite(c.ctrl(:)))
        fail(sprintf('curve %d: ctrl must be a finite [n x 3] array', k));
    end
    check_knots(c.knots, c.degree, n, sprintf('curve %d', k));
    check_weights(c.weights(:), [n 1], sprintf('curve %d', k));
end

function check_surface(s, k)
    need(s, {'type'}, sprintf('surface %d', k));
    if strcmp(s.type, 'plane')
        need(s, {'origin', 'normal'}, sprintf('surface %d', k));
        if numel(s.origin) ~= 3 || numel(s.normal) ~= 3 || norm(s.normal) == 0 || ...
                ~all(isfinite([s.origin(:); s.normal(:)]))
            fail(sprintf('surface %d: plane needs a finite origin and a non-zero normal', k));
        end
    elseif strcmp(s.type, 'bspline')
        need(s, {'degree', 'ctrl', 'knots', 'weights'}, sprintf('surface %d', k));
        sz = size(s.ctrl);
        if numel(sz) ~= 3 || sz(3) ~= 3 || ~all(isfinite(s.ctrl(:)))
            fail(sprintf('surface %d: ctrl must be a finite [nu x nv x 3] array', k));
        end
        if ~iscell(s.knots) || numel(s.knots) ~= 2 || numel(s.degree) ~= 2
            fail(sprintf('surface %d: knots must be {ku, kv} and degree [du dv]', k));
        end
        check_knots(s.knots{1}, s.degree(1), sz(1), sprintf('surface %d u', k));
        check_knots(s.knots{2}, s.degree(2), sz(2), sprintf('surface %d v', k));
        check_weights(s.weights, sz(1:2), sprintf('surface %d', k));
    else
        fail(sprintf('surface %d: type must be plane or bspline', k));
    end
end

function check_knots(t, p, n, label)
    if ~isscalar(p) || p < 1 || p ~= fix(p) || n < p + 1 || numel(t) ~= n + p + 1 || ...
            ~all(isfinite(t(:))) || any(diff(t(:)) < 0) || t(end) <= t(1)
        fail(sprintf('%s: inconsistent degree, control point count and knot vector', label));
    end
    if any(t(1:p + 1) ~= t(1)) || any(t(end - p:end) ~= t(end))
        fail(sprintf('%s: knot vector must be clamped (end knots of multiplicity degree+1)', label));
    end
    if any(histc_mult(t) > p + 1) || any(histc_mult(t(2:end - 1)) > p)
        fail(sprintf('%s: interior knot multiplicity exceeds the degree', label));
    end
end

function check_weights(w, sz, label)
    if ~isempty(w) && (~isequal(size(w), sz) || ~all(isfinite(w(:))) || any(w(:) <= 0))
        fail(sprintf('%s: weights must be empty or positive with the size of the control net', label));
    end
end

function need(s, names, label)
for i = 1:numel(names)
    if ~isfield(s, names{i})
        error('mwecmass:step:invalid', 'write_step: %s lacks field %s', label, names{i});
    end
end
end

function ok = isindex(v, n)
ok = isnumeric(v) & v == fix(v) & v >= 1 & v <= n;
end

function m = histc_mult(t)
[~, ~, j] = unique(t(:));
m = accumarray(j, 1);
end

function [a, b] = loop_ends(edges, L)
v = vertcat(edges(abs(L)).vertices);
fwd = L(:) > 0;
a = v(:, 1);
b = v(:, 2);
a(~fwd) = v(~fwd, 2);
b(~fwd) = v(~fwd, 1);
a = a(:)';
b = b(:)';
end

function [pos, neg] = edge_uses(brep, S)
ne = numel(brep.edges);
pos = zeros(1, ne);
neg = zeros(1, ne);
for i = 1:numel(S)
    f = brep.faces(abs(S(i)));
    flip = sign(S(i));
    for j = 1:numel(f.loops)
        for e = f.loops{j}(:)' * flip
            if e > 0
                pos(e) = pos(e) + 1;
            else
                neg(-e) = neg(-e) + 1;
            end
        end
    end
end
end
