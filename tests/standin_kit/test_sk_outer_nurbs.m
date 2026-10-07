function test_sk_outer_nurbs()
%TEST_SK_OUTER_NURBS  Stand-in F1 and F3 on the fixture decks: S1/S1b fields, exactness, seams, poles, orientation.
%   Reads geo.analytic (stand-in marker of F1): J1, which merges the real F1 and F3, deletes or
%   rewrites this test.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
s1_fields = {'name', 'source', 'type', 'flips', 'surf', 'outward', 'exact', 'z_of_u', 'u_range', ...
    'z_range', 'offset_kind', 'pole', 'c0_u', 'c0_v', 'seam_u0', 'seam_u1', 'seam_v0', 'seam_v1', ...
    'visible', 'fit', 'swap_uv'};
for name = {'cylinder', 'box'}
    evalc('model = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''tests'', ''standins'', ''fixtures'', [name{1} ''.ms2'']));');
    geo = mwecmass.solid.outer_nurbs(model);
    check(isequal(sort(fieldnames(geo))', sort({'hull_name', 'outer', 'z_range', 'analytic', 'flat'})), 'S1b fields');
    check(isequal(sort(fieldnames(geo.outer))', sort(s1_fields)), 'S1 fields');
    check(strcmp(geo.hull_name, name{1}) && isequal({geo.outer.name}, model.visible_surfs), 'names in visible order');
    check(isempty(geo.flat), 'flat is empty');
    check(isequal([geo.outer.visible], 1:numel(geo.outer)) && all(cellfun(@isempty, {geo.outer.fit})) && ...
        ~any([geo.outer.swap_uv]), 'visible = k, fit empty, swap_uv false');
    % F1 with cache and opts (the general-path inputs) ignores them
    check(isequaln(geo, mwecmass.solid.outer_nurbs(model, [], struct())) && ...
        isequaln(geo, mwecmass.solid.outer_nurbs(model, struct('unused', 1), struct('t_min', 0.1, 'max_passes', 3, 'force_general', true))), ...
        '%s: outer_nurbs(model, cache, opts) differs from outer_nurbs(model)', name{1});
    fx = geo.analytic;
    hl = sti_closed_form('hull', fx);
    worst_eval = 0;
    for k = 1:numel(geo.outer)
        p = geo.outer(k);
        s = p.surf;
        [nu, nv, ~] = size(s.ctrl);
        check(all(all(s.ctrl(:, :, 3) == s.ctrl(:, 1, 3))) && p.z_of_u, '%s: every control row has one z', p.name);
        check(isequal(p.z_range, [s.ctrl(1, 1, 3), s.ctrl(end, 1, 3)]), '%s: z_range', p.name);
        check(p.exact, '%s: exact', p.name);
        rows = {squeeze(s.ctrl(1, :, :)), squeeze(s.ctrl(end, :, :))};
        for e = 1:2
            check(p.pole(e) == all(all(rows{e} == rows{e}(1, :))), '%s: pole flag %d', p.name, e);
        end
        for c = p.c0_u
            check(sum(s.knots{1} == c) == s.degree(1), '%s: c0_u knot multiplicity', p.name);
        end
        % one parameter point per span: S_u x S_v against the direction from the hull centroid
        u = 0.5 * (s.knots{1}(1:end - 1) + s.knots{1}(2:end));
        u = u(diff(s.knots{1}) > 0)';
        [S, Su, Sv] = mwecmass.solid.eval_bspline_surface(s, u, 0.5 * ones(size(u)));
        n = cross(Su, Sv, 2);
        dirn = sign(sum(n .* (S - hl.centroid), 2));
        check(all(dirn == (2 * p.outward - 1)), '%s: outward flag', p.name);
        % exactness against MS2Parser: same point on the deck's surface (the oracle parameter is
        % found from the NURBS point: the angle about the axis for the RevSurf, the same (u, v)
        % for the bilinear RuledSurf)
        [uu, vv] = meshgrid(linspace(0, 1, 7), linspace(0, 1, 5));
        P = mwecmass.solid.eval_bspline_surface(s, uu(:), vv(:));
        for j = 1:numel(uu)
            if strcmp(fx.kind, 'cylinder')
                phi = atan2(abs(P(j, 2)), abs(P(j, 1)));
                Q = model.eval_surface(p.name, uu(j), phi / (pi / 2));
            else
                Q = model.eval_surface(p.name, uu(j), vv(j));
            end
            worst_eval = max(worst_eval, norm(P(j, :) - Q));
        end
    end
    fprintf('%s: %d patches, largest |NURBS - MS2Parser| %.3e m\n', name{1}, numel(geo.outer), worst_eval);
    % both are the same exact map evaluated with a few dozen rounded operations on coordinates <= 3 m
    check(worst_eval <= 64 * eps * 3, '%s: NURBS differs from MS2Parser beyond rounding', name{1});

    % seams: the boundary rows of the two patches are bitwise equal (either direction)
    fields = {'seam_v0', 'seam_u1', 'seam_v1', 'seam_u0'};
    n_seam = 0;
    for k = 1:numel(geo.outer)
        for b = 1:4
            nb = geo.outer(k).(fields{b});
            if isempty(nb)
                check(b == 2 && geo.outer(k).pole(2) || b == 4 && geo.outer(k).pole(1), ...
                    '%s: empty seam %d that is not a pole', geo.outer(k).name, b);
                continue
            end
            [ra, wa] = boundary(geo.outer(k).surf, b);
            [rb, wb] = boundary(geo.outer(nb(1)).surf, nb(2));
            back = geo.outer(nb(1)).(fields{nb(2)});
            check(isequal(back, [k b]), '%s: seam %d is not mutual', geo.outer(k).name, b);
            same = isequal(ra, rb) && isequal(wa, wb);
            rev = isequal(ra, flipud(rb)) && isequal(wa, flipud(wb));
            check(same || rev, '%s: seam %d rows differ', geo.outer(k).name, b);
            n_seam = n_seam + 1;
        end
    end
    fprintf('%s: %d seam boundaries bitwise equal\n', name{1}, n_seam);
end

evalc('c1 = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''Input'', ''C1.ms2''));');
try
    mwecmass.solid.outer_nurbs(c1);
    error('test_sk_outer_nurbs: C1 did not error');
catch err
    check(strcmp(err.identifier, 'mwecmass:standin:NotAnalytic'), 'C1: %s', err.identifier);
end
fprintf('C1 deck: outer_nurbs stand-in errors mwecmass:standin:NotAnalytic\n');
end

function [r, w] = boundary(s, b)
% boundary 1 = v0, 2 = u1, 3 = v1, 4 = u0 (MS2Parser EdgeSnake numbering)
W = s.weights;
if isempty(W)
    W = ones(size(s.ctrl, 1), size(s.ctrl, 2));
end
switch b
    case 1, r = squeeze(s.ctrl(:, 1, :)); w = W(:, 1);
    case 2, r = squeeze(s.ctrl(end, :, :)); w = W(end, :)';
    case 3, r = squeeze(s.ctrl(:, end, :)); w = W(:, end);
    case 4, r = squeeze(s.ctrl(1, :, :)); w = W(1, :)';
end
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'standins'), '-end');
addpath(fullfile(root, 'tests', 'standins', 'fixtures'), '-end');
end

function check(cond, varargin)
if ~cond
    error('test_sk_outer_nurbs:fail', varargin{:});
end
end
