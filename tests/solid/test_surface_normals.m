function test_surface_normals
%TEST_SURFACE_NORMALS Normals from mwecmass.solid.surface_normals on the C1 hull.

    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    if exist('OCTAVE_VERSION', 'builtin')
        warning('off', 'Octave:shadowed-function');
        addpath(fullfile(root, 'tests', 'octave_shims'));
    end
    addpath(fullfile(root, 'src'));

    model = mwecmass.geometry.MS2Parser.parse(fullfile(root, 'Input', 'C1.ms2'));
    cache = mwecmass.geometry.precompute_boundary_cache(model, 100);
    [~, grids] = mwecmass.solid.outer_rows(model, cache, []);
    zr = grids.z_range;
    names = model.visible_surfs;
    delta = 1e-3;

    tic;
    [~, orient] = mwecmass.solid.surface_normals(model, cache, names{1}, 0.5, 0.5);
    fprintf('orientation of %d patches decided in %.1f s, signs %s\n', ...
            numel(names), toc, mat2str(orient.sign'));

    zs = zr(1) + [0.03 0.2 0.4 0.6 0.8 0.97]' * (zr(2) - zr(1));
    rows = mwecmass.solid.outer_rows(model, cache, zs, grids);
    ult = cellfun(@(nm) ultimate_source(cache, nm), names, 'UniformOutput', false);

    normals = cell(numel(rows), 1);
    d_unit = 0; d_rev = 0; d_horiz = 0; max_step = 0; max_seam_mirror = 0; max_seam_other = 0;
    n_pts = 0; n_limit = 0; n_horiz = 0;
    for k = 1:numel(rows)
        r = rows(k);
        [n, ~, lim] = mwecmass.solid.surface_normals(model, cache, r.patch, r.u, r.v, orient);
        n_pts = n_pts + r.n;
        n_limit = n_limit + nnz(lim);
        d_unit = max(d_unit, max(abs(sqrt(sum(n.^2, 2)) - 1)));

        nx = [2:r.n, 1];
        jump = atan2(sqrt(sum(cross(n, n(nx, :), 2).^2, 2)), sum(n .* n(nx, :), 2));
        max_step = max(max_step, max(jump));

        nh = n(:, 1:2);
        sel = sqrt(sum(nh.^2, 2)) > 0.1;
        dirs = nh(sel, :) ./ sqrt(sum(nh(sel, :).^2, 2));
        base = r.pts(sel, 1:2);
        check(~any(inside(base + delta * dirs, r.pts)), ...
              'z = %g: horizontal step along +n stays inside the section', zs(k));
        check(all(inside(base - delta * dirs, r.pts)), ...
              'z = %g: horizontal step along -n leaves the section', zs(k));
        n_horiz = n_horiz + nnz(sel);
        normals{k} = n;

        % independent oracle (C1 only): surface1 and its mirrors are revolved about x = 0, y = +-1
        fam1 = strncmp(names(r.patch), 'surface1', 8);
        for i = find(fam1(:))'
            c = [0, sign(r.pts(i, 2))];
            radial = r.pts(i, 1:2) - c;
            if norm(radial) > 0
                d_rev = max(d_rev, abs(n(i, 1) * radial(2) - n(i, 2) * radial(1)) / norm(radial));
            end
        end

        for i = find(r.seam)'
            n1 = n(i, :);
            n2 = mwecmass.solid.surface_normals(model, cache, r.patch2(i), r.u2(i), r.v2(i), orient);
            ang = atan2(norm(cross(n1, n2)), dot(n1, n2));
            if strcmp(ult{r.patch(i)}, ult{r.patch2(i)})
                max_seam_mirror = max(max_seam_mirror, ang);
            else
                max_seam_other = max(max_seam_other, ang);
            end
        end
    end

    out_ok = 0;
    for k = 1:numel(rows)
        r = rows(k);
        n = normals{k};
        for i = round(linspace(1, r.n, 8))
            out_ok = out_ok + step_test(model, cache, grids, r.pts(i, :), n(i, :), delta);
        end
    end
    for p = 1:numel(names)
        for u = [0 1]
            P = model.eval_surface(names{p}, u, 0.5);
            [n, ~, lim] = mwecmass.solid.surface_normals(model, cache, p, u, 0.5, orient);
            out_ok = out_ok + step_test(model, cache, grids, P, n, delta);
            n_limit = n_limit + nnz(lim);
        end
    end

    h = 1e-6;
    d_fd = 0;
    for p = [1 2 3 4]
        for uv = [0.2 0.3; 0.5 0.5; 0.8 0.9]'
            [~, Su, Sv] = model.eval_surface_with_derivs(names{p}, uv(1), uv(2));
            Su_fd = (model.eval_surface(names{p}, uv(1) + h, uv(2)) - ...
                     model.eval_surface(names{p}, uv(1) - h, uv(2))) / (2 * h);
            Sv_fd = (model.eval_surface(names{p}, uv(1), uv(2) + h) - ...
                     model.eval_surface(names{p}, uv(1), uv(2) - h)) / (2 * h);
            d_fd = max([d_fd, norm(Su - Su_fd) / norm(Su), norm(Sv - Sv_fd) / norm(Sv)]);
        end
    end

    fprintf('%d section points, %d normals taken as limits (poles)\n', n_pts, n_limit);
    fprintf('largest | |n| - 1 |                                   %.3e\n', d_unit);
    fprintf('horizontal outward step (%g m) passed at %d points\n', delta, n_horiz);
    fprintf('3-D outward step along +n, section at the stepped height: %d points outside, all passed\n', out_ok);
    fprintf('largest angle between normals of consecutive section points: %.3e rad\n', max_step);
    fprintf('largest angle jump at a seam, mirror-image patches:           %.3e rad\n', max_seam_mirror);
    fprintf('largest angle jump at a seam, surface1 to surface2:           %.3e rad\n', max_seam_other);
    fprintf('revolution-plane oracle (surface1 family), largest deviation: %.3e\n', d_rev);
    fprintf('analytic S_u, S_v vs central differences (h = 1e-6), relative: %.3e\n', d_fd);

    check(d_unit <= 1e-12, 'normals are not unit length to 1e-12');
    check(d_rev <= 1e-12, 'normal leaves the plane through the revolution axis');
    % reflection in a symmetry plane maps the outward normal onto itself at a point of that plane
    check(max_seam_mirror <= 1e-12, 'normals of mirror-image patches differ at their seam');
    check_flat_patches(root);
    fprintf('test_surface_normals passed\n');
end

function check_flat_patches(root)
% Synthetic body revolved about a vertical axis: a flat horizontal deck (z = 1), a
% thin cap (z in [0.98, 1], about 2.3 deg from horizontal), a vertical wall and a flat bottom.
    model = mwecmass.geometry.MS2Parser.parse( ...
                fullfile(root, 'tests', 'solid', 'fixtures', 'capped_cylinder.ms2'));
    cache = mwecmass.geometry.precompute_boundary_cache(model, 60);
    [~, grids] = mwecmass.solid.outer_rows(model, cache, []);
    names = model.visible_surfs;
    [~, orient] = mwecmass.solid.surface_normals(model, cache, names{1}, 0.5, 0.5);
    delta = 1e-2;   % above the 60-point section sagitta (1.4e-3 m) of the unit circle
    n_pts = 0;
    for p = 1:numel(names)
        for u = [0.1 0.5 0.9]
            for v = [0.05 0.25 0.5 0.9]
                P = model.eval_surface(names{p}, u, v);
                n = mwecmass.solid.surface_normals(model, cache, p, u, v, orient);
                check(abs(norm(n) - 1) <= 1e-12, 'flat-patch test: %s normal not unit', names{p});
                step_test(model, cache, grids, P, n, delta);
                radial = P(1:2) / norm(P(1:2));
                switch names{p}
                    case {'deck', 'cap'}
                        check(n(3) > 0 && dot(n(1:2), radial) >= 0, ...
                              'flat-patch test: %s normal does not point up and outward', names{p});
                    case 'bottom'
                        check(n(3) < 0, 'flat-patch test: bottom normal does not point down');
                    case 'wall'
                        check(abs(n(3)) <= 1e-12 && dot(n(1:2), radial) > 0.999999, ...
                              'flat-patch test: wall normal not horizontal outward');
                end
                n_pts = n_pts + 1;
            end
        end
    end
    fprintf('flat/thin-patch body: %d points on 4 patches, outward by the 3-D step test\n', n_pts);
end

function ok = step_test(model, cache, grids, p, n, delta)
    q_out = p + delta * n;
    r_out = mwecmass.solid.outer_rows(model, cache, q_out(3), grids);
    outside = r_out.degenerate || ~inside(q_out(1:2), r_out.pts);
    check(outside, 'step of %g m along +n from (%g, %g, %g) stays inside the hull', ...
          delta, p(1), p(2), p(3));
    q_in = p - delta * n;
    r_in = mwecmass.solid.outer_rows(model, cache, q_in(3), grids);
    if ~r_in.degenerate
        check(inside(q_in(1:2), r_in.pts), ...
              'step of %g m against n from (%g, %g, %g) leaves the hull', delta, p(1), p(2), p(3));
    end
    ok = 1;
end

function in = inside(Q, poly)
% Even-odd crossing test of the rows of Q against the polygon poly (xy columns).
    x = poly(:, 1)'; y = poly(:, 2)';
    xn = [x(2:end), x(1)]; yn = [y(2:end), y(1)];
    in = false(size(Q, 1), 1);
    for i = 1:size(Q, 1)
        straddle = (y > Q(i, 2)) ~= (yn > Q(i, 2));
        xc = x + (Q(i, 2) - y) .* (xn - x) ./ (yn - y);
        in(i) = mod(nnz(straddle & (Q(i, 1) < xc)), 2) == 1;
    end
end

function u = ultimate_source(cache, name)
    u = name;
    for m = 1:numel(cache.mirrors)
        if strcmp(cache.mirrors(m).name, name)
            u = cache.mirrors(m).ultimate_source;
        end
    end
end

function check(cond, varargin)
    if ~cond
        error('test_surface_normals:fail', varargin{:});
    end
end
