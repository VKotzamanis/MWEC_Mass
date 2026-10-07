function test_outer_rows
%TEST_OUTER_ROWS Sections from mwecmass.solid.outer_rows on the C1 hull.

    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    if exist('OCTAVE_VERSION', 'builtin')
        warning('off', 'Octave:shadowed-function');
        addpath(fullfile(root, 'tests', 'octave_shims'));
    end
    addpath(fullfile(root, 'src'));

    model = mwecmass.geometry.MS2Parser.parse(fullfile(root, 'Input', 'C1.ms2'));
    n_u = 100;
    cache = mwecmass.geometry.precompute_boundary_cache(model, n_u);
    [~, grids] = mwecmass.solid.outer_rows(model, cache, []);
    zr = grids.z_range;
    fprintf('hull z range [%.6f, %.6f] m\n', zr);

    % strip edges of C1 in the body frame (AGENTS.md 2.5), interior of the z range
    z_edges = [-2.706, -1.619, -0.531, 0.556];
    zs = [zr(1) + (1:25)' / 26 * (zr(2) - zr(1)); z_edges'];
    rows = mwecmass.solid.outer_rows(model, cache, zs, grids);
    check(numel(rows) == numel(zs), 'one row per height');
    fprintf('%d heights\n', numel(zs));

    ends = mwecmass.solid.outer_rows(model, cache, zr(:), grids);
    check(all([ends.degenerate]) && all([ends.n] == 0), 'range ends are degenerate sections');

    worst = struct('area', 0, 'ixx', 0, 'iyy', 0, 'area_zk', 0, 'zres_ref', 0, ...
                   'sym', 0, 'gap', 0, 'mind', Inf, 'zrow', 0, 'shoe', 0);
    for k = 1:numel(zs)
        r = rows(k);
        P = r.pts;
        N = r.n;
        check(~r.degenerate && N >= 3, 'z = %g: section is empty', zs(k));
        nx = [2:N, 1];

        s_area = 0.5 * sum(P(:, 1) .* P(nx, 2) - P(nx, 1) .* P(:, 2));
        check(s_area > 0, 'z = %g: section is not counter-clockwise', zs(k));
        worst.shoe = max(worst.shoe, abs(s_area - r.area) / r.area);

        seg = sqrt(sum((P(nx, :) - P).^2, 2));
        check(seg(N) <= max(seg(1:N-1)), 'z = %g: closing chord longer than every other', zs(k));
        D = sqrt((P(:, 1) - P(:, 1)').^2 + (P(:, 2) - P(:, 2)').^2);
        D(1:N+1:end) = Inf;
        mind = min(D(:));
        check(mind > r.seam_gap, 'z = %g: repeated point (min distance %g, seam gap %g)', ...
              zs(k), mind, r.seam_gap);
        worst.mind = min(worst.mind, mind);
        worst.gap = max(worst.gap, r.seam_gap);
        check(any(r.seam) && all(r.patch2(r.seam) > 0 & r.patch2(r.seam) ~= r.patch(r.seam)) && ...
              all(r.patch2(~r.seam) == 0), 'z = %g: seam bookkeeping', zs(k));

        worst.zrow = max(worst.zrow, max(abs(P(:, 3) - zs(k))));

        for flip = [-1 1 1; 1 -1 1]'
            M = P .* flip';
            dm = sqrt((M(:, 1) - P(:, 1)').^2 + (M(:, 2) - P(:, 2)').^2);
            worst.sym = max(worst.sym, max(min(dm, [], 2)));
            check(max(min(dm, [], 2)) <= 1e-12, ...
                  'z = %g: mirror symmetry differs by %g m', zs(k), max(min(dm, [], 2)));
        end

        ref = mwecmass.geometry.extract_isocurve_at_z(model, zs(k), n_u, cache);
        [A0, Ixx0, Iyy0] = mwecmass.hydrostatics.waterplane_properties(ref);
        [A1, Ixx1, Iyy1] = mwecmass.hydrostatics.waterplane_properties(P);
        worst.area_zk = max(worst.area_zk, abs(A1 - A0) / A0);
        worst.zres_ref = max(worst.zres_ref, max(abs(ref(:, 3) - zs(k))));

        z_ref = mean(ref(:, 3));
        rr = mwecmass.solid.outer_rows(model, cache, z_ref, grids);
        [A2, Ixx2, Iyy2] = mwecmass.hydrostatics.waterplane_properties(rr.pts);
        worst.area = max(worst.area, abs(A2 - A0) / A0);
        worst.ixx = max(worst.ixx, abs(Ixx2 - Ixx0) / Ixx0);
        worst.iyy = max(worst.iyy, abs(Iyy2 - Iyy0) / Iyy0);
    end

    fprintf('seam gap (largest over all heights)          %.3e m\n', worst.gap);
    fprintf('smallest distance between two points         %.3e m\n', worst.mind);
    fprintf('shoelace vs row.area, relative               %.3e\n', worst.shoe);
    fprintf('|z - z_k| of section points                  %.3e m\n', worst.zrow);
    fprintf('mirror symmetry error x->-x, y->-y (point to nearest point) %.3e m\n', worst.sym);
    fprintf('vs extract_isocurve_at_z at z_k: area diff   %.3e (its points miss z_k by up to %.3e m)\n', ...
            worst.area_zk, worst.zres_ref);
    fprintf('vs extract_isocurve_at_z at its own z: area %.3e, Ixx %.3e, Iyy %.3e (relative)\n', ...
            worst.area, worst.ixx, worst.iyy);
    check(worst.shoe <= 4 * eps, 'row.area is the shoelace area');
    check(worst.zrow <= 16 * eps(max(1, max(abs(zs)))), ...
          'section points lie at z_k to 16 ulp (the root tolerance of outer_rows)');
    % 1e-12 m: the same exact point reached from two patches differs only by z/curve rounding (~1e-15)
    check(worst.gap <= 1e-12, 'seam gap above 1e-12 m');
    check(max([worst.area, worst.ixx, worst.iyy]) <= 1e-12, ...
          'same evaluator, same height: sections differ by more than 1e-12');

    zq = zs([3 9 15 21 26]);
    d_eval = 0;
    for k = 1:numel(zq)
        r = mwecmass.solid.outer_rows(model, cache, zq(k), grids);
        names = model.visible_surfs;
        for i = [1:11:r.n, find(r.seam)']
            S = model.eval_surface(names{r.patch(i)}, r.u(i), r.v(i));
            d_eval = max(d_eval, norm(S - r.pts(i, :)));
            if r.seam(i)
                S2 = model.eval_surface(names{r.patch2(i)}, r.u2(i), r.v2(i));
                d_eval = max(d_eval, norm(S2 - r.pts(i, :)));
            end
        end
    end
    fprintf('MS2Parser.eval_surface at the stored (patch,u,v) vs point: %.3e m\n', d_eval);
    check(d_eval <= 1e-12, 'stored (patch,u,v) does not reproduce the point');

    % any surface type outside RuledSurf/RevSurf takes the eval_surface path
    nc = 21;
    c_fast = mwecmass.geometry.precompute_boundary_cache(model, nc);
    c_gen = c_fast;
    c_gen.data = containers.Map(keys(c_fast.data), values(c_fast.data));
    for s = 1:numel(c_fast.sources)
        d = c_gen.data(c_fast.sources{s});
        d.type = 'DevSurf';
        c_gen.data(c_fast.sources{s}) = d;
    end
    zt = [-3.0, -1.0, 0.7];
    rf = mwecmass.solid.outer_rows(model, c_fast, zt);
    rg = mwecmass.solid.outer_rows(model, c_gen, zt);
    dg = 0;
    for k = 1:numel(zt)
        check(rf(k).n == rg(k).n, 'generic path point count');
        dg = max(dg, max(abs(rf(k).pts(:) - rg(k).pts(:))));
    end
    fprintf('eval_surface path vs fast path, n_u = %d: largest point difference %.3e m\n', nc, dg);
    check(dg <= 1e-12, 'eval_surface path differs from the cached path');

    zt300 = zr(1) + (1:300)' / 301 * (zr(2) - zr(1));
    tic;
    c100 = mwecmass.geometry.precompute_boundary_cache(model, n_u);
    r300 = mwecmass.solid.outer_rows(model, c100, zt300);
    t300 = toc;
    check(~any([r300.degenerate]), '300 heights all non-degenerate');
    fprintf('time to build rows for 300 heights (n_u = %d, grids included): %.1f s\n', n_u, t300);
    fprintf('test_outer_rows passed\n');
end

function check(cond, varargin)
    if ~cond
        error('test_outer_rows:fail', varargin{:});
    end
end
