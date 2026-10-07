function rec = ms2_parser_record(model, opts)
%MS2_PARSER_RECORD Evaluate a parsed deck through the public MS2Parser API and keep every output.
%   rec = ms2_parser_record(model, opts) with model = MS2Parser.parse(deck) and opts fields
%     sg      surface parameter values (u and v grid)
%     tg      curve parameter values
%     n_u     boundary-cache sample counts, one iso-z run per entry
%     n_z     iso-z heights per entry of n_u (the first run also gets the edge heights)
%   Each call is stored as a struct with label, ok, out (cell of outputs) and, when it threw,
%   the error identifier and message, so a later run can be compared call by call.

    rec = struct('labels', {{}}, 'calls', {{}});
    rec.model = model_summary(model);

    names = model.entities.keys();
    z_min = inf;
    z_max = -inf;

    for i = 1:numel(names)
        nm = names{i};
        e = model.entities(nm);
        if e.is_surface
            continue;
        end
        t = opts.tg(:);
        rec = add(rec, ['eval_point ' nm], @() model.eval_point(nm), 1);
        rec = add(rec, ['eval_any_point ' nm], @() model.eval_any_point(nm), 1);
        rec = add(rec, ['eval_curve ' nm], @() model.eval_curve(nm, t), 1);
        rec = add(rec, ['eval_snake ' nm], @() model.eval_snake(nm, t), 1);
        rec = add(rec, ['eval_curve_or_snake ' nm], @() model.eval_curve_or_snake(nm, t), 1);
        rec = add(rec, ['eval_curve_with_deriv ' nm], @() model.eval_curve_with_deriv(nm, t), 2);
        rec = add(rec, ['eval_snake_with_deriv ' nm], @() model.eval_snake_with_deriv(nm, t), 2);
        rec = add(rec, ['eval_curve scalar ' nm], @() model.eval_curve(nm, 0.37), 1);
        rec = add(rec, ['eval_snake scalar ' nm], @() model.eval_snake(nm, 0.37), 1);
    end

    n_s = numel(opts.sg);
    rec = add_type_evaluators(rec, model, names, opts);

    for i = 1:numel(model.visible_surfs)
        nm = model.visible_surfs{i};
        pts = zeros(n_s * n_s, 3);
        S = pts; Su = pts; Sv = pts;
        failed = '';
        k = 0;
        for a = 1:n_s
            for b = 1:n_s
                k = k + 1;
                try
                    pts(k, :) = model.eval_surface(nm, opts.sg(a), opts.sg(b));
                    [S(k, :), Su(k, :), Sv(k, :)] = ...
                        model.eval_surface_with_derivs(nm, opts.sg(a), opts.sg(b));
                catch err
                    failed = sprintf('%s | %s', err.identifier, err.message);
                    break;
                end
            end
            if ~isempty(failed), break; end
        end
        rec = push(rec, ['surface grid ' nm], struct('ok', isempty(failed), 'out', {{pts, S, Su, Sv}}, ...
                   'id', failed, 'msg', ''));
        if isempty(failed)
            z_min = min(z_min, min(pts(:, 3)));
            z_max = max(z_max, max(pts(:, 3)));
        end
        ug = opts.sg(1:3:end);
        rec = add(rec, ['eval_surface_grid ' nm], @() model.eval_surface_grid(nm, ug, ug), 1);
    end
    for i = 1:numel(names)
        nm = names{i};
        e = model.entities(nm);
        if e.is_surface && ~any(strcmp(nm, model.visible_surfs))
            rec = add(rec, ['eval_surface hidden ' nm], @() model.eval_surface(nm, 0.3, 0.6), 1);
            rec = add(rec, ['eval_surface_with_derivs hidden ' nm], ...
                      @() model.eval_surface_with_derivs(nm, 0.3, 0.6), 3);
        end
    end
    rec = add(rec, 'eval_surface missing', @() model.eval_surface('no_such_surface', 0.3, 0.6), 1);
    rec = add(rec, 'eval_surface on a curve', @() model.eval_surface('bc1', 0.3, 0.6), 1);
    rec = add(rec, 'eval_surface_with_derivs missing', ...
              @() model.eval_surface_with_derivs('no_such_surface', 0.3, 0.6), 3);
    rec = add(rec, 'eval_curve missing', @() model.eval_curve('no_such_curve', 0.5), 1);
    rec = add(rec, 'eval_snake missing', @() model.eval_snake('no_such_curve', 0.5), 1);
    rec = add(rec, 'eval_curve_or_snake missing', @() model.eval_curve_or_snake('no_such_curve', 0.5), 1);
    rec = add(rec, 'eval_curve_with_deriv missing', @() model.eval_curve_with_deriv('no_such_curve', 0.5), 2);
    rec = add(rec, 'eval_snake_with_deriv missing', @() model.eval_snake_with_deriv('no_such_curve', 0.5), 2);
    rec = add(rec, 'eval_point missing', @() model.eval_point('no_such_point'), 1);
    rec = add(rec, 'eval_any_point missing', @() model.eval_any_point('no_such_point'), 1);

    rec = add(rec, 'classify_visible_surfaces', @() model.classify_visible_surfaces(), 1);
    for i = 1:numel(model.visible_surfs)
        nm = model.visible_surfs{i};
        rec = add(rec, ['get_required_entities ' nm], @() model.get_required_entities(nm), 1);
        rec = add(rec, ['resolve_mirror_chain ' nm], @() model.resolve_mirror_chain(nm), 2);
        rec = add(rec, ['get_boundary_entities ' nm], @() model.get_boundary_entities(nm), 1);
    end

    if isfinite(z_min)
        for r = 1:numel(opts.n_u)
            n_u = opts.n_u(r);
            z = linspace(z_min + 1e-4, z_max - 1e-4, opts.n_z(r));
            if r == 1
                z = [z, z_min, z_max, z_min - 0.05, z_max + 0.05, 0.5 * (z_min + z_max)];
            end
            [cache, rec] = record_cache(model, n_u, rec);
            if isempty(cache), continue; end
            rec.heights{r} = z;
            for k = 1:numel(z)
                rec = add(rec, sprintf('extract_isocurve_at_z n_u=%d z=%.17g', n_u, z(k)), ...
                          @() mwecmass.geometry.extract_isocurve_at_z(model, z(k), n_u, cache), 1);
            end
        end
    end
end

function rec = add_type_evaluators(rec, model, names, opts)
% The public evaluators named after an entity type take the entity struct of model.entities.
    curve_fns = struct( ...
        'BCurve',     {{'eval_bcurve', 'eval_bcurve_deriv'}}, ...
        'Conic',      {{'eval_conic', 'eval_conic_deriv'}}, ...
        'CopyCurve',  {{'eval_copy_curve', 'eval_copy_curve_deriv'}}, ...
        'Line',       {{'eval_line', 'eval_line_deriv'}}, ...
        'BSubCurve',  {{'eval_bsub_curve', 'eval_bsub_curve_deriv'}}, ...
        'Arc',        {{'eval_arc', 'eval_arc_deriv'}}, ...
        'PolyCurve2', {{'eval_polycurve2', 'eval_polycurve2_deriv'}}, ...
        'ProjCurve',  {{'eval_proj_curve', 'eval_proj_curve_deriv'}}, ...
        'EdgeSnake',  {{'eval_edge_snake', 'eval_edge_snake_deriv'}}, ...
        'BSubSnake',  {{'eval_bsub_snake', 'eval_bsub_snake_deriv'}});
    surf_fns = struct( ...
        'RuledSurf', {{'eval_ruled_surf', 'eval_ruled_surf_derivs'}}, ...
        'RevSurf',   {{'eval_rev_surf', 'eval_rev_surf_derivs'}}, ...
        'BLoftSurf', {{'eval_bloft_surf', 'eval_bloft_surf_derivs'}}, ...
        'DevSurf',   {{'eval_dev_surf', 'eval_dev_surf_derivs'}}, ...
        'MirrSurf',  {{'eval_mirr_surf', 'eval_mirr_surf_derivs'}});
    t = opts.tg(:);
    uv = [0.3, 0.6; 1/3, 2/3; 0, 1; 1 - 1e-11, 0.25];
    for i = 1:numel(names)
        nm = names{i};
        e = model.entities(nm);
        if isfield(curve_fns, e.type)
            fn = curve_fns.(e.type);
            rec = add(rec, [fn{1} ' entity ' nm], @() feval(fn{1}, model, e, t), 1);
            rec = add(rec, [fn{2} ' entity ' nm], @() feval(fn{2}, model, e, t), 2);
            rec = add(rec, [fn{1} ' entity scalar ' nm], @() feval(fn{1}, model, e, 0.37), 1);
        elseif isfield(surf_fns, e.type)
            fn = surf_fns.(e.type);
            for k = 1:size(uv, 1)
                u = uv(k, 1);
                v = uv(k, 2);
                rec = add(rec, sprintf('%s entity %s uv%d', fn{1}, nm, k), @() feval(fn{1}, model, e, u, v), 1);
                rec = add(rec, sprintf('%s entity %s uv%d', fn{2}, nm, k), @() feval(fn{2}, model, e, u, v), 3);
            end
            if strcmp(e.type, 'BLoftSurf')
                rec = add(rec, ['eval_bloft_surf_at_u entity ' nm], ...
                          @() model.eval_bloft_surf_at_u(e, 0.3, linspace(0, 1, 11)), 1);
            end
        end
    end
end

function [cache, rec] = record_cache(model, n_u, rec)
    cache = [];
    try
        cache = mwecmass.geometry.precompute_boundary_cache(model, n_u);
    catch err
        rec = push(rec, sprintf('precompute_boundary_cache n_u=%d', n_u), ...
                   struct('ok', false, 'out', {{}}, 'id', err.identifier, 'msg', err.message));
        return;
    end
    k = cache.data.keys();
    v = cache.data.values();
    flat = struct('sources', {cache.sources}, 'mirrors', cache.mirrors, ...
                  'u_samples', cache.u_samples, 'data_keys', {k}, 'data_values', {v});
    rec = push(rec, sprintf('precompute_boundary_cache n_u=%d', n_u), ...
               struct('ok', true, 'out', {{flat}}, 'id', '', 'msg', ''));
end

function s = model_summary(model)
    k = model.entities.keys();
    v = model.entities.values();
    [~, fname, fext] = fileparts(model.filename);
    s = struct('keys', {k}, 'values', {v}, 'visible_surfs', {model.visible_surfs}, ...
               'units', model.units, 'extents', model.extents, 'filename', [fname fext], ...
               'file_symmetry', {model.file_symmetry});
end

function rec = add(rec, label, fn, n_out)
    out = cell(1, n_out);
    try
        [out{:}] = fn();
        r = struct('ok', true, 'out', {out}, 'id', '', 'msg', '');
    catch err
        r = struct('ok', false, 'out', {{}}, 'id', err.identifier, 'msg', err.message);
    end
    rec = push(rec, label, r);
end

function rec = push(rec, label, r)
    rec.labels{end + 1} = label;
    rec.calls{end + 1} = r;
end
