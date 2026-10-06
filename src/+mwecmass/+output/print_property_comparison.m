function print_property_comparison(results, final_props, fids)
%PRINT_PROPERTY_COMPARISON Print the 2D-surrogate / 3D-truth comparison table.
    if nargin < 3, fids = 1; end
    config = results.config;

    mwecmass.output.emit(fids, '\n┌──────────────────────────────────────────────────────────────────┐\n');
    mwecmass.output.emit(fids, '│  PROPERTY COMPARISON: 2D Surrogate / 3D Ground Truth \n');
    mwecmass.output.emit(fids, '├──────────────────────────────────────────────────────────────────┤\n');

    x_final = results.stage2_3d.x_optimal;
    try
        p2d = mwecmass.hydrostatics.properties_2d(x_final, config);
    catch ME
        % Failures degrade gracefully (NaN in the table, warning instead of error).
        warning('mwecmass:output:propertyComparisonSkipped', ...
            '2D-surrogate recompute failed, showing dashes for those columns: %s', ME.message);
        p2d = [];
    end

    p3d = final_props;

    mwecmass.output.emit(fids, '│  %-22s %12s %12s %8s\n', ...
        'Quantity', '2D Surr.', '3D Truth', 'Err(%)');
    mwecmass.output.emit(fids, '│  %s\n', repmat('─', 1, 57));

    if ~isempty(p2d)
        m2 = p2d.mass_total;
    else
        m2 = NaN;
    end
    m3 = p3d.mass_total;
    print_row(fids, 'Mass [kg]', m2, m3, '%.1f');

    if ~isempty(p2d)
        b2 = p2d.mass_buoyant_force;
    else
        b2 = NaN;
    end
    print_row(fids, 'Buoyancy [kg]', b2, p3d.mass_buoyant_force, '%.1f');

    if ~isempty(p2d)
        v2 = p2d.V_sub;
    else
        v2 = NaN;
    end
    print_row(fids, 'V_sub [m^3]', v2, p3d.V_sub, '%.4f');

    if ~isempty(p2d)
        gm2 = p2d.GM;
        gm2_raw = p2d.GM_uncorrected;
    else
        gm2 = NaN;
        gm2_raw = NaN;
    end
    print_row(fids, 'GM [m]', gm2, p3d.GM_L, '%.4f');
    print_row(fids, 'GM_raw [m]', gm2_raw, p3d.GM_L, '%.4f');

    if ~isempty(p2d)
        cg2 = p2d.CG_total(3);
    else
        cg2 = NaN;
    end
    print_row(fids, 'CG_z [m]', cg2, p3d.CG_total(3), '%+.4f');

    if ~isempty(p2d)
        cb2 = p2d.CB(3);
    else
        cb2 = NaN;
    end
    print_row(fids, 'CB_z [m]', cb2, p3d.CB(3), '%+.4f');

    if ~isempty(p2d)
        iyy2 = p2d.Iyy;
    else
        iyy2 = NaN;
    end
    print_row(fids, 'Iyy [kg*m^2]', iyy2, p3d.Iyy, '%.1f');

    if ~isempty(p2d)
        aw2 = p2d.Aw;
    else
        aw2 = NaN;
    end
    print_row(fids, 'Aw [m^2]', aw2, p3d.Aw, '%.4f');

    mwecmass.output.emit(fids, '│  %s\n', repmat('─', 1, 57));
    if ~isempty(p2d)
        th2 = p2d.periods.heave;
        tp2 = p2d.periods.pitch;
    else
        th2 = NaN;
        tp2 = NaN;
    end
    print_row(fids, 'T_heave [s]', th2, p3d.periods.heave, '%.3f');
    print_row(fids, 'T_pitch [s]', tp2, p3d.periods.pitch, '%.3f');

    if ~isinf(p3d.periods.surge) && ~isnan(p3d.periods.surge)
        if ~isempty(p2d)
            ts2 = p2d.periods.surge;
        else
            ts2 = NaN;
        end
        print_row(fids, 'T_surge [s]', ts2, p3d.periods.surge, '%.3f');
    else
        mwecmass.output.emit(fids, '│  %-22s %12s %12s %8s\n', 'T_surge [s]', 'Inf', 'Inf', '—');
    end

    mwecmass.output.emit(fids, '│  %s\n', repmat('─', 1, 57));
    print_row(fids, 'A33 [kg]', NaN, p3d.A33, '%.1f');
    print_row(fids, 'A55 [kg*m^2]', NaN, p3d.A55, '%.1f');

    mwecmass.output.emit(fids, '└──────────────────────────────────────────────────────────────────┘\n');
end

function print_row(fids, label, val_2d, val_3d, fmt)
    fmt_val = fmt;
    err_str = '—';

    if isfinite(val_2d) && isfinite(val_3d) && abs(val_3d) > 1e-12
        err = 100 * (val_2d - val_3d) / val_3d;
        err_str = sprintf('%+.1f', err);
    end

    s2d = fmt_or_dash(val_2d, fmt_val);
    s3d = fmt_or_dash(val_3d, fmt_val);

    mwecmass.output.emit(fids, '│  %-22s %12s %12s %8s\n', label, s2d, s3d, err_str);
end

function s = fmt_or_dash(val, fmt)
    if isnan(val) || isinf(val)
        s = '—';
    else
        s = sprintf(fmt, val);
    end
end
