function v_root = bspline_param_at_z(knots, z_ctrl, degree, z_wl)
%BSPLINE_PARAM_AT_Z Find the B-spline parameter at a target elevation.
%   Inputs define the scalar spline; z_wl is in m. The first bracketed crossing is bisected.

    spans = unique(knots);
    n_spans = length(spans) - 1;
    v_root = NaN;

    z_at_spans = zeros(length(spans), 1);
    for j = 1:length(spans)
        z_at_spans(j) = eval_scalar_bspline( ...
                            knots, z_ctrl, degree, spans(j));
    end

    f_spans = z_at_spans - z_wl;

    % Count sign changes before bisection; non-monotone profiles use the first.
    n_crossings = sum(f_spans(1:end-1) .* f_spans(2:end) <= 0);
    if n_crossings > 1
        warning('mwecmass:geometry:MultipleZCrossings', ...
            ['bspline_param_at_z: %d sign changes found for z_wl=%.4f. ' ...
             'Only the first is returned. Check for non-monotone ' ...
             'BLoftSurf z-profile.'], n_crossings, z_wl);
    end

    for s = 1:n_spans
        if f_spans(s) * f_spans(s + 1) <= 0
            a = spans(s);
            b = spans(s + 1);
            fa = f_spans(s);

            for iter = 1:50
                m = (a + b) / 2;
                fm = eval_scalar_bspline( ...
                         knots, z_ctrl, degree, m) - z_wl;
                if abs(fm) < 1e-12, break; end
                if fa * fm < 0
                    b = m;
                else
                    a = m; fa = fm;
                end
            end

            v_root = (a + b) / 2;
            return;
        end
    end
end

function z = eval_scalar_bspline(knots, z_ctrl, degree, t)
%EVAL_SCALAR_BSPLINE Evaluate a scalar B-spline at parameter t.
%   knots, z_ctrl, and degree define the spline; t is in [0,1]. Returns the
%   scalar value z(t), clamping endpoint evaluations to the end control values.
    if t >= 1 - 1e-10
        z = z_ctrl(end);
        return;
    end
    if t <= 1e-10
        z = z_ctrl(1);
        return;
    end
    n_ctrl = length(z_ctrl);
    basis = mwecmass.geometry.MS2Parser.bspline_basis_all( ...
                knots, degree, min(max(t, 0), 1 - 1e-12), n_ctrl);
    z = basis * z_ctrl(:);
end
