function [A_inner, Iyy_inner, Ixx_inner] = inner_properties_at_z( ...
        config, z_level, t_steel, max_slope_factor, n_u, ...
        A_outer, Iyy_outer, Ixx_outer)
%INNER_PROPERTIES_AT_Z Compute inner (void) area and second moments at z.
%   Extracts and orders the outer contour, offsets it inward by the
%   slope-corrected steel thickness, and evaluates polygon properties.
%   A_inner is [m^2], Iyy_inner=∫∫x^2 dA and Ixx_inner=∫∫y^2 dA [m^4].
%   See docs/METHODS_ENGINE.md#thin-shell-slope-offset

    wl_pts = mwecmass.geometry.extract_isocurve_at_z( ...
                 config.ms2_model, z_level, n_u, config.boundary_cache);

    if isempty(wl_pts) || size(wl_pts, 1) < 3
        error('mwecmass:thin_shell:NoIsocurveAtZ', ...
              ['extract_isocurve_at_z returned <3 contour points at ', ...
               'z=%g where Aw_table reports A_outer=%g (>0).  ', ...
               'Check that config.boundary_cache and config.Aw_table* ', ...
               'were built for the same ms2_model.'], z_level, A_outer);
    end

    [~, ~, ~, pts_ord] = mwecmass.hydrostatics.waterplane_properties(wl_pts);
    x_poly = pts_ord(:, 1);
    y_poly = pts_ord(:, 2);

    x_cl = [x_poly; x_poly(1)];
    y_cl = [y_poly; y_poly(1)];
    P_outer = sum(sqrt(diff(x_cl).^2 + diff(y_cl).^2));

    cos_alpha = mwecmass.realise.thin_shell.hull_slope_cos_at_z( ...
                    config, z_level, A_outer, P_outer);
    cos_alpha = max(cos_alpha, 1.0 / max_slope_factor);
    offset_dist = t_steel / cos_alpha;

    % Use polygon_properties on the miter polygon for BOTH area and second moments
    % (consistent: same boundary for mass and inertia integrals).
    [x_off, y_off] = mwecmass.internal.offset_polygon( ...
                          x_poly, y_poly, offset_dist);
    if length(x_off) >= 3
        [geom_in, iner_in, ~] = mwecmass.hydrostatics.polygon_properties(x_off, y_off);
        A_miter = geom_in(1);
        % Validity guard: collapsed polygon (A<=0) or acute-corner overshoot
        % (A_miter >= A_outer occurs when interior angles are acute and the
        % miter vertex moves past the outer boundary — physically invalid).
        if A_miter <= 1e-10 || A_miter >= A_outer
            A_inner = 0;  Iyy_inner = 0;  Ixx_inner = 0;
            return;
        end
        A_inner   = A_miter;
        Iyy_inner = abs(iner_in(2));   % ∫∫ x² dA about origin
        Ixx_inner = abs(iner_in(1));   % ∫∫ y² dA about origin
        if ~isfinite(Iyy_inner)
            Iyy_inner = Iyy_outer * (A_inner / A_outer)^2;
        end
        if ~isfinite(Ixx_inner)
            Ixx_inner = Ixx_outer * (A_inner / A_outer)^2;
        end
    else
        A_inner = 0;  Iyy_inner = 0;  Ixx_inner = 0;
        return;
    end

    Iyy_inner = min(Iyy_inner, Iyy_outer);
    Ixx_inner = min(Ixx_inner, Ixx_outer);
end
