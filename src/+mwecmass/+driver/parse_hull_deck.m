function hull = parse_hull_deck(ms2_file)
%PARSE_HULL_DECK Parse an MS2 deck and compute geometry/profile mass properties.
% Syntax: hull = mwecmass.driver.parse_hull_deck(ms2_file).
% Input: ms2_file char path. Output hull contains the parser, midplane profile [M x 2] m,
% actual sampled z extents [m], volume [m^3], centroid [1 x 3] m, and second volume moments [m^5].
% Surface evaluation is used for extents; divergence-theorem integration supplies enclosed properties.

    if ~exist(ms2_file, 'file')
        error('WEC:FileNotFound', 'MS2 file not found: %s', ms2_file);
    end

    fprintf('  Parsing geometry: %s\n', ms2_file);
    hull.ms2_model = mwecmass.geometry.MS2Parser.parse(ms2_file);
    hull.profile   = mwecmass.geometry.extract_midplane_profile(hull.ms2_model);

    fprintf('  Visible surfaces: %d\n', length(hull.ms2_model.visible_surfs));
    fprintf('  Profile vertices: %d\n', size(hull.profile, 1));

    % Hull z-extents from surface evaluation
    %
    %  model.extents comes from the .ms2 file header — the bounding box
    %  of ALL entities including B-spline control points.  B-splines
    %  APPROXIMATE control points (don't pass through interior ones),
    %  so control points can lie well outside the actual hull surface.
    %  For C0, a control point at z = -3.0 gives hull_z_min = -3.0,
    %  but the surface only reaches z ≈ -2.5.  This wastes a density
    %  strip on empty space below the keel.
    %
    %  Sample all visible surfaces on a coarse grid and extract
    %  the actual min/max z.  Cost: ~0.05 s (negligible at config time).
    [hull_z_lo, hull_z_hi] = compute_surface_z_range(hull.ms2_model);
    hull.extents_z = [hull_z_lo, hull_z_hi];
    fprintf('  Surface z-range: [%.4f, %.4f] m (file header: [%.4f, %.4f])\n', ...
            hull_z_lo, hull_z_hi, ...
            min(hull.ms2_model.extents(3), hull.ms2_model.extents(6)), ...
            max(hull.ms2_model.extents(3), hull.ms2_model.extents(6)));
    hp_full = mwecmass.hydrostatics.compute_hull(hull.ms2_model, ...
                  struct('n_quad', 20, 'verbose', false));
    hull.volume      = hp_full.volume;
    hull.centroid    = hp_full.centroid;

    % Second volume moments for the HAMS inertia tensor
    hull.int_x2 = hp_full.int_x2;
    hull.int_y2 = hp_full.int_y2;
    hull.int_z2 = hp_full.int_z2;

    fprintf('  Hull volume: %.4f m^3 (parametric divergence theorem)\n', hull.volume);
end

function [z_min, z_max] = compute_surface_z_range(parser, n_sample)
%COMPUTE_SURFACE_Z_RANGE Sample every visible surface on an n_sample×n_sample
% parametric grid and return sampled z-extents. Sampling avoids using .ms2
% B-spline control-point bounds, which may include empty space.

    if nargin < 2, n_sample = 50; end

    z_all = zeros(n_sample * n_sample * length(parser.visible_surfs), 1);
    idx = 0;
    t_s = linspace(0, 1, n_sample);

    for s = 1:length(parser.visible_surfs)
        sname = parser.visible_surfs{s};
        for ui = 1:n_sample
            for vi = 1:n_sample
                pt = parser.eval_surface(sname, t_s(ui), t_s(vi));
                idx = idx + 1;
                z_all(idx) = pt(3);
            end
        end
    end

    z_all = z_all(1:idx);
    z_min = min(z_all);
    z_max = max(z_all);
end
