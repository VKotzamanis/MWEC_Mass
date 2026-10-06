function props = properties_2d(x, config)
%PROPERTIES_2D Strip-theory (2D midplane-extrusion) surrogate for WEC physical properties.
% See docs/METHODS_ENGINE.md#coupled-modal-energy-share.
% Inputs: x = [vertical_shift (m); ballast densities (kg/m^3)]; config from build_config.
% Outputs include hydrostatic geometry/mass properties, stiffness matrices, and periods [s].
% Coupled heave/pitch periods are named by kinetic-energy share; per-axis values are retained.
% Sign/frame conventions: z=0 is the free surface, z up, and positive vertical_shift moves the
% hull up. GM = KM - CG_z. The free-floating model has K_pto = 0 and K_total = K_hydro.
try

%% Section 1 - unpack design vector
%  vertical_shift [m] applied to the 2D profile; positive moves the hull up, so the body-frame
%  waterline moves downward and physical draft decreases.
%  densities_at_nodes [kg/m^3] density at each ballast node (N values).
%  shifted_profile [m x 2] midplane profile in the waterline frame (z=0 is the free surface).
props.vertical_shift = x(1);
densities_at_nodes   = x(2:end)';
shifted_profile      = config.profile + [0, props.vertical_shift];

%% Section 2 - draft-dependent effective width
%  Y-span table stores transverse hull width per z-level (built at configuration time from
%  waterplane polygon slicing). eff_w_floor guards divide-by-zero as the waterplane shrinks.
z_wp_body = -props.vertical_shift;
has_yspan = ~isempty(config.y_span_z_levels) && ...
            ~isempty(config.y_span_table);

if has_yspan
    eff_w = interp1(config.y_span_z_levels, config.y_span_table, ...
                    z_wp_body, 'linear', 'extrap');
    eff_w = max(eff_w, config.eff_w_floor);
else
    eff_w = config.effective_width;
end

%% Section 3 - hydrostatic geometry
%  sub_area [m^2] submerged cross-section area; cb_2d [m x 2] centroid (x,z); wp_width [m]
%  waterplane chord length. V_sub here is a first estimate, overwritten by strip integration below.
try
    [props.sub_area, cb_2d, props.wp_width] = ...
        compute_submerged_properties(shifted_profile);
catch
    props.sub_area = 0;  cb_2d = [0, 0];  props.wp_width = 0;
end

k_vol = config.k_vol;

props.V_sub              = props.sub_area * eff_w * k_vol;
props.Aw                 = props.wp_width * eff_w;
props.effective_width_used = eff_w;
props.mass_buoyant_force = props.V_sub * config.RHO_WATER;
props.CB                 = [cb_2d(1), 0, cb_2d(2)];

%  Waterplane second moment of area about the pitch axis (y-axis); BM = I_wp/V_sub.
try
    wl_pts = mwecmass.geometry.find_waterline_intersections(shifted_profile, 0);

    if ~isempty(wl_pts) && size(wl_pts, 1) >= 2
        x_coords     = wl_pts(:, 1);
        wp_width_calc = max(x_coords) - min(x_coords);
        wp_centroid_x = mean(x_coords);

        I_wp_local = eff_w * wp_width_calc^3 / 12;
        d_x        = wp_centroid_x - props.CB(1);
        I_wp_yy    = I_wp_local + props.Aw * d_x^2;
    else
        I_wp_yy = 0;
    end
catch
    I_wp_yy = 0;
end

% Metacentric height before CG is known (KM = KB + BM)
if props.V_sub > 1e-6
    BM = I_wp_yy / props.V_sub;
else
    BM = 0;
end
KM = props.CB(3) + BM;

%% Section 4 - mass properties via strip integration
%  num_strips horizontal strips from keel to deck; each strip is (profile chord) x (Y-span) x dz.
%  k_vol multiplies every strip volume (PID-driven correction). V_sub is re-computed strip-by-
%  strip and overwrites the Section 3 estimate (accounts for per-strip width variation).
num_strips  = config.n_density_strips;
hull_z_min  = min(shifted_profile(:,2));
hull_z_max  = max(shifted_profile(:,2));
strip_height = (hull_z_max - hull_z_min) / num_strips;
z_strips    = linspace(hull_z_min + strip_height/2, ...
                       hull_z_max - strip_height/2, num_strips)';

% Pre-compute deduplicated shell lookup table (z_levels may have duplicates at strip
% boundaries from linspace overlap in the thin-shell solve). interp1 needs unique sorted points.
if ~isempty(config.shell)
    [z_shell_unique, ia_shell] = unique(config.shell.z_levels, 'stable');
    sf_shell_unique = config.shell.shell_fraction(ia_shell);
end

total_mass   = 0;
cg_numerator = [0, 0];
V_sub_strips = 0;

for i = 1:num_strips
    z_cur = z_strips(i);

    try
        isects = mwecmass.geometry.find_waterline_intersections(shifted_profile, z_cur);
    catch
        continue;
    end
    if isempty(isects) || size(isects, 1) < 2
        continue;
    end

    x_coords    = isects(:, 1);
    strip_width = max(x_coords) - min(x_coords);
    strip_area  = strip_width * strip_height;

    % Per-strip Y-span: hull can taper toward keel or deck
    z_body = z_cur - props.vertical_shift;
    if has_yspan
        strip_eff_w = interp1(config.y_span_z_levels, config.y_span_table, ...
                              z_body, 'linear', 'extrap');
        strip_eff_w = max(strip_eff_w, config.eff_w_floor);
    else
        strip_eff_w = config.effective_width;
    end

    strip_vol = strip_area * strip_eff_w * k_vol;

    % Only strips below z = 0 contribute to displaced volume
    if z_cur <= 0
        V_sub_strips = V_sub_strips + strip_vol;
    end

    % Density at this strip (piecewise-constant when constructability layout is used, else
    % interpolated from the density-node profile).
    z_orig  = z_cur - props.vertical_shift;
    if config.enable_constructability && ~isempty(config.strip_edges)
        bin = discretize(z_orig, config.strip_edges);
        if isnan(bin)
            if z_orig <= config.strip_edges(1)
                bin = 1;
            else
                bin = length(densities_at_nodes);
            end
        end
        rho_raw = densities_at_nodes(bin);
    else
        rho_raw = interp1(config.density_nodes_z, densities_at_nodes, ...
                          z_orig, 'linear', 'extrap');
    end
    rho     = max(config.ballast_density_bounds(1), ...
              min(config.ballast_density_bounds(2), rho_raw));

    % Composite shell branching: rho applies only to the core region when shell is enabled;
    % shell mass uses fixed rho_shell. Reduces to mass = vol*rho when shell is disabled.
    if ~isempty(config.shell)
        sf = interp1(z_shell_unique, sf_shell_unique, ...
                     z_body, 'linear', 'extrap');
        sf = max(0, min(1, sf));

        shell_vol  = strip_vol * sf;
        core_vol   = strip_vol * (1 - sf);
        strip_mass = shell_vol * config.shell.rho_shell + core_vol * rho;
    else
        strip_mass = strip_vol * rho;
    end

    total_mass   = total_mass + strip_mass;
    cg_numerator = cg_numerator + strip_mass * [mean(x_coords), z_cur];
end

% Overwrite Section 3 estimate with strip-integrated displaced volume (more accurate for
% non-prismatic hulls where per-strip width differs from the single waterplane width).
if V_sub_strips > 0
    props.V_sub              = V_sub_strips;
    props.mass_buoyant_force = V_sub_strips * config.RHO_WATER;
end

props.mass_total = total_mass;

if props.mass_total < 1
    warning('mwecmass:hydrostatics:ZeroMass', ...
            'Mass near zero (%.3f kg). Check profile/densities.', props.mass_total);
end

if props.mass_total > 1e-6
    cg_2d          = cg_numerator / props.mass_total;
    % 2D surrogate keeps a non-zero CG_total(1) (cg_2d(1)); the 3D path (properties_3d.m) and
    % the realised path (build_realised_props.m) both force CG_total(1) = 0. This surrogate is
    % not the mode of record (see docs/METHODS_ENGINE.md) and is internally consistent: the
    % 2D Iyy parallel-axis term uses this same props.CG_total.
    props.CG_total = [cg_2d(1), 0, cg_2d(2)];
else
    props.CG_total = [0, 0, 0];
end

%% Section 5 - metacentric height, raw and corrected
%  GM_uncorrected: strip-theory prediction with no PID bias. GM: k_gm applied to CG_z before
%  KM subtraction (CG-depth bias, not a direct GM scale, preserves GM = KM - CG_z).
GM_raw             = KM - props.CG_total(3);
props.GM_uncorrected = GM_raw;

CG_z_corrected = props.CG_total(3) * config.k_gm;
props.GM       = KM - CG_z_corrected;

%% Section 6 - pitch moment of inertia
%  Iyy about the CG via parallel-axis theorem over the same strip decomposition as Section 4
%  (re-iterated because CG was unknown on the first pass). Only Iyy is needed for the 2D surrogate.
Iyy_total = 0;

if props.mass_total > 1e-6
    for i = 1:num_strips
        z_cur = z_strips(i);

        try
            isects = mwecmass.geometry.find_waterline_intersections( ...
                         shifted_profile, z_cur);
        catch
            continue;
        end
        if isempty(isects) || size(isects, 1) < 2
            continue;
        end

        x_coords    = isects(:, 1);
        strip_width = max(x_coords) - min(x_coords);
        strip_area  = strip_width * strip_height;

        z_body = z_cur - props.vertical_shift;
        if has_yspan
            strip_eff_w = interp1(config.y_span_z_levels, config.y_span_table, ...
                                  z_body, 'linear', 'extrap');
            strip_eff_w = max(strip_eff_w, config.eff_w_floor);
        else
            strip_eff_w = config.effective_width;
        end

        z_orig  = z_cur - props.vertical_shift;
        if config.enable_constructability && ~isempty(config.strip_edges)
            bin = discretize(z_orig, config.strip_edges);
            if isnan(bin)
                if z_orig <= config.strip_edges(1)
                    bin = 1;
                else
                    bin = length(densities_at_nodes);
                end
            end
            rho_raw = densities_at_nodes(bin);
        else
            rho_raw = interp1(config.density_nodes_z, densities_at_nodes, ...
                              z_orig, 'linear', 'extrap');
        end
        rho     = max(config.ballast_density_bounds(1), ...
                  min(config.ballast_density_bounds(2), rho_raw));

        % Composite shell branching (same as Section 4)
        strip_vol_iyy = strip_area * strip_eff_w * k_vol;
        if ~isempty(config.shell)
            sf = interp1(z_shell_unique, sf_shell_unique, ...
                         z_body, 'linear', 'extrap');
            sf = max(0, min(1, sf));
            strip_mass = strip_vol_iyy * sf * config.shell.rho_shell + ...
                         strip_vol_iyy * (1 - sf) * rho;
        else
            strip_mass = strip_area * strip_eff_w * k_vol * rho;
        end

        % Parallel-axis: d is vector from CG to strip centroid
        d = [mean(x_coords), 0, z_cur] - props.CG_total;
        Iyy_total = Iyy_total + strip_mass * (d(1)^2 + d(3)^2);
    end
end

props.Iyy = Iyy_total;

%% Section 7 - hydrostatic stiffness
%  K33 (heave) = rho_water * g * Aw. K55 (pitch) = m * g * GM (corrected GM, so the 2D
%  surrogate's period predictions include the PID bias). K11 (surge) = 0, free-floating.
K33_hydro = config.RHO_WATER * config.G * props.Aw;

% K55_hydro = M_total·g·GM (total mass at every site; equals rho_w·V_sub·g·GM at flotation balance).
if props.GM > 0
    K55_hydro = props.mass_total * config.G * props.GM;
else
    K55_hydro = 0;
end

%% Section 8 - WAMIT added mass (full 3x3 matrices)
%  Interpolates WAMIT .1 file data at the current draft: diagonal terms (A11,A33,A55) for
%  uncoupled periods, full 3x3 (A_full,B_full) for coupled eigenvalue analysis. A11/A33/A55/
%  A_full are added mass at omega->infinity (hydro_table.added_mass_inf, build_config.m); B_full is
%  radiation damping, band-averaged (hydro_table.radiation_damping_band_avg) — not added mass, despite being computed
%  in the same call. The frequency qualifier is dropped at this consumer; documentation-only note.
try
    % Pass the 2D strip-theory CG (world frame, z-component) so A55 references the right axis;
    % matches properties_3d and the realisation solvers (one consistent transform path).
    % NOTE: this is the RAW props.CG_total(3), not the k_gm-biased CG_z_corrected used above for
    % GM/K55 -- the A55 congruence transform and the pitch stiffness reference two different CG
    % conventions inside this same function. Both choices are as written; documentation-only
    % note, no arithmetic change.
    [props.A11, props.A33, props.A55, props.A_full, props.B_full] = ...
        mwecmass.bem.interpolate_at_draft( ...
            props.vertical_shift, config, props.CG_total(3));
catch ME
    warning('mwecmass:hydrostatics:WAMITFailed', ...
            'WAMIT interpolation failed: %s. Using zero added mass.', ME.message);
    props.A11 = 0;  props.A33 = 0;  props.A55 = 0;
    props.A_full = zeros(3);  props.B_full = zeros(3);
end

%% Section 9 - total stiffness (no PTO block -- see file header)
%  Free-floating body: no PTO stiffness contribution. K_pto is the fixed zero matrix; K_total
%  equals K_hydro. Surge has zero hydrostatic restoring, so K11_total = 0 -> T_surge = Inf.
K11_total = 0;
K33_total = K33_hydro;
K55_total = K55_hydro;

props.K_hydro = diag([0, K33_hydro, K55_hydro]);
props.K_pto   = diag([0, 0, 0]);
props.K_total = diag([K11_total, K33_total, K55_total]);

%% Section 10 - uncoupled natural periods
%  T = 2*pi*sqrt(M_virtual/K), M_virtual = M_physical + A_added. Pitch uses Iyy, not hull mass.
%  K ~ 0 -> period set to Inf (free mode, no restoring).
M11_virtual = props.mass_total + props.A11;
M33_virtual = props.mass_total + props.A33;
M55_virtual = Iyy_total        + props.A55;

if K11_total > 1e-6
    props.periods.surge = 2*pi * sqrt(M11_virtual / K11_total);
else
    props.periods.surge = inf;
end

if K33_total > 1e-6
    props.periods.heave = 2*pi * sqrt(M33_virtual / K33_total);
else
    props.periods.heave = inf;
end

if K55_total > 1e-6
    props.periods.pitch = 2*pi * sqrt(M55_virtual / K55_total);
else
    props.periods.pitch = inf;
end

%% Section 11 - coupled eigenvalue analysis; heave and pitch periods named by share
%  props.periods.heave and props.periods.pitch are the eigen-periods of the coupled 3-DOF
%  (1 surge, 2 heave, 3 pitch) undamped problem, each named by the DOF holding the largest share
%  of that mode's kinetic energy. The per-axis values from Section 10 are kept beside them as
%  .heave_uncoupled/.pitch_uncoupled for comparison. The coupled solve uses the full mass matrix.
M_phys  = diag([props.mass_total, props.mass_total, Iyy_total]);   % [kg, kg, kg m^2]
M_total = M_phys + props.A_full;   % [kg, kg m, kg m^2] physical + omega->inf added mass
K_full  = props.K_total;           % [N/m, N/m, N m/rad] hydrostatic only, K_pto = 0

props.periods.heave_uncoupled = props.periods.heave;   % [s] per-axis value
props.periods.pitch_uncoupled = props.periods.pitch;   % [s] per-axis value

[T_coupled, share, Phi] = ...
    mwecmass.hydrostatics.coupled_periods_by_share(M_total, K_full);

props.periods.heave         = T_coupled.heave;   % [s] heave-dominated coupled mode
props.periods.pitch         = T_coupled.pitch;   % [s] pitch-dominated coupled mode
props.participation_factors = share;             % [%] kinetic-energy share, rows surge/heave/pitch
props.coupled_periods       = [T_coupled.surge; T_coupled.heave; T_coupled.pitch];   % [s]
props.coupled_modes         = Phi;               % columns ordered surge, heave, pitch
props.surge_per_pitch       = Phi(1,3) / Phi(3,3);   % [m/rad] surge per unit pitch in the pitch mode

%% Section 12 - mass balance
%  mass_discrepancy > 0 -> hull heavier than displaced water; < 0 -> lighter (would pop up).
%  The fmincon equality constraint (owned by +optim) drives this toward zero.
props.mass_discrepancy  = props.mass_total - props.mass_buoyant_force;

%% Section 13 - store density profile
props.densities_at_nodes = densities_at_nodes;

catch ME
    %% Graceful degradation: preserve fields computed before the error, fill the rest with
    %  safe defaults so the caller (fmincon, via +optim) does not receive NaN/Inf.
    if ~exist('props', 'var'), props = struct(); end

    defaults = struct( ...
        'vertical_shift',    x(1), ...
        'mass_total',        0, ...
        'CG_total',          [0 0 0], ...
        'CB',                [0 0 0], ...
        'GM',                0, ...
        'GM_uncorrected',    0, ...
        'V_sub',             0, ...
        'Aw',                0, ...
        'sub_area',          0, ...
        'wp_width',          0, ...
        'mass_buoyant_force', 0, ...
        'A11',               0, ...
        'A33',               0, ...
        'A55',               0, ...
        'A_full',            zeros(3), ...
        'B_full',            zeros(3), ...
        'Iyy',               0, ...
        'mass_discrepancy',  0, ...
        'densities_at_nodes', x(2:end)');

    fnames = fieldnames(defaults);
    for i = 1:numel(fnames)
        if ~isfield(props, fnames{i})
            props.(fnames{i}) = defaults.(fnames{i});
        end
    end

    if ~isfield(props, 'periods')
        props.periods.surge = inf;
        props.periods.heave = inf;
        props.periods.pitch = inf;
    end

    warning('mwecmass:hydrostatics:PartialFailure', ...
            'Partial failure: %s. Returning partial results.', ME.message);
end

end

function [sub_area, cb, wp_width] = compute_submerged_properties(profile)
%COMPUTE_SUBMERGED_PROPERTIES 2D submerged polygon area, centroid, and waterplane chord width.
% Clips the midplane profile polygon below z=0, then computes the submerged area, its
% centroid, and the chord length at z=0 (waterplane). Returns zeros if the hull is entirely
% above water or the clipped polygon degenerates (fewer than 3 vertices).
% Inputs:  profile  [m x 2] double, midplane profile polygon (x, z), z=0 at the free surface
% Outputs: sub_area [m^2] double, signed area of the submerged polygon (0 if none)
%          cb       [1x2] double, centroid [x, z] of the submerged polygon (m)
%          wp_width [m]   double, chord length at z=0 (0 if the waterline does not intersect)
try
    if min(profile(:,2)) >= 0
        sub_area = 0;  cb = [0, 0];  wp_width = 0;
        return;
    end

    sub_poly = mwecmass.geometry.clip_polygon_at_z(profile, 0, 'below');

    if isempty(sub_poly) || size(sub_poly, 1) < 3
        sub_area = 0;  cb = [0, 0];  wp_width = 0;
        return;
    end

    [geom, ~, ~] = mwecmass.hydrostatics.polygon_properties(sub_poly(:,1), sub_poly(:,2));
    sub_area     = geom(1);
    cb           = [geom(2), geom(3)];

    wl_pts = mwecmass.geometry.find_waterline_intersections(profile, 0);
    if ~isempty(wl_pts) && size(wl_pts, 1) >= 2
        wp_width = max(wl_pts(:,1)) - min(wl_pts(:,1));
    else
        wp_width = 0;
    end
catch
    sub_area = 0;  cb = [0, 0];  wp_width = 0;
end
end
