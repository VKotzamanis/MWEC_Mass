function props = calculate_2d_properties(x, config)
% CALCULATE_2D_PROPERTIES  Strip-theory surrogate for WEC physical properties.
%
%   props = CALCULATE_2D_PROPERTIES(x, config)
%
%   This is the FAST (O(ms)) model used by Stage 1 (fmincon inner loop).
%   It approximates the 3D hull as a 2D midplane profile extruded through
%   a draft-dependent effective width (Y-span table).  The 2D surrogate
%   systematically over-predicts submerged volume because rectangular
%   extrusion ignores hull curvature.  The PID correction factors
%   config.k_vol and config.k_gm compensate for this bias.
%
%   COMPUTATION SEQUENCE
%     1. Shift profile vertically by draft (x(1))
%     2. Look up effective width at waterplane from Y-span table
%     3. Compute 2D hydrostatic properties (submerged area, Aw, I_wp)
%     4. Strip integration for mass, CG, and submerged volume
%        — k_vol scales each strip volume (corrects extrusion overestimate)
%     5. GM calculation: raw (uncorrected) and biased (k_gm on CG_z)
%     6. Pitch moment of inertia via parallel-axis theorem
%     7. Hydrostatic + PTO stiffness → uncoupled natural periods
%     8. Coupled eigenvalue analysis (if full 3×3 WAMIT matrices available)
%
%   CORRECTION FACTORS (applied here, driven by PID loop in WEC_Main_Optimizer)
%     k_vol : multiplicative on strip volume.  Corrects the ~12% overestimate
%             from rectangular extrusion of a curved hull cross-section.
%     k_gm  : multiplicative on CG_z.  Corrects the CG depth mismatch
%             between 2D strips and 3D divergence-theorem integration.
%             k_gm > 1 → deeper CG → higher GM (3D CG below 2D CG).
%             k_gm < 1 → shallower CG → lower GM (3D CG above 2D CG).
%
%   INPUTS
%     x      : [1+N × 1] optimisation vector [draft; rho_1 ... rho_N]
%     config : struct from WEC_Configuration_Builder
%
%   OUTPUT
%     props  : struct — see field list at end of file
%
%   See also: calculate_3d_properties, run_2d_optimizer, WEC_Configuration_Builder

try

%% §1  UNPACK DESIGN VECTOR  ──────────────────────────────────────────
%
%  draft           [m]  vertical shift applied to 2D profile.
%                       Positive = hull moves up → deeper draft.
%  densities_at_nodes [kg/m³]  density at each ballast node (N values).
%  shifted_profile [m × 2]    midplane profile in the waterline frame
%                              (z = 0 is the free surface).

props.vertical_shift = x(1);
densities_at_nodes   = x(2:end)';
shifted_profile      = config.profile + [0, props.vertical_shift];

%% §2  DRAFT-DEPENDENT EFFECTIVE WIDTH  ────────────────────────────────
%
%  The Y-span table stores the transverse hull width at each z-level,
%  computed during configuration from waterplane polygon slicing.
%  At the waterplane (z_wp_body), this gives the effective extrusion
%  width for waterplane area (Aw) and second moment of area (I_wp).
%
%  Per-strip widths (used in §4) also interpolate from this table,
%  so that tapered hulls get narrower strips near the keel.
%
%  eff_w_floor prevents division-by-zero when the waterplane shrinks
%  to a sliver at extreme drafts.

z_wp_body = -props.vertical_shift;                        % waterplane in body frame
has_yspan = ~isempty(config.y_span_z_levels) && ...
            ~isempty(config.y_span_table);

if has_yspan
    eff_w = interp1(config.y_span_z_levels, config.y_span_table, ...
                    z_wp_body, 'linear', 'extrap');
    eff_w = max(eff_w, config.eff_w_floor);
else
    eff_w = config.effective_width;                        % constant-width fallback
end

%% §3  HYDROSTATIC GEOMETRY  ──────────────────────────────────────────
%
%  sub_area  [m²]    submerged cross-section area (2D polygon below z = 0)
%  cb_2d     [m × 2] centroid of submerged polygon (x, z)
%  wp_width  [m]     waterplane chord length
%
%  V_sub here is a first estimate (area × width × k_vol).
%  It gets overwritten by the more accurate strip-integrated value in §4.

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

%  Waterplane second moment of area about the pitch axis (y-axis).
%  Uses the polygon-based width + parallel-axis shift to CB_x.
%  I_wp_yy drives the metacentric radius BM = I_wp / V_sub.
try
    wl_pts = WEC_Core_Functions.find_waterline_intersections(shifted_profile, 0);

    if ~isempty(wl_pts) && size(wl_pts, 1) >= 2
        x_coords     = wl_pts(:, 1);
        wp_width_calc = max(x_coords) - min(x_coords);
        wp_centroid_x = mean(x_coords);

        I_wp_local = eff_w * wp_width_calc^3 / 12;          % rectangle about own centroid
        d_x        = wp_centroid_x - props.CB(1);            % shift CB → waterplane centroid
        I_wp_yy    = I_wp_local + props.Aw * d_x^2;         % parallel-axis theorem
    else
        I_wp_yy = 0;
    end
catch
    I_wp_yy = 0;
end

% Metacentric height BEFORE CG is known (KM = KB + BM)
if props.V_sub > 1e-6
    BM = I_wp_yy / props.V_sub;
else
    BM = 0;
end
KM = props.CB(3) + BM;

%% §4  MASS PROPERTIES — STRIP INTEGRATION  ───────────────────────────
%
%  The hull is sliced into num_strips horizontal strips from keel to deck.
%  Each strip is a rectangle: (profile chord) × (Y-span at that z) × dz.
%
%  k_vol multiplies every strip volume.  This is where the PID-driven
%  correction enters the mass/buoyancy calculation — the solver sees a
%  corrected displaced volume while the density field stays physical.
%
%  Submerged volume is re-computed strip-by-strip (V_sub_strips) and
%  overwrites the §3 estimate because strip integration accounts for
%  per-strip width variation that the single-width estimate misses.
%
%  VARIABLE DICTIONARY (loop-scope)
%    z_strips      [m × num_strips]  strip midpoint elevations (shifted frame)
%    strip_height  [m]               uniform strip thickness
%    strip_eff_w   [m]               Y-span at this strip's body-frame z
%    strip_vol     [m³]              strip volume (chord × width × dz × k_vol)
%    rho           [kg/m³]           clamped interpolated density at strip z
%    strip_mass    [kg]              strip_vol × rho

num_strips  = config.n_density_strips;
hull_z_min  = min(shifted_profile(:,2));
hull_z_max  = max(shifted_profile(:,2));
strip_height = (hull_z_max - hull_z_min) / num_strips;
z_strips    = linspace(hull_z_min + strip_height/2, ...
                       hull_z_max - strip_height/2, num_strips)';

% Pre-compute deduplicated shell lookup table (z_levels may have
% duplicates at strip boundaries from linspace overlap in
% WEC_Shell_Offset.compute).  interp1 requires unique sorted points.
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
        isects = WEC_Core_Functions.find_waterline_intersections(shifted_profile, z_cur);
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

    % Density at this strip (piecewise-constant or interpolated)
    %
    %  When constructability is enabled, use piecewise-constant lookup
    %  to match the realised hull's discrete density per strip.
    %  See calculate_3d_properties §3 header for full rationale.
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

    % --- Composite shell branching ---
    %  When shell is enabled, rho applies only to the core region.
    %  The shell volume fraction at this z is interpolated from the
    %  pre-computed lookup table.  Shell mass uses fixed rho_shell.
    %  When shell is disabled, this reduces to the original: mass = vol × rho.
    if ~isempty(config.shell)
        sf = interp1(z_shell_unique, sf_shell_unique, ...
                     z_body, 'linear', 'extrap');
        sf = max(0, min(1, sf));   % clamp to [0, 1]

        shell_vol  = strip_vol * sf;
        core_vol   = strip_vol * (1 - sf);
        strip_mass = shell_vol * config.shell.rho_shell + core_vol * rho;
    else
        strip_mass = strip_vol * rho;
    end

    total_mass   = total_mass + strip_mass;
    cg_numerator = cg_numerator + strip_mass * [mean(x_coords), z_cur];
end

% Overwrite §3 estimate with strip-integrated displaced volume.
% This is more accurate for non-prismatic hulls where per-strip
% width differs from the single waterplane width.
if V_sub_strips > 0
    props.V_sub              = V_sub_strips;
    props.mass_buoyant_force = V_sub_strips * config.RHO_WATER;
end

props.mass_total = total_mass;

if props.mass_total < 1
    warning('calculate_2d_properties:ZeroMass', ...
            'Mass near zero (%.3f kg). Check profile/densities.', props.mass_total);
end

if props.mass_total > 1e-6
    cg_2d          = cg_numerator / props.mass_total;
    props.CG_total = [cg_2d(1), 0, cg_2d(2)];
else
    props.CG_total = [0, 0, 0];
end

%% §5  METACENTRIC HEIGHT (RAW AND CORRECTED)  ────────────────────────
%
%  GM_uncorrected is the strip-theory prediction with no PID bias.
%  GM (corrected) applies k_gm to the CG_z before subtracting from KM.
%
%  WHY bias CG instead of scaling GM directly?
%    The physical mechanism is a depth error: the 2D strips place the
%    CG at a different depth than the 3D divergence-theorem integral.
%    Scaling CG_z by k_gm preserves the correct relationship
%    GM = KM − CG_z  and avoids masking whether the error comes from
%    buoyancy (KB) or mass distribution (CG).

GM_raw             = KM - props.CG_total(3);
props.GM_uncorrected = GM_raw;

CG_z_corrected = props.CG_total(3) * config.k_gm;
props.GM       = KM - CG_z_corrected;

%% §6  PITCH MOMENT OF INERTIA  ──────────────────────────────────────
%
%  Iyy about the CG, computed by parallel-axis theorem over the same
%  strip decomposition used in §4.  Re-iterates over strips because
%  the CG was not known during the first pass.
%
%  Only Iyy (pitch) is needed for the 2D surrogate — Ixx and Izz are
%  computed in the 3D model where the full mesh is available.

Iyy_total = 0;

if props.mass_total > 1e-6
    for i = 1:num_strips
        z_cur = z_strips(i);

        try
            isects = WEC_Core_Functions.find_waterline_intersections( ...
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

        % --- Composite shell branching (same as §4) ---
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

%% §7  HYDROSTATIC STIFFNESS  ─────────────────────────────────────────
%
%  K33 (heave) = ρ g Aw          — restoring force per unit heave displacement
%  K55 (pitch) = m g GM           — restoring moment per unit pitch angle
%                                   Uses CORRECTED GM so the 2D surrogate's
%                                   period predictions include the PID bias.
%  K11 (surge) = 0 in free-floating mode (no hydrostatic restoring in surge).

K33_hydro = config.RHO_WATER * config.G * props.Aw;

if props.GM > 0
    K55_hydro = props.mass_total * config.G * props.GM;
else
    K55_hydro = 0;
end

%% §8  WAMIT ADDED MASS (FULL 3×3 MATRICES)  ──────────────────────────
%
%  Interpolates WAMIT .1 file data at the current draft.
%  Returns diagonal terms (A11, A33, A55) for uncoupled periods AND
%  the full 3×3 matrix (A_full, B_full) for coupled eigenvalue analysis.
%  A13 ≠ 0 → surge–pitch coupling (asymmetric hull about the waterplane).

try
    % Pass the 2D strip-theory CG (world frame, z-component) so A55
    % references the right axis.  Matches calculate_3d_properties and
    % the realisation solvers — one consistent transform path.
    [props.A11, props.A33, props.A55, ~, props.A_full, props.B_full] = ...
        WEC_Core_Functions.interpolate_wamit_added_mass( ...
            props.vertical_shift, config, props.CG_total(3));
catch ME
    warning('calculate_2d_properties:WAMITFailed', ...
            'WAMIT interpolation failed: %s. Using zero added mass.', ME.message);
    props.A11 = 0;  props.A33 = 0;  props.A55 = 0;
    props.A_full = zeros(3);  props.B_full = zeros(3);
end

%% §9  PTO STIFFNESS  ─────────────────────────────────────────────────
%
%  In PTO-augmented mode, the Power Take-Off mechanism adds extra
%  restoring stiffness to each DOF.  The PTO is modelled as a pair of
%  symmetric linear dashpot–springs at angle pto_angle_deg from vertical.
%
%  K_pto (base) is scaled to the hydrostatic heave stiffness so that
%  PTO effects remain proportional to hull size.
%
%  In free-floating mode (enable_pto_effects = 0), K_pto = 0 for all DOF.
%  This is the default during mass-distribution optimisation; PTO coupling
%  is studied later in the kinematic analysis suite.

if config.enable_pto_effects == 1
    K_pto = K33_hydro;
    a_rad = deg2rad(config.pto_angle_deg);

    K11_pto = 2 * K_pto * sin(a_rad)^2;
    K33_pto = 2 * K_pto * cos(a_rad)^2;
    K55_pto = 2 * K_pto * (config.ht * cos(a_rad) - config.bt * sin(a_rad))^2;
else
    K11_pto = 0;  K33_pto = 0;  K55_pto = 0;
end

%% §10  TOTAL STIFFNESS  ──────────────────────────────────────────────
%
%  Surge has zero hydrostatic restoring — only PTO contributes.
%  In free-floating mode K11 = 0 → T_surge = Inf (physically correct:
%  a floating body drifts freely in surge).

K11_total = K11_pto;
K33_total = K33_hydro + K33_pto;
K55_total = K55_hydro + K55_pto;

props.K_hydro = diag([0, K33_hydro, K55_hydro]);
props.K_pto   = diag([K11_pto, K33_pto, K55_pto]);
props.K_total = diag([K11_total, K33_total, K55_total]);

%% §11  UNCOUPLED NATURAL PERIODS  ────────────────────────────────────
%
%  T = 2π √(M_virtual / K)   where  M_virtual = M_physical + A_added
%
%  For pitch: M_physical = Iyy (not hull mass).
%  When K ≈ 0 the period is set to Inf (free mode, no restoring).

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

%% §12  COUPLED EIGENVALUE ANALYSIS  ──────────────────────────────────
%
%  When full 3×3 WAMIT matrices are available, solve the generalised
%  eigenvalue problem:
%
%      K φ = ω² M φ    →    (M⁻¹ K) φ = ω² φ
%
%  This captures surge–pitch coupling (A13 ≠ 0) which shifts the
%  natural periods relative to the uncoupled estimates above.
%
%  The 2D version uses M⁻¹K (standard eigenvalue) rather than
%  eig(K, M) (generalised) — acceptable here because M is always
%  well-conditioned in the 2D surrogate.  The 3D model uses the
%  generalised form for robustness.

if isfield(props, 'A_full') && ~isempty(props.A_full) && any(props.A_full(:) ~= 0)
    try
        M_phys     = diag([props.mass_total, props.mass_total, Iyy_total]);
        M_total    = M_phys + props.A_full;
        K_full     = props.K_total;

        if det(M_total) > 1e-12
            [V_eig, D_eig] = eig(M_total \ K_full);
            omega_sq        = diag(D_eig);

            valid = omega_sq > 1e-6;
            if any(valid)
                omega_n = sqrt(omega_sq(valid));
                T_n     = 2*pi ./ omega_n;

                [T_sorted, si]          = sort(T_n, 'descend');
                props.coupled_periods   = T_sorted;
                props.coupled_modes     = V_eig(:, valid);
                props.coupled_modes     = props.coupled_modes(:, si);
                props.participation_factors = ...
                    compute_participation_factors(props.coupled_modes, M_total);
            end
        end
    catch
        % Coupled analysis failed — uncoupled periods already set above.
    end
end

%% §13  MASS BALANCE  ─────────────────────────────────────────────────
%  mass_discrepancy > 0  → hull is heavier than its displaced water
%  mass_discrepancy < 0  → hull is lighter (would pop up)
%  The fmincon equality constraint drives this toward zero.

props.mass_discrepancy  = props.mass_total - props.mass_buoyant_force;

%% §14  STORE DENSITY PROFILE  ────────────────────────────────────────
props.densities_at_nodes = densities_at_nodes;

catch ME
    %% GRACEFUL DEGRADATION ────────────────────────────────────────────
    %  Preserve any fields that were successfully computed before the
    %  error, and fill missing fields with safe defaults so the solver
    %  can continue (returning NaN/Inf would crash fmincon).

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

    warning('calculate_2d_properties:PartialFailure', ...
            'Partial failure: %s. Returning partial results.', ME.message);
end

end  % calculate_2d_properties


%% ═════════════════════════════════════════════════════════════════════
%%  LOCAL HELPER: SUBMERGED 2D GEOMETRY
%% ═════════════════════════════════════════════════════════════════════

function [sub_area, cb, wp_width] = compute_submerged_properties(profile)
% COMPUTE_SUBMERGED_PROPERTIES  2D submerged polygon area, centroid, Aw chord.
%
%   Clips the midplane profile polygon below z = 0, then computes:
%     sub_area  [m²]    signed area of submerged polygon
%     cb        [1×2]   centroid [x, z] of submerged polygon
%     wp_width  [m]     chord length at z = 0 (waterplane)
%
%   Returns zeros if the hull is entirely above water or the polygon
%   degenerates (fewer than 3 vertices after clipping).

try
    if min(profile(:,2)) >= 0
        sub_area = 0;  cb = [0, 0];  wp_width = 0;
        return;
    end

    sub_poly = WEC_Core_Functions.clipPolygon(profile, 0, 'below');

    if isempty(sub_poly) || size(sub_poly, 1) < 3
        sub_area = 0;  cb = [0, 0];  wp_width = 0;
        return;
    end

    [geom, ~, ~] = WEC_Core_Functions.polygeom(sub_poly(:,1), sub_poly(:,2));
    sub_area     = geom(1);
    cb           = [geom(2), geom(3)];

    wl_pts = WEC_Core_Functions.find_waterline_intersections(profile, 0);
    if ~isempty(wl_pts) && size(wl_pts, 1) >= 2
        wp_width = max(wl_pts(:,1)) - min(wl_pts(:,1));
    else
        wp_width = 0;
    end
catch
    sub_area = 0;  cb = [0, 0];  wp_width = 0;
end
end


%% ═════════════════════════════════════════════════════════════════════
%%  LOCAL HELPER: MODAL PARTICIPATION FACTORS
%% ═════════════════════════════════════════════════════════════════════

function PF = compute_participation_factors(modes, M)
% COMPUTE_PARTICIPATION_FACTORS  DOF contribution to each coupled mode.
%
%   PF(i,j) = φ_i² M_ii / (φ' M φ) × 100     [%]
%
%   Rows: surge (1), heave (2), pitch (3).   Columns: modes sorted by
%   descending period.  Each column sums to ~100%.
%
%   INPUTS
%     modes : [3×N]  mode shapes (columns = eigenvectors)
%     M     : [3×3]  total mass matrix (physical + added)
%
%   OUTPUT
%     PF    : [3×N]  participation factors in percent

try
    n_modes = size(modes, 2);
    PF      = zeros(3, n_modes);

    for j = 1:n_modes
        phi        = modes(:, j);
        modal_mass = phi' * M * phi;

        if abs(modal_mass) > 1e-12
            for i = 1:3
                PF(i, j) = (phi(i)^2 * M(i,i)) / modal_mass * 100;
            end
        end
    end
catch
    PF = zeros(3, size(modes, 2));
end
end