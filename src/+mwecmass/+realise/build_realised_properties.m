function realised = build_realised_properties(final_props, steel_data, config)
%BUILD_REALISED_PROPERTIES Rebuild hydrostatics and dynamics from realised geometry.
% steel_data supplies feasible as-built mass, CG, inertia, and hydrostatics; config supplies tables.
% Added-mass coefficients are interpolated at the realised draft/CG and coupled periods are recomputed.
% If steel_data is infeasible or missing its flag, final_props is returned unchanged.
% See docs/METHODS_ENGINE.md#realise-property-rebuild
    if ~isstruct(steel_data) || ~isfield(steel_data, 'feasible') || ~steel_data.feasible
        realised = final_props;
        return;
    end

    realised = struct();

    % Body-frame placement
    realised.vertical_shift = steel_data.vertical_shift;
    realised.draft          = steel_data.draft;

    % As-built hydrostatics with table interpolation where required
    z_wl_body = -steel_data.vertical_shift;
    z_sub_top = min(z_wl_body, config.hull_z_max);

    realised.Aw      = steel_data.Aw;
    realised.I_wp_yy = steel_data.I_wp_yy;
    if isfield(config, 'I_wp_xx_table') && ~isempty(config.I_wp_xx_table)
        realised.I_wp_xx = max(0, interp1(config.Aw_table_z, ...
                                           config.I_wp_xx_table, ...
                                           z_sub_top, 'linear', 0));
    else
        realised.I_wp_xx = realised.I_wp_yy;   % symmetric-hull fallback
    end
    realised.V_sub = steel_data.V_sub;
    realised.CB    = [0, 0, steel_data.CB_z_world];

    % Wetted surface area (handle missing S_wet_table)
    if isfield(config, 'S_wet_table') && ~isempty(config.S_wet_table)
        realised.A_sub = max(0, interp1(config.Aw_table_z, ...
                                         config.S_wet_table, ...
                                         z_sub_top, 'linear', 0));
    elseif isfield(final_props, 'A_sub')
        realised.A_sub = final_props.A_sub;
    else
        realised.A_sub = 0;
    end

    realised.mass_buoyant_force = config.RHO_WATER * realised.V_sub;
    realised.KM                 = steel_data.KM_world;

    % ── Mass / CG / Inertia tensor (full 3×3) ─────────────────
    realised.mass_total = steel_data.M_total;
    realised.CG_total   = [0, 0, steel_data.CG_z_world];

    realised.Ixx = steel_data.Ixx_about_cg;
    realised.Iyy = steel_data.Iyy_about_cg;
    realised.Izz = steel_data.Izz_about_cg;
    realised.Inertia_Tensor = diag([realised.Ixx, realised.Iyy, realised.Izz]);

    realised.GM_L              = steel_data.GM_realised;
    realised.mass_discrepancy  = realised.mass_total - realised.mass_buoyant_force;

    % Loud warnings on physically suspect realised values
    if realised.mass_total <= 0
        error('mwecmass:realise:RealisedZeroMass', ...
              'Realised mass_total = %g (<= 0).', realised.mass_total);
    end
    if realised.GM_L <= 0
        warning('mwecmass:realise:RealisedNegativeGM', ...
                ['Realised GM_L = %.4f m is non-positive — the steel-fill hull is ', ...
                 'STATICALLY UNSTABLE.  The optimiser-targeted GM may not have been ', ...
                 'achievable under the steel/air mass partition + t_min constraint.'], ...
                realised.GM_L);
    end
    if realised.Iyy <= 0
        warning('mwecmass:realise:RealisedZeroIyy', ...
                'Realised Iyy_about_cg = %g — pitch period will be undefined.', ...
                realised.Iyy);
    end

    % ── Hydrostatic stiffness ─────────────────────────────────
    K33_hydro = config.RHO_WATER * config.G * realised.Aw;
    % K55_hydro = M_total·g·GM (total mass at every site; equals rho_w·V_sub·g·GM at flotation balance).
    if realised.GM_L > 0
        K55_hydro = realised.mass_total * config.G * realised.GM_L;
    else
        K55_hydro = 0;
    end
    realised.K_hydro = diag([0, K33_hydro, K55_hydro]);

    % ── Added mass at the realised draft, referenced to the
    %    REALISED CG in a single call.  The delta congruence
    %    transform (formerly inline here) now lives inside
    %    mwecmass.bem.interpolate_at_draft.  See
    %    that function for the math.
    try
        [A11, A33, A55, A_full, B_full] = ...
            mwecmass.bem.interpolate_at_draft( ...
                steel_data.vertical_shift, config, realised.CG_total(3));
    catch ME
        error('mwecmass:realise:RealisedInterpolationFailed', ...
              ['mwecmass.bem.interpolate_at_draft failed for the realised draft ', ...
               'vs=%.4f m: %s. The hydro cache must cover the steel draft range.'], ...
               steel_data.vertical_shift, ME.message);
    end

    realised.A11    = A11;
    realised.A33    = A33;     % translation-invariant
    realised.A55    = A55;
    realised.A_full = A_full;
    realised.B_full = B_full;

    % ── PTO stiffness: the suite models a free-floating body, so K_pto is the constant zero matrix,
    % matching properties_2d.m/properties_3d.m, and K_total equals K_hydro exactly.
    realised.K_pto   = zeros(3, 3);
    realised.K_total = realised.K_hydro;

    % ── Uncoupled natural periods ────────────────────────────
    M11_virtual = realised.mass_total + A11;
    M33_virtual = realised.mass_total + A33;  %#ok<NASGU>
    M55_virtual = realised.Iyy        + A55;  %#ok<NASGU>

    % Free-floating body: no PTO stiffness contribution. K11_total is the surge hydrostatic
    % stiffness, which is zero by definition for a free-floating body, so K11_total = 0 ->
    % periods.surge = Inf always. Kept as an explicit
    % conditional (matching properties_2d.m/properties_3d.m's own style) rather than collapsed
    % to a bare `inf`, so a future re-introduction of surge stiffness only needs K11_total
    % redefined here.
    K11_total = 0;

    if K11_total > 1e-6
        realised.periods.surge = 2 * pi * sqrt(M11_virtual / K11_total);
    else
        realised.periods.surge = inf;
    end
    % Per-axis heave/pitch already computed by the steel solver — use those directly so we don't
    % re-derive and risk a tiny numerical drift. They are kept as the logged comparison; the
    % periods the objective reads come from the coupled eigenproblem below.
    realised.periods.heave_uncoupled = steel_data.T_heave_realised;   % [s] per-axis
    realised.periods.pitch_uncoupled = steel_data.T_pitch_realised;   % [s] per-axis

    % ── Coupled eigenproblem; heave and pitch named by share ─
    %  realised.periods.heave/.pitch are the eigen-periods of the coupled 3-DOF
    %  (1 surge, 2 heave, 3 pitch) undamped problem, each named by the DOF holding the largest
    %  share of that mode's kinetic energy. The stiffness-determinant guard that used to sit
    %  here -- false by construction, since K11_total = 0 makes that determinant exactly zero --
    %  and its 100*eye(3) placeholder are gone: the coupled problem is solved on every call.
    M_phys  = diag([realised.mass_total, realised.mass_total, realised.Iyy]);   % [kg, kg, kg m^2]
    M_total_3x3 = M_phys + A_full;   % [kg, kg m, kg m^2] physical + omega->inf added mass

    [T_coupled, share, Phi] = ...
        mwecmass.hydrostatics.coupled_periods_by_share(M_total_3x3, realised.K_total);

    realised.periods.heave = T_coupled.heave;   % [s] heave-dominated coupled mode
    realised.periods.pitch = T_coupled.pitch;   % [s] pitch-dominated coupled mode
    realised.participation_factors = share;     % [%] kinetic-energy share, rows surge/heave/pitch
    realised.coupled_periods = [T_coupled.surge; T_coupled.heave; T_coupled.pitch];   % [s]
    realised.coupled_modes = Phi;               % columns ordered surge, heave, pitch
    realised.surge_per_pitch = Phi(1,3) / Phi(3,3);   % [m/rad] surge per unit pitch in the pitch mode

    % ── 6×6 mass matrices (centre and origin) ────────────────
    realised.MassMatrix_CG = mwecmass.hydrostatics.build_mass_matrix( ...
        realised.mass_total, [0, 0, 0], realised.Inertia_Tensor);
    realised.MassMatrix_Origin = mwecmass.hydrostatics.build_mass_matrix( ...
        realised.mass_total, realised.CG_total, realised.Inertia_Tensor);

    % ── Realised per-strip mass partition ─────────────────────
    % If steel_data carries per-strip realisation arrays (UHPC path
    % via solve, or the modular-precast strip extraction in the steel
    % path), expose them on final_props.realised_strips so
    % diagnostics, visualisation and IO can read the AS-BUILT state
    % directly without recomputing.

    realised.fill_method             = 'steel_fill';   % overwritten below (line ~317) from the caller's own tag when supplied
    realised.density_profile_source = 'realised_partition';

    rs = struct();
    rs.fill_method = realised.fill_method;
    rs.t_offset = [];
    rs.is_solid = [];
    rs.is_wall  = [];
    rs.z_lo     = [];
    rs.z_hi     = [];
    rs.V_uhpc   = [];
    rs.V_void   = [];
    rs.mass_uhpc = [];
    rs.mass_void = [];
    rs.uhpc_volume_fraction = [];
    rs.contours_outer = {};
    rs.contours_inner = {};

    if isfield(steel_data, 't_offset_strip') && ~isempty(steel_data.t_offset_strip)
        rs.t_offset = steel_data.t_offset_strip(:);
    end
    if isfield(steel_data, 'is_solid_strip') && ~isempty(steel_data.is_solid_strip)
        rs.is_solid = steel_data.is_solid_strip(:);
    end
    if isfield(steel_data, 'wall_strip_idx') && ~isempty(steel_data.wall_strip_idx)
        N_rs = max([length(rs.t_offset), length(rs.is_solid), 0]);
        if N_rs > 0
            rs.is_wall = false(N_rs, 1);
            rs.is_wall(steel_data.wall_strip_idx) = true;
        end
    end
    if isfield(steel_data, 'strip_z_lo')
        rs.z_lo = steel_data.strip_z_lo(:);
    end
    if isfield(steel_data, 'strip_z_hi')
        rs.z_hi = steel_data.strip_z_hi(:);
    end
    if isfield(steel_data, 'strip_V_UHPC')
        rs.V_uhpc = steel_data.strip_V_UHPC(:);
    end
    if isfield(steel_data, 'strip_V_void')
        rs.V_void = steel_data.strip_V_void(:);
    end
    if isfield(steel_data, 'strip_mass_UHPC')
        rs.mass_uhpc = steel_data.strip_mass_UHPC(:);
    end
    if isfield(steel_data, 'strip_mass_void')
        rs.mass_void = steel_data.strip_mass_void(:);
    end
    if ~isempty(rs.V_uhpc) && ~isempty(rs.V_void)
        Vt = rs.V_uhpc + rs.V_void;
        Vt(Vt <= 0) = 1;
        rs.uhpc_volume_fraction = rs.V_uhpc ./ Vt;
    end
    if isfield(steel_data, 'contours_outer')
        rs.contours_outer = steel_data.contours_outer;
    end
    if isfield(steel_data, 'contours_inner')
        rs.contours_inner = steel_data.contours_inner;
    end
    realised.realised_strips = rs;

    % Backward-compat carry-overs (cross_section is the hull profile,
    % not material; densities_at_nodes is the optimiser's per-node rho —
    % we keep both for the figure functions).
    if isfield(final_props, 'cross_section')
        realised.cross_section = final_props.cross_section;
    else
        realised.cross_section = [];
    end
    if isfield(final_props, 'densities_at_nodes')
        realised.densities_at_nodes = final_props.densities_at_nodes;
    else
        realised.densities_at_nodes = [];
    end

    % Realised per-strip equivalent density (set by either pipeline):
    %   Steel:  steel_data.strip_rho_eff  (computed in solve)
    %   UHPC :  cstr.strip_rho_eff        (computed in
    %           mwecmass.realise.modular_precast.extract_strip_geometry)
    % Exposed at the top level of realised so visualisation
    % functions can render the AS-BUILT cake layers without
    % digging into solver-specific sub-structs.
    if isfield(steel_data, 'strip_rho_eff') && ~isempty(steel_data.strip_rho_eff)
        realised.realised_strip_density = steel_data.strip_rho_eff(:);
    else
        realised.realised_strip_density = [];
    end
    if isfield(steel_data, 'strip_edges') && ~isempty(steel_data.strip_edges)
        realised.realised_strip_edges = steel_data.strip_edges(:);
    elseif isfield(config, 'strip_edges')
        realised.realised_strip_edges = config.strip_edges(:);
    else
        realised.realised_strip_edges = [];
    end
    if isfield(steel_data, 'fill_method') && ~isempty(steel_data.fill_method)
        realised.fill_method = steel_data.fill_method;   % `steel_data` here is the generic second argument -- for the modular-precast caller it is actually `cstr` (extract_strip_geometry.m), whose fill_method is 'uhpc_fill'
    else
        realised.fill_method = '';
    end

    % `components`: REPLACE the optimiser's per-strip rho with the
    % realised UHPC/void partition where available, so any consumer
    % that reads final_props.components sees as-built data.
    if ~isempty(rs.mass_uhpc)
        N = length(rs.mass_uhpc);
        comps = repmat(struct('density', 0, 'z_level', 0, ...
                              'V_uhpc', 0, 'V_void', 0, ...
                              'mass_uhpc', 0, 'mass_void', 0, ...
                              'is_solid', false, 'is_wall', false, ...
                              't_offset', NaN), N, 1);
        for ii = 1:N
            if ~isempty(rs.z_lo) && ~isempty(rs.z_hi)
                z_mid = 0.5 * (rs.z_lo(ii) + rs.z_hi(ii));
            else
                z_mid = 0;
            end
            Vti = 0;
            if ~isempty(rs.V_uhpc), Vti = Vti + rs.V_uhpc(ii); end
            if ~isempty(rs.V_void), Vti = Vti + rs.V_void(ii); end
            Mti = rs.mass_uhpc(ii);
            if ~isempty(rs.mass_void), Mti = Mti + rs.mass_void(ii); end
            if Vti > 1e-12
                comps(ii).density = Mti / Vti;
            end
            comps(ii).z_level = z_mid;
            if ~isempty(rs.V_uhpc), comps(ii).V_uhpc = rs.V_uhpc(ii); end
            if ~isempty(rs.V_void), comps(ii).V_void = rs.V_void(ii); end
            comps(ii).mass_uhpc = rs.mass_uhpc(ii);
            if ~isempty(rs.mass_void), comps(ii).mass_void = rs.mass_void(ii); end
            if ~isempty(rs.is_solid), comps(ii).is_solid = rs.is_solid(ii); end
            if ~isempty(rs.is_wall),  comps(ii).is_wall  = rs.is_wall(ii);  end
            if ~isempty(rs.t_offset), comps(ii).t_offset = rs.t_offset(ii); end
        end
        realised.components = comps;
    elseif isfield(final_props, 'components')
        realised.components = final_props.components;
    else
        realised.components = struct('density', {}, 'z_level', {});
    end
end
