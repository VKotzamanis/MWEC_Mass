function result = stage1_screen_draft(vs_k, config)
%STAGE1_SCREEN_DRAFT Probe one selected draft without fmincon.
% Builds a flotation-balanced density vector, pinning the modular-precast wall strip when enabled,
% then evaluates 3-D properties once. Tier 1 checks only mass balance and GM; full constraints belong
% to Tier 2. The result contains the candidate, objective value, properties, and feasibility flag.
    N        = config.num_ballast_sections;
    V_strips = sum(config.strip_V);   % total hull strip volume from config

    % Probe with initial densities to get V_sub at this draft
    x_probe     = [vs_k, config.initial_densities];
    props_probe = mwecmass.hydrostatics.properties_3d(x_probe, config);
    V_sub       = props_probe.V_sub;

    if V_sub < 1e-6 || V_strips < 1e-6
        result.feasible = false;
        result.fval     = Inf;
        result.x        = x_probe;
        result.props    = props_probe;
        return
    end

    % Build mass-balanced density vector
    if config.enable_constructability && ~isempty(config.wall_strip_index)
        % Wall strip pinned; distribute remaining required mass to platform
        w         = config.wall_strip_index;
        mass_wall = config.strip_V(w) * config.constructability_rho_hull;
        V_free    = V_strips - config.strip_V(w);
        mass_need = config.RHO_WATER * V_sub - mass_wall;
        rho_free  = mass_need / max(V_free, eps);
        rho_free  = max(config.ballast_density_bounds(1), ...
                    min(config.ballast_density_bounds(2), rho_free));

        rho_bal      = repmat(rho_free, 1, N);
        rho_bal(w)   = config.constructability_rho_hull;
    else
        rho_unif = config.RHO_WATER * V_sub / V_strips;
        rho_unif = max(config.ballast_density_bounds(1), ...
                   min(config.ballast_density_bounds(2), rho_unif));
        rho_bal  = repmat(rho_unif, 1, N);
    end

    x_bal  = [vs_k, rho_bal];
    props  = mwecmass.hydrostatics.properties_3d(x_bal, config);
    fval   = mwecmass.optim.stage2_objective(x_bal, config);

    mass_err    = abs(props.mass_total - props.mass_buoyant_force) / max(props.mass_total, 1);
    is_feasible = mass_err < 0.05 && props.GM_L >= config.gm_min;

    result.feasible = is_feasible;
    result.fval     = fval;
    result.x        = x_bal;
    result.props    = props;
end
