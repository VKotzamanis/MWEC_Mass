function rho_eq = compute_strip_equivalent_density(config, densities_core)
%COMPUTE_STRIP_EQUIVALENT_DENSITY Compute equivalent density per strip from core densities.
% densities_core is N-by-1; uses config.shell_density, config.shell.V_shell,
% config.shell.V_core to return volume-weighted equivalent rho_eq [kg/m^3].
    N       = length(densities_core);
    rho_eq  = densities_core(:);
    if isempty(config.shell), return; end
    V_shell = config.shell.V_shell;
    V_core  = config.shell.V_core;
    V_total = V_shell + V_core;
    for i = 1:N
        if V_total(i) > 1e-12
            rho_eq(i) = (config.shell_density * V_shell(i) + ...
                         densities_core(i)    * V_core(i)) / V_total(i);
        end
    end
end
