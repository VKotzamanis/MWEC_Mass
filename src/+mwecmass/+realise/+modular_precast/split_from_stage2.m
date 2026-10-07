function split = split_from_stage2(rho_stage2, V, rho_uhpc, rho_air, wall)
%SPLIT_FROM_STAGE2  UHPC and air volume of every module from its Stage-2 density (AGENTS 3 item 4.1).
%
%   split = mwecmass.realise.modular_precast.split_from_stage2(rho_stage2, V, rho_uhpc, rho_air, wall)
%
%   rho_stage2 [N x 1] Stage-2 bulk densities (x_opt(2:end)) [kg/m^3]; V [N x 1] module volumes of
%   the exact outer body [m^3]; rho_uhpc, rho_air [kg/m^3]; wall: index of the solid wall module
%   ([] when there is none).
%
%   V_uhpc_target = V (rho - rho_air) / (rho_uhpc - rho_air), V_air_target = V - V_uhpc_target, so
%   rho_uhpc V_uhpc_target + rho_air V_air_target = rho V for every module.
%   solid: the wall module and every module at rho_uhpc, the Stage-2 upper bound (full section).
%   k_star: the lowest module that is not solid, the ballast module ([] when every module is
%   solid). Every module below it is solid. hollow: the modules above k_star that are not solid;
%   each holds its V_uhpc_target in a shell of thickness t_i >= t_min around air.

rho_stage2 = rho_stage2(:);
V = V(:);
N = numel(V);
solid = rho_stage2 >= rho_uhpc;
solid(wall) = true;
k_star = find(~solid, 1);
hollow = false(N, 1);
hollow(k_star + 1:end) = ~solid(k_star + 1:end);
V_uhpc_target = V .* (rho_stage2 - rho_air) / (rho_uhpc - rho_air);
split = struct('rho_stage2', rho_stage2, 'V', V, 'V_uhpc_target', V_uhpc_target, ...
    'V_air_target', V - V_uhpc_target, 'solid', solid, 'k_star', k_star, 'hollow', hollow, ...
    'wall', wall);
end
