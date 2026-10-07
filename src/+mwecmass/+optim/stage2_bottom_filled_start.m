function x0 = stage2_bottom_filled_start(vertical_shift, config, lb, ub)
%STAGE2_BOTTOM_FILLED_START Stage-2 start with the mass in the lowest modules.
% Every module starts at its lower density bound (lb(2:end)); modules are then raised to their
% upper bound from the keel up (module 1 is the lowest) until the mass equals the displaced mass
% rho_w * V_sub at vertical_shift, the last one partly. A pinned module (lb = ub, the UHPC wall
% module) keeps its density and is skipped. If the lower bounds already carry more than the
% displaced mass, or even all modules at their upper bound fall short of it, the vector is returned
% as far as it got and flotation is left to the solver. x0 is [vertical_shift, rho_1 ... rho_N].
    rho = lb(2:end);
    rho = rho(:)';
    props = mwecmass.hydrostatics.properties_3d([vertical_shift, rho], config);
    deficit = props.mass_buoyant_force - props.mass_total;
    for i = 1:numel(rho)
        if deficit <= 0
            break;
        end
        room = config.strip_V(i) * (ub(1 + i) - lb(1 + i));
        if room <= 0
            continue;
        end
        if deficit >= room
            rho(i) = ub(1 + i);
            deficit = deficit - room;
        else
            rho(i) = lb(1 + i) + deficit / config.strip_V(i);
            deficit = 0;
        end
    end
    x0 = [vertical_shift, rho];
end
