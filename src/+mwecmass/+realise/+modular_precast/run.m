function [results, final_props] = run(config, x_opt, opt_results)
%RUN Execute modular-precast UHPC/void realisation after Stage 2.
% x_opt is the Stage-2 design vector and opt_results contains un-realised Final3D properties.
% The solve/extraction result updates results.constructability and results.stage2_3d.properties;
% results.Final3D remains the optimiser result. Outputs are results and realised final_props.
    % Reuse the optimiser design and pre-realisation properties supplied by the caller.
    x_opt_3d = x_opt;
    final_props_optimiser = opt_results.Final3D;
    final_props = final_props_optimiser;

    if config.enable_constructability
        fprintf('\n╔══════════════════════════════════════════════════╗\n');
        fprintf('║ CONSTRUCTABILITY POST-PROCESSING                 ║\n');
        fprintf('╚══════════════════════════════════════════════════╝\n\n');

        constructability = mwecmass.realise.modular_precast.solve_and_extract( ...
            config, x_opt_3d, final_props_optimiser);
        % The constructability view writes two images, the midplane elevation and the per-strip
        % plan view; the call is made when either is wanted. A configuration carrying no
        % output options -- a caller that invokes this realisation directly rather than through
        % mwecmass.driver.run -- keeps the unconditional draw.
        if ~isfield(config, 'output') || config.output.save.stage3.precast_midplane || ...
                config.output.save.stage3.precast_strips
            mwecmass.output.figures.plot_modular_precast(constructability, config);
        end
        final_props = mwecmass.realise.modular_precast.build_realised_properties( ...
            final_props_optimiser, constructability, config);

        fprintf('\n  Constructability realization complete.\n');
    else
        constructability = [];
    end

    % Preserve Final3D and replace only constructability and realised Stage-2 properties.
    results = opt_results;
    results.constructability = constructability;
    results.stage2_3d.properties = final_props;
end
