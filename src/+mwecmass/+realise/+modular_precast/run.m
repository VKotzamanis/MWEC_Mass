function [results, final_props] = run(config, x_opt, opt_results)
%RUN  Modular-precast Stage 3 (contract F14): realise the Stage-2 solution in UHPC modules.
%
%   [results, final_props] = mwecmass.realise.modular_precast.run(config, x_opt, opt_results)
%
%   x_opt = [vs; rho_1 ... rho_N] and opt_results.Final3D are the Stage-2 solution; it is split,
%   built and checked by solve_and_extract. results = opt_results with results.stage3 (S8) and
%   results.stage2_3d.properties = final_props; results.Final3D stays the Stage-2 design.
%   final_props always describes the realised design, accepted or failed. The figures
%   (plot_modular_precast) and the STEP files (mwecmass.output.step.export_stage3 into
%   Output/modular_precast/step) are produced for every status; config.output switches them off
%   (stage3.precast_midplane, stage3.precast_strips, stage3.step), and a config without output
%   options produces both.

fprintf('\n╔══════════════════════════════════════════════════╗\n');
fprintf('║ STAGE 3: MODULAR PRECAST (UHPC) REALISATION      ║\n');
fprintf('╚══════════════════════════════════════════════════╝\n');
[realised, final_props] = mwecmass.realise.modular_precast.solve_and_extract(config, x_opt, ...
    opt_results.Final3D);

has_out = isfield(config, 'output');
if ~has_out || config.output.save.stage3.precast_midplane || config.output.save.stage3.precast_strips
    mwecmass.output.figures.plot_modular_precast(realised, config);
end
if ~has_out || config.output.save.stage3.step
    realised.step_files = mwecmass.output.step.export_stage3(realised, ...
        mwecmass.output.output_dir('modular_precast'));
end

results = opt_results;
results.stage3 = realised;
results.stage2_3d.properties = final_props;
end
