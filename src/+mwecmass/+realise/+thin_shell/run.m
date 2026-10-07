function [results, final_props] = run(config, x_opt, opt_results)
%RUN  Thin-shell Stage 3 (contract F14): realise the Stage-2 design on the exact geometry.
%
%   [results, final_props] = mwecmass.realise.thin_shell.run(config, x_opt, opt_results)
%
%   x_opt = [vs; rho_1..N] and opt_results.Final3D are the Stage-2 solution. The realised design,
%   or the closest design when Stage 3 fails (status 'failed', per-metric report), is stored as
%   results.stage3 (S8), drawn (mwecmass.output.figures.plot_steel_solve) and exported as STEP
%   (mwecmass.output.step.export_stage3), whatever its status; final_props describes it and
%   never the Stage-2 properties. results.Final3D keeps the Stage-2 properties.
%   Switches (config.output, when present): out.save.stage3.steel_solve_log (report log),
%   out.save.stage3.steel_solve (figure), out.save.stage3.step (STEP files); a caller without
%   config.output gets the report on the console, the figure and the STEP files.

has_output = isfield(config, 'output');
fids = 1;
if has_output && config.output.save.stage3.steel_solve_log
    out = config.output;
    fid = mwecmass.output.open_log(out, 'thin_shell', 'steel_solve');
    log_cleanup = onCleanup(@() close_fid_if_open(fid));
    if out.console_echo
        fids = [1 fid];
    else
        fids = fid;
    end
end
realised = mwecmass.realise.thin_shell.solve(config, x_opt, opt_results.Final3D, fids);
clear log_cleanup
if ~has_output || config.output.save.stage3.steel_solve
    mwecmass.output.figures.plot_steel_solve(realised, config);
end
if ~has_output || config.output.save.stage3.step
    realised.step_files = mwecmass.output.step.export_stage3(realised, ...
        mwecmass.output.output_dir('thin_shell'));
end

final_props = realised.props;
final_props.stage3_status = realised.status;
final_props.stage3_check = realised.check;
results = opt_results;
results.stage3 = realised;
results.stage2_3d.properties = final_props;
end

function close_fid_if_open(fid)
if ~isempty(fid) && fid > 0 && ~isempty(fopen(fid))
    mwecmass.output.close_log(fid);
end
end
