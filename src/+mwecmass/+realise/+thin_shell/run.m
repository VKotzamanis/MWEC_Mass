function [results, final_props] = run(config, x_opt, opt_results)
%RUN Perform thin-shell material realisation after Stage 2.
%   [results,final_props] = run(config,x_opt,opt_results) solves for steel
%   thickness/fill when enabled, rebuilds realised properties when feasible,
%   and stores steel_data plus realised stage2_3d.properties. x_opt is
%   [1+N x 1] ([vertical shift; densities]); infeasible solves retain the
%   optimiser properties as a documented fallback.

    % Preserve the optimiser record for the infeasible-solve fallback.
    x_opt_3d = x_opt;
    final_props_optimiser = opt_results.Final3D;
    final_props = final_props_optimiser;
    steel_data  = [];

    if isfield(config, 'enable_steel_solve') && config.enable_steel_solve
        % The steel-solve console text is written to Output/thin_shell/steel_solve.log whenever
        % config.output (WEC_Output_Options) is present and out.save.stage3.steel_solve_log is
        % true, echoed to stdout too when out.console_echo is true. A caller that
        % invokes this realisation directly with no config.output keeps the unconditional
        % console-only print (solve.m defaults fids to 1), matching the figure gate immediately
        % below.
        if isfield(config, 'output') && config.output.save.stage3.steel_solve_log
            out = config.output;
            fid_steel = mwecmass.output.open_log(out, 'thin_shell', 'steel_solve');
            % onCleanup ensures cleanup on all exit paths; the helper is idempotent against
            % the explicit close_log call (fopen returns '' for closed files).
            steel_log_cleanup = onCleanup(@() close_fid_if_open(fid_steel));
            if out.console_echo
                fids_steel = [1 fid_steel];
            else
                fids_steel = fid_steel;
            end
            steel_data = mwecmass.realise.thin_shell.solve(config, x_opt_3d, final_props_optimiser, [], fids_steel);
            mwecmass.output.close_log(fid_steel);
        else
            steel_data = mwecmass.realise.thin_shell.solve(config, x_opt_3d, final_props_optimiser);
        end
        % The steel-fill diagnostic writes every configured export format when it is wanted. A
        % configuration carrying no output options -- a caller that
        % invokes this realisation directly rather than through mwecmass.driver.run -- keeps the
        % unconditional draw.
        if ~isfield(config, 'output') || config.output.save.stage3.steel_solve
            mwecmass.output.figures.plot_steel_solve(config, steel_data);
        end
        if isfield(steel_data, 'feasible') && steel_data.feasible
            final_props = mwecmass.realise.build_realised_properties( ...
                final_props_optimiser, steel_data, config);
        else
            warning('mwecmass:thin_shell:SteelInfeasibleNoSwap', ...
                'Steel solve infeasible — final_props NOT updated, falling back to optimiser.');
        end
    end

    % results = opt_results: opt_results already holds every result field except the two this
    %   realisation type changes (results.steel_data and results.stage2_3d.properties, both
    %   overwritten immediately below). results.Final3D is left untouched: it must stay at the
    %   pre-realisation value, which is what the optimisation step put there.
    results = opt_results;
    results.steel_data = steel_data;
    results.stage2_3d.properties = final_props;
end

function close_fid_if_open(fid)
%CLOSE_FID_IF_OPEN Close a file identifier only if it is still open; safe against prior close_log calls.
    if ~isempty(fid) && fid > 0 && ~isempty(fopen(fid))
        mwecmass.output.close_log(fid);
    end
end
