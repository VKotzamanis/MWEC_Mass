function action = hams_enrichment_action(config)
%HAMS_ENRICHMENT_ACTION Select live HAMS-MREL enrichment versus cached interpolation.
% config.run_HAMS_MREL must be true and config.hams_dir non-empty to return 'run'.
% Otherwise return 'skip'; off-node hydrostatics remain available through cache interpolation.
% See docs/METHODS_ENGINE.md#optim-pid-correction
    if config.run_HAMS_MREL && ~isempty(config.hams_dir)
        action = 'run';
    else
        action = 'skip';
    end
end
