function [results, final_props] = load_results(file)
%LOAD_RESULTS Load the results/final_props pair from a mass-suite export .mat file.
    L = load(file, 'results', 'final_props');
    results = L.results;
    final_props = L.final_props;
    if isfield(results, 'mode') && ~isfield(results, 'realisation_type')
        results.realisation_type = results.mode;
    end
    if isfield(final_props, 'realisation_mode') && ~isfield(final_props, 'fill_method')
        final_props.fill_method = final_props.realisation_mode;
    end
    if isfield(results, 'steel_data') && isstruct(results.steel_data) ...
            && isfield(results.steel_data, 'realisation_mode') && ~isfield(results.steel_data, 'fill_method')
        results.steel_data.fill_method = results.steel_data.realisation_mode;
    end
end
