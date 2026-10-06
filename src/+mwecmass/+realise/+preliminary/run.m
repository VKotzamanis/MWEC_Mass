function [results, final_props] = run(~, ~, opt_results)
%RUN Pass through Stage-2 results for the preliminary (no-realisation) mode.
% config and x_opt are accepted for dispatch parity but unused. results equals opt_results and
% final_props equals opt_results.Final3D.
    % opt_results already holds the state this type returns (see header). Nothing to add.
    results     = opt_results;
    final_props = opt_results.Final3D;   % the un-realised Stage-2 result
end
