function [negate_x, negate_y, quadrant_str] = detect_source_quadrant(parser, sources)
%DETECT_SOURCE_QUADRANT Infer coordinate reflections needed to place source surfaces in Q1.
% Syntax: [negate_x,negate_y,quadrant_str] = detect_source_quadrant(parser,sources).
% Inputs: parser geometry parser; sources cell array of surface names.
% Outputs: logical reflection flags and 'Q1'...'Q4', based on sampled mean x/y signs.

    negate_x = false;
    negate_y = false;
    quadrant_str = 'Q1';

    x_sum = 0;
    y_sum = 0;
    n_pts = 0;
    % Mean sign unreliable for axis-spanning sources; track extrema to flag them.
    x_min = Inf; x_max = -Inf;
    y_min = Inf; y_max = -Inf;

    u_probe = [0.25, 0.5, 0.75];
    v_probe = [0.25, 0.5, 0.75];

    for s = 1:length(sources)
        for ui = 1:length(u_probe)
            for vi = 1:length(v_probe)
                try
                    pt = parser.eval_surface(sources{s}, u_probe(ui), v_probe(vi));
                    x_sum = x_sum + pt(1);
                    y_sum = y_sum + pt(2);
                    x_min = min(x_min, pt(1)); x_max = max(x_max, pt(1));
                    y_min = min(y_min, pt(2)); y_max = max(y_max, pt(2));
                    n_pts = n_pts + 1;
                catch
                    % skip evaluation failures
                end
            end
        end
    end

    if n_pts == 0
        % Fail loud; silently defaulting to Q1 previously masked parser failures.
        error('mwecmass:mesh:SourceQuadrantUndetermined', ...
              ['Could not evaluate any source surface points (%d source(s), %d probes each) -- ' ...
               'quadrant is undetermined.'], length(sources), numel(u_probe) * numel(v_probe));
    end

    if (x_min < 0 && x_max > 0) || (y_min < 0 && y_max > 0)
        warning('mwecmass:mesh:SourceQuadrantAmbiguous', ...
                ['Source surface samples span x=0 or y=0; the mean-sign quadrant ' ...
                 'classification below may be unreliable for this source.']);
    end

    x_avg = x_sum / n_pts;
    y_avg = y_sum / n_pts;

    if x_avg < 0, negate_x = true; end
    if y_avg < 0, negate_y = true; end

    if     negate_x && negate_y, quadrant_str = 'Q3';
    elseif negate_x,             quadrant_str = 'Q2';
    elseif negate_y,             quadrant_str = 'Q4';
    end
end
