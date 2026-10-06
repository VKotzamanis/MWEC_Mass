function [open_edges, all_edges] = detect_open_edges(parser, n_samples, tol)
%DETECT_OPEN_EDGES Find visible-surface boundary edges not shared by another surface.
%   Samples are [N x 3] coordinates in m; n_samples and tol control matching.

    if nargin < 2 || isempty(n_samples), n_samples = 50; end
    if nargin < 3 || isempty(tol),       tol = 1e-4;     end

    t_sample   = linspace(0, 1, n_samples);
    surf_names = parser.visible_surfs;
    n_surfs    = length(surf_names);

    % ── Sample all 4 edges of each visible surface ──────────
    all_edges = {};
    for s = 1:n_surfs
        sname = surf_names{s};
        for edge_idx = 1:4
            pts = zeros(n_samples, 3);
            for k = 1:n_samples
                t = t_sample(k);
                switch edge_idx
                    case 1  % v = 0, u varies
                        pts(k,:) = parser.eval_surface(sname, t, 0);
                    case 2  % u = 1, v varies
                        pts(k,:) = parser.eval_surface(sname, 1, t);
                    case 3  % v = 1, u varies
                        pts(k,:) = parser.eval_surface(sname, t, 1);
                    case 4  % u = 0, v varies
                        pts(k,:) = parser.eval_surface(sname, 0, t);
                end
            end

            edge = struct();
            edge.surface  = sname;
            edge.edge_idx = edge_idx;
            edge.pts      = pts;
            edge.z_mean   = mean(pts(:, 3));
            edge.matched  = false;
            all_edges{end+1} = edge; %#ok<AGROW>
        end
    end

    % ── Match edges pairwise ────────────────────────────────
    %  Two edges match if their sampled points coincide (up to
    %  tolerance) in either forward or reversed order.
    n_edges = length(all_edges);

    for a = 1:n_edges
        if all_edges{a}.matched, continue; end
        for b = a+1:n_edges
            if all_edges{b}.matched, continue; end
            % Don't match edges from the same surface
            if strcmp(all_edges{a}.surface, all_edges{b}.surface)
                continue;
            end

            pts_a = all_edges{a}.pts;
            pts_b = all_edges{b}.pts;

            dist_fwd = max(vecnorm(pts_a - pts_b, 2, 2));
            dist_rev = max(vecnorm(pts_a - flipud(pts_b), 2, 2));

            if min(dist_fwd, dist_rev) < tol
                all_edges{a}.matched = true; %#ok<AGROW> -- existing elements
                all_edges{b}.matched = true; %#ok<AGROW> -- existing elements
                break;  % edge a is matched, move on
            end
        end
    end

    % ── Collect unmatched ───────────────────────────────────
    open_edges = {};
    for i = 1:n_edges
        if ~all_edges{i}.matched
            open_edges{end+1} = all_edges{i}; %#ok<AGROW>
        end
    end
end
