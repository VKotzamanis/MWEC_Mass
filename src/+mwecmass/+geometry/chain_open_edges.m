function groups = chain_open_edges(open_edges, tol)
%CHAIN_OPEN_EDGES Group and chain open boundary edges into loops.
%   Edge points are [N x 3] in m; groups contain loop points and mean elevation.

    if nargin < 2, tol = 1e-3; end

    n = length(open_edges);
    if n == 0
        groups = {};
        return;
    end

    % ── Group edges by z-level ──────────────────────────────
    z_vals   = cellfun(@(e) e.z_mean, open_edges);
    assigned = false(n, 1);
    z_groups = {};

    for i = 1:n
        if assigned(i), continue; end
        members = i;
        assigned(i) = true;
        for j = i+1:n
            if ~assigned(j) && abs(z_vals(j) - z_vals(i)) < 0.05
                members(end+1) = j; %#ok<AGROW>
                assigned(j) = true;
            end
        end
        z_groups{end+1} = members; %#ok<AGROW>
    end

    % ── Chain each group into a loop ────────────────────────
    groups = {};
    for g = 1:length(z_groups)
        idx   = z_groups{g};
        edges = open_edges(idx);
        n_e   = length(edges);

        % Greedy endpoint-matching chain
        used     = false(n_e, 1);
        order    = zeros(n_e, 1);
        rev_flag = false(n_e, 1);

        order(1) = 1;
        used(1)  = true;

        for k = 2:n_e
            prev = order(k-1);
            if rev_flag(prev)
                current_end = edges{prev}.pts(1, :);
            else
                current_end = edges{prev}.pts(end, :);
            end

            best_dist = inf;
            best_idx  = 0;
            best_rev  = false;

            for j = 1:n_e
                if used(j), continue; end
                d_fwd = norm(current_end - edges{j}.pts(1, :));
                d_rev = norm(current_end - edges{j}.pts(end, :));

                if d_fwd < best_dist
                    best_dist = d_fwd;
                    best_idx  = j;
                    best_rev  = false;
                end
                if d_rev < best_dist
                    best_dist = d_rev;
                    best_idx  = j;
                    best_rev  = true;
                end
            end

            if best_dist > tol
                % Raise error if no edge within tol to avoid silently corrupting cap area with spurious bridge.
                error('mwecmass:geometry:chain_open_edges:NoCandidateWithinTolerance', ...
                    ['chain_open_edges: at z~=%.4f, could not extend the open-edge chain ' ...
                     'within tol=%.4g m after %d of %d edges (nearest available gap = %.4g m). ' ...
                     'The z-group likely contains more than one disjoint boundary loop.'], ...
                    mean(z_vals(idx)), tol, k-1, n_e, best_dist);
            end

            order(k)    = best_idx;
            used(best_idx) = true;
            rev_flag(k) = best_rev;
        end

        % Concatenate points in chain order
        loop = [];
        for k = 1:n_e
            pts = edges{order(k)}.pts;
            if rev_flag(k)
                pts = flipud(pts);
            end
            if k > 1
                pts = pts(2:end, :);  % avoid duplicate junction point
            end
            loop = [loop; pts]; %#ok<AGROW>
        end

        % Unclosed loop gives wrong cap area; raise error rather than silently return invalid polyline.
        if size(loop,1) > 1
            if norm(loop(1,:) - loop(end,:)) < tol
                loop = loop(1:end-1, :);
            else
                error('mwecmass:geometry:chain_open_edges:UnclosedLoop', ...
                    ['chain_open_edges: at z~=%.4f, the %d chained edge(s) do not close into ' ...
                     'a loop (endpoint gap = %.4g m > tol = %.4g m).'], ...
                    mean(z_vals(idx)), n_e, norm(loop(1,:) - loop(end,:)), tol);
            end
        end

        grp = struct();
        grp.loop_pts = loop;
        grp.z_mean   = mean(loop(:, 3));
        grp.z_range  = max(loop(:,3)) - min(loop(:,3));
        groups{end+1} = grp; %#ok<AGROW>
    end
end
