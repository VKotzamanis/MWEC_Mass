function profile = extract_midplane_profile(parser)
%EXTRACT_MIDPLANE_PROFILE Collect the hull boundary profile near y=0.
%   Returns [N x 2] [x,z] coordinates in m.

    try
        n_edge_pts = 200;
        t_edge     = linspace(0, 1, n_edge_pts)';
        all_xz     = zeros(0, 2);

        for s = 1:length(parser.visible_surfs)
            sname = parser.visible_surfs{s};
            e = parser.entities(sname);

            switch e.type
                case 'RevSurf'
                    % Profile curve IS the y=0 silhouette (at start angle).
                    pts = parser.eval_curve(e.params.profile, t_edge);
                    all_xz = [all_xz; pts(:, [1, 3])]; %#ok<AGROW>

                case 'BLoftSurf'
                    % Evaluate edges at u=0 and u=1 (boundary curves)
                    for ie = [1, 3]  % edges 1 and 3: v=0,1 with u varying
                        for it = 1:n_edge_pts
                            if ie == 1
                                pt = parser.eval_surface(sname, t_edge(it), 0);
                            else
                                pt = parser.eval_surface(sname, t_edge(it), 1);
                            end
                            if abs(pt(2)) < 0.05  % near y=0
                                all_xz(end+1, :) = [pt(1), pt(3)]; %#ok<AGROW>
                            end
                        end
                    end
                    % Also check edges at v=0 and v=1
                    for ie = [2, 4]
                        for it = 1:n_edge_pts
                            if ie == 2
                                pt = parser.eval_surface(sname, 1, t_edge(it));
                            else
                                pt = parser.eval_surface(sname, 0, t_edge(it));
                            end
                            if abs(pt(2)) < 0.05
                                all_xz(end+1, :) = [pt(1), pt(3)]; %#ok<AGROW>
                            end
                        end
                    end

                case 'RuledSurf'
                    pts1 = parser.eval_curve(e.params.curve1, t_edge);
                    pts2 = parser.eval_curve(e.params.curve2, t_edge);
                    near1 = abs(pts1(:,2)) < 0.05;
                    near2 = abs(pts2(:,2)) < 0.05;
                    all_xz = [all_xz; pts1(near1, [1,3]); pts2(near2, [1,3])]; %#ok<AGROW>

                case 'DevSurf'
                    pts_s = parser.eval_snake(e.params.snake, t_edge);
                    pts_c = parser.eval_curve(e.params.curve, t_edge);
                    near_s = abs(pts_s(:,2)) < 0.05;
                    near_c = abs(pts_c(:,2)) < 0.05;
                    all_xz = [all_xz; pts_s(near_s, [1,3]); pts_c(near_c, [1,3])]; %#ok<AGROW>

                case 'MirrSurf'
                    % MirrSurf about Y=0: the source boundary AT y=0
                    % is the profile.  Already captured by the source
                    % surface's edge evaluation.  Skip to avoid
                    % duplicating points.
                    continue;
            end
        end

        if size(all_xz, 1) < 3
            warning('mwecmass:geometry:ProfileEmpty', ...
                    'Profile extraction found < 3 points');
            profile = all_xz;
            return;
        end

        % Mirror across x = 0 to get the full symmetric profile.
        % The source surfaces (RevSurf, DevSurf, etc.) produce
        % x >= 0 points only.  MirrSurf(Y=0) creates the x < 0
        % half of the 3D hull, but doesn't appear in the y=0
        % silhouette because MirrSurf flips y, not x.  For the
        % XZ profile view, bilateral symmetry means the left half
        % is a mirror of the right half about x = 0.
        mirrored = all_xz;
        mirrored(:, 1) = -mirrored(:, 1);
        % Remove mirrored points that are at x ≈ 0 (would duplicate)
        not_on_axis = abs(mirrored(:, 1)) > 1e-6;
        all_xz = [all_xz; mirrored(not_on_axis, :)];

        % Remove duplicates
        all_xz = unique(round(all_xz * 1e6) / 1e6, 'rows', 'stable');

        % Sort by angle from centroid to form proper polygon
        cx = mean(all_xz(:, 1));
        cz = mean(all_xz(:, 2));
        angles = atan2(all_xz(:, 2) - cz, all_xz(:, 1) - cx);
        [~, order] = sort(angles);
        profile = all_xz(order, :);

    catch ME
        % Fail loud; a silent fallback profile would propagate into downstream GM/waterplane calculations.
        newME = MException('mwecmass:geometry:ProfileFailed', ...
            'extract_midplane_profile: profile extraction failed for this hull deck: %s', ...
            ME.message);
        newME = addCause(newME, ME);
        throw(newME);
    end
end
