function cache = precompute_boundary_cache(parser, n_u)
%PRECOMPUTE_BOUNDARY_CACHE Cache source-surface boundaries for iso-z queries.
%   n_u sets the parameter grid; cached coordinates and axis data are in m.

    topo = parser.classify_visible_surfaces();
    cache.sources   = topo.sources;
    cache.mirrors   = topo.mirrors;
    cache.u_samples = linspace(0, 1, n_u)';
    cache.data      = containers.Map();

    for s = 1:length(topo.sources)
        sname = topo.sources{s};
        e = parser.entities(sname);
        d = struct();
        d.type = e.type;

        switch e.type
            case {'DevSurf', 'RuledSurf'}
                if strcmp(e.type, 'DevSurf')
                    d.pts_1 = parser.eval_snake(e.params.snake, cache.u_samples);
                    d.pts_2 = parser.eval_curve(e.params.curve, cache.u_samples);
                else
                    d.pts_1 = parser.eval_curve(e.params.curve1, cache.u_samples);
                    d.pts_2 = parser.eval_curve(e.params.curve2, cache.u_samples);
                end
                d.z_boundary_1 = d.pts_1(:, 3);
                d.z_boundary_2 = d.pts_2(:, 3);

            case 'RevSurf'
                profile_pts = parser.eval_curve_or_snake( ...
                                  e.params.profile, cache.u_samples);
                d.profile_pts = profile_pts;
                d.z_profile   = profile_pts(:, 3);
                d.profile_name = e.params.profile;

                axis_ent = parser.entities(e.params.axis);
                d.axis_start = parser.eval_any_point(axis_ent.params.pt_start);
                d.axis_end   = parser.eval_any_point(axis_ent.params.pt_end);
                d.axis_vec   = d.axis_end - d.axis_start;
                d.axis_dir   = d.axis_vec / norm(d.axis_vec);
                d.angle_start = e.params.angle_start;
                d.angle_end   = e.params.angle_end;

            case 'BLoftSurf'
                sec_names = e.params.section_names;
                n_sec = length(sec_names);
                d.degree = e.params.degree;
                d.n_sections = n_sec;
                d.section_pts = zeros(n_sec, n_u, 3);
                d.z_sections  = zeros(n_sec, n_u);

                for k = 1:n_sec
                    pts_k = parser.eval_curve_or_snake( ...
                                sec_names{k}, cache.u_samples);
                    d.section_pts(k, :, :) = pts_k;
                    d.z_sections(k, :) = pts_k(:, 3)';
                end

                d.knots_v = mwecmass.geometry.MS2Parser.make_clamped_knots(n_sec, d.degree);
        end

        cache.data(sname) = d;
    end
end
