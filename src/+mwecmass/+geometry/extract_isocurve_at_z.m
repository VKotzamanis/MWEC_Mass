function pts = extract_isocurve_at_z(parser, z_wl, n_u, cache)
%EXTRACT_ISOCURVE_AT_Z Compute a cached visible-surface iso-z contour.
%   z_wl and returned [N x 3] points are in m; n_u matches cache.u_samples.

    pts = [];
    source_pts = containers.Map();

    for s = 1:length(cache.sources)
        sname = cache.sources{s};
        d = cache.data(sname);

        switch d.type
            case {'DevSurf', 'RuledSurf'}
                src_pts = mwecmass.geometry.isocurve_devsurf( ...
                              z_wl, d.pts_1, d.pts_2);

            case 'RevSurf'
                src_pts = mwecmass.geometry.isocurve_revsurf( ...
                              parser, z_wl, n_u, d);

            case 'BLoftSurf'
                src_pts = mwecmass.geometry.isocurve_bloftsurf( ...
                              z_wl, n_u, d);
            otherwise
                src_pts = [];
        end

        source_pts(sname) = src_pts;
        if ~isempty(src_pts)
            pts = [pts; src_pts]; %#ok<AGROW>
        end
    end

    % Mirror contributions: coordinate flip on source contour
    for m = 1:length(cache.mirrors)
        mirr = cache.mirrors(m);
        ult = mirr.ultimate_source;

        if source_pts.isKey(ult) && ~isempty(source_pts(ult))
            mirr_pts = source_pts(ult);
            for fi = 1:length(mirr.effective_flips)
                if strcmp(mirr.effective_flips{fi}, 'Y')
                    mirr_pts(:, 2) = -mirr_pts(:, 2);
                else
                    mirr_pts(:, 1) = -mirr_pts(:, 1);
                end
            end
            pts = [pts; mirr_pts]; %#ok<AGROW>
        end
    end
end
