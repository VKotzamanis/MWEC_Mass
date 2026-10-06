function [verts, panels, surf_ids, is_cap, wl_verts] = ...
        trim_at_waterline(verts, panels, surf_ids, is_cap, z_wl)
%TRIM_AT_WATERLINE Keep the submerged mesh and clip panels crossing z = z_wl.
% Syntax: [verts,panels,surf_ids,is_cap,wl_verts] = trim_at_waterline(...,z_wl).
% Inputs: mesh arrays (vertices in m, panels as indices, ids/cap flags) and waterline z [m].
% Outputs retain array conventions; wl_verts lists generated waterline vertices. A 1e-10 m band
% classifies panels already on the waterline as submerged; straddling panels are clipped linearly.

    n_panels    = size(panels, 1);
    new_verts   = verts;
    new_panels  = zeros(0, 4);
    new_sids    = zeros(0, 1);
    new_caps    = false(0, 1);
    wl_verts    = [];

    for p = 1:n_panels
        vidx = panels(p, :);
        z    = verts(vidx, 3);

        if all(z <= z_wl + 1e-10)
            % Fully submerged — keep
            new_panels(end+1, :) = vidx; %#ok<AGROW>
            new_sids(end+1, 1)   = surf_ids(p); %#ok<AGROW>
            new_caps(end+1, 1)   = is_cap(p); %#ok<AGROW>

        elseif all(z > z_wl - 1e-10)
            % Fully above water — discard
            continue;

        else
            % Straddling — split at waterline
            [sub_p, sub_v, sub_wl] = ...
                mwecmass.mesh.split_panel_at_z(verts, vidx, z_wl);
            if ~isempty(sub_p)
                vo = size(new_verts, 1);
                new_verts = [new_verts; sub_v]; %#ok<AGROW>
                for sp = 1:size(sub_p, 1)
                    new_panels(end+1, :) = sub_p(sp, :) + vo; %#ok<AGROW>
                    new_sids(end+1, 1)   = surf_ids(p); %#ok<AGROW>
                    new_caps(end+1, 1)   = is_cap(p); %#ok<AGROW>
                end
                wl_verts = [wl_verts; sub_wl + vo]; %#ok<AGROW>
            end
        end
    end

    verts    = new_verts;
    panels   = new_panels;
    surf_ids = new_sids;
    is_cap   = new_caps;
end
