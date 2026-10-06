function clipped = clip_profile_to_z_range(profile, z_lo, z_hi)
%CLIP_PROFILE_TO_Z_RANGE Clip a 2D profile polygon to a [z_lo, z_hi] band.
% profile is n-by-2 [x, z] vertices; z_lo, z_hi in [m]. Returns clipped
% polygon or empty array if fewer than 3 vertices remain.
    poly = profile;
    poly = mwecmass.output.figures.clip_to_halfplane(poly, z_lo, +1);
    poly = mwecmass.output.figures.clip_to_halfplane(poly, z_hi, -1);
    if size(poly, 1) < 3
        clipped = [];
    else
        clipped = poly;
    end
end
