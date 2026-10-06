function [xc, zc] = clip_z(x, z, z_cut, side)
%CLIP_Z Clip polygon vertices against the horizontal line z = z_cut.
% side='below' retains z <= z_cut; side='above' retains z >= z_cut.
% Inputs: x,z matching vertex vectors; z_cut scalar; side character vector.
% Outputs: xc,zc clipped vertices, or empty arrays when fewer than three inputs remain.
    keep_below = strcmp(side, 'below');
    n = length(x);
    if n < 3,  xc = [];  zc = [];  return;  end
    xc = zeros(2*n, 1);
    zc = zeros(2*n, 1);
    cnt = 0;
    for i = 1:n
        j  = mod(i, n) + 1;
        zi = z(i);  zj = z(j);
        if keep_below
            in_i = (zi <= z_cut);
            in_j = (zj <= z_cut);
        else
            in_i = (zi >= z_cut);
            in_j = (zj >= z_cut);
        end
        if in_i
            cnt = cnt + 1;
            xc(cnt) = x(i);  zc(cnt) = zi;
        end
        if in_i ~= in_j && abs(zj - zi) > 1e-14
            t    = (z_cut - zi) / (zj - zi);
            cnt  = cnt + 1;
            xc(cnt) = x(i) + t * (x(j) - x(i));
            zc(cnt) = z_cut;
        end
    end
    xc = xc(1:cnt);
    zc = zc(1:cnt);
end
