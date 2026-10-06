function cmap = cividis_map(n)
%CIVIDIS_MAP Return an n-row Cividis colormap array.
% cmap = cividis_map(n) returns n-by-3 RGB values via interpolation on
% nine control points. Cividis is a blue-yellow colormap, approximately
% colorblind-safe.
    ctrl = [0.000, 0.135, 0.305; 0.000, 0.205, 0.380; 0.123, 0.263, 0.406;
            0.253, 0.318, 0.420; 0.365, 0.373, 0.432; 0.479, 0.435, 0.440;
            0.608, 0.518, 0.430; 0.770, 0.640, 0.380; 0.995, 0.906, 0.144];
    cmap = flipud(interp1(linspace(0, 1, 9), ctrl, linspace(0, 1, n)));
end
