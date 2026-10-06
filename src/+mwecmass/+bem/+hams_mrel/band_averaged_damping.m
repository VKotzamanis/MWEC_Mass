function B_avg = band_averaged_damping(hams_data, T_min, T_max)
%BAND_AVERAGED_DAMPING Average radiation damping over a period band.
% hams_data supplies omega [rad/s] and B [6x6xN]; T_min and T_max are [s].
% Returns a 6x6 mean over the inclusive band, or the nearest sample if empty.
% See docs/METHODS_ENGINE.md#bem-frequency-band.

    omega = hams_data.omega;
    T = 2 * pi ./ omega;
    mask = (T >= T_min) & (T <= T_max);

    if ~any(mask)
        T_center = 0.5 * (T_min + T_max);
        [~, closest] = min(abs(T - T_center));
        mask(closest) = true;
    end

    B_avg = mean(hams_data.B(:, :, mask), 3);
end
