function c = hams_constants()
%HAMS_CONSTANTS  Return constants used to dimensionalise HAMS output.
%   c.RHO_WATER [kg/m^3] and c.G [m/s^2] match HAMS's internal normalisation;
%   keep this density separate from physical seawater density in hydrostatics.
%   See docs/METHODS_ENGINE.md#bem-hams-normalization.
    c.RHO_WATER = 1000.0;   % [kg/m^3] HAMS internal rho -- DO NOT change to 1025
    c.G         = 9.80665;  % [m/s^2]  HAMS internal g -- matches WavDynMods.f90
end
