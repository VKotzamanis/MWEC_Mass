function phi = range_penalty(r, k_amp)
%RANGE_PENALTY C1 piecewise-quadratic penalty for a normalized target error.
% r is dimensionless (actual-target)/half_range; k_amp is dimensionless curvature outside |r|=1.
% Inside the band phi=r^2. Outside, delta=abs(r)-1 and phi=1+2*delta+k_amp*delta^2.
% See docs/METHODS_ENGINE.md#optim-stage2-formulation
if r < -1
    delta = -r - 1;     % positive
    phi = 1 + 2*delta + k_amp * delta^2;
elseif r > 1
    delta = r - 1;      % positive
    phi = 1 + 2*delta + k_amp * delta^2;
else
    phi = r^2;
end
end
