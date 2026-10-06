function [V_shell, A_outer, debug] = compute_perpendicular_shell_volume( ...
                                          parser, z_lo, z_hi, t, n_z, b_cache)
%COMPUTE_PERPENDICULAR_SHELL_VOLUME Compute volume and outer lateral area
% of an axisymmetric shell with inward normal thickness t.
% The shell spans [z_lo,z_hi]; the optional debug output exposes samples.
% The shell is an inward normal offset of thickness t over [z_lo,z_hi].
% The optional debug output exposes the sampled profile and offset radius.
%   FORMULA (leading order in t)
%     For profile radius r(z), the inward radial offset is
%     t·sqrt(1+r'(z)^2):
%
%         r_inner(z) = max(0, r(z) − t·sqrt(1 + r'(z)²))
%
%     V_shell = π · ∫_{z_lo}^{z_hi} (r²(z) − r_inner²(z)) dz
%     A_outer = ∫_{z_lo}^{z_hi} 2π·r(z)·sqrt(1+r'(z)²) dz   (lateral)
%
%     Clamping r_inner at zero makes a slice solid when the offset
%     exceeds its available radial radius.
%
%   The profile formula is a leading-order normal-offset approximation and
%   avoids second derivatives of r(z).
%
%   AXISYMMETRIC ASSUMPTION
%     compute_rmin_at_z supplies the inscribed-circle radius. For a
%     non-axisymmetric cross-section, the resulting volume is a lower bound.
%
%   INPUTS
%     parser   — mwecmass.geometry.MS2Parser object
%     z_lo, z_hi — strip z bounds [m] (z_lo < z_hi)
%     t        — perpendicular shell thickness [m] (t ≥ 0)
%     n_z      — number of z-samples (uses 100 if []).
%     b_cache  — boundary cache from precompute_boundary_cache
%
%   OUTPUTS
%     V_shell  — shell volume [m³]
%     A_outer  — outer lateral surface area [m²]
%     debug    — struct with z_samples, r_samples, rp_samples,
%                L_samples, r_inner_samples, solid_fraction
%                (fraction of samples with r_inner = 0)

    if z_hi <= z_lo
        error('mwecmass:hydrostatics:BadStripBounds', ...
              'z_hi (%.4f) <= z_lo (%.4f)', z_hi, z_lo);
    end
    if t < 0
        error('mwecmass:hydrostatics:NegativeThickness', ...
              't (%.4f) < 0', t);
    end
    if nargin < 5 || isempty(n_z), n_z = 100; end
    if ~isscalar(n_z) || n_z ~= round(n_z) || n_z < 2
        error('mwecmass:hydrostatics:BadNumZSamples', ...
            'compute_perpendicular_shell_volume: n_z (%s) must be an integer >= 2 -- the central-difference slope estimate needs at least 2 z-samples.', ...
            mat2str(n_z));
    end
    if nargin < 6 || isempty(b_cache)
        b_cache = mwecmass.geometry.precompute_boundary_cache(parser, 100);
    end

    z_samples = linspace(z_lo, z_hi, n_z)';
    r_samples = zeros(n_z, 1);

    for k = 1:n_z
        [r_k, ~] = mwecmass.geometry.compute_rmin_at_z( ...
                        parser, z_samples(k), 100, b_cache);
        r_samples(k) = r_k;
    end

    % Central differences provide the profile slope at each sampled z.
    rp_samples = zeros(n_z, 1);
    for k = 1:n_z
        if k == 1
            dz = z_samples(2) - z_samples(1);
            rp_samples(k) = (r_samples(2) - r_samples(1)) / dz;
        elseif k == n_z
            dz = z_samples(n_z) - z_samples(n_z-1);
            rp_samples(k) = (r_samples(n_z) - r_samples(n_z-1)) / dz;
        else
            dz = z_samples(k+1) - z_samples(k-1);
            rp_samples(k) = (r_samples(k+1) - r_samples(k-1)) / dz;
        end
    end
    % Dome-tip correction: zero rp at samples where r < t (the
    % surface has effectively reached its closing point).  Mirrors
    % realisation's own strip extraction.
    rp_samples(r_samples < t) = 0;

    L_samples       = sqrt(1 + rp_samples.^2);
    r_inner_samples = max(0, r_samples - t .* L_samples);

    % Volume integrand: π · (r² - r_inner²)
    integrand_V = pi .* (r_samples.^2 - r_inner_samples.^2);
    V_shell     = trapz(z_samples, integrand_V);

    % Lateral surface area (no caps)
    integrand_A = 2*pi .* r_samples .* L_samples;
    A_outer     = trapz(z_samples, integrand_A);

    if nargout >= 3
        debug.z_samples       = z_samples;
        debug.r_samples       = r_samples;
        debug.rp_samples      = rp_samples;
        debug.L_samples       = L_samples;
        debug.r_inner_samples = r_inner_samples;
        debug.solid_fraction  = sum(r_inner_samples == 0) / n_z;
    end
end
