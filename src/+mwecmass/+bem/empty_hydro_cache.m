function cache = empty_hydro_cache(config)
%EMPTY_HYDRO_CACHE  Create an empty hydro_cache/hydro_table struct in the field's terms.
    %
    %   cache = mwecmass.bem.empty_hydro_cache()
    %   cache = mwecmass.bem.empty_hydro_cache(config)

    % config is optional; allows populating solver_params.depth when water_depth is present.
    if nargin < 1, config = struct(); end

    cache.drafts = []; % [m], body-frame vertical shift, signed z positive up.
    cache.z_cg = []; % [m], coefficient reference-point z, positive up.
    cache.submerged_volume = []; % [m^3], positive displaced volume.
    cache.displaced_mass = []; % [kg], positive Archimedean displaced mass.
    cache.omega = []; % [rad/s], positive circular wave frequency.
    cache.added_mass_inf = {}; % [kg, kg m, kg m^2], body-origin reference.
    cache.radiation_damping_band_avg = {}; % [N s/m and rotational equivalents], body-origin reference.
    cache.added_mass_omega = {}; % [kg, kg m, kg m^2], body-origin reference.
    cache.radiation_damping_omega = {}; % [N s/m and rotational equivalents], body-origin reference.
    cache.exciting_force_omega = {}; % [N and N m], body-origin reference.
    cache.period_band = [4.0, 16.0]; % [s], positive periods used for damping average.
    cache.solver_params = mwecmass.bem.hams_mrel.default_hams_params(config);
    cache.water_depth = NaN; % [m], positive down, unknown until a solver run supplies it.
    cache.ulen = NaN; % [m], positive reference length, unknown until a solver run supplies it.
    cache.timestamp = datestr(now); %#ok<DATST,TNOW1> -- cache schema
    cache.ms2_file = '';
    cache.ms2_date = '';
end
