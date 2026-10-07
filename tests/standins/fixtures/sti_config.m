function config = sti_config(name)
%STI_CONFIG  config fields the Stage-3 code reads, for a stand-in fixture ('cylinder' | 'box').
%
%   config = sti_config(name)
%
%   RHO_WATER, G, mass_acceptable_pct: the author inputs of WEC_User_Input.m (1025 kg/m^3,
%   9.80665 m/s^2, 10 %). strip_edges: fixture module edges (body frame). profile: the y = 0
%   section of the fixture, [x z] counter-clockwise, from its closed form. ms2_model,
%   boundary_cache, hull_z_min, hull_z_max as build_config sets them. One-draft hydro table:
%   hydro_drafts, hydro_z_cg, added_mass_diagonal, added_mass_full, radiation_damping_full are
%   test data chosen symmetric positive definite (with a surge-pitch coupling A15), not BEM
%   results. config.hull_solid is not set: tests put the S1b there (contract section 3).

fx = sti_closed_form('fixture', name);
config = struct();
config.RHO_WATER = 1025;
config.G = 9.80665;
config.mass_acceptable_pct = 10;
evalc('config.ms2_model = mwecmass.geometry.MS2Parser.parse(fx.deck);');
config.boundary_cache = mwecmass.geometry.precompute_boundary_cache(config.ms2_model, 21);
config.hull_z_min = fx.z(1);
config.hull_z_max = fx.z(2);
if strcmp(fx.kind, 'cylinder')
    config.strip_edges = [-3; -2; -1; 0; 1];
    x = [-fx.R fx.R];
    A = [16000 0 -6000; 0 4000 0; -6000 0 12000];
    B = diag([500 800 300]);
    z_cg = -1.6;
else
    config.strip_edges = [-2.5; -1.5; -0.5; 0.5];
    x = fx.x;
    A = [9000 0 -2000; 0 6000 0; -2000 0 4000];
    B = diag([300 600 200]);
    z_cg = -1.2;
end
config.profile = [x(1) fx.z(1); x(2) fx.z(1); x(2) fx.z(2); x(1) fx.z(2)];
config.hydro_drafts = 0.5;
config.hydro_z_cg = z_cg;
config.added_mass_diagonal = diag(A)';
config.added_mass_full = {A};
config.radiation_damping_full = {B};
end
