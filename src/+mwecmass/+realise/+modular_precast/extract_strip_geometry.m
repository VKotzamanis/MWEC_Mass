function modules = extract_strip_geometry(bp, design, rho_stage2, rho_floor)
%EXTRACT_STRIP_GEOMETRY  Per-module record of a realised body (S8 modules).
%
%   modules = mwecmass.realise.modular_precast.extract_strip_geometry(bp, design, rho_stage2, rho_floor)
%
%   bp: S6 of the realised body; design: S3; rho_stage2, rho_floor [N x 1] the Stage-2 density and
%   the Stage-2 density floor of every module [kg/m^3] (NaN where none is set).
%   modules(i): z_lo, z_hi (module edges, body frame), t (shell thickness; NaN when the module has
%   no air), h_ballast (ballast height from the module bottom, OD6 option ii), V and V_<region>
%   for every region of bp, mass, rho_eff = mass / V, rho_stage2, rho_floor, CG_world =
%   CG_body + [0 0 vs].

e = design.edges(:);
N = numel(e) - 1;
names = fieldnames(bp.regions)';
modules = struct('z_lo', num2cell(e(1:end - 1)), 'z_hi', num2cell(e(2:end)));
for i = 1:N
    m = bp.modules(i);
    modules(i).t = NaN;
    if m.V_air > 0
        modules(i).t = design.t(i);
    end
    modules(i).h_ballast = min(max(design.z_ballast - e(i), 0), e(i + 1) - e(i));
    modules(i).V = m.V;
    for r = 1:numel(names)
        modules(i).(['V_' names{r}]) = m.(['V_' names{r}]);
    end
    modules(i).mass = m.mass;
    modules(i).rho_eff = m.rho_eff;
    modules(i).rho_stage2 = rho_stage2(i);
    modules(i).rho_floor = rho_floor(i);
    modules(i).CG_world = m.CG_body + [0 0 design.vs];
end
end
