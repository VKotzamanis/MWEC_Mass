function fl = density_floors(model, cache, geo, edges, t_min, rho_solid, rho_air, solid_modules, opts)
%DENSITY_FLOORS Lowest buildable density of each module: a t_min shell around air, no ballast.
%   fl = DENSITY_FLOORS(model, cache, geo, edges, t_min, rho_solid, rho_air, solid_modules)
%   fl = DENSITY_FLOORS(..., opts)
%
%   model: MS2Parser model of the deck, cache: its boundary cache, geo: the hull's outer patches
%   (mwecmass.solid.outer_nurbs), edges [N+1 x 1] module edges (body frame, m), t_min the shell
%   thickness [m], rho_solid and rho_air [kg/m^3], solid_modules the indices of the modules that
%   are solid throughout (the UHPC wall module; [] for thin shell).
%
%   The hollow modules get the shell of design thickness t_min that Stage 3 builds for t_i = t_min
%   (inner surface = outer surface offset along the exact normal by t_min + eps_fit/2, eps_fit =
%   0.01 t_min, fitted by mwecmass.solid.offset_surface), assembled into one body by
%   mwecmass.solid.build_body with z_ballast at the keel (no ballast); region volumes come from
%   mwecmass.solid.body_properties. With V_solid and V_air the volumes of the shell and of the
%   void in a module and V = V_solid + V_air,
%       rho_min = (rho_solid V_solid + rho_air V_air) / V.
%   A solid module has rho_min = rho_solid, V_air = 0.
%
%   opts (optional): mode ('modular_precast' | 'thin_shell'; default 'modular_precast' when
%   solid_modules is not empty, else 'thin_shell'), n_gauss (forwarded to body_properties) and
%   max_passes (forwarded to offset_surface). The first and last edge must equal the hull's z
%   range to 16 ulp (the bound of the outer rows); they are replaced by the exact z range.
%
%   fl: rho_min, V, V_solid, V_air [N x 1] and fit, the fit report (S2r) of the inner set.

    if nargin < 9 || isempty(opts), opts = struct(); end
    edges = edges(:);
    N = numel(edges) - 1;
    solid_modules = solid_modules(:)';
    if isfield(opts, 'mode')
        mode = opts.mode;
    elseif isempty(solid_modules)
        mode = 'thin_shell';
    else
        mode = 'modular_precast';
    end
    switch mode
        case 'modular_precast'
            solid_region = 'uhpc';
        case 'thin_shell'
            solid_region = 'shell';
        otherwise
            error('mwecmass:driver:BadMode', 'density_floors: mode ''%s'' is neither modular_precast nor thin_shell.', mode);
    end
    hollow = setdiff(1:N, solid_modules);
    if isempty(hollow)
        error('mwecmass:driver:NoHollowModule', 'density_floors: every module is solid, there is no shell to size.');
    end

    z_hull = geo.z_range;
    ends = edges([1 end])';
    if any(abs(ends - z_hull) > 16 * eps(max(1, abs(z_hull))))
        error('mwecmass:driver:EdgesNotOnHull', ...
              'density_floors: outer module edges [%.17g %.17g] m differ from the hull z extent [%.17g %.17g] m.', ...
              ends, z_hull);
    end
    edges([1 end]) = z_hull;

    t = repmat(t_min, N, 1);
    t(solid_modules) = NaN;
    fit_opts = struct('t_min', t_min);
    if isfield(opts, 'max_passes'), fit_opts.max_passes = opts.max_passes; end
    z_hollow = [edges(min(hollow)), edges(max(hollow) + 1)];
    [inner, fit] = mwecmass.solid.offset_surface(model, cache, geo, t_min, z_hollow, fit_opts);

    design = struct('mode', mode, 'edges', edges, 'vs', 0, 't', t, ...
                    'z_ballast', z_hull(1), 'solid_modules', solid_modules);
    body = mwecmass.solid.build_body(geo, design, inner);

    rho = struct(solid_region, rho_solid, 'air', rho_air);
    if strcmp(mode, 'thin_shell')
        rho.ballast = rho_solid;   % the body has no ballast: its volume is 0 in every module
    end
    bp_opts = struct();
    if isfield(opts, 'n_gauss'), bp_opts.n_gauss = opts.n_gauss; end
    bp = mwecmass.solid.body_properties(body, rho, bp_opts);

    V = zeros(N, 1);
    V_solid = zeros(N, 1);
    V_air = zeros(N, 1);
    for i = 1:N
        V(i) = bp.modules(i).V;
        V_solid(i) = bp.modules(i).(['V_' solid_region]);
        V_air(i) = bp.modules(i).V_air;
    end
    fl = struct('rho_min', (rho_solid * V_solid + rho_air * V_air) ./ V, ...
                'V', V, 'V_solid', V_solid, 'V_air', V_air, 'fit', fit);
end
