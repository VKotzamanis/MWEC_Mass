function bp = body_properties(body, rho, opts) %#ok<INUSD>
%BODY_PROPERTIES  Stand-in of contract F6: S6 mass properties of a stand-in body in closed form.
%
%   bp = mwecmass.solid.body_properties(body, rho, opts)
%
%   body: S4 from the build_body stand-in (body.analytic set, else mwecmass:standin:NotAnalytic).
%   rho: densities by region [kg/m^3] (uhpc, air for modular precast; ballast, shell, air for thin
%   shell). V, S = int x dV, J = int x x' dV per region and module come from sti_closed_form
%   'regions' (prisms; source in its help); I = sum over regions of rho (trace(J) E - J), and
%   I_cg = I_origin - M (|c|^2 E - c c') (parallel-axis theorem). Body frame, about the origin.
%   quad is empty: no quadrature is used.

if ~isstruct(body) || ~isfield(body, 'analytic') || isempty(body.analytic)
    error('mwecmass:standin:NotAnalytic', 'body_properties stand-in: body.analytic is empty (not a stand-in body)');
end
if exist('sti_closed_form', 'file') ~= 2
    addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'fixtures'), '-end');
end
a = body.analytic;
reg = sti_closed_form('regions', a.fixture, body.design, a.d);
names = setdiff(fieldnames(reg), {'V_module'}, 'stable');
N = numel(body.design.edges) - 1;
modules = struct('V', cell(1, N));
M = 0;
Sm = zeros(1, 3);
I0 = zeros(3);
for r = 1:numel(names)
    bp.regions.(names{r}) = reg.(names{r});
end
for i = 1:N
    V = 0;
    mass = 0;
    S = zeros(1, 3);
    for r = 1:numel(names)
        R = reg.(names{r});
        modules(i).(['V_' names{r}]) = R.V(i);
        V = V + R.V(i);
        mass = mass + rho.(names{r}) * R.V(i);
        S = S + rho.(names{r}) * R.S(i, :);
        J = R.J(:, :, i);
        I0 = I0 + rho.(names{r}) * (trace(J) * eye(3) - J);
    end
    modules(i).V = V;
    modules(i).mass = mass;
    modules(i).rho_eff = mass / V;
    modules(i).CG_body = S / mass;
    M = M + mass;
    Sm = Sm + S;
end
c = Sm / M;
bp.modules = modules;
bp.total = struct('mass', M, 'CG_body', c, 'I_origin', I0, 'I_cg', I0 - M * (dot(c, c) * eye(3) - c' * c));
bp.quad = [];
end
