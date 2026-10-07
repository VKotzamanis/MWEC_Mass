function bp = body_properties(body, rho, opts)
%BODY_PROPERTIES  Volume, first and second moments and mass properties of a B-rep body (contract F6, S6).
%
%   bp = mwecmass.solid.body_properties(body, rho, opts)
%
%   body: S4 (build_body). rho: densities by region [kg/m^3]: uhpc, air (modular_precast) or
%   ballast, shell, air (thin_shell); for any other body.design.mode the regions are the fields
%   of rho (hydrostatics_at_draft integrates the submerged outer pieces that way). opts.n_gauss: Gauss-Legendre points per knot span in u and
%   in v (default 8: exact for every face up to bicubic, whose integrands have degree 14 in each
%   parameter; convergent for rational faces).
%
%   Divergence theorem with fields (f, 0, 0), so only the x component of the outward normal
%   enters: V = surf x n_x, int x = surf x^2/2 n_x, int y = surf x y n_x, int z = surf x z n_x,
%   int x^2 = surf x^3/3 n_x, int y^2 = surf x y^2 n_x, int z^2 = surf x z^2 n_x, int x y =
%   surf x^2 y/2 n_x, int x z = surf x^2 z/2 n_x, int y z = surf x y z n_x, with n dA = S_u x S_v
%   du dv on each face (face normal: surface normal, reversed if not same_sense). Plane faces and
%   constant-z faces have n_x = 0 and add exactly 0, so only lateral faces are integrated. Region
%   r of module i sums the faces of module i with inside = r and subtracts those with outside = r.
%
%   bp.regions.<r>: V [N x 1] (m^3), S [N x 3] = int x dV (m^4), J [3 x 3 x N] = int x x' dV
%   (m^5); bp.modules(i): V, V_<r>, mass, rho_eff = mass/V, CG_body; bp.total: mass, CG_body,
%   I_origin, I_cg (I = int rho (|x|^2 E - x x') dV, body frame); bp.quad = n_gauss.

if nargin < 3 || isempty(opts)
    opts = struct();
end
ng = 8;
if isfield(opts, 'n_gauss') && ~isempty(opts.n_gauss)
    ng = opts.n_gauss;
end
switch body.design.mode
    case 'modular_precast'
        names = {'uhpc', 'air'};
    case 'thin_shell'
        names = {'ballast', 'shell', 'air'};
    otherwise
        names = fieldnames(rho)';
end
N = numel(body.design.edges) - 1;
for r = 1:numel(names)
    if ~isfield(rho, names{r})
        error('mwecmass:solid:MissingDensity', 'body_properties: rho.%s is required', names{r});
    end
    reg.(names{r}) = zeros(N, 10);
end
B = body.brep;
for k = 1:numel(B.faces)
    f = B.faces(k);
    s = B.surfaces{f.surface};
    if ~strcmp(s.type, 'bspline') || all(reshape(s.ctrl(:, :, 3), [], 1) == s.ctrl(1, 1, 3))
        continue
    end
    m = face_moments(s, ng);
    if ~f.same_sense
        m = -m;
    end
    i = f.module(1);
    if isfield(reg, f.inside)
        reg.(f.inside)(i, :) = reg.(f.inside)(i, :) + m;
    end
    if isfield(reg, f.outside)
        reg.(f.outside)(i, :) = reg.(f.outside)(i, :) - m;
    end
end

modules = struct('V', cell(1, N));
M = 0;
Sm = zeros(1, 3);
I0 = zeros(3);
for r = 1:numel(names)
    R = reg.(names{r});
    J = zeros(3, 3, N);
    for i = 1:N
        J(:, :, i) = sym3(R(i, 5:10));
    end
    bp.regions.(names{r}) = struct('V', R(:, 1), 'S', R(:, 2:4), 'J', J);
end
for i = 1:N
    V = 0;
    mass = 0;
    S = zeros(1, 3);
    for r = 1:numel(names)
        R = bp.regions.(names{r});
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
bp.quad = ng;
end

function J = sym3(q)
% [int x^2, int y^2, int z^2, int xy, int xz, int yz] as a symmetric matrix
J = [q(1) q(4) q(5); q(4) q(2) q(6); q(5) q(6) q(3)];
end

function m = face_moments(s, ng)
% [V, int x, int y, int z, int x^2, int y^2, int z^2, int xy, int xz, int yz] of one face
[xg, wg] = gauss_legendre(ng);
ku = unique(s.knots{1});
kv = unique(s.knots{2});
au = ku(1:end - 1);
bu = ku(2:end);
av = kv(1:end - 1);
bv = kv(2:end);
uu = (au(:)' + bu(:)') / 2 + (bu(:)' - au(:)') / 2 .* xg;
vv = (av(:)' + bv(:)') / 2 + (bv(:)' - av(:)') / 2 .* xg;
wu = (bu(:)' - au(:)') / 2 .* wg;
wv = (bv(:)' - av(:)') / 2 .* wg;
[U, V] = ndgrid(uu(:), vv(:));
[Wu, Wv] = ndgrid(wu(:), wv(:));
[P, Su, Sv] = mwecmass.solid.eval_bspline_surface(s, U(:), V(:));
nx = Su(:, 2) .* Sv(:, 3) - Su(:, 3) .* Sv(:, 2);
w = Wu(:) .* Wv(:) .* nx;
x = P(:, 1);
y = P(:, 2);
z = P(:, 3);
F = [x, x.^2 / 2, x .* y, x .* z, x.^3 / 3, x .* y.^2, x .* z.^2, x.^2 .* y / 2, x.^2 .* z / 2, x .* y .* z];
m = w' * F;
end

function [x, w] = gauss_legendre(n)
% Golub-Welsch: nodes and weights of n-point Gauss-Legendre on [-1, 1] (columns)
persistent cache
if isempty(cache)
    cache = cell(1, 64);
end
if ~isempty(cache{n})
    x = cache{n}(:, 1);
    w = cache{n}(:, 2);
    return
end
b = (1:n - 1) ./ sqrt(4 * (1:n - 1).^2 - 1);
[Vv, D] = eig(diag(b, 1) + diag(b, -1));
[x, i] = sort(diag(D));
w = 2 * Vv(1, i)'.^2;
cache{n} = [x w];
end
