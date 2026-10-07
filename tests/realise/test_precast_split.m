function test_precast_split()
%TEST_PRECAST_SPLIT  Split of the Stage-2 densities into UHPC and air (AGENTS section 3 item 4.1).
%   No kernel: hand-made volumes and densities.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'src'));
f = @mwecmass.realise.modular_precast.split_from_stage2;
rho_u = 2500;
rho_a = 1.2;
V = [3.1; 2.7; 4.4; 1.9; 5.3];

% C1-like: module 1 at the upper bound, wall module on top
s = f([2500; 1430; 369.4; 369.4; 2500], V, rho_u, rho_a, 5);
check(isequal(s.solid', [true false false false true]) && s.k_star == 2 && ...
    isequal(s.hollow', [false false true true false]), 'C1-like roles');
rho = [2500; 1430; 369.4; 369.4; 2500];
% two products and a sum per module: within 4 rounding steps of rho V
lhs = rho_u * s.V_uhpc_target + rho_a * s.V_air_target;
check(all(abs(lhs - rho .* V) <= 4 * eps * rho_u * V), 'mass of the split');
check(isequal(s.V_uhpc_target, V .* (rho - rho_a) / (rho_u - rho_a)) && ...
    isequal(s.V_air_target, V - s.V_uhpc_target), 'split formula');
check(all(s.V_uhpc_target(s.solid & rho == rho_u) == V(s.solid & rho == rho_u)), 'solid modules all UHPC');
fprintf('C1-like split: V_uhpc %s, V_air %s, mass error %s kg\n', mat2str(s.V_uhpc_target', 6), ...
    mat2str(s.V_air_target', 6), mat2str((lhs - rho .* V)', 3));

% a module above k* at the upper bound stays solid; modules below k* are solid
s = f([2500; 1200; 2500; 400; 2500], V, rho_u, rho_a, 5);
check(isequal(find(s.solid)', [1 3 5]) && s.k_star == 2 && isequal(find(s.hollow)', 4), 'non-monotone roles');

% wall at the bottom, no other solid module
s = f([2500; 900; 700; 500; 300], V, rho_u, rho_a, 1);
check(isequal(find(s.solid)', 1) && s.k_star == 2 && isequal(find(s.hollow)', [3 4 5]), 'bottom wall');

% no wall module
s = f([1500; 900; 700; 500; 300], V, rho_u, rho_a, []);
check(~any(s.solid) && s.k_star == 1 && isequal(find(s.hollow)', 2:5) && isempty(s.wall), 'no wall');

% every module solid: no ballast module
s = f(rho_u * ones(5, 1), V, rho_u, rho_a, 5);
check(all(s.solid) && isempty(s.k_star) && ~any(s.hollow), 'all solid');
end

function check(cond, msg)
if ~cond
    error('test_precast_split:fail', '%s', msg);
end
end
