% run_eta_power_stiffness_compatible.m
%   Stage-1 joint placement + damping + geometry optimisation under
%   eta_power, with FULL passive-feasible stiffness compatibility from
%   the symmetric-tether architecture.
%
%   ============================================================
%   ARCHITECTURE & ASSUMPTIONS (defend in peer review)
%   ============================================================
%   The two-tether PTO contributes BOTH stiffness and damping with the
%   same matrix multipliers tied to (alpha, ell):
%       K_PTO_matrix = 2 K_PTO [cos^2 a   0    -ell cos a;
%                                   0   sin^2 a    0     ;
%                               -ell cos a 0    ell^2    ]
%       B_PTO_matrix = (same structure, B_PTO instead of K_PTO)
%   Surge has K_hyd,11 = 0 but an INDEPENDENT mooring spring K_moor >= 0
%   provides the surge restoring force (consistent with the existing
%   script's '[Surge mooring]' headline).
%
%   Therefore the three stiffness-matching equations are:
%       (1) 2 K_PTO cos^2 a  + K_moor = omega_n,1^2 (M + A_11(omega_n,1))    ≡ S_1
%       (2) 2 K_PTO sin^2 a            = omega_n,2^2 (M + A_22(omega_n,2)) - K_hyd,22 ≡ S_2
%       (3) 2 K_PTO ell^2              = omega_n,3^2 (Iyy + A_33(omega_n,3)) - K_hyd,33 ≡ S_3
%
%   Solving:
%       K_PTO  = S_2 / (2 sin^2 a)
%       ell    = sin(a) sqrt(S_3 / S_2)
%       K_moor = S_1 - S_2 cot^2 a
%
%   Free design variables (5): omega_n,1, omega_n,2, omega_n,3, alpha, B_PTO.
%
%   Passive constraints:
%     S_2 >= 0       -> omega_n,2 >= omega_n_hydro,h  (heave T_n <= hydrostatic)
%     S_3 >= 0       -> omega_n,3 >= omega_n_hydro,p  (pitch T_n <= hydrostatic)
%     K_PTO >= 0     (implied by S_2 >= 0 with sin^2 a > 0)
%     K_moor >= 0    -> S_2 cot^2 a <= S_1            (couples alpha, surge, heave)
%     B_PTO >= 1
%     alpha in [30 deg, 60 deg]
%     ell <= b_t,max sin(a) + h_t,max cos(a)  (hull)
%
%   Objective: same eta_power as run_eta_power_joint.m
%     eta_power = sum_k nu_k * 4 beta_k/(1+beta_k)^2 * integral Lorentzian / denom
%   Optimiser: fmincon SQP with 16 multistart restarts.

clear; clc;
addpath(pwd);

%% PHASE 0 -- baseline from script
climate = fullfile(pwd, 'WIS_Output_WAM', 'ST84040_climate_grid.mat');
props   = fullfile(pwd, 'WEC_Mass_UHPC.mat');
opts    = struct('plot', false, 'verbose', false);

fprintf('================================================================\n');
fprintf('  STIFFNESS-COMPATIBLE PASSIVE ETA_POWER OPTIMISATION\n');
fprintf('================================================================\n');
res = run_MWEC_tuning(climate, props, opts);

bem     = res.bem;
closure = res.closure;
M       = closure.M;
Iyy     = closure.Iyy;
M_diag  = [M; M; Iyy];
nu      = [2; 1; 2];

omega_grid = res.climate.omega(:);
S_ew       = res.climate.S_ew(:);
omega_L    = res.climate.omega_L;
omega_H    = res.climate.omega_H;
mask       = (omega_grid >= omega_L) & (omega_grid <= omega_H);
SF         = S_ew ./ max(omega_grid, 1e-9).^3;
denom_int  = sum(nu) * trapz(omega_grid(mask), SF(mask));

% Hydrostatic K and Tn (passive feasibility bounds)
fp_S = load(props);  fp = fp_S.final_props;
K22_hyd = fp.K_hydro(2,2);
K33_hyd = fp.K_hydro(3,3);
om_hydro_h = 2*pi / bem.T_n_hydrostatic(2);
om_hydro_p = 2*pi / bem.T_n_hydrostatic(3);

% Hull
cs = fp.cross_section;
bt_max = max(abs(cs(:,1)));
ht_max = max(cs(:,2));
ht_min = min(cs(:,2));

% Bundle for objective
data.bem = bem; data.M_diag = M_diag; data.nu = nu;
data.omega = omega_grid; data.SF = SF; data.mask = mask; data.denom = denom_int;
data.M = M; data.Iyy = Iyy;
data.K22_hyd = K22_hyd; data.K33_hyd = K33_hyd;
data.bt_max = bt_max;  data.ht_max = ht_max;

fprintf('\n  Hydrostatic omega_n: heave %.4f r/s (T %.3f s),  pitch %.4f r/s (T %.3f s)\n', ...
    om_hydro_h, 2*pi/om_hydro_h, om_hydro_p, 2*pi/om_hydro_p);
fprintf('  90%% band: T in [%.2f, %.2f] s\n', 2*pi/omega_H, 2*pi/omega_L);
fprintf('  T_p = %.3f s\n', 2*pi/res.climate.IEC.omega_p);
fprintf('  K_hyd_22 = %.4g N/m,  K_hyd_33 = %.4g N.m/rad\n', K22_hyd, K33_hyd);
fprintf('  Hull: |b_t|<=%.3f m, h_t in [%.3f, %.3f] m\n', bt_max, ht_min, ht_max);
fprintf('  Falnes weights: nu = [%g, %g, %g]\n', nu);

% Baseline reference values
om_n_baseline = res.omega_n_target(:);
T_n_baseline  = res.T_n_target_s(:);
fprintf('\n  Baseline script placement: T_n = [%.3f, %.3f, %.3f] s\n', T_n_baseline);
fprintf('  Baseline eta_F (script) = %.4f\n', res.eta_F);

%% PHASE 1 -- bounds and seeds
% Heave omega_n >= hydrostatic (passive)
% Pitch omega_n >= hydrostatic (passive)
% Surge omega_n free in 90% band
lb = [omega_L;     om_hydro_h;  om_hydro_p;  deg2rad(30); 1   ];
ub = [omega_H;     omega_H;     omega_H;     deg2rad(60); 1e5 ];

rng(42);
nstarts = 16;
seeds = zeros(nstarts, 5);
% Hand-picked
seeds(1,:)  = [0.5;   om_hydro_h*1.001;  om_hydro_p*1.001;  deg2rad(45);   500];
seeds(2,:)  = [0.7;   om_hydro_h*1.001;  om_hydro_p*1.001;  deg2rad(45);  2000];
seeds(3,:)  = [1.0;   om_hydro_h*1.001;  om_hydro_p*1.001;  deg2rad(45);  5000];
seeds(4,:)  = [0.5;   1.5;               om_hydro_p*1.001;  deg2rad(45);  3000];
seeds(5,:)  = [0.5;   om_hydro_h*1.001;  1.7;               deg2rad(45);  1000];
seeds(6,:)  = [0.5;   om_hydro_h*1.001;  om_hydro_p*1.001;  deg2rad(35);   500];
seeds(7,:)  = [0.5;   om_hydro_h*1.001;  om_hydro_p*1.001;  deg2rad(55);   500];
seeds(8,:)  = [0.7;   1.5;               1.7;               deg2rad(50);  2000];
for s = 9:nstarts
    seeds(s,1) = lb(1) + (ub(1)-lb(1))*rand;
    seeds(s,2) = lb(2) + (ub(2)-lb(2))*rand;
    seeds(s,3) = lb(3) + (ub(3)-lb(3))*rand;
    seeds(s,4) = lb(4) + (ub(4)-lb(4))*rand;
    seeds(s,5) = exp(log(lb(5)) + (log(ub(5))-log(lb(5)))*rand);
end

%% PHASE 2 -- optimise
fprintf('\n  Running fmincon SQP with %d multistart restarts ...\n', nstarts);
obj      = @(x) -eta_power_compat(x, data);
nonlcon  = @(x) feasibility_compat(x, data);
optsfm   = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'TolFun', 1e-10, 'TolX', 1e-8, 'MaxIterations', 1500, ...
    'MaxFunctionEvaluations', 8000, 'ConstraintTolerance', 1e-8);

bestE = -inf; bestX = nan(5,1);
results = nan(nstarts, 9);   % [eta, om1, om2, om3, a_deg, B_PTO, K_moor, K_PTO, ell]
fprintf('\n  start  Tn_s    Tn_h    Tn_p    a(d)   B_PTO     K_PTO     ell    eta_pwr\n');
for s = 1:nstarts
    x0 = seeds(s,:)';
    try
        [x_s, neg_e_s, exf] = fmincon(obj, x0, [], [], [], [], lb, ub, nonlcon, optsfm);
        e_s = -neg_e_s;
        % Recompute derived quantities at x_s
        [~, deriv] = eta_power_compat(x_s, data);
        Tn = 2*pi ./ x_s(1:3);
        results(s,:) = [e_s, x_s(1), x_s(2), x_s(3), rad2deg(x_s(4)), x_s(5), ...
                        deriv.K_moor, deriv.K_PTO, deriv.ell];
        fprintf('   %2d  %5.2f   %5.2f   %5.2f   %4.1f  %8.2f  %8.2f  %.3f  %.4f\n', s, ...
            Tn(1), Tn(2), Tn(3), rad2deg(x_s(4)), x_s(5), deriv.K_PTO, deriv.ell, e_s);
        if e_s > bestE
            bestE = e_s;  bestX = x_s;
        end
    catch ME
        fprintf('   %2d  FAILED: %s\n', s, ME.message);
    end
end

%% PHASE 3 -- diagnostics at the optimum
[~, deriv_opt] = eta_power_compat(bestX, data);
om_n_opt = bestX(1:3);
T_n_opt  = 2*pi ./ om_n_opt;
a_opt    = bestX(4);
B_PTO_opt = bestX(5);
K_PTO_opt = deriv_opt.K_PTO;
ell_opt  = deriv_opt.ell;
K_moor_opt = deriv_opt.K_moor;
g_opt    = [2*cos(a_opt)^2; 2*sin(a_opt)^2; 2*ell_opt^2];
B_kk_opt = B_PTO_opt * g_opt;
K_kk_opt = K_PTO_opt * g_opt;

B_rad_opt = max([bem.B_kk_func(1, om_n_opt(1));
                 bem.B_kk_func(2, om_n_opt(2));
                 bem.B_kk_func(3, om_n_opt(3))], 0);
A_at_opt  = max([bem.A_kk_func(1, om_n_opt(1));
                 bem.A_kk_func(2, om_n_opt(2));
                 bem.A_kk_func(3, om_n_opt(3))], 0);
beta_opt  = B_kk_opt ./ max(B_rad_opt, eps);
eta_match = 4*beta_opt ./ (1+beta_opt).^2;
Q_opt     = om_n_opt .* (M_diag + A_at_opt) ./ max(B_rad_opt + B_kk_opt, eps);
B13_opt   = -2*B_PTO_opt*ell_opt*cos(a_opt);
K13_opt   = -2*K_PTO_opt*ell_opt*cos(a_opt);

% Per-mode contribution to eta_power
contrib = zeros(3,1);
for k = 1:3
    det = omega_grid./om_n_opt(k) - om_n_opt(k)./omega_grid;
    Pk  = 1 ./ (1 + Q_opt(k)^2 .* det.^2);
    contrib(k) = nu(k) * eta_match(k) * trapz(omega_grid(mask), Pk(mask) .* SF(mask)) / denom_int;
end

fprintf('\n================================================================\n');
fprintf('  GLOBAL OPTIMUM  (passive-feasible, stiffness-compatible)\n');
fprintf('================================================================\n');
fprintf('  Placed T_n (s) : surge %.4f, heave %.4f, pitch %.4f\n', T_n_opt);
fprintf('  omega_n (r/s)  : %.4f, %.4f, %.4f\n', om_n_opt);
fprintf('  alpha          : %.3f deg\n', rad2deg(a_opt));
fprintf('  ell            : %.4f m  (derived)\n', ell_opt);
fprintf('  K_PTO          : %.4g N/m  (per-tether scalar, derived)\n', K_PTO_opt);
fprintf('  B_PTO          : %.4g N.s/m (per-tether scalar)\n', B_PTO_opt);
fprintf('  K_moor (surge) : %.4g N/m  (passive spring, derived)\n', K_moor_opt);
fprintf('  eta_power*     : %.4f\n', bestE);

fprintf('\n  Per-mode diagnostics:\n');
fprintf('    Mode    B_PTO,kk     B_rad        beta      eta_match   Q_k    contrib\n');
mode_lbl = {'surge';'heave';'pitch'};
for k = 1:3
    fprintf('    %-5s  %10.4g   %10.4g  %.4f    %.4f    %6.2f   %.4f\n', mode_lbl{k}, ...
        B_kk_opt(k), B_rad_opt(k), beta_opt(k), eta_match(k), Q_opt(k), contrib(k));
end

fprintf('\n  Stiffness diagnostics (passive realisability check):\n');
fprintf('    K_PTO,11 = %.4g  (from matrix)\n', K_kk_opt(1));
fprintf('    K_PTO,22 = %.4g  (matches heave req)\n', K_kk_opt(2));
fprintf('    K_PTO,33 = %.4g  (matches pitch req)\n', K_kk_opt(3));
fprintf('    K_moor   = %.4g  (= S_1 - K_PTO,11, must be >= 0)\n', K_moor_opt);
fprintf('    K_total,11 = K_PTO,11 + K_moor = %.4g\n', K_kk_opt(1) + K_moor_opt);

fprintf('\n  Stage-2 gap check (should be near zero by construction):\n');
S1 = om_n_opt(1)^2 * (M + A_at_opt(1));
S2 = om_n_opt(2)^2 * (M + A_at_opt(2)) - K22_hyd;
S3 = om_n_opt(3)^2 * (Iyy + A_at_opt(3)) - K33_hyd;
fprintf('    Heave:  required S_2 = %.4g,  delivered 2 K_PTO sin^2 a = %.4g\n', S2, 2*K_PTO_opt*sin(a_opt)^2);
fprintf('    Pitch:  required S_3 = %.4g,  delivered 2 K_PTO ell^2   = %.4g\n', S3, 2*K_PTO_opt*ell_opt^2);
fprintf('    Surge:  required S_1 = %.4g,  delivered 2 K_PTO cos^2 a + K_moor = %.4g\n', ...
    S1, 2*K_PTO_opt*cos(a_opt)^2 + K_moor_opt);

fprintf('\n  Hull feasibility at the optimum:\n');
ell_max_a = bt_max*sin(a_opt) + ht_max*cos(a_opt);
fprintf('    ell* = %.4f m,  feasible max = %.4f m,  margin %.1f%%\n', ...
    ell_opt, ell_max_a, 100*(1 - ell_opt/ell_max_a));
fprintf('    h_t (m)    b_t (m)    feasible?\n');
ht_values = [0, 1.0, -1.0, ht_max, -2.26, ht_min];
for ht = ht_values
    bt = (ell_opt - ht*cos(a_opt)) / sin(a_opt);
    flag = 'YES';
    if abs(bt) > bt_max
        flag = 'NO (|b_t|>half-beam)';
    elseif (ht > ht_max) || (ht < ht_min)
        flag = 'NO (h_t outside hull)';
    end
    fprintf('     %+6.3f    %+7.3f    %s\n', ht, bt, flag);
end

fprintf('\n  Off-diagonal entries (diagnostic):\n');
fprintf('    B_PTO,13 = %.4g  (vs B_PTO,11 = %.4g, B_PTO,33 = %.4g)\n', ...
    B13_opt, B_kk_opt(1), B_kk_opt(3));
fprintf('    K_PTO,13 = %.4g  (vs K_PTO,11 = %.4g, K_PTO,33 = %.4g)\n', ...
    K13_opt, K_kk_opt(1), K_kk_opt(3));

%% PHASE 4 -- side-by-side three-way comparison
fprintf('\n================================================================\n');
fprintf('  THREE-WAY COMPARISON\n');
fprintf('================================================================\n');
fprintf('                              | Initial script | Free eta_power | Passive-comp\n');
fprintf('  ----------------------------+----------------+----------------+----------------\n');
fprintf('  T_n surge        (s)        | %14.3f | %14.3f | %14.3f\n', ...
    T_n_baseline(1), 5.918, T_n_opt(1));
fprintf('  T_n heave        (s)        | %14.3f | %14.3f | %14.3f\n', ...
    T_n_baseline(2), 9.116, T_n_opt(2));
fprintf('  T_n pitch        (s)        | %14.3f | %14.3f | %14.3f\n', ...
    T_n_baseline(3), 5.712, T_n_opt(3));
fprintf('  alpha            (deg)      | %14s | %14.2f | %14.2f\n', ...
    'n/a', 49.12, rad2deg(a_opt));
fprintf('  ell              (m)        | %14s | %14.4f | %14.4f\n', ...
    'n/a', 0.0168, ell_opt);
fprintf('  B_PTO            (N.s/m)    | %14.1f | %14.1f | %14.1f\n', ...
    closure.B_PTO_scalar, 4247.5, B_PTO_opt);
fprintf('  K_PTO needed             (N/m) heave|  %12.4g | %14s | %14.4g\n', ...
    -39530, 'n/a (free)', K_kk_opt(2));
fprintf('  K_PTO,22 required        (N/m)      |    (negative)  |   (negative)   | %14.4g\n', ...
    K_kk_opt(2));
fprintf('  K_PTO,33 required (N.m/rad) pitch   |  %12.4g | %14s | %14.4g\n', ...
    -27840, '(negative)', K_kk_opt(3));
fprintf('  K_moor surge          (N/m)         | %14.1f | %14s | %14.4g\n', ...
    6005.4, 'n/a', K_moor_opt);
fprintf('  Passive realisable?                 | %14s | %14s | %14s\n', ...
    'NO (K_PTO<0)', 'NO (K_PTO<0)', 'YES');
fprintf('  Script eta_F                        | %14.4f | %14s | %14s\n', ...
    res.eta_F, '-', '-');
fprintf('  eta_power (Lorentzian + matching)   | %14.4f | %14.4f | %14.4f\n', ...
    0.0401, 0.0545, bestE);
fprintf('  (baseline ref: same placement and matrix as script, ~upper bound at that placement)\n');

%% Save
save('eta_power_compat_results.mat', 'bestE','bestX','results', ...
    'om_n_opt','T_n_opt','a_opt','ell_opt','B_PTO_opt','K_PTO_opt','K_moor_opt', ...
    'beta_opt','eta_match','Q_opt','contrib','B_kk_opt','K_kk_opt','B_rad_opt','A_at_opt');
fprintf('\n  Saved to eta_power_compat_results.mat\n');

%% =====================================================
%% LOCAL FUNCTIONS
%% =====================================================

function [eta, deriv] = eta_power_compat(x, d)
    % Derive (K_PTO, ell, K_moor) from placement under stiffness compatibility
    om_n  = x(1:3);
    a     = x(4);
    B_PTO = x(5);
    A_at = max([d.bem.A_kk_func(1, om_n(1));
                d.bem.A_kk_func(2, om_n(2));
                d.bem.A_kk_func(3, om_n(3))], 0);
    S1 = om_n(1)^2 * (d.M   + A_at(1));
    S2 = om_n(2)^2 * (d.M   + A_at(2)) - d.K22_hyd;
    S3 = om_n(3)^2 * (d.Iyy + A_at(3)) - d.K33_hyd;
    % Derived geometry & passive stiffness
    K_PTO = S2 / (2*sin(a)^2);
    if S2 > 0 && S3 >= 0
        ell = sin(a) * sqrt(max(S3, 0)/S2);
    else
        % Infeasible region; objective penalty
        eta = -1e3;
        deriv = struct('K_PTO', NaN, 'ell', NaN, 'K_moor', NaN);
        return;
    end
    K_moor = S1 - 2*K_PTO*cos(a)^2;
    % B_PTO matrix diagonal
    g = [2*cos(a)^2; 2*sin(a)^2; 2*ell^2];
    B_PTO_kk = B_PTO * g;
    B_rad = max([d.bem.B_kk_func(1, om_n(1));
                 d.bem.B_kk_func(2, om_n(2));
                 d.bem.B_kk_func(3, om_n(3))], 0);
    beta = B_PTO_kk ./ max(B_rad, eps);
    eta_match = 4*beta ./ (1+beta).^2;
    Q = om_n .* (d.M_diag + A_at) ./ max(B_rad + B_PTO_kk, eps);
    eta = 0;
    for k = 1:3
        det = d.omega./om_n(k) - om_n(k)./d.omega;
        Pk  = 1 ./ (1 + Q(k)^2 .* det.^2);
        eta = eta + d.nu(k) * eta_match(k) * trapz(d.omega(d.mask), Pk(d.mask) .* d.SF(d.mask));
    end
    eta = eta / d.denom;
    deriv = struct('K_PTO', K_PTO, 'ell', ell, 'K_moor', K_moor, ...
                   'S1', S1, 'S2', S2, 'S3', S3);
end

function [c, ceq] = feasibility_compat(x, d)
    om_n  = x(1:3);
    a     = x(4);
    A_at = max([d.bem.A_kk_func(1, om_n(1));
                d.bem.A_kk_func(2, om_n(2));
                d.bem.A_kk_func(3, om_n(3))], 0);
    S1 = om_n(1)^2 * (d.M   + A_at(1));
    S2 = om_n(2)^2 * (d.M   + A_at(2)) - d.K22_hyd;
    S3 = om_n(3)^2 * (d.Iyy + A_at(3)) - d.K33_hyd;
    % (a) S_2 >= 0 (heave passive); redundant with lb on om_n,2 but keep
    c1 = -S2;
    % (b) S_3 >= 0 (pitch passive)
    c2 = -S3;
    % (c) K_moor >= 0:  S_1 - S_2 cot^2 a >= 0
    c3 = S2*cos(a)^2/max(sin(a)^2,eps) - S1;
    % (d) ell <= ell_max
    if S2 > 0 && S3 >= 0
        ell = sin(a)*sqrt(max(S3,0)/S2);
        ell_max = d.bt_max*sin(a) + d.ht_max*cos(a);
        c4 = ell - ell_max;
    else
        c4 = -1;  % trivially feasible if S2/S3 infeasible (the c1/c2 handle that)
    end
    c   = [c1; c2; c3; c4];
    ceq = [];
end
