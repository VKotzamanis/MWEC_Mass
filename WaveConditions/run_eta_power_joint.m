% run_eta_power_joint.m
%   Joint Stage-1 optimisation of placement + damping + geometry under a
%   corrected broadband absorbed-power figure of merit (eta_power).
%
%   ============================================================
%   ASSUMPTIONS (state up front; defend in peer review)
%   ============================================================
%   A1. SDOF Lorentzian retained per mode:
%         P_k(omega) = 1 / (1 + Q_k^2 (omega/omega_n,k - omega_n,k/omega)^2)
%       same as the script's existing placement -- maintains parity.
%       Off-diagonal B_PTO,13 NOT in the objective; reported diagnostically.
%       Full 3x3 coupled RAO is a Stage-2 calculation.
%
%   A2. Per-mode resonant matching efficiency multiplies the Lorentzian:
%         eta_match,k = 4 beta_k / (1 + beta_k)^2,  beta_k = B_PTO,kk / B_kk^rad(omega_n,k)
%       Canonical Falnes/Evans factor for SDOF resistance-matched absorption.
%
%   A3. Damping matrix structure from the linearised symmetric-tether derivation:
%         B_PTO,11 = 2 B_PTO cos^2 a
%         B_PTO,22 = 2 B_PTO sin^2 a
%         B_PTO,33 = 2 B_PTO ell^2
%         B_PTO,13 = -2 B_PTO ell cos a   (off-diagonal; diagnostic only)
%
%   A4. Hull feasibility (hard nonlinear inequality):
%         ell <= b_t_max sin(a) + h_t_max cos(a)
%       guarantees at least one (h_t, b_t) on the hull envelope satisfies
%       ell = h_t cos a + b_t sin a with b_t in [-b_t_max, b_t_max],
%       h_t in [h_t_min, h_t_max].  Above-CG attachment assumed (ell>0).
%
%   A5. Placement zones (inherited from script's spectral partition):
%         surge omega_n in [omega_L, omega_partition]   (long periods)
%         heave omega_n in [omega_L, omega_H]           (full 90% band)
%         pitch omega_n in [omega_partition, omega_H]   (short periods)
%
%   A6. Mooring angle bounds:  alpha in [30 deg, 60 deg].
%   A7. B_PTO bounds:  [1, 1e5] N.s/m  (loose; expect interior optimum now).
%
%   A8. No K_PTO optimisation: the placed omega_n are design targets to be
%       realised by Stage-2 mooring.  Stiffness budget is reported but not
%       optimised here (consistent with the Stage-1 deliverable).
%
%   A9. Falnes weights nu = [2; 1; 2] (axisymmetric deep-water).
%
%   A10. Frozen BEM cache (interpolants from run_MWEC_tuning's output).
%        No re-running of BEM.
%
%   Objective:
%     eta_power = sum_k nu_k * eta_match,k * integral( P_k S_ew / omega^3 ) / denom
%   normalised by  denom = (sum nu) * integral S_ew / omega^3  over the 90% band.
%
%   Optimiser:  fmincon SQP with 16 multistart restarts.

clear; clc;
addpath(pwd);

%% ============================================================
%% PHASE 0 -- BASELINE FROM run_MWEC_tuning (for context only)
%% ============================================================
climate = fullfile(pwd, 'WIS_Output_WAM', 'ST84040_climate_grid.mat');
props   = fullfile(pwd, 'WEC_Mass_UHPC.mat');
opts    = struct('plot', false, 'verbose', false);

fprintf('================================================================\n');
fprintf('  ETA_POWER JOINT OPTIMISATION  (placement + damping + geometry)\n');
fprintf('================================================================\n');
fprintf('  Baseline call to run_MWEC_tuning ...\n');
res = run_MWEC_tuning(climate, props, opts);

bem     = res.bem;
closure = res.closure;
M_diag  = [closure.M; closure.M; closure.Iyy];
nu      = [2; 1; 2];

omega_grid = res.climate.omega(:);
S_ew       = res.climate.S_ew(:);
omega_L    = res.climate.omega_L;
omega_H    = res.climate.omega_H;
omega_part = res.climate.spectral_info.partition_omega;
mask       = (omega_grid >= omega_L) & (omega_grid <= omega_H);

SF        = S_ew ./ max(omega_grid, 1e-9).^3;
denom_int = sum(nu) * trapz(omega_grid(mask), SF(mask));

% Hull bounds
fp_S = load(props);  fp = fp_S.final_props;
cs = fp.cross_section;
bt_max = max(abs(cs(:,1)));         % half-beam
ht_max = max(cs(:,2));              % top
ht_min = min(cs(:,2));              % keel

% Hydrostatic refs
T_n_hydro_h = bem.T_n_hydrostatic(2);
T_n_hydro_p = bem.T_n_hydrostatic(3);
T_p = 2*pi / res.climate.IEC.omega_p;

% Bundle data
data.bem = bem; data.M_diag = M_diag; data.nu = nu;
data.omega = omega_grid; data.SF = SF; data.mask = mask; data.denom = denom_int;

fprintf('\n  Climate         : 90%% band T in [%.2f, %.2f] s,  T_p = %.2f s\n', ...
    2*pi/omega_H, 2*pi/omega_L, T_p);
fprintf('  Partition omega : %.4f r/s (T = %.3f s)\n', omega_part, 2*pi/omega_part);
fprintf('  Hydrostatic T_n : heave %.3f s,  pitch %.3f s\n', T_n_hydro_h, T_n_hydro_p);
fprintf('  Hull bounds     : |b_t| <= %.3f m,  h_t in [%.3f, %.3f] m\n', bt_max, ht_min, ht_max);
fprintf('  Falnes weights  : nu = [%g, %g, %g]\n', nu);

% Baseline at the script's existing placement (with B_PTO_scalar applied
% uniformly per the current code, NOT the matrix structure)
om_n_baseline = res.omega_n_target(:);
T_n_baseline  = res.T_n_target_s(:);
fprintf('\n  Baseline placement (current pipeline): T_n = [%.3f, %.3f, %.3f] s\n', T_n_baseline);
fprintf('  Baseline eta_F (script) : %.4f\n', res.eta_F);

% Equivalent baseline eta_power if we treat the script's scalar as the
% matrix scalar with alpha=45deg, ell=k_gyr -- documented as an upper-bound
% reference, not an apples-to-apples comparison.
x_baseline_ref = [om_n_baseline; deg2rad(45); 0.0625; closure.B_PTO_scalar];
eta_baseline_ref = eta_power_obj(x_baseline_ref, data);
fprintf('  Reference eta_power at baseline placement w/ Case-B geometry: %.4f\n', eta_baseline_ref);

%% ============================================================
%% PHASE 1 -- BOUNDS, MULTISTART SEEDS
%% ============================================================
lb = [omega_L;     omega_L;     omega_part; deg2rad(30); 0.01;  1   ];
ub = [omega_part;  omega_H;     omega_H;    deg2rad(60); 1.85;  1e5 ];

% Hand-picked seeds (8) + random (8)
rng(42);
nstarts = 16;
seeds = zeros(nstarts, 6);
seeds(1, :) = [0.50; 0.55; 1.10; deg2rad(45); 0.30; 500  ];   % near baseline
seeds(2, :) = [0.50; 2*pi/T_n_hydro_h; 2*pi/T_n_hydro_p; deg2rad(45); 0.50; 1000];   % at hydrostatic
seeds(3, :) = [0.70; 0.80; 1.40; deg2rad(45); 0.50; 300  ];   % clustered
seeds(4, :) = [0.50; 1.00; 1.30; deg2rad(40); 1.00; 2000 ];
seeds(5, :) = [0.45; 0.85; 1.50; deg2rad(50); 0.20; 800  ];
seeds(6, :) = [0.55; 1.20; 1.40; deg2rad(35); 0.80; 200  ];
seeds(7, :) = [0.45; 0.65; 1.80; deg2rad(55); 0.05; 100  ];
seeds(8, :) = [0.85; 1.30; 1.30; deg2rad(45); 0.60; 1500 ];
for s = 9:nstarts
    seeds(s,1) = lb(1) + (ub(1)-lb(1))*rand;
    seeds(s,2) = lb(2) + (ub(2)-lb(2))*rand;
    seeds(s,3) = lb(3) + (ub(3)-lb(3))*rand;
    seeds(s,4) = lb(4) + (ub(4)-lb(4))*rand;
    seeds(s,5) = lb(5) + (ub(5)-lb(5))*rand;
    seeds(s,6) = exp(log(lb(6)) + (log(ub(6))-log(lb(6)))*rand);
end
% Clamp ell at each seed to hull feasibility margin
for s = 1:nstarts
    aa = seeds(s,4);
    ell_max = bt_max*sin(aa) + ht_max*cos(aa);
    if seeds(s,5) > 0.95*ell_max
        seeds(s,5) = 0.95*ell_max;
    end
end

%% ============================================================
%% PHASE 2 -- OPTIMISE
%% ============================================================
fprintf('\n  Running fmincon SQP with %d multistart restarts ...\n', nstarts);
obj_neg  = @(x) -eta_power_obj(x, data);
nonlcon  = @(x) hull_constraint(x, bt_max, ht_max);
opts_fmin = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'TolFun', 1e-10, 'TolX', 1e-8, 'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000);

bestE = -inf; bestX = nan(6,1);
results = nan(nstarts, 8);   % [eta, om1, om2, om3, a, ell, B_PTO, exitflag]
fprintf('\n  start    Tn_s   Tn_h   Tn_p    a(deg)  ell    B_PTO    eta_power\n');
for s = 1:nstarts
    x0 = seeds(s,:)';
    try
        [x_s, neg_e_s, ex] = fmincon(obj_neg, x0, [], [], [], [], lb, ub, nonlcon, opts_fmin);
        e_s = -neg_e_s;
        Tn = 2*pi ./ x_s(1:3);
        results(s,:) = [e_s, x_s(1), x_s(2), x_s(3), x_s(4), x_s(5), x_s(6), ex];
        fprintf('   %2d  %6.2f %6.2f %6.2f  %5.1f  %.3f  %8.2f   %.4f\n', s, ...
            Tn(1), Tn(2), Tn(3), rad2deg(x_s(4)), x_s(5), x_s(6), e_s);
        if e_s > bestE
            bestE = e_s; bestX = x_s;
        end
    catch ME
        fprintf('   %2d  FAILED: %s\n', s, ME.message);
    end
end

%% ============================================================
%% PHASE 3 -- DIAGNOSTICS AT GLOBAL OPTIMUM
%% ============================================================
om_n_opt = bestX(1:3); a_opt = bestX(4); ell_opt = bestX(5); B_PTO_opt = bestX(6);
T_n_opt = 2*pi ./ om_n_opt;
g_opt = [2*cos(a_opt)^2; 2*sin(a_opt)^2; 2*ell_opt^2];
B_PTO_kk = B_PTO_opt * g_opt;
B_rad_opt = max([bem.B_kk_func(1, om_n_opt(1));
                 bem.B_kk_func(2, om_n_opt(2));
                 bem.B_kk_func(3, om_n_opt(3))], 0);
A_at_opt  = max([bem.A_kk_func(1, om_n_opt(1));
                 bem.A_kk_func(2, om_n_opt(2));
                 bem.A_kk_func(3, om_n_opt(3))], 0);
beta_opt  = B_PTO_kk ./ max(B_rad_opt, eps);
eta_match = 4*beta_opt ./ (1+beta_opt).^2;
Q_opt     = om_n_opt .* (M_diag + A_at_opt) ./ max(B_rad_opt + B_PTO_kk, eps);
B13_opt   = -2*B_PTO_opt*ell_opt*cos(a_opt);

% Per-mode contribution to eta_power
contrib = zeros(3,1);
for k = 1:3
    det = omega_grid./om_n_opt(k) - om_n_opt(k)./omega_grid;
    Pk  = 1 ./ (1 + Q_opt(k)^2 .* det.^2);
    contrib(k) = nu(k) * eta_match(k) * trapz(omega_grid(mask), Pk(mask) .* SF(mask)) / denom_int;
end

fprintf('\n================================================================\n');
fprintf('  GLOBAL OPTIMUM\n');
fprintf('================================================================\n');
fprintf('  Placed T_n (s)     : surge %.4f, heave %.4f, pitch %.4f\n', T_n_opt);
fprintf('  Placed omega (r/s) : %.4f, %.4f, %.4f\n', om_n_opt);
fprintf('  Geometry           : alpha = %.3f deg,  ell = %.4f m\n', rad2deg(a_opt), ell_opt);
fprintf('  Damping            : B_PTO = %.3f N.s/m\n', B_PTO_opt);
fprintf('  eta_power*         : %.4f\n', bestE);
fprintf('\n  Per-mode diagnostics:\n');
fprintf('    Mode    B_PTO,kk     B_rad        beta      eta_match   Q_k    contrib\n');
mode_lbl = {'surge';'heave';'pitch'};
for k = 1:3
    fprintf('    %-5s  %10.4g   %10.4g  %.4f    %.4f   %.3f   %.4f\n', mode_lbl{k}, ...
        B_PTO_kk(k), B_rad_opt(k), beta_opt(k), eta_match(k), Q_opt(k), contrib(k));
end

fprintf('\n  Stage-2 gap (placed vs hydrostatic):\n');
fprintf('    Heave: target %.3f s,  hydrostatic %.3f s,  gap %+.3f s\n', ...
    T_n_opt(2), T_n_hydro_h, T_n_opt(2)-T_n_hydro_h);
fprintf('    Pitch: target %.3f s,  hydrostatic %.3f s,  gap %+.3f s\n', ...
    T_n_opt(3), T_n_hydro_p, T_n_opt(3)-T_n_hydro_p);

fprintf('\n  Off-diagonal coupling (diagnostic):\n');
fprintf('    B_PTO,13 = %.4g N.s/m   (B_PTO,11 = %.4g, B_PTO,33 = %.4g)\n', ...
    B13_opt, B_PTO_kk(1), B_PTO_kk(3));
fprintf('    |B13|/sqrt(B11*B33) = %.4f  (rank-1 sub-block: =1 means singular)\n', ...
    abs(B13_opt)/sqrt(max(B_PTO_kk(1)*B_PTO_kk(3), eps)));

fprintf('\n  Hull-feasibility at the optimum:\n');
ell_max_a = bt_max*sin(a_opt) + ht_max*cos(a_opt);
fprintf('    ell* = %.4f m;  feasible max at this alpha = %.4f m  (margin %.1f%%)\n', ...
    ell_opt, ell_max_a, 100*(1 - ell_opt/ell_max_a));
fprintf('    h_t (m)   b_t (m)    feasible?\n');
ht_values = [0, 1.0, -1.0, ht_max, -2.26, ht_min];
for ht = ht_values
    bt = (ell_opt - ht*cos(a_opt)) / sin(a_opt);
    flag = 'YES';
    if abs(bt) > bt_max
        flag = 'NO (|b_t|>half-beam)';
    elseif (ht > ht_max) || (ht < ht_min)
        flag = 'NO (h_t off hull)';
    end
    fprintf('     %+6.3f    %+7.3f    %s\n', ht, bt, flag);
end

fprintf('\n  Stiffness implied by placement (informational only):\n');
% Heave: K_22_total = ω_n,2^2 * (M + A_22(omega_n,2))
K22_total_needed = om_n_opt(2)^2 * (closure.M + A_at_opt(2));
K22_hydro = fp.K_hydro(2,2);
K22_PTO_needed = K22_total_needed - K22_hydro;
% Pitch
K33_total_needed = om_n_opt(3)^2 * (closure.Iyy + A_at_opt(3));
K33_hydro = fp.K_hydro(3,3);
K33_PTO_needed = K33_total_needed - K33_hydro;
% Surge
K11_total_needed = om_n_opt(1)^2 * (closure.M + A_at_opt(1));
fprintf('    Heave K_PTO,22 needed = %.4g N/m  (hydrostatic = %.4g)\n', K22_PTO_needed, K22_hydro);
fprintf('    Pitch K_PTO,33 needed = %.4g N.m/rad  (hydrostatic = %.4g)\n', K33_PTO_needed, K33_hydro);
fprintf('    Surge K_PTO,11 needed = %.4g N/m  (hydrostatic = 0)\n', K11_total_needed);

%% Save
save('eta_power_joint_results.mat', 'bestE', 'bestX', 'results', ...
    'om_n_opt', 'T_n_opt', 'a_opt', 'ell_opt', 'B_PTO_opt', ...
    'beta_opt', 'eta_match', 'Q_opt', 'contrib', 'B_PTO_kk', 'B_rad_opt', 'B13_opt');
fprintf('\n  Saved to eta_power_joint_results.mat\n');

%% ============================================================
%% LOCAL FUNCTIONS
%% ============================================================
function eta = eta_power_obj(x, d)
    om_n = x(1:3); a = x(4); ell = x(5); B_PTO = x(6);
    g = [2*cos(a)^2; 2*sin(a)^2; 2*ell^2];
    B_rad = max([d.bem.B_kk_func(1, om_n(1));
                 d.bem.B_kk_func(2, om_n(2));
                 d.bem.B_kk_func(3, om_n(3))], 0);
    A_at = max([d.bem.A_kk_func(1, om_n(1));
                d.bem.A_kk_func(2, om_n(2));
                d.bem.A_kk_func(3, om_n(3))], 0);
    B_PTO_kk = B_PTO * g;
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
end

function [c, ceq] = hull_constraint(x, bt_max, ht_max)
    a = x(4); ell = x(5);
    ell_max = bt_max*sin(a) + ht_max*cos(a);
    c   = ell - ell_max;      % must be <= 0
    ceq = [];
end
