%% debug_CG_transform.m — verify the single-step CG congruence vs final_props
%
%   Phase-0 gate for the WEC_GM v5.0 refactor.
%
%   We want to confirm that
%       A_55,CG_world  =  (T'_world · A_3x3_O · T_world)(3,3)
%   with
%       T_world = [1 0 -z_cg_world ; 0 1 0 ; 0 0 1]
%   reproduces  final_props.A_full(3,3)  to within ~1 %.
%
%   If it does not, we dump the production chain step-by-step so the bug
%   localises to:  (a) wrong cache index lookup,
%                  (b) wrong reference frame for z_cg,
%                  (c) interp_cache_at_vs versus interpolate_3x3_matrix path,
%                  (d) extraction index mismatch.

clear; clc;

this_dir     = fileparts(mfilename('fullpath'));
project_root = fileparts(this_dir);
mass_mat     = fullfile(project_root, 'WEC_Mass_UHPC.mat');

addpath(project_root);

fprintf('\n=== Phase 0  CG-transform debug ====================================\n');
fprintf('  Mass mat : %s\n', mass_mat);

S           = load(mass_mat);
cfg         = S.results.config;
fp          = S.final_props;
cache       = cfg.hydro_cache;

vs          = fp.vertical_shift;
cg_z_world  = fp.CG_total(3);
cg_z_body   = cg_z_world - vs;

fprintf('  vs               = %+.4f m\n', vs);
fprintf('  cg_z_world       = %+.4f m   (= final_props.CG_total(3))\n', cg_z_world);
fprintf('  cg_z_body        = %+.4f m   (= cg_z_world - vs)\n',         cg_z_body);
fprintf('  final A_full(3,3)= %.4f kg*m^2   (= reference target)\n',    fp.A_full(3,3));
fprintf('  final A_full(1,1)= %.4f kg      (surge-surge at CG)\n',       fp.A_full(1,1));
fprintf('  final A_full(1,3)= %.4f kg*m    (surge-pitch at CG)\n',       fp.A_full(1,3));
fprintf('  final A55        = %.4f kg*m^2\n',                            fp.A55);

%% ----- Step 1: identify the bracketing cache drafts -----
[drafts_sorted, ord] = sort(cache.drafts(:));
fprintf('\n  Cache drafts (sorted)        : %s\n', mat2str(drafts_sorted', 4));
fprintf('  Cache z_cg (HAMS) sorted     : %s\n', mat2str(cache.z_cg(ord)', 4));

k_hi = find(drafts_sorted >= vs, 1, 'first');
k_lo = k_hi - 1;
vs_lo = drafts_sorted(k_lo);   vs_hi = drafts_sorted(k_hi);
t_interp = (vs - vs_lo) / (vs_hi - vs_lo);

idx_lo = ord(k_lo);
idx_hi = ord(k_hi);

fprintf('\n  vs = %.4f  brackets  vs_lo=%.4f (idx %d)  vs_hi=%.4f (idx %d)  t=%.3f\n', ...
        vs, vs_lo, idx_lo, vs_hi, idx_hi, t_interp);

z_cg_hams_lo = cache.z_cg(idx_lo);
z_cg_hams_hi = cache.z_cg(idx_hi);
z_cg_hams_interp = (1-t_interp)*z_cg_hams_lo + t_interp*z_cg_hams_hi;
fprintf('  HAMS-CG bracketing values    : %.4f  ->  %.4f  (interpolated %.4f)\n', ...
        z_cg_hams_lo, z_cg_hams_hi, z_cg_hams_interp);

dz_delta = cg_z_world - z_cg_hams_interp;
fprintf('  Δz = z_cg_world - z_cg_hams  : %+.4f m\n', dz_delta);

%% ----- Step 2: take raw cache.A_inf at the two bracket drafts -----
A_inf_lo = cache.A_inf{idx_lo};        % 6x6 at HAMS reference
A_inf_hi = cache.A_inf{idx_hi};
A_inf_interp = (1-t_interp)*A_inf_lo + t_interp*A_inf_hi;

fprintf('\n  cache.A_inf{lo}(5,5) = %.3e   cache.A_inf{hi}(5,5) = %.3e   interp = %.3e kg*m^2\n', ...
        A_inf_lo(5,5), A_inf_hi(5,5), A_inf_interp(5,5));
fprintf('  cache.A_inf{lo}(1,1) = %.3e   cache.A_inf{hi}(1,1) = %.3e   interp = %.3e kg\n', ...
        A_inf_lo(1,1), A_inf_hi(1,1), A_inf_interp(1,1));
fprintf('  cache.A_inf{lo}(1,5) = %+.3e  cache.A_inf{hi}(1,5) = %+.3e  interp = %+.3e kg*m\n', ...
        A_inf_lo(1,5), A_inf_hi(1,5), A_inf_interp(1,5));

A_3x3_O = A_inf_interp([1 3 5], [1 3 5]);
A_3x3_O = 0.5 * (A_3x3_O + A_3x3_O');     % symmetrise

%% ----- Test A: single-step transform with z_cg = cg_z_world -----
T_world  = [1, 0, -cg_z_world; 0, 1, 0; 0, 0, 1];
A_CG_test_world = T_world' * A_3x3_O * T_world;
A_CG_test_world = 0.5 * (A_CG_test_world + A_CG_test_world');

fprintf('\n  --- Test A: T(z_cg_world) one-step ---\n');
fprintf('     A_CG(3,3) = %.4f   (vs target %.4f)   ratio = %.3f\n', ...
        A_CG_test_world(3,3), fp.A_full(3,3), A_CG_test_world(3,3)/fp.A_full(3,3));
fprintf('     A_CG(1,1) = %.4f   (vs target %.4f)   ratio = %.3f\n', ...
        A_CG_test_world(1,1), fp.A_full(1,1), A_CG_test_world(1,1)/fp.A_full(1,1));
fprintf('     A_CG(1,3) = %+.4f  (vs target %+.4f)\n', ...
        A_CG_test_world(1,3), fp.A_full(1,3));

%% ----- Test B: single-step transform with z_cg = cg_z_body -----
T_body  = [1, 0, -cg_z_body; 0, 1, 0; 0, 0, 1];
A_CG_test_body = T_body' * A_3x3_O * T_body;
A_CG_test_body = 0.5 * (A_CG_test_body + A_CG_test_body');

fprintf('\n  --- Test B: T(z_cg_body) one-step ---\n');
fprintf('     A_CG(3,3) = %.4f   (vs target %.4f)   ratio = %.3f\n', ...
        A_CG_test_body(3,3), fp.A_full(3,3), A_CG_test_body(3,3)/fp.A_full(3,3));
fprintf('     A_CG(1,3) = %+.4f  (vs target %+.4f)\n', ...
        A_CG_test_body(1,3), fp.A_full(1,3));

%% ----- Test C: full PRODUCTION chain (the one final_props was built with) -----
% Step C1: at each cache draft, apply T(z_cg_hams_i) to A_inf_origin
A_at_hams_cg = cell(numel(cache.drafts), 1);
for i = 1:numel(cache.drafts)
    z_i = cache.z_cg(i);
    Ti  = [1, 0, -z_i; 0, 1, 0; 0, 0, 1];
    Ai  = cache.A_inf{i}([1 3 5], [1 3 5]);
    Ai  = 0.5 * (Ai + Ai');
    A_at_hams_cg{i} = 0.5 * ((Ti' * Ai * Ti) + (Ti' * Ai * Ti)');
end

% Step C2: linear interp in draft (same as interpolate_3x3_matrix)
A_at_hams_cg_interp = (1-t_interp)*A_at_hams_cg{idx_lo} + t_interp*A_at_hams_cg{idx_hi};

% Step C3: delta transform from HAMS-CG to target
T_delta = [1, 0, -dz_delta; 0, 1, 0; 0, 0, 1];
A_CG_chain = T_delta' * A_at_hams_cg_interp * T_delta;
A_CG_chain = 0.5 * (A_CG_chain + A_CG_chain');

fprintf('\n  --- Test C: TWO-STEP production chain (T(z_cg_hams) then T(dz)) ---\n');
fprintf('     A_at_hams_cg_interp(3,3) = %.4f kg*m^2  (intermediate, at z_cg_hams)\n', ...
        A_at_hams_cg_interp(3,3));
fprintf('     after delta T(dz=%+.4f) :  A_CG(3,3) = %.4f   target %.4f   ratio %.3f\n', ...
        dz_delta, A_CG_chain(3,3), fp.A_full(3,3), A_CG_chain(3,3)/fp.A_full(3,3));
fprintf('     A_CG(1,3) = %+.4f   target %+.4f\n', A_CG_chain(1,3), fp.A_full(1,3));

%% ----- Test D: check whether WEC_Core_Functions reproduces final_props -----
% Re-run interpolate_wamit_added_mass with vs and cg_z_world target.
fprintf('\n  --- Test D: re-run WEC_Core_Functions.interpolate_wamit_added_mass ---\n');
try
    [A11_d, A33_d, A55_d, ~, A_full_d, B_full_d] = ...
        WEC_Core_Functions.interpolate_wamit_added_mass(vs, cfg, cg_z_world);
    fprintf('     production A11/A33/A55 = %.4f / %.4f / %.4f\n', A11_d, A33_d, A55_d);
    fprintf('     production A_full(3,3) = %.4f   target %.4f   ratio %.3f\n', ...
            A_full_d(3,3), fp.A_full(3,3), A_full_d(3,3)/fp.A_full(3,3));
    fprintf('     production A_full(1,3) = %+.4f   target %+.4f\n', ...
            A_full_d(1,3), fp.A_full(1,3));
catch ME
    fprintf('     Error calling production interp: %s\n', ME.message);
end

%% ----- Test E: also compare A(omega_max) and final_props.A55  -----
% final_props.A55 should be the infinite-freq pitch added mass at the
% actual UHPC CG.  The last omega in cache is ~3.14 rad/s, which is NOT
% omega → infty.  But the cache also has A_inf (the actual ω → ∞ limit).
fprintf('\n  --- Test E: A55(ω_inf) extrapolation vs final_props.A55 ---\n');
A_omega_max_origin = (1-t_interp)*cache.A{idx_lo}(:,:,end) + t_interp*cache.A{idx_hi}(:,:,end);
A_omega_max_3x3 = A_omega_max_origin([1 3 5], [1 3 5]);
A_omega_max_3x3 = 0.5*(A_omega_max_3x3 + A_omega_max_3x3');
A_omega_max_CG = T_world' * A_omega_max_3x3 * T_world;
fprintf('     A55(ω_max)_CG using T_world  : %.4f kg*m^2\n', A_omega_max_CG(3,3));
A_omega_max_CG_body = T_body' * A_omega_max_3x3 * T_body;
fprintf('     A55(ω_max)_CG using T_body   : %.4f kg*m^2\n', A_omega_max_CG_body(3,3));
fprintf('     final_props.A55              : %.4f kg*m^2  ← target\n', fp.A55);

%% ----- Summary verdict -----
err_world = abs(A_CG_test_world(3,3) - fp.A_full(3,3)) / fp.A_full(3,3) * 100;
err_body  = abs(A_CG_test_body(3,3)  - fp.A_full(3,3)) / fp.A_full(3,3) * 100;
err_chain = abs(A_CG_chain(3,3)      - fp.A_full(3,3)) / fp.A_full(3,3) * 100;

fprintf('\n========================  VERDICT  ========================\n');
fprintf('  Target A_CG(3,3)            = %.4f kg*m^2\n', fp.A_full(3,3));
fprintf('  Test A  (T_world one-step ) = %.4f  err %.2f%%\n', A_CG_test_world(3,3), err_world);
fprintf('  Test B  (T_body  one-step ) = %.4f  err %.2f%%\n', A_CG_test_body(3,3),  err_body);
fprintf('  Test C  (production chain ) = %.4f  err %.2f%%\n', A_CG_chain(3,3),      err_chain);
fprintf('===========================================================\n');

if err_world < 2
    fprintf('  CONCLUSION: T(z_cg_world) is the correct one-step transform.\n');
    fprintf('              My WEC_GM v4 code is using the right formula.\n');
    fprintf('              The 13× discrepancy I quoted earlier came from\n');
    fprintf('              comparing UHPC final_props against STEEL-cache-index-1\n');
    fprintf('              values — apples to oranges.  No code bug.\n');
elseif err_body < 2
    fprintf('  CONCLUSION: T(z_cg_BODY) is correct — current code is WRONG.\n');
    fprintf('              Fix interpolate_pitch_BEM_at_design to use cg_z_body.\n');
elseif err_chain < 2
    fprintf('  CONCLUSION: Two-step chain reproduces final_props.  My one-step\n');
    fprintf('              must be subtly different — possibly different symmetrisation\n');
    fprintf('              or interpolation ordering.  Inspect Test C vs Test A.\n');
else
    fprintf('  CONCLUSION: No tested form matches.  Deeper investigation needed.\n');
end
