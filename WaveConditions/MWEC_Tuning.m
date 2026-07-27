%% MWEC_Tuning.m  —  Stage-1 natural-period placement pipeline  v2 (coupled)
%
%  Deliverable (the ONLY output): the three design natural periods
%       (T_n^surge, T_n^heave, T_n^pitch)   per region per closure.
%  K_PTO is a free Stage-2 variable — never an output of this pipeline.
%
%  Architecture (mirrors WEC_GM):
%    MWEC_Tuning.m          — driver; all hardcoded values in §0.
%    MWEC_Tuning_Kernels.m  — computational static methods (coupled engine).
%    MWEC_Tuning_Plots.m    — figure-generation static methods.
%
%  Methodology (v2 — see MWEC_TUNING_FRAMEWORK_PLAN.md):
%    T1  Climate spectrum S_ew + IEC moments + 90% band + partition.
%    T2  BEM at realised CG: full 3x3 A,B (off-diagonals A15,B15 retained)
%        + complex Fe.  Per-region body (own mass + cache).
%    T3  Per-mode PTO closure (Evans or 40%-critical).
%    T4  Coupled 3x3 placement: maximise eta_Falnes over the three periods.
%        Surge-pitch separation emerges from the coupled objective; the band
%        partition is the structural backstop.  No Picard, no Lorentzian.
%    T5  Per-mode capture width CW_k(T) vs the combined-dipole Falnes ceiling.
%
%  Idealised diagonal per-mode PTO: the only surge-pitch coupling in-model is
%  the intrinsic hydrodynamic A15,B15.  Tether geometry is a Stage-2 concern.

clear; clc; close all;

%% ============================================================================
%%  §0  CONFIGURATION  (all hardcoded values live here)
%% ============================================================================

cfg = struct();

%% §0.1  Paths
this_dir            = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
cfg.this_dir        = this_dir;
cfg.body_dir        = this_dir;                       % where region body files live
cfg.climate_dir     = fullfile(this_dir, 'WIS_Output_WAM');
cfg.out_dir         = fullfile(this_dir, 'WaveConditions_Results');
cfg.fig_dir         = fullfile(cfg.out_dir, 'figs');
if ~exist(cfg.out_dir, 'dir'), mkdir(cfg.out_dir); end
if ~exist(cfg.fig_dir, 'dir'), mkdir(cfg.fig_dir); end

%% §0.2  CASE TABLE  (EDIT HERE — one row per region, body = mass + cache)
%   region    : label used in output filenames and figures
%   station   : WIS climate-grid id  -> WIS_Output_WAM/ST<...>_climate_grid.mat
%   mass      : per-region mass file (final_props + results.config)
%   cache     : per-region WAMIT cache ('' -> read hydro_table from mass file)
cfg.cases = struct( ...
    'region',  {'Pacific',                     'NA',                          'SouthPass'}, ...
    'station', {'ST84040',                     'ST63044',                     'ST73135'}, ...
    'mass',    {'Pacific_WEC_UHPC_MASS.mat',   'NA_WEC_UHPC_MASS.mat',        'SouthPass_WEC_UHPC_MASS.mat'}, ...
    'cache',   {'Pacific_C1_wamit_cache.mat',  'NA_C1_wamit_cache.mat',       'SouthPass_C1_wamit_cache.mat'});
cfg.cases_to_run    = {'Pacific', 'NA', 'SouthPass'};   % subset selector

%% §0.3  Climate band
cfg.bandwidth.pct_low      = 0.05;
cfg.bandwidth.pct_high     = 0.95;
cfg.bandwidth.T_min_s      = 3.0;            % device-physics floor (s); 0 disables

%% §0.4  Stage-1 closure mode  (choose ONE; or compare both via §0.5)
%   'evans'            B_PTO,k = B_kk(omega_n,k)                  (resistance match)
%   'critical_damping' B_PTO,k = zeta*2*sqrt(K_total,k*(M_k+A_inf,k))  per mode
cfg.closure.mode           = 'evans';
cfg.closure.zeta           = 0.40;           % only for critical_damping

%% §0.5  Run both closures back-to-back for the cross-pipeline comparison?
cfg.compare_closures       = true;
cfg.closure_list           = {'evans', 'critical_damping'};

%% §0.6  Optimisation
cfg.opt.ga.PopulationSize  = 80;
cfg.opt.ga.MaxGenerations  = 120;
cfg.opt.fminsearch.TolX    = 1e-5;
cfg.opt.fmincon.TolFun     = 1e-10;
cfg.opt.verbose            = true;

%% §0.7  Absorption diagnostic
cfg.absorption.T_grid_s    = (2:0.1:14).';
cfg.absorption.L_ref       = [];             % [] = auto (cross_section beam)

%% §0.8  Validation gate tolerances
cfg.gate.m0_pct_max        = 5.0;            % G1
cfg.gate.T_pitch_pct_max   = 2.0;            % G2c
cfg.gate.fp_tol            = 1e-6;           % G2b
cfg.gate.fp_maxit          = 50;
cfg.gate.l_falnes_pct_max  = 1.0;            % G5a

%% §0.9  Plot output
cfg.plot.enable            = true;
cfg.plot.formats           = {'pdf', 'png', 'fig'};

%% §0.10 Random seed (reproducibility)
cfg.rng_seed               = 0;

%% ============================================================================
%%  §1  RUN  (loop over closures, then over cases)
%% ============================================================================

sep = repmat('=', 1, 72);
fprintf('\n%s\n  MWEC_Tuning v2 — Stage-1 coupled placement\n%s\n', sep, sep);

if cfg.compare_closures
    closures = cfg.closure_list;
else
    closures = {cfg.closure.mode};
end

% master container: results.(closureTag).(region)
master = struct();

for ic = 1:numel(closures)
    closure_mode = closures{ic};
    closure_cfg  = struct('mode', closure_mode, 'zeta', cfg.closure.zeta);
    ctag = matlab.lang.makeValidName(closure_mode);

    fprintf('\n%s\n  CLOSURE: %s', sep, closure_mode);
    if strcmpi(closure_mode, 'critical_damping'), fprintf('  (zeta = %.2f)', cfg.closure.zeta); end
    fprintf('\n%s\n', sep);

    summary = struct();  sidx = 0;

    for k = 1:numel(cfg.cases)
        C = cfg.cases(k);
        if ~any(strcmp(cfg.cases_to_run, C.region)), continue; end

        fprintf('\n%s\n  CASE: %s  (station %s)  body: %s\n%s\n', ...
                sep, C.region, C.station, C.mass, sep);
        t_case = tic;   % sidx incremented only once a case actually completes (§7)

        %% --- §1.1  Load per-region body (mass + cache) -------------------
        mass_path = fullfile(cfg.body_dir, C.mass);
        if ~isfile(mass_path)
            warning('MWEC_Tuning:noMass', 'Body file missing: %s — skipping case %s.', mass_path, C.region);
            continue;
        end
        S = load(mass_path);
        assert(isfield(S, 'final_props'), 'Mass file %s missing final_props.', C.mass);
        final_props = S.final_props;
        assert(isfield(S, 'results') && isfield(S.results, 'config'), 'Mass file %s missing results.config.', C.mass);
        config = S.results.config;
        assert(isfield(config, 'RHO_WATER') && isfield(config, 'G'), ...
            'Mass-file %s config missing RHO_WATER / G.', C.mass);
        rho   = config.RHO_WATER;
        g_acc = config.G;

        % hydro_table: external cache if present, else embedded in mass file
        hydro_table = local_load_hydro(cfg.body_dir, C.cache, config, C.region);
        fprintf('  [§1] rho=%.2f kg/m^3  g=%.5f m/s^2\n', rho, g_acc);

        %% --- §2  T2 — BEM at realised CG (per region) --------------------
        fprintf('  [§2] T2 — BEM at realised CG (full 3x3 + coupling diagnostic)\n');
        bem = MWEC_Tuning_Kernels.bem_at_design(final_props, hydro_table, cfg.gate, true);
        [L_ref, L_ref_src] = MWEC_Tuning_Kernels.resolve_L_ref(final_props, cfg.absorption.L_ref);
        fprintf('       L_ref=%.3f m  (%s)\n', L_ref, L_ref_src);

        %% --- §3  T1 — Climate spectrum + band + partition ----------------
        cg_path = fullfile(cfg.climate_dir, [C.station '_climate_grid.mat']);
        if ~isfile(cg_path)
            warning('MWEC_Tuning:noClimate', 'Climate file missing: %s — skipping.', cg_path);
            continue;
        end
        cg = load(cg_path);
        assert(isfield(cg, 'climateGrid'), 'Climate grid %s missing climateGrid.', cg_path);
        climateGrid = cg.climateGrid;
        omega = climateGrid.omega(:);

        fprintf('  [§3] T1 — S_ew + IEC moments\n');
        [F_ew, S_ew, sew_meta] = MWEC_Tuning_Kernels.build_S_ew( ...
            climateGrid, omega, rho, g_acc, cfg.gate.m0_pct_max);
        IEC  = MWEC_Tuning_Kernels.compute_iec_moments(omega, S_ew);
        band = MWEC_Tuning_Kernels.energy_band(omega, S_ew, ...
            cfg.bandwidth.pct_low, cfg.bandwidth.pct_high, cfg.bandwidth.T_min_s);
        part = MWEC_Tuning_Kernels.band_partition(omega, S_ew, band.omega_L, band.omega_H);
        omega_bounds = [band.omega_L,        part.partition_omega; ...
                        band.omega_L,        band.omega_H; ...
                        part.partition_omega, band.omega_H];
        fprintf('       m0=%.5f m^2 (G1 err %.2f%%)  Hm0=%.3f m  Te=%.2f s  Tp=%.2f s  eps=%.3f\n', ...
                sew_meta.m0_Sew, sew_meta.m0_err_pct, IEC.Hm0, IEC.Te, IEC.Tp, IEC.eps_bw);
        fprintf('       T band [%.2f, %.2f] s   partition T=%.2f s   %s\n', ...
                band.T_L, band.T_H, 2*pi/part.partition_omega, ternary(part.bimodal, 'BIMODAL', 'unimodal'));
        if bem.omega_BEM(1) > band.omega_L || bem.omega_BEM(end) < band.omega_H
            warning('MWEC_Tuning:G2d', 'G2d WARN: BEM omega [%.3f,%.3f] does not cover band [%.3f,%.3f].', ...
                    bem.omega_BEM(1), bem.omega_BEM(end), band.omega_L, band.omega_H);
        end

        %% --- §4  T3/T4 — Closure rule + coupled placement ----------------
        fprintf('  [§4] T3 — closure: %s\n', closure_cfg.mode);
        fprintf('       T4 — coupled placement (maximise eta_Falnes over the 3 periods)\n');
        rng(cfg.rng_seed, 'twister');
        [omega_n_opt, eta_opt, place_info] = MWEC_Tuning_Kernels.place_natural_periods( ...
            bem, closure_cfg, omega, S_ew, band, omega_bounds, rho, g_acc, cfg.opt, 'freq');
        res = MWEC_Tuning_Kernels.evaluate_placement( ...
            omega_n_opt, bem, closure_cfg, omega, S_ew, band, rho, g_acc, 'freq');
        T_n_opt = 2*pi ./ omega_n_opt;
        CWR     = res.CW / max(L_ref, eps) * 100;

        fprintf('       eta_Falnes = %.4f   (%% of band-capturable energy)\n', eta_opt);
        fprintf('       T_n design = [%.3f %.3f %.3f] s   (surge heave pitch)\n', T_n_opt(1), T_n_opt(2), T_n_opt(3));
        fprintf('       <P_abs>    = [%.1f %.1f %.1f] W   total %.1f W\n', ...
                res.P_abs_per_mode(1), res.P_abs_per_mode(2), res.P_abs_per_mode(3), res.P_abs_total);
        fprintf('       CW / CWR   = %.3f m / %.1f %%   (L_ref=%.3f m)\n', res.CW, CWR, L_ref);
        fprintf('       B_PTO/mode = [%.3e %.3e %.3e]   surge-pitch: %s\n', ...
                res.closure.B_PTO_per_mode(1), res.closure.B_PTO_per_mode(2), res.closure.B_PTO_per_mode(3), ...
                ternary(place_info.sep_binding, 'partition BINDING', 'separated (partition slack)'));

        %% --- §5  T5 — Capture-width diagnostic ---------------------------
        fprintf('  [§5] T5 — per-mode capture width vs Falnes ceiling\n');
        absorption = MWEC_Tuning_Kernels.absorption_vs_T( ...
            cfg.absorption.T_grid_s, bem, res.closure, omega_n_opt, rho, g_acc, L_ref, cfg.gate, 'freq');
        fprintf('       max Falnes-ceiling utilisation overrun = %.2f%% (gate %.1f%%)\n', ...
                absorption.max_viol_pct, cfg.gate.l_falnes_pct_max);

        %% --- §6  A_inf sensitivity ---------------------------------------
        fprintf('  [§6] A_inf sensitivity — re-place with A(omega) -> A_inf\n');
        rng(cfg.rng_seed, 'twister');
        A_inf_comp = MWEC_Tuning_Kernels.sensitivity_A_inf( ...
            bem, closure_cfg, omega, S_ew, band, omega_bounds, rho, g_acc, cfg.opt, L_ref, cfg.gate);
        if isfield(A_inf_comp, 'error')
            fprintf('       SENSITIVITY FAILED: %s\n', A_inf_comp.error);
        else
            fprintf('                       %-12s | %-12s | %-9s\n', 'A(omega)', 'A_inf', 'Delta');
            lbls = {'T_n,surge (s)', 'T_n,heave (s)', 'T_n,pitch (s)'};
            for m = 1:3
                fprintf('       %-15s %12.3f | %12.3f | %+8.2f%%\n', lbls{m}, ...
                        T_n_opt(m), A_inf_comp.T_n_placed_s(m), 100*(A_inf_comp.T_n_placed_s(m)/T_n_opt(m)-1));
            end
            fprintf('       %-15s %12.4f | %12.4f | %+8.2f%%\n', 'eta_Falnes', ...
                    eta_opt, A_inf_comp.eta_Falnes, 100*(A_inf_comp.eta_Falnes/eta_opt-1));
            fprintf('       %-15s %12.1f | %12.1f | %+8.2f%%\n', '<P_abs> tot (W)', ...
                    res.P_abs_total, A_inf_comp.P_abs_total, 100*(A_inf_comp.P_abs_total/res.P_abs_total-1));
        end

        %% --- §7  Pack + save ---------------------------------------------
        results = struct();
        results.region            = C.region;
        results.station_id        = C.station;
        results.closure_mode      = closure_cfg.mode;
        results.T_n_target_s      = T_n_opt;           % THE deliverable
        results.omega_n_target    = omega_n_opt;
        results.eta_Falnes        = eta_opt;
        results.absorbed_power_W       = res.P_abs_per_mode;
        results.absorbed_power_W_total = res.P_abs_total;
        results.CW                = res.CW;
        results.CWR               = CWR;
        results.closure           = res.closure;
        results.place_info        = place_info;
        results.absorption_vs_T   = absorption;
        results.A_inf_comparison  = A_inf_comp;
        results.L_ref             = L_ref;
        results.L_ref_source      = L_ref_src;
        results.bem               = bem;
        results.climate           = struct('omega', omega, 'S_ew', S_ew, 'F_ew', F_ew, ...
                                            'band', band, 'partition', part, ...
                                            'omega_bounds', omega_bounds, 'm0_meta', sew_meta, 'IEC', IEC);
        results.cfg               = cfg;
        results.meta              = struct('climate_grid_path', cg_path, 'mass_path', mass_path, ...
                                           'rho', rho, 'g_acc', g_acc, 't_elapsed_s', toc(t_case));

        out_path = fullfile(cfg.out_dir, sprintf('Tuning_results_%s_%s.mat', ctag, C.region));
        save(out_path, 'results', '-v7');
        fprintf('       Saved: %s\n', out_path);
        master.(ctag).(C.region) = results;

        sidx = sidx + 1;   % advance only for a completed case (no holes on skip)
        summary(sidx).region        = C.region;
        summary(sidx).station       = C.station;
        summary(sidx).T_n_target_s  = T_n_opt;
        summary(sidx).eta_Falnes    = eta_opt;
        summary(sidx).CW            = res.CW;
        summary(sidx).CWR           = CWR;
        summary(sidx).P_abs_total_W = res.P_abs_total;
        summary(sidx).sep_binding   = place_info.sep_binding;
        if isfield(A_inf_comp, 'error')
            summary(sidx).A_inf_eta        = NaN;
            summary(sidx).A_inf_eta_dpct   = NaN;
            summary(sidx).A_inf_Pabs_dpct  = NaN;
        else
            summary(sidx).A_inf_eta        = A_inf_comp.eta_Falnes;
            summary(sidx).A_inf_eta_dpct   = 100*(A_inf_comp.eta_Falnes/eta_opt-1);
            summary(sidx).A_inf_Pabs_dpct  = 100*(A_inf_comp.P_abs_total/res.P_abs_total-1);
        end

        %% --- §8  Per-case figures ----------------------------------------
        if cfg.plot.enable
            try
                MWEC_Tuning_Plots.fig_t1_climate(C.region, ctag, results, cfg);
                MWEC_Tuning_Plots.fig_t2_bem(C.region, ctag, results, cfg);
                MWEC_Tuning_Plots.fig_t4_placement(C.region, ctag, results, cfg);
                MWEC_Tuning_Plots.fig_t5_absorption(C.region, ctag, results, cfg);
                MWEC_Tuning_Plots.fig_t6_scatter(C.region, ctag, results, climateGrid, cfg);
                MWEC_Tuning_Plots.fig_c1_interference(C.region, ctag, results, cfg);
            catch ME_p
                warning('MWEC_Tuning:plotFail', 'Per-case plotting failed (%s): %s', C.region, ME_p.message);
            end
        end

        fprintf('  Case %s [%s] done (%.1f s)\n', C.region, closure_cfg.mode, toc(t_case));
    end

    %% --- Cross-case summary for this closure ----------------------------
    fprintf('\n%s\n  CROSS-CASE SUMMARY  (closure = %s)\n%s\n', sep, closure_mode, sep);
    fprintf('  %-10s | %9s | %8s | %8s | %8s | %7s | %10s\n', ...
            'Region', 'eta_Faln', 'T_surge', 'T_heave', 'T_pitch', 'CWR %', '<P_abs> W');
    fprintf('  %s\n', repmat('-', 1, 78));
    for s = 1:numel(summary)
        if ~isfield(summary(s), 'region') || isempty(summary(s).region), continue; end
        r = summary(s);
        fprintf('  %-10s | %9.4f | %8.3f | %8.3f | %8.3f | %7.1f | %10.1f\n', ...
                r.region, r.eta_Falnes, r.T_n_target_s(1), r.T_n_target_s(2), r.T_n_target_s(3), r.CWR, r.P_abs_total_W);
    end
    fprintf('%s\n', sep);

    master_summary.(ctag) = summary; %#ok<SAGROW>
end

%% ============================================================================
%%  §2  SAVE + CROSS-PIPELINE FIGURES
%% ============================================================================

save(fullfile(cfg.out_dir, 'Tuning_summary.mat'), 'master', 'master_summary', 'cfg', '-v7');
fprintf('\n  Saved cross-pipeline summary: %s\n', fullfile(cfg.out_dir, 'Tuning_summary.mat'));

% Evans vs Critical-damping comparison (CWR metric) — only if both ran
local_evans_vs_crit(master_summary, cfg, sep);

if cfg.plot.enable
    try
        MWEC_Tuning_Plots.fig_t7_cross(master_summary, cfg);
        MWEC_Tuning_Plots.fig_t8_closure_compare(master_summary, cfg);
    catch ME_p
        warning('MWEC_Tuning:plotFail', 'Cross-pipeline plotting failed: %s', ME_p.message);
    end
end

fprintf('\n%s\n  MWEC_Tuning v2 complete\n%s\n\n', sep, sep);


%% ============================================================================
%%  Local utilities (driver-only)
%% ============================================================================
function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

function hydro_table = local_load_hydro(body_dir, cache_name, config, region)
%LOCAL_LOAD_HYDRO  Resolve hydro_table from external cache or embedded copy.
    if ~isempty(cache_name)
        cache_path = fullfile(body_dir, cache_name);
        if isfile(cache_path)
            WC = load(cache_path);
            if isfield(WC, 'hydro_table') && isstruct(WC.hydro_table)
                hydro_table = WC.hydro_table;  return;
            elseif isfield(WC, 'results') && isfield(WC.results, 'config') && isfield(WC.results.config, 'hydro_table')
                hydro_table = WC.results.config.hydro_table;  return;
            else
                error('MWEC_Tuning:badCache', 'Cache %s missing hydro_table.', cache_name);
            end
        end
    end
    % Fall back to embedded copy in the mass file's config
    if isfield(config, 'hydro_table')
        hydro_table = config.hydro_table;
    elseif isfield(config, 'hydro_cache')
        hydro_table = config.hydro_cache;
    else
        error('MWEC_Tuning:noHydro', ...
              'No external cache (%s) and no embedded hydro_table for region %s.', cache_name, region);
    end
end

function local_evans_vs_crit(master_summary, cfg, sep)
%LOCAL_EVANS_VS_CRIT  Console comparison of the two closures, keyed on CWR.
    if ~(isfield(master_summary, 'evans') && isfield(master_summary, 'critical_damping'))
        return;   % both pipelines must have run
    end
    SE = master_summary.evans;  SC = master_summary.critical_damping;
    fprintf('\n%s\n  EVANS vs CRITICAL-DAMPING   (deliverable = design periods; metric = CWR)\n%s\n', sep, sep);
    fprintf('  %-10s | %-19s | %-19s | %7s | %7s | %9s\n', ...
            'Region', 'T_n Evans (s)', 'T_n Crit (s)', 'CWR_E %', 'CWR_C %', 'd(C-E)pp');
    fprintf('  %s\n', repmat('-', 1, 86));
    cwrE = []; cwrC = [];
    for k = 1:numel(cfg.cases)
        reg = cfg.cases(k).region;
        if ~any(strcmp(cfg.cases_to_run, reg)), continue; end
        e = local_find_region(SE, reg);  c = local_find_region(SC, reg);
        if isempty(e) || isempty(c), continue; end
        cwrE(end+1) = e.CWR;  cwrC(end+1) = c.CWR; %#ok<AGROW>
        fprintf('  %-10s | %5.2f/%5.2f/%5.2f | %5.2f/%5.2f/%5.2f | %7.1f | %7.1f | %+9.1f\n', ...
                reg, e.T_n_target_s(1), e.T_n_target_s(2), e.T_n_target_s(3), ...
                c.T_n_target_s(1), c.T_n_target_s(2), c.T_n_target_s(3), ...
                e.CWR, c.CWR, c.CWR - e.CWR);
    end
    fprintf('  %s\n', repmat('-', 1, 86));
    if ~isempty(cwrE)
        fprintf('  %-10s | %19s | %19s | %7.1f | %7.1f | %+9.1f\n', ...
                'mean', '', '', mean(cwrE), mean(cwrC), mean(cwrC) - mean(cwrE));
        names = {'Evans', 'Critical-damping'};
        [~, wi] = max([mean(cwrE), mean(cwrC)]);
        fprintf('  Higher mean CWR: %s   (d(C-E) in percentage points)\n', names{wi});
    end
    fprintf('%s\n', sep);
end

function e = local_find_region(summary, region)
%LOCAL_FIND_REGION  Return the summary entry for a region, or [] if absent.
    e = [];
    for i = 1:numel(summary)
        if isfield(summary(i), 'region') && ~isempty(summary(i).region) && strcmp(summary(i).region, region)
            e = summary(i);  return;
        end
    end
end
