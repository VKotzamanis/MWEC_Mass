%% test_tuning_multistation.m  —  smoke test for MWEC_Tuning v2 (coupled)
%
%  Sequence:
%    1. Run MWEC_Tuning.m (all configured closures x cases).
%    2. Load Tuning_summary.mat (master, master_summary).
%    3. Apply sanity gates per closure per region:
%         - eta_Falnes in (0, 1]              (bounded capture-efficiency)
%         - T_n_target_s finite, inside the placement band
%         - surge resonance < pitch resonance (partition: T_surge > T_pitch)
%         - absorbed_power_W_total > 0,  CW > 0
%         - bem.gate_G2c_pass == true        (CG-transform regression)
%         - absorption_vs_T.max_viol_pct within the Falnes-ceiling gate
%    4. Re-run once to confirm reproducibility under cfg.rng_seed.
%
%  NOTE: requires the per-region body files (Pacific/NA/SouthPass _WEC_UHPC_MASS
%  + _C0_wamit_cache).  If absent, MWEC_Tuning warns-and-skips and this test
%  reports "no cases ran" rather than failing.

clear; clc;
bar = repmat('=', 1, 72);
fprintf('\n%s\n  Smoke test: MWEC_Tuning v2\n%s\n\n', bar, bar);

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
out_dir  = fullfile(this_dir, 'WaveConditions_Results');

%% --- 1. Run -------------------------------------------------------------
run(fullfile(this_dir, 'MWEC_Tuning.m'));

%% --- 2. Load ------------------------------------------------------------
summary_path = fullfile(out_dir, 'Tuning_summary.mat');
assert(isfile(summary_path), 'Smoke test: expected %s after run.', summary_path);
Sld = load(summary_path);
master         = Sld.master;
master_summary = Sld.master_summary;

closures = fieldnames(master);
if isempty(closures) || all(structfun(@(x) isempty(fieldnames(x)), master))
    fprintf('\n  No cases ran (body files likely absent).  Test inconclusive — supply\n');
    fprintf('  Pacific/NA/SouthPass _WEC_UHPC_MASS + _C0_wamit_cache and re-run.\n%s\n', bar);
    return;
end

%% --- 3. Sanity gates ----------------------------------------------------
fprintf('\n%s\n  Sanity gates\n%s\n', repmat('-', 1, 72), repmat('-', 1, 72));
n_pass = 0; n_fail = 0;
for ic = 1:numel(closures)
    ctag = closures{ic};
    regions = fieldnames(master.(ctag));
    for ir = 1:numel(regions)
        r = master.(ctag).(regions{ir});
        issues = {};
        Tn = r.T_n_target_s;  om_n = r.omega_n_target;
        band = r.climate.band;

        if ~isfinite(r.eta_Falnes) || r.eta_Falnes <= 0 || r.eta_Falnes > 1 + 1e-9
            issues{end+1} = sprintf('eta_Falnes=%.3f outside (0,1]', r.eta_Falnes); %#ok<*AGROW>
        end
        if any(~isfinite(Tn))
            issues{end+1} = 'NaN/Inf in T_n_target_s';
        end
        if any(om_n < band.omega_L - 1e-6) || any(om_n > band.omega_H + 1e-6)
            issues{end+1} = sprintf('omega_n outside band [%.3f,%.3f]', band.omega_L, band.omega_H);
        end
        if ~(om_n(1) < om_n(3))   % surge resonance below pitch resonance (partition)
            issues{end+1} = sprintf('surge/pitch not separated: om=[%.3f %.3f %.3f]', om_n);
        end
        if r.absorbed_power_W_total <= 0 || ~isfinite(r.absorbed_power_W_total)
            issues{end+1} = sprintf('<P_abs>=%g W', r.absorbed_power_W_total);
        end
        if r.CW <= 0 || ~isfinite(r.CW)
            issues{end+1} = sprintf('CW=%g m', r.CW);
        end
        if isfield(r.bem, 'gate_G2c_pass') && ~r.bem.gate_G2c_pass
            issues{end+1} = 'G2c: pitch period mismatch';
        end
        if isfield(r.absorption_vs_T, 'max_viol_pct') && ...
                r.absorption_vs_T.max_viol_pct > r.cfg.gate.l_falnes_pct_max + 1e-6
            issues{end+1} = sprintf('G5a: CW exceeds ceiling by %.2f%%', r.absorption_vs_T.max_viol_pct);
        end

        if isempty(issues)
            fprintf('    [%s] %-10s : PASS  (eta=%.4f, CWR=%.1f%%, <P>=%.0f W)\n', ...
                ctag, regions{ir}, r.eta_Falnes, r.CWR, r.absorbed_power_W_total);
            n_pass = n_pass + 1;
        else
            fprintf('    [%s] %-10s : FAIL  %s\n', ctag, regions{ir}, strjoin(issues, ' | '));
            n_fail = n_fail + 1;
        end
    end
end
fprintf('\n  %d/%d (closure x region) cases passed.\n', n_pass, n_pass + n_fail);

%% --- 4. Reproducibility -------------------------------------------------
fprintf('\n%s\n  Reproducibility check (re-run, same cfg.rng_seed)\n%s\n', repmat('-', 1, 72), repmat('-', 1, 72));
snap = struct();
for ic = 1:numel(closures)
    ctag = closures{ic};
    regions = fieldnames(master.(ctag));
    for ir = 1:numel(regions)
        r = master.(ctag).(regions{ir});
        snap.(ctag).(regions{ir}) = struct('T', r.T_n_target_s, 'eta', r.eta_Falnes, 'P', r.absorbed_power_W_total);
    end
end

run(fullfile(this_dir, 'MWEC_Tuning.m'));
S2 = load(summary_path);  master2 = S2.master;

max_diff = 0; n_repro = 0; n_drift = 0;
for ic = 1:numel(closures)
    ctag = closures{ic};
    regions = fieldnames(master.(ctag));
    for ir = 1:numel(regions)
        a = snap.(ctag).(regions{ir});
        b = master2.(ctag).(regions{ir});
        d = max([max(abs(a.T - b.T_n_target_s)), abs(a.eta - b.eta_Falnes)]);
        if d > max_diff, max_diff = d; end
        if d < 1e-6
            fprintf('    [%s] %-10s : reproducible (max diff %.2e)\n', ctag, regions{ir}, d);
            n_repro = n_repro + 1;
        else
            fprintf('    [%s] %-10s : DRIFT %.2e\n', ctag, regions{ir}, d);
            n_drift = n_drift + 1;
        end
    end
end
fprintf('\n  %d/%d reproducible.  Max headline drift: %.2e\n%s\n', n_repro, n_repro + n_drift, max_diff, bar);
