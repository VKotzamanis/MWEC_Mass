function check = check_against_stage2(props, stage2, pct, tol_eq, rho_water)
%CHECK_AGAINST_STAGE2  Stand-in of contract F10: Stage-3 metrics and equalities against Stage 2.
%
%   check = mwecmass.realise.check_against_stage2(props, stage2, pct, tol_eq, rho_water)
%
%   props must come from the evaluate_realised stand-in (props.analytic set), else
%   mwecmass:standin:NotAnalytic. Formulas of contract F10: for Z_CG = CG_total(3), GM = GM_L,
%   T_heave = periods.heave, T_pitch = periods.pitch: rel_dev = |X3 - X2| / |X2|, limit = pct/100,
%   pass = rel_dev <= limit (AGENTS section 3 items 4.4 and 16); equalities flotation
%   M / (rho_w V_sub) - 1 and GM / GM2 - 1, pass = |residual| <= tol_eq.

if ~isstruct(props) || ~isfield(props, 'analytic') || isempty(props.analytic)
    error('mwecmass:standin:NotAnalytic', 'check_against_stage2 stand-in: props.analytic is empty (not stand-in output)');
end
names = {'Z_CG', 'GM', 'T_heave', 'T_pitch'};
vals = [props.CG_total(3), props.GM_L, props.periods.heave, props.periods.pitch];
ref = [stage2.Z_CG, stage2.GM, stage2.T_heave, stage2.T_pitch];
metrics = struct('name', names, 'value', num2cell(vals), 'stage2', num2cell(ref), ...
    'rel_dev', num2cell(abs(vals - ref) ./ abs(ref)), 'limit', pct / 100, 'pass', false);
for k = 1:numel(metrics)
    metrics(k).pass = metrics(k).rel_dev <= metrics(k).limit;
end
res = [props.mass_total / (rho_water * props.V_sub) - 1, props.GM_L / stage2.GM - 1];
equalities = struct('name', {'flotation', 'GM'}, 'residual', num2cell(res), 'tol', tol_eq, 'pass', false);
for k = 1:numel(equalities)
    equalities(k).pass = abs(equalities(k).residual) <= tol_eq;
end
failed = [names(~[metrics.pass]), {equalities(~[equalities.pass]).name}];
parts = {};
for k = find(~[metrics.pass])
    parts{end + 1} = sprintf('%s deviates %.3g %% from Stage 2 (limit %.3g %%)', names{k}, ...
        100 * metrics(k).rel_dev, pct); %#ok<AGROW>
end
for k = find(~[equalities.pass])
    parts{end + 1} = sprintf('%s residual %.3g exceeds %.3g', equalities(k).name, ...
        equalities(k).residual, tol_eq); %#ok<AGROW>
end
check = struct('metrics', {metrics}, 'equalities', {equalities}, 'pass', isempty(failed), ...
    'failed', {failed}, 'reason', strjoin(parts, '; '));
end
