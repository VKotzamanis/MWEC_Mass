function check = check_against_stage2(props, stage2, pct, tol_eq, rho_water)
%CHECK_AGAINST_STAGE2  Compare a realised Stage-3 design with its Stage-2 solution (contract F10).
%
%   check = mwecmass.realise.check_against_stage2(props, stage2, pct, tol_eq, rho_water)
%
%   props: final_props of the realised design (evaluate_realised). stage2: the S8 stage2 struct
%   (Z_CG = Final3D.CG_total(3) in the world frame, GM = Final3D.GM_L, coupled T_heave and
%   T_pitch). pct: in.pid.mass_acceptable_pct [%]. tol_eq: the constraint tolerance the Stage-3
%   equalities are held to. rho_water [kg/m^3].
%
%   metrics(k), k = Z_CG, GM, T_heave, T_pitch: value (realised), stage2, rel_dev =
%   |X3 - X2| / |X2|, limit = pct/100, pass = rel_dev <= limit (AGENTS section 3 items 4.4, 16).
%   equalities(k): flotation, residual M / (rho_w V_sub) - 1 (as Stage 2), and GM, residual
%   GM / GM2 - 1; tol = tol_eq, pass = |residual| <= tol.
%   pass: every metric passes and flotation holds. This is the owner's acceptance rule (AGENTS
%   section 3 item 32): mass balance satisfied, and Z_CG, GM and the periods within
%   mass_acceptable_pct. The GM equality is a constraint of the Stage-3 optimisation (item 5);
%   its row is reported here for the solver's closest-fail rule and does not enter pass.
%   failed: names of the metrics and equalities in pass that fail; reason: one sentence each.

names = {'Z_CG', 'GM', 'T_heave', 'T_pitch'};
vals = [props.CG_total(3), props.GM_L, props.periods.heave, props.periods.pitch];
ref = [stage2.Z_CG, stage2.GM, stage2.T_heave, stage2.T_pitch];
rel = abs(vals - ref) ./ abs(ref);
metrics = struct('name', names, 'value', num2cell(vals), 'stage2', num2cell(ref), ...
    'rel_dev', num2cell(rel), 'limit', pct / 100, 'pass', num2cell(rel <= pct / 100));

res = [props.mass_total / (rho_water * props.V_sub) - 1, props.GM_L / stage2.GM - 1];
equalities = struct('name', {'flotation', 'GM'}, 'residual', num2cell(res), 'tol', tol_eq, ...
    'pass', num2cell(abs(res) <= tol_eq));

parts = {};
failed = {};
for k = find(~[metrics.pass])
    failed{end + 1} = names{k}; %#ok<AGROW>
    parts{end + 1} = sprintf('%s = %.6g deviates %.3g %% from Stage 2 (%.6g; limit %.3g %%)', ...
        names{k}, vals(k), 100 * rel(k), ref(k), pct); %#ok<AGROW>
end
if ~equalities(1).pass
    failed{end + 1} = 'flotation';
    parts{end + 1} = sprintf('flotation residual M/(rho_w V_sub) - 1 = %.3g exceeds %.3g', ...
        res(1), tol_eq);
end
check = struct('metrics', {metrics}, 'equalities', {equalities}, 'pass', isempty(failed), ...
    'failed', {failed}, 'reason', strjoin(parts, '; '));
end
