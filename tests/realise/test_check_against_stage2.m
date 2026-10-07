function test_check_against_stage2()
%TEST_CHECK_AGAINST_STAGE2  Contract F10 on hand-made property sets (no kernel).
%   Values are binary fractions, so every relative deviation and residual below is exact and the
%   limits are tested at equality.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'src'));
f10 = @mwecmass.realise.check_against_stage2;
s2 = struct('vs', 0.5, 'rho', [2500; 1000], 'mass', 2048, 'Z_CG', -2, 'GM', 0.25, ...
    'T_heave', 4, 'T_pitch', 8);
rho_w = 1024;
props = struct('CG_total', [0 0 -2.25], 'GM_L', 0.25, 'periods', struct('heave', 4.5, 'pitch', 7), ...
    'mass_total', 2048, 'V_sub', 2);

% rel_dev: Z_CG 0.25/2, GM 0, T_heave 0.5/4, T_pitch 1/8 (all 0.125 or 0); limit 12.5 % holds at equality
ck = f10(props, s2, 12.5, 1e-6, rho_w);
names = {ck.metrics.name};
check(isequal(names, {'Z_CG', 'GM', 'T_heave', 'T_pitch'}), 'metric names');
check(isequal([ck.metrics.rel_dev], [0.125 0 0.125 0.125]), 'rel_dev');
check(isequal([ck.metrics.value], [-2.25 0.25 4.5 7]) && isequal([ck.metrics.stage2], [-2 0.25 4 8]), 'values');
check(all([ck.metrics.limit] == 0.125) && all([ck.metrics.pass]), 'limit inclusive');
check(isequal({ck.equalities.name}, {'flotation', 'GM'}) && isequal([ck.equalities.residual], [0 0]), ...
    'equality residuals');
check(ck.pass && isempty(ck.failed) && isempty(ck.reason), 'accepted case');
fprintf('accepted case: rel_dev %s, residuals %s\n', mat2str([ck.metrics.rel_dev]), ...
    mat2str([ck.equalities.residual]));

% one ulp below the limit fails only the metrics at 0.125
ck = f10(props, s2, 12.5 * (1 - eps), 1e-6, rho_w);
check(isequal(ck.failed, {'Z_CG', 'T_heave', 'T_pitch'}) && ~ck.pass, 'failed metrics list');
check(~isempty(strfind(ck.reason, 'T_pitch')) && isempty(strfind(ck.reason, 'flotation')), 'reason text');

% flotation residual 2^-10 against tol 2^-11 fails; at tol 2^-10 it holds (|r| <= tol)
p = props;
p.mass_total = 2050;
ck = f10(p, s2, 50, 2^-11, rho_w);
check(ck.equalities(1).residual == 2^-10 && ~ck.equalities(1).pass, 'flotation residual');
check(~ck.pass && isequal(ck.failed, {'flotation'}) && ~isempty(strfind(ck.reason, 'flotation')), ...
    'flotation enters pass');
ck = f10(p, s2, 50, 2^-10, rho_w);
check(ck.pass && ck.equalities(1).pass, 'flotation at tolerance');

% the GM equality is reported but is not part of pass (AGENTS section 3 item 32)
p = props;
p.GM_L = 0.25 * (1 + 2^-4);
ck = f10(p, s2, 12.5, 1e-6, rho_w);
check(ck.equalities(2).residual == 2^-4 && ~ck.equalities(2).pass, 'GM equality residual');
check(ck.pass && isempty(ck.failed), 'GM equality outside pass');

% a non-finite realised value never passes
p = props;
p.periods.pitch = Inf;
ck = f10(p, s2, 12.5, 1e-6, rho_w);
check(~ck.pass && isequal(ck.failed, {'T_pitch'}), 'Inf period fails');
p.periods.pitch = NaN;
ck = f10(p, s2, 12.5, 1e-6, rho_w);
check(~ck.pass && isequal(ck.failed, {'T_pitch'}), 'NaN period fails');
fprintf('failure cases reported: %s\n', ck.reason);
end

function check(cond, msg)
if ~cond
    error('test_check_against_stage2:fail', '%s', msg);
end
end
