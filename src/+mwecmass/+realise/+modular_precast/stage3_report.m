function stage3_report(realised, notes)
%STAGE3_REPORT  Print the realised modular-precast design and its Stage-3 report.
%
%   mwecmass.realise.modular_precast.stage3_report(realised, notes)
%
%   realised: S8 (modular precast); notes: cellstr of the split and escalation notes. Prints the
%   modules (shell thickness, ballast height from the module bottom, UHPC and air volumes against
%   the split targets, densities), the ballast level and mass, every metric of the F10 check
%   (value, Stage-2 value, deviation, limit, pass), the equality residuals (flotation enters the
%   check; GM is the solver's equality, reported beside it), every escalation step run (exit flag,
%   iterations, objective, largest equality residual), the notes and the status with its reason.

r = realised;
fprintf('      realised modules (body frame), escalation %s, vs = %.6f m, draft %.6f m:\n', ...
    r.escalation, r.vs, r.draft);
fprintf('      %-3s %9s %10s %12s %12s %12s %10s %10s\n', 'mod', 't[mm]', 'h_ball[m]', 'V_uhpc[m3]', ...
    'target[m3]', 'V_air[m3]', 'rho_eff', 'rho2');
for i = 1:numel(r.modules)
    m = r.modules(i);
    fprintf('      %-3d %9.3f %10.5f %12.6f %12.6f %12.6f %10.3f %10.3f\n', i, 1000 * m.t, m.h_ballast, ...
        m.V_uhpc, r.V_uhpc_target(i), m.V_air, m.rho_eff, m.rho_stage2);
end
fprintf('      z_ballast = %.6f m (body), mass %.3f kg, displaced mass %.3f kg\n', r.design.z_ballast, ...
    r.props.mass_total, r.props.mass_buoyant_force);
fprintf('      %-8s %12s %12s %10s %8s  %s\n', 'metric', 'Stage 2', 'realised', 'dev[%]', 'limit', 'pass');
for m = r.check.metrics
    fprintf('      %-8s %12.6f %12.6f %10.4f %8.2f  %d\n', m.name, m.stage2, m.value, 100 * m.rel_dev, ...
        100 * m.limit, m.pass);
end
for q = r.check.equalities
    role = 'in the check';
    if strcmp(q.name, 'GM')
        role = 'solver equality, not in the check';
    end
    fprintf('      %-9s equality residual %11.3g (tol %.3g, %s), pass %d\n', q.name, q.residual, q.tol, ...
        role, q.pass);
end
fprintf('      %-12s %9s %11s %14s %16s\n', 'step', 'exitflag', 'iterations', 'objective', 'max|eq residual|');
for s = r.solver
    fprintf('      %-12s %9g %11d %14.6g %16.3g\n', s.step, s.exitflag, s.iterations, s.fval, ...
        s.max_eq_violation);
end
for k = 1:numel(notes)
    fprintf('      note: %s\n', notes{k});
end
fprintf('      Stage 3 status: %s%s\n', r.status, ...
    mwecmass.internal.ternary(isempty(r.reason), '', [' (' r.reason ')']));
end
