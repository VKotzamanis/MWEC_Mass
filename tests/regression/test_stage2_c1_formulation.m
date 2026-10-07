function test_stage2_c1_formulation()
%TEST_STAGE2_C1_FORMULATION Stage 1 and Stage 2 on C1 (no realisation): check the Stage-2 formulation
%   at the returned point. Per mode: both starts are logged, the kept start has the lower objective
%   among the starts that end within the solver's constraint tolerance, the result lies inside the
%   mode's bounds (upper bound = solid density of the mode's material), the inequality vector holds the
%   GM floor and the density-ratio pairs only, and when the solver reports success (exitflag > 0) the
%   GM floor and the flotation equality hold at the returned point to that tolerance. The values
%   printed are Octave's; they are not compared with MATLAB v1.0 numbers.
%   Environment variable TESTS_STAGE2_MODES: comma separated modes, default thin_shell,modular_precast.
%   The first run of a mode builds the geometry cache (Output/cache) and takes long; tests/run_tests.m
%   runs this file only with MWEC_REGRESSION=1.
  root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  addpath(fullfile(root, 'tests', 'octave_shims'), fullfile(root, 'src'), root);
  warning('off', 'Octave:shadowed-function');
  selected = getenv('TESTS_STAGE2_MODES');
  if isempty(selected), selected = 'thin_shell,modular_precast'; end
  modes = strsplit(selected, ',');
  constraint_tol = 1e-8;   % ConstraintTolerance set in mwecmass.optim.run

  for m = 1:numel(modes)
    mode = modes{m};
    in = WEC_User_Input();
    in.materials.realisation_type = mode;
    if strcmp(mode, 'modular_precast')
      in.materials.modular_precast.n_sub = 101;   % the n_sub = 100 grid hits a roundoff-sized area in the circle-based floors (tools/baseline_run.m)
    end
    hydro_table = mwecmass.bem.load_hydro_cache(fullfile(root, 'tests', 'fixtures', 'C1_wamit_cache_v5.mat'));
    mesh_sizing = struct();
    [mesh_sizing.mesh_Nu, mesh_sizing.mesh_Nv, ~] = mwecmass.bem.wamit.restore_mesh_sizing_from_cache(hydro_table, in);
    t0 = tic;
    config = mwecmass.driver.build_config(in, hydro_table, mesh_sizing);
    out = WEC_Output_Options();
    out.console_echo = false;
    out.save = set_all(out.save, false);
    config.output = out;
    t_config = toc(t0);

    t0 = tic;
    [opt_results, x_opt] = mwecmass.optim.run(config);
    t_optim = toc(t0);
    s2 = opt_results.stage2_3d;

    % starts
    starts = s2.starts;
    if numel(starts) ~= 2 || sum([starts.kept]) ~= 1
      error('%s: expected 2 starts with one kept, got %d with %d kept', mode, numel(starts), sum([starts.kept]));
    end
    feasible = [starts.constrviolation] <= constraint_tol;
    fvals = [starts.fval];
    if any(feasible)
      cand = find(feasible);
      [~, j] = min(fvals(cand));
      expected_kept = cand(j);
    else
      [~, expected_kept] = min([starts.constrviolation]);
    end
    kept = find([starts.kept]);
    if kept ~= expected_kept || s2.fval ~= starts(kept).fval || ~isequal(s2.x_optimal(:), starts(kept).x(:))
      error('%s: kept start %d does not follow the rule (expected %d)', mode, kept, expected_kept);
    end
    if ~isequal(x_opt(:), s2.x_optimal(:))
      error('%s: second output of optim.run is not the Stage-2 optimum', mode);
    end

    % the bottom-filled start is the one the function builds at the Stage-1 draft
    [lb, ub] = mwecmass.optim.stage2_bounds(config);
    x1 = opt_results.stage1_2d.x_optimal;
    x_bf = mwecmass.optim.stage2_bottom_filled_start(x1(1), config, lb, ub);
    if ~isequal(starts(2).x0(:), x_bf(:)) || ~isequal(starts(1).x0(:), x1(:))
      error('%s: logged start vectors differ from the Stage-1 result and the bottom-filled start', mode);
    end
    p_bf = mwecmass.hydrostatics.properties_3d(x_bf, config);

    % bounds of the mode
    x = s2.x_optimal(:)';
    if strcmp(mode, 'thin_shell'), solid = config.rho_ballast; else, solid = config.constructability_rho_hull; end
    ub_modules = ub(2:end);
    free = ub_modules > lb(2:end);
    if any(ub_modules(free) ~= solid)
      error('%s: upper bound %s, solid density of the mode %g', mode, mat2str(unique(ub_modules(free))), solid);
    end
    if any(x < lb - 1e-9 * abs(lb)) || any(x > ub + 1e-9 * abs(ub))
      error('%s: result outside the bounds', mode);
    end

    % the formulation at the returned point
    [c, ceq] = mwecmass.optim.stage2_constraints(x, config);
    w = config.wall_strip_index;
    N = config.num_ballast_sections;
    n_pairs = 0;
    for i = 1:N-1
      if isempty(w) || (i ~= w && i + 1 ~= w), n_pairs = n_pairs + 1; end
    end
    if numel(c) ~= 1 + n_pairs
      error('%s: %d inequality entries, expected 1 GM + %d ratio pairs', mode, numel(c), n_pairs);
    end
    props = mwecmass.hydrostatics.properties_3d(x, config);
    violation = max([0; c(:); abs(ceq)]);
    if abs(violation - s2.output.constrviolation) > 1e-12
      error('%s: constraint violation %.3g differs from the solver''s %.3g', mode, violation, s2.output.constrviolation);
    end
    fprintf('\n%s: config %.0f s, Stage 1 + 2 %.0f s\n', mode, t_config, t_optim);
    fprintf('  Stage 1: vs = %.4f m, rho = [%s]\n', x1(1), sprintf(' %.1f', x1(2:end)));
    for k = 1:2
      fprintf('  start %-15s f(x0) = %-10.5g fval = %-10.5g exitflag = %2d violation = %.3g%s\n', ...
              starts(k).label, starts(k).f0, starts(k).fval, starts(k).exitflag, ...
              starts(k).constrviolation, mwecmass.internal.ternary(starts(k).kept, '  kept', ''));
    end
    fprintf('  bottom-filled x0: mass %.3f kg, displaced %.3f kg\n', p_bf.mass_total, p_bf.mass_buoyant_force);
    fprintf('  optimum: vs = %.4f m, rho = [%s]\n', x(1), sprintf(' %.1f', x(2:end)));
    fprintf('  M = %.2f kg, buoyant mass %.2f kg, ceq = %.3e, GM = %.4f m (gm_min %.3f), c = [%s]\n', ...
            props.mass_total, props.mass_buoyant_force, ceq, props.GM_L, config.gm_min, sprintf(' %.3g', c));
    fprintf('  T_heave = %.3f s, T_pitch = %.3f s, exitflag = %d, iterations = %d\n', ...
            props.periods.heave, props.periods.pitch, s2.exitflag, s2.output.iterations);
    if s2.exitflag > 0
      if c(1) > constraint_tol || abs(ceq) > constraint_tol || max(c) > constraint_tol
        error('%s: exitflag %d but max(c) = %.3g, GM floor c = %.3g, |ceq| = %.3g exceed %.0e', ...
              mode, s2.exitflag, max(c), c(1), abs(ceq), constraint_tol);
      end
    else
      fprintf('  exitflag %d: constraints not asserted (violation %.3g)\n', s2.exitflag, violation);
    end
  end
end

function s = set_all(s, value)
  names = fieldnames(s);
  for k = 1:numel(names)
    if isstruct(s.(names{k}))
      s.(names{k}) = set_all(s.(names{k}), value);
    else
      s.(names{k}) = value;
    end
  end
end
