function test_stage2_bottom_filled_start()
%TEST_STAGE2_BOTTOM_FILLED_START The bottom-filled start carries the displaced mass, with the
%   modules at their upper bound from the keel up, one partly filled, the rest at their lower
%   bound, and a pinned module untouched. Stand-in cylinder deck, both modes.
%   The mass identity holds to rounding: the partly filled module takes exactly the remaining
%   deficit, so the sum of N products differs from the target by a few ulp of the total.
  repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  addpath(fullfile(repo_root, 'tests', 'optim'));
  modes = {'modular_precast', 'thin_shell'};
  for m = 1:numel(modes)
    config = stage2_test_config(modes{m});
    [lb, ub] = mwecmass.optim.stage2_bounds(config);
    N = config.num_ballast_sections;
    vs = -0.5;
    x0 = mwecmass.optim.stage2_bottom_filled_start(vs, config, lb, ub);
    if x0(1) ~= vs || numel(x0) ~= N + 1 || ~isrow(x0)
      error('%s: x0 is not [vs, rho_1..rho_N]', modes{m});
    end
    rho = x0(2:end);
    if any(rho < lb(2:end)) || any(rho > ub(2:end))
      error('%s: x0 outside the bounds', modes{m});
    end

    props = mwecmass.hydrostatics.properties_3d(x0, config);
    rel = abs(props.mass_total - props.mass_buoyant_force) / props.mass_buoyant_force;
    if rel > 64 * eps
      error('%s: mass %.10g differs from displaced mass %.10g (relative %.2e)', ...
            modes{m}, props.mass_total, props.mass_buoyant_force, rel);
    end

    % state per module: 2 at the upper bound, 1 between, 0 at the lower bound; pinned modules skipped
    free = find(ub(2:end) > lb(2:end));
    state = zeros(1, numel(free));
    for j = 1:numel(free)
      i = free(j);
      if rho(i) == ub(1 + i), state(j) = 2; elseif rho(i) > lb(1 + i), state(j) = 1; end
    end
    if any(diff(state) > 0) || sum(state == 1) > 1
      error('%s: module states from the keel up %s are not filled / one partial / empty', ...
            modes{m}, mat2str(state));
    end
    pinned = find(ub(2:end) == lb(2:end));
    if any(rho(pinned) ~= lb(1 + pinned))
      error('%s: a pinned module changed', modes{m});
    end
    fprintf('%s: vs = %.2f, states from the keel %s, mass %.4f kg = displaced %.4f kg (relative %.1e)\n', ...
            modes{m}, vs, mat2str(state), props.mass_total, props.mass_buoyant_force, rel);

    % Lower bounds already heavier than the displaced mass: the start is the lower bounds.
    lb_heavy = lb;
    lb_heavy(2:end) = max(lb(2:end), min(ub(2:end), 1500));
    x_heavy = mwecmass.optim.stage2_bottom_filled_start(vs, config, lb_heavy, ub);
    p_heavy = mwecmass.hydrostatics.properties_3d([vs, lb_heavy(2:end)], config);
    if p_heavy.mass_total <= p_heavy.mass_buoyant_force
      error('%s: raised lower bounds are not heavier than the displaced mass', modes{m});
    end
    if ~isequal(x_heavy, [vs, lb_heavy(2:end)])
      error('%s: lower bounds heavier than the displaced mass must return the lower bounds', modes{m});
    end
    fprintf('%s: lower bounds %.1f kg > displaced %.1f kg return the lower bounds\n', ...
            modes{m}, p_heavy.mass_total, p_heavy.mass_buoyant_force);
  end
end
