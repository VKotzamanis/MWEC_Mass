function test_stage2_constraints()
%TEST_STAGE2_CONSTRAINTS Stage-2 constraints are the GM floor and the adjacent density ratio, and
%   nothing else: no monotonic-density and no minimum-mass constraint, in either copy.
%   Stand-in cylinder deck, both realisation modes. Equality of c and ceq with the formulas written
%   out here is by identity of two code paths (same operations, same order), so isequal.
  repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  addpath(fullfile(repo_root, 'tests', 'optim'));
  modes = {'modular_precast', 'thin_shell'};
  for m = 1:numel(modes)
    config = stage2_test_config(modes{m});
    N = config.num_ballast_sections;
    % Density rising with height: the removed monotonic constraint would have been violated here.
    x = [0.3, linspace(300, 900, N)];
    densities = x(2:end);
    w = config.wall_strip_index;
    pairs = [];
    for i = 1:N-1
      if isempty(w) || (i ~= w && i + 1 ~= w), pairs(end+1) = i; end %#ok<AGROW>
    end
    props = mwecmass.hydrostatics.properties_3d(x, config);
    expected = 1 - props.GM_L / config.gm_min;
    for i = pairs
      expected(end+1, 1) = densities(i) / (densities(i+1) + 1) - config.max_density_ratio; %#ok<AGROW>
    end
    [c, ceq] = mwecmass.optim.stage2_constraints(x, config);
    if ~isequal(c, expected)
      error('%s: c differs from GM floor + density-ratio pairs', modes{m});
    end
    if ~isequal(ceq, props.mass_total / props.mass_buoyant_force - 1)
      error('%s: ceq is not mass / buoyant force - 1', modes{m});
    end

    % The fallback for a failed evaluation has the size of the normal output.
    broken = rmfield(config, 'max_density_ratio');
    [c_fb, ceq_fb] = mwecmass.optim.stage2_constraints(x, broken);
    if numel(c_fb) ~= numel(c) || any(c_fb ~= 1) || ceq_fb ~= 1
      error('%s: fallback output has %d entries, normal %d', modes{m}, numel(c_fb), numel(c));
    end
    fprintf('%s: %d inequality entries (1 GM + %d ratio pairs), ceq = %.6g, fallback size equal\n', ...
            modes{m}, numel(c), numel(pairs), ceq);
  end

  % check_3d_convergence no longer judges density order.
  config = stage2_test_config('thin_shell');
  x = [0.3, linspace(300, 900, config.num_ballast_sections)];
  out = struct('constrviolation', 0, 'firstorderopt', 0, 'iterations', 1);
  [converged, qm] = mwecmass.optim.check_3d_convergence(x, config, 1, out);
  if isfield(qm, 'monotonic')
    error('quality_metrics still has a monotonic field');
  end
  if converged ~= (qm.mass_balance && qm.GM_satisfied && qm.fmincon_optimal)
    error('converged is not the conjunction of mass balance, GM and solver status');
  end
  fprintf('check_3d_convergence: no monotonic field, converged = %d from mass balance %d, GM %d, solver %d\n', ...
          converged, qm.mass_balance, qm.GM_satisfied, qm.fmincon_optimal);

  % The deleted names are gone from the sources.
  names = 'c_mono|c_monotonic|c_mass_min|config\.m_min_constructability|geo\.m_min_constructability|Monotonic density';
  hits = scan_sources(fullfile(repo_root, 'src'), names);
  if ~isempty(hits)
    error('deleted names still referenced: %s', strjoin(hits, ', '));
  end
  fprintf('no reference to the deleted constraints in src/\n');
end

function hits = scan_sources(folder, pattern)
  hits = {};
  entries = dir(folder);
  for k = 1:numel(entries)
    name = entries(k).name;
    if any(strcmp(name, {'.', '..'})), continue; end
    full = fullfile(folder, name);
    if entries(k).isdir
      hits = [hits, scan_sources(full, pattern)]; %#ok<AGROW>
    elseif ~isempty(regexp(name, '\.m$', 'once'))
      if ~isempty(regexp(fileread(full), pattern, 'once')), hits{end+1} = full; end %#ok<AGROW>
    end
  end
end
