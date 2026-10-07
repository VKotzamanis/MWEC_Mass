function test_stage2_bounds()
%TEST_STAGE2_BOUNDS The upper density bound of a mode is the solid density of its material input,
%   rho_air is 1.2 in both modes, and the thin-shell minimum thickness is the one-inch input.
%   Expected values are read from the author inputs (WEC_User_Input), not from the config.
  repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  addpath(fullfile(repo_root, 'tests', 'optim'));
  in = WEC_User_Input();

  % thin shell: steel ballast, no wall module
  config = stage2_test_config('thin_shell');
  [lb, ub] = mwecmass.optim.stage2_bounds(config);
  N = config.num_ballast_sections;
  expect_ub = in.materials.thin_shell.rho_ballast;
  if ~isequal(ub(2:end), expect_ub * ones(1, N)) || expect_ub ~= 7500
    error('thin shell: ub = %s, expected %g for every module', mat2str(ub(2:end)), expect_ub);
  end
  if ~isequal(lb(2:end), in.bounds.ballast_density_bounds(1) * ones(1, N))
    error('thin shell: lb = %s', mat2str(lb(2:end)));
  end
  if config.rho_air ~= 1.2 || config.rho_air ~= in.materials.thin_shell.rho_air
    error('thin shell: config.rho_air = %g', config.rho_air);
  end
  if config.rho_ballast ~= in.materials.thin_shell.rho_ballast
    error('thin shell: config.rho_ballast = %g', config.rho_ballast);
  end
  if in.materials.thin_shell.t_min ~= 0.0254 || in.materials.thin_shell.t_init ~= 0.0254 ...
      || config.steel_t_min ~= 0.0254 || config.steel_t_init ~= 0.0254
    error('thin shell: t_min = %g, t_init = %g (config %g, %g)', in.materials.thin_shell.t_min, ...
          in.materials.thin_shell.t_init, config.steel_t_min, config.steel_t_init);
  end
  fprintf('thin_shell: ub = %g, lb = %g kg/m^3, rho_air = %g, t_min = t_init = %g m\n', ...
          ub(2), lb(2), config.rho_air, config.steel_t_min);

  % modular precast: UHPC, wall module pinned
  config = stage2_test_config('modular_precast');
  [lb, ub] = mwecmass.optim.stage2_bounds(config);
  w = config.wall_strip_index;
  rho_uhpc = in.materials.modular_precast.rho_hull;
  free = setdiff(1:N, w);
  if ~isequal(ub(1 + free), rho_uhpc * ones(1, numel(free))) || rho_uhpc ~= 2500
    error('modular precast: ub = %s, expected %g', mat2str(ub(2:end)), rho_uhpc);
  end
  if lb(1 + w) ~= rho_uhpc || ub(1 + w) ~= rho_uhpc
    error('modular precast: wall module %d not pinned to %g', w, rho_uhpc);
  end
  if any(lb(1 + free) < config.per_strip_density_lb(free)) ...
      || any(lb(1 + free) < in.bounds.ballast_density_bounds(1))
    error('modular precast: lb %s below the floors', mat2str(lb(2:end)));
  end
  if config.constructability_rho_air ~= 1.2 ...
      || config.constructability_rho_air ~= in.materials.modular_precast.rho_air
    error('modular precast: config.constructability_rho_air = %g', config.constructability_rho_air);
  end
  fprintf('modular_precast: ub = %g, wall module %d pinned, rho_air = %g kg/m^3\n', ...
          ub(1 + free(1)), w, config.constructability_rho_air);

  % preliminary has no material: the input's own upper bound stays
  config = stage2_test_config('preliminary');
  if ~isequal(config.ballast_density_bounds, in.bounds.ballast_density_bounds)
    error('preliminary: bounds = %s', mat2str(config.ballast_density_bounds));
  end
  fprintf('preliminary: bounds = %s kg/m^3 (the input)\n', mat2str(config.ballast_density_bounds));

  % Every UHPC path reads the precast air density, every thin-shell path the thin-shell one.
  expect_reads = {
    fullfile('+realise', '+modular_precast', 'solve_and_extract.m'), 'config.constructability_rho_air';
    fullfile('+realise', '+modular_precast', 'solve_and_extract.m'), 'config.constructability_rho_hull';
    fullfile('+realise', '+thin_shell', 'solve.m'), 'config.rho_air';
    fullfile('+realise', '+thin_shell', 'solve.m'), 'config.rho_ballast'};
  for k = 1:size(expect_reads, 1)
    text = fileread(fullfile(repo_root, 'src', '+mwecmass', expect_reads{k, 1}));
    if isempty(strfind(text, expect_reads{k, 2}))
      error('%s does not read %s', expect_reads{k, 1}, expect_reads{k, 2});
    end
  end
  fprintf('realiser sources read the named material inputs\n');
end
