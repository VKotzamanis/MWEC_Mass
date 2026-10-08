function test_density_floors()
%TEST_DENSITY_FLOORS Stage-2 density floors (contract F8): the density of a module built with a
%   t_min shell around air and no ballast, in both modes, on the cylinder fixture deck.
%   Oracle (independent of the kernel and its stand-ins): the cylinder of radius R between z0 and z1
%   read from the deck points K, P1, T; a shell of design thickness t_min is built at d = t_min +
%   0.005 t_min (eps_fit = 0.01 t_min), so the void of a hollow module is the cylinder of radius
%   R - d between max(edge_i, z0 + d) and min(edge_i+1, z1 - d).
%   Exact by construction (asserted for every kernel): V_solid + V_air = V per module, the solid
%   modules are solid (rho_min equal to rho_solid bitwise, also for a module volume V where
%   (rho_solid V) / V differs from rho_solid), the floor lies in (rho_air, rho_solid], build_config
%   stores max(lower input bound, floor) per module and the same floors in the geometry products.
%   The sum of the module volumes equal to the hull volume and the agreement with the closed form
%   are asserted at rounding when the kernel is the closed-form stand-in (geo.analytic) and
%   printed for the real kernel, whose rational faces are only approximated by quadrature.
  repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
  deck = fullfile(repo_root, 'tests', 'standins', 'fixtures', 'cylinder.ms2');
  saved_path = path();
  addpath(fullfile(repo_root, 'tests', 'standins'), '-end');   % kernel stand-ins until the real one is merged
  cleanup = onCleanup(@() path(saved_path));

  evalc('model = mwecmass.geometry.MS2Parser.parse(deck);');
  K = model.eval_point('K');
  P1 = model.eval_point('P1');
  T = model.eval_point('T');
  R = P1(1);
  z0 = K(3);
  z1 = T(3);

  modes = {'modular_precast', 'thin_shell'};
  for m = 1:numel(modes)
    mode = modes{m};
    in = deck_input(deck, mode, '');
    if strcmp(mode, 'modular_precast')
      mat = in.materials.modular_precast;
      t_min = mat.t_min;
      rho_solid = mat.rho_hull;
    else
      mat = in.materials.thin_shell;
      t_min = mat.t_min;
      rho_solid = mat.rho_shell;
    end
    rho_air = mat.rho_air;
    d = t_min + 0.005 * t_min;

    evalc('config = mwecmass.driver.build_config(in, [], struct());');
    edges = config.strip_edges(:);
    N = numel(edges) - 1;
    solid = config.wall_strip_index;
    if isempty(config.hull_solid) || isempty(config.density_floors)
      error('%s: build_config left hull_solid or density_floors empty', mode);
    end
    if strcmp(mode, 'thin_shell') && ~isempty(solid)
      error('thin_shell: wall module %d set', solid);
    end
    fl = config.density_floors;
    analytic = isfield(config.hull_solid, 'analytic') && ~isempty(config.hull_solid.analytic);

    % closed form
    V_oracle = pi * R^2 * diff(edges);
    z_void = [max(edges(1:N), z0 + d), min(edges(2:N+1), z1 - d)];
    V_air_oracle = pi * (R - d)^2 * max(0, diff(z_void, 1, 2));
    V_air_oracle(solid) = 0;
    V_solid_oracle = V_oracle - V_air_oracle;
    rho_oracle = (rho_solid * V_solid_oracle + rho_air * V_air_oracle) ./ V_oracle;

    % exact by construction
    if numel(fl.rho_min) ~= N || any(size(fl.V) ~= [N 1])
      error('%s: floors have %d entries for %d modules', mode, numel(fl.rho_min), N);
    end
    closure = abs(fl.V_solid + fl.V_air - fl.V) ./ fl.V;
    closure_bound = 8 * eps;   % two region integrals added to the module volume of the same faces
    if any(closure > closure_bound)
      error('%s: V_solid + V_air - V = %s of V exceeds %.1e', mode, mat2str(closure, 3), closure_bound);
    end
    hull_V = abs(sum(fl.V) - pi * R^2 * (z1 - z0)) / (pi * R^2 * (z1 - z0));
    if analytic && hull_V > 8 * eps
      error('%s: module volumes sum to %.3g relative off the hull volume', mode, hull_V);
    end
    if any(fl.rho_min(solid) ~= rho_solid) || any(fl.V_air(solid) ~= 0)
      error('%s: a solid module is not solid (rho_min %s, V_air %s)', mode, ...
            mat2str(fl.rho_min(solid)), mat2str(fl.V_air(solid)));
    end
    hollow = setdiff(1:N, solid);
    if any(fl.rho_min(hollow) <= rho_air) || any(fl.rho_min(hollow) > rho_solid)
      error('%s: hollow floors %s outside (%g, %g]', mode, mat2str(fl.rho_min(hollow)', 6), rho_air, rho_solid);
    end
    if any(fl.V_air(hollow) <= 0)
      error('%s: a hollow module has no void: V_air = %s', mode, mat2str(fl.V_air', 5));
    end
    recomposed = (rho_solid * fl.V_solid + rho_air * fl.V_air) ./ fl.V;
    if ~isequal(recomposed(hollow), fl.rho_min(hollow))
      error('%s: a hollow rho_min is not (rho_solid V_solid + rho_air V_air) / V', mode);
    end
    lb = mwecmass.optim.stage2_bounds(config);
    expect_lb = max(in.bounds.ballast_density_bounds(1), fl.rho_min(:))';
    if ~isequal(config.per_strip_density_lb, expect_lb)
      error('%s: per_strip_density_lb %s is not max(lower input bound, floor) %s', mode, ...
            mat2str(config.per_strip_density_lb, 6), mat2str(expect_lb, 6));
    end
    if ~isequal(lb(2:end), expect_lb) && isempty(solid)
      error('%s: Stage-2 lower bounds %s differ from the floors', mode, mat2str(lb(2:end), 6));
    end

    % closed form
    dev_V = max(abs(fl.V - V_oracle) ./ V_oracle);
    dev_air = max(abs(fl.V_air - V_air_oracle) ./ max(V_oracle, eps));
    dev_rho = max(abs(fl.rho_min - rho_oracle) ./ rho_oracle);
    if analytic
      closed_form_bound = 16 * eps;   % a handful of products and sums on each side
      if dev_V > closed_form_bound || dev_air > closed_form_bound || dev_rho > closed_form_bound
        error('%s: floors differ from the closed form by V %.2e, V_air %.2e, rho %.2e (bound %.1e)', ...
              mode, dev_V, dev_air, dev_rho, closed_form_bound);
      end
    end

    fprintf('%s: kernel %s, t_min %.4f m, d %.5f m, modules %d, solid modules %s\n', mode, ...
            ternary(analytic, 'stand-in', 'real'), t_min, d, N, mat2str(solid));
    fprintf('  floors [kg/m3] %s, closed form %s\n', mat2str(fl.rho_min', 6), mat2str(rho_oracle', 6));
    fprintf('  V closure max %.2e, hull volume relative %.2e, closed form: V %.2e, V_air %.2e, rho %.2e\n', ...
            max(closure), hull_V, dev_V, dev_air, dev_rho);
  end

  % preliminary: no material, no shell, no floors
  in = deck_input(deck, 'preliminary', '');
  evalc('config = mwecmass.driver.build_config(in, [], struct());');
  if ~isempty(config.per_strip_density_lb) || ~isempty(config.hull_solid) || ~isempty(config.density_floors)
    error('preliminary: floors or hull_solid present');
  end
  fprintf('preliminary: no floors, no hull_solid\n');

  % refusals
  evalc('config = mwecmass.driver.build_config(deck_input(deck, ''thin_shell'', ''''), [], struct());');
  args = {config.ms2_model, config.boundary_cache, config.hull_solid};
  e = config.strip_edges(:);
  ballast = struct('rho_ballast', 7500);
  expect_error(@() mwecmass.driver.density_floors(args{:}, [e(1) + 0.1; e(2:end)], 0.0254, 7500, 1.2, [], ballast), ...
               'mwecmass:driver:EdgesNotOnHull');
  expect_error(@() mwecmass.driver.density_floors(args{:}, e, 0.0254, 7500, 1.2, []), ...
               'mwecmass:driver:MissingBallastDensity');
  expect_error(@() mwecmass.driver.density_floors(args{:}, e, 0.0254, 7500, 1.2, 1:numel(e) - 1), ...
               'mwecmass:driver:NoHollowModule');
  expect_error(@() mwecmass.driver.density_floors(args{:}, e, 0.0254, 7500, 1.2, [], struct('mode', 'preliminary')), ...
               'mwecmass:driver:BadMode');
  fprintf('refusals: EdgesNotOnHull, MissingBallastDensity, NoHollowModule, BadMode\n');

  % a solid module whose volume V makes (rho_solid V) / V differ from rho_solid
  rho_solid = 2500;
  top = e(end);
  found = false;
  for k = 1:200
    edge_k = top - 0.01 * (1 + 0.1337 * k) - 0.9;
    edges_k = [e(1); linspace(e(1), edge_k, 4)(2:end)'; top];
    fl_k = mwecmass.driver.density_floors(args{:}, edges_k, 0.0254, rho_solid, 1.2, numel(edges_k) - 1, ...
                                          struct('mode', 'modular_precast'));
    V_w = fl_k.V(end);
    if (rho_solid * V_w) / V_w ~= rho_solid
      found = true;
      if fl_k.rho_min(end) ~= rho_solid || fl_k.V_solid(end) ~= V_w || fl_k.V_air(end) ~= 0
        error('solid module with V = %.17g: rho_min %.17g, V_solid %.17g, V_air %.3g', ...
              V_w, fl_k.rho_min(end), fl_k.V_solid(end), fl_k.V_air(end));
      end
      fprintf('solid module V = %.17g m^3: (rho V)/V - rho = %.3g, rho_min - rho = 0\n', ...
              V_w, (rho_solid * V_w) / V_w - rho_solid);
      break;
    end
  end
  if ~found
    error('no candidate wall module volume made (rho_solid V) / V differ from rho_solid');
  end
end

function expect_error(fun, id)
  try
    fun();
  catch err
    if ~strcmp(err.identifier, id)
      error('expected %s, got %s: %s', id, err.identifier, err.message);
    end
    return;
  end
  error('expected error %s, none raised', id);
end

function out = ternary(c, a, b)
  if c, out = a; else, out = b; end
end
