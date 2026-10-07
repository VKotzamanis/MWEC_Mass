function [x, fval, exitflag, output, lambda] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options)
%FMINCON test shim over Octave's core sqp (GNU Octave 8.4 has no fmincon).
%   [x, fval, exitflag, output, lambda] = fmincon(fun, x0, A, b, Aeq, beq, lb, ub, nonlcon, options)
%   Trailing arguments may be omitted; [] means "none". x0 may be row or column; fun, nonlcon and
%   the OutputFcn see x in the shape of x0, and x is returned in that shape.
%
%   Problem mapping to sqp(x0, f, g, h, lb, ub, maxiter, tol):
%     equalities   g(x) = [ceq(x); Aeq*x - beq] = 0
%     inequalities h(x) = -[c(x); A*x - b] >= 0
%     bounds       lb, ub passed to sqp (infinite entries dropped)
%   x0 is first projected onto [lb, ub], as MATLAB's sqp does; it may violate every other
%   constraint (sqp starts from an infeasible point).
%
%   Options read (case-insensitive, all optional): MaxIterations (default 400), OptimalityTolerance
%   and ConstraintTolerance (default 1e-6; sqp has a single tolerance, set to the smaller of the
%   two), OutputFcn (function handle or cell of handles). Other options (Algorithm, Display,
%   StepTolerance, MaxFunctionEvaluations, FiniteDifferenceStepSize, ScaleProblem, ...) are accepted
%   and ignored: sqp has no function-evaluation cap, uses its own forward-difference step
%   sqrt(eps), and prints nothing. Gradients are always finite differences.
%
%   exitflag from sqp's info:
%     101 (terminated normally)                  ->  1
%     104 (step too small) or 102 (BFGS stalled) ->  2 if constrviolation <= ConstraintTolerance,
%                                                     otherwise -2 (no feasible point found)
%     103 (maximum iterations reached)           ->  0
%     OutputFcn requested stop at 'init'         -> -1
%
%   output: iterations (sqp iterations that took a step), funcCount (sqp objective evaluations;
%   finite-difference evaluations not counted), constrviolation (max violation at x, recomputed
%   here), stepsize (NaN), firstorderopt (NaN: sqp does not report it), algorithm, message, info
%   (the raw sqp code).
%   lambda: eqlin, eqnonlin, ineqlin, ineqnonlin, lower, upper, in MATLAB's sign convention.
%   When sqp returns no multipliers (its last QP subproblem failed), lambda holds empty and zero fields.
%
%   OutputFcn is called with state 'init' (iteration 0, at x0) and 'done' (at the solution).
%   Per-iteration ('iter') calls are not available because sqp has no callback.
  if nargin < 2, error('fmincon:nargin', 'fmincon needs at least fun and x0.'); end
  if nargin < 3, A = []; end
  if nargin < 4, b = []; end
  if nargin < 5, Aeq = []; end
  if nargin < 6, beq = []; end
  if nargin < 7, lb = []; end
  if nargin < 8, ub = []; end
  if nargin < 9, nonlcon = []; end
  if nargin < 10 || isempty(options), options = struct(); end

  x0_shape = size(x0);
  n = numel(x0);
  xc = x0(:);
  lbc = full_bound(lb, n, -Inf);
  ubc = full_bound(ub, n, Inf);
  if any(lbc > ubc)
    error('fmincon:badBounds', 'Some lower bound exceeds its upper bound.');
  end
  xc = min(max(xc, lbc), ubc);

  maxit = get_opt(options, 'MaxIterations', 400);
  tol = min(get_opt(options, 'OptimalityTolerance', 1e-6), ...
            get_opt(options, 'ConstraintTolerance', 1e-6));
  ctol = get_opt(options, 'ConstraintTolerance', 1e-6);
  outfcns = get_opt(options, 'OutputFcn', {});
  if isa(outfcns, 'function_handle'), outfcns = {outfcns}; end

  cache = containers.Map();
  shp = @(v) reshape(v, x0_shape);
  f_sq = @(v) cached_objective(cache, fun, shp, v);
  nl = @(v) cached_nonlcon(cache, nonlcon, shp, v);

  has_eq = ~isempty(Aeq) || ~isempty(nonlcon);
  has_in = ~isempty(A) || ~isempty(nonlcon);
  g_fn = [];
  h_fn = [];
  if has_eq
    g_fn = @(v) eq_values(nl, Aeq, beq, v);
  end
  if has_in
    h_fn = @(v) in_values(nl, A, b, v);
  end

  fval0 = f_sq(xc);
  ov = struct('iteration', 0, 'funccount', 0, 'fval', fval0, 'constrviolation', ...
              violation(nl, A, b, Aeq, beq, lbc, ubc, xc), 'stepsize', NaN, ...
              'firstorderopt', NaN, 'lssteplength', NaN, 'directionalderivative', NaN, ...
              'procedure', '');
  if call_output(outfcns, shp(xc), ov, 'init')
    x = shp(xc); fval = fval0; exitflag = -1;
    output = make_output(0, 0, ov.constrviolation, -1, 'Stopped by OutputFcn.');
    lambda = empty_lambda(n);
    return;
  end

  finite_lb = any(isfinite(lbc));
  finite_ub = any(isfinite(ubc));
  if finite_lb || finite_ub
    lb_arg = lbc; ub_arg = ubc;
  else
    lb_arg = []; ub_arg = [];
  end
  [xs, fval, info, iter, nf, lam] = sqp(xc, f_sq, g_fn, h_fn, lb_arg, ub_arg, maxit + 1, tol);
  iterations = iter - 1;

  cv = violation(nl, A, b, Aeq, beq, lbc, ubc, xs);
  switch info
    case 101
      exitflag = 1; msg = 'Local minimum found (sqp terminated normally).';
    case {102, 104}
      if cv <= ctol
        exitflag = 2; msg = 'Step size below tolerance; constraints satisfied.';
      else
        exitflag = -2; msg = 'No feasible point found; step size below tolerance.';
      end
    case 103
      exitflag = 0; msg = 'Maximum number of iterations reached.';
    otherwise
      error('fmincon:unknownInfo', 'Unexpected sqp info %d.', info);
  end
  x = shp(xs);
  output = make_output(iterations, nf, cv, info, msg);
  lambda = split_lambda(lam, nl, xs, Aeq, A, lbc, ubc);

  ov.iteration = iterations; ov.funccount = nf; ov.fval = fval; ov.constrviolation = cv;
  call_output(outfcns, x, ov, 'done');
end

function v = full_bound(bd, n, fill)
  if isempty(bd)
    v = fill * ones(n, 1);
  else
    v = bd(:);
    if numel(v) ~= n
      error('fmincon:badBounds', 'Bound vectors must have one entry per variable.');
    end
  end
end

function v = get_opt(options, name, default)
  v = default;
  names = fieldnames(options);
  k = find(strcmpi(names, name), 1);
  if ~isempty(k) && ~isempty(options.(names{k}))
    v = options.(names{k});
  end
end

function f = cached_objective(cache, fun, shp, v)
  if isKey(cache, 'f')
    s = cache('f');
    if isequal(s.x, v)
      f = s.val;
      return;
    end
  end
  f = fun(shp(v));
  if ~isscalar(f) || ~isreal(f)
    error('fmincon:badObjective', 'The objective must return a real scalar.');
  end
  cache('f') = struct('x', v, 'val', f);
end

function [c, ceq] = cached_nonlcon(cache, nonlcon, shp, v)
  if isempty(nonlcon)
    c = zeros(0, 1); ceq = zeros(0, 1);
    return;
  end
  if isKey(cache, 'n')
    s = cache('n');
    if isequal(s.x, v)
      c = s.c; ceq = s.ceq;
      return;
    end
  end
  [c, ceq] = nonlcon(shp(v));
  c = c(:); ceq = ceq(:);
  if isempty(c), c = zeros(0, 1); end
  if isempty(ceq), ceq = zeros(0, 1); end
  cache('n') = struct('x', v, 'c', c, 'ceq', ceq);
end

function g = eq_values(nl, Aeq, beq, v)
  [~, ceq] = nl(v);
  g = ceq;
  if ~isempty(Aeq), g = [g; Aeq * v - beq(:)]; end
end

function h = in_values(nl, A, b, v)
  c = nl(v);
  h = -c;
  if ~isempty(A), h = [h; -(A * v - b(:))]; end
end

function cv = violation(nl, A, b, Aeq, beq, lbc, ubc, v)
  [c, ceq] = nl(v);
  parts = [abs(ceq); max(c, 0); max(lbc - v, 0); max(v - ubc, 0)];
  if ~isempty(A), parts = [parts; max(A * v - b(:), 0)]; end
  if ~isempty(Aeq), parts = [parts; abs(Aeq * v - beq(:))]; end
  if isempty(parts), cv = 0; else cv = max(parts); end
end

function out = make_output(iterations, nf, cv, info, msg)
  out = struct('iterations', iterations, 'funcCount', nf, 'constrviolation', cv, ...
               'stepsize', NaN, 'firstorderopt', NaN, 'algorithm', 'sqp (Octave core sqp)', ...
               'message', msg, 'info', info);
end

function lam = empty_lambda(n)
  lam = struct('eqlin', [], 'eqnonlin', [], 'ineqlin', [], 'ineqnonlin', [], ...
               'lower', zeros(n, 1), 'upper', zeros(n, 1));
end

function lam = split_lambda(l, nl, xs, Aeq, A, lbc, ubc)
% sqp stacks multipliers as [g; h; lower bounds; upper bounds]; g = [ceq; Aeq rows],
% h = -[c; A rows]. sqp's equality sign is opposite to MATLAB's; inequalities agree.
  n = numel(xs);
  [c, ceq] = nl(xs);
  lo_idx = find(isfinite(lbc));
  up_idx = find(isfinite(ubc));
  counts = [numel(ceq), rows(Aeq), numel(c), rows(A), numel(lo_idx), numel(up_idx)];
  l = l(:);
  if isempty(l)
    % sqp returns no multipliers when its last QP subproblem failed (qp gives a 0x0 lambda).
    lam = empty_lambda(n);
    return;
  end
  if numel(l) ~= sum(counts)
    error('fmincon:lambdaSize', 'sqp returned %d multipliers for %d constraints.', numel(l), sum(counts));
  end
  stop_at = cumsum(counts);
  seg = @(i) l(stop_at(i) - counts(i) + (1:counts(i)));
  lam = empty_lambda(n);
  lam.eqnonlin = -seg(1);
  lam.eqlin = -seg(2);
  lam.ineqnonlin = seg(3);
  lam.ineqlin = seg(4);
  lam.lower(lo_idx) = seg(5);
  lam.upper(up_idx) = seg(6);
end

function stop = call_output(outfcns, x, ov, state)
  stop = false;
  for k = 1:numel(outfcns)
    stop = stop || logical(outfcns{k}(x, ov, state));
  end
end
