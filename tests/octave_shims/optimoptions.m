function opts = optimoptions(solver, varargin)
%OPTIMOPTIONS shim: returns a plain struct of the given name/value pairs.
%   opts = optimoptions('fmincon', 'Name', value, ...) and optimoptions(opts, 'Name', value, ...).
%   The solver name is stored in opts.SolverName. Names are kept as written; the fmincon shim
%   reads them case-insensitively. No defaults are filled in.
  if isstruct(solver)
    opts = solver;
  else
    opts = struct('SolverName', solver);
  end
  if mod(numel(varargin), 2) ~= 0
    error('optimoptions:badPairs', 'Options must be name/value pairs.');
  end
  for k = 1:2:numel(varargin)
    opts.(varargin{k}) = varargin{k+1};
  end
end
