function out = logical(varargin)
%LOGICAL shim adding the MATLAB form logical.empty(m, n) (export_results.m); logical(x) is the built-in.
%   See double.m.
  if nargin == 0
    out = struct('empty', @(varargin) false(varargin{:}));
  else
    out = builtin('logical', varargin{:});
  end
end
