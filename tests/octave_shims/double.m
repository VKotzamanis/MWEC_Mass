function out = double(varargin)
%DOUBLE shim adding the MATLAB form double.empty(m, n) (export_results.m); double(x) is the built-in.
%   Octave 8.4 evaluates double.empty as a call double() and fails. With no arguments this shim
%   returns a struct whose field empty builds the empty array, so double.empty(0, 0) works.
  if nargin == 0
    out = struct('empty', @(varargin) zeros(varargin{:}));
  else
    out = builtin('double', varargin{:});
  end
end
