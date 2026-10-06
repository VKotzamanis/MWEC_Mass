function tf = issorted(x, varargin)
%ISSORTED shim: adds the MATLAB modes 'strictascend' and 'strictdescend' that Octave 8.4 lacks.
%   Every other call form is passed to Octave's built-in issorted unchanged.
  if numel(varargin) == 1 && ischar(varargin{1}) && any(strcmpi(varargin{1}, {'strictascend', 'strictdescend'}))
    if ~isvector(x) && ~isempty(x)
      error('issorted:notVector', 'The strict modes are supported for vectors only.');
    end
    d = diff(x(:));
    if strcmpi(varargin{1}, 'strictascend')
      tf = all(d > 0);
    else
      tf = all(d < 0);
    end
  else
    tf = builtin('issorted', x, varargin{:});
  end
end
