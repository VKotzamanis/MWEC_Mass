function r = startsWith(str, pattern, varargin)
%STARTSWITH shim: Octave 8's startsWith breaks on whitespace-only patterns (cellstr strips them).
  if ischar(str), str = {str}; end
  if ischar(pattern), pattern = {pattern}; end
  ic = false;
  if numel(varargin) >= 2, ic = logical(varargin{2}); end
  r = false(size(str));
  for i = 1:numel(str)
    s = str{i}; if isempty(s), s = ''; end
    for j = 1:numel(pattern)
      pj = pattern{j}; n = numel(pj);
      if numel(s) >= n
        if ic, m = strcmpi(s(1:n), pj); else m = strcmp(s(1:n), pj); end
        if n == 0, m = true; end
        r(i) = r(i) || m;
      end
    end
  end
end
