function r = endsWith(str, pattern, varargin)
%ENDSWITH shim: see startsWith shim.
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
        t = s(end-n+1:end); if n == 0, t = ''; end
        if ic, m = strcmpi(t, pj); else m = strcmp(t, pj); end
        if n == 0, m = true; end
        r(i) = r(i) || m;
      end
    end
  end
end
