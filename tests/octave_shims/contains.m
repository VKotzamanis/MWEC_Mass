function tf = contains(str, pattern, varargin)
%CONTAINS shim for MATLAB contains: true where str has any pattern as a substring.
%   tf = contains(str, pattern) and contains(str, pattern, 'IgnoreCase', tf).
%   str is char, cellstr or string array; pattern is char or cellstr. tf has the shape of str.
  ignore_case = false;
  if numel(varargin) >= 2 && strcmpi(varargin{1}, 'IgnoreCase')
    ignore_case = logical(varargin{2});
  end
  was_char = ischar(str);
  if was_char, str = {str}; end
  if ischar(pattern), pattern = {pattern}; end
  tf = false(size(str));
  for i = 1:numel(str)
    s = str{i};
    for j = 1:numel(pattern)
      p = pattern{j};
      if ignore_case
        found = ~isempty(strfind(lower(s), lower(p)));
      else
        found = ~isempty(strfind(s, p));
      end
      if isempty(p), found = true; end
      tf(i) = tf(i) || found;
    end
  end
end
