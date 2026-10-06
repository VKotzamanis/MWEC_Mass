function s = datetime(when, varargin)
%DATETIME shim for the only form the pipeline calls: char(datetime('now', 'Format', fmt)).
%   Returns the formatted time as a char row, so char() of the result is the identity.
%   fmt uses MATLAB tokens (yyyy, MM month, MMM month name, dd, HH, mm minute, ss); the shim
%   translates them to datestr tokens. Without 'Format' the MATLAB default 'dd-MMM-yyyy HH:mm:ss'
%   is used. Any other first argument raises an error.
  if nargin < 1 || ~ischar(when) || ~strcmpi(when, 'now')
    error('datetime:unsupported', 'datetime shim supports only datetime(''now'', ...).');
  end
  fmt = 'dd-MMM-yyyy HH:mm:ss';
  for k = 1:2:numel(varargin)
    if strcmpi(varargin{k}, 'Format')
      fmt = varargin{k+1};
    else
      error('datetime:unsupported', 'datetime shim supports only the Format option.');
    end
  end
  s = datestr(now, translate_format(fmt));
end

function out = translate_format(fmt)
% MATLAB datetime tokens to Octave datestr tokens (month and minute letters swap case).
  out = '';
  k = 1;
  while k <= numel(fmt)
    c = fmt(k);
    j = k;
    while j < numel(fmt) && fmt(j+1) == c
      j = j + 1;
    end
    run = j - k + 1;
    switch c
      case 'M'
        out = [out repmat('m', 1, run)]; %#ok<AGROW>
      case 'm'
        out = [out repmat('M', 1, run)]; %#ok<AGROW>
      case 's'
        out = [out repmat('S', 1, run)]; %#ok<AGROW>
      otherwise
        out = [out repmat(c, 1, run)]; %#ok<AGROW>
    end
    k = j + 1;
  end
end
