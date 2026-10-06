function bin = discretize(x, edges)
%DISCRETIZE shim for MATLAB discretize(x, edges) with a numeric edge vector.
%   bin(i) is the index k with edges(k) <= x(i) < edges(k+1); the last bin also takes
%   x == edges(end). Values outside the edges give NaN, as in MATLAB. edges must be increasing.
  edges = edges(:);
  if numel(edges) < 2 || any(diff(edges) <= 0)
    error('discretize:badEdges', 'edges must be an increasing vector with at least two entries.');
  end
  bin = NaN(size(x));
  n = numel(edges) - 1;
  for i = 1:numel(x)
    xi = x(i);
    if isnan(xi) || xi < edges(1) || xi > edges(end)
      continue;
    end
    if xi == edges(end)
      bin(i) = n;
    else
      bin(i) = find(edges <= xi, 1, 'last');
    end
  end
end
