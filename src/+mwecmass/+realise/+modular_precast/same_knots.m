function same = same_knots(a, b)
%SAME_KNOTS True when two inner-surface sets have the same pieces, degrees and knot vectors (the
%structure that F2 opts.knots_from keeps).
same = numel(a.patches) == numel(b.patches);
for p = 1:numel(a.patches)
    if ~same
        return
    end
    same = isequal(a.patches(p).surf.degree, b.patches(p).surf.degree) && ...
        isequal(a.patches(p).surf.knots, b.patches(p).surf.knots);
end
end
