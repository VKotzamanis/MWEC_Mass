function r = ternary(cond, t, f)
%TERNARY Return t when cond is true, otherwise f.
% cond is a logical scalar; t and f may be any values of compatible usage.
    if cond, r = t; else, r = f; end
end
