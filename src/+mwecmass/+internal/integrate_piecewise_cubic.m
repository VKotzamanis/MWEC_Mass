function I = integrate_piecewise_cubic(pp, a, b)
%INTEGRATE_PIECEWISE_CUBIC Integrate a piecewise cubic on [a,b].
% pp.breaks and pp.coefs follow MATLAB's piecewise-polynomial format;
% each overlapping interval is integrated analytically.

    I = 0;
    breaks = pp.breaks(:)';
    coefs  = pp.coefs;
    for k = 1:length(breaks)-1
        t_lo = max(a, breaks(k))   - breaks(k);
        t_hi = min(b, breaks(k+1)) - breaks(k);
        if t_hi <= t_lo + 1e-14, continue; end
        c = coefs(k, :);
        I = I + c(1)*(t_hi^4 - t_lo^4)/4 ...
              + c(2)*(t_hi^3 - t_lo^3)/3 ...
              + c(3)*(t_hi^2 - t_lo^2)/2 ...
              + c(4)*(t_hi   - t_lo);
    end
end
