function [V_env, V_mat, V_void, V_ballast, V_shell] = strip_partition_volumes( ...
        z_grid, A_outer, A_inner, z_ballast, strip_edges)
%STRIP_PARTITION_VOLUMES Partition each strip into envelope, material, and void.
%   Inputs are z_grid [m], section areas [m^2], ballast elevation z_ballast [m],
%   and strip_edges [m]. Outputs are [N x 1] volumes [m^3]. Below z_ballast the
%   whole section is material; above it material is the jacket annulus.
%   V_mat=V_ballast+V_shell and V_env=V_mat+V_void. A ballast breakpoint is included
%   only when it lies inside a strip, preventing out-of-range double counting.
%   See docs/METHODS_ENGINE.md#thin-shell-strip-partition
%   Outputs V_env,V_mat,V_void,V_ballast,V_shell are [N x 1] volumes [m^3] with
%   V_env=V_mat+V_void and V_mat=V_ballast+V_shell.

    z_grid  = z_grid(:);
    A_outer = A_outer(:);
    A_inner = A_inner(:);
    se      = strip_edges(:);

    assert(numel(A_outer) == numel(z_grid) && ...
           numel(A_inner) == numel(z_grid), ...
        'mwecmass:thin_shell:strip_partition_volumes:SizeMismatch', ...
        'z_grid, A_outer and A_inner must be the same length.');
    assert(numel(se) >= 2, ...
        'mwecmass:thin_shell:strip_partition_volumes:BadEdges', ...
        'strip_edges needs at least 2 entries.');
    % z_ballast = NaN silently makes BOTH the below- and above-masks
    % all-false, so both trapz branches are skipped and the method
    % returns V_mat = V_void = 0 with V_env nonzero — closure
    % violated with no diagnostic.  +/-Inf is handled sensibly
    % (all-material / all-jacket-plus-void) and stays accepted.
    assert(isscalar(z_ballast) && ~isnan(z_ballast), ...
        'mwecmass:thin_shell:strip_partition_volumes:BadZBallast', ...
        'z_ballast must be a non-NaN scalar.');
    % Descending edges send every strip down the degenerate-strip
    % `continue`, returning all zeros with no error; interleaved
    % edges silently over-count.
    assert(all(diff(se) >= 0), ...
        'mwecmass:thin_shell:strip_partition_volumes:NonMonotonicEdges', ...
        'strip_edges must be non-decreasing.');
    % interp1(..., 'linear', 0) extrapolates area to zero, so an edge
    % beyond the grid fabricates a phantom triangular volume.
    assert(se(1) >= z_grid(1) - 1e-9 && se(end) <= z_grid(end) + 1e-9, ...
        'mwecmass:thin_shell:strip_partition_volumes:EdgesOutsideGrid', ...
        ['strip_edges [%.6f, %.6f] fall outside z_grid [%.6f, %.6f]; ' ...
         'linear extrapolation to zero would fabricate volume.'], ...
        se(1), se(end), z_grid(1), z_grid(end));
    assert(all(diff(z_grid) > 0), ...
        'mwecmass:thin_shell:strip_partition_volumes:NonMonotonicGrid', ...
        'z_grid must be strictly increasing.');

    N         = numel(se) - 1;
    V_env     = zeros(N, 1);
    V_mat     = zeros(N, 1);
    V_void    = zeros(N, 1);
    V_ballast = zeros(N, 1); % below-z_ballast solid volume (V_mat's below-z_ballast term)
    V_shell   = zeros(N, 1); % above-z_ballast solid volume (V_mat's above-z_ballast term)

    for i = 1:N
        z_lo = se(i);
        z_hi = se(i + 1);
        if z_hi <= z_lo + 1e-12   % [m] degenerate-strip tolerance
            continue;
        end

        bp = [z_lo; z_hi; z_grid(z_grid > z_lo & z_grid < z_hi)];
        if z_ballast > z_lo && z_ballast < z_hi
            bp = [bp; z_ballast];   %#ok<AGROW>  N is small
        end
        bp = unique(sort(bp));
        if numel(bp) < 2
            continue;
        end

        Ao = interp1(z_grid, A_outer, bp, 'linear', 0);
        Ai = interp1(z_grid, A_inner, bp, 'linear', 0);
        % Two distinct degeneracies, deliberately handled differently.
        %
        % (a) GENUINE INVERSION -- positive outer area with the inner
        %     area at or beyond it, i.e. negative wall thickness.  The
        %     section is solid.  This matches
        %     the modular-precast strip extraction, which resolves the same
        %     degeneracy the same way.  Reporting it as void instead
        %     yields a plausible-looking but physically impossible
        %     rho_eff that PASSES the [rho_air, rho_hull] invariant
        %     test, so it must be loud.
        %
        % (b) ZERO-AREA HULL TIPS -- A_outer == A_inner == 0 at the keel
        %     and deck samples, where "Ai >= Ao" is satisfied by 0 >= 0.
        %     Benign and present on every real run (2 of 300 samples on
        %     both C1 results), so warning on it would fire every time
        %     and camouflage case (a).  Clamped silently instead.
        inverted = (Ao > 0) & (Ai >= Ao);
        if any(inverted)
            warning('mwecmass:thin_shell:strip_partition_volumes:InvertedSection', ...
                ['Negative wall thickness at %d of %d samples in strip %d ' ...
                 '(A_inner >= A_outer > 0). Treating those samples as solid.'], ...
                sum(inverted), numel(inverted), i);
        end
        Ai(inverted) = 0;
        Ai = max(0, min(Ai, Ao));   % [m^2] silent clamp: case (b)
        Aj = max(0, Ao - Ai);

        V_env(i) = trapz(bp, Ao);

        below = bp <= z_ballast;
        above = bp >= z_ballast;
        if sum(below) >= 2
            V_below_i    = trapz(bp(below), Ao(below));
            V_mat(i)     = V_mat(i) + V_below_i;
            V_ballast(i) = V_ballast(i) + V_below_i;   % same term, also kept split
        end
        if sum(above) >= 2
            V_above_i  = trapz(bp(above), Aj(above));
            V_mat(i)   = V_mat(i)   + V_above_i;
            V_shell(i) = V_shell(i) + V_above_i;   % same term, also kept split
            V_void(i)  = V_void(i)  + trapz(bp(above), Ai(above));
        end
    end
end
