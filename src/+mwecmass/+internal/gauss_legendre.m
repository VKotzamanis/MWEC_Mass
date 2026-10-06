function [x, w] = gauss_legendre(n)
%GAUSS_LEGENDRE Return n-point Gauss-Legendre nodes and weights on [0,1].
% [x,w] = gauss_legendre(n) uses the Golub-Welsch Jacobi eigensystem;
% the rule is exact for polynomials through degree 2n-1.

    if n == 1
        x = 0.5;
        w = 1.0;
        return;
    end

    % Jacobi matrix sub-diagonal for Legendre polynomials
    k    = 1:n-1;
    beta = k ./ sqrt(4*k.^2 - 1);
    J    = diag(beta, 1) + diag(beta, -1);

    [V, D] = eig(J);
    x_ref  = diag(D);            % nodes on [−1, +1]
    w_ref  = 2 * V(1, :)'.^2;    % weights on [−1, +1]

    % Transform [−1, +1] → [0, 1]
    x = (x_ref + 1) / 2;
    w = w_ref / 2;

    % Sort ascending
    [x, idx] = sort(x);
    w = w(idx);
end
