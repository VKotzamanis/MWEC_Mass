function [periods, share, Phi, lambda] = coupled_periods_by_share(M_total, K_full)
%COUPLED_PERIODS_BY_SHARE Natural periods of the coupled 3-DOF undamped problem, each mode named
% by the DOF holding the largest share of its modal kinetic energy.
% See docs/METHODS_ENGINE.md#coupled-modal-energy-share.
% Inputs: M_total [3x3] symmetric positive-definite mass [kg, kg m, kg m^2];
%         K_full [3x3] symmetric stiffness [N/m, N/m, N m/rad], DOFs surge/heave/pitch.
% Outputs: periods [s], kinetic-energy share [%], mode shapes Phi, and lambda [rad^2/s^2],
%          all reordered to surge/heave/pitch. Shares use the full mass matrix, including
%          off-diagonal coupling, and columns sum to 100%. Non-positive lambda gives Inf.
% A positive restoring stiffness must have a mode with at least 50% share in its DOF;
% otherwise the function raises mwecmass:hydrostatics:AmbiguousModeLabel.

n_dof = 3;   % surge, heave, pitch
validateattributes(M_total, {'double'}, {'size', [3 3], 'finite'}, mfilename, 'M_total');
validateattributes(K_full,  {'double'}, {'size', [3 3], 'finite'}, mfilename, 'K_full');

[Phi, Lam] = eig(M_total \ K_full);

% M_total is symmetric positive definite and K_full symmetric, so the pencil (K_full, M_total) is
% symmetric-definite and both its eigenvalues and its eigenvectors are real. eig() is called on
% the non-symmetric product M_total\K_full, so round-off can leave imaginary parts of order eps;
% real() drops those and nothing else.
lambda = real(diag(Lam));   % [rad^2/s^2] omega_n^2 per mode
Phi    = real(Phi);

share = zeros(n_dof, n_dof);   % [%] kinetic-energy share, rows = DOF, columns = mode
for j = 1:n_dof
    phi  = Phi(:, j);
    Mphi = M_total * phi;   % [kg m; kg m; kg m^2/rad] generalised momentum per unit modal rate
    share(:, j) = 100 * (phi .* Mphi) / (phi.' * Mphi);
end

[s_max, dom] = max(share, [], 1);   % dominant DOF of each mode and its share [%]

T = inf(n_dof, 1);   % [s] natural period per mode; Inf for a mode with no restoring stiffness
positive_lambda    = lambda > 0;
T(positive_lambda) = 2*pi ./ sqrt(lambda(positive_lambda));

% First match, because only one mode can be dominated by a given DOF whenever the labelling is
% unambiguous at all, and the ambiguous case is the error below. A column dominated by DOF k with
% lambda <= 0 is an unrestored mode and yields Inf, which is the per-axis convention.
j_heave = find(dom == 2, 1);
j_pitch = find(dom == 3, 1);

periods = struct('surge', inf, 'heave', inf, 'pitch', inf);   % [s]
if ~isempty(j_heave)
    periods.heave = T(j_heave);
end
if ~isempty(j_pitch)
    periods.pitch = T(j_pitch);
end

% A DOF with a positive diagonal stiffness must own a mode; one with K(k,k) <= 0 has no restored
% mode at all (K33 = rho*g*Aw = 0 out of the water; K55 = m*g*GM = 0 when GM <= 0), so its Inf
% period above is the answer, not a failure.
dof_of_slot = [2, 3];               % heave, pitch
j_of_slot   = {j_heave, j_pitch};
for s = 1:numel(dof_of_slot)
    k   = dof_of_slot(s);
    j_k = j_of_slot{s};
    if K_full(k, k) > 0 && (isempty(j_k) || s_max(j_k) < 50)
        error('mwecmass:hydrostatics:AmbiguousModeLabel', ...
              ['DOF %d (1 surge, 2 heave, 3 pitch) has diagonal stiffness %.6g > 0, so a ' ...
               'restored mode must exist, but no coupled mode holds at least 50%% of its ' ...
               'modal kinetic energy in that DOF -- the modes cannot be named. ' ...
               'Kinetic-energy share [%%], rows surge/heave/pitch, columns in eig order: %s'], ...
              k, K_full(k, k), mat2str(share));
    end
end

% Columns reordered to (surge, heave, pitch). A slot with no named mode takes the lowest-numbered
% unclaimed column, so the reordering is a permutation in every case.
order = zeros(1, n_dof);
if ~isempty(j_heave)
    order(2) = j_heave;
end
if ~isempty(j_pitch)
    order(3) = j_pitch;
end
order(order == 0) = setdiff(1:n_dof, order(order > 0));

share  = share(:, order);
Phi    = Phi(:, order);
lambda = lambda(order);

end
