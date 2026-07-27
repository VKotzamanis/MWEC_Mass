classdef MWEC_Tuning_Kernels
%MWEC_TUNING_KERNELS  Computational static methods for the MWEC Stage-1
%                     natural-period placement pipeline (v2 — coupled engine).
%
%   v2 methodology (see MWEC_TUNING_FRAMEWORK_PLAN.md):
%     - One engine, two closures (Evans / 40%-critical), selected by cfg.
%     - Dynamics: full 3x3 COUPLED frequency-domain solve  x = Z^-1 Fe, so the
%       surge-pitch hydrodynamic coupling (A15,B15) — and the destructive
%       interference it carries — is in the absorbed-power number, not assumed.
%     - Objective: dimensionless eta_Falnes in [0,1] = captured spectral energy
%       divided by the Falnes-limited capturable energy.  The surge+pitch dipole
%       is counted ONCE (nu=2 combined), heave monopole nu=1, ceiling 3 g/omega^2.
%       Same argmax as raw <P_abs> (denominator is placement-independent).
%     - Decision variables: the three natural frequencies (= design natural
%       periods, the ONLY deliverable).  K_total is a silent internal
%       intermediate; PTO stiffness is never an output.  Idealised diagonal
%       per-mode PTO (only hull coupling A15,B15 acts in-model).
%     - Surge-pitch separation: emergent from the coupled objective, with the
%       band partition (surge low-omega / pitch high-omega) as a structural
%       backstop.
%
%   Method inventory by step:
%     T1  build_S_ew, compute_iec_moments
%     T2  energy_band, band_partition, bem_at_design (full-matrix interpolants)
%     T3  compute_closure         per-mode B_PTO (Evans / critical_damping)
%     T4  coupled_rao             3x3 Z^-1 Fe steady-state solve
%         evaluate_placement      eta_Falnes + <P_abs> + capture width at a placement
%         place_natural_periods   GA -> SQP max eta_Falnes (no Picard, no FWHM)
%     T5  absorption_vs_T         per-mode CW(T) vs Falnes ceiling (combined dipole)
%     Sens sensitivity_A_inf      re-place with A(omega) -> A_inf
%     Util resolve_L_ref, clamp_x, interp_matrix3, interp_fe3

    methods (Static)

        %% =================================================================
        %%  T1 — Climate spectrum
        %% =================================================================

        function [F_ew, S_ew, meta] = build_S_ew(climateGrid, omega_out, rho, g_acc, m0_pct_gate)
        %BUILD_S_EW  Probability-weighted empirical climate spectrum.
        %   Empirical-only contract (schema v2.1.0): a cell with p>0 but no WAM
        %   spectrum hard-errors.  F_ew is deep-water energy-flux density
        %   J(omega) = rho*g^2/(2*omega) * S_ew(omega)  [W/m per rad/s].
        %   Gate G1: |m0(S_ew) - sum(p*Hs^2)/16| / target <= m0_pct_gate.
            required = {'probability_grid','source_grid','S_omega_grid', ...
                        'omega','Hs_centers','Te_centers'};
            for k = 1:numel(required)
                if ~isfield(climateGrid, required{k})
                    error('MWEC_Tuning_Kernels:build_S_ew:missingField', ...
                          'climateGrid missing required field ''%s''.', required{k});
                end
            end

            prob     = climateGrid.probability_grid;
            src      = climateGrid.source_grid;
            Sgrid    = climateGrid.S_omega_grid;
            omega_cg = climateGrid.omega(:).';
            Hs_c     = climateGrid.Hs_centers(:);

            if nargin < 2 || isempty(omega_out)
                omega_out = omega_cg(:);
            end
            omega_out = omega_out(:);
            Nf = numel(omega_out);

            [N_H, N_T] = size(prob);
            S_ew = zeros(Nf, 1);
            p_wam = 0;  n_wam = 0;

            for i = 1:N_H
                for j = 1:N_T
                    p_ij = prob(i,j);
                    if p_ij <= 0, continue; end
                    src_ij = src(i,j);
                    if src_ij == "wam2d"
                        S_cg = squeeze(Sgrid(i, j, :));
                        S_cg = S_cg(:);
                        if numel(omega_cg) > 1
                            S_cell = interp1(omega_cg, S_cg, omega_out, 'pchip', 0);
                        else
                            S_cell = zeros(Nf, 1);
                        end
                        S_cell = max(S_cell, 0);
                        p_wam = p_wam + p_ij;
                        n_wam = n_wam + 1;
                        S_ew = S_ew + p_ij * S_cell;
                    else
                        error('MWEC_Tuning_Kernels:build_S_ew:noSpectrum', ...
                              ['Cell (i=%d, j=%d) has p=%.4f but source=%s. ' ...
                               'Schema v2.1.0 is empirical-only.'], ...
                              i, j, p_ij, char(src_ij));
                    end
                end
            end

            omega_safe = max(omega_out, 1e-10);
            F_ew = (rho * g_acc^2 / 2) .* S_ew ./ omega_safe;   % [W/m per rad/s]

            % G1 — m0 budget gate
            m0_Sew    = trapz(omega_out, S_ew);
            m0_target = sum(sum(prob .* repmat(Hs_c.^2, 1, N_T))) / 16;
            if m0_target > 0
                m0_err_pct = 100 * abs(m0_Sew - m0_target) / m0_target;
            else
                m0_err_pct = NaN;
            end
            if isfinite(m0_err_pct) && m0_err_pct > m0_pct_gate
                error('MWEC_Tuning_Kernels:build_S_ew:m0Gate', ...
                      ['G1 FAIL: m0(S_ew)=%.5f m^2 vs target=%.5f m^2 ' ...
                       '(err %.2f%% > %.1f%% gate).  Check omega_out range ' ...
                       '[%.3f, %.3f] rad/s or climate-grid construction.'], ...
                      m0_Sew, m0_target, m0_err_pct, m0_pct_gate, ...
                      omega_out(1), omega_out(end));
            end

            meta = struct( ...
                'n_cells_wam',    n_wam, ...
                'p_wam_fraction', p_wam, ...
                'm0_Sew',         m0_Sew, ...
                'm0_target',      m0_target, ...
                'm0_err_pct',     m0_err_pct, ...
                'omega_out',      omega_out, ...
                'gate_G1_pass',   true);
        end


        function IEC = compute_iec_moments(omega, S_ew)
        %COMPUTE_IEC_MOMENTS  Omega-form IEC TS 62600-101 spectral moments.
            omega = omega(:);  S_ew = S_ew(:);
            m0  = trapz(omega, S_ew);
            m_1 = trapz(omega, S_ew ./ max(omega, 1e-12));
            m1  = trapz(omega, S_ew .* omega);
            m2  = trapz(omega, S_ew .* omega.^2);
            [~, ip] = max(S_ew);
            omega_p = omega(ip);
            IEC = struct( ...
                'm0',      m0, ...
                'm_1',     m_1, ...
                'm1',      m1, ...
                'm2',      m2, ...
                'Hm0',     4*sqrt(max(m0, 0)), ...
                'Te',      2*pi * m_1 / max(m0, eps), ...
                'T01',     2*pi * m0  / max(m1, eps), ...
                'T02',     2*pi * sqrt(max(m0, 0)/max(m2, eps)), ...
                'Tp',      2*pi / max(omega_p, eps), ...
                'omega_p', omega_p, ...
                'eps_bw',  sqrt(max(0, 1 - m1.^2/(max(m0, eps)*max(m2, eps)))));
        end


        %% =================================================================
        %%  T2 — Band, partition, BEM at realised CG
        %% =================================================================

        function band = energy_band(omega, S_ew, pct_low, pct_high, T_min_s)
        %ENERGY_BAND  CDF-based central energy band with T_min device-physics floor.
            if pct_low <= 0 || pct_high >= 1 || pct_low >= pct_high
                error('MWEC_Tuning_Kernels:energy_band:badPct', ...
                      '0 < pct_low < pct_high < 1 required; got %.3f, %.3f.', pct_low, pct_high);
            end
            omega = omega(:);  S_ew = S_ew(:);
            if any(omega <= 0)
                error('MWEC_Tuning_Kernels:energy_band:negFreq', 'omega must be strictly positive.');
            end
            m0 = trapz(omega, S_ew);
            if m0 <= 0
                error('MWEC_Tuning_Kernels:energy_band:zeroEnergy', 'S_ew integrates to zero.');
            end
            c = cumtrapz(omega, S_ew) / m0;
            [c_u, ia] = unique(c, 'last');
            omega_u   = omega(ia);
            omega_L = interp1(c_u, omega_u, pct_low,  'linear', 'extrap');
            omega_H = interp1(c_u, omega_u, pct_high, 'linear', 'extrap');
            omega_L = max(omega_L, omega(1));
            omega_H = min(omega_H, omega(end));

            omega_H_raw = omega_H;
            if T_min_s > 0
                omega_max_phys = 2*pi / T_min_s;
                if omega_H > omega_max_phys
                    omega_H = omega_max_phys;
                end
            end

            band = struct( ...
                'omega_L',     omega_L, ...
                'omega_H',     omega_H, ...
                'omega_H_raw', omega_H_raw, ...
                'B90',         omega_H - omega_L, ...
                'T_H',         2*pi / omega_L, ...   % longest period (low omega)
                'T_L',         2*pi / omega_H, ...   % shortest period (high omega)
                'pct_low',     pct_low, ...
                'pct_high',    pct_high, ...
                'T_min_s',     T_min_s);
        end


        function info = band_partition(omega, S_ew, omega_L, omega_H)
        %BAND_PARTITION  Energy-weighted centroid + optional bimodal valley.
        %   Partition serves as the surge/pitch separation backstop: surge is
        %   bounded to [omega_L, partition], pitch to [partition, omega_H].
            omega = omega(:);  S_ew = S_ew(:);
            mask  = (omega >= omega_L) & (omega <= omega_H);
            om_b  = omega(mask);  Sb = S_ew(mask);
            m0_b  = trapz(om_b, Sb);
            if m0_b <= 0
                omega_c = (omega_L + omega_H) / 2;
            else
                omega_c = trapz(om_b, Sb .* om_b) / m0_b;
            end

            N_b   = numel(om_b);
            sigma = max(0.05 * N_b, 2);
            hw    = min(ceil(2.5 * sigma), floor((N_b - 1) / 2));
            kk    = (-hw:hw)';
            kern  = exp(-0.5 * (kk / sigma).^2);
            kern  = kern / sum(kern);
            Sb_sm = conv(Sb, kern, 'same');

            min_height = 0.10 * max(Sb_sm);
            is_max = false(N_b, 1);
            for i = 2:N_b-1
                if Sb_sm(i) > Sb_sm(i-1) && Sb_sm(i) > Sb_sm(i+1) && Sb_sm(i) >= min_height
                    is_max(i) = true;
                end
            end
            pk_idx  = find(is_max);
            pk_vals = Sb_sm(pk_idx);

            bimodal      = false;
            valley_ratio = NaN;
            valley_omega = NaN;
            peak_omegas  = [];

            if numel(pk_idx) >= 2
                [~, order] = sort(pk_vals, 'descend');
                idx_a = pk_idx(order(1));
                idx_b = pk_idx(order(2));
                if idx_a > idx_b, [idx_a, idx_b] = deal(idx_b, idx_a); end
                [val_val, val_rel] = min(Sb_sm(idx_a:idx_b));
                val_abs = idx_a + val_rel - 1;
                ratio   = val_val / max(pk_vals(order(1)), pk_vals(order(2)));
                if ratio < 0.60
                    bimodal      = true;
                    valley_ratio = ratio;
                    valley_omega = om_b(val_abs);
                    peak_omegas  = [om_b(idx_a); om_b(idx_b)];
                end
            end

            if bimodal,  partition_omega = valley_omega;
            else,        partition_omega = omega_c;        end

            info = struct( ...
                'omega_c',         omega_c, ...
                'partition_omega', partition_omega, ...
                'bimodal',         bimodal, ...
                'valley_ratio',    valley_ratio, ...
                'valley_omega',    valley_omega, ...
                'peak_omegas',     peak_omegas);
        end


        function bem = bem_at_design(final_props, hydro_table, gates, verbose)
        %BEM_AT_DESIGN  WAMIT cache -> realised draft + full 3x3 CG-frame matrices.
        %
        %   Canonical CG congruence (do NOT deviate):
        %       T = [1 0 -z_cg; 0 1 0; 0 0 1]
        %       A_CG = T'*A_O*T,  B_CG = T'*B_O*T,  Fe5_CG = Fe5 - z_cg*Fe1
        %   z_cg = final_props.CG_total(3) (world frame).
        %
        %   Stores the FULL 3x3 A_3DOF, B_3DOF (off-diagonals A15,B15 retained —
        %   the coupled engine consumes them) plus complex Fe_3DOF.
        %
        %   Gates: G2a draft in BEM range; G2b fixed-point T_n converged;
        %          G2c T_n^pitch matches final_props.periods.pitch to tol.
            req_h = {'drafts', 'omega', 'A', 'B'};
            for k = 1:numel(req_h)
                if ~isfield(hydro_table, req_h{k})
                    error('MWEC_Tuning_Kernels:bem_at_design:badHydroTable', ...
                          'hydro_table missing required field ''%s''.', req_h{k});
                end
            end
            omega_BEM = hydro_table.omega(:);
            N_omega   = numel(omega_BEM);
            drafts    = hydro_table.drafts(:);
            n_drafts  = numel(drafts);

            req_fp = {'mass_total', 'Iyy', 'draft', 'CG_total', 'K_hydro'};
            for k = 1:numel(req_fp)
                if ~isfield(final_props, req_fp{k})
                    error('MWEC_Tuning_Kernels:bem_at_design:badFinalProps', ...
                          'final_props missing required field ''%s''.', req_fp{k});
                end
            end
            M_body   = final_props.mass_total;
            Iyy      = final_props.Iyy;
            draft    = final_props.draft;
            z_cg     = final_props.CG_total(3);
            % K_hydro convention: 3x3 matrix in [surge,heave,pitch] order, so
            % (2,2)=heave, (3,3)=pitch (matches WEC_GM fp.K_hydro(2,2)).  Assert
            % it explicitly — a 6x6 or vector store would silently mis-map modes.
            if ~isequal(size(final_props.K_hydro), [3, 3])
                error('MWEC_Tuning_Kernels:bem_at_design:badKhydro', ...
                      ['final_props.K_hydro must be the 3x3 [surge,heave,pitch] ' ...
                       'reduced hydrostatic-stiffness matrix; got size %s.'], ...
                      mat2str(size(final_props.K_hydro)));
            end
            K_hydro  = diag(final_props.K_hydro);
            K_hydro  = K_hydro(:);

            % G2a — draft within BEM range
            [drafts_sorted, idx_sort] = sort(drafts);
            if draft < drafts_sorted(1) - 1e-6 || draft > drafts_sorted(end) + 1e-6
                error('MWEC_Tuning_Kernels:bem_at_design:G2a', ...
                      ['G2a FAIL: realised draft %.4f m outside BEM range ' ...
                       '[%.4f, %.4f] m.'], draft, drafts_sorted(1), drafts_sorted(end));
            end

            A_cell_sorted = hydro_table.A(idx_sort);
            B_cell_sorted = hydro_table.B(idx_sort);

            if n_drafts == 1
                j_lo = 1; j_hi = 1; w = 0;
                A_6_at_draft = A_cell_sorted{1};
                B_6_at_draft = B_cell_sorted{1};
            else
                j_hi = find(drafts_sorted >= draft, 1, 'first');
                if isempty(j_hi) || j_hi == 1
                    j_lo = 1; j_hi = min(2, n_drafts);
                else
                    j_lo = j_hi - 1;
                end
                if j_lo == j_hi
                    w = 1.0;
                else
                    w = (draft - drafts_sorted(j_lo)) / (drafts_sorted(j_hi) - drafts_sorted(j_lo));
                end
                A_6_at_draft = (1 - w) * A_cell_sorted{j_lo} + w * A_cell_sorted{j_hi};
                B_6_at_draft = (1 - w) * B_cell_sorted{j_lo} + w * B_cell_sorted{j_hi};
            end

            idx_3dof = [1, 3, 5];
            T = [1, 0, -z_cg; 0, 1, 0; 0, 0, 1];
            A_3DOF = zeros(3, 3, N_omega);
            B_3DOF = zeros(3, 3, N_omega);
            for k = 1:N_omega
                AO = A_6_at_draft(idx_3dof, idx_3dof, k);  AO = 0.5*(AO + AO.');
                BO = B_6_at_draft(idx_3dof, idx_3dof, k);  BO = 0.5*(BO + BO.');
                ACG = T.' * AO * T;   A_3DOF(:, :, k) = 0.5*(ACG + ACG.');
                BCG = T.' * BO * T;   B_3DOF(:, :, k) = 0.5*(BCG + BCG.');
            end

            A_diag = zeros(3, N_omega);  B_diag = zeros(3, N_omega);
            for k = 1:N_omega
                A_diag(:, k) = diag(A_3DOF(:, :, k));
                B_diag(:, k) = diag(B_3DOF(:, :, k));
            end

            % A_inf at CG (full 3x3).  Prefer omega->inf WAMIT solution; else
            % use highest-omega BEM sample as proxy.  CG congruence applied.
            A_inf_3DOF_CG = zeros(3, 3);
            if isfield(hydro_table, 'A_inf') && ~isempty(hydro_table.A_inf)
                A_inf_cell_sorted = hydro_table.A_inf(idx_sort);
                if n_drafts == 1
                    A_inf_6 = A_inf_cell_sorted{1};
                else
                    A_inf_6 = (1 - w) * A_inf_cell_sorted{j_lo} + w * A_inf_cell_sorted{j_hi};
                end
                AOi = A_inf_6(idx_3dof, idx_3dof);  AOi = 0.5*(AOi + AOi.');
                tmp = T.' * AOi * T;  A_inf_3DOF_CG = 0.5*(tmp + tmp.');
                A_inf_source = 'hydro_table.A_inf (proper omega->inf)';
            else
                A_inf_3DOF_CG = A_3DOF(:, :, end);
                A_inf_source  = 'highest-omega BEM sample (proxy)';
            end
            A_inf_diag = diag(A_inf_3DOF_CG);

            % Complex excitation force (CG frame): Fe5_CG = Fe5 - z_cg*Fe1.
            % The inter-mode PHASE carries the surge-pitch interference, so Fe
            % must be complex (magnitude-only would be physically wrong).
            Fe_3DOF = nan(3, N_omega);
            Fe_meta = struct('present', false, 'cg_transfer_applied', false);
            if isfield(hydro_table, 'Fe') && ~isempty(hydro_table.Fe)
                Fe_cell_sorted = hydro_table.Fe(idx_sort);
                if n_drafts == 1
                    Fe_6_at_draft = Fe_cell_sorted{1};
                else
                    Fe_6_at_draft = (1 - w) * Fe_cell_sorted{j_lo} + w * Fe_cell_sorted{j_hi};
                end
                if isequal(size(Fe_6_at_draft), [6, N_omega])
                    Fe1 = Fe_6_at_draft(1, :);  Fe3 = Fe_6_at_draft(3, :);  Fe5 = Fe_6_at_draft(5, :);
                    Fe_3DOF = [Fe1; Fe3; Fe5 - z_cg .* Fe1];
                    Fe_meta.present = true;
                    Fe_meta.cg_transfer_applied = true;
                    Fe_meta.z_cg = z_cg;
                end
            end
            if ~Fe_meta.present
                error('MWEC_Tuning_Kernels:bem_at_design:noFe', ...
                      ['Excitation force Fe (complex, 6 x N_omega) is required ' ...
                       'for the coupled engine and is missing from the cache.']);
            end
            % Fe must carry phase: the surge-pitch interference lives in the
            % relative phase of Fe1 and Fe5.  A magnitude-only (real) Fe would
            % silently destroy the physics the coupled engine exists to capture.
            if max(abs(imag(Fe_3DOF(:)))) <= 1e-9 * max(abs(real(Fe_3DOF(:))))
                error('MWEC_Tuning_Kernels:bem_at_design:realFe', ...
                      ['Excitation force Fe is real-only (no phase).  The coupled ' ...
                       'engine needs the complex Fe (with inter-mode phase) from ' ...
                       'WAMIT; a magnitude-only cache is not usable.']);
            end
            Fe_phase_warn = false(3, 1);
            for k = 1:3
                if any(abs(diff(unwrap(angle(Fe_3DOF(k, :))))) > pi/2)
                    Fe_phase_warn(k) = true;
                end
            end

            M_diag = [M_body; M_body; Iyy];

            % Hydrostatic T_n via fixed-point (heave, pitch) — diagnostic/figure
            % markers and the G2c regression test only.  Not used by placement.
            fp_tol = gates.fp_tol;  fp_maxit = gates.fp_maxit;
            T_n_hydro = nan(3, 1);  fp_iters = zeros(3, 1);
            for k = 2:3
                Ck = K_hydro(k);  Mk = M_diag(k);
                if Ck <= 0, continue; end
                omega_guess = sqrt(Ck / (Mk + A_diag(k, 1)));
                converged = false;  omega_new = omega_guess;
                for it = 1:fp_maxit
                    A_at = interp1(omega_BEM, A_diag(k, :), omega_guess, 'pchip', NaN);
                    if ~isfinite(A_at), A_at = A_diag(k, 1); end
                    omega_new = sqrt(Ck / (Mk + A_at));
                    if abs(omega_new - omega_guess) < fp_tol * omega_guess
                        converged = true; break;
                    end
                    omega_guess = omega_new;
                end
                if ~converged
                    error('MWEC_Tuning_Kernels:bem_at_design:G2b', ...
                          'G2b FAIL: fixed-point T_n (mode %d) not converged in %d its.', k, fp_maxit);
                end
                fp_iters(k) = it;  T_n_hydro(k) = 2*pi / omega_new;
            end

            % G2c — pitch period regression (catches CG-transform corruption).
            % Compare BETTER of {fixed-point A(om_n), asymptotic A_inf}: the
            % FP-vs-asymptotic gap is a legitimate 5-10%, the bug is 30-40%.
            gate_G2c_pass = true;
            T_pitch_err_pct = NaN;  T_pitch_err_inf_pct = NaN;  T_pitch_asymptotic = NaN;
            if K_hydro(3) > 0 && (M_diag(3) + A_inf_diag(3)) > 0
                T_pitch_asymptotic = 2*pi / sqrt(K_hydro(3) / (M_diag(3) + A_inf_diag(3)));
            end
            if isfield(final_props, 'periods') && isstruct(final_props.periods) && ...
               isfield(final_props.periods, 'pitch') && isfinite(final_props.periods.pitch) && ...
               final_props.periods.pitch > 0
                T_ref = final_props.periods.pitch;
                if isfinite(T_n_hydro(3)) && T_n_hydro(3) > 0
                    T_pitch_err_pct = 100 * abs(T_n_hydro(3) - T_ref) / T_ref;
                end
                if isfinite(T_pitch_asymptotic) && T_pitch_asymptotic > 0
                    T_pitch_err_inf_pct = 100 * abs(T_pitch_asymptotic - T_ref) / T_ref;
                end
                best_err = min([T_pitch_err_pct, T_pitch_err_inf_pct]);
                if best_err > gates.T_pitch_pct_max
                    error('MWEC_Tuning_Kernels:bem_at_design:G2c', ...
                          ['G2c FAIL: neither pitch period matches mass file.\n' ...
                           '  ref=%.4f s | FP=%.4f s (%.2f%%) | asymp=%.4f s (%.2f%%) | gate %.1f%%.\n' ...
                           '  A mismatch under both comparisons indicates a CG-transform error.'], ...
                          T_ref, T_n_hydro(3), T_pitch_err_pct, T_pitch_asymptotic, ...
                          T_pitch_err_inf_pct, gates.T_pitch_pct_max);
                end
            end

            % Per-mode diagonal interpolants (B_PTO closure + K_total)
            A_kk_func = @(k, omega_q) interp1(omega_BEM, A_diag(k, :), omega_q(:), 'pchip', NaN);
            B_kk_func = @(k, omega_q) interp1(omega_BEM, B_diag(k, :), omega_q(:), 'pchip', NaN);

            % Surge-pitch coupling diagnostic: |A15|, |B15| normalised by the
            % geometric mean of the diagonals.  ~1 => rank-deficient shared
            % dipole block (the regime where surge & pitch genuinely compete).
            A15 = squeeze(A_3DOF(1, 3, :)).';   B15 = squeeze(B_3DOF(1, 3, :)).';
            denomA = sqrt(max(A_diag(1, :), 0) .* max(A_diag(3, :), 0));
            denomB = sqrt(max(B_diag(1, :), 0) .* max(B_diag(3, :), 0));
            coupling = struct( ...
                'A15', A15, 'B15', B15, ...
                'A15_norm', A15 ./ max(denomA, eps), ...
                'B15_norm', B15 ./ max(denomB, eps), ...
                'A15_norm_max', max(abs(A15 ./ max(denomA, eps))), ...
                'B15_norm_max', max(abs(B15 ./ max(denomB, eps))));

            bem = struct( ...
                'omega_BEM',       omega_BEM, ...
                'A_3DOF',          A_3DOF, ...
                'B_3DOF',          B_3DOF, ...
                'A_diag',          A_diag, ...
                'B_diag',          B_diag, ...
                'M_diag',          M_diag, ...
                'K_hydro',         K_hydro, ...
                'T_n_hydrostatic', T_n_hydro, ...
                'fp_iters',        fp_iters, ...
                'draft',           draft, ...
                'z_cg',            z_cg, ...
                'Fe_3DOF',         Fe_3DOF, ...
                'Fe_meta',         Fe_meta, ...
                'Fe_phase_warn',   Fe_phase_warn, ...
                'A_kk_func',       A_kk_func, ...
                'B_kk_func',       B_kk_func, ...
                'A_inf_3DOF_CG',   A_inf_3DOF_CG, ...
                'A_inf_diag',      A_inf_diag, ...
                'A_inf_source',    A_inf_source, ...
                'coupling',        coupling, ...
                'gate_G2a_pass',   true, ...
                'gate_G2b_pass',   true, ...
                'gate_G2c_pass',   gate_G2c_pass, ...
                'T_pitch_err_pct', T_pitch_err_pct, ...
                'T_pitch_err_inf_pct', T_pitch_err_inf_pct, ...
                'T_pitch_asymptotic',  T_pitch_asymptotic);

            if verbose
                fprintf('   [BEM] omega [%.3f, %.3f] rad/s  N=%d  draft=%.4f m  z_cg=%.4f m\n', ...
                        omega_BEM(1), omega_BEM(end), N_omega, draft, z_cg);
                fprintf('   [BEM] M=%.1f kg  Iyy=%.1f kg.m^2  A_inf diag=[%.2e %.2e %.2e]\n', ...
                        M_body, Iyy, A_inf_diag(1), A_inf_diag(2), A_inf_diag(3));
                fprintf('   [BEM] T_n_hydro: heave=%.3f s  pitch=%.3f s\n', T_n_hydro(2), T_n_hydro(3));
                fprintf('   [BEM] surge-pitch coupling: max|A15|/sqrt(A11.A33)=%.3f  max|B15|/sqrt(B11.B33)=%.3f\n', ...
                        coupling.A15_norm_max, coupling.B15_norm_max);
                if any(Fe_phase_warn)
                    fprintf('   [BEM] WARNING: Fe phase jumps > pi/2 in modes %s\n', mat2str(find(Fe_phase_warn)));
                end
            end
        end


        %% =================================================================
        %%  T3 — Stage-1 closure (per-mode PTO damping)
        %% =================================================================

        function clo = compute_closure(omega_n, bem, closure_cfg, A_mode)
        %COMPUTE_CLOSURE  Per-mode B_PTO and the silent K_total at a placement.
        %
        %   omega_n      [3x1]  candidate natural frequencies (surge,heave,pitch)
        %   closure_cfg  .mode = 'evans' | 'critical_damping', .zeta
        %   A_mode       'freq' (A(omega)) | 'inf' (A_inf) — for the sensitivity run
        %
        %   K_total,k = (M_k + A_kk(omega_n,k)) * omega_n,k^2   (silent; never output)
        %   Evans            : B_PTO,k = B_kk(omega_n,k)        resistance matching
        %   critical_damping : B_PTO,k = zeta*2*sqrt(K_total,k*(M_k+A_inf,k))
        %                      (defined for surge because K_total,k uses the placed
        %                      stiffness, not a hydrostatic one)
        %
        %   Gate G3: B_PTO,k finite and > 0 for all modes.
            if nargin < 4 || isempty(A_mode), A_mode = 'freq'; end
            omega_n = omega_n(:);
            M_diag  = bem.M_diag(:);

            A_n = zeros(3, 1);
            for k = 1:3
                if strcmpi(A_mode, 'inf')
                    A_n(k) = bem.A_inf_diag(k);
                else
                    A_n(k) = bem.A_kk_func(k, omega_n(k));
                    if ~isfinite(A_n(k)), A_n(k) = bem.A_inf_diag(k); end
                end
            end
            K_total = (M_diag + A_n) .* omega_n.^2;        % [silent internal]

            B_PTO = zeros(3, 1);
            switch lower(closure_cfg.mode)
                case 'evans'
                    for k = 1:3
                        Bk = bem.B_kk_func(k, omega_n(k));
                        if ~isfinite(Bk) || Bk <= 0
                            error('MWEC_Tuning_Kernels:compute_closure:G3', ...
                                  ['G3 FAIL (Evans): B_%d%d(omega_n=%.4f)=%g <= 0; ' ...
                                   'radiation resistance must be positive for matching.'], ...
                                  k, k, omega_n(k), Bk);
                        end
                        B_PTO(k) = Bk;                      % resistance matching
                    end
                    label = 'Evans (per-mode resistance matching B_PTO,k=B_kk(omega_n,k))';
                case 'critical_damping'
                    A_inf = bem.A_inf_diag(:);
                    for k = 1:3
                        arg = K_total(k) * (M_diag(k) + A_inf(k));
                        if ~(arg > 0)
                            error('MWEC_Tuning_Kernels:compute_closure:G3crit', ...
                                  ['G3 FAIL (crit): K_total(%d)*(M+A_inf)=%g <= 0.'], k, arg);
                        end
                        B_PTO(k) = closure_cfg.zeta * 2 * sqrt(arg);
                    end
                    label = sprintf('%g%% critical damping per mode (B_PTO,k=zeta*2*sqrt(K_total,k*(M_k+A_inf,k)))', ...
                                    100 * closure_cfg.zeta);
                otherwise
                    error('MWEC_Tuning_Kernels:compute_closure:badMode', ...
                          'Unknown closure mode ''%s'' (use ''evans'' or ''critical_damping'').', closure_cfg.mode);
            end

            % NOTE: K_total is the silent internal intermediate used above to set
            % the per-mode resonance/stiffness; per the contract it is the free
            % Stage-2 variable and is NOT exposed in the returned struct.
            clo = struct( ...
                'mode',           lower(closure_cfg.mode), ...
                'label',          label, ...
                'zeta',           closure_cfg.zeta, ...
                'omega_n',        omega_n, ...
                'A_at_n',         A_n, ...
                'B_PTO_per_mode', B_PTO);
        end


        %% =================================================================
        %%  T4 — Coupled response + placement
        %% =================================================================

        function R = coupled_rao(omega_q, omega_n, bem, B_PTO_per_mode, A_mode)
        %COUPLED_RAO  Full 3x3 steady-state frequency-domain solve x = Z^-1 Fe.
        %
        %   Z(omega) = -omega^2 (M + A(omega)) + i*omega (B(omega) + B_PTO) + K_total
        %   with the FULL 3x3 A,B (off-diagonals A15,B15 active).  B_PTO and
        %   K_total are diagonal (idealised independent per-mode PTO).
        %
        %   K_total realises resonance at omega_n,k (silent internal):
        %       K_total,k = (M_k + A_kk(omega_n,k)) * omega_n,k^2.
        %
        %   Returns X [3 x Nq] complex displacement RAO (per unit wave amplitude).
            if nargin < 5 || isempty(A_mode), A_mode = 'freq'; end
            omega_q = omega_q(:);  Nq = numel(omega_q);
            M_diag  = bem.M_diag(:);
            B_PTO   = B_PTO_per_mode(:);

            % K_total from the placement (diagonal A at omega_n)
            A_n = zeros(3, 1);
            for k = 1:3
                if strcmpi(A_mode, 'inf')
                    A_n(k) = bem.A_inf_diag(k);
                else
                    A_n(k) = bem.A_kk_func(k, omega_n(k));
                    if ~isfinite(A_n(k)), A_n(k) = bem.A_inf_diag(k); end
                end
            end
            K_total = (M_diag + A_n) .* omega_n(:).^2;
            Kmat    = diag(K_total);
            Mmat    = diag(M_diag);
            Bptomat = diag(B_PTO);

            % Interpolate full A,B and complex Fe onto omega_q
            if strcmpi(A_mode, 'inf')
                A_q = repmat(bem.A_inf_3DOF_CG, 1, 1, Nq);
            else
                A_q = MWEC_Tuning_Kernels.interp_matrix3(omega_q, bem.omega_BEM, bem.A_3DOF);
            end
            B_q  = MWEC_Tuning_Kernels.interp_matrix3(omega_q, bem.omega_BEM, bem.B_3DOF);
            Fe_q = MWEC_Tuning_Kernels.interp_fe3(omega_q, bem.omega_BEM, bem.Fe_3DOF);  % 3 x Nq

            X = complex(zeros(3, Nq));
            singular = false;
            for q = 1:Nq
                om = omega_q(q);
                Aq = A_q(:, :, q);  Bq = B_q(:, :, q);
                Z = -om^2 * (Mmat + Aq) + 1i * om * (Bq + Bptomat) + Kmat;
                if rcond(Z) < 1e-14 || any(~isfinite(Z(:)))
                    singular = true;
                    X(:, q) = 0;
                else
                    X(:, q) = Z \ Fe_q(:, q);
                end
            end

            R = struct('omega', omega_q, 'X', X, 'K_total', K_total, ...
                       'A_at_n', A_n, 'singular', singular);
        end


        function res = evaluate_placement(omega_n, bem, closure_cfg, omega, S_ew, band, rho, g_acc, A_mode)
        %EVALUATE_PLACEMENT  eta_Falnes + <P_abs> + capture width at a placement.
        %
        %   eta_Falnes = INT_band sum_k p_k S domega / INT_band CW_max F_ew domega
        %     p_k(omega) = B_PTO,k * omega^2 * |X_k|^2      (absorbed power density)
        %     CW_max(omega) = 3 g / omega^2                 (heave nu=1 + dipole nu=2)
        %   Falnes denominator closed form: (3*rho*g^3/2) INT_band S/omega^3 domega.
        %   Dimensionless in [0,1]; same argmax as raw <P_abs> (denominator is
        %   placement-independent).  <P_abs> integrated over the full omega grid.
            if nargin < 9 || isempty(A_mode), A_mode = 'freq'; end
            omega_n = omega_n(:);
            omega   = omega(:);  S_ew = S_ew(:);

            clo = MWEC_Tuning_Kernels.compute_closure(omega_n, bem, closure_cfg, A_mode);
            R   = MWEC_Tuning_Kernels.coupled_rao(omega, omega_n, bem, clo.B_PTO_per_mode, A_mode);

            B_PTO = clo.B_PTO_per_mode(:);
            p = zeros(numel(omega), 3);                    % power density per mode
            for k = 1:3
                p(:, k) = B_PTO(k) .* omega.^2 .* abs(R.X(k, :)).'.^2;
            end

            % <P_abs> over the full spectrum
            P_per_mode = zeros(3, 1);
            for k = 1:3
                P_per_mode(k) = trapz(omega, p(:, k) .* S_ew);
            end
            P_total = sum(P_per_mode);

            % eta_Falnes over the placement band
            mask = (omega >= band.omega_L) & (omega <= band.omega_H);
            om_b = omega(mask);  S_b = S_ew(mask);
            num  = trapz(om_b, sum(p(mask, :), 2) .* S_b);
            den  = (3 * rho * g_acc^3 / 2) * trapz(om_b, S_b ./ om_b.^3);
            if den <= 0
                eta_Falnes = 0;
            else
                eta_Falnes = num / den;
            end

            % Capture width (band-integrated flux as reference)
            F_ew    = (rho * g_acc^2 / 2) .* S_ew ./ max(omega, 1e-10);
            J_total = trapz(omega, F_ew);                  % wave power per crest width [W/m]
            CW      = P_total / max(J_total, eps);         % [m]

            res = struct( ...
                'omega_n',        omega_n, ...
                'eta_Falnes',     eta_Falnes, ...
                'eta_num',        num, ...
                'eta_den',        den, ...
                'P_abs_per_mode', P_per_mode, ...   % over full spectrum
                'P_abs_total',    P_total, ...      % over full spectrum (reported headline)
                'P_abs_band',     num, ...          % over 90% band = exactly what eta optimises
                'CW',             CW, ...
                'J_total',        J_total, ...
                'closure',        clo, ...
                'singular',       R.singular);
        end


        function eta = eta_objective(omega_n, bem, closure_cfg, om_b, S_b, rho, g_acc, A_mode)
        %ETA_OBJECTIVE  Band-only eta_Falnes — the fast objective for the optimiser.
        %   om_b, S_b are the band-restricted omega and S_ew (precomputed by the
        %   caller, so each evaluation solves the coupled RAO only over the
        %   placement band).  Returns the scalar eta_Falnes in [0,1].
            omega_n = omega_n(:);
            clo = MWEC_Tuning_Kernels.compute_closure(omega_n, bem, closure_cfg, A_mode);
            R   = MWEC_Tuning_Kernels.coupled_rao(om_b, omega_n, bem, clo.B_PTO_per_mode, A_mode);
            B_PTO = clo.B_PTO_per_mode(:);
            psum = zeros(numel(om_b), 1);
            for k = 1:3
                psum = psum + B_PTO(k) .* om_b.^2 .* abs(R.X(k, :)).'.^2 .* S_b;
            end
            num = trapz(om_b, psum);
            den = (3 * rho * g_acc^3 / 2) * trapz(om_b, S_b ./ om_b.^3);
            if den <= 0, eta = 0; else, eta = num / den; end
        end


        function [omega_n_opt, eta_opt, info] = place_natural_periods( ...
                bem, closure_cfg, omega, S_ew, band, omega_bounds, rho, g_acc, opt_cfg, A_mode)
        %PLACE_NATURAL_PERIODS  GA -> SQP maximisation of eta_Falnes.
        %
        %   Decision variables: omega_n [3x1] (the design natural frequencies).
        %   omega_bounds [3x2] per-mode [lb,ub]: surge low sub-band, heave full,
        %   pitch high sub-band — the partition is the surge/pitch separation
        %   backstop.  No FWHM constraint, no Picard: the coupled objective is a
        %   direct function of omega_n, and separation emerges from it.
            if nargin < 10 || isempty(A_mode), A_mode = 'freq'; end
            verbose = opt_cfg.verbose;

            lb = max(omega_bounds(:, 1), band.omega_L);
            ub = min(omega_bounds(:, 2), band.omega_H);
            lb = lb(:);  ub = ub(:);

            % Band-restricted grid (the objective integrates over [omega_L,omega_H])
            omega = omega(:);  S_ew = S_ew(:);
            mask = (omega >= band.omega_L) & (omega <= band.omega_H);
            om_b = omega(mask);  S_b = S_ew(mask);

            obj_n = @(xv) -MWEC_Tuning_Kernels.eta_objective( ...
                MWEC_Tuning_Kernels.clamp_x(xv, lb, ub), ...
                bem, closure_cfg, om_b, S_b, rho, g_acc, A_mode);

            x0 = (lb + ub) / 2;
            x_best = x0;  eta_best = -obj_n(x0);

            ga_available = license('test', 'GADS_Toolbox') && exist('ga', 'file') == 2;
            if ga_available
                if verbose, fprintf('     Stage 1 (GA): pop=%d gen=%d...\n', ...
                        opt_cfg.ga.PopulationSize, opt_cfg.ga.MaxGenerations); end
                opts_ga = optimoptions('ga', ...
                    'PopulationSize',          opt_cfg.ga.PopulationSize, ...
                    'MaxGenerations',          opt_cfg.ga.MaxGenerations, ...
                    'FunctionTolerance',       1e-8, ...
                    'InitialPopulationMatrix', x0', ...
                    'UseParallel',             false, ...
                    'Display',                 'off');
                try
                    [x_ga, fval_ga] = ga(@(xv) obj_n(xv(:)), 3, [], [], [], [], lb, ub, [], opts_ga);
                    if -fval_ga > eta_best, x_best = x_ga(:); eta_best = -fval_ga; end
                catch ME
                    warning('MWEC_Tuning_Kernels:place:gaFailed', ...
                            'GA failed (%s); using multi-start.', ME.message);
                    ga_available = false;
                end
            end
            if ~ga_available
                if verbose, fprintf('     Stage 1 (multi-start fminsearch)\n'); end
                span = ub - lb;
                starts = {x0, lb + [0.25; 0.50; 0.80].*span, ...
                              lb + [0.15; 0.45; 0.85].*span, ...
                              lb + [0.35; 0.55; 0.75].*span};
                fmin_opts = optimset('Display', 'off', 'TolX', opt_cfg.fminsearch.TolX, ...
                                     'TolFun', 1e-7, 'MaxIter', 600);
                for s = 1:numel(starts)
                    try
                        [x_s, f_s] = fminsearch(@(xv) obj_n(xv(:)), starts{s}, fmin_opts);
                        if -f_s > eta_best, x_best = MWEC_Tuning_Kernels.clamp_x(x_s, lb, ub); eta_best = -f_s; end
                    catch
                        continue
                    end
                end
            end

            stage1_eta = eta_best;
            x_best = MWEC_Tuning_Kernels.clamp_x(x_best, lb, ub);

            sqp_available = license('test', 'Optimization_Toolbox') && exist('fmincon', 'file') == 2;
            if sqp_available
                if verbose, fprintf('     Stage 2 (SQP)\n'); end
                opts_sqp = optimoptions('fmincon', 'Algorithm', 'sqp', ...
                    'MaxIterations', 800, 'FunctionTolerance', opt_cfg.fmincon.TolFun, ...
                    'StepTolerance', 1e-11, 'Display', 'off');
                try
                    [x_sqp, fval_sqp] = fmincon(@(xv) obj_n(xv(:)), x_best, [], [], [], [], lb, ub, [], opts_sqp);
                    if -fval_sqp > eta_best, x_best = x_sqp(:); eta_best = -fval_sqp; end
                catch ME
                    warning('MWEC_Tuning_Kernels:place:sqpFailed', ...
                            'fmincon failed (%s); keeping Stage-1.', ME.message);
                    sqp_available = false;
                end
            end

            omega_n_opt = MWEC_Tuning_Kernels.clamp_x(x_best, lb, ub);
            eta_opt     = eta_best;

            % Separation backstop status: does pitch sit on the partition floor
            % (lb of pitch == ub of surge == partition)?  Binding => the coupled
            % objective wanted them closer than the partition allows.
            tol_bind = 1e-3;
            surge_at_ub = abs(omega_n_opt(1) - ub(1)) < tol_bind * max(ub(1), 1);
            pitch_at_lb = abs(omega_n_opt(3) - lb(3)) < tol_bind * max(lb(3), 1);
            sep_binding = surge_at_ub || pitch_at_lb;

            info = struct( ...
                'omega_n_opt',   omega_n_opt, ...
                'T_n_opt',       2*pi ./ omega_n_opt, ...
                'eta_opt',       eta_opt, ...
                'stage1_eta',    stage1_eta, ...
                'ga_available',  ga_available, ...
                'sqp_available', sqp_available, ...
                'lb',            lb, 'ub', ub, ...
                'sep_binding',   sep_binding, ...
                'omega_bounds',  omega_bounds);
        end


        %% =================================================================
        %%  T5 — Capture-width diagnostic (per mode vs Falnes ceiling)
        %% =================================================================

        function tbl = absorption_vs_T(T_grid, bem, closure_at_opt, omega_n_opt, rho, g_acc, L_ref, gates, A_mode)
        %ABSORPTION_VS_T  Per-mode capture width CW_k(T) [m] vs Falnes ceiling.
        %
        %   CW_k(omega) = 2 B_PTO,k omega^3 |X_k(omega)|^2 / (rho g^2)   (coupled X)
        %   Ceilings (Falnes 2002 Table 6.1, deep water):
        %       heave (monopole)      : g/omega^2
        %       surge+pitch (dipole)  : 2 g/omega^2   COMBINED — counted once
        %   G5a: CW_heave <= g/omega^2 ; CW_surge + CW_pitch <= 2 g/omega^2.
            if nargin < 9 || isempty(A_mode), A_mode = 'freq'; end
            T_grid     = T_grid(:);
            omega_grid = 2*pi ./ T_grid;
            B_PTO      = closure_at_opt.B_PTO_per_mode(:);

            R = MWEC_Tuning_Kernels.coupled_rao(omega_grid, omega_n_opt, bem, B_PTO, A_mode);
            X = R.X;                                        % 3 x N

            CW = zeros(numel(T_grid), 3);
            for k = 1:3
                CW(:, k) = 2 .* B_PTO(k) .* omega_grid.^3 .* abs(X(k, :)).'.^2 ./ max(rho * g_acc^2, eps);
            end

            ceil_heave  = g_acc ./ omega_grid.^2;           % nu=1
            ceil_dipole = 2 * g_acc ./ omega_grid.^2;       % nu=2 combined

            viol_heave  = max(0, CW(:, 2) - ceil_heave);
            viol_dipole = max(0, (CW(:, 1) + CW(:, 3)) - ceil_dipole);
            max_viol_pct = 100 * max([ max(viol_heave  ./ max(ceil_heave,  eps)); ...
                                       max(viol_dipole ./ max(ceil_dipole, eps)) ]);
            if max_viol_pct > gates.l_falnes_pct_max
                warning('MWEC_Tuning_Kernels:absorption:G5a', ...
                        ['G5a WARN: CW exceeds the Falnes ceiling by up to %.2f%% ' ...
                         '(gate %.1f%%).  Numerical artefact near an irregular frequency?'], ...
                        max_viol_pct, gates.l_falnes_pct_max);
            end

            tbl = struct( ...
                'T',              T_grid, ...
                'omega',          omega_grid, ...
                'CW_surge',       CW(:, 1), ...
                'CW_heave',       CW(:, 2), ...
                'CW_pitch',       CW(:, 3), ...
                'CW_dipole_sum',  CW(:, 1) + CW(:, 3), ...
                'ceil_heave',     ceil_heave, ...
                'ceil_dipole',    ceil_dipole, ...
                'CWR_surge',      CW(:, 1) ./ max(L_ref, eps) * 100, ...
                'CWR_heave',      CW(:, 2) ./ max(L_ref, eps) * 100, ...
                'CWR_pitch',      CW(:, 3) ./ max(L_ref, eps) * 100, ...
                'B_PTO_per_mode', B_PTO, ...
                'L_ref',          L_ref, ...
                'max_viol_pct',   max_viol_pct);
        end


        %% =================================================================
        %%  Sensitivity — re-place with A(omega) -> A_inf
        %% =================================================================

        function comp = sensitivity_A_inf(bem, closure_cfg, omega, S_ew, band, omega_bounds, ...
                                          rho, g_acc, opt_cfg, L_ref, gates)
        %SENSITIVITY_A_INF  Re-run placement with A(omega) replaced by A_inf.
        %   B(omega) kept frequency-dependent so the comparison isolates the
        %   A-frequency-dependence effect (broadband-robustness finding).
            try
                [omega_n_inf, eta_inf, info_inf] = MWEC_Tuning_Kernels.place_natural_periods( ...
                    bem, closure_cfg, omega, S_ew, band, omega_bounds, rho, g_acc, opt_cfg, 'inf');
                res_inf = MWEC_Tuning_Kernels.evaluate_placement( ...
                    omega_n_inf, bem, closure_cfg, omega, S_ew, band, rho, g_acc, 'inf');
            catch ME
                comp = struct('error', sprintf('A_inf placement failed: %s', ME.message));
                return;
            end
            comp = struct( ...
                'A_inf_source',  bem.A_inf_source, ...
                'omega_n',       omega_n_inf, ...
                'T_n_placed_s',  2*pi ./ omega_n_inf, ...
                'eta_Falnes',    eta_inf, ...
                'P_abs_per_mode',res_inf.P_abs_per_mode, ...
                'P_abs_total',   res_inf.P_abs_total, ...
                'CW',            res_inf.CW, ...
                'CWR',           res_inf.CW / max(L_ref, eps) * 100, ...
                'sep_binding',   info_inf.sep_binding);
        end


        %% =================================================================
        %%  Utilities
        %% =================================================================

        function [L_ref, src] = resolve_L_ref(final_props, override)
        %RESOLVE_L_REF  Characteristic width (raft beam) for CWR normalisation.
            if ~isempty(override) && isfinite(override) && override > 0
                L_ref = double(override);  src = 'cfg override (manual)';  return;
            end
            if isfield(final_props, 'cross_section')
                cs = final_props.cross_section;
                L_ref = max(cs(:, 1)) - min(cs(:, 1));
                src   = 'final_props.cross_section col 1 (beam y-extent)';
            else
                error('MWEC_Tuning_Kernels:Lref:missing', ...
                      'final_props.cross_section missing; provide cfg.absorption.L_ref override.');
            end
        end


        function xc = clamp_x(xv, lb, ub)
            xv = xv(:);
            xc = min(max(xv, lb(:)), ub(:));
        end


        function Mq = interp_matrix3(omega_q, omega_BEM, M3)
        %INTERP_MATRIX3  pchip-interpolate a real 3x3xN array onto omega_q.
        %   Returns 3 x 3 x numel(omega_q).
            omega_q = omega_q(:);  Nq = numel(omega_q);
            Mq = zeros(3, 3, Nq);
            for a = 1:3
                for b = 1:3
                    seq = squeeze(M3(a, b, :));
                    Mq(a, b, :) = interp1(omega_BEM, seq, omega_q, 'pchip', 'extrap');
                end
            end
        end


        function Feq = interp_fe3(omega_q, omega_BEM, Fe3)
        %INTERP_FE3  pchip-interpolate complex 3xN excitation onto omega_q.
        %   Real and imaginary parts interpolated separately (preserves phase).
            omega_q = omega_q(:).';  Nq = numel(omega_q);
            Feq = complex(zeros(3, Nq));
            for k = 1:3
                re = interp1(omega_BEM, real(Fe3(k, :)), omega_q, 'pchip', 0);
                im = interp1(omega_BEM, imag(Fe3(k, :)), omega_q, 'pchip', 0);
                Feq(k, :) = re + 1i * im;
            end
        end

    end
end
