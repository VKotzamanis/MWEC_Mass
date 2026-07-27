classdef WEC_File_IO
    % WEC_FILE_IO  WAMIT data reader with denormalisation and coordinate transform.
    %
    %   This static class handles the full pipeline from a WAMIT .1 file
    %   (dimensionless coefficients at the body origin) to physical added-mass
    %   and radiation-damping matrices referenced to the centre of gravity:
    %
    %     .1 file  →  parse  →  denormalise  →  transform to CG  →  verify
    %
    %   METHOD INVENTORY
    %   ────────────────────────────────────────────────────────────────────
    %   PUBLIC
    %     read_wamit_data                    — top-level reader (entry point)
    %     denormalize_coefficients           — WAMIT dimensionless → SI units
    %     transform_hydrodynamic_matrices_3x3 — body-origin → CG reference
    %     verify_transformation_3x3          — analytic spot-checks on transform
    %
    %   PRIVATE
    %     print_matrix_3x3   — formatted 3×3 console dump
    %     print_pass_fail    — PASS / FAIL label from error vs tolerance
    %   ────────────────────────────────────────────────────────────────────
    %
    %   COORDINATE TRANSFORMATION (body origin → CG)
    %
    %     WAMIT outputs A and B at the body origin (0,0,0).  The optimiser
    %     needs them at the CG (0, 0, Z_CG).
    %
    %     Kinematic relationship for planar surge–heave–pitch:
    %       q_origin = T · q_CG
    %
    %       T = [ 1   0  −Z_CG ]       T⁻¹ = [ 1   0   Z_CG ]
    %           [ 0   1   0    ]              [ 0   1   0    ]
    %           [ 0   0   1    ]              [ 0   0   1    ]
    %
    %     Bilinear-form transformation:
    %       M_CG = T' · M_origin · T
    %
    %   RADIATION DAMPING FILTERING
    %     Only periods in [4, 16] s are averaged for B_rad.
    %     This covers the MWEC operating envelope:
    %       pitch ∼3 s, heave ∼7 s, surge ∼15 s.
    %     [VERIFY] T_min = 4 s and T_max = 16 s are currently hardcoded
    %       inside read_wamit_data.  Consider moving to the driver as
    %       config.brad_period_range if different hull families need
    %       different bands.
    %
    %   See also: WEC_Configuration_Builder, WEC_Core_Functions
    %
    %   Author:  WEC Optimisation Team
    %   Version: 4.0 — fixed congruence transform (T' M T, not T⁻ᵀ M T⁻¹)
    
    methods (Static)

        %% ═════════════════════════════════════════════════════════════
        %%  TOP-LEVEL READER
        %% ═════════════════════════════════════════════════════════════

        function wamit_data = read_wamit_data(file_prefix, rho, g, L, Z_CG)
            % READ_WAMIT_DATA  Read a WAMIT .1 file and return CG-referenced
            %   added-mass and radiation-damping matrices.
            %
            %   wamit_data = WEC_File_IO.read_wamit_data(prefix, rho, g, L, Z_CG)
            %
            %   PIPELINE
            %     1. Parse .1 file → raw dimensionless coefficient tensors
            %     2. Denormalise → physical SI units
            %     3. Transform from body origin to CG
            %     4. Verify transform with analytic spot-checks
            %     5. Extract infinite-frequency A and average B over [4,16] s
            %
            %   VARIABLE DICTIONARY (scope: this function)
            %     periods          [s]      column 1 of .1 file (0 = ω → ∞)
            %     i_modes, j_modes [-]      WAMIT mode indices (1,3,5)
            %     added_mass_coeff [-]      col 4, dimensionless Ā_ij
            %     radiation_damping[-]      col 5, dimensionless B̄_ij
            %     mode_map         [-]      {1→1, 3→2, 5→3} compact index
            %     A_raw, B_raw     [-]      [3×3×N] raw coefficient tensors
            %     A_denorm, B_denorm [SI]   after denormalisation
            %     A_cg, B_cg       [SI]    after coordinate transform to CG
            %     inf_freq_idx     [-]      index of period ≤ 0 (WAMIT convention)
            %     T_min, T_max     [s]      radiation-damping averaging band
            %     B_rad_avg        [SI]     mean B over [T_min, T_max]
            %
            %   INPUTS
            %     file_prefix : char   path without extension (e.g. 'data/case1')
            %     rho         : [kg/m³]  water density
            %     g           : [m/s²]   gravitational acceleration
            %     L           : [m]      WAMIT characteristic length
            %     Z_CG        : [m]      CG z-coordinate (negative if below origin)
            %
            %   OUTPUT
            %     wamit_data : struct with fields
            %       .A           [3×3]      infinite-freq added mass at CG
            %       .B_rad_avg   [3×3]      averaged radiation damping at CG
            %       .A_full      [3×3×N]    all-frequency A at CG
            %       .B_full      [3×3×N]    all-frequency B at CG
            %       .periods     [N×1 s]    unique periods from .1 file
            %       .periods_used[M×1 s]    periods inside [T_min, T_max]
            
            try
                filename = [file_prefix, '.1'];
                
                % Open file with error checking
                fid = fopen(filename, 'r');
                if fid == -1
                    error('Cannot open WAMIT file: %s', filename);
                end
                
                % Read header and data
                header = fgetl(fid); %#ok<NASGU>
                data = textscan(fid, '%f %d %d %f %f');
                fclose(fid);
                
                periods = data{1};
                i_modes = data{2};
                j_modes = data{3};
                added_mass_coeff = data{4};
                radiation_damping = data{5};
                
                % Extract modes of interest [surge=1, heave=3, pitch=5]
                modes_of_interest = [1, 3, 5];
                mode_map = containers.Map({1, 3, 5}, {1, 2, 3});
                
                unique_periods = unique(periods);
                n_periods = length(unique_periods);
                
                % Preallocate matrices
                A_raw = zeros(3, 3, n_periods);
                B_raw = zeros(3, 3, n_periods);
                
                % Populate raw coefficient matrices
                for idx = 1:length(periods)
                    period = periods(idx);
                    i_mode = i_modes(idx);
                    j_mode = j_modes(idx);
                    
                    if ismember(i_mode, modes_of_interest) && ismember(j_mode, modes_of_interest)
                        period_idx = find(unique_periods == period, 1);
                        i_idx = mode_map(i_mode);
                        j_idx = mode_map(j_mode);
                        
                        A_raw(i_idx, j_idx, period_idx) = added_mass_coeff(idx);
                        
                        if ~isnan(radiation_damping(idx))
                            B_raw(i_idx, j_idx, period_idx) = radiation_damping(idx);
                        end
                    end
                end
                
                % Denormalize coefficients
                [A_denorm, B_denorm] = WEC_File_IO.denormalize_coefficients(...
                    A_raw, B_raw, unique_periods, rho, g, L);
                
                % Apply CORRECTED coordinate transformation from body origin to CG
                fprintf('  Applying CORRECTED coordinate transformation: Body Origin -> CG (Z_CG = %.4f m)\n', Z_CG);
                [A_cg, B_cg] = WEC_File_IO.transform_hydrodynamic_matrices_3x3(...
                    A_denorm, B_denorm, Z_CG);
                
                % Verify transformation
                WEC_File_IO.verify_transformation_3x3(...
                    A_denorm, A_cg, B_denorm, B_cg, Z_CG);
                
                % WHY period ≤ 0 for infinite frequency?
                %   WAMIT convention: period = 0 (or negative) in the .1
                %   file represents ω → ∞.  The added mass at infinite
                %   frequency (A∞) is the constant that appears in the
                %   Cummins equation:  m_eff = m + A∞.  It is always the
                %   first (or last) record in the file.
                inf_freq_idx = find(unique_periods <= 0, 1);
                if isempty(inf_freq_idx)
                    warning('WEC_File_IO:NoInfFreq', 'No infinite frequency data, using highest frequency (shortest period)');
                    [~, inf_freq_idx] = min(unique_periods(unique_periods > 0));
                end
                A = A_cg(:, :, inf_freq_idx);
                
                % WHY filter to [4, 16] s?
                %   Radiation damping is frequency-dependent.  Averaging
                %   over all frequencies would dilute the physically
                %   relevant values with high-frequency (short-period)
                %   data where B → 0, and low-frequency (long-period) data
                %   where B also → 0.  The [4, 16] s band brackets the
                %   three DOF natural periods of interest:
                %     pitch ∼ 3 s,  heave ∼ 7 s,  surge ∼ 15 s
                %   with 1 s margin on each side.
                %
                % [VERIFY] T_min and T_max are hardcoded here.  If a new
                %   hull family has significantly different natural periods,
                %   these should move to config.brad_period_range.
                T_min = 4.0;   % [s]
                T_max = 16.0;  % [s]
                
                period_mask = (unique_periods >= T_min) & (unique_periods <= T_max);
                
                if sum(period_mask) == 0
                    warning('WEC_File_IO:NoPeriodInRange', ...
                            'No periods in %.1f-%.1fs range. Using closest to 7.5s.', T_min, T_max);
                    positive_periods = unique_periods(unique_periods > 0);
                    [~, closest_idx] = min(abs(positive_periods - 7.5));
                    period_mask = (unique_periods == positive_periods(closest_idx));
                end
                
                periods_filtered = unique_periods(period_mask);
                B_filtered = B_cg(:, :, period_mask);
                
                % Average radiation damping over filtered range
                B_rad_avg = mean(B_filtered, 3);
                
                % Symmetrize the averaged matrix
                B_rad_avg = 0.5 * (B_rad_avg + B_rad_avg');
                
                fprintf('  Radiation damping averaged over %d periods: [%.1fs, %.1fs]\n', ...
                        sum(period_mask), min(periods_filtered), max(periods_filtered));
                
                % Diagnostic output: show the matrices
                fprintf('\n  === DIAGNOSTIC: CG-Referenced Matrices ===\n');
                fprintf('  A_inf (3x3 at CG):\n');
                WEC_File_IO.print_matrix_3x3(A, '    ');
                fprintf('  B_rad_avg (3x3 at CG):\n');
                WEC_File_IO.print_matrix_3x3(B_rad_avg, '    ');
                
                % Return processed data with FULL matrices
                wamit_data = struct('A', A, ...
                                    'B_rad_avg', B_rad_avg, ...
                                    'A_full', A_cg, ...
                                    'B_full', B_cg, ...
                                    'periods', unique_periods, ...
                                    'periods_used', periods_filtered);
                
                fprintf('  SUCCESS: WAMIT data processed with CORRECTED transformations\n\n');
                
            catch ME
                warning('WEC_File_IO:ReadFailed', ...
                        'WAMIT data read failed: %s. Returning default values.', ME.message);
                wamit_data = struct('A', zeros(3,3), ...
                                    'B_rad_avg', zeros(3,3), ...
                                    'A_full', zeros(3,3,1), ...
                                    'B_full', zeros(3,3,1), ...
                                    'periods', [], ...
                                    'periods_used', []);
            end
        end
        
        %% ═════════════════════════════════════════════════════════════
        %%  DENORMALISATION
        %% ═════════════════════════════════════════════════════════════

        function [A_denorm, B_denorm] = denormalize_coefficients(A_raw, B_raw, periods, rho, g, L)
            % DENORMALIZE_COEFFICIENTS  Convert WAMIT dimensionless coefficients
            %   to physical SI units.
            %
            %   [A_denorm, B_denorm] = WEC_File_IO.denormalize_coefficients(...)
            %
            %   WAMIT NORMALISATION CONVENTION
            %     Ā_ij = A_ij / (ρ · L^k(i,j))
            %     B̄_ij = B_ij / (ρ · ω · L^k(i,j))
            %
            %   where k(i,j) is the length-power exponent matrix for the
            %   compact mode set [surge=1, heave=3, pitch=5] → [1,2,3]:
            %
            %     k = [ 3  3  4 ]   surge–surge  surge–heave  surge–pitch
            %         [ 3  3  4 ]   heave–surge  heave–heave  heave–pitch
            %         [ 4  4  5 ]   pitch–surge  pitch–heave  pitch–pitch
            %
            %   WHY this exponent pattern?
            %     Translational DOFs (surge, heave) scale as L³ (volume),
            %     rotational DOFs (pitch) scale as L⁴ (volume × lever).
            %     Off-diagonal coupling between translation and rotation
            %     scales as L^((3+4)/2) but WAMIT rounds to the larger
            %     exponent, so surge–pitch = 4 and pitch–pitch = 5.
            %
            %   WHY B = 0 at period = 0?
            %     Period = 0 means ω → ∞.  At infinite frequency the
            %     free-surface boundary condition becomes a rigid wall;
            %     no waves are radiated, so damping is exactly zero.
            %
            %   INPUTS
            %     A_raw, B_raw : [3×3×N]  dimensionless coefficients
            %     periods      : [N×1 s]  wave periods (0 = ω → ∞)
            %     rho, g, L    : normalisation constants
            
            try
                n_periods = size(A_raw, 3);
                
                % Power exponents for each matrix element
                % For modes [surge=1, heave=3, pitch=5] mapped to [1,2,3]
                k_matrix = [3, 3, 4;
                            3, 3, 4;
                            4, 4, 5];
                
                A_denorm = zeros(size(A_raw));
                B_denorm = zeros(size(B_raw));
                
                for p_idx = 1:n_periods
                    period = periods(p_idx);
                    
                    % Angular frequency (handle infinite frequency case)
                    if period <= 0
                        omega = inf;  % Will be handled specially for B
                    else
                        omega = 2 * pi / period;
                    end
                    
                    for ii = 1:3
                        for jj = 1:3
                            k = k_matrix(ii, jj);
                            scale_A = rho * L^k;
                            
                            A_denorm(ii, jj, p_idx) = A_raw(ii, jj, p_idx) * scale_A;
                            
                            % B is zero at infinite frequency (period = 0)
                            if period > 0
                                scale_B = rho * omega * L^k;
                                B_denorm(ii, jj, p_idx) = B_raw(ii, jj, p_idx) * scale_B;
                            else
                                B_denorm(ii, jj, p_idx) = 0;
                            end
                        end
                    end
                end
                
            catch ME
                warning('WEC_File_IO:DenormFailed', ...
                        'Denormalization failed: %s. Returning raw values.', ME.message);
                A_denorm = A_raw;
                B_denorm = B_raw;
            end
        end
        
        %% ═════════════════════════════════════════════════════════════
        %%  COORDINATE TRANSFORMATION (body origin → CG)
        %% ═════════════════════════════════════════════════════════════

        function [A_cg, B_cg] = transform_hydrodynamic_matrices_3x3(A_origin, B_origin, Z_CG)
            % TRANSFORM_HYDRODYNAMIC_MATRICES_3X3  Shift the reference
            %   point of added-mass and damping matrices from the body
            %   origin (0,0,0) to the centre of gravity (0, 0, Z_CG).
            %
            %   [A_cg, B_cg] = WEC_File_IO.transform_hydrodynamic_matrices_3x3(...)
            %
            %   DERIVATION (surge–heave–pitch planar sub-system)
            %
            %     q_origin = T · q_CG
            %
            %       T = [ 1   0  −Z_CG ]       T⁻¹ = [ 1   0   Z_CG ]
            %           [ 0   1   0    ]              [ 0   1   0    ]
            %           [ 0   0   1    ]              [ 0   0   1    ]
            %
            %     Energy invariance:
            %       E = ½ q_o' M_o q_o = ½ (T q_CG)' M_o (T q_CG)
            %         = ½ q_CG' (T' M_o T) q_CG
            %
            %     We have M_o and want M_CG.  Since q_o = T q_CG:
            %       M_CG = T' · M_o · T
            %
            %   WHY symmetrise both before and after the transform?
            %     WAMIT output can have O(ε) asymmetry from finite-panel
            %     resolution.  Pre-symmetrising ensures the bilinear-form
            %     identity M_CG = M_CG' holds exactly rather than to O(ε²).
            %     Post-symmetrising catches any residual from the matrix
            %     multiplications.
            %
            %   INPUTS
            %     A_origin, B_origin : [3×3×N]  matrices at body origin
            %     Z_CG               : [m]      CG z-coordinate (typically < 0)
            %
            %   OUTPUTS
            %     A_cg, B_cg : [3×3×N]  matrices at CG
            
            try
                [n_i, n_j, n_periods] = size(A_origin);
                
                % CORRECTED Transformation matrix
                % Maps CG velocities to body origin velocities: q_origin = T * q_CG
                T = [1,  0, -Z_CG;
                     0,  1,  0;
                     0,  0,  1];
                
                % Inverse transformation: q_CG = Tinv * q_origin
                Tinv = [1,  0, Z_CG;
                        0,  1,  0;
                        0,  0,  1];
                
                % Verify inverse
                T_check = T * Tinv;
                if max(max(abs(T_check - eye(3)))) > 1e-12
                    warning('WEC_File_IO:InverseError', 'T * Tinv != I');
                end
                
                A_cg = zeros(n_i, n_j, n_periods);
                B_cg = zeros(n_i, n_j, n_periods);
                
                for p_idx = 1:n_periods
                    % Symmetrize input matrices (WAMIT may have small asymmetries)
                    A_symm = 0.5 * (A_origin(:, :, p_idx) + A_origin(:, :, p_idx)');
                    B_symm = 0.5 * (B_origin(:, :, p_idx) + B_origin(:, :, p_idx)');
                    
                    % Congruence transform: M_CG = T' * M_origin * T
                    %   From energy invariance: E = q_O' M_O q_O
                    %     = (T q_CG)' M_O (T q_CG) = q_CG' (T' M_O T) q_CG
                    %   Therefore M_CG = T' M_O T.
                    %   T has −z_G in position (1,3);  T⁻¹ has +z_G.
                    %   Using T⁻¹ here would flip the sign of the 2·z_G·A₁₅
                    %   cross-term in A₅₅ — an error of 4·z_G·A₁₅.
                    A_transformed = T' * A_symm * T;
                    B_transformed = T' * B_symm * T;
                    
                    % Re-symmetrize (should be exact, but enforce for numerics)
                    A_cg(:, :, p_idx) = 0.5 * (A_transformed + A_transformed');
                    B_cg(:, :, p_idx) = 0.5 * (B_transformed + B_transformed');
                end
                
            catch ME
                warning('WEC_File_IO:TransformFailed', ...
                        'Coordinate transformation failed: %s. Returning input matrices.', ME.message);
                A_cg = A_origin;
                B_cg = B_origin;
            end
        end
        
        %% ═════════════════════════════════════════════════════════════
        %%  TRANSFORMATION VERIFICATION
        %% ═════════════════════════════════════════════════════════════

        function verify_transformation_3x3(A_origin, A_cg, B_origin, B_cg, Z_CG)
            % VERIFY_TRANSFORMATION_3X3  Analytic spot-checks on the
            %   body-origin → CG coordinate transform.
            %
            %   WEC_File_IO.verify_transformation_3x3(A_o, A_cg, B_o, B_cg, Z_CG)
            %
            %   WHY verify analytically?
            %     The transform is a single matrix multiply that is easy
            %     to get wrong by sign or index.  Expanding H' M H
            %     element-by-element gives closed-form relations:
            %       A_cg(1,1) = A_o(1,1)                             [invariant]
            %       A_cg(2,2) = A_o(2,2)                             [invariant]
            %       A_cg(1,3) = A_o(1,3) − Z_CG · A_o(1,1)
            %       A_cg(3,3) = Z_CG² · A_o(1,1) − 2·Z_CG · A_o(1,3) + A_o(3,3)
            %     If any of these fail, the T matrix or the bilinear formula
            %     is wrong.  Printed to console for every WAMIT read so the
            %     user sees PASS/FAIL before the optimiser starts.
            %
            %   Additional checks:
            %     - Matrix symmetry |A − A'| < ε
            %     - Positive semi-definiteness  min(eig(A)) ≥ 0
            
            try
                fprintf('\n  ===============================================\n');
                fprintf('  |   TRANSFORMATION VERIFICATION (v4.0)       |\n');
                fprintf('  ===============================================\n\n');
                
                A_origin_inf = A_origin(:, :, end);
                A_cg_inf = A_cg(:, :, end);
                A_o = 0.5 * (A_origin_inf + A_origin_inf');
                
                fprintf('  Z_CG = %.4f m\n\n', Z_CG);
                
                % A(1,1): surge-surge [INVARIANT]
                expected_11 = A_o(1,1);
                actual_11 = A_cg_inf(1,1);
                error_11 = abs(actual_11 - expected_11);
                fprintf('  A(1,1) Surge-Surge (invariant):\n');
                fprintf('    Original:  %12.3f kg\n', A_o(1,1));
                fprintf('    Actual:    %12.3f kg\n', actual_11);
                fprintf('    Error:     %12.6f', error_11);
                WEC_File_IO.print_pass_fail(error_11, 1e-6);
                
                % A(2,2): heave-heave [INVARIANT]
                expected_22 = A_o(2,2);
                actual_22 = A_cg_inf(2,2);
                error_22 = abs(actual_22 - expected_22);
                fprintf('  A(2,2) Heave-Heave (invariant):\n');
                fprintf('    Original:  %12.3f kg\n', A_o(2,2));
                fprintf('    Actual:    %12.3f kg\n', actual_22);
                fprintf('    Error:     %12.6f', error_22);
                WEC_File_IO.print_pass_fail(error_22, 1e-6);
                
                % A(1,3): surge-pitch coupling
                expected_13 = A_o(1,3) - Z_CG * A_o(1,1);
                actual_13 = A_cg_inf(1,3);
                error_13 = abs(actual_13 - expected_13);
                fprintf('  A(1,3) Surge-Pitch coupling:\n');
                fprintf('    Formula: A_o(1,3) - Z_CG*A_o(1,1)\n');
                fprintf('    Original:  %12.3f kg*m\n', A_o(1,3));
                fprintf('    Expected:  %12.3f kg*m\n', expected_13);
                fprintf('    Actual:    %12.3f kg*m\n', actual_13);
                fprintf('    Error:     %12.6f', error_13);
                WEC_File_IO.print_pass_fail(error_13, 1e-3);
                
                % A(3,3): pitch-pitch
                expected_33 = Z_CG^2 * A_o(1,1) - 2*Z_CG * A_o(1,3) + A_o(3,3);
                actual_33 = A_cg_inf(3,3);
                error_33 = abs(actual_33 - expected_33);
                fprintf('  A(3,3) Pitch-Pitch:\n');
                fprintf('    Formula: Z_CG^2*A_o(1,1) - 2*Z_CG*A_o(1,3) + A_o(3,3)\n');
                fprintf('    Original:  %12.3f kg*m^2\n', A_o(3,3));
                fprintf('    Expected:  %12.3f kg*m^2\n', expected_33);
                fprintf('    Actual:    %12.3f kg*m^2\n', actual_33);
                fprintf('    Error:     %12.6f', error_33);
                WEC_File_IO.print_pass_fail(error_33, 1e-3);
                
                % Symmetry
                A_symmetry_error = max(max(abs(A_cg_inf - A_cg_inf')));
                fprintf('  Matrix Symmetry:\n');
                fprintf('    A_cg symmetry error: %.2e', A_symmetry_error);
                WEC_File_IO.print_pass_fail(A_symmetry_error, 1e-12);
                
                % Positive semi-definiteness
                A_eigenvals = eig(A_cg_inf);
                min_eigenval = min(A_eigenvals);
                fprintf('  Positive Semi-Definiteness:\n');
                fprintf('    Minimum eigenvalue: %12.3f', min_eigenval);
                if min_eigenval > -1e-6
                    fprintf('  PASS\n');
                else
                    fprintf('  FAIL\n');
                end
                
                fprintf('\n  ===============================================\n');
                all_pass = (error_11 < 1e-6) && (error_22 < 1e-6) && ...
                           (error_13 < 1e-3) && (error_33 < 1e-3) && ...
                           (A_symmetry_error < 1e-12) && (min_eigenval > -1e-6);
                
                if all_pass
                    fprintf('  VERIFICATION RESULT: ALL CHECKS PASSED\n');
                else
                    fprintf('  VERIFICATION RESULT: SOME CHECKS FAILED\n');
                end
                fprintf('  ===============================================\n\n');
                
            catch ME
                warning('WEC_File_IO:VerificationFailed', ...
                        'Verification failed: %s', ME.message);
            end
        end
        
    end
    
    %% ═════════════════════════════════════════════════════════════
    %%  PRIVATE FORMATTING HELPERS
    %% ═════════════════════════════════════════════════════════════

    methods (Static)  % Convention Bridge transforms (public — called by WEC_Visualization)
        function print_matrix_3x3(M, prefix)
            % Print a 3x3 matrix with formatting
            for ii = 1:3
                fprintf('%s[', prefix);
                for jj = 1:3
                    fprintf('%12.3f', M(ii,jj));
                    if jj < 3
                        fprintf(', ');
                    end
                end
                fprintf(']\n');
            end
        end
        
        function print_pass_fail(error, tol)
            % Print PASS/FAIL based on error tolerance
            if error < tol
                fprintf('  PASS\n\n');
            else
                fprintf('  FAIL (tol: %.0e)\n\n', tol);
            end
        end


        function [A_2D, B_2D, Fe_2D] = bem_to_2d_cg(A_3x3_O, B_3x3_O, z_G, Fe_3x1_O)
            % BEM_TO_2D_CG  Convention Bridge doc v1.1, Steps 2-3.
            %
            %   Single authoritative implementation.  WEC_Visualization and
            %   WEC_Main_Optimizer both delegate here.
            %
            %   INPUTS
            %     A_3x3_O  -- [3x3]  added mass at BEM origin, [surge,heave,pitch]
            %     B_3x3_O  -- [3x3]  radiation damping at BEM origin
            %     z_G      -- [m]    CG z-coordinate in BEM frame (signed; <0 below WL)
            %     Fe_3x1_O -- [3x1]  (optional) excitation force at origin
            %
            %   OUTPUTS
            %     A_2D, B_2D -- [3x3]  at CG, CCW-positive pitch (2D EOM indices 1,2,3)
            %     Fe_2D      -- [3x1]  at CG, CCW-positive pitch ([] if not supplied)
            %
            %   STEP 2 -- Congruence transform origin to CG:
            %     T    = [1  0  -z_G;  0  1  0;  0  0  1]
            %     A_CG = T' * A_O * T          (T, NOT T^-1)
            %
            %   STEP 3 -- Pitch sign change BEM to 2D EOM:
            %     S    = diag(1, 1, -1)
            %     A_2D = S * A_CG * S
            %
            %   Net effect on key entries:
            %     A11_2D = A11_O
            %     A22_2D = A22_O
            %     A33_2D = z_G^2*A11_O - 2*z_G*A13_O + A33_O
            %     A13_2D = z_G*A11_O - A13_O
            %     Fe3_2D = z_G*Fe1_O - Fe5_O

            T = [1, 0, -z_G; 0, 1, 0; 0, 0, 1];
            S = diag([1, 1, -1]);

            A_symm = 0.5 * (A_3x3_O + A_3x3_O');
            B_symm = 0.5 * (B_3x3_O + B_3x3_O');

            A_CG = T' * A_symm * T;
            B_CG = T' * B_symm * T;

            A_2D = S * A_CG * S;
            B_2D = S * B_CG * S;

            if nargin >= 4 && ~isempty(Fe_3x1_O)
                Fe_2D = S * (T' * Fe_3x1_O(:));
            else
                Fe_2D = [];
            end
        end

    end
end