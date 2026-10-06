function [A_cg, B_cg] = transform_to_cg(A_origin, B_origin, Z_CG)
%TRANSFORM_TO_CG  Shift 3x3 added-mass and damping matrices to the CG.
%   Inputs are [3x3xN] surge/heave/pitch matrices at the body origin and
%   Z_CG [m] (positive up). Uses q_origin=T*q_CG and M_CG=T'*M_origin*T;
%   each slice is symmetrised before and after transfer.
%   See docs/METHODS_ENGINE.md#bem-cg-congruence.
    try
        [n_i, n_j, n_periods] = size(A_origin);
        
        % Map CG velocities to body-origin velocities: q_origin = T * q_CG.
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
            warning('mwecmass:bem:InverseError', 'T * Tinv != I');
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
        warning('mwecmass:bem:TransformFailed', ...
                'Coordinate transformation failed: %s. Returning input matrices.', ME.message);
        A_cg = A_origin;
        B_cg = B_origin;
    end
end
