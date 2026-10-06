function M_interp = interpolate_matrix(drafts, M_cell, draft_target)
%INTERPOLATE_MATRIX  Element-wise linear interpolation of 3x3 matrices stored in a cell array.
    % drafts and M_cell define source cases; draft_target is extrapolated linearly.
    % Elements are stacked in row-major order, interpolated in one [N×9] call,
    % then reshaped/transposed to recover the [3×3] matrix.

    try
        n_cases = length(drafts);

        % Stack all matrix elements: rows = drafts, cols = (i,j)
        elem_mat = zeros(n_cases, 9);
        col = 1;
        for ii = 1:3
            for jj = 1:3
                for k = 1:n_cases
                    elem_mat(k, col) = M_cell{k}(ii, jj);
                end
                col = col + 1;
            end
        end

        % Single vectorised interp1 call across all 9 elements
        row_interp = interp1(drafts, elem_mat, draft_target, 'linear', 'extrap');

        % Reshape and symmetrize
        M_interp = reshape(row_interp, 3, 3)';
        M_interp = 0.5 * (M_interp + M_interp');

    catch
        M_interp = zeros(3, 3);
    end
end
