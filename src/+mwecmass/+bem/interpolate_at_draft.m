function [A11, A33, A55, A_full, B_full] = interpolate_at_draft(vertical_shift, config, z_cg_target)
%INTERPOLATE_AT_DRAFT  Interpolate cached added mass and damping by vertical shift.
% Inputs are vertical_shift [m], a config from build_config, and optional target CG [m].
% Outputs A11/A33/A55 and A_full are added mass at omega->infinity; B_full is band-averaged
% radiation damping.  Scalar and full matrices are 3x3; missing full fields fall back to
% diagonal/zero matrices.  Drafts outside the cached range are clamped.
% Matrices are stored at the HAMS-CG in config.hydro_z_cg; when z_cg_target is supplied,
% mwecmass.bem.transform_to_cg applies M_target = T' * M * T once.  Omitted targets retain
% HAMS-CG values.  See docs/METHODS_ENGINE.md#73-hydrodynamic-cache-interpolation-and-cg-transfer.  This free-floating-
% body interface has no PTO stiffness output.

    if nargin < 3, z_cg_target = []; end

    try
        if ~isfield(config, 'hydro_drafts') || ~isfield(config, 'added_mass_diagonal')
            error('Hydrodynamic coefficient table not found in config');
        end

        drafts = config.hydro_drafts;
        A_data = config.added_mass_diagonal;

        % Clamp to data range
        draft_clamped = max(min(vertical_shift, max(drafts)), min(drafts));

        % Handle single-draft case
        if length(drafts) == 1 %#ok<ISCL> -- scalar-cache comparison
            A11 = A_data(1, 1);
            A33 = A_data(1, 2);
            A55 = A_data(1, 3);

            if isfield(config, 'added_mass_full') && ~isempty(config.added_mass_full)
                A_full = config.added_mass_full{1};
            else
                A_full = diag([A11, A33, A55]);
            end

            if isfield(config, 'radiation_damping_full') && ~isempty(config.radiation_damping_full)
                B_full = config.radiation_damping_full{1};
            else
                B_full = zeros(3, 3);
            end
        else
            A11 = interp1(drafts, A_data(:,1), draft_clamped, 'linear', 'extrap');
            A33 = interp1(drafts, A_data(:,2), draft_clamped, 'linear', 'extrap');
            A55 = interp1(drafts, A_data(:,3), draft_clamped, 'linear', 'extrap');

            if isfield(config, 'added_mass_full') && ~isempty(config.added_mass_full)
                A_full = mwecmass.bem.interpolate_matrix(...
                    drafts, config.added_mass_full, draft_clamped);
            else
                A_full = diag([A11, A33, A55]);
            end

            if isfield(config, 'radiation_damping_full') && ~isempty(config.radiation_damping_full)
                B_full = mwecmass.bem.interpolate_matrix(...
                    drafts, config.radiation_damping_full, draft_clamped);
            else
                B_full = zeros(3, 3);
            end
        end

        % Ensure non-negative diagonal terms
        A11 = max(A11, 0);
        A33 = max(A33, 0);
        A55 = max(A55, 0);

        for ii = 1:3
            A_full(ii, ii) = max(A_full(ii, ii), 0);
            B_full(ii, ii) = max(B_full(ii, ii), 0);
        end

        % Optional delta transform: HAMS-CG → z_cg_target.
        %  Apply the centralized congruence transform.
        if ~isempty(z_cg_target) && ...
                isfield(config, 'hydro_z_cg') && ...
                ~isempty(config.hydro_z_cg)
            if length(drafts) == 1 %#ok<ISCL> -- scalar-cache comparison
                z_cg_hams = config.hydro_z_cg(1);
            else
                z_cg_hams = interp1(drafts, config.hydro_z_cg, ...
                                    draft_clamped, 'linear', 'extrap');
            end
            dz = z_cg_target - z_cg_hams;
            if abs(dz) > 1e-4
                % mwecmass.bem.transform_to_cg
                % applies M_target = T'·M·T where T encodes the
                % CG shift specified by its third argument.
                % Composes correctly: applying T(z_cg_hams) (during
                % rebuild_config_hydro) followed by T(dz) here is
                % algebraically identical to T(z_cg_target) applied
                % to the origin-frame matrices, so the net result
                % is M at z_cg_target irrespective of the HAMS-CG
                % intermediate.
                [A_full, B_full] = ...
                    mwecmass.bem.transform_to_cg( ...
                        A_full, B_full, dz);
                % Re-extract scalar coefficients after the transform.
                A11 = max(0, A_full(1, 1));
                A33 = max(0, A_full(2, 2));
                A55 = max(0, A_full(3, 3));
            end
        end

    catch ME
        warning('mwecmass:bem:CoefficientInterpolationFailed', ...
                'Hydrodynamic coefficient interpolation failed: %s. Returning zeros.', ME.message);
        A11 = 0; A33 = 0; A55 = 0;
        A_full = zeros(3, 3);
        B_full = zeros(3, 3);
    end
end
