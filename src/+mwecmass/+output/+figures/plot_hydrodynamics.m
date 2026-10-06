function plot_hydrodynamics(hydro_cache, vertical_shift, final_props, config)
%PLOT_HYDRODYNAMICS Plot frequency-dependent added mass A(ω), radiation damping B(ω), RAO, and A_inf linearity check from HAMS cache.
% Inputs: hydro_cache (added_mass_omega, radiation_damping_omega, exciting_force_omega, added_mass_inf, omega, period_band, drafts); vertical_shift [m]; final_props (mass_total, periods, CG, CB).

    style = mwecmass.output.figures.presentation_style(config);
    try

    dof6     = [1 3 5];
    dof_num  = {'1','3','5'};
    dof_name = {'Surge','Heave','Pitch'};

    % ── Find closest cached draft ─────────────────────────────
    [~, idx] = min(abs(hydro_cache.drafts - vertical_shift));
    fprintf('  Plotting hydrodynamics at vs = %+.4f m (cache idx %d)\n', ...
            hydro_cache.drafts(idx), idx);

    A_full = hydro_cache.added_mass_omega{idx};
    B_full = hydro_cache.radiation_damping_omega{idx};
    omega  = hydro_cache.omega;

    if isempty(A_full) || isempty(omega)
        warning('mwecmass:figures:NoHydroData', ...
                'No frequency-dependent data at vs=%+.4f.', vertical_shift);
        return;
    end

    T  = 2*pi ./ omega;
    Nf = length(omega);
    A_3x3 = A_full(dof6, dof6, :);
    B_3x3 = B_full(dof6, dof6, :);
    A_inf  = hydro_cache.added_mass_inf{idx}(dof6, dof6);
    % Guard against period_band missing in older cache files.
    if isfield(hydro_cache, 'period_band') && ~isempty(hydro_cache.period_band)
        T_band = hydro_cache.period_band;
    else
        T_band = [4, 25];   % safe default: 4-25 s operational band
        warning('mwecmass:figures:NoTBand', ...
            'hydro_cache.period_band missing — using default [4, 25] s.');
    end

    % ══════════════════════════════════════════════════════════
    %  FIGURE 1: A(w) and B(w) vs Period
    % ══════════════════════════════════════════════════════════
    fig1 = mwecmass.output.figures.new_figure(style, 'double_column');
    set(fig1, 'Name', 'Hydrodynamic Coefficients');
    layout1 = tiledlayout(fig1, 2, 3);

    A_units = {'[kg]','[kg]','[kg$\cdot$m$^2$]'};
    B_units = {'[kg/s]','[kg/s]','[kg$\cdot$m$^2$/s]'};

    for d = 1:3
        ax = nexttile(layout1, d);  hold(ax, 'on');
        Av = squeeze(A_3x3(d,d,:));
        Bv = squeeze(B_3x3(d,d,:));
        Ai = A_inf(d,d);

        yyaxis(ax, 'left');
        hA = plot(ax, T, Av, '-', 'Color', style.color.dof(d,:));
        mwecmass.output.figures.style_line(hA, style, 'curve');
        hAi = yline(ax, Ai, '--', 'Color', style.color.reference, ...
              'HandleVisibility', 'off');
        mwecmass.output.figures.style_line(hAi, style, 'reference');
        ylabel(ax, sprintf('$A_{%s%s}$ %s', dof_num{d}, dof_num{d}, A_units{d}));
        set(ax, 'YColor', style.color.dof(d,:));

        yyaxis(ax, 'right');
        hB = plot(ax, T, Bv, '-', 'Color', style.color.dof_secondary(d,:));
        mwecmass.output.figures.style_line(hB, style, 'curve');
        ylabel(ax, sprintf('$B_{%s%s}$ %s', dof_num{d}, dof_num{d}, B_units{d}));
        set(ax, 'YColor', style.color.dof_secondary(d,:));

        shade_T_band(ax, T_band, style);
        xlabel(ax, '$T$ [s]');
        title(ax, sprintf('%s -- $A_\\infty$ = %.0f', dof_name{d}, Ai));
        lg = legend(ax, {'$A(\omega)$','$B(\omega)$'}, 'Location', 'best');
        mwecmass.output.figures.style_legend(lg, style);
        % apply_axes_style is called on each yyaxis side in turn, after every text object on the
        % axes exists: a yyaxis axes exposes only the active side's YLabel (verified against
        % MathWorks doc "Modify Properties of Charts with Two y-Axes"), and the axes-level
        % properties (font, grid, title, xlabel) this helper also sets must not be applied before
        % title/xlabel/legend are created, or a later title()/xlabel() call can reset them.
        yyaxis(ax, 'left');  mwecmass.output.figures.apply_axes_style(ax, style);
        yyaxis(ax, 'right'); mwecmass.output.figures.apply_axes_style(ax, style);
    end

    pairs = [1 2; 1 3; 2 3];
    for k = 1:3
        ax = nexttile(layout1, 3+k);  hold(ax, 'on');
        ii = pairs(k,1);  jj = pairs(k,2);
        Av = squeeze(A_3x3(ii,jj,:));
        Bv = squeeze(B_3x3(ii,jj,:));
        Ai = A_inf(ii,jj);

        yyaxis(ax, 'left');
        hA = plot(ax, T, Av, '-', 'Color', style.color.dof_pair(k,:));
        mwecmass.output.figures.style_line(hA, style, 'curve');
        hAi = yline(ax, Ai, '--', 'Color', style.color.reference, ...
              'HandleVisibility', 'off');
        mwecmass.output.figures.style_line(hAi, style, 'reference');
        ylabel(ax, sprintf('$A_{%s%s}$', dof_num{ii}, dof_num{jj}));
        set(ax, 'YColor', style.color.dof_pair(k,:));

        yyaxis(ax, 'right');
        hB = plot(ax, T, Bv, '-', 'Color', style.color.dof_pair_secondary(k,:));
        mwecmass.output.figures.style_line(hB, style, 'curve');
        ylabel(ax, sprintf('$B_{%s%s}$', dof_num{ii}, dof_num{jj}));
        set(ax, 'YColor', style.color.dof_pair_secondary(k,:));

        shade_T_band(ax, T_band, style);
        xlabel(ax, '$T$ [s]');
        title(ax, sprintf('%s--%s -- $A_\\infty$ = %.1f', ...
              dof_name{ii}, dof_name{jj}, Ai));
        lg = legend(ax, {'$A(\omega)$','$B(\omega)$'}, 'Location', 'best');
        mwecmass.output.figures.style_legend(lg, style);
        yyaxis(ax, 'left');  mwecmass.output.figures.apply_axes_style(ax, style);
        yyaxis(ax, 'right'); mwecmass.output.figures.apply_axes_style(ax, style);
    end

    title(layout1, sprintf( ...
        'Hydrodynamic Coefficients at Origin ($v_s$ = %+.3f m)', vertical_shift));
    mwecmass.output.figures.apply_layout_style(layout1, style);

    saved1 = mwecmass.output.figures.export_figure(fig1, 'WEC_HydroCoeffs', config);
    fprintf('  Figure saved: %s\n', saved1{1});

    % ══════════════════════════════════════════════════════════
    %  FIGURE 2: RAO
    % ══════════════════════════════════════════════════════════
    has_Fe = isfield(hydro_cache, 'exciting_force_omega') && length(hydro_cache.exciting_force_omega) >= idx ...
             && ~isempty(hydro_cache.exciting_force_omega{idx});

    fig2 = mwecmass.output.figures.new_figure(style, 'double_column');
    set(fig2, 'Name', 'Response Amplitude Operators');
    layout2 = tiledlayout(fig2, 2, 3);

    mass   = final_props.mass_total;
    Iyy_cg = final_props.Iyy;
    z_G    = final_props.CG_total(3);
    rho_w  = config.RHO_WATER;
    g_acc  = config.G;

    Iyy_O = Iyy_cg + mass * z_G^2;

    M_O = [mass,       0,    mass*z_G;
           0,          mass, 0;
           mass*z_G,   0,    Iyy_O   ];

    assert(Iyy_O >= Iyy_cg - 1e-6, ...
        'Steiner produced Iyy_O < Iyy_CG — sign bug in z_G');

    C33 = rho_w * g_acc * final_props.Aw;
    z_B = final_props.CB(3);
    C55 = rho_w*g_acc*final_props.I_wp_yy ...
        + rho_w*g_acc*final_props.V_sub*z_B ...
        - mass*g_acc*z_G;
    C55 = max(C55, 0);

    C55_check = mass * g_acc * final_props.GM_L;
    if abs(C55 - C55_check) / max(abs(C55), 1) > 0.01
        fprintf('    WARNING: C55 formula (%.1f) vs mg·GM (%.1f) differ by %.1f%%\n', ...
                C55, C55_check, abs(C55-C55_check)/max(abs(C55),1)*100);
    end

    C_O = diag([0, C33, C55]);

    if isfield(config, 'B_visc_diag') && ~isempty(config.B_visc_diag)
        B_visc = diag(config.B_visc_diag);
    else
        B_visc = zeros(3);
    end
    has_visc = any(diag(B_visc) > 0);

    if has_Fe
        Fe_3xM = hydro_cache.exciting_force_omega{idx}(dof6, :);
    end

    RAO_complex = zeros(3, Nf);
    for m_idx = 1:Nf
        w = omega(m_idx);
        A_o = 0.5*(A_3x3(:,:,m_idx) + A_3x3(:,:,m_idx)');
        B_o = 0.5*(B_3x3(:,:,m_idx) + B_3x3(:,:,m_idx)');
        Z = -w^2*(M_O + A_o) + 1i*w*(B_o + B_visc) + C_O;

        if has_Fe
            if abs(det(Z)) > 1e-30
                X_O = Z \ Fe_3xM(:, m_idx);
                RAO_complex(:, m_idx) = [X_O(1) + z_G * X_O(3);
                                         X_O(2);
                                         X_O(3)];
            end
        else
            if abs(det(Z)) > 1e-30
                H = Z \ eye(3);
                RAO_complex(1, m_idx) = H(1,1) + z_G * H(3,1);
                RAO_complex(2, m_idx) = H(2,2);
                RAO_complex(3, m_idx) = H(3,3);
            end
        end
    end

    RAO   = abs(RAO_complex);
    PHASE = angle(RAO_complex) * (180/pi);

    if has_Fe
        RAO(3,:) = RAO(3,:) * (180/pi);
    end

    if has_Fe
        rao_titles = {'Surge RAO at CG $|X_1^{CG}/A|$', ...
                      'Heave RAO $|X_3/A|$', ...
                      'Pitch RAO $|X_5/A|$'};
        rao_ylabels = {'$|X_1^{CG}/A|$ [m/m]', ...
                       '$|X_3/A|$ [m/m]', ...
                       '$|X_5/A|$ [deg/m]'};
        phase_titles = {'Surge phase $\angle X_1^{CG}$', ...
                        'Heave phase $\angle X_3$', ...
                        'Pitch phase $\angle X_5$'};
    else
        rao_titles = {'Surge $|H_{11}^{CG}|$ (no $F_e$)', ...
                      'Heave $|H_{33}|$ (no $F_e$)', ...
                      'Pitch $|H_{55}|$ (no $F_e$)'};
        rao_ylabels = {'$|H_{11}^{CG}|$ [m/N]', ...
                       '$|H_{33}|$ [m/N]', ...
                       '$|H_{55}|$ [rad/(N$\cdot$m)]'};
        phase_titles = {'Surge phase $\angle H_{11}^{CG}$', ...
                        'Heave phase $\angle H_{33}$', ...
                        'Pitch phase $\angle H_{55}$'};
    end

    T_natural = [final_props.periods.surge, final_props.periods.heave, final_props.periods.pitch];
    T_nat_lbl = {'$T_{surge}$', '$T_{heave}$', '$T_{pitch}$'};

    for d = 1:3
        ax = nexttile(layout2, d);  hold(ax, 'on');
        h = plot(ax, T, RAO(d,:), '-', 'Color', style.color.dof(d,:));
        mwecmass.output.figures.style_line(h, style, 'curve');
        if isfinite(T_natural(d)) && T_natural(d) > 0 && T_natural(d) < max(T)*1.5
            hT = xline(ax, T_natural(d), ':', 'Color', style.color.reference, ...
                  'Label', sprintf('%s = %.1f s', T_nat_lbl{d}, T_natural(d)), ...
                  'LabelVerticalAlignment', 'top');
            mwecmass.output.figures.style_line(hT, style, 'reference');
            mwecmass.output.figures.style_text(hT, style, 'label');
        end
        xlabel(ax, '$T$ [s]');
        ylabel(ax, rao_ylabels{d});
        title(ax, rao_titles{d});
        mwecmass.output.figures.apply_axes_style(ax, style);
    end

    for d = 1:3
        ax = nexttile(layout2, d + 3);  hold(ax, 'on');
        h = plot(ax, T, PHASE(d,:), '-', 'Color', style.color.dof(d,:));
        mwecmass.output.figures.style_line(h, style, 'curve');
        if isfinite(T_natural(d)) && T_natural(d) > 0 && T_natural(d) < max(T)*1.5
            hT = xline(ax, T_natural(d), ':', 'Color', style.color.reference);
            mwecmass.output.figures.style_line(hT, style, 'reference');
        end
        h0  = yline(ax, 0, '-', 'Color', style.color.reference);
        hm  = yline(ax, -90, '--', 'Color', style.color.reference);
        hp  = yline(ax, 90, '--', 'Color', style.color.reference);
        mwecmass.output.figures.style_line(h0, style, 'reference');
        mwecmass.output.figures.style_line(hm, style, 'reference');
        mwecmass.output.figures.style_line(hp, style, 'reference');
        ylim(ax, [-180 180]);
        set(ax, 'YTick', [-180 -90 0 90 180]);
        xlabel(ax, '$T$ [s]');
        ylabel(ax, 'Phase [deg]');
        title(ax, phase_titles{d});
        mwecmass.output.figures.apply_axes_style(ax, style);
    end

    if has_visc
        visc_str = sprintf('$B_{visc}$ = [%.0f, %.0f, %.0f]', ...
            config.B_visc_diag(1), config.B_visc_diag(2), config.B_visc_diag(3));
    else
        visc_str = 'Potential flow only ($B_{visc} = 0$)';
    end

    if has_Fe, sgt = 'RAO (from HAMS excitation force)';
    else,      sgt = 'Frequency Response (no $F_e$ data)';
    end
    title(layout2, sprintf('%s -- $v_s$ = %+.3f m, $m$ = %.0f kg -- %s', ...
            sgt, vertical_shift, mass, visc_str));
    mwecmass.output.figures.apply_layout_style(layout2, style);

    saved2 = mwecmass.output.figures.export_figure(fig2, 'WEC_RAO', config);
    fprintf('  Figure saved: %s\n', saved2{1});

    % ══════════════════════════════════════════════════════════
    %  FIGURE 3: A_inf diagonals vs vertical_shift
    % ══════════════════════════════════════════════════════════
    N_cache = length(hydro_cache.drafts);
    if N_cache < 2
        fprintf('  Skipping A_inf vs draft plot (need >= 2 cache entries)\n');
        return;
    end

    drafts_all = hydro_cache.drafts(:);
    A_inf_diag = zeros(N_cache, 3);
    valid = true(N_cache, 1);
    for kk = 1:N_cache
        Ak = hydro_cache.added_mass_inf{kk};
        if isempty(Ak) || max(abs(Ak(:))) < 1e-6
            valid(kk) = false; continue;
        end
        A_inf_diag(kk,:) = [Ak(1,1), Ak(3,3), Ak(5,5)];
    end

    drafts_v = drafts_all(valid);
    A_vals_v = A_inf_diag(valid,:);
    if size(A_vals_v,1) < 2; return; end

    fig3 = mwecmass.output.figures.new_figure(style, 'double_column');
    set(fig3, 'Name', 'A_inf vs Draft');
    layout3 = tiledlayout(fig3, 1, 3);
    [drafts_s, si] = sort(drafts_v);
    A_s = A_vals_v(si,:);

    A_labels = {'$A_{11}^{\infty}$ [kg]', '$A_{33}^{\infty}$ [kg]', ...
                '$A_{55}^{\infty}$ [kg$\cdot$m$^2$]'};
    for d = 1:3
        ax = nexttile(layout3, d);  hold(ax, 'on');
        % style_line's marker branch applies style.marker.size (Marker 'o' is not 'none').
        h = plot(ax, drafts_s, A_s(:,d), 'o', 'Color', style.color.dof(d,:), ...
             'MarkerFaceColor', style.color.dof(d,:));
        mwecmass.output.figures.style_line(h, style, 'curve');

        p = polyfit(drafts_s, A_s(:,d), 1);
        xf = linspace(min(drafts_s), max(drafts_s), 100);
        % The fit line reuses the curve's own DOF colour (role dof) rather than an arithmetic
        % darkened shade of it; the two are distinguished by line style ('--' vs '-') instead.
        hfit = plot(ax, xf, polyval(p, xf), '--', 'Color', style.color.dof(d,:));
        mwecmass.output.figures.style_line(hfit, style, 'reference');

        A_pred = polyval(p, drafts_s);
        SS_res = sum((A_s(:,d) - A_pred).^2);
        SS_tot = sum((A_s(:,d) - mean(A_s(:,d))).^2);
        if SS_tot > 1e-12; R2 = 1 - SS_res/SS_tot; else; R2 = 1; end

        hx = xline(ax, vertical_shift, ':', 'Color', style.color.reference, 'Label', '$v_s^*$');
        mwecmass.output.figures.style_line(hx, style, 'reference');
        mwecmass.output.figures.style_text(hx, style, 'label');

        xlabel(ax, 'Vertical Shift [m]');
        ylabel(ax, A_labels{d});
        title(ax, sprintf('%s -- $R^2$ = %.4f (slope = %.1f/m)', ...
              dof_name{d}, R2, p(1)));
        lg = legend(ax, {'HAMS', 'Linear fit'}, 'Location', 'best');
        mwecmass.output.figures.style_legend(lg, style);
        mwecmass.output.figures.apply_axes_style(ax, style);
    end

    title(layout3, sprintf( ...
        '$A^{\\infty}$ vs Vertical Shift (%d points -- linearity check)', ...
        length(drafts_v)));
    mwecmass.output.figures.apply_layout_style(layout3, style);

    saved3 = mwecmass.output.figures.export_figure(fig3, 'WEC_Ainf_vs_Draft', config);
    fprintf('  Figure saved: %s\n', saved3{1});

    catch ME
        warning('mwecmass:figures:HydroPlotFailed', ...
                'Hydrodynamic plot failed: %s', ME.message);
    end
end

function shade_T_band(ax, T_band, style)
%SHADE_T_BAND Shade the operational-period band region on a plot.
    yl = ylim(ax);
    fill(ax, [T_band(1) T_band(2) T_band(2) T_band(1)], ...
         [yl(1) yl(1) yl(2) yl(2)], ...
         style.color.band, 'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
end
