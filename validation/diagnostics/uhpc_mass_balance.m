function uhpc_mass_balance(result_file, output_folder)
%UHPC_MASS_BALANCE Four-panel UHPC mass-balance figure (closure, distribution, build-up, ledger).
% Input: result_file (path to .mat holding results/final_props); output_folder optional
% (default Output/diagnostics/). Writes UHPC_Mass_Balance.png/.fig, or skips if result file
% is not modular_precast type. z positive upward; draft positive downward.
% Units: z, d [m]; A [m^2]; M [kg]; rho [kg/m^3].

    if nargin < 2 || isempty(output_folder)
        output_folder = mwecmass.output.output_dir('diagnostics');
    end
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end

    [results, ~] = mwecmass.output.load_results(result_file);
    if ~strcmp(results.realisation_type, 'modular_precast')
        fprintf('uhpc_mass_balance skipped: realisation_type = ''%s'' (modular_precast only).\n', ...
                results.realisation_type);
        return;
    end

    r   = results;
    cfg = r.config;
    ct  = r.constructability;

    rho_w    = cfg.RHO_WATER;             % [kg/m^3]
    rho_uhpc = ct.rho_uhpc;               % [kg/m^3] UHPC
    rho_air  = ct.rho_air;                % [kg/m^3]
    zmin     = cfg.hull_z_min;
    zmax     = cfg.hull_z_max;

    zg  = ct.z_grid(:);                   % [m] body frame
    Ao  = ct.A_outer_grid(:);             % [m^2]
    Ai  = ct.A_inner_grid(:);             % [m^2]
    Aj  = max(0, Ao - Ai);                % [m^2] jacket annulus
    zbal  = ct.z_ballast;                      % [m]
    vs  = ct.vertical_shift;              % [m]
    zwl = -vs;                            % [m] waterline, body frame
    se  = ct.strip_edges(:);
    Ns  = numel(se) - 1;
    wIdx = ct.wall_strip_idx;

    % EFFECTIVE material split, respecting the z_ballast rule that integrate_split
    % applies: below z_ballast the WHOLE section is UHPC, whatever A_inner says.
    % The strip-aware grid can carry A_inner > 0 below z_ballast, because a strip is
    % flagged solid only when its ENTIRE span lies below z_ballast.  Shading by
    % A_inner alone therefore draws void where the mass model counts material.
    below  = zg <= zbal;
    A_mat  = Aj;   A_mat(below) = Ao(below);   % [m^2] UHPC area
    A_vd   = Ai;   A_vd(below)  = 0;           % [m^2] void area

    fig = figure('Position', [50 50 1500 950], 'Color', 'w', ...
                 'Visible', 'off');
    tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf(['UHPC realisation - mass balance   ' ...
          '(\\rho_{UHPC}=%.0f, \\rho_{void}=%.2f kg/m^3, t=%.1f mm, N=%d strips)'], ...
          rho_uhpc, rho_air, ct.t_min*1000, Ns), ...
          'FontWeight', 'bold', 'FontSize', 13);

    C_uhpc = [0.35 0.42 0.55];
    C_void = [0.93 0.93 0.90];
    C_wl   = [0.10 0.45 0.80];

    %% ---- (a) CLOSURE: buoyancy curve vs realised mass ---------------------
    ax = nexttile(tl, 1); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    % Buoyancy as a function of draft.  draft = -(zmin + vs) and z_wl = -vs,
    % so sweeping z_wl over the hull height sweeps the draft.
    zwl_s  = linspace(zmin, zmax, 400)';
    Vsub_s = max(0, interp1(cfg.Aw_table_z, cfg.V_sub_table, zwl_s, 'linear', 0));
    d_s    = zwl_s - zmin;                          % [m] draft, positive
    Mbuoy_s = rho_w * Vsub_s;                       % [kg]

    plot(ax, d_s, Mbuoy_s/1e3, '-', 'Color', C_wl, 'LineWidth', 2.0, ...
         'DisplayName', '\rho_w V_{sub}(d)   (buoyancy)');
    yline(ax, ct.M_total/1e3, '--', 'Color', [0.75 0.25 0.15], 'LineWidth', 1.8, ...
          'Label', sprintf('M_{total} = %.0f kg', ct.M_total), ...
          'LabelHorizontalAlignment', 'left', 'FontSize', 9, ...
          'DisplayName', 'M_{total} (realised)');
    d_star = ct.draft;
    plot(ax, d_star, ct.M_total/1e3, 'o', 'MarkerSize', 9, 'LineWidth', 1.8, ...
         'MarkerFaceColor', 'w', 'Color', [0.75 0.25 0.15], ...
         'DisplayName', sprintf('equilibrium: d* = %.3f m', d_star));

    xlabel(ax, 'draft  d  [m]'); ylabel(ax, 'mass  [tonne]');
    title(ax, '(a)  Closure: mass = displaced mass', 'FontSize', 11);
    legend(ax, 'Location', 'northwest', 'FontSize', 9);
    xlim(ax, [0, zmax - zmin]);
    text(ax, 0.98, 0.06, sprintf('residual |ceq| = %.1e', ...
         abs(ct.mass_balance_error_pct/100)), 'Units', 'normalized', ...
         'HorizontalAlignment', 'right', 'FontSize', 9, 'Color', [0.3 0.3 0.3]);

    %% ---- (b) DISTRIBUTION: where the material sits ------------------------
    ax = nexttile(tl, 2); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    fill(ax, [Ao; flipud(A_vd)], [zg; flipud(zg)], C_uhpc, ...
         'EdgeColor', 'none', 'FaceAlpha', 0.85, 'DisplayName', 'UHPC');
    fill(ax, [A_vd; zeros(size(A_vd))], [zg; flipud(zg)], C_void, ...
         'EdgeColor', [0.6 0.6 0.6], 'LineWidth', 0.4, 'DisplayName', 'void');
    plot(ax, Ao, zg, 'k-', 'LineWidth', 1.4, 'DisplayName', 'A_{outer}(z)');
    plot(ax, Ai, zg, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.1, ...
         'DisplayName', 'A_{inner}(z)  (geometry only)');

    yline(ax, zbal,  '-',  'Color', [0.85 0.35 0.10], 'LineWidth', 1.8, ...
          'Label', sprintf('z_{ballast} = %.3f m', zbal), 'FontSize', 9, ...
          'LabelHorizontalAlignment', 'right', 'DisplayName', 'z_{ballast}');
    yline(ax, zwl, '--', 'Color', C_wl, 'LineWidth', 1.8, ...
          'Label', sprintf('waterline = %.3f m', zwl), 'FontSize', 9, ...
          'LabelHorizontalAlignment', 'left', 'DisplayName', 'waterline');
    for k = 2:Ns
        yline(ax, se(k), ':', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.7, ...
              'HandleVisibility', 'off');
    end

    xlabel(ax, 'sectional area  [m^2]'); ylabel(ax, 'z  (body frame)  [m]');
    title(ax, '(b)  Material distribution: solid below z_{ballast}, hollow above', ...
          'FontSize', 11);
    legend(ax, 'Location', 'east', 'FontSize', 7.5);
    ylim(ax, [zmin, zmax]);
    text(ax, 0.97, 0.055, ...
         sprintf('note: A_{inner}>0 between %.3f and %.3f m,\nbut that band is below z_{ballast} \\Rightarrow counted solid', ...
                 se(2), zbal), 'Units','normalized','HorizontalAlignment','right', ...
         'FontSize',7.5,'Color',[0.35 0.35 0.35]);

    %% ---- (c) BUILD-UP: cumulative mass vs cumulative displacement ---------
    ax = nexttile(tl, 3); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    % A_mat / A_vd already carry the z_ballast rule (computed once, above).
    dM_dz  = rho_uhpc*A_mat + rho_air*A_vd;                     % [kg/m]
    M_cum  = cumtrapz(zg, dM_dz);                                % [kg]

    Vsub_c = max(0, interp1(cfg.Aw_table_z, cfg.V_sub_table, zg, 'linear', 0));
    Mb_cum = rho_w * Vsub_c;                                     % [kg]

    plot(ax, M_cum/1e3, zg, '-', 'Color', C_uhpc, 'LineWidth', 2.2, ...
         'DisplayName', 'cumulative hull mass  M(z)');
    plot(ax, Mb_cum/1e3, zg, '-', 'Color', C_wl, 'LineWidth', 1.6, ...
         'DisplayName', 'cumulative displaced mass  \rho_w V_{sub}(z)');
    yline(ax, zwl, '--', 'Color', C_wl, 'LineWidth', 1.5, 'HandleVisibility','off');
    yline(ax, zbal,  '-',  'Color', [0.85 0.35 0.10], 'LineWidth', 1.5, 'HandleVisibility','off');
    plot(ax, ct.M_total/1e3, zwl, 'o', 'MarkerSize', 9, 'LineWidth', 1.8, ...
         'MarkerFaceColor', 'w', 'Color', [0.75 0.25 0.15], ...
         'DisplayName', 'balance point at the waterline');
    yline(ax, ct.CG_z_body, ':', 'Color', [0.2 0.6 0.2], 'LineWidth', 1.6, ...
          'Label', sprintf('CG_z = %.3f m', ct.CG_z_body), 'FontSize', 9, ...
          'HandleVisibility','off');
    xlim(ax, [0, 1.05*max(ct.M_total, max(Mb_cum))/1e3]);

    xlabel(ax, 'cumulative mass from keel  [tonne]');
    ylabel(ax, 'z  (body frame)  [m]');
    title(ax, '(c)  Build-up from the keel', 'FontSize', 11);
    legend(ax, 'Location', 'southeast', 'FontSize', 8);
    ylim(ax, [zmin, zmax]);

    %% ---- (d) LEDGER: per-strip mass split + rho comparison ----------------
    ax = nexttile(tl, 4); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

    rho_target = r.stage2_3d.x_optimal(:);  rho_target = rho_target(2:end);  % [kg/m^3]
    M_u = zeros(Ns,1);  M_v = zeros(Ns,1);  V_o = zeros(Ns,1);

    for k = 1:Ns
        zlo = se(k);  zhi = se(k+1);
        % Include z_ballast as a breakpoint ONLY when it lies inside this strip.
        % Adding it unconditionally extends the integration past the strip edge
        % and double-counts volume (this is the bug that made sum(V) = 2x V_hull).
        bp = [zlo; zhi; zg(zg > zlo & zg < zhi)];
        if zbal > zlo && zbal < zhi, bp = [bp; zbal]; end %#ok<AGROW> -- bp is rebuilt from [zlo; zhi; ...] fresh each k iteration; this conditionally appends at most one breakpoint, not an accumulating loop.
        bp = unique(sort(bp));
        Ao_b = interp1(zg, Ao, bp, 'linear', 0);
        Ai_b = interp1(zg, Ai, bp, 'linear', 0);
        Aj_b = max(0, Ao_b - Ai_b);
        V_o(k) = trapz(bp, Ao_b);
        lo = bp <= zbal;  hi = bp >= zbal;
        Vu = 0;  Vv = 0;
        if sum(lo) >= 2, Vu = Vu + trapz(bp(lo), Ao_b(lo)); end
        if sum(hi) >= 2
            Vu = Vu + trapz(bp(hi), Aj_b(hi));
            Vv = Vv + trapz(bp(hi), Ai_b(hi));
        end
        M_u(k) = rho_uhpc * Vu;
        M_v(k) = rho_air * Vv;
    end
    rho_real = (M_u + M_v) ./ max(V_o, eps);   % [kg/m^3] realised effective

    hb = bar(ax, 1:Ns, [M_u, M_v]/1e3, 0.62, 'stacked');
    hb(1).FaceColor = C_uhpc;  hb(1).DisplayName = 'UHPC mass';
    hb(2).FaceColor = [0.85 0.75 0.45];  hb(2).DisplayName = 'void mass';
    xlabel(ax, 'strip index  (1 = keel)'); ylabel(ax, 'strip mass  [tonne]');
    title(ax, '(d)  Per-strip ledger: \Sigma = M_{total}', 'FontSize', 11);

    yyaxis(ax, 'right');
    plot(ax, 1:Ns, rho_target, 's--', 'MarkerSize', 8, 'LineWidth', 1.5, ...
         'Color', [0.75 0.25 0.15], 'MarkerFaceColor', 'w', ...
         'DisplayName', '\rho_i  optimiser target');
    plot(ax, 1:Ns, rho_real, 'd-', 'MarkerSize', 8, 'LineWidth', 1.8, ...
         'Color', [0.15 0.35 0.15], 'MarkerFaceColor', [0.6 0.85 0.6], ...
         'DisplayName', '\rho_{eff}  as realised');
    ylabel(ax, '\rho  [kg/m^3]');
    ax.YAxis(2).Color = [0.15 0.35 0.15];
    yyaxis(ax, 'left');
    xticks(ax, 1:Ns);
    xtl = arrayfun(@(k) sprintf('%d', k), 1:Ns, 'UniformOutput', false);
    xtl{wIdx} = sprintf('%d\n(wall)', wIdx);
    xticklabels(ax, xtl);
    legend(ax, 'Location', 'northoutside', 'Orientation', 'horizontal', 'FontSize', 7.5);
    text(ax, 0.02, 0.955, sprintf('void mass \\leq %.1f kg per strip - not resolvable at this scale', ...
         max(M_v)), 'Units','normalized','FontSize',7.5,'Color',[0.4 0.4 0.4]);

    txt = sprintf(['M_{UHPC} = %.1f kg   M_{void} = %.1f kg   M_{total} = %.1f kg\n' ...
                   'V_{UHPC} = %.3f m^3 (%.1f%% of envelope)   V_{void} = %.3f m^3\n' ...
                   '\\rho_w V_{sub} = %.1f kg   -> residual %.1e %%'], ...
        ct.M_uhpc, ct.M_air, ct.M_total, ...
        ct.V_uhpc, 100*ct.V_uhpc/ct.V_hull, ct.V_air, ...
        rho_w*ct.V_sub, ct.mass_balance_error_pct);
    annotation(fig, 'textbox', [0.52 0.005 0.46 0.055], 'String', txt, ...
        'FontSize', 8.5, 'EdgeColor', [0.8 0.8 0.8], 'BackgroundColor', [0.98 0.98 0.98], ...
        'VerticalAlignment', 'middle');

    png = fullfile(output_folder, 'UHPC_Mass_Balance.png');
    exportgraphics(fig, png, 'Resolution', 200);
    savefig(fig, fullfile(output_folder, 'UHPC_Mass_Balance.fig'));
    close(fig);   % 16 GB machine: never leave figures open

    fprintf('--- per-strip ledger ---\n');
    fprintf('strip   V_out[m3]  M_UHPC[kg]  M_void[kg]  rho_eff  rho_target\n');
    for k = 1:Ns
        fprintf('%4d %10.4f %11.1f %11.2f %9.1f %11.1f\n', ...
                k, V_o(k), M_u(k), M_v(k), rho_real(k), rho_target(k));
    end
    fprintf('sum  %10.4f %11.1f %11.2f\n', sum(V_o), sum(M_u), sum(M_v));
    fprintf('CHECK sum(M_UHPC)+sum(M_void) = %.4f kg  vs M_total = %.4f kg  (delta %.3e)\n', ...
            sum(M_u)+sum(M_v), ct.M_total, sum(M_u)+sum(M_v)-ct.M_total);
    fprintf('CHECK rho_w*V_sub = %.4f kg  vs M_total = %.4f kg  (delta %.3e)\n', ...
            rho_w*ct.V_sub, ct.M_total, rho_w*ct.V_sub - ct.M_total);
    fprintf('Written: %s\n', png);
end
