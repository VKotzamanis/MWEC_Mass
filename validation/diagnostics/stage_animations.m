function stage_animations(result_file, output_folder, which, frames)
%STAGE_ANIMATIONS One animation per optimisation procedure in the mass-distribution pipeline.
% Inputs: result_file (path to .mat with results/final_props); output_folder (default
% Output/diagnostics/); which (default 'all'; one of '1'|'2'|'2b'|'2c'|'3'); frames (default [],
% all iterates; or subset for GIF 2c, or scalar for one 300-dpi PNG preview). Builds five GIFs
% (stages 1, 2, 2b, 2c, 3) from result, applies presentation style, writes animated output.

    if nargin < 2 || isempty(output_folder)
        output_folder = mwecmass.output.output_dir('diagnostics');
    end
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end
    if nargin < 3 || isempty(which), which = 'all'; end
    if nargin < 4, frames = []; end

    [r, ~] = mwecmass.output.load_results(result_file);
    cfg = reconstruct_live_config(r.config);   % restores ms2_model, boundary_cache.data (GIF 3)
    st  = mwecmass.output.figures.presentation_style(cfg);
    st.frames  = frames(:)';   % [] = every iterate; otherwise a subset (GIF 2c)
    st.is_demo = contains(upper(result_file), 'DEMO');
    st.suffix  = '';
    if st.is_demo
        st.suffix = '_DEMO';
        fprintf('DEMO source: %s  (outputs suffixed _DEMO and stamped)\n', result_file);
    end
    fprintf('source: %s   stage1_mode = %s\n', result_file, cfg.stage1_mode);

    outdir = output_folder;

    % Ordered closed CCW silhouette, body frame (PUBLIC method)
    prof = mwecmass.output.figures.build_silhouette_profile(cfg, 300);
    fprintf('silhouette: %d vertices, x in [%.3f, %.3f], z in [%.3f, %.3f]\n', ...
        size(prof,1), min(prof(:,1)), max(prof(:,1)), min(prof(:,2)), max(prof(:,2)));

    sel = @(t) strcmp(which,'all') || strcmp(which,t);

    if sel('1'),  gif1_stage1(r, cfg, prof, st, outdir);        end
    if sel('2'),  gif2_stage2(r, cfg, prof, st, outdir);        end
    if sel('2b'), gif2b_stage2_uhpc(r, cfg, prof, st, outdir);  end
    if sel('2c'), gif2c_stage2_pair(r, cfg, prof, st, outdir);  end
    if sel('3'),  gif3_stage3(r, cfg, prof, st, outdir);        end
    fprintf('\nGIFS DONE (%s)\n', which);
end


%% =====================================================================
function gif1_stage1(r, cfg, prof, st, outdir)
% STAGE 1 - draft screening: one frame per cached draft.
    if ~isfield(r.stage1_2d,'convergence_data') || ...
            ~isfield(r.stage1_2d.convergence_data,'sweep')
        fprintf('GIF 1 skipped: no sweep data.\n'); return
    end
    sw = r.stage1_2d.convergence_data.sweep;
    N  = numel(sw.vs);
    fprintf('\n=== GIF 1  Stage 1 draft screening (%d drafts) ===\n', N);

    GM = arrayfun(@(k) sw.props{k}.GM_L, 1:N);
    Th = arrayfun(@(k) sw.props{k}.periods.heave, 1:N);
    Tp = arrayfun(@(k) sw.props{k}.periods.pitch, 1:N);
    Th(~isfinite(Th)) = NaN;  Tp(~isfinite(Tp)) = NaN;
    fv   = sw.fval(:)';
    feas = logical(sw.feasible);
    [~, best] = min(fv + (~feas)*1e12);

    % Tier-2 shortlist: ranked on |T_heave - goal| (see run_stage1_sweep)
    Thr = Th;  Thr(isnan(Thr)) = Inf;
    [~, rk] = sort(abs(Thr - cfg.T_heave_goal));
    topK = rk(1:min(cfg.n_sweep_refine, N));

    % Fixed colour scale across all frames so strip colours do not jump
    allr = cell2mat(cellfun(@(x) x(2:end)', sw.x, 'UniformOutput', false)');
    clim_rho = [min(allr(:)), max(allr(:))];
    cmap = mwecmass.output.cividis_map(256);

    gif = fullfile(outdir,'Stage1_Draft_Screening.gif');
    if exist(gif,'file'), delete(gif); end

    for k = 1:N
        vs  = sw.vs(k);
        rho = sw.x{k}(:);  rho = rho(2:end);

        fig = new_fig(st);
        tl  = tiledlayout(fig,3,2);
        ttl = sprintf('Stage 1: Draft Screening  --  draft %d of %d,  $v_s$ = %+.3f m,  %s', ...
              k, N, vs, mwecmass.internal.ternary(feas(k),'feasible','infeasible'));
        if k == best, ttl = [ttl '  $\Leftarrow$ selected warm start']; end %#ok<AGROW> -- ttl is rebuilt via sprintf fresh each k iteration; this appends a fixed suffix once, not an accumulating loop.
        title(tl, ttl);
        mwecmass.output.figures.apply_layout_style(tl, st);

        ax = nexttile(tl,1,[3 1]);
        draw_density_section(ax, prof, vs, cfg, rho, cmap, clim_rho, st);
        title(ax, sprintf('$d$ = %.3f m,  $V_{sub}$ = %.2f m$^3$,  $GM$ = %+.3f m', ...
              abs(cfg.hull_z_min+vs), sw.props{k}.V_sub, GM(k)));

        % objective (log: the penalty guard is 4 decades above the feasible pair)
        ax = nexttile(tl,2); hold(ax,'on');
        plot(ax, sw.vs, fv,'o-','Color',[0.45 0.45 0.45],'LineWidth',st.line_width.boundary, ...
             'MarkerSize',5,'DisplayName','objective f');
        set(ax,'YScale','log');
        plot(ax, sw.vs(topK), fv(topK),'s','MarkerSize',11,'LineWidth',1.4, ...
             'Color',st.fill_palette.waterline,'DisplayName','Tier-2 shortlist');
        plot(ax, vs, fv(k),'o','MarkerSize',12,'LineWidth',2.0, ...
             'Color',st.fill_palette.ballast_level,'DisplayName','current');
        ylabel(ax,'$f$ [--]');
        lg = legend(ax,'Location','southeast');
        mwecmass.output.figures.style_legend(lg, st);
        title(ax,sprintf('$f$ = %.3f   (guard = %g)',fv(k),cfg.penalty_guard));
        mwecmass.output.figures.apply_axes_style(ax, st);

        % GM
        ax = nexttile(tl,4); hold(ax,'on');
        yl = [min([GM cfg.gm_min])-0.15, max([GM cfg.gm_min])+0.15];
        patch(ax,[min(sw.vs) max(sw.vs) max(sw.vs) min(sw.vs)], ...
                 [yl(1) yl(1) cfg.gm_min cfg.gm_min], st.fill_palette.wall_boundary, ...
                 'FaceAlpha',0.10,'EdgeColor','none','DisplayName','GM < GM_{min}');
        plot(ax, sw.vs, GM,'o-','Color',[0.15 0.45 0.20],'LineWidth',st.line_width.main, ...
             'MarkerSize',5,'DisplayName','GM');
        plot(ax, vs, GM(k),'o','MarkerSize',12,'LineWidth',2.0,'Color',st.fill_palette.ballast_level, ...
             'HandleVisibility','off');
        ylim(ax,yl);
        ylabel(ax,'$GM$ [m]');
        lg = legend(ax,'Location','northwest');
        mwecmass.output.figures.style_legend(lg, st);
        mwecmass.output.figures.apply_axes_style(ax, st);

        % periods
        ax = nexttile(tl,6); hold(ax,'on');
        xb = [min(sw.vs) max(sw.vs)];
        band(ax, xb, cfg.T_heave_range, st.fill_palette.waterline, 'T_2 band');
        band(ax, xb, cfg.T_pitch_range, st.fill_palette.ballast_level,  'T_3 band');
        plot(ax, sw.vs, Th,'o-','Color',st.fill_palette.waterline,'LineWidth',st.line_width.main, ...
             'MarkerSize',5,'DisplayName','T_2 heave');
        plot(ax, sw.vs, Tp,'d-','Color',st.fill_palette.ballast_level,'LineWidth',st.line_width.main, ...
             'MarkerSize',5,'DisplayName','T_3 pitch (NaN = \infty)');
        plot(ax, vs, Th(k),'o','MarkerSize',12,'LineWidth',2.0,'Color',st.fill_palette.ballast_level, ...
             'HandleVisibility','off');
        xlabel(ax,'$v_s$ [m]');
        ylabel(ax,'$T$ [s]');
        lg = legend(ax,'Location','northeast');
        mwecmass.output.figures.style_legend(lg, st);
        mwecmass.output.figures.apply_axes_style(ax, st);

        dt = 1.0;  if k == best, dt = 3.0; end
        append_frame(fig, gif, k==1, dt);
        fprintf('  frame %2d/%d  vs=%+.3f  f=%9.3f  GM=%+.3f  feas=%d\n', ...
                k,N,vs,fv(k),GM(k),feas(k));
    end
    fprintf('  -> %s\n', gif);
end


%% =====================================================================
function gif2_stage2(r, cfg, prof, st, outdir)
% STAGE 2 - mass-distribution optimisation: one frame per logged SQP iterate.
    if ~isfield(r.stage2_3d,'trajectory') || isempty(r.stage2_3d.trajectory)
        fprintf('\nGIF 2 skipped: results.stage2_3d.trajectory absent.\n'); return
    end
    t = r.stage2_3d.trajectory;  n = t.n;
    fprintf('\n=== GIF 2  Stage 2 mass distribution (%d iterates) ===\n', n);

    se = cfg.strip_edges(:);  Ns = numel(se)-1;
    X  = t.x;   rho_all = X(2:end,:);
    clim_rho = [min(rho_all(:)), max(rho_all(:))];
    cmap = mwecmass.output.cividis_map(256);

    Th = t.T_heave; Th(~isfinite(Th)) = NaN;
    Tp = t.T_pitch; Tp(~isfinite(Tp)) = NaN;
    Ml = pad([t.mass_total t.mass_buoy], 0.10);
    Gl = pad([t.GM cfg.gm_min], 0.25);
    Tl = pad([Th Tp cfg.T_heave_goal cfg.T_pitch_goal], 0.12);

    gif = fullfile(outdir,'Stage2_Mass_Distribution.gif');
    if exist(gif,'file'), delete(gif); end

    for k = 1:n
        vs = t.vs(k);  rho = X(2:end,k);

        fig = new_fig(st);
        tl  = tiledlayout(fig,3,3);
        title(tl, sprintf(['Stage 2: Mass-Distribution Optimisation  --  ' ...
              'logged iterate %d of %d,  $f$ = %.4f'], k, n, t.fval(k)));
        mwecmass.output.figures.apply_layout_style(tl, st);

        ax = nexttile(tl,1,[3 1]);
        draw_density_section(ax, prof, vs, cfg, rho, cmap, clim_rho, st);
        title(ax, sprintf('$d$ = %.3f m,  $M$ = %.0f kg', t.draft(k), t.mass_total(k)));

        % rho(z) staircase with the constructability floor
        ax = nexttile(tl,2,[3 1]); hold(ax,'on');
        for i = 1:Ns
            patch(ax,[0 rho(i) rho(i) 0],[se(i) se(i) se(i+1) se(i+1)], ...
                  rho_rgb(rho(i),clim_rho,cmap),'EdgeColor',st.fill_palette.boundary, ...
                  'LineWidth',0.7,'FaceAlpha',0.85,'HandleVisibility','off');
            h_slbl = text(ax, rho(i)+0.02*clim_rho(2), 0.5*(se(i)+se(i+1)), ...
                 sprintf('S%d  %.0f', i, rho(i)), 'VerticalAlignment','middle');
            mwecmass.output.figures.style_text(h_slbl, st, 'annotation');
        end
        if ~isempty(cfg.per_strip_density_lb)
            lb = cfg.per_strip_density_lb(:);
            zz = reshape([se(1:end-1), se(2:end)]',[],1);
            rr = reshape([lb, lb]',[],1);
            plot(ax, rr, zz,'--','Color',st.fill_palette.wall_boundary,'LineWidth',1.4, ...
                 'DisplayName','\rho_{min} constructability');
            lg = legend(ax,'Location','southeast');
            mwecmass.output.figures.style_legend(lg, st);
        end
        xlim(ax,[0 1.30*clim_rho(2)]); ylim(ax,[min(se) max(se)]);
        xlabel(ax,'$\rho_i$ [kg/m$^3$]');
        ylabel(ax,'$z$ [m]  (body frame)');
        title(ax,'Ballast density per strip');
        mwecmass.output.figures.apply_axes_style(ax, st);

        trace_panel(nexttile(tl,3), 0:n-1, k, {t.mass_total, t.mass_buoy}, ...
            {'M_{total}','\rho_w V_{sub}'}, {'o-','s--'}, ...
            {[0.35 0.42 0.55], st.fill_palette.waterline}, Ml, '$M$ [kg]', ...
            sprintf('Mass balance:  %+.2e', t.mass_total(k)/max(t.mass_buoy(k),eps)-1), ...
            [], [], st);

        trace_panel(nexttile(tl,6), 0:n-1, k, {t.GM}, {'GM'}, {'o-'}, ...
            {[0.15 0.45 0.20]}, Gl, '$GM$ [m]', ...
            sprintf('$GM$ = %.4f m', t.GM(k)), cfg.gm_min, st.fill_palette.wall_boundary, st);

        ax = nexttile(tl,9); hold(ax,'on');
        band(ax,[-0.5 n-0.5], cfg.T_heave_range, st.fill_palette.waterline,'T_2 band');
        band(ax,[-0.5 n-0.5], cfg.T_pitch_range, st.fill_palette.ballast_level, 'T_3 band');
        plot(ax,0:k-1,Th(1:k),'o-','Color',st.fill_palette.waterline,'LineWidth',st.line_width.main, ...
             'MarkerSize',4,'DisplayName','T_2');
        plot(ax,0:k-1,Tp(1:k),'d-','Color',st.fill_palette.ballast_level,'LineWidth',st.line_width.main, ...
             'MarkerSize',4,'DisplayName','T_3');
        xlim(ax,[-0.5 n-0.5]); ylim(ax,Tl);
        xlabel(ax,'SQP iterate');
        ylabel(ax,'$T$ [s]');
        lg = legend(ax,'Location','southeast','NumColumns',2);
        mwecmass.output.figures.style_legend(lg, st);
        title(ax,sprintf('$T_2$ = %.3f s,  $T_3$ = %.3f s   (goals %.2f / %.2f)', ...
              Th(k),Tp(k),cfg.T_heave_goal,cfg.T_pitch_goal));
        mwecmass.output.figures.apply_axes_style(ax, st);

        dt = 0.55;  if k==n, dt = 3.0; end
        append_frame(fig, gif, k==1, dt);
        if k==1 || k==n || mod(k,5)==0
            fprintf('  frame %2d/%d  vs=%+.4f  M=%8.1f  GM=%.4f  T2=%.3f  f=%.4f\n', ...
                    k,n,vs,t.mass_total(k),t.GM(k),Th(k),t.fval(k));
        end
    end
    fprintf('  -> %s\n', gif);
end


%% =====================================================================
function gif2b_stage2_uhpc(r, cfg, prof, st, outdir)
% STAGE 2, RENDERED AS A UHPC PARTITION.
%
% Every Stage-2 iterate carries a per-strip effective density rho_i, but no
% UHPC geometry -- the realisation has not run yet.  Mass equivalence against
% a UHPC + air partition converts rho_i into a drawable volume fraction with
% no extra solve:
%
%     rho_i*V_i = rho_UHPC*V_u + rho_air*(V_i - V_u)
%   =>  f_i = V_u/V_i = (rho_i - rho_air) / (rho_UHPC - rho_air)
%
% This is the same relation behind the modular-precast realisation.scale_contour,
% whose s = sqrt((rho_UHPC - rho_eff)/(rho_UHPC - rho_air)) has s^2 equal to
% the void fraction.
%
% CAVEAT, stated on the figure: these frames show the UHPC distribution
% IMPLIED by each iterate's densities, not a solved realisation.  Only the
% final Stage-3 solve fixes an actual (t_i, z_ballast).
%
% Units: V [m^3], rho [kg/m^3], z/x/t [m].

    if ~isfield(r.stage2_3d,'trajectory') || isempty(r.stage2_3d.trajectory)
        fprintf('\nGIF 2b skipped: no Stage-2 trajectory.\n'); return
    end
    t = r.stage2_3d.trajectory;  n = t.n;
    fprintf('\n=== GIF 2b  Stage 2 as UHPC distribution (%d iterates) ===\n', n);

    se    = cfg.strip_edges(:);  Ns = numel(se)-1;
    Vs    = cfg.strip_V(:);                      % canonical Aw-trapz strip volumes
    rho_u = cfg.constructability_rho_hull;
    rho_v = cfg.constructability_rho_air;
    t_min = cfg.constructability_t_min;
    wIdx  = cfg.wall_strip_index;
    lb    = cfg.per_strip_density_lb(:);         % constructability floor per strip

    % t_min offset polygon, computed ONCE in body frame (a per-frame vs shift
    % is a pure z-translation, so the offset never needs recomputing).
    P_in_body = inner_offset(prof(:,1), prof(:,2), t_min);

    X = t.x;  RHO = X(2:end,:);
    F = min(1, max(0, (RHO - rho_v) ./ (rho_u - rho_v)));    % UHPC volume fraction
    Vu = F .* Vs;   Vv = max(0, Vs - Vu);
    Vmax = max(Vs);

    Th = t.T_heave; Th(~isfinite(Th)) = NaN;
    Tp = t.T_pitch; Tp(~isfinite(Tp)) = NaN;
    Ml = pad([t.mass_total t.mass_buoy], 0.10);
    Gl = pad([t.GM cfg.gm_min], 0.25);
    Tl = pad([Th Tp cfg.T_heave_goal cfg.T_pitch_goal], 0.12);

    gif = fullfile(outdir,'Stage2_UHPC_Distribution.gif');
    if exist(gif,'file'), delete(gif); end

    for k = 1:n
        vs = t.vs(k);

        fig = new_fig(st);
        tl  = tiledlayout(fig,3,3);
        title(tl, sprintf(['Stage 2: Implied UHPC Distribution  --  ' ...
              'logged iterate %d of %d,  $f$ = %.4f'], k, n, t.fval(k)));
        mwecmass.output.figures.apply_layout_style(tl, st);

        ax = nexttile(tl,1,[3 1]);
        draw_uhpc_from_fraction(ax, prof, P_in_body, vs, cfg, F(:,k), ...
                                Vu(:,k), Vv(:,k), lb, RHO(:,k), t_min, wIdx, st);
        title(ax, sprintf(['$d$ = %.3f m,  $V_{UHPC}$ = %.3f m$^3$,  ' ...
              '$V_{void}$ = %.3f m$^3$'], t.draft(k), sum(Vu(:,k)), sum(Vv(:,k))));

        % per-strip UHPC / void volume, stacked horizontally
        ax = nexttile(tl,2,[3 1]); hold(ax,'on');
        for i = 1:Ns
            zc = 0.5*(se(i)+se(i+1));  hh = 0.80*(se(i+1)-se(i));
            patch(ax,[0 Vu(i,k) Vu(i,k) 0],[zc-hh/2 zc-hh/2 zc+hh/2 zc+hh/2], ...
                  st.fill_palette.jacket_material,'EdgeColor',st.fill_palette.boundary,'LineWidth',0.7, ...
                  'HandleVisibility',mwecmass.internal.ternary(i==1,'on','off'),'DisplayName','V_{UHPC}');
            patch(ax,[Vu(i,k) Vs(i) Vs(i) Vu(i,k)],[zc-hh/2 zc-hh/2 zc+hh/2 zc+hh/2], ...
                  st.fill_palette.void,'EdgeColor',st.fill_palette.inner_boundary,'LineStyle','--','LineWidth',0.9, ...
                  'HandleVisibility',mwecmass.internal.ternary(i==1,'on','off'),'DisplayName','V_{void}');
            % A literal '%' character (from the sprintf '%%' escape) starts a LaTeX comment under
            % the math interpreter style_text now forces, so it is escaped to '\%' first.
            lbl_pct = strrep(sprintf('S%d  %.0f%%', i, 100*F(i,k)), '%', '\%');
            h_ulbl = text(ax, Vs(i)+0.03*Vmax, zc, lbl_pct, 'VerticalAlignment','middle');
            mwecmass.output.figures.style_text(h_ulbl, st, 'annotation');
        end
        xlim(ax,[0 1.32*Vmax]); ylim(ax,[min(se) max(se)]);
        xlabel(ax,'strip volume [m$^3$]');
        ylabel(ax,'$z$ [m]  (body frame)');
        title(ax,'UHPC vs void per strip');
        lg = legend(ax,'Location','southeast');
        mwecmass.output.figures.style_legend(lg, st);
        mwecmass.output.figures.apply_axes_style(ax, st);

        trace_panel(nexttile(tl,3), 0:n-1, k, {t.mass_total, t.mass_buoy}, ...
            {'M_{total}','\rho_w V_{sub}'}, {'o-','s--'}, ...
            {[0.35 0.42 0.55], st.fill_palette.waterline}, Ml, '$M$ [kg]', ...
            sprintf('Mass balance:  %+.2e', t.mass_total(k)/max(t.mass_buoy(k),eps)-1), ...
            [], [], st);

        trace_panel(nexttile(tl,6), 0:n-1, k, {t.GM}, {'GM'}, {'o-'}, ...
            {[0.15 0.45 0.20]}, Gl, '$GM$ [m]', ...
            sprintf('$GM$ = %.4f m', t.GM(k)), cfg.gm_min, st.fill_palette.wall_boundary, st);

        ax = nexttile(tl,9); hold(ax,'on');
        band(ax,[-0.5 n-0.5], cfg.T_heave_range, st.fill_palette.waterline,'T_2 band');
        band(ax,[-0.5 n-0.5], cfg.T_pitch_range, st.fill_palette.ballast_level, 'T_3 band');
        plot(ax,0:k-1,Th(1:k),'o-','Color',st.fill_palette.waterline,'LineWidth',st.line_width.main, ...
             'MarkerSize',4,'DisplayName','T_2');
        plot(ax,0:k-1,Tp(1:k),'d-','Color',st.fill_palette.ballast_level,'LineWidth',st.line_width.main, ...
             'MarkerSize',4,'DisplayName','T_3');
        xlim(ax,[-0.5 n-0.5]); ylim(ax,Tl);
        xlabel(ax,'SQP iterate');
        ylabel(ax,'$T$ [s]');
        lg = legend(ax,'Location','southeast','NumColumns',2);
        mwecmass.output.figures.style_legend(lg, st);
        title(ax,sprintf('$T_2$ = %.3f s,  $T_3$ = %.3f s   (goals %.2f / %.2f)', ...
              Th(k),Tp(k),cfg.T_heave_goal,cfg.T_pitch_goal));
        mwecmass.output.figures.apply_axes_style(ax, st);

        dt = 0.55;  if k==n, dt = 3.0; end
        append_frame(fig, gif, k==1, dt);
        if k==1 || k==n || mod(k,5)==0
            fprintf('  frame %2d/%d  f_UHPC = [', k, n);
            fprintf('%.2f ', F(:,k));  fprintf(']  V_UHPC=%.2f m^3\n', sum(Vu(:,k)));
        end
    end
    fprintf('  -> %s\n', gif);
end


function gif2c_stage2_pair(r, cfg, prof, st, outdir)
% STAGE 2, TWO PANELS ONLY -- the presentation figure.
%   (a) uniform per-segment ballast density (the optimiser's design variable)
%   (b) the UHPC / air partition it implies (the same state, read as material)
% The panels are butted together and share the vertical axis: (b) carries no
% Z label and no Z tick labels, so the pair reads as one section.
%
% CENTRE OF GRAVITY / CENTRE OF BUOYANCY  (panel a only)
%   Z_CG  trajectory.CG_z, i.e. props.CG_total(3) -- ALREADY WORLD FRAME
%         (mwecmass.optim.run.m:475, commented as such).
%   Z_CB  NOT stored per iterate.  Reconstructed from the precomputed hull
%         table that properties_3d.m:81 itself uses:
%             z_wl_body = -vs                       (stage_animations.m header)
%             CB_body   = interp1(cfg.Aw_table_z, cfg.CB_z_table, z_wl_body)
%             CB_world  = CB_body + vs
%   VERIFIED before use: CB_world + I_wp_yy/V_sub - CG_z reproduces
%         trajectory.GM to 0.0e+00 m on all 34 iterates -- the identity
%         GM_L = KM - CG_z with KM = CB_z + BM (properties_3d.m
%         :137,217).  A body-frame misreading of CG_z would show as a 1.075 m
%         residual, so the check discriminates.
%
% Units : Z, X, CB, CG, BM [m];  rho [kg/m^3];  T_2, T_3 [s].
% Sign  : world = body + vs, waterline at Z = 0, Z positive up.
%
% LAYOUT.  Axes Positions are computed from the true data aspect ratio so the
% drawn box EQUALS the Position rectangle.  Without this, `axis equal` shrinks
% the drawn box inside a wider Position and the southoutside colorbar -- which
% tracks Position, not the drawn box -- runs wider than the panel it labels.
% This is also why tiledlayout is not used here: it owns Position.

    if ~isfield(r.stage2_3d,'trajectory') || isempty(r.stage2_3d.trajectory)
        fprintf('\nGIF 2c skipped: no Stage-2 trajectory.\n'); return
    end
    t = r.stage2_3d.trajectory;  n = t.n;

    ks = 1:n;
    if isfield(st,'frames') && ~isempty(st.frames)
        ks = st.frames(ks_valid(st.frames, n));
    end
    preview = isscalar(ks);
    fprintf('\n=== GIF 2c  Stage 2 density + UHPC pair (%d of %d iterates) ===\n', ...
            numel(ks), n);

    % Presentation-only sizes for this one collided-panel figure: no out.style role
    % names a "large" variant of a font, so these stay local literals here, the same
    % way the small-multiple sizes of plot_modular_precast's per-strip panels do.
    FS_TITLE_PRES  = 18;   % pt, the 3-line figure title
    FS_LABEL_PRES  = 18;   % pt, panel titles and axis labels
    FS_AXIS_PRES   = 18;   % pt, tick numbers
    FS_ANNOT_PRES  = 15;   % pt, in-panel strip labels
    FS_LEGEND_PRES = 15;   % pt, the legend of panel (b)

    Vs    = cfg.strip_V(:);
    rho_u = cfg.constructability_rho_hull;
    rho_v = cfg.constructability_rho_air;
    t_min = cfg.constructability_t_min;
    wIdx  = cfg.wall_strip_index;
    lb    = cfg.per_strip_density_lb(:);

    P_in_body = inner_offset(prof(:,1), prof(:,2), t_min);

    RHO = t.x(2:end,:);
    F   = min(1, max(0, (RHO - rho_v) ./ (rho_u - rho_v)));
    Vu  = F .* Vs;   Vv = max(0, Vs - Vu);
    clim_rho = [min(RHO(:)), max(RHO(:))];      % fixed over the trajectory
    cmap = mwecmass.output.cividis_map(256);

    % --- CG and CB in the world frame, one value per iterate --------------
    vs_all    = t.vs(:);
    z_wl_body = -vs_all;                                            % [m]
    CB_world  = interp1(cfg.Aw_table_z, cfg.CB_z_table, z_wl_body, ...
                        'linear','extrap') + vs_all;                % [m]
    CG_world  = t.CG_z(:);                                          % [m]

    % --- fixed axis limits over the WHOLE trajectory ----------------------
    % Per-frame limits let the axes creep while the hull sits still; fixed
    % limits make the hull sink and rise against the waterline, which is the
    % point of the animation.
    px = prof(:,1);  pz_body = prof(:,2);
    xr = [min(px)-0.20, max(px)+0.20];
    zr = [min(pz_body)+min(vs_all)-0.15, max(pz_body)+max(vs_all)+0.25];
    xt = tick_vector(xr, 0.5);
    yt = tick_vector(zr, 0.5);

    % --- geometry of the two collided panels ------------------------------
    % Bands sized for 18 pt type: TOP holds the 3-line figure title plus the
    % panel titles; BOT holds the X label, then the colorbar (with its own
    % ticks and caption) under (a) and the two-row legend under (b).
    FIG_W = 1200;  FIG_H = 950;                    % [px]
    TOP   = mwecmass.internal.ternary(st.is_demo, 0.210, 0.175);     % title band
    BOT   = 0.200;                                 % X label + colorbar / legend
    GAP   = 0.012;                                 % the collision
    Hax   = 1 - TOP - BOT;
    Wax   = Hax * (FIG_H/FIG_W) * (diff(xr)/diff(zr));   % => drawn box == Position
    X0    = 0.5 - (2*Wax + GAP)/2 + 0.030;         % nudged right for the Z label
    pos_a = [X0,             BOT, Wax, Hax];
    pos_b = [X0 + Wax + GAP, BOT, Wax, Hax];

    gif    = fullfile(outdir, ['Stage2_Density_and_UHPC' st.suffix '.gif']);
    pngdir = fullfile(outdir, ['Stage2_Density_and_UHPC' st.suffix '_frames']);
    if ~preview
        if exist(gif,'file'), delete(gif); end
        if ~exist(pngdir,'dir'), mkdir(pngdir); end
    end

    first = true;
    for k = ks
        vs = t.vs(k);

        fig = figure('Color','w','Position',[40 40 FIG_W FIG_H], ...
                     'Visible','off','Renderer','painters', ...
                     'DefaultAxesFontName',st.font_name);

        % ---- figure title: iterate counter + the two natural periods ------
        ttl = {sprintf('Iteration %d of %d', k, n), ...
               sprintf('Natural Period in Heave, $T_2$ = %.3f s', t.T_heave(k)), ...
               sprintf('Natural Period in Pitch, $T_3$ = %.3f s', t.T_pitch(k))};
        if st.is_demo
            ttl = [{'DEMONSTRATION -- cold-start variant, not the production result'}, ttl]; %#ok<AGROW>
        end
        h_ttl2c = annotation(fig,'textbox', ...
            [X0, 1-TOP+0.050, 2*Wax+GAP, TOP-0.055], 'String', ttl, ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'EdgeColor','none','FitBoxToText','off');
        % FS_TITLE_PRES has no out.style role (18 pt, larger than font_size.title, for this one
        % collided-panel presentation figure -- see the comment above); style_text is still the
        % one place that reads the font name and math interpreter, with only the size locally
        % overridden (generalised for this figure's presentation scale).
        title2c_style = st;
        title2c_style.font_size.title = FS_TITLE_PRES;
        mwecmass.output.figures.style_text(h_ttl2c, title2c_style, 'title');

        % ---- (a) uniform segment density ---------------------------------
        oa = struct();
        oa.wl_style  = '-';                       % solid waterline
        oa.xlabel    = '$X$ (m)';
        oa.ylabel    = 'Distance from SWL, $Z$ (m)';
        oa.cb_label  = 'Segment''s Uniform Density, $\rho$ [kg/m$^3$]';
        oa.xlim = xr;  oa.ylim = zr;  oa.xtick = xt;  oa.ytick = yt;
        oa.fs_label  = FS_LABEL_PRES;
        oa.fs_axis   = FS_AXIS_PRES;
        oa.fs_annot  = FS_ANNOT_PRES;
        oa.fs_cb     = FS_AXIS_PRES;
        oa.fs_title  = FS_LABEL_PRES;
        % 90% of the panel width, centred: the bar stays inside the panel it
        % labels (it used to overrun it) and the end tick labels, which are
        % centred on the bar ends, no longer spill into panel (b)'s legend.
        oa.cb_pos    = [pos_a(1)+0.05*pos_a(3), 0.098, 0.90*pos_a(3), 0.024];
        oa.z_cg      = CG_world(k);
        oa.z_cb      = CB_world(k);

        ax_a = axes('Parent',fig,'Position',pos_a,'FontName',st.font_name);
        draw_density_section(ax_a, prof, vs, cfg, RHO(:,k), cmap, clim_rho, st, oa);
        % apply_axes_style, called inside draw_density_section above with oa.fs_title routed
        % into its call_style, already set this Title's font/size/math interpreter; setting the
        % string only (no style args) here does not disturb them.
        title(ax_a,'Uniform Segment Density');

        % ---- (b) implied UHPC distribution -------------------------------
        ob = struct();
        ob.wl_style     = '-';
        ob.xlabel       = '$X$ (m)';
        ob.show_ylabel  = false;                  % collided: no Z label
        ob.ytick_labels = false;                  % collided: no Z numbers
        ob.xlim = xr;  ob.ylim = zr;  ob.xtick = xt;  ob.ytick = yt;
        ob.fs_label     = FS_LABEL_PRES;
        ob.fs_axis      = FS_AXIS_PRES;
        ob.fs_annot     = FS_ANNOT_PRES;
        ob.fs_legend    = FS_LEGEND_PRES;
        ob.fs_title     = FS_LABEL_PRES;
        ob.void_word    = 'Air';                  % Void -> Air throughout
        ob.shell_name   = 'UHPC Thickness';
        ob.show_wall_bound  = false;              % red wall line removed
        ob.show_implied_note = false;             % "implied by rho_i" removed
        ob.show_flags   = false;                  % drawing-caveat tag overflows
        ob.wl_label_offset = 0.12;                % [m] keep labels off the SWL
        ob.legend_items = {'wl','wall','shell','air'};
        ob.legend_cols  = 2;
        ob.legend_pos   = [pos_b(1), 0.022, pos_b(3), 0.082];

        ax_b = axes('Parent',fig,'Position',pos_b,'FontName',st.font_name);
        draw_uhpc_from_fraction(ax_b, prof, P_in_body, vs, cfg, F(:,k), ...
                                Vu(:,k), Vv(:,k), lb, RHO(:,k), t_min, wIdx, st, ob);
        % apply_axes_style, called inside draw_uhpc_from_fraction above with ob.fs_title routed
        % into its call_style, already set this Title's font/size/math interpreter; setting the
        % string only (no style args) here does not disturb them.
        title(ax_b,'Distribution of UHPC');

        % ---- export ------------------------------------------------------
        if preview
            png = fullfile(outdir, sprintf('Stage2_Density_and_UHPC%s_iter%02d.png', ...
                                           st.suffix, k));
            im = print(fig,'-RGBImage','-r300');   % 300 dpi
            close(fig);
            imwrite(im, png);
            fprintf('  iterate %d: %d x %d px  Z_CG = %+.4f m  Z_CB = %+.4f m\n', ...
                    k, size(im,2), size(im,1), CG_world(k), CB_world(k));
            fprintf('  -> %s\n', png);
        else
            dt = 0.55;  if k==ks(end), dt = 3.0; end
            append_frame(fig, gif, first, dt, ...
                         fullfile(pngdir, sprintf('frame_%03d.png', k)));
            first = false;
            if k==ks(1) || k==ks(end) || mod(k,10)==0
                fprintf('  frame %2d/%d  V_UHPC=%.3f m^3  Z_CG=%+.4f  Z_CB=%+.4f\n', ...
                        k, n, sum(Vu(:,k)), CG_world(k), CB_world(k));
            end
        end
    end
    if ~preview
        fprintf('  -> %s\n', gif);
        fprintf('  -> %s  (300 dpi frames)\n', pngdir);
    end
end


function idx = ks_valid(f, n)
    idx = find(f >= 1 & f <= n);
    if isempty(idx), error('stage_animations:badFrame', ...
        'requested frame(s) outside 1..%d', n); end
end


function tv = tick_vector(lim, step)
% Deterministic ticks, so the two collided panels cannot disagree.
    tv = ceil(lim(1)/step)*step : step : floor(lim(2)/step)*step;
end


function draw_uhpc_from_fraction(ax, prof, P_in_body, vs, cfg, f, ~, ~, ...
                                 lb, rho_i, t_min, wIdx, st, opts)
% Draw the UHPC/void partition implied by per-strip volume fractions f.
%
% OPTS (optional; defaults reproduce the original figure, so GIF 2b is
% unaffected):
%   wl_style        waterline style ('--' default, '-' presentation)
%   void_word       'Void' default; 'Air' for the presentation figure
%   shell_name      legend name of the shell patch ('UHPC Shell' default)
%   show_wall_bound  draw the red wall-boundary line (true default)
%   show_implied_note  draw the "implied by rho_i" caveat (true default)
%   legend_items    subset/order of {'wall','ballast','shell','void'|'air','wl','wb'}
%   legend_cols / legend_pos / fs_legend
%   xlabel / ylabel / xlim / ylim / xtick / ytick / ytick_labels
%   fs_label / fs_axis / fs_annot
% Same visual grammar as mwecmass.output.figures.plot_modular_precast FIGURE 1:
%   STEP A  t_min offset of the FULL profile (precomputed, body frame)
%   STEP B  clip to strip, then scale HORIZONTALLY about the void's x-centroid
%           until its 2D area equals A_strip*(1-f_i).  z preserved.
% Scaling is clamped so a >= t_min horizontal rim always survives; when the
% clamp binds, the strip is annotated -- that iterate wants a thinner wall
% than t_min allows, i.e. it is not realisable as drawn.
    if nargin < 14 || ~isstruct(opts), opts = struct(); end
    wl_style   = get_or(opts,'wl_style','--');
    void_word  = get_or(opts,'void_word','Void');
    fs_label   = get_or(opts,'fs_label', st.font_size.axes); %#ok<NASGU> -- kept for the opts
    % contract this function documents (a caller may still set ob.fs_label independently); the
    % xlabel/ylabel calls that read it directly are gone (apply_axes_style's pre-existing
    % call_style.font_size.axes = fs_axis now styles them instead, see the comment there).
    fs_axis    = get_or(opts,'fs_axis',  st.font_size.axes);
    fs_annot   = get_or(opts,'fs_annot', st.font_size.annotation);
    fs_legend  = get_or(opts,'fs_legend', 9);
    % Defaults to the normal title size (GIF 2b's call, no opts), overridden to FS_LABEL_PRES by
    % GIF 2c's opts -- kept independent of fs_label so this ax's Title does not silently follow
    % the axis-label size (kept independent for this presentation-scale case).
    fs_title   = get_or(opts,'fs_title', st.font_size.title);
    % style_text's 'annotation' role reads font_size.annotation; fs_annot-scaled locally so the
    % in-panel strip labels below keep tracking this function's own fs_annot override.
    annot_style = st;
    annot_style.font_size.annotation = fs_annot;

    hold(ax,'on');  axis(ax,'equal');
    se = cfg.strip_edges(:);  Ns = numel(se)-1;
    px = prof(:,1);  pz = prof(:,2) + vs;

    for i = 1:Ns
        zlo = se(i)+vs;  zhi = se(i+1)+vs;
        out = clip_poly(px, pz, zlo, zhi);
        if size(out,1) < 3, continue; end

        solid = f(i) >= 0.999;
        if i == wIdx,     fc = st.fill_palette.solid_material;
        elseif solid,     fc = st.fill_palette.ballast_material;
        else,             fc = st.fill_palette.jacket_material;
        end
        patch(ax, out(:,1), out(:,2), fc,'EdgeColor','none','HandleVisibility','off');

        flag = '';
        if ~solid && ~isempty(P_in_body)
            inn = clip_poly(P_in_body(:,1), P_in_body(:,2)+vs, zlo, zhi);
            if size(inn,1) >= 3
                A_strip = polyarea(out(:,1), out(:,2));
                A_void  = polyarea(inn(:,1), inn(:,2));
                A_tgt   = A_strip * (1 - f(i));
                if A_void > 1e-9 && A_tgt > 1e-9
                    sh   = A_tgt / A_void;
                    xmax_out = max(abs(out(:,1)));
                    xmax_in  = max(abs(inn(:,1)));
                    sh_cap   = max(0, (xmax_out - t_min) / max(xmax_in, eps));
                    % Clamp binding is a DRAWING limit, not a feasibility
                    % verdict: horizontal scaling of the t_min offset cannot
                    % always reach the target void area even for a strip that
                    % is perfectly realisable (the real solver varies t per
                    % strip instead).  Label it as a drawing caveat.  The
                    % genuine feasibility signal is rho_i < rho_min, below.
                    if sh > sh_cap
                        sh = sh_cap;
                        % Drawing caveat only (see above), so it is suppressed
                        % on the presentation figure via show_flags=false.
                        if get_or(opts,'show_flags',true)
                            flag = sprintf('  (%s clamped)', lower(void_word));
                        end
                    end
                    cx = mean(inn(:,1));
                    inn(:,1) = (inn(:,1)-cx)*sh + cx;
                end
                patch(ax, inn(:,1), inn(:,2), st.fill_palette.void, ...
                      'EdgeColor',st.fill_palette.inner_boundary,'LineStyle','--','LineWidth',1.0, ...
                      'HandleVisibility','off');
                mwecmass.output.figures.draw_hatch_strips(ax, inn(:,1), inn(:,2), ...
                      st.hatch_spacing, st.fill_palette.hatch);
            end
        end
        if rho_i(i) < lb(i) - 1e-6, flag = [flag '  ($\rho<\rho_{min}$)']; end %#ok<AGROW> -- flag is reset to '' at the top of each i iteration (line 678); this appends a fixed suffix once, not an accumulating loop.

        if i == wIdx
            lbl = sprintf('S%d: 100%% UHPC (wall)', i);
        elseif solid
            lbl = sprintf('S%d: 100%% UHPC (solid ballast)', i);
        else
            lbl = sprintf('S%d: %.0f%% UHPC / %.0f%% %s', ...
                          i, 100*f(i), 100*(1-f(i)), void_word);
        end
        % A strip straddling the waterline puts its label on top of the
        % waterline; push it clear rather than let the label's white box
        % erase a line the figure is meant to show.
        z_lbl = 0.5*(zlo+zhi);
        wl_off = get_or(opts,'wl_label_offset',0);
        if wl_off > 0 && abs(z_lbl) < wl_off, z_lbl = z_lbl - 1.6*wl_off; end
        % lbl's literal '%' characters (from '100%%'/'%.0f%%') would start a LaTeX comment
        % under the math interpreter style_text now forces, so they are escaped first; flag's
        % own '$\rho<\rho_{min}$' is already valid LaTeX math and is left unescaped.
        h_striplbl = text(ax, 0, z_lbl, [strrep(lbl,'%','\%') flag], ...
             'HorizontalAlignment','center', 'BackgroundColor','w','Margin',0.5);
        mwecmass.output.figures.style_text(h_striplbl, annot_style, 'annotation');
    end

    h_outer = plot(ax, px, pz,'-','Color',st.fill_palette.boundary);
    mwecmass.output.figures.style_line(h_outer, st, 'boundary');
    xr = get_or(opts,'xlim',[min(px)-0.2, max(px)+0.2]);
    zr = get_or(opts,'ylim',[min(pz)-0.15, max(pz)+0.25]);
    h.wl = plot(ax, xr,[0 0], wl_style,'Color',st.fill_palette.waterline, ...
                'DisplayName','Waterline');
    mwecmass.output.figures.style_line(h.wl, st, 'reference');
    if get_or(opts,'show_wall_bound',true)
        h.wb = plot(ax, xr,[se(wIdx)+vs se(wIdx)+vs],'-','Color',st.fill_palette.wall_boundary, ...
                    'DisplayName','Wall boundary');
        mwecmass.output.figures.style_line(h.wb, st, 'boundary');
    end
    h.wall    = patch(ax,nan,nan,st.fill_palette.solid_material,'DisplayName','UHPC Solid Wall');
    h.ballast = patch(ax,nan,nan,st.fill_palette.ballast_material,'DisplayName','UHPC Solid Ballast');
    h.shell   = patch(ax,nan,nan,st.fill_palette.jacket_material, ...
                      'DisplayName',get_or(opts,'shell_name','UHPC Shell'));
    h.void    = patch(ax,nan,nan,st.fill_palette.void,'EdgeColor',st.fill_palette.inner_boundary, ...
                      'LineStyle','--','DisplayName', ...
                      mwecmass.internal.ternary(strcmpi(void_word,'Void'),'Void (Air)',void_word));
    h.air     = h.void;                     % alias, so 'air' selects the same handle

    xlim(ax,xr);  ylim(ax,zr);
    if ~isempty(get_or(opts,'xtick',[])), xticks(ax, opts.xtick); end
    if ~isempty(get_or(opts,'ytick',[])), yticks(ax, opts.ytick); end
    % apply_axes_style below (call_style.font_size.axes = fs_axis, pre-existing) now styles
    % these labels, in place of the direct font/size/math-interpreter literals this pair used to
    % carry; fs_axis equals fs_label at every call site today (both default to
    % st.font_size.axes, and GIF 2c's opts set both to FS_AXIS_PRES/FS_LABEL_PRES = 18), so this
    % is not a value change.
    xlabel(ax, get_or(opts,'xlabel','$x$ [m]'));
    % show_ylabel, not an empty ylabel string: get_or() falls back to the
    % default on '' , which silently restored the label on the collided panel.
    if get_or(opts,'show_ylabel',true)
        ylabel(ax, get_or(opts,'ylabel','$z$ [m]   (world frame, waterline = 0)'));
    end

    items = get_or(opts,'legend_items',{'wall','ballast','shell','void','wl','wb'});
    hv = gobjects(0);
    for q = 1:numel(items)
        if isfield(h, items{q}), hv(end+1) = h.(items{q}); end %#ok<AGROW>
    end
    lg = legend(ax, hv,'Location','southoutside', ...
                'NumColumns',get_or(opts,'legend_cols',3));
    if get_or(opts,'show_implied_note',true)
        h_note = text(ax, 0.02, 0.015, 'implied by $\rho_i$ -- not a solved realisation', ...
             'Units','normalized','Color',[0.35 0.35 0.35]);
        mwecmass.output.figures.style_text(h_note, st, 'annotation');
    end
    % apply_axes_style reads style.font_size.tick_label (not .axes), so both are overridden
    % here to keep the fs_axis/fs_legend override this dual-use function relies on (GIF 2c
    % passes 'presentation' sizes via opts; GIF 2b's default fs_axis/fs_legend equal the
    % current tick_label/legend roles, so that call is unaffected).
    call_style = st;
    call_style.font_size.tick_label = fs_axis;
    call_style.font_size.axes       = fs_axis;
    call_style.font_size.legend     = fs_legend;
    call_style.font_size.title      = fs_title;
    mwecmass.output.figures.apply_axes_style(ax, call_style);
    mwecmass.output.figures.style_legend(lg, call_style);
    if ~get_or(opts,'ytick_labels',true)
        ax.YTickLabel = [];      % collided pair: (b) borrows (a)'s Z numbers
    end
    lp = get_or(opts,'legend_pos',[]);
    if ~isempty(lp), lg.Position = lp; end
end


%% =====================================================================
function gif3_stage3(r, cfg, prof, st, outdir)
% STAGE 3 - UHPC realisation, rendered the way visualize FIGURE 1 does it.
    h = [];
    if isfield(r,'constructability') && isstruct(r.constructability) && ...
            isfield(r.constructability,'iter_history')
        h = r.constructability.iter_history;
    end
    if isempty(h) || isempty(h.iter)
        fprintf('\nGIF 3 skipped: no UHPC iteration history.\n'); return
    end

    se   = cfg.strip_edges(:);  Ns = numel(se)-1;
    wIdx = cfg.wall_strip_index;
    is_sol = false(Ns,1);  is_sol(wIdx) = true;
    n_zg = get_or(cfg,'uhpc_n_z_grid',cfg.steel_n_z_grid);
    slf  = get_or(cfg,'uhpc_max_slope_factor',cfg.steel_max_slope_factor);

    keep = true(1,numel(h.iter));
    for k=2:numel(h.iter)
        keep(k) = ~(abs(h.vs(k)-h.vs(k-1))<1e-12 && abs(h.z_ballast(k)-h.z_ballast(k-1))<1e-12);
    end
    idx = find(keep);  nf = numel(idx);
    fprintf('\n=== GIF 3  Stage 3 UHPC realisation (%d of %d records) ===\n', nf, numel(h.iter));

    G = cell(nf,1);
    for j=1:nf
        [g,nz] = mwecmass.realise.modular_precast.build_geometry_grid( ...
            cfg, se, h.t_strip(:,idx(j)), is_sol, n_zg, slf);
        if nz > 0.05*n_zg, error('stage_animations:Degenerate','frame %d',j); end
        G{j} = g;
    end

    Ml = pad([h.M_total(idx) h.M_buoy(idx)], 0.10);
    Gl = pad([h.GM(idx) cfg.gm_min], 0.25);
    Tl = pad([h.T_heave(idx) h.T_pitch(idx) cfg.T_heave_goal cfg.T_pitch_goal], 0.12);

    gif = fullfile(outdir,'Stage3_UHPC_Realisation.gif');
    if exist(gif,'file'), delete(gif); end

    for j=1:nf
        k = idx(j);  g = G{j};
        fig = new_fig(st);
        tl  = tiledlayout(fig,3,2);
        title(tl, sprintf(['Stage 3: UHPC Realisation  --  SQP iterate %d of %d,' ...
              '  $\\Phi$ = %.4f'], h.iter(k), max(h.iter), h.phi(k)));
        mwecmass.output.figures.apply_layout_style(tl, st);

        ax = nexttile(tl,1,[3 1]);
        draw_uhpc_section(ax, prof, h.vs(k), h.z_ballast(k), cfg, g, is_sol, ...
                          h.t_strip(:,k), h.V_uhpc(k), h.V_air(k), st);

        trace_panel(nexttile(tl,2), h.iter(idx), j, ...
            {h.M_total(idx), h.M_buoy(idx)}, {'M_{total}','\rho_w V_{sub}'}, ...
            {'o-','s--'}, {[0.35 0.42 0.55], st.fill_palette.waterline}, Ml, '$M$ [kg]', ...
            sprintf('Mass balance:  $c_{eq}$ = %+.2e', h.ceq(k)), [], [], st);

        trace_panel(nexttile(tl,4), h.iter(idx), j, {h.GM(idx)}, {'GM'}, {'o-'}, ...
            {[0.15 0.45 0.20]}, Gl, '$GM$ [m]', ...
            sprintf('$GM$ = %.4f m', h.GM(k)), cfg.gm_min, st.fill_palette.wall_boundary, st);

        ax = nexttile(tl,6); hold(ax,'on');
        xb = [-0.4 max(h.iter)+0.4];
        band(ax, xb, cfg.T_heave_range, st.fill_palette.waterline,'T_2 band');
        band(ax, xb, cfg.T_pitch_range, st.fill_palette.ballast_level, 'T_3 band');
        plot(ax,h.iter(idx(1:j)),h.T_heave(idx(1:j)),'o-','Color',st.fill_palette.waterline, ...
             'LineWidth',st.line_width.main,'MarkerSize',4.5,'DisplayName','T_2');
        plot(ax,h.iter(idx(1:j)),h.T_pitch(idx(1:j)),'d-','Color',st.fill_palette.ballast_level, ...
             'LineWidth',st.line_width.main,'MarkerSize',4.5,'DisplayName','T_3');
        xlim(ax,xb); ylim(ax,Tl);
        xlabel(ax,'SQP iterate');
        ylabel(ax,'$T$ [s]');
        lg = legend(ax,'Location','southeast','NumColumns',2);
        mwecmass.output.figures.style_legend(lg, st);
        title(ax,sprintf('$T_2$ = %.3f s,  $T_3$ = %.3f s   (goals %.2f / %.2f)', ...
              h.T_heave(k),h.T_pitch(k),cfg.T_heave_goal,cfg.T_pitch_goal));
        mwecmass.output.figures.apply_axes_style(ax, st);

        dt = 0.85;  if j==1, dt=1.2; elseif j==nf, dt=3.0; end
        append_frame(fig, gif, j==1, dt);
        fprintf('  frame %d/%d  iter=%d  vs=%+.4f  zbal=%+.4f  M=%.1f  T2=%.3f\n', ...
                j,nf,h.iter(k),h.vs(k),h.z_ballast(k),h.M_total(k),h.T_heave(k));
    end
    fprintf('  -> %s\n', gif);
end


%% =====================================================================
%%  SECTION RENDERERS
%% =====================================================================
function draw_density_section(ax, prof, vs, cfg, rho, cmap, clim_rho, st, opts)
% Mirrors mwecmass.output.figures.plot_equivalent_density_2d: smooth silhouette shifted
% to the world frame, strips filled by scalar CData through the colormap at
% FaceAlpha 0.85, strip boundary lines drawn from actual profile crossings,
% outer boundary and waterline on top.  axis equal, as in the class.
%
% OPTS (all optional; the defaults reproduce the original figure exactly, so
% GIFs 1 and 2 are unaffected):
%   wl_style  waterline line style ('--' default, '-' for the presentation)
%   xlabel / ylabel / cb_label   label strings (LaTeX)
%   xlim / ylim / xtick / ytick  fixed limits and ticks ([] = auto per frame)
%   fs_label / fs_axis / fs_annot  font sizes [pt]
%   cb_pos    explicit colorbar Position, so it cannot outrun the panel width
%   z_cg / z_cb  world-frame [m] centre of gravity / buoyancy; drawn as a red
%             and a blue bullet at X = 0 when finite.
    if nargin < 9 || ~isstruct(opts), opts = struct(); end
    wl_style = get_or(opts,'wl_style','--');
    fs_label = get_or(opts,'fs_label', st.font_size.axes);
    fs_axis  = get_or(opts,'fs_axis',  st.font_size.axes);
    fs_annot = get_or(opts,'fs_annot', st.font_size.annotation);
    % Defaults to the normal title size (GIFs 1/2, no opts); GIF 2c's opts override it to
    % FS_LABEL_PRES, independent of fs_label -- see draw_uhpc_from_fraction's own fs_title.
    fs_title = get_or(opts,'fs_title', st.font_size.title);

    % Centre-of-gravity / centre-of-buoyancy marker colours: out.style names no role for
    % these two markers (plot_steel_solve's own CG/CB markers are plain 'r'/'b' short-name
    % colours, unparameterised), so they stay local literals here, unchanged from before.
    CG_RED = [0.85 0.10 0.10];
    CB_BLUE = [0.00 0.20 0.75];

    hold(ax,'on');  axis(ax,'equal');
    px = prof(:,1);  pz = prof(:,2) + vs;
    bounds = cfg.strip_edges(:) + vs;

    colormap(ax, cmap);  clim(ax, clim_rho);
    for i = 1:numel(rho)
        slab = clip_poly(px, pz, bounds(i), bounds(i+1));
        if size(slab,1) < 3, continue; end
        patch(ax, slab(:,1), slab(:,2), rho(i), ...
              'EdgeColor','none','FaceAlpha',0.85,'HandleVisibility','off');
    end
    for i = 1:numel(bounds)
        xc = crossings(px, pz, bounds(i));
        if numel(xc) >= 2
            h_bound = plot(ax,[min(xc) max(xc)],[bounds(i) bounds(i)],'-', ...
                 'Color',st.fill_palette.boundary,'HandleVisibility','off');
            mwecmass.output.figures.style_line(h_bound, st, 'boundary');
        end
    end
    h_outer = plot(ax, px, pz,'-','Color',st.fill_palette.boundary);
    mwecmass.output.figures.style_line(h_outer, st, 'boundary');
    xr = get_or(opts,'xlim',[min(px)-0.2, max(px)+0.2]);
    zr = get_or(opts,'ylim',[min(pz)-0.15, max(pz)+0.25]);
    h_wl = plot(ax, xr,[0 0], wl_style,'Color',st.fill_palette.waterline);
    mwecmass.output.figures.style_line(h_wl, st, 'reference');

    % Centre of gravity (red) and centre of buoyancy (blue), world frame, at
    % X = 0 on the axis of revolution.  White edge so they stay legible over
    % the dark end of the cividis ramp.
    z_cg = get_or(opts,'z_cg',[]);
    z_cb = get_or(opts,'z_cb',[]);
    dx   = 0.08*diff(xr);
    % fs_annot-scaled locally (GIF 2c's presentation call overrides it via opts); style_text's
    % 'annotation' role reads font_size.annotation, so it is overridden here rather than
    % silently dropping to the un-scaled default. The MarkerSize=11 on
    % these two point markers is left as a literal: it is this figure's deliberate CG/CB
    % emphasis size, distinct from the uniform style.marker.size style_line's marker branch
    % would apply.
    annot_style = st;
    annot_style.font_size.annotation = fs_annot;
    if ~isempty(z_cg) && isfinite(z_cg)
        plot(ax, 0, z_cg,'o','MarkerSize',11,'MarkerFaceColor',CG_RED, ...
             'MarkerEdgeColor','w','LineWidth',1.2,'HandleVisibility','off');
        h_cglbl = text(ax, -dx, z_cg,'$Z_{CG}$', ...
             'HorizontalAlignment','right','VerticalAlignment','middle','Color',CG_RED, ...
             'BackgroundColor','w','Margin',0.5);
        mwecmass.output.figures.style_text(h_cglbl, annot_style, 'annotation');
    end
    if ~isempty(z_cb) && isfinite(z_cb)
        plot(ax, 0, z_cb,'o','MarkerSize',11,'MarkerFaceColor',CB_BLUE, ...
             'MarkerEdgeColor','w','LineWidth',1.2,'HandleVisibility','off');
        h_cblbl = text(ax, dx, z_cb,'$Z_{CB}$', ...
             'HorizontalAlignment','left','VerticalAlignment','middle','Color',CB_BLUE, ...
             'BackgroundColor','w','Margin',0.5);
        mwecmass.output.figures.style_text(h_cblbl, annot_style, 'annotation');
    end

    xlim(ax,xr);  ylim(ax,zr);
    if ~isempty(get_or(opts,'xtick',[])), xticks(ax, opts.xtick); end
    if ~isempty(get_or(opts,'ytick',[])), yticks(ax, opts.ytick); end
    % apply_axes_style reads style.font_size.axes for XLabel/YLabel and style.font_size.tick_label
    % for the axes itself; fs_label here is the same "axes" role under this function's own name,
    % so xlabel/ylabel go through it unstyled and call_style below carries both fields.
    xlabel(ax, get_or(opts,'xlabel','$x$ [m]'));
    ylabel(ax, get_or(opts,'ylabel','$z$ [m]   (world frame, waterline = 0)'));
    cb = colorbar(ax,'Location','southoutside');
    cb.Label.String = get_or(opts,'cb_label','$\rho$ [kg/m$^3$]');
    % Only when asked: a colorbar takes its tick size from the axes at
    % creation, so setting it unconditionally would silently enlarge the
    % ticks in GIFs 1 and 2, which never pass this option.
    fs_cb = get_or(opts,'fs_cb',[]);
    call_style = st;
    call_style.font_size.tick_label = fs_axis;
    call_style.font_size.axes       = fs_label;
    call_style.font_size.title      = fs_title;
    % style_colorbar now reaches cb.Label too; it has no dedicated
    % font_size.label field, so its 'label' role falls back to font_size.annotation unless
    % overridden here -- done to keep this call's fs_label sizing (the "axes" role under this
    % function's own opts name) rather than silently drop to the annotation size.
    call_style.font_size.annotation = fs_label;
    if ~isempty(fs_cb), call_style.font_size.colorbar = fs_cb; end
    mwecmass.output.figures.apply_axes_style(ax, call_style);
    mwecmass.output.figures.style_colorbar(cb, call_style);
    % AFTER apply_axes_style: the colorbar tracks ax.Position, which `axis equal`
    % does not shrink here because Position already carries the data aspect.
    % Pinning it makes that guarantee explicit rather than incidental.
    cb_pos = get_or(opts,'cb_pos',[]);
    if ~isempty(cb_pos), cb.Position = cb_pos; end
end


function draw_uhpc_section(ax, prof, vs, zbal, cfg, g, is_sol, t_strip, Vu, Vv, st)
% Mirrors mwecmass.output.figures.plot_modular_precast FIGURE 1.
%   STEP A  offset the FULL profile inward by t_strip(i) (never a strip-clipped
%           piece -- offsetting the clipped piece's fake horizontal edges is
%           what produces shelves at strip boundaries), cached per thickness.
%   STEP B  clip the cached inner polygon to the strip, then contract it
%           HORIZONTALLY about its x-centroid until its 2D area equals
%           A_strip * V_void(i)/V_total(i).  z preserved -> no shelves.
    hold(ax,'on');  axis(ax,'equal');
    se = cfg.strip_edges(:);  Ns = numel(se)-1;
    px = prof(:,1);  pz = prof(:,2) + vs;

    zg = g.z(:);  Ao = g.A_outer(:);  Ai = g.A_inner(:);
    Ai(zg <= zbal) = 0;                     % z_ballast rule used by integrate_split

    % per-strip volumes for the STEP-B area target, and for the labels
    Vu_i = zeros(Ns,1);  Vv_i = zeros(Ns,1);  solid_i = false(Ns,1);
    for i = 1:Ns
        m = zg >= se(i) & zg <= se(i+1);
        if sum(m) >= 2
            Vu_i(i) = trapz(zg(m), max(Ao(m)-Ai(m),0));
            Vv_i(i) = trapz(zg(m), Ai(m));
        end
        solid_i(i) = is_sol(i) || se(i+1) <= zbal;
    end

    % STEP A: cache one inner polygon per unique thickness
    keys = {};  polys = {};
    for i = 1:Ns
        if solid_i(i), continue; end
        ti = t_strip(i);
        if ~isfinite(ti) || ti <= 0, continue; end
        kk = sprintf('%.6f', ti);
        if any(strcmp(keys,kk)), continue; end
        keys{end+1} = kk;                               %#ok<AGROW>
        polys{end+1} = inner_offset(px, pz, ti);        %#ok<AGROW>
    end

    for i = 1:Ns
        zlo = se(i)+vs;  zhi = se(i+1)+vs;
        out = clip_poly(px, pz, zlo, zhi);
        if size(out,1) < 3, continue; end

        if is_sol(i),            fc = st.fill_palette.solid_material;
        elseif solid_i(i),       fc = st.fill_palette.ballast_material;
        else,                    fc = st.fill_palette.jacket_material;
        end
        patch(ax, out(:,1), out(:,2), fc,'EdgeColor','none','FaceAlpha',1.0, ...
              'HandleVisibility','off');

        if ~solid_i(i)
            kk = sprintf('%.6f', t_strip(i));
            q  = find(strcmp(keys,kk),1);
            if ~isempty(q) && ~isempty(polys{q})
                P = polys{q};
                inn = clip_poly(P(:,1), P(:,2), zlo, zhi);
                if size(inn,1) >= 3
                    Vt = Vu_i(i) + Vv_i(i);
                    A_strip = polyarea(out(:,1), out(:,2));
                    A_void  = polyarea(inn(:,1), inn(:,2));
                    if Vt > 1e-12 && A_strip > 1e-9 && A_void > 1e-9
                        A_tgt = A_strip * (Vv_i(i)/Vt);
                        if A_void > A_tgt
                            sh = A_tgt / A_void;   cx = mean(inn(:,1));
                            inn(:,1) = (inn(:,1)-cx)*sh + cx;
                        end
                    end
                    patch(ax, inn(:,1), inn(:,2), st.fill_palette.void, ...
                          'EdgeColor',st.fill_palette.inner_boundary,'LineStyle','--','LineWidth',1.0, ...
                          'HandleVisibility','off');
                    mwecmass.output.figures.draw_hatch_strips(ax, inn(:,1), inn(:,2), ...
                          st.hatch_spacing, st.fill_palette.hatch);
                end
            end
        end

        % strip label, matching visualize's wording
        if is_sol(i)
            lbl = sprintf('S%d: 100%% UHPC (wall)', i);
        elseif solid_i(i)
            lbl = sprintf('S%d: 100%% UHPC (solid ballast)', i);
        else
            f = 100*Vu_i(i)/max(Vu_i(i)+Vv_i(i),eps);
            lbl = sprintf('S%d: %.0f%% UHPC / %.0f%% Void   t = %.2f in', ...
                          i, f, 100-f, t_strip(i)/0.0254);
        end
        % lbl's literal '%' characters (from '100%%'/'%.0f%%') would start a LaTeX comment under
        % the math interpreter style_text now forces, so they are escaped first.
        h_uhpclbl = text(ax, 0, 0.5*(zlo+zhi), strrep(lbl,'%','\%'), ...
             'HorizontalAlignment','center','BackgroundColor','w','Margin',0.5);
        mwecmass.output.figures.style_text(h_uhpclbl, st, 'annotation');
    end

    h_outer = plot(ax, px, pz,'-','Color',st.fill_palette.boundary);
    mwecmass.output.figures.style_line(h_outer, st, 'boundary');
    xr = [min(px)-0.2, max(px)+0.2];
    h_wl = plot(ax, xr,[0 0],'--','Color',st.fill_palette.waterline, ...
                'DisplayName','Waterline');
    mwecmass.output.figures.style_line(h_wl, st, 'reference');
    h_zbal = plot(ax, xr,[zbal+vs zbal+vs],'-','Color',st.fill_palette.ballast_level, ...
                'DisplayName','z_{ballast}');
    mwecmass.output.figures.style_line(h_zbal, st, 'reference');
    h_wb = plot(ax, xr,[se(cfg.wall_strip_index)+vs se(cfg.wall_strip_index)+vs], ...
                '-','Color',st.fill_palette.wall_boundary,'DisplayName','Wall boundary');
    mwecmass.output.figures.style_line(h_wb, st, 'boundary');

    hS = patch(ax,nan,nan,st.fill_palette.solid_material,'DisplayName','UHPC Solid Wall');
    hF = patch(ax,nan,nan,st.fill_palette.ballast_material,'DisplayName','UHPC Solid Ballast');
    hU = patch(ax,nan,nan,st.fill_palette.jacket_material,'DisplayName','UHPC Shell');
    hV = patch(ax,nan,nan,st.fill_palette.void,'EdgeColor',st.fill_palette.inner_boundary,'LineStyle','--', ...
               'DisplayName','Void (Air)');

    xlim(ax,xr);  ylim(ax,[min(pz)-0.15 max(pz)+0.25]);
    xlabel(ax,'$x$ [m]');
    ylabel(ax,'$z$ [m]   (world frame, waterline = 0)');
    title(ax,sprintf('$d$ = %.3f m,  $V_{UHPC}$ = %.3f m$^3$,  $V_{void}$ = %.3f m$^3$', ...
          abs(cfg.hull_z_min+vs), Vu, Vv));
    lg = legend(ax,[hS hF hU hV h_wl h_zbal h_wb],'Location','southoutside','NumColumns',4);
    mwecmass.output.figures.style_legend(lg, st);
    mwecmass.output.figures.apply_axes_style(ax, st);
end


%% =====================================================================
%%  MIRRORED / SHARED HELPERS
%% =====================================================================
function c = rho_rgb(rho, clim_rho, cmap)
    f = max(0, min(1, (rho-clim_rho(1))/max(diff(clim_rho),eps)));
    c = cmap(1+round(f*(size(cmap,1)-1)), :);
end

function P = clip_poly(x, z, z_lo, z_hi)
% Strip clip via the class's PUBLIC Sutherland-Hodgman helper.
    [xa, za] = mwecmass.internal.clip_z(x(:), z(:), z_lo, 'above');
    if numel(xa) < 3, P = zeros(0,2); return; end
    [xb, zb] = mwecmass.internal.clip_z(xa, za, z_hi, 'below');
    if numel(xb) < 3, P = zeros(0,2); else, P = [xb(:), zb(:)]; end
end

function xc = crossings(x, z, z0)
    xc = [];  n = numel(x);
    for j = 1:n
        j2 = mod(j,n)+1;
        z1 = z(j);  z2 = z(j2);
        if (z1-z0)*(z2-z0) <= 0 && abs(z2-z1) > 1e-12
            t = (z0-z1)/(z2-z1);
            xc(end+1) = x(j) + t*(x(j2)-x(j));   %#ok<AGROW>
        end
    end
end

function P = inner_offset(px, pz, t)
% Inward perpendicular offset of the FULL profile by t, via polyshape
% Minkowski erosion (the path visualize prefers).  Returns the largest
% resulting region, or empty if the profile erodes away entirely.
    P = zeros(0,2);
    try
        ps = polyshape(px(:), pz(:), 'Simplify', true, 'KeepCollinearPoints', false);
        pe = polybuffer(ps, -abs(t));
        if pe.NumRegions < 1, return; end
        R = regions(pe);
        [~, b] = max(area(R));
        v = R(b).Vertices;
        v = v(all(isfinite(v),2), :);
        if size(v,1) >= 3, P = v; end
    catch
        P = zeros(0,2);
    end
end

function band(ax, xb, rng, col, name)
    patch(ax,[xb fliplr(xb)],[rng(1) rng(1) rng(2) rng(2)], col, ...
          'FaceAlpha',0.09,'EdgeColor','none','DisplayName',name);
end

function trace_panel(ax, xv, k, series, names, styles, cols, yl, ylab, ttl, floor_v, floor_c, st)
    hold(ax,'on');
    if ~isempty(floor_v)
        patch(ax,[min(xv)-0.5 max(xv)+0.5 max(xv)+0.5 min(xv)-0.5], ...
                 [yl(1) yl(1) floor_v floor_v], floor_c, ...
                 'FaceAlpha',0.10,'EdgeColor','none','DisplayName','below floor');
        yline(ax, floor_v,'--','Color',floor_c,'LineWidth',1.4,'HandleVisibility','off');
    end
    for s = 1:numel(series)
        v = series{s};
        h_s = plot(ax, xv(1:k), v(1:k), styles{s},'Color',cols{s}, ...
             'MarkerSize',4.5,'DisplayName',names{s});
        mwecmass.output.figures.style_line(h_s, st, 'curve');
    end
    xlim(ax,[min(xv)-0.5 max(xv)+0.5]);  ylim(ax,yl);
    ylabel(ax, ylab);
    title(ax, ttl);
    lg = legend(ax,'Location','southeast');
    mwecmass.output.figures.style_legend(lg, st);
    mwecmass.output.figures.apply_axes_style(ax, st);
end

function fig = new_fig(st)
    fig = figure('Color','w','Position',[40 40 1500 820], ...
                 'Visible','off','Renderer','painters', ...
                 'DefaultAxesFontName',st.font_name);
end

function append_frame(fig, gif, is_first, delay, png_path)
% Without png_path: the original 90-dpi path (GIFs 1, 2, 2b, 3 -- unchanged).
% With png_path: render at 300 dpi, keep that full-resolution frame as a PNG,
% and downscale a copy for the GIF.  A 300-dpi frame of this figure is
% ~3875 x 2500 px; 34 of those in one GIF is >100 MB and unusable in a
% presentation, so the animation is capped at GIF_W px wide while the PNGs
% keep the full 300 dpi.
    GIF_W = 1600;                      % [px] animation width cap
    if nargin < 5 || isempty(png_path)
        im = print(fig,'-RGBImage','-r90');
        close(fig);                    % 16 GB machine: one figure at a time
    else
        im = print(fig,'-RGBImage','-r300');
        close(fig);
        imwrite(im, png_path);
        if size(im,2) > GIF_W
            im = imresize(im, [NaN GIF_W], 'lanczos3');
        end
    end
    [A,map] = rgb2ind(im,256);
    if is_first
        imwrite(A,map,gif,'gif','LoopCount',Inf,'DelayTime',delay);
    else
        imwrite(A,map,gif,'gif','WriteMode','append','DelayTime',delay);
    end
end

function L = pad(v, f)
    v = v(isfinite(v));
    L = [min(v) max(v)];  L = L + f*max(diff(L),eps)*[-1 1];
end

function v = get_or(s,f,d)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
