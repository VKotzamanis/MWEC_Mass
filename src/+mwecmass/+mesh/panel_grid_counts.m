function [Nu, Nv, report] = panel_grid_counts(parser, panel_size, opts)
%PANEL_GRID_COUNTS Estimate parametric grid counts from sampled physical surface lengths.
% Syntax: [Nu,Nv,report] = panel_grid_counts(parser,panel_size,opts).
% Inputs: parser geometry parser; panel_size [m] target edge; opts probe_n/n_min/verbose.
% Outputs: Nu,Nv scalar counts (at least n_min) and report with sampled lengths/aspect data.

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'probe_n'), opts.probe_n = 30;   end
    if ~isfield(opts, 'n_min'),   opts.n_min   = 3;    end
    if ~isfield(opts, 'verbose'), opts.verbose = true; end

    assert(panel_size > 0, 'panel_size must be positive');

    topo    = parser.classify_visible_surfaces();
    sources = topo.sources;
    n_src   = numel(sources);

    assert(n_src >= 1, ...
        'mwecmass:mesh:NoSources', ...
        'No source (non-mirror) surfaces found in geometry.');

    u_probe = linspace(0, 1, opts.probe_n);
    v_probe = linspace(0, 1, opts.probe_n);

    S_u_all = zeros(n_src, 1);
    S_v_all = zeros(n_src, 1);

    for k = 1:n_src
        sname = sources{k};
        G = parser.eval_surface_grid(sname, u_probe, v_probe);
        % G : [probe_n × probe_n × 3]

        % Iso-v curves are parameterized by u → arc length in u.
        % Polyline length: sum of chord lengths along i-axis.
        dU = diff(G, 1, 1);                          % [(N-1) × N × 3]
        len_u_seg = sqrt(sum(dU.^2, 3));             % [(N-1) × N]
        S_u_per_v = sum(len_u_seg, 1);               % [1 × N]
        S_u_all(k) = mean(S_u_per_v);

        dV = diff(G, 1, 2);                          % [N × (N-1) × 3]
        len_v_seg = sqrt(sum(dV.^2, 3));             % [N × (N-1)]
        S_v_per_u = sum(len_v_seg, 2);               % [N × 1]
        S_v_all(k) = mean(S_v_per_u);
    end

    Nu = max(opts.n_min, ceil(max(S_u_all) / panel_size));
    Nv = max(opts.n_min, ceil(max(S_v_all) / panel_size));

    % Per-surface realised panel sizes and aspect ratio.
    panel_u_per_surf = S_u_all / Nu;
    panel_v_per_surf = S_v_all / Nv;
    % Aspect ratio defined as max(Δu,Δv) / min(Δu,Δv) so values
    % are always ≥ 1 (1 = isotropic panel).
    aspect_per_surf  = max(panel_u_per_surf, panel_v_per_surf) ./ ...
                       max(min(panel_u_per_surf, panel_v_per_surf), eps);

    % Waterplane lid target edge — matched to panel_size so the WP
    % lid panels are comparable to hull panels at the waterline.
    % This satisfies the irregular-frequency-removal rule
    % ("wp_target_edge must be COMPARABLE to the hull panel edge").
    % Average hull panel edge ≈ panel_size by construction of the
    % sizer (mean Δu and mean Δv both targeted at panel_size).
    wp_target_edge = panel_size;

    report.panel_size        = panel_size;
    report.sources           = sources;
    report.S_u               = S_u_all;
    report.S_v               = S_v_all;
    report.panel_u_per_surf  = panel_u_per_surf;
    report.panel_v_per_surf  = panel_v_per_surf;
    report.aspect_per_surf   = aspect_per_surf;
    report.Nu                = Nu;
    report.Nv                = Nv;
    report.wp_target_edge    = wp_target_edge;

    if opts.verbose
        fprintf('\n========== MESH SIZING ==========\n');
        fprintf('  Target panel size : %.3f m\n', panel_size);
        fprintf('  Source surfaces   : %d\n', n_src);
        fprintf('  %-20s %10s %10s %10s %10s %8s\n', ...
                'surface', 'S_u [m]', 'S_v [m]', ...
                'Δu [m]', 'Δv [m]', 'aspect');
        for k = 1:n_src
            fprintf('  %-20s %10.3f %10.3f %10.3f %10.3f %8.2f\n', ...
                    sources{k}, S_u_all(k), S_v_all(k), ...
                    panel_u_per_surf(k), panel_v_per_surf(k), ...
                    aspect_per_surf(k));
        end
        fprintf('  ---------------------------------------------------------------------\n');
        fprintf('  Chosen Nu = %d   (max S_u = %.3f m / panel_size)\n', ...
                Nu, max(S_u_all));
        fprintf('  Chosen Nv = %d   (max S_v = %.3f m / panel_size)\n', ...
                Nv, max(S_v_all));
        fprintf('  WP target edge    = %.3f m   (= panel_size, IRFR-matched)\n', ...
                wp_target_edge);
        fprintf('  Total panels per source surf : %d\n', (Nu-1)*(Nv-1));
        fprintf('==================================\n\n');
    end
end
