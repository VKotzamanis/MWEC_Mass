classdef WEC_Mesh_Sizing
    % WEC_MESH_SIZING  Pick BEM panelizer (Nu, Nv) from a target panel size.
    %
    %   Given an MS2 parser and a target physical panel edge length
    %   panel_size [m], compute global (Nu, Nv) such that the mean panel
    %   side length in each parametric direction is approximately
    %   panel_size on the geometrically dominant surface, and the mean
    %   panel aspect ratio is ~ 1.
    %
    %   STRATEGY
    %   ────────────────────────────────────────────────────────────────
    %   For each source (non-mirror) surface, evaluate the parametric
    %   map r(u,v) on a coarse probe grid.  Compute:
    %     S_u(surf) = mean over v of the polyline arc length of an
    %                 iso-v curve (parameterized by u).
    %     S_v(surf) = mean over u of the polyline arc length of an
    %                 iso-u curve (parameterized by v).
    %
    %   These are the average physical lengths "along u" and "along v"
    %   on the surface.  Choose:
    %     Nu = ceil(max_s S_u(s) / panel_size)
    %     Nv = ceil(max_s S_v(s) / panel_size)
    %
    %   With uniform parametric spacing this gives a mean panel side
    %   length on each surface s of:
    %     <Δ_u>_s ≈ S_u(s)/Nu,    <Δ_v>_s ≈ S_v(s)/Nv
    %   so the worst-case surface gets ~ panel_size; others get smaller
    %   panels (still mean aspect ratio ≈ 1 because both N's scale
    %   together).
    %
    %   WHY GLOBAL Nu/Nv
    %   ────────────────────────────────────────────────────────────────
    %   WEC_Panelizer.generate uses a single u_grid / v_grid across every
    %   source surface to guarantee watertight junctions at shared edges
    %   (matched parametric subdivision → bit-identical 3D vertices →
    %   exact merge).  See the comment block at WEC_Panelizer.m §3
    %   "WHY this guarantees watertight junctions".  A strict per-surface
    %   (Nu, Nv) would require a parallel refactor — out of scope here.
    %
    %   USAGE
    %     [Nu, Nv, report] = WEC_Mesh_Sizing.compute(parser, panel_size)
    %     [Nu, Nv, report] = WEC_Mesh_Sizing.compute(parser, panel_size, opts)
    %
    %   OPTIONS (struct, all optional)
    %     .probe_n    — probe grid size in each parametric direction
    %                   used for arc length estimation (default: 30)
    %     .n_min      — floor on Nu/Nv (default: 3)
    %     .verbose    — print per-surface report (default: true)
    %
    %   OUTPUT
    %     Nu, Nv      — chosen panel grid density
    %     report      — struct with fields:
    %                     .panel_size, .sources, .S_u, .S_v,
    %                     .panel_u_per_surf, .panel_v_per_surf,
    %                     .aspect_per_surf, .Nu, .Nv

    methods (Static)

        function [Nu, Nv, report] = compute(parser, panel_size, opts)
            if nargin < 3, opts = struct(); end
            if ~isfield(opts, 'probe_n'), opts.probe_n = 30;   end
            if ~isfield(opts, 'n_min'),   opts.n_min   = 3;    end
            if ~isfield(opts, 'verbose'), opts.verbose = true; end

            assert(panel_size > 0, 'panel_size must be positive');

            topo    = parser.classify_visible_surfaces();
            sources = topo.sources;
            n_src   = numel(sources);

            assert(n_src >= 1, ...
                'WEC_Mesh_Sizing:NoSources', ...
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
            % This satisfies the IRFR rule documented in WEC_Driver §2
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

    end
end
