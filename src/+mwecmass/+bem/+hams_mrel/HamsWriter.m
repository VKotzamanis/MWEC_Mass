classdef HamsWriter
%HAMSWRITER Static-method writers for the four fixed-format HAMS-MREL solver input files.
% Provides stateless writers for fixed-format HAMS input files.
    methods (Static)
        function matrix_6x6(fid, M)
        %MATRIX_6X6 Write a 6x6 matrix in HAMS Hydrostatic.in row format.
        % Fortran expects six fields using 2x,E12.5; the explicit two-space
        % prefix and %12.5E produce the required 14-character fields.
            for i = 1:6
                fprintf(fid, '  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E  %12.5E\n', ...
                    M(i,1), M(i,2), M(i,3), M(i,4), M(i,5), M(i,6));
            end
        end

        function control_file(filepath, params)
        %CONTROL_FILE Write a HAMS ControlFile.in solver-control file.
        % PARAMS contains depth [m], frequency/heading settings, reference
        % body center [3x1 m], body length [m], solver flags, and field points
        % [N×3 m]. Negative frequency/heading counts select auto-ranging.

            % CRITICAL: label widths match Fortran READ field widths in InputFiles.f90.
            % Verified against write_hams_controlfile.m (unit-tested).
            % Wrong widths cause HAMS to silently read 0 for IRSP/NTHREAD/ISOL.
            fid = fopen(filepath, 'w');
            if fid == -1
                error('Cannot open file for writing: %s', filepath);
            end

            % Lines 1-2: consumed by read(*) x2
            fprintf(fid, '   --------------HAMS Control file---------------\n');
            fprintf(fid, '\n');

            % Waterdepth  [14x,f30.15]
            fprintf(fid, '%-14s%30.15f\n', '   Waterdepth ', params.depth);
            fprintf(fid, '\n');
            fprintf(fid, '   #Start Definition of Wave Frequencies\n');

            % SYBO  [27x,i16]
            fprintf(fid, '%-27s%16d\n', '    0_inf_frequency_limits ', params.zero_inf_limits);
            % INFT  [25x,i16]
            fprintf(fid, '%-25s%16d\n', '    Input_frequency_type ', params.input_freq_type);
            % OUFT  [25x,i17]
            fprintf(fid, '%-25s%17d\n', '    Output_frequency_type', params.output_freq_type);
            % NPET  [26x,i16]
            fprintf(fid, '%-26s%16d\n', '    Number_of_frequencies ', params.n_frequencies);
            if params.n_frequencies < 0
                % WK1  [27x,f30.15]
                fprintf(fid, '%-27s%30.15f\n', '    Minimum_frequency_Wmin ', params.min_frequency);
                % DWK  [19x,f30.15]
                fprintf(fid, '%-19s%30.15f\n', '    Frequency_step ', params.freq_step);
            else
                freq_str = sprintf('%.6f ', params.frequency_list);
                fprintf(fid, '    %s\n', strtrim(freq_str));
            end
            fprintf(fid, '   #End Definition of Wave Frequencies\n');
            fprintf(fid, '\n');

            % Wave headings — always use auto-range (n_headings < 0)
            fprintf(fid, '   #Start Definition of Wave Headings\n');
            if params.n_headings < 0
                % NBETA  [23x,i16]
                fprintf(fid, '%-23s%16d\n', '    Number_of_headings ', params.n_headings);
                % BETA1  [20x,f30.15]
                fprintf(fid, '%-20s%30.15f\n', '    Minimum_heading ', params.min_heading);
                % DBETA  [17x,f30.15]
                fprintf(fid, '%-17s%30.15f\n', '    Heading_step ', params.heading_step);
            else
                % Positive n_headings: emit as auto-range with first value
                fprintf(fid, '%-23s%16d\n', '    Number_of_headings ', -1);
                h0 = params.heading_values(1);
                fprintf(fid, '%-20s%30.15f\n', '    Minimum_heading ', h0);
                fprintf(fid, '%-17s%30.15f\n', '    Heading_step ', 90.0);
            end
            fprintf(fid, '   #End Definition of Wave Headings\n');
            fprintf(fid, '\n');

            % XR   [28x,3f12.3]
            fprintf(fid, '%-28s%12.3f%12.3f%12.3f\n', ...
                '    Reference_body_center   ', ...
                params.ref_body_center(1), params.ref_body_center(2), params.ref_body_center(3));
            % REFL  [26x,f30.15]
            fprintf(fid, '%-26s%30.15f\n', '    Reference_body_length ', params.ref_body_length);
            % ISOL  [26x,i16]  — label EXACTLY 26 chars
            fprintf(fid, '%-26s%16d\n', '    Wave_diffrac_solution ', params.wave_diffrac_soln);
            % IRSP  [23x,i16]  — label EXACTLY 23 chars
            fprintf(fid, '%-23s%16d\n', '    If_remove_irr_freq ', params.remove_irr_freq);
            % NTHREAD [23x,i16]  — label EXACTLY 23 chars
            fprintf(fid, '%-23s%16d\n', '    Number of threads  ', params.n_threads);
            fprintf(fid, '\n');
            fprintf(fid, '   #Start Definition of Pressure and/or Elevation (PE)\n');

            % NFP  [27x,i16]
            fprintf(fid, '%-27s%16d\n', '    Number_of_field_points ', params.n_field_points);
            for i = 1:params.n_field_points
                fprintf(fid, '    %f    %f    %f    Global_coords_point_%d\n', ...
                    params.field_points(i,1), params.field_points(i,2), ...
                    params.field_points(i,3), i);
            end
            % Three lines required by InputFiles.f90 after NFP (lines 190,191,193)
            fprintf(fid, '   #End Definition of Pressure and/or Elevation\n');
            fprintf(fid, '   ----------End HAMS Control file----------\n');
            fprintf(fid, '\n');

            fclose(fid);
            fprintf('  Wrote %s\n', filepath);
        end

        function hydrostatic_file(filepath, CG, M_6x6, B_ext_lin, B_ext_quad, C_hydro, K_ext)
        %HYDROSTATIC_FILE Write HAMS Hydrostatic.in mass, damping, and stiffness.
        % CG is [3x1 m]; all matrix inputs are [6x6] (SI units), comprising
        % body mass, external linear/quadratic damping, hydrostatic restoring,
        % and external restoring (PTO/mooring).

            fid = fopen(filepath, 'w');
            if fid == -1
                error('Cannot open file for writing: %s', filepath);
            end

            % CG — Fortran reads with format similar to matrices.
            %   Use E12.5 fields (12 chars, 5 decimal digits) with
            %   2-space prefix per value = 14 chars per field.
            %   %.15E produced 21-char numbers → Fortran parse failure.
            fprintf(fid, ' Center of Gravity:\n');
            fprintf(fid, '  %12.5E  %12.5E  %12.5E\n', CG(1), CG(2), CG(3));

            % Mass matrix
            fprintf(fid, ' Body Mass Matrix:\n');
            mwecmass.bem.hams_mrel.HamsWriter.matrix_6x6(fid, M_6x6);

            % External linear damping
            fprintf(fid, ' External Linear Damping Matrix:\n');
            mwecmass.bem.hams_mrel.HamsWriter.matrix_6x6(fid, B_ext_lin);

            % External quadratic damping
            fprintf(fid, ' External Quadratic Damping Matrix:\n');
            mwecmass.bem.hams_mrel.HamsWriter.matrix_6x6(fid, B_ext_quad);

            % Hydrostatic restoring
            fprintf(fid, ' Hydrostatic Restoring Matrix:\n');
            mwecmass.bem.hams_mrel.HamsWriter.matrix_6x6(fid, C_hydro);

            % External restoring
            fprintf(fid, ' External Restoring Matrix:\n');
            mwecmass.bem.hams_mrel.HamsWriter.matrix_6x6(fid, K_ext);

            fclose(fid);
            fprintf('  Wrote %s\n', filepath);
        end

        function pnl_file(filepath, nodes, panels, panel_nverts, x_sym, y_sym)
        %PNL_FILE Write a HAMS .pnl mesh (nodes, connectivity, symmetry flags).
        % nodes is [N×3], panels is [P×4], and panel_nverts is 3 or 4 per row;
        % x_sym and y_sym select HAMS symmetry conventions.

            n_nodes = size(nodes, 1);
            n_panels = size(panels, 1);

            fid = fopen(filepath, 'w');
            if fid == -1
                error('Cannot open file for writing: %s', filepath);
            end

            fprintf(fid, '    --------------Hull Mesh File---------------\n');
            fprintf(fid, ' \n');
            fprintf(fid, '    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry\n');
            fprintf(fid, '    %8d    %8d       %5d       %5d\n', ...
                n_panels, n_nodes, x_sym, y_sym);
            fprintf(fid, ' \n');

            % Node coordinates
            fprintf(fid, '    # Start Definition of Node Coordinates     ! node_number   x   y   z\n');
            for i = 1:n_nodes
                fprintf(fid, ' %4d    %14.6f    %14.6f    %14.6f\n', ...
                    i, nodes(i,1), nodes(i,2), nodes(i,3));
            end
            fprintf(fid, '    # End Definition of Node Coordinates\n');
            fprintf(fid, ' \n');

            % Panel connectivity
            fprintf(fid, '  # Start Definition of Node Relations   ! panel_number  number_of_vertices   Vertex1_ID   Vertex2_ID   Vertex3_ID   (Vertex4_ID\n');
            for i = 1:n_panels
                nv = panel_nverts(i);
                if nv == 4
                    fprintf(fid, ' %4d    %d    %6d    %6d    %6d    %6d\n', ...
                        i, 4, panels(i,1), panels(i,2), panels(i,3), panels(i,4));
                else
                    fprintf(fid, ' %4d    %d    %6d    %6d    %6d\n', ...
                        i, 3, panels(i,1), panels(i,2), panels(i,3));
                end
            end
            fprintf(fid, '    # End Definition of Node Relations\n');
            fprintf(fid, ' \n');
            fprintf(fid, '    --------------End Hull Mesh File---------------\n');

            fclose(fid);
            fprintf('  Wrote %s (%d nodes, %d panels)\n', filepath, n_nodes, n_panels);
        end
    end
end
