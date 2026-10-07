function test_baseline_vs_matlab()
%TEST_BASELINE_VS_MATLAB Print the recorded Octave baseline next to the MATLAB v1.0 values.
%   Reads only tests/baseline/*.json; no pipeline run. Octave's sqp is not MATLAB's, so the two
%   differ by construction and nothing is gated here except that every printed value exists.
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    matlab_ref = jsondecode(fileread(fullfile(root, 'tests', 'baseline', 'matlab_v1_reference.json')));
    files = {'octave_v1_baseline.json', 'full'; 'octave_v1_baseline_fast.json', 'fast'};
    scalars = {'vertical_shift', 'draft', 'mass_total', 'CG_total_z', 'GM_L', ...
               'T_heave_coupled', 'T_pitch_coupled', 'stage2_exitflag'};
    stage3_keys = {'z_ballast', 'exitflag', 'feasible'};
    for mode = {'modular_precast', 'thin_shell'}
        s2 = matlab_ref.(mode{1}).stage2;
        fprintf('%s, MATLAB v1.0 Stage-2 design (results.Final3D), the reference for Stage 3:\n', mode{1});
        for name = fieldnames(s2)'
            fprintf('  %-22s %16.10g\n', name{1}, s2.(name{1}));
        end
        fprintf('\n');
    end
    for f = 1:size(files, 1)
        recorded = jsondecode(fileread(fullfile(root, 'tests', 'baseline', files{f, 1})));
        for mode = {'modular_precast', 'thin_shell'}
            mode = mode{1};
            if ~isfield(recorded, mode), continue; end
            o = recorded.(mode);
            m = matlab_ref.(mode);
            fprintf('%s, preset %s: Octave run %.0f s, Octave %s\n', mode, files{f, 2}, ...
                    recorded.(['seconds_' mode]), recorded.info.octave);
            fprintf('  %-22s %16s %16s %14s\n', 'quantity', 'Octave', 'MATLAB v1.0', 'Octave-MATLAB');
            for k = 1:numel(scalars)
                print_row(scalars{k}, o.(scalars{k}), m.(scalars{k}));
            end
            for k = 1:numel(stage3_keys)
                print_row(['stage3.' stage3_keys{k}], o.stage3.(stage3_keys{k}), m.stage3.(stage3_keys{k}));
            end
            for name = {'stage1_x', 'stage2_x', 'realised_strip_density'}
                a = o.(name{1})(:); b = m.(name{1})(:);
                fprintf('  %-22s max |Octave - MATLAB| = %.6g\n', name{1}, max(abs(a - b)));
                fprintf('    Octave %s\n    MATLAB %s\n', mat2str(a', 8), mat2str(b', 8));
            end
            fprintf('\n');
        end
    end
end

function print_row(name, a, b)
    fprintf('  %-22s %16.8g %16.8g %14.6g\n', name, a, b, a - b);
end
