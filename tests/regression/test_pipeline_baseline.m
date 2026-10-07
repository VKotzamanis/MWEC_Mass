function test_pipeline_baseline()
%TEST_PIPELINE_BASELINE Rerun the Octave pipeline on C1 and require the baseline numbers exactly.
%   The pipeline runs through mwecmass.driver.run (tools/baseline_run.m) in the modes that the
%   selected baseline file holds. The code is deterministic, so every recorded number must
%   reproduce to the last bit: the gate is equality of jsonencode(summary), which prints each
%   double with the shortest text that reads back to the same double, against the member text
%   stored in the baseline file. There is no tolerance. The MATLAB v1.0 values are printed
%   alongside as information only: Octave's sqp is not MATLAB's, so those differ by construction.
%   The Octave run is long (see the seconds_<mode> members of the baseline files; most of it is
%   the hull-section tables in build_config, where Octave's containers.Map lookups dominate).
%   Environment variables:
%     TESTS_BASELINE_PRESET  'fast' (default): baseline_run preset with coarser z-grids, compared
%                            with octave_v1_baseline_fast.json; 'full': the author inputs,
%                            compared with octave_v1_baseline.json.
%     TESTS_BASELINE_MODES   comma separated modes to run, for example thin_shell; the default is
%                            every mode that the selected baseline file holds.
%   tests/run_tests.m runs this test only when MWEC_REGRESSION=1.
    root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(fullfile(root, 'tests', 'octave_shims'), fullfile(root, 'tools'));
    preset = getenv('TESTS_BASELINE_PRESET');
    if isempty(preset), preset = 'fast'; end
    suffix = '';
    if ~strcmp(preset, 'full'), suffix = ['_' preset]; end
    % A struct, not a containers.Map: baseline_run clears functions, which breaks Map objects.
    members = read_members(fullfile(root, 'tests', 'baseline', ['octave_v1_baseline' suffix '.json']));
    matlab_ref = jsondecode(fileread(fullfile(root, 'tests', 'baseline', 'matlab_v1_reference.json')));

    modes = intersect({'modular_precast', 'thin_shell'}, fieldnames(members)');
    selected = getenv('TESTS_BASELINE_MODES');
    if ~isempty(selected)
        modes = strsplit(selected, ',');
    end
    mismatches = {};
    for k = 1:numel(modes)
        mode = modes{k};
        if ~isfield(members, mode)
            error('test_pipeline_baseline:noBaseline', 'no baseline entry for %s', mode);
        end
        t0 = tic;
        summary = baseline_run(mode, false, preset);
        fprintf('%s (%s): pipeline run %.0f s (baseline run %s s)\n', mode, preset, toc(t0), members.(['seconds_' mode]));
        fresh_text = jsonencode(summary);
        if ~strcmp(fresh_text, members.(mode))
            found = compare(mode, jsondecode(members.(mode)), jsondecode(fresh_text));
            if isempty(found)
                found = {[mode ': texts differ but decode to the same doubles']};
            end
            mismatches = [mismatches, found]; %#ok<AGROW>
        end
        print_reference(mode, jsondecode(fresh_text), matlab_ref.(mode));
    end
    if ~isempty(mismatches)
        error('test_pipeline_baseline:differs', 'numbers differ from the Octave baseline:\n  %s', ...
              strjoin(mismatches, sprintf('\n  ')));
    end
    fprintf('all baseline numbers identical\n');
end

function members = read_members(path)
% One JSON member per line, written by tools/write_octave_baseline.m.
    members = struct();
    lines = strsplit(fileread(path), sprintf('\n'));
    for k = 1:numel(lines)
        tok = regexp(lines{k}, '^"([^"]+)": (.*?),?$', 'tokens', 'once');
        if ~isempty(tok), members.(tok{1}) = tok{2}; end
    end
end

function bad = compare(prefix, expected, fresh)
    bad = {};
    names = fieldnames(expected);
    for k = 1:numel(names)
        n = names{k};
        label = [prefix '.' n];
        if ~isfield(fresh, n)
            bad{end+1} = [label ' missing']; %#ok<AGROW>
        elseif isstruct(expected.(n))
            bad = [bad, compare(label, expected.(n), fresh.(n))]; %#ok<AGROW>
        elseif ~isequal(expected.(n), fresh.(n))
            if isnumeric(expected.(n)) && isnumeric(fresh.(n)) && isequal(size(expected.(n)), size(fresh.(n)))
                detail = sprintf('max abs diff %.3e', max(abs(expected.(n)(:) - fresh.(n)(:))));
            else
                detail = 'different value';
            end
            bad{end+1} = sprintf('%s: %s', label, detail); %#ok<AGROW>
        end
    end
end

function print_reference(mode, fresh, ref)
    fprintf('%s, Octave vs MATLAB v1.0 (information, not asserted):\n', mode);
    keys = {'mass_total', 'CG_total_z', 'GM_L', 'T_heave_coupled', 'T_pitch_coupled', ...
            'vertical_shift', 'stage2_exitflag'};
    for k = 1:numel(keys)
        fprintf('  %-16s Octave % .8g   MATLAB % .8g\n', keys{k}, fresh.(keys{k}), ref.(keys{k}));
    end
    for name = {'stage1_x', 'stage2_x'}
        fprintf('  %-16s max |Octave - MATLAB| = %.3e\n', name{1}, ...
                max(abs(fresh.(name{1})(:) - ref.(name{1})(:))));
    end
    for name = {'z_ballast', 'exitflag', 'feasible'}
        if isfield(fresh.stage3, name{1}) && isfield(ref.stage3, name{1})
            fprintf('  stage3.%-9s Octave % .8g   MATLAB % .8g\n', name{1}, ...
                    fresh.stage3.(name{1}), ref.stage3.(name{1}));
        end
    end
end
