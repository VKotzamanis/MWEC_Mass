function test_thin_shell_deletions()
%TEST_THIN_SHELL_DELETIONS  Plan deletions D2, D4, D6, D7, D17 on the thin-shell path.
%   The thin-shell realisation holds none of the deleted names and the replaced helpers are gone.
%   Outside it, the names may remain only in files other tasks own at this base (contract section
%   4 and the deletion register): they are listed in `pending` with their owner and printed; any
%   other file holding one fails the test. J2 and T11 shrink the list as those files change.

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
names = {'compute_rmin_at_z', 'hull_slope_cos_at_z', 'inner_properties_at_z', 'max_slope_factor', ...
    't_min_active', 'thin_shell.build_geometry_grid', 'strip_partition_volumes', 'steel_n_z_grid', ...
    'cos_alpha', '0.95*t_max', '0.9*t_max'};
pending = {
    'src/+mwecmass/+driver/build_config.m', 'T4b (D1, D3), then the T7 forwards at J2'
    'WEC_User_Input.m', 'T5, then the T7 inputs at J2'
    'src/+mwecmass/+geometry/compute_rmin_at_z.m', 'D2 file, deleted once T4b, T5 and T6 removed its callers'
    'src/+mwecmass/+hydrostatics/compute_perpendicular_shell_volume.m', 'T4b (D1)'
    'src/+mwecmass/+output/+figures/plot_steel_solve.m', 'T8 (D8)'
    'src/+mwecmass/+output/export_schema.m', 'T5 (S8 replaces steel_data)'
    'src/+mwecmass/+realise/+modular_precast/build_geometry_grid.m', 'T5, T6'
    'src/+mwecmass/+realise/+modular_precast/extract_strip_geometry.m', 'T5 (D5)'
    'src/+mwecmass/+realise/+modular_precast/solve.m', 'T6'
    'src/+mwecmass/+realise/+modular_precast/solve_and_extract.m', 'T5'
    'validation/diagnostics/stage_animations.m', 'T8 (D8)'
    'tools/baseline_run.m', 'T5'
    };
for f = {'build_geometry_grid', 'hull_slope_cos_at_z', 'inner_properties_at_z', 'strip_partition_volumes'}
    p = fullfile(root, 'src', '+mwecmass', '+realise', '+thin_shell', [f{1} '.m']);
    check(exist(p, 'file') ~= 2, 'thin_shell/%s.m still exists', f{1});
end
files = [list_m(fullfile(root, 'src')), list_m(fullfile(root, 'validation')), list_m(fullfile(root, 'tools')), ...
    {fullfile(root, 'WEC_User_Input.m'), fullfile(root, 'WEC_Output_Options.m')}];
n_pending = 0;
for k = 1:numel(files)
    txt = fileread(files{k});
    hits = names(cellfun(@(n) ~isempty(strfind(txt, n)), names));
    if isempty(hits)
        continue
    end
    rel = strrep(files{k}(numel(root) + 2:end), filesep, '/');
    row = find(strcmp(pending(:, 1), rel), 1);
    check(~isempty(row), '%s holds %s', rel, strjoin(hits, ', '));
    check(isempty(strfind(rel, '+thin_shell/')), '%s holds %s', rel, strjoin(hits, ', '));
    n_pending = n_pending + 1;
    fprintf('pending (%s): %s holds %s\n', pending{row, 2}, rel, strjoin(hits, ', '));
end
fprintf('%d files outside the thin-shell path still hold a deleted name\n', n_pending);
end

function files = list_m(folder)
files = {};
entries = dir(folder);
for k = 1:numel(entries)
    e = entries(k);
    if any(strcmp(e.name, {'.', '..'}))
        continue
    end
    p = fullfile(folder, e.name);
    if e.isdir
        files = [files, list_m(p)]; %#ok<AGROW>
    elseif numel(e.name) > 2 && strcmp(e.name(end - 1:end), '.m')
        files{end + 1} = p; %#ok<AGROW>
    end
end
end

function check(cond, varargin)
if ~cond
    error('test_thin_shell_deletions:fail', varargin{:});
end
end
