% run_buildClimateGrid.m
% Driver: parse a WAM 2-D spectra file, build the (Hs, Te) climate grid, and
% write <station>_climate_grid.mat per STABILITY_HANDOFF_PLAN.md §4.
%
% Usage patterns
%   1) Interactive: just run the script; you will be prompted for the WAM file,
%      site metadata, and water depth.
%   2) Scripted: define the variables in the CONFIG block below before running,
%      and the prompts wilSTl be skipped.
%
% Outputs
%   <wam_dir>/<station>_climate_grid.mat   (contains top-level var 'climateGrid')
%   prints a one-line summary to the console
%
% Standards: see STABILITY_HANDOFF_PLAN.md §10.

clear; clc;

% =====================================================================
% USER-EDITABLE BIN WIDTHS
% Change these to alter the climate-grid resolution and re-run.
% Smaller widths -> more cells, fewer records/cell, finer plots.
% =====================================================================
dbinHs = 0.5;    % Hs bin width [m]   (default 0.5)
dbinTe = 1.0;    % Te bin width [s]   (default 1.0)
% =====================================================================

% ----- CONFIG (override these variables to skip prompts) ---------------------
% wam_file        = '/path/to/station.dat';
% station_id      = 'ST63044';
% region          = 'North Atlantic';
% lat             = 41.96;   lon = -67.31;
% water_depth_m   = 75;
% out_dir         = '';      % default: same folder as wam_file
%
% Optional grid overrides (bin widths: edit dbinHs / dbinTe at the top of this file):
% Hs_envelope_pctl = 98;
% omega_min = 0.30; omega_max = 6.00; omega_step = 0.05;
% gamma_default = 3.3; min_records_for_fit = 10;
% -----------------------------------------------------------------------------

helper_dir = fullfile(fileparts(mfilename('fullpath')), 'MATLAB Script for Plotting');
if exist(helper_dir, 'dir')
    addpath(helper_dir);
else
    error('run_buildClimateGrid:helpers', ...
          'Helper folder not found: %s', helper_dir);
end

if ~exist('wam_file','var') || isempty(wam_file)
    [fn, fp] = uigetfile( ...
        {'*.wam;*.dat;*.spc;*.txt','WAM 2-D spectra or WIS OneLine file (*.wam, *.dat, *.spc, *.txt)'; ...
         '*.wam','WAM 2-D spectra (*.wam)'; ...
         '*.*',  'All files'}, ...
        'Select wave-conditions input file(s) — hold Ctrl/Shift for multi-select', ...
        'MultiSelect', 'on');
    if isequal(fn, 0); fprintf('Cancelled.\n'); return; end
    if ischar(fn)
        wam_file = fullfile(fp, fn);
    else
        wam_file = cellfun(@(x) fullfile(fp, x), fn, 'UniformOutput', false);
    end
end

% Normalise wam_file into a cellstr list of paths
if ischar(wam_file) || isstring(wam_file)
    wam_files = {char(wam_file)};
elseif iscell(wam_file)
    wam_files = cellfun(@char, wam_file(:)', 'UniformOutput', false);
else
    error('run_buildClimateGrid:wam','wam_file must be a path string or cellstr of paths');
end
wam_files = sort(wam_files);   % deterministic ordering (e.g. by month suffix)
for ii = 1:numel(wam_files)
    if ~exist(wam_files{ii}, 'file')
        error('run_buildClimateGrid:wam','Input file not found: %s', wam_files{ii});
    end
end
wam_file_primary = wam_files{1};   % used for station-id default and out_dir

if ~exist('station_id','var') || isempty(station_id)
    [~, bases] = cellfun(@fileparts, wam_files, 'UniformOutput', false);
    % Common alphanumeric prefix across all selected basenames (e.g. "ST84040")
    common_prefix = bases{1};
    for ii = 2:numel(bases)
        n = min(numel(common_prefix), numel(bases{ii}));
        m = find(common_prefix(1:n) ~= bases{ii}(1:n), 1, 'first');
        if isempty(m), m = n + 1; end
        common_prefix = common_prefix(1:m-1);
    end
    common_prefix = regexprep(common_prefix, '[_\W]+$', '');  % trim trailing _ or punct
    if isempty(common_prefix)
        default_id = regexprep(bases{1}, '\W', '_');
    else
        default_id = regexprep(common_prefix, '\W', '_');
    end
    station_id = input(sprintf('station_id [%s]: ', default_id), 's');
    if isempty(station_id), station_id = default_id; end
end
if ~exist('region','var'),        region = input('region (e.g. North Atlantic): ', 's'); end
if ~exist('lat','var'),           lat = readScalarPrompt('lat (deg)', NaN); end
if ~exist('lon','var'),           lon = readScalarPrompt('lon (deg)', NaN); end
if ~exist('water_depth_m','var'), water_depth_m = readScalarPrompt('water_depth_m', NaN); end
if ~exist('out_dir','var') || isempty(out_dir)
    out_dir = fileparts(wam_file_primary);
end

fprintf('\n[1/3] Parsing %d wave file(s) ...\n', numel(wam_files));
records = [];
for ii = 1:numel(wam_files)
    fpath_ii = wam_files{ii};
    fmt = detectWaveFormat(fpath_ii);
    [~, base_ii, ext_ii] = fileparts(fpath_ii);
    fprintf('       (%d/%d) %s%s  [%s]\n', ii, numel(wam_files), base_ii, ext_ii, fmt);
    switch fmt
        case 'oneline'
            parsed  = parseWisOneline(fpath_ii);
            rec_ii  = onelineToRecords(parsed);
        case 'wam'
            rec_ii  = parseWamSpectra(fpath_ii);
            if iscell(rec_ii); rec_ii = [rec_ii{:}]; end
        otherwise
            error('run_buildClimateGrid:fmt','Unrecognised input format for %s', fpath_ii);
    end
    if isempty(records)
        records = rec_ii(:);
    else
        records = [records; rec_ii(:)]; %#ok<AGROW>
    end
end
fprintf('       %d records parsed across %d file(s).\n', numel(records), numel(wam_files));

opts = struct();
opts.station_id    = station_id;
opts.region        = region;
opts.lat           = lat;
opts.lon           = lon;
opts.water_depth_m = water_depth_m;
opts.dH = dbinHs;
opts.dT = dbinTe;
optsFromBase = {'Hs_envelope_pctl','Hs_max','Te_max', ...
                'omega_min','omega_max','omega_step', ...
                'gamma_default','min_records_for_fit'};
for k = 1:numel(optsFromBase)
    name = optsFromBase{k};
    if evalin('base', sprintf('exist(''%s'',''var'')', name))
        opts.(name) = evalin('base', name);
    end
end

fprintf('[2/4] Building climate grid (Hs x Te) ...\n');
climateGrid = buildClimateGrid(records, opts);
fprintf('       grid: %d x %d cells; Hs_P98 = %.2f m; %d operational cells.\n', ...
        numel(climateGrid.Hs_centers), numel(climateGrid.Te_centers), ...
        climateGrid.Hs_P98, nnz(climateGrid.operational_mask));

out_path = fullfile(out_dir, sprintf('%s_climate_grid.mat', station_id));
fprintf('[3/4] Writing %s ...\n', out_path);
report = writeClimateGridMat(climateGrid, out_path);
fprintf('       %s\n', report.summary);

fprintf('[4/4] Rendering diagnostic plots ...\n');
try
    plot_opts = struct();
    if exist('plot_visible','var') && ~isempty(plot_visible)
        plot_opts.visible = plot_visible;
    end
    plotClimateGrid(climateGrid, plot_opts);
catch ME_plot
    warning('run_buildClimateGrid:plot', ...
        'plotClimateGrid failed (grid was written OK): %s', ME_plot.message);
end

fprintf('\nDone.\n');

% --- local helpers -----------------------------------------------------------

function v = readScalarPrompt(label, default)
    txt = input(sprintf('%s [%g]: ', label, default), 's');
    if isempty(txt)
        v = default;
    else
        v = str2double(txt);
        if isnan(v); v = default; end
    end
end

function fmt = detectWaveFormat(filename)
% Heuristic: WIS OneLine files are pure-numeric one-record-per-row (>= 20 cols);
% WAM 2-D spectra files have header lines with embedded text or short numeric
% lines indicating (ML, KL) band counts.
    fid = fopen(filename, 'r');
    if fid < 0, error('detectWaveFormat:open','cannot open %s', filename); end
    cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
    line1 = '';
    for tries = 1:5
        L = fgetl(fid);
        if ~ischar(L); break; end
        if ~isempty(strtrim(L)); line1 = L; break; end
    end
    if isempty(line1)
        fmt = 'unknown'; return
    end
    tokens = strsplit(strtrim(line1));
    nums = str2double(tokens);
    if numel(tokens) >= 20 && all(~isnan(nums))
        fmt = 'oneline';
    else
        fmt = 'wam';
    end
end

function records = onelineToRecords(parsed)
% Adapter: WIS OneLine parsed struct (.All.* column arrays) -> record struct
% array compatible with buildClimateGrid. No directional spectrum is supplied;
% buildClimateGrid will fall back to gamma_default and use WaveDir for theta.
    A = parsed.All;
    N = numel(A.Hm0);
    if N == 0
        records = struct([]); return
    end
    template = struct( ...
        'Hm0', NaN, 'Te', NaN, 'Tp', NaN, 'WaveDir', NaN, ...
        'WindSpeed', NaN, 'WindDir', NaN, ...
        'Hm0WS', 0, 'TpfWS', 0, 'TeWS', 0, ...
        'Hm0SW', 0, 'TpfSW', 0, 'TeSW', 0, ...
        'f', [], 'theta', [], 'E', []);
    records = repmat(template, N, 1);
    for k = 1:N
        records(k).Hm0       = A.Hm0(k);
        records(k).Te        = A.Te(k);
        records(k).Tp        = A.Tp(k);
        records(k).WaveDir   = safeIdx(A, 'WaveDir',   k, NaN);
        records(k).WindSpeed = safeIdx(A, 'WindSpeed', k, NaN);
        records(k).WindDir   = safeIdx(A, 'WindDir',   k, NaN);
        records(k).Hm0WS     = safeIdx(A, 'Hm0WS', k, 0);
        records(k).TpfWS     = safeIdx(A, 'TpfWS', k, 0);
        records(k).TeWS      = safeIdx(A, 'TeWS',  k, 0);
        records(k).Hm0SW     = safeIdx(A, 'Hm0SW', k, 0);
        records(k).TpfSW     = safeIdx(A, 'TpfSW', k, 0);
        records(k).TeSW      = safeIdx(A, 'TeSW',  k, 0);
    end
end

function v = safeIdx(S, field, k, default)
    if isfield(S, field) && numel(S.(field)) >= k && isfinite(S.(field)(k))
        v = S.(field)(k);
    else
        v = default;
    end
end
