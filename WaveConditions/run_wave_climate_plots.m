% run_wave_climate_plots.m
% Driver: render the station-only wave-climate figures (W1..W4) for every
% <station>_climate_grid.mat found under WIS_Output_WAM, and print a
% cross-station resource table.
%
% These figures describe the wave resource alone.  Nothing here needs a body
% file, a WAMIT cache, or a tuning run — point it at the climate grids and go.
% For the WEC-coupled figures (T1..T8) use MWEC_Tuning.m instead.
%
% Usage
%   1) Interactive: run the script; it picks up every grid in WIS_Output_WAM.
%   2) Scripted: set the CONFIG variables below before running.
%
% Outputs
%   <fig_dir>/fig_w1_scatter_<station>.{pdf,png,fig}
%   <fig_dir>/fig_w2_spectrum_<station>.{pdf,png,fig}
%   <fig_dir>/fig_w3_energy_<station>.{pdf,png,fig}
%   <fig_dir>/fig_w6_surface_<station>.{pdf,png,fig}
%   <fig_dir>/fig_w4_stations.{pdf,png,fig}
%   <fig_dir>/fig_w5_histograms.{pdf,png,fig}
%
% ----- CONFIG (override before running to skip the defaults) -----------------
% climate_dir = 'WIS_Output_WAM';    % where the *_climate_grid.mat files live
% fig_dir     = 'WaveConditions_Results/WaveClimate';
% stations    = {'ST63044','ST84040'};   % subset; default is all found
% formats     = {'pdf','png'};           % default {'pdf','png','fig'}
% sigma_factor = 0.75;   % W6 smoothing width, in climate-grid bin widths
% station_labels = struct('ST63044', 'North Atlantic, NH', ...
%                         'ST84040', 'Santa Barbara, CA');   % W5 column headers
% -----------------------------------------------------------------------------

this_dir = fileparts(mfilename('fullpath'));
if isempty(this_dir), this_dir = pwd; end
addpath(this_dir);

if ~exist('climate_dir', 'var') || isempty(climate_dir)
    climate_dir = fullfile(this_dir, 'WIS_Output_WAM');
end
if ~exist('fig_dir', 'var') || isempty(fig_dir)
    fig_dir = fullfile(this_dir, 'WaveConditions_Results', 'WaveClimate');
end
if ~isfolder(fig_dir), mkdir(fig_dir); end

cfg = struct('fig_dir', fig_dir);
if exist('formats', 'var') && ~isempty(formats)
    cfg.plot.formats = formats;
end
if exist('station_labels', 'var') && ~isempty(station_labels)
    cfg.station_labels = station_labels;
end

surf_opts = struct();
if exist('sigma_factor', 'var') && ~isempty(sigma_factor)
    surf_opts.sigma_factor = sigma_factor;
end

bar_line = repmat('=', 1, 78);
fprintf('\n%s\n  Wave-climate figures (station data only)\n%s\n', bar_line, bar_line);
fprintf('  Climate grids : %s\n  Figures       : %s\n\n', climate_dir, fig_dir);

listing = dir(fullfile(climate_dir, '*_climate_grid.mat'));
if isempty(listing)
    error('run_wave_climate_plots:noGrids', ...
          'No *_climate_grid.mat found in %s.  Build one with run_buildClimateGrid.m.', climate_dir);
end

if exist('stations', 'var') && ~isempty(stations)
    wanted = cellfun(@(x) [x '_climate_grid.mat'], cellstr(stations), 'UniformOutput', false);
    listing = listing(ismember({listing.name}, wanted));
    if isempty(listing)
        error('run_wave_climate_plots:noMatch', ...
              'None of the requested stations were found in %s.', climate_dir);
    end
end

grids = {};  rows = cell(0, 11);
for i = 1:numel(listing)
    path_i = fullfile(listing(i).folder, listing(i).name);
    L = load(path_i);
    if ~isfield(L, 'climateGrid')
        warning('run_wave_climate_plots:schema', ...
                'Skipping %s — no ''climateGrid'' variable.', listing(i).name);
        continue
    end
    cg = L.climateGrid;
    [sid, region, depth] = MWEC_WaveClimate_Plots.ident(cg);

    fprintf('  [%d/%d] %s (%s) ... ', i, numel(listing), sid, region);
    k = MWEC_WaveClimate_Plots.derive(cg);

    MWEC_WaveClimate_Plots.fig_w1_scatter(cg, cfg);
    MWEC_WaveClimate_Plots.fig_w2_spectrum(cg, cfg);
    MWEC_WaveClimate_Plots.fig_w3_energy(cg, cfg);
    MWEC_WaveClimate_Plots.fig_w6_energy_surface(cg, cfg, surf_opts);
    fprintf('W1 W2 W3 W6 done\n');

    q = MWEC_WaveClimate_Plots.energy_density(cg, surf_opts);
    fprintf('          surface volume %.4f kW/m vs sum(p*J) %.4f kW/m  (%.2f%% apart)\n', ...
            q.volume, q.P_wave, 100*abs(q.volume - q.P_wave)/max(q.P_wave, eps));

    if k.m0_err_pct > 1
        warning('run_wave_climate_plots:m0', ...
                '%s: m0(S_ew)=%.5f m^2 vs Hs budget %.5f m^2 (%.2f%% apart).', ...
                sid, k.IEC.m0, k.m0_target, k.m0_err_pct);
    end

    grids{end+1} = cg; %#ok<SAGROW>
    rows(end+1, :) = {sid, region, depth, k.IEC.Hm0, k.IEC.Te, k.IEC.Tp, ...
                      k.IEC.eps_bw, k.P_wave/1000, k.band.T_L, k.band.T_H, ...
                      k.n_cells_occupied}; %#ok<SAGROW>
end

if ~isempty(grids)
    fprintf('\n  Marginal histograms (all stations) ... ');
    MWEC_WaveClimate_Plots.fig_w5_histograms(grids, cfg);
    fprintf('W5 done\n');
end
if numel(grids) > 1
    fprintf('  Cross-station comparison ... ');
    MWEC_WaveClimate_Plots.fig_w4_stations(grids, cfg);
    fprintf('W4 done\n');
end

% --- resource table ----------------------------------------------------------
fprintf('\n%s\n  Resource summary (derived from the climate grids)\n%s\n', ...
        repmat('-', 1, 78), repmat('-', 1, 78));
fprintf('  %-9s %-14s %6s %7s %6s %6s %6s %9s %14s %6s\n', ...
        'station', 'region', 'h[m]', 'Hm0[m]', 'Te[s]', 'Tp[s]', 'eps', 'P[kW/m]', '90% band[s]', 'cells');
for i = 1:size(rows, 1)
    fprintf('  %-9s %-14s %6.0f %7.2f %6.2f %6.2f %6.3f %9.2f %6.1f-%-7.1f %6d\n', ...
            rows{i,1}, rows{i,2}, rows{i,3}, rows{i,4}, rows{i,5}, rows{i,6}, ...
            rows{i,7}, rows{i,8}, rows{i,9}, rows{i,10}, rows{i,11});
end

summary_path = fullfile(fig_dir, 'WaveClimate_summary.mat');
wave_climate_summary = cell2table(rows, 'VariableNames', ...
    {'station','region','depth_m','Hm0_m','Te_s','Tp_s','eps_bw','P_wave_kWm', ...
     'band_T_L_s','band_T_H_s','cells_occupied'});
save(summary_path, 'wave_climate_summary');

fprintf('\n  Summary table : %s\n', summary_path);
fprintf('  Figures       : %s\n%s\n\n', fig_dir, bar_line);
