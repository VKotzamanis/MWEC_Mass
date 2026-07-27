% _run_gm_test.m - headless harness for WEC_GM.m
% Runs WEC_GM, then saves every open figure to _run_artifacts/ as PNG.

this_dir = fileparts(mfilename('fullpath'));
artdir   = fullfile(this_dir, '_run_artifacts');
if ~exist(artdir, 'dir'); mkdir(artdir); end
setenv('WEC_GM_ARTDIR', artdir);   % survives WEC_GM.m's `clear`
setenv('WEC_GM_THIS', this_dir);

% Force figures off-screen so headless run is OK
set(groot, 'DefaultFigureVisible', 'off');

cd(this_dir);

diary(fullfile(artdir, 'gm_log.txt'));
diary on;

try
    run(fullfile(this_dir, 'WEC_GM.m'));
    ok = true;
    err_msg = '';
catch ME
    ok = false;
    err_msg = sprintf('%s\n%s', ME.message, getReport(ME, 'extended'));
    fprintf(2, '\n[HARNESS] WEC_GM.m FAILED:\n%s\n', err_msg);
end

artdir   = getenv('WEC_GM_ARTDIR');
this_dir = getenv('WEC_GM_THIS');

% Save all open figures
figs = findall(groot, 'Type', 'figure');
fprintf('\n[HARNESS] %d figure(s) open. Saving to %s\n', numel(figs), artdir);
for k = 1:numel(figs)
    f = figs(k);
    name = get(f, 'Name');
    if isempty(name), name = sprintf('fig%d', k); end
    safe = regexprep(name, '[^A-Za-z0-9_]', '_');
    out  = fullfile(artdir, sprintf('%02d_%s.png', k, safe));
    try
        exportgraphics(f, out, 'Resolution', 150);
        fprintf('  saved %s\n', out);
    catch
        saveas(f, out);
        fprintf('  saved (legacy) %s\n', out);
    end
end

diary off;

if ~ok
    fprintf(2, '[HARNESS] exit non-zero\n');
    exit(1);
end
fprintf('[HARNESS] done.\n');
exit(0);
