function test_step_sheet()
%TEST_STEP_SHEET  Open sheet of the four bicubic side faces: surfaces only, no volume.

addpath(fileparts(mfilename('fullpath')));
b = stp_new();
[b, lateral] = stp_bulged_box(b);
b.bodies(1) = struct('name', 'open_sheet', 'kind', 'sheet', 'shells', {{lateral}});
file = [tempname() '.step'];
cleanup = onCleanup(@() delete(file));
mwecmass.output.step.write_step(b, file);
r = stp_check(file);

fprintf('sheet: volumes %d, surfaces %d, open mesh edges (the sheet boundary) %d\n', ...
    r.n_volumes, r.n_surfaces, r.open_edges);
if r.n_volumes ~= 0 || r.n_surfaces ~= 4
    error('test_step_sheet: expected 0 volumes and 4 surfaces');
end
end
