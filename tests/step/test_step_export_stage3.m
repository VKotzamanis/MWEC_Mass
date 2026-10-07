function test_step_export_stage3()
%TEST_STEP_EXPORT_STAGE3  export_stage3 (F11) on the stand-in bodies (box, cylinder), both modes.
%   Per file: gmsh/OpenCASCADE import by tests/step_check.py (METRE, solid count, no open mesh edge on a
%   solid, a sheet's boundary as the only open edges), imported volumes against the kernel volumes
%   (S6 of the stand-in F6) and against sti_closed_form, which is the independent oracle. The box is
%   planar, so its volumes are asserted to the bound of test_step_cube (OCC rounding only); the
%   cylinder's rational faces are compared and printed. The offset distance of the closed form is
%   d = t + 0.01 t_min / 2 (contract section 0).

root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
setup(root);
% {fixture, mode, Stage-2 strip densities, t, z_ballast, solid_modules, vs, label}
cases = {
    'box', 'modular_precast', [NaN; 600; 500], [0.08; 0.08; 0.08], -2.3, [], 0, 'box precast, ballast in module 1'
    'box', 'modular_precast', [NaN; 600; 500], [0.08; 0.08; NaN], -2.3, 3, 0, 'box precast, solid wall module on top'
    'box', 'modular_precast', [NaN; 600; 500], NaN(3, 1), -2.5, 1:3, 0, 'box precast, every module solid'
    'box', 'thin_shell', [NaN; 500; 300], 0.0254 * ones(3, 1), -2.3, [], 0.5, 'box thin shell, ballast in module 1'
    'box', 'thin_shell', [NaN; 500; 300], 0.0254 * ones(3, 1), -1.5, [], 0.5, 'box thin shell, ballast at a module edge'
    'box', 'thin_shell', [NaN; 500; 300], 0.0254 * ones(3, 1), -2.6, [], 0.5, 'box thin shell, no ballast'
    'cylinder', 'modular_precast', [NaN; 250; 200; 150], 0.1 * ones(4, 1), -2.5, [], 0.5, 'cylinder precast, ballast in module 1'
    'cylinder', 'thin_shell', [NaN; 400; 250; 200], 0.0254 * ones(4, 1), -2.85, [], 0.5, 'cylinder thin shell, ballast in module 1'
    };
worst_box = 0;
for c = 1:size(cases, 1)
    [name, mode, rho2, t, zb, solid, vs, label] = cases{c, :};
    config = sti_config(name);
    config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
    [~, stage2] = sti_stage2(config, vs, rho2);
    design = struct('mode', mode, 'edges', config.strip_edges, 'vs', vs, 't', t, 'z_ballast', zb, ...
        'solid_modules', solid);
    precast = strcmp(mode, 'modular_precast');
    if precast
        rho = struct('uhpc', 2500, 'air', 1.2);
    else
        rho = struct('ballast', 7500, 'shell', 7850, 'air', 1.2);
    end
    t_min = t(find(isfinite(t), 1));
    if isempty(t_min)
        t_min = 0.08;
    end
    r = sti_realised(config, design, rho, stage2, struct('t_min', t_min));
    dir_out = tempname();
    cleanup = onCleanup(@() rmdir_if(dir_out)); %#ok<NASGU>
    files = mwecmass.output.step.export_stage3(r, dir_out);

    N = numel(r.modules);
    hull = r.hull_name;
    if precast
        want = [arrayfun(@(i) sprintf('%s_UHPC_module_%d', hull, i), 1:N, 'UniformOutput', false), ...
            {[hull '_UHPC_all']}];
    else
        want = {};
        if zb > config.hull_z_min
            want{end + 1} = [hull '_STEEL_ballast'];
        end
        want{end + 1} = [hull '_STEEL_shell'];
        if numel(want) == 2
            want{end + 1} = [hull '_STEEL_all'];
        end
    end
    check(isequal({files.name}, want), '%s: file names %s', label, strjoin({files.name}, ', '));
    for k = 1:numel(files)
        check(exist(files(k).path, 'file') == 2 && strcmp(files(k).path, fullfile(dir_out, 'step', [files(k).name '.step'])), ...
            '%s: path of %s', label, files(k).name);
    end
    listed = dir(fullfile(dir_out, 'step', '*.step'));
    check(numel(listed) == numel(files), '%s: %d files on disk, %d returned', label, numel(listed), numel(files));
    check(isempty(dir(fullfile(dir_out, '*.step'))), '%s: a file was written outside the step folder', label);

    fx = sti_closed_form('fixture', config.hull_solid);
    d = t + 0.01 * t_min / 2;
    d(solid) = NaN;
    reg = sti_closed_form('regions', fx, design, d);
    parts = cell(1, numel(files));
    for k = 1:numel(files)
        parts{k} = stp_check(files(k).path);
        p = parts{k};
        check(strcmp(p.declared_length_unit, 'METRE'), '%s: %s units %s', label, files(k).name, p.declared_length_unit);
        check(p.open_edges_volumes == 0, '%s: %s has open edges on a solid', label, files(k).name);
        text = fileread(files(k).path);
        for q = 1:numel(files(k).bodies)
            tag = sprintf('PRODUCT(''%s''', files(k).bodies{q});
            check(~isempty(strfind(text, tag)), '%s: body %s unnamed in %s', ...
                label, files(k).bodies{q}, files(k).name);
        end
        check(numel(regexp(text, '=PRODUCT\(', 'match')) == numel(files(k).bodies), '%s: product count in %s', label, files(k).name);
    end

    if precast
        for i = 1:N
            worst_box = max(worst_box, compare(label, files(i).name, parts{i}, 1, r.modules(i).V_uhpc, reg.uhpc.V(i), strcmp(name, 'box')));
        end
        worst_box = max(worst_box, compare(label, files(N + 1).name, parts{N + 1}, 1, sum([r.modules.V_uhpc]), sum(reg.uhpc.V), strcmp(name, 'box')));
        check(isequal(files(N + 1).bodies, {want{N + 1}}), '%s: fused file bodies', label);
        for i = 1:N
            check(parts{i}.open_edges_sheets == 0 && parts{N + 1}.open_edges_sheets == 0, '%s: free sheet in a solid file', label);
        end
    else
        Vb = sum([r.modules.V_ballast]);
        Vb0 = sum(reg.ballast.V);
        kshell = find(strcmp({files.name}, [hull '_STEEL_shell']));
        sheet = parts{kshell};
        % with ballast the sheet ends at the junction curve; without, it is the closed hull surface
        check(sheet.n_volumes == 0 && (sheet.open_edges_sheets > 0) == (zb > config.hull_z_min), ...
            '%s: sheet boundary (%d open sheet edges)', label, sheet.open_edges_sheets);
        if zb > config.hull_z_min
            kb = 1;
            worst_box = max(worst_box, compare(label, files(kb).name, parts{kb}, 1, Vb, Vb0, strcmp(name, 'box')));
            check(parts{kb}.open_edges_sheets == 0, '%s: ballast file holds a sheet', label);
            ka = find(strcmp({files.name}, [hull '_STEEL_all']));
            both = parts{ka};
            worst_box = max(worst_box, compare(label, files(ka).name, both, 1, Vb, Vb0, strcmp(name, 'box')));
            check(both.open_edges_sheets > 0 && both.n_surfaces == parts{kb}.n_surfaces + sheet.n_surfaces, ...
                '%s: combined file is not the two bodies', label);
            check(isequal(files(ka).bodies, {files(kb).bodies{1}, files(kshell).bodies{1}}), '%s: combined bodies', label);
            % the junction curve is written once: edges of the combined file = ballast + sheet - junction
            brep = r.body.brep;
            eb = edges_of(brep, brep.bodies(1).shells);
            es = edges_of(brep, brep.bodies(2).shells);
            junction = intersect(eb, es);
            n_all = numel(regexp(fileread(files(ka).path), '=EDGE_CURVE\(', 'match'));
            check(~isempty(junction) && n_all == numel(union(eb, es)), '%s: %d EDGE_CURVE in the combined file, %d expected', ...
                label, n_all, numel(union(eb, es)));
            fprintf('  junction edges shared by ballast and sheet: %d; EDGE_CURVE in ballast, shell, combined file: %d, %d, %d\n', ...
                numel(junction), count_edges(files(kb).path), count_edges(files(kshell).path), n_all);
        else
            check(numel(files) == 1, '%s: no ballast, one file', label);
        end
    end
    fprintf('%-52s files %s\n', label, strjoin({files.name}, ', '));
    clear cleanup
    rmdir_if(dir_out);
end
fprintf('largest relative |OCC - kernel| and |OCC - closed form| on the box files: %.3e\n', worst_box);

% stale files of the mode's naming are removed, other files stay
d = tempname();
cleanup = onCleanup(@() rmdir_if(d));
mkdir(fullfile(d, 'step'));
config = sti_config('box');
config.hull_solid = mwecmass.solid.outer_nurbs(config.ms2_model);
[~, stage2] = sti_stage2(config, 0, [NaN; 600; 500]);
design = struct('mode', 'modular_precast', 'edges', config.strip_edges, 'vs', 0, 't', 0.08 * ones(3, 1), ...
    'z_ballast', -2.3, 'solid_modules', []);
r = sti_realised(config, design, struct('uhpc', 2500, 'air', 1.2), stage2, struct('t_min', 0.08));
for f = {'box_UHPC_module_7.step', 'box_UHPC_all.step', 'box_STEEL_all.step', 'notes.txt'}
    fid = fopen(fullfile(d, 'step', f{1}), 'w');
    fclose(fid);
end
mwecmass.output.step.export_stage3(r, d);
check(exist(fullfile(d, 'step', 'box_UHPC_module_7.step'), 'file') == 0 && ...
    exist(fullfile(d, 'step', 'notes.txt'), 'file') == 2 && exist(fullfile(d, 'step', 'box_STEEL_all.step'), 'file') == 2, ...
    'stale files');
sz = dir(fullfile(d, 'step', 'box_UHPC_all.step'));
check(sz.bytes > 0, 'box_UHPC_all.step was not rewritten');

% a combined export whose sheet boundary is not the ballast's top outline errors, and so do bad inputs
design = struct('mode', 'thin_shell', 'edges', config.strip_edges, 'vs', 0.5, 't', 0.0254 * ones(3, 1), ...
    'z_ballast', -2.3, 'solid_modules', []);
[~, stage2] = sti_stage2(config, 0.5, [NaN; 500; 300]);
r = sti_realised(config, design, struct('ballast', 7500, 'shell', 7850, 'air', 1.2), stage2, struct('t_min', 0.0254));
bad = r;
sh = bad.body.brep.bodies(2).shells{1};
bad.body.brep.bodies(2).shells{1} = sh(1:end - 1);
expect_error(@() mwecmass.output.step.export_stage3(bad, d), 'mwecmass:step:JunctionNotShared');
% loops given as a column cell must give the same junction check and the same files
col = r;
for i = 1:numel(col.body.brep.faces)
    col.body.brep.faces(i).loops = col.body.brep.faces(i).loops(:);
end
fr = mwecmass.output.step.export_stage3(r, d);
text_row = arrayfun(@(f) data_section(f.path), fr, 'UniformOutput', false);
fc = mwecmass.output.step.export_stage3(col, d);
check(isequal({fr.name}, {fc.name}), 'column loops: file names');
for k = 1:numel(fc)
    check(strcmp(text_row{k}, data_section(fc(k).path)), 'column loops: %s differs', fc(k).name);
end
check(isequal(edges_of(r.body.brep, r.body.brep.bodies(1).shells), edges_of(col.body.brep, col.body.brep.bodies(1).shells)), ...
    'column loops: edges of the ballast solid');
bad = r;
bad.mode = 'preliminary';
expect_error(@() mwecmass.output.step.export_stage3(bad, d), 'mwecmass:step:BadMode');
bad = r;
bad.body = [];
expect_error(@() mwecmass.output.step.export_stage3(bad, d), 'mwecmass:step:NoBody');
fprintf('error paths: JunctionNotShared, BadMode, NoBody\n');
end

function worst = compare(label, file, res, n, V_kernel, V_closed, asserted)
% worst: relative OCC deviations of a planar (asserted) file, else 0
occ = res.occ_volumes(1);
check(res.n_volumes == n, '%s: %s imported %d solids, expected %d', label, file, res.n_volumes, n);
e1 = abs(occ - V_kernel) / V_kernel;
e2 = abs(occ - V_closed) / V_closed;
fprintf('  %-34s OCC %.12g  kernel %.12g  closed form %.12g  rel(OCC-kernel) %.2e  rel(OCC-closed) %.2e\n', ...
    file, occ, V_kernel, V_closed, e1, e2);
worst = 0;
if asserted
    % planar faces are exact: only OCC rounding remains (bound of test_step_cube)
    worst = max(e1, e2);
    check(e1 < 1e-9 && e2 < 1e-9, '%s: %s volume differs from the box closed form', label, file);
end
end

function s = data_section(file)
% the header carries a time stamp; the DATA section is what the loops decide
s = fileread(file);
s = s(strfind(s, 'DATA;'):end);
end

function n = count_edges(file)
n = numel(regexp(fileread(file), '=EDGE_CURVE\(', 'match'));
end

function e = edges_of(brep, shells)
e = [];
for q = 1:numel(shells)
    for i = abs(shells{q}(:)')
        for j = 1:numel(brep.faces(i).loops)
            e = union(e, abs(brep.faces(i).loops{j}(:)'));
        end
    end
end
end

function expect_error(fn, id)
try
    fn();
catch err
    if strcmp(err.identifier, id)
        return
    end
    error('test_step_export_stage3:fail', 'expected %s, got %s: %s', id, err.identifier, err.message);
end
error('test_step_export_stage3:fail', 'expected error %s, none raised', id);
end

function rmdir_if(d)
if exist(d, 'dir')
    rmdir(d, 's');
end
end

function setup(root)
if exist('OCTAVE_VERSION', 'builtin')
    warning('off', 'Octave:shadowed-function');
    addpath(fullfile(root, 'tests', 'octave_shims'));
end
addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'tests', 'step'));
addpath(fullfile(root, 'tests', 'standins'), '-end');
addpath(fullfile(root, 'tests', 'standins', 'fixtures'), '-end');
end

function check(cond, varargin)
if ~cond
    error('test_step_export_stage3:fail', varargin{:});
end
end
