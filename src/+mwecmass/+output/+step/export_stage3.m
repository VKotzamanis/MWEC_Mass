function files = export_stage3(realised, out_dir)
%EXPORT_STAGE3  Write the STEP files of a realised Stage-3 design (contract F11).
%
%   files = mwecmass.output.step.export_stage3(realised, out_dir)
%
%   realised: results.stage3 with its body (S4); the body's B-rep carries the bodies of the mode.
%   out_dir: the type folder, mwecmass.output.output_dir(type); the files go into its subfolder
%   'step' (created when missing). Files already there that follow this mode's naming are
%   removed first, so a design with fewer modules or no ballast leaves no stale file.
%
%   Modular precast: <hull>_UHPC_module_<i>.step, one solid per module, and <hull>_UHPC_all.step,
%   one solid whose outer shell is the hull surface and whose void shells are the cavities
%   (S4 shells.all_outer, shells.all_voids), so the fused solid has no joint faces.
%   Thin shell: <hull>_STEEL_ballast.step (solid), <hull>_STEEL_shell.step (sheet from z_ballast
%   up) and, when both bodies exist, <hull>_STEEL_all.step holding both. The two bodies share the
%   junction curve at z_ballast by edge index, so the file writes it once.
%
%   files: struct array with name (file stem), path and bodies (cellstr of body names in it).

if ~isstruct(realised) || ~isfield(realised, 'body') || isempty(realised.body) || ...
        ~isfield(realised.body, 'brep')
    error('mwecmass:step:NoBody', 'export_stage3: realised has no body with a B-rep');
end
brep = realised.body.brep;
hull = realised.hull_name;

switch realised.mode
    case 'modular_precast'
        stale = {[hull '_UHPC_module_*.step'], [hull '_UHPC_all.step']};
        parts = brep.bodies;
        fused = struct('name', [hull '_UHPC_all'], 'kind', 'solid', ...
            'shells', {[{realised.body.shells.all_outer}, realised.body.shells.all_voids(:)']});
        extra = {fused};
    case 'thin_shell'
        stale = {[hull '_STEEL_ballast.step'], [hull '_STEEL_shell.step'], [hull '_STEEL_all.step']};
        parts = brep.bodies;
        extra = {};
        if numel(parts) == 2
            check_junction(brep, parts);
            extra = {parts};
        end
    otherwise
        error('mwecmass:step:BadMode', 'export_stage3: no STEP export for mode %s', realised.mode);
end

step_dir = fullfile(out_dir, 'step');
if ~exist(step_dir, 'dir')
    [ok, msg] = mkdir(step_dir);
    if ~ok
        error('mwecmass:step:io', 'export_stage3: cannot create %s: %s', step_dir, msg);
    end
end
for k = 1:numel(stale)
    old = dir(fullfile(step_dir, stale{k}));
    for q = 1:numel(old)
        delete(fullfile(step_dir, old(q).name));
    end
end

files = struct('name', {}, 'path', {}, 'bodies', {});
for k = 1:numel(parts)
    files(end + 1) = write_file(brep, parts(k), parts(k).name, step_dir); %#ok<AGROW>
end
if strcmp(realised.mode, 'modular_precast')
    files(end + 1) = write_file(brep, extra{1}, extra{1}.name, step_dir);
elseif ~isempty(extra)
    files(end + 1) = write_file(brep, extra{1}, [hull '_STEEL_all'], step_dir);
end
end

function f = write_file(brep, bodies, stem, step_dir)
brep.bodies = bodies;
path = fullfile(step_dir, [stem '.step']);
mwecmass.output.step.write_step(brep, path);
f = struct('name', stem, 'path', path, 'bodies', {{bodies.name}});
end

function check_junction(brep, parts)
% The sheet's free edges (used by one face of the sheet) must be exactly the edges the sheet has
% in common with the ballast solid: the junction curve at z_ballast, one edge index in both bodies.
kinds = {parts.kind};
solid = parts(strcmp(kinds, 'solid'));
sheet = parts(strcmp(kinds, 'sheet'));
if numel(solid) ~= 1 || numel(sheet) ~= 1
    error('mwecmass:step:JunctionNotShared', 'export_stage3: expected one ballast solid and one shell sheet');
end
sheet_edges = [];
for q = 1:numel(sheet.shells)
    sheet_edges = [sheet_edges, edge_list(brep, sheet.shells{q})]; %#ok<AGROW>
end
solid_edges = [];
for q = 1:numel(solid.shells)
    solid_edges = [solid_edges, edge_list(brep, solid.shells{q})]; %#ok<AGROW>
end
[u, ~, j] = unique(sheet_edges);
free = u(accumarray(j(:), 1)' == 1);
common = intersect(unique(sheet_edges), unique(solid_edges));
if isempty(free) || ~isequal(sort(free), sort(common))
    error('mwecmass:step:JunctionNotShared', ...
        'export_stage3: the sheet boundary is not the curve shared with the ballast solid');
end
end

function e = edge_list(brep, signed_faces)
e = [];
for i = abs(signed_faces(:)')
    for j = 1:numel(brep.faces(i).loops)
        e = [e, abs(brep.faces(i).loops{j}(:)')]; %#ok<AGROW>
    end
end
end
