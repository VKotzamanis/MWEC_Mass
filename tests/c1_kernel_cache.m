function [geo, inner] = c1_kernel_cache(t, z_range)
%C1_KERNEL_CACHE  C1 outer NURBS (F1) and inner sets (F2) for tests, built once per kernel version.
%
%   geo = c1_kernel_cache()
%   [geo, inner] = c1_kernel_cache(t, z_range)
%
%   geo = mwecmass.solid.outer_nurbs(model) of Input/C1.ms2; inner = mwecmass.solid.offset_surface(
%   model, [], geo, t, z_range, struct('t_min', t)). Results are saved in one file per key under
%   /home/user/geomcache (outside the repository; MWEC_GEOMCACHE overrides the folder) and loaded on
%   every later call with the same key: the SHA-256 of Input/C1.ms2, every file in
%   src/+mwecmass/+solid, MS2Parser.m (F1 and F2 evaluate the deck through it) and this file. Any
%   change to one of them gives a new key, so the products are rebuilt once. An inner set is looked
%   up by t and z_range compared bitwise and added to the file the first time it is asked for.
%   The caller puts src on the path (and the Octave shims first, as the tests do).

root = fileparts(fileparts(mfilename('fullpath')));
folder = getenv('MWEC_GEOMCACHE');
if isempty(folder)
    folder = '/home/user/geomcache';
end
key = kernel_key(root);
file = fullfile(folder, ['c1_' key '.mat']);
model = [];
if exist(file, 'file')
    S = load(file);
    geo = S.geo;
    sets = S.sets;
    fprintf('c1_kernel_cache: loaded %s\n', file);
else
    model = parse_c1(root);
    tic;
    geo = mwecmass.solid.outer_nurbs(model);
    sets = struct('t', {}, 'z_range', {}, 'inner', {});
    fprintf('c1_kernel_cache: built F1 of C1 in %.1f s (key %s)\n', toc, key);
    store(folder, file, geo, sets);
end
if nargin == 0
    inner = [];
    return
end
k = find(arrayfun(@(s) isequal(s.t, t) && isequal(s.z_range, z_range), sets), 1);
if ~isempty(k)
    inner = sets(k).inner;
    return
end
if isempty(model)
    model = parse_c1(root);
end
tic;
inner = mwecmass.solid.offset_surface(model, [], geo, t, z_range, struct('t_min', t));
fprintf('c1_kernel_cache: built F2 of C1 at t = %.17g over %s in %.1f s\n', t, mat2str(z_range, 17), toc);
sets(end + 1) = struct('t', t, 'z_range', z_range, 'inner', inner);
store(folder, file, geo, sets);
end

function model = parse_c1(root)
evalc('model = mwecmass.geometry.MS2Parser.parse(fullfile(root, ''Input'', ''C1.ms2''));');
end

function store(folder, file, geo, sets)
% written to a temporary name in the same folder and renamed, so a stopped run leaves no partial file
if ~exist(folder, 'dir')
    mkdir(folder);
end
tmp = [file '.part'];
save(tmp, 'geo', 'sets', '-v7');
movefile(tmp, file, 'f');
end

function key = kernel_key(root)
files = [{fullfile(root, 'Input', 'C1.ms2'), fullfile(root, 'src', '+mwecmass', '+geometry', 'MS2Parser.m'), ...
          [mfilename('fullpath') '.m']}, sort(list_files(fullfile(root, 'src', '+mwecmass', '+solid')))];
parts = cell(1, numel(files));
for k = 1:numel(files)
    fid = fopen(files{k}, 'r');
    bytes = fread(fid, Inf, 'uint8=>uint8')';
    fclose(fid);
    parts{k} = [files{k}(numel(root) + 2:end) ':' sha256(bytes) ';'];
end
key = sha256(uint8([parts{:}]));
end

function files = list_files(folder)
d = dir(folder);
files = {};
for k = 1:numel(d)
    if any(strcmp(d(k).name, {'.', '..'}))
        continue
    end
    f = fullfile(folder, d(k).name);
    if d(k).isdir
        files = [files, list_files(f)]; %#ok<AGROW>
    else
        files{end + 1} = f; %#ok<AGROW>
    end
end
end

function h = sha256(bytes)
if exist('OCTAVE_VERSION', 'builtin')
    h = hash('sha256', char(bytes));
else
    md = java.security.MessageDigest.getInstance('SHA-256');
    v = typecast(md.digest(bytes), 'uint8');
    h = lower(reshape(dec2hex(v, 2)', 1, []));
end
end
