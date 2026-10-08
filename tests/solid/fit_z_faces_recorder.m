function dir_out = fit_z_faces_recorder()
%FIT_Z_FACES_RECORDER  Test-only recording stub of mwecmass.solid.fit_z_faces (helper, not a test).
%
%   d = fit_z_faces_recorder()
%
%   Writes, in a new temporary folder d, a package function +mwecmass/+solid/fit_z_faces.m that
%   stores its inputs in the global FIT_Z_FACES_RECORD (fields faces, general, opts) and raises
%   fit_z_faces_stub:reached. The caller adds d to the front of the path (a package function in an
%   earlier path entry wins, so the stub also shadows the real function after J1), calls F1 or F2,
%   expects that error, reads the record, and removes d with rmpath and rmdir.

dir_out = tempname();
pkg = fullfile(dir_out, '+mwecmass', '+solid');
mkdir(pkg);
src = {
    'function [faces, flat, rep] = fit_z_faces(model, cache, faces, general, opts) %#ok<INUSL,STOUT>'
    'global FIT_Z_FACES_RECORD'
    'FIT_Z_FACES_RECORD = struct(''faces'', {faces}, ''general'', {general}, ''opts'', {opts});'
    'error(''fit_z_faces_stub:reached'', ''recording stub of fit_z_faces reached'');'
    'end'};
fid = fopen(fullfile(pkg, 'fit_z_faces.m'), 'w');
fprintf(fid, '%s\n', src{:});
fclose(fid);
end
