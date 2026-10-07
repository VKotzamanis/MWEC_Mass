function test_step_units()
%TEST_STEP_UNITS  The file declares METRE and a unit-aware import returns metres.
%
%   A copy with the length unit relabelled MILLI shows that the check discriminates: the same
%   numbers import 1000 times smaller.

addpath(fileparts(mfilename('fullpath')));
size_xyz = [2 3 5];
b = stp_new();
[b, shell] = stp_box(b, [1 -2 0.5], [1 -2 0.5] + size_xyz);
b.bodies(1) = struct('name', 'box', 'kind', 'solid', 'shells', {{shell}});
file = [tempname() '.step'];
milli = [tempname() '.step'];
cleanup = onCleanup(@() cellfun(@delete, {file, milli}));
mwecmass.output.step.write_step(b, file);

text = fileread(file);
if isempty(strfind(text, 'SI_UNIT($,.METRE.)')) || ~isempty(strfind(text, '.MILLI.'))
    error('test_step_units: the file must declare SI_UNIT($,.METRE.) with no prefix');
end
fid = fopen(milli, 'w');
fwrite(fid, strrep(text, 'SI_UNIT($,.METRE.)', 'SI_UNIT(.MILLI.,.METRE.)'));
fclose(fid);

r = stp_check(file);
rm = stp_check(milli);
lo = [1 -2 0.5];
expected = [lo lo + size_xyz];
fprintf('units: declared %s, OCC bbox %s\n', r.declared_length_unit, mat2str(r.bbox_occ(:)', 9));
fprintf('units: mesh bbox %s, max |mesh bbox - expected| %.3g\n', mat2str(r.bbox_mesh(:)', 9), ...
    max(abs(r.bbox_mesh(:)' - expected)));
fprintf('units: relabelled MILLI: declared %s, OCC bbox %s\n', rm.declared_length_unit, mat2str(rm.bbox_occ(:)', 9));
if ~strcmp(r.declared_length_unit, 'METRE')
    error('test_step_units: declared length unit is %s, not METRE', r.declared_length_unit);
end
% OpenCASCADE widens its bounding box by the shape tolerance, the stated file uncertainty
% (1e-7 m); 4 eps of the largest coordinate covers the rounding of that sum
unc = 1e-7 + 4 * eps(max(abs(expected)));
if max(abs(r.bbox_occ(:)' - expected)) > unc
    error('test_step_units: imported bounding box is not in metres');
end
% 1e-12: mesh nodes at the corners are the file's CARTESIAN_POINT values
if max(abs(r.bbox_mesh(:)' - expected)) > 1e-12
    error('test_step_units: mesh bounding box differs from the box in metres');
end
if max(abs(rm.bbox_mesh(:)' - expected / 1000)) > 1e-12
    error('test_step_units: the MILLI copy did not import 1000 times smaller');
end
end
