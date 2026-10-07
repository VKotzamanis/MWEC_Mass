function d_close = void_closing_distance(model, cache, geo, z_range) %#ok<INUSL>
%VOID_CLOSING_DISTANCE  Stand-in of contract F2b: offset distance at which the fixture's void closes.
%
%   d_close = mwecmass.solid.void_closing_distance(model, cache, geo, z_range)
%
%   Closed form (contract section 3, Stand-in kit SK): cylinder min(R, H/2), box half its smallest
%   dimension; opposite offset layers of a prism meet at half the distance between the faces.
%   z_range must lie in the hull; the closed form does not depend on it. Requires geo.analytic,
%   else mwecmass:standin:NotAnalytic.

if ~isstruct(geo) || ~isfield(geo, 'analytic') || isempty(geo.analytic)
    error('mwecmass:standin:NotAnalytic', 'void_closing_distance stand-in: geo.analytic is empty (not a stand-in fixture)');
end
if exist('sti_closed_form', 'file') ~= 2
    addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'fixtures'), '-end');
end
zr = sort(z_range(:)');
if numel(zr) ~= 2 || zr(1) < geo.analytic.z(1) || zr(2) > geo.analytic.z(2)
    error('mwecmass:solid:ZOutside', 'void_closing_distance stand-in: z_range outside the hull');
end
d_close = sti_closed_form('d_close', geo.analytic);
end
