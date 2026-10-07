function hs = hydrostatics_at_draft(geo, vs, opts) %#ok<INUSD>
%HYDROSTATICS_AT_DRAFT  Stand-in of contract F7: S7 hydrostatics of a stand-in fixture in closed form.
%
%   hs = mwecmass.solid.hydrostatics_at_draft(geo, vs, opts)
%
%   geo.analytic must be set (stand-in F1), else mwecmass:standin:NotAnalytic. Numbers from
%   sti_closed_form 'hydrostatics' (vertical prism: V_sub = A (zw - z0), CB at mid-depth, I_wp about
%   the centre of flotation by the parallel-axis theorem, S_wet = A + perimeter * depth; source in
%   its help); the waterline loop is the slice_bspline_surface stand-in of geo.outer at z = -vs
%   (partial submersion only). Waterline at or above z_max: 'full'; at or below z_min: 'none'.

if ~isstruct(geo) || ~isfield(geo, 'analytic') || isempty(geo.analytic)
    error('mwecmass:standin:NotAnalytic', 'hydrostatics_at_draft stand-in: geo.analytic is empty (not a stand-in fixture)');
end
if exist('sti_closed_form', 'file') ~= 2
    addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'fixtures'), '-end');
end
hs = sti_closed_form('hydrostatics', geo.analytic, vs);
hs.waterline = [];
if strcmp(hs.submersion, 'partial')
    hs.waterline = mwecmass.solid.slice_bspline_surface(geo.outer, -vs);
end
end
