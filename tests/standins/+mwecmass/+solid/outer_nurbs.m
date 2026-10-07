function geo = outer_nurbs(model)
%OUTER_NURBS  Stand-in of contract F1: S1b of a stand-in fixture deck, with geo.analytic set.
%
%   geo = mwecmass.solid.outer_nurbs(model)
%
%   model is MS2Parser.parse of tests/standins/fixtures/cylinder.ms2 or box.ms2. The patches are
%   built from the closed form of the fixture (sti_closed_form 'patches'): the cylinder is the
%   PolyCurve2 of three Lines (degree 1, C0 knots 1/3 and 2/3 at the parser's parameter map)
%   revolved by the rational quadratic 90-degree arc (weights 1, sqrt(2)/2, 1), the box faces are
%   bilinear RuledSurf patches; source: contract F1 (conversions of Line, PolyCurve2, RevSurf,
%   RuledSurf and mirrors). Any other deck errors mwecmass:standin:NotAnalytic.

standin_path();
fx = sti_closed_form('fixture', model);
outer = sti_closed_form('patches', fx, 0);
geo = struct('hull_name', fx.name, 'outer', outer, 'z_range', fx.z, 'analytic', fx);
end

function standin_path()
if exist('sti_closed_form', 'file') ~= 2
    addpath(fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'fixtures'), '-end');
end
end
