function unc = declared_uncertainty(brep)
%DECLARED_UNCERTAINTY  Distance accuracy stated in the file, in metres (default 1e-7).

unc = 1e-7;
if isfield(brep, 'uncertainty') && ~isempty(brep.uncertainty)
    unc = brep.uncertainty;
    if ~isscalar(unc) || ~isfinite(unc) || unc <= 0
        error('mwecmass:step:invalid', 'write_step: uncertainty must be a positive finite scalar');
    end
end
end
