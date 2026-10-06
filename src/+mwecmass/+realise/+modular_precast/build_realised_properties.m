function realised = build_realised_properties(final_props, cstr, config)
%BUILD_REALISED_PROPERTIES Rebuild hydrostatics and dynamics from realised UHPC geometry.
% cstr supplies feasible as-built mass, CG, inertia, and hydrostatics; config supplies interpolation tables.
% Added-mass coefficients and coupled periods are recomputed. Infeasible input follows the shared fallback.
% See docs/METHODS_ENGINE.md#realise-property-rebuild
    % cstr originates from mwecmass.realise.modular_precast.solve() via
    % extract_strip_geometry(), so it carries all steel_data fields
    % required by Shell_Offset's version PLUS per-strip realisation
    % arrays consumed by build_realised_properties.
    realised = mwecmass.realise.build_realised_properties(final_props, cstr, config);
    realised.fill_method = 'uhpc_fill';
    if isfield(realised, 'realised_strips') && isstruct(realised.realised_strips)
        realised.realised_strips.fill_method = 'uhpc_fill';
    end
end
