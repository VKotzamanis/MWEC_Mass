function style_line(h, style, role)
%STYLE_LINE  Apply the shared width for a curve, boundary, or reference line, and the shared
%marker size to any handle that has a marker.
% Inputs: h is a line/patch/quiver/etc. handle or handle array; style is out.style; role is curve,
%         boundary, or reference.
% Outputs: none; LineWidth (pt) is updated on every handle. MarkerSize (pt, from
%          style.marker.size) is additionally applied to any handle that both exposes a Marker
%          property and has it set to something other than 'none' -- a plotted-point marker, not
%          a uniformly-sized trace/emphasis marker pair, which some diagnostic animations set
%          directly and deliberately at two different sizes (see
%          validation/diagnostics/stage_animations.m, not routed through this helper).
    valid_roles = {'curve', 'boundary', 'reference'};
    if ~any(strcmp(role, valid_roles))
        error('mwecmass:figures:UnknownLineRole', ...
              'role must be curve, boundary, or reference.');
    end
    set(h, 'LineWidth', style.line_width.(role));
    for k = 1:numel(h)
        hk = h(k);
        if isprop(hk, 'Marker') && ~strcmp(get(hk, 'Marker'), 'none')
            set(hk, 'MarkerSize', style.marker.size);
        end
    end
end
