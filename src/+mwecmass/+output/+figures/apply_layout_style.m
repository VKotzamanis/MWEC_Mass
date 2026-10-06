function apply_layout_style(tl, style)
%APPLY_LAYOUT_STYLE Apply shared tile spacing, padding, and title/subtitle styling to a tiledlayout.
% tl is a tiledlayout handle; style is out.style. Updates TileSpacing,
% Padding, and tl.Title and tl.Subtitle text objects in place.
    set(tl, 'TileSpacing', style.layout.tile_spacing, 'Padding', style.layout.padding);
    mwecmass.output.figures.style_text(tl.Title, style, 'title');
    mwecmass.output.figures.style_text(tl.Subtitle, style, 'label');
end
