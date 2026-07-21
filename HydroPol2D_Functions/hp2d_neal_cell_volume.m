function V = hp2d_neal_cell_volume(h, River_Width, River_Depth, Resolution)
%HP2D_NEAL_CELL_VOLUME Cell storage volume for Neal-style subgrid cells [m3].
%   h is the representative water depth above channel bed [m].

h = max(h, 0);

V = Resolution.^2 .* h;

has_channel = River_Width > 0 & River_Depth > 0;
if ~any(has_channel(:))
    return;
end

inbank = has_channel & h <= River_Depth;
overbank = has_channel & h > River_Depth;

V(inbank) = Resolution .* River_Width(inbank) .* h(inbank);
V(overbank) = Resolution .* River_Width(overbank) .* River_Depth(overbank) + ...
    Resolution.^2 .* (h(overbank) - River_Depth(overbank));

V(~isfinite(V)) = 0;

end
