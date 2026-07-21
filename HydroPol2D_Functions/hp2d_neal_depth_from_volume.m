function h = hp2d_neal_depth_from_volume(V, River_Width, River_Depth, Resolution)
%HP2D_NEAL_DEPTH_FROM_VOLUME Invert Neal-style cell storage to depth [m].

V = max(V, 0);
h = V ./ (Resolution.^2);

has_channel = River_Width > 0 & River_Depth > 0;
if ~any(has_channel(:))
    h(~isfinite(h)) = 0;
    return;
end

V_bank = Resolution .* River_Width .* River_Depth;
inbank = has_channel & V <= V_bank;
overbank = has_channel & V > V_bank;

h(inbank) = V(inbank) ./ max(Resolution .* River_Width(inbank), eps);
h(overbank) = River_Depth(overbank) + ...
    (V(overbank) - V_bank(overbank)) ./ (Resolution.^2);

h(~isfinite(h)) = 0;
h = max(h, 0);

end
