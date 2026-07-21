function A = hp2d_neal_cell_area(h, River_Width, River_Depth, Resolution)
%HP2D_NEAL_CELL_AREA Active water-surface area for Neal-style subgrid cells.
%   h is the representative water depth above channel bed [m].

A = Resolution.^2 .* ones(size(h), 'like', h);

has_channel = River_Width > 0 & River_Depth > 0;
inbank = has_channel & h <= River_Depth;

A(inbank) = River_Width(inbank) .* Resolution;
A(~isfinite(A)) = 0;

end
