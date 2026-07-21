function d_surface_mm = hp2d_neal_surface_depth(d_rep_mm, River_Depth)
%HP2D_NEAL_SURFACE_DEPTH Convert Neal representative depth to flood depth [mm].

d_surface_mm = d_rep_mm;

has_channel = River_Depth > 0;
d_surface_mm(has_channel) = max(d_rep_mm(has_channel) - 1000 .* River_Depth(has_channel), 0);
d_surface_mm(~isfinite(d_surface_mm)) = 0;

end
