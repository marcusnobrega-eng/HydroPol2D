function h = hp2d_neal_apply_volume_change(h, delta_volume_m3, River_Width, River_Depth, Resolution)
%HP2D_NEAL_APPLY_VOLUME_CHANGE Conservatively add volume to Neal cells.

V = hp2d_neal_cell_volume(h, River_Width, River_Depth, Resolution);
V = max(V + delta_volume_m3, 0);
h = hp2d_neal_depth_from_volume(V, River_Width, River_Depth, Resolution);

end
