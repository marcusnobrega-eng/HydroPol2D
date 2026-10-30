function [surface_depth_mm, channel_depth_mm] = hp2d_d8_surface_depth(D8, d_t_mm)
%HP2D_D8_SURFACE_DEPTH Map depths of a D8 local-inertial run [mm].
%   d_t_mm is the cell volume over the cell area. surface_depth_mm is the
%   water above the terrain, max(stage - terrain, 0): zero while a river
%   stays in its banks. channel_depth_mm is the water above the channel bed
%   on river cells (stage - datum) and zero elsewhere. These are the
%   surface_depth_map_m and channel_depth_map_m of HydroPol2D-Python.

volume_m3 = double(gather(d_t_mm(:))) / 1000 * D8.cell_area;
volume_m3(~isfinite(volume_m3)) = 0;
[stage_m, channel_m] = hp2d_d8_cell_stage(D8.Cells, volume_m3);
surface_depth_mm = reshape(1000 * max(stage_m - D8.Cells.terrain, 0), D8.ny, D8.nx);
channel_depth_mm = reshape(1000 * channel_m, D8.ny, D8.nx);
surface_depth_mm(~D8.valid_grid) = NaN;
channel_depth_mm(~D8.valid_grid) = NaN;

end
