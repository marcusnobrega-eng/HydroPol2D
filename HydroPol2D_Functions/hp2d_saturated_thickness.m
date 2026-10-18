function sat = hp2d_saturated_thickness(column, head_m)
%HP2D_SATURATED_THICKNESS Thickness of each slot below the water table, [N x L] (m).
depth = hp2d_water_table_depth(column, head_m);
sat = max(column.z_bot_m - max(depth, column.z_top_m), 0);
end
