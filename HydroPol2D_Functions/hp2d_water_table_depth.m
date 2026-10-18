function depth = hp2d_water_table_depth(column, head_m)
%HP2D_WATER_TABLE_DEPTH Water-table depth below ground (m), clipped to [0, column base].
depth = min(max(column.ground_m - head_m(:), 0), column.base_depth_m);
end
