function unsat = hp2d_unsaturated_thickness(column, head_m)
%HP2D_UNSATURATED_THICKNESS Thickness of each slot above the water table, [N x L] (m).
depth = hp2d_water_table_depth(column, head_m);
unsat = max(min(depth, column.z_bot_m) - column.z_top_m, 0);
end
