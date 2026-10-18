function sy = hp2d_effective_specific_yield(column, head_m)
%HP2D_EFFECTIVE_SPECIFIC_YIELD theta_s - theta_retained of the slot just below the table.
%   Used for time-step limits; the deepest slot's value where the table is at the base.
depth = hp2d_water_table_depth(column, head_m);
inside = column.z_top_m <= depth & depth < column.z_bot_m;
yield_ = max(column.theta_s - column.theta_retained, 0);
sy = sum(yield_.*inside, 2);
none = ~any(inside, 2);
sy(none) = yield_(none, end);
end
