function storage = hp2d_column_storage_mm(column, head_m, water_mm)
%HP2D_COLUMN_STORAGE_MM Water above theta_r in the whole column (mm over the cell), [N x 1].
%   Unsaturated water plus (theta_s - theta_r) times the saturated thickness.
sat = hp2d_saturated_thickness(column, head_m);
storage = sum(water_mm, 2) + sum((column.theta_s - column.theta_r).*sat*1000, 2);
storage(~column.valid) = 0;
end
