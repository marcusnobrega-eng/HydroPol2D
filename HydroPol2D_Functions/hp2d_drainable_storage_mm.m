function drainable = hp2d_drainable_storage_mm(column, head_m)
%HP2D_DRAINABLE_STORAGE_MM Water a falling table can release before the column base (mm).
sat = hp2d_saturated_thickness(column, head_m);
drainable = sum(max(column.theta_s - column.theta_retained, 0).*sat*1000, 2);
end
