function T = hp2d_layered_transmissivity(column, head_m)
%HP2D_LAYERED_TRANSMISSIVITY T = sum over slots of Kh * saturated thickness (m2/s), [N x 1].
T = sum(column.kh_m_s.*hp2d_saturated_thickness(column, head_m), 2);
T(~column.valid) = 0;
end
