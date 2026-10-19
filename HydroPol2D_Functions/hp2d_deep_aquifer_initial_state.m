function [storage_m3, discharge_m3_s, subbasin_ids] = hp2d_deep_aquifer_initial_state( ...
    subbasin_id, cell_area_m2, response_time_days, recharge_fraction, ...
    initial_baseflow_mm_day, mean_recharge_mm_year)
%HP2D_DEEP_AQUIFER_INITIAL_STATE Linear-reservoir starting state per subbasin.
%   Each subbasin starts at S0 = q0 x area x tau, with q0 the deep baseflow per unit
%   area (HydroPol2D-Python parity):
%     initial_baseflow_mm_day  q0 given, e.g. from a gauge (estimate-deep-baseflow);
%     mean_recharge_mm_year    steady state, q0 = recharge_fraction x mean recharge.
%   Give at most one (the other NaN); both NaN gives empty reservoirs. cell_area_m2 is a
%   scalar or one value per cell of subbasin_id. Rows follow the sorted positive
%   subbasin IDs, the order used by hp2d_deep_aquifer_step.
ids = round(double(subbasin_id(:)));
active = isfinite(ids) & ids > 0;
subbasin_ids = unique(ids(active));
has_baseflow = isfinite(initial_baseflow_mm_day);
has_recharge = isfinite(mean_recharge_mm_year);
assert(~(has_baseflow && has_recharge), ...
    'Set deep_aquifer_initial_baseflow_mm_day or deep_aquifer_mean_recharge_mm_year, not both.');
if has_baseflow
    assert(initial_baseflow_mm_day >= 0, 'deep_aquifer_initial_baseflow_mm_day must be nonnegative.');
    q0_m_s = initial_baseflow_mm_day/1000/86400;
elseif has_recharge
    assert(mean_recharge_mm_year >= 0, 'deep_aquifer_mean_recharge_mm_year must be nonnegative.');
    q0_m_s = recharge_fraction*mean_recharge_mm_year/1000/(365.25*86400);
else
    q0_m_s = 0;
end
[~, compact] = ismember(ids(active), subbasin_ids);
if isscalar(cell_area_m2)
    area_m2 = accumarray(compact, 1, [numel(subbasin_ids) 1])*double(cell_area_m2);
else
    cell_area = double(cell_area_m2(:));
    area_m2 = accumarray(compact, cell_area(active), [numel(subbasin_ids) 1]);
end
discharge_m3_s = q0_m_s.*area_m2;
storage_m3 = discharge_m3_s.*response_time_days.*86400;
end
