function [state, shallow_recharge_rate_m_s, surface_addition_mm] = ...
    hp2d_deep_aquifer_step(state, recharge_rate_m_s, subbasin_id, ...
    channel_weight, recharge_fraction, response_time_days, cell_area_m2, dt_s)
%HP2D_DEEP_AQUIFER_STEP Exact subbasin linear reservoirs with channel return.

assert(recharge_fraction > 0 && recharge_fraction <= 1, ...
    'deep_recharge_fraction must be in (0,1].');
assert(response_time_days > 0 && dt_s > 0, ...
    'Deep-aquifer response time and timestep must be positive.');

if ~isfield(state,'deep_aquifer_compact_index') || ...
        ~isequal(size(state.deep_aquifer_compact_index),size(subbasin_id))
    ids = round(double(subbasin_id));
    active = isfinite(ids) & ids > 0;
    unique_ids = unique(ids(active));
    n = numel(unique_ids);
    assert(n > 0,'Deep-aquifer subbasin map contains no positive IDs.');
    compact = zeros(size(ids));
    weights = zeros(size(channel_weight));
    for k = 1:n
        cells = active & ids == unique_ids(k);
        compact(cells) = k;
        total_weight = sum(channel_weight(cells),'omitnan');
        assert(total_weight > 0, ...
            'Deep-aquifer subbasin %g has no positive channel weight.', unique_ids(k));
        weights(cells) = max(channel_weight(cells),0) ./ total_weight;
    end
    state.deep_aquifer_compact_index = compact;
    state.deep_aquifer_channel_weight = weights;
    state.deep_aquifer_subbasin_ids = unique_ids(:);
else
    compact = state.deep_aquifer_compact_index;
    weights = state.deep_aquifer_channel_weight;
    active = compact > 0;
    n = numel(state.deep_aquifer_subbasin_ids);
end

if ~isfield(state,'deep_aquifer_storage_m3') || ...
        numel(state.deep_aquifer_storage_m3) ~= n
    state.deep_aquifer_storage_m3 = zeros(n,1);
    state.cumulative_deep_aquifer_recharge_m3 = zeros(n,1);
    state.cumulative_deep_aquifer_discharge_m3 = zeros(n,1);
end

deep_depth_m = max(recharge_rate_m_s,0) .* dt_s .* recharge_fraction;
input_m3 = accumarray(compact(active), ...
    deep_depth_m(active) .* cell_area_m2, [n 1], @sum, 0);
input_rate_m3_s = input_m3 ./ dt_s;
tau_s = response_time_days .* 86400;
decay = exp(-dt_s ./ tau_s);
storage_old = max(state.deep_aquifer_storage_m3,0);
storage_new = storage_old .* decay + input_rate_m3_s .* tau_s .* (1-decay);
discharge_m3 = max(input_m3 - (storage_new-storage_old),0);

surface_addition_m3 = zeros(size(compact));
for k = 1:n
    cells = compact == k;
    surface_addition_m3(cells) = discharge_m3(k) .* weights(cells);
end

state.deep_aquifer_storage_m3 = storage_new;
state.deep_aquifer_discharge_m3_s = discharge_m3 ./ dt_s;
state.cumulative_deep_aquifer_recharge_m3 = ...
    state.cumulative_deep_aquifer_recharge_m3 + input_m3;
state.cumulative_deep_aquifer_discharge_m3 = ...
    state.cumulative_deep_aquifer_discharge_m3 + discharge_m3;
shallow_recharge_rate_m_s = recharge_rate_m_s .* (1-recharge_fraction);
surface_addition_mm = surface_addition_m3 ./ cell_area_m2 .* 1000;
surface_addition_mm(~active) = 0;
end
