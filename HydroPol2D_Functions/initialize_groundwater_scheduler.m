function [flags, GW_States] = initialize_groundwater_scheduler(flags, GW_States, template, idx_nan)
%INITIALIZE_GROUNDWATER_SCHEDULER Defaults and state for async GW coupling.

if nargin < 4 || isempty(idx_nan)
    idx_nan = false(size(template));
end

flags = local_default_flag(flags, 'flag_groundwater_async', 1);
flags = local_default_flag(flags, 'flag_capillary_rise', 1);
flags = local_default_flag(flags, 'groundwater_target_dt_min', 1440);
flags = local_default_flag(flags, 'groundwater_min_dt_min', 1);
flags = local_default_flag(flags, 'groundwater_max_head_change_m', 0.25);
flags = local_default_flag(flags, 'groundwater_courant', 0.25);

if ~isfield(GW_States, 'cell_area_m2') || isempty(GW_States.cell_area_m2)
    GW_States.cell_area_m2 = ones(size(template), 'like', template);
end

state_fields = { ...
    'pending_net_exchange_m', ...
    'pending_recharge_m', ...
    'pending_capillary_m', ...
    'last_surface_recharge_rate_m_s', ...
    'last_capillary_rate_m_s', ...
    'last_net_exchange_rate_m_s', ...
    'last_groundwater_recharge_rate_m_s', ...
    'last_q_exf_m_s', ...
    'last_q_river_m_s'};

for i = 1:numel(state_fields)
    name = state_fields{i};
    if ~isfield(GW_States, name) || isempty(GW_States.(name)) || ...
            ~isequal(size(GW_States.(name)), size(template))
        value = zeros(size(template), 'like', template);
        value(idx_nan) = nan;
        GW_States.(name) = value;
    end
end

scalar_defaults = { ...
    'elapsed_s', 0; ...
    'last_groundwater_dt_s', 0; ...
    'n_updates', 0; ...
    'total_recharge_m3', 0; ...
    'total_capillary_m3', 0; ...
    'total_exfiltration_m3', 0};

for i = 1:size(scalar_defaults, 1)
    name = scalar_defaults{i, 1};
    if ~isfield(GW_States, name) || isempty(GW_States.(name)) || ...
            ~isfinite(GW_States.(name))
        GW_States.(name) = scalar_defaults{i, 2};
    end
end

GW_States.target_dt_s = max(double(flags.groundwater_target_dt_min) * 60, eps);
GW_States.min_dt_s = max(double(flags.groundwater_min_dt_min) * 60, eps);
GW_States.max_head_change_m = max(double(flags.groundwater_max_head_change_m), eps);
GW_States.groundwater_courant = max(min(double(flags.groundwater_courant), 1), eps);
end

function flags = local_default_flag(flags, name, default_value)
if ~isfield(flags, name) || isempty(flags.(name)) || ~isfinite(double(flags.(name)))
    flags.(name) = default_value;
end
end
