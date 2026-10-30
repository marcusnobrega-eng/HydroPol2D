function [due, target_dt_s] = hp2d_groundwater_lateral_due(elapsed_s, target_dt_s, stable_dt_s, ...
    suggested_dt_s, min_dt_s, head_change_m, max_head_change_m, is_final_step)
%HP2D_GROUNDWATER_LATERAL_DUE Whether the lateral groundwater update runs now.
%   The interval is max(min(target, stable, suggested), min_dt): the
%   configured target (groundwater_target_dt_s, default one day), shortened
%   by the Boussinesq stability limit and the adaptive suggestion, but never
%   below groundwater_min_dt_s. The update also runs when the head has moved
%   by max_head_change_m since the last one, and on the last surface step.
%   This is groundwater_update_due in HydroPol2D-Python
%   (modules/groundwater.py), used by its legacy and layered schemes.
%   Pass [] or inf for a limit that does not apply.

if isempty(stable_dt_s), stable_dt_s = inf; end
if isempty(suggested_dt_s), suggested_dt_s = inf; end
if isempty(head_change_m), head_change_m = 0; end
if isempty(max_head_change_m), max_head_change_m = inf; end
target_dt_s = max(min([target_dt_s, stable_dt_s, suggested_dt_s]), min_dt_s);
largest_change = max(abs(head_change_m(:)), [], 'omitnan');
if isempty(largest_change), largest_change = 0; end
due = elapsed_s >= target_dt_s || largest_change >= max_head_change_m || is_final_step;
end
