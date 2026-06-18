function dt_stable_s = estimate_groundwater_stable_timestep(Soil_Properties, BC_States, z_bed, dx, dy, flags, idx_nan)
%ESTIMATE_GROUNDWATER_STABLE_TIMESTEP Conservative explicit GW timestep cap.

target_s = 1440 * 60;
if isfield(flags, 'groundwater_target_dt_min') && isfinite(double(flags.groundwater_target_dt_min))
    target_s = max(double(flags.groundwater_target_dt_min) * 60, eps);
end

if ~isfield(flags, 'flag_baseflow') || flags.flag_baseflow ~= 1
    dt_stable_s = target_s;
    return
end

courant = 0.25;
if isfield(flags, 'groundwater_courant') && isfinite(double(flags.groundwater_courant))
    courant = max(min(double(flags.groundwater_courant), 1), eps);
end

Sy = Soil_Properties.Sy;
Sy = max(Sy, 1e-6);
K = Soil_Properties.ksat_gw ./ 1000 ./ 3600;
K(~isfinite(K)) = nan;

saturated_thickness = max(BC_States.h_t - z_bed, 1e-4);
saturated_thickness(idx_nan) = nan;

diffusivity = K .* saturated_thickness ./ Sy;
max_diffusivity = max(diffusivity(:), [], 'omitnan');

if isempty(max_diffusivity) || ~isfinite(max_diffusivity) || max_diffusivity <= 0
    dt_stable_s = target_s;
    return
end

grid_length = min(dx, dy);
dt_stable_s = courant * grid_length^2 / max_diffusivity / 4;
dt_stable_s = min(max(dt_stable_s, 1), target_s);
end
