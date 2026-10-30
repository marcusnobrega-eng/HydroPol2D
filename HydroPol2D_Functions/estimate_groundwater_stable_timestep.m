function dt_stable_s = estimate_groundwater_stable_timestep(Soil_Properties, BC_States, z_bed, dx, dy, idx_nan, courant)
%ESTIMATE_GROUNDWATER_STABLE_TIMESTEP Conservative explicit GW timestep cap.
%   courant * min(dx,dy)^2 / max(K b / Sy) / 4, at least 1 s
%   (estimate_boussinesq_stable_dt_s in HydroPol2D-Python). courant defaults
%   to 0.25.

Sy = Soil_Properties.Sy;
Sy = max(Sy, 1e-6);
K = Soil_Properties.ksat_gw ./ 1000 ./ 3600;
K(~isfinite(K)) = nan;

saturated_thickness = max(BC_States.h_t - z_bed, 1e-4);
saturated_thickness(idx_nan) = nan;

diffusivity = K .* saturated_thickness ./ Sy;
max_diffusivity = max(diffusivity(:), [], 'omitnan');

if isempty(max_diffusivity) || ~isfinite(max_diffusivity) || max_diffusivity <= 0
    dt_stable_s = inf;
    return
end

if nargin < 7 || isempty(courant), courant = 0.25; end
grid_length = min(dx, dy);
dt_stable_s = courant * grid_length^2 / max_diffusivity / 4;
dt_stable_s = max(dt_stable_s, 1);
end
