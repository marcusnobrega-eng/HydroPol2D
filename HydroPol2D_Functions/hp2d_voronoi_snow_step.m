function [snow, rain_m_s, melt_volume_m3, sublimation_volume_m3] = hp2d_voronoi_snow_step(snow, precipitation_m_s, meteorology, dt_s)
%HP2D_VORONOI_SNOW_STEP Snow on every HRU of a mesh, before the canopy.
%   precipitation_m_s is the gross precipitation of each cell [m/s].
%   meteorology holds temperature_c, minimum_temperature_c, wind_speed_m_s
%   and latitude_deg per cell, and a uniform day_of_year. The pack of each
%   HRU advances with Snow_Model_Function, as on the grid (Snow_Module.m)
%   and as HydroPol2D-Python's mesh runs snow_step on its HRU pieces.
%   rain_m_s [n x h] is the liquid part, which goes on to the canopy;
%   melt_volume_m3 [n x 1] leaves the pack on the ground and goes straight
%   to the surface. sublimation_volume_m3 leaves the model.

h = size(snow.swe_mm, 2);
on_hrus = @(value) repmat(double(value(:)), 1, h);
P_mm = on_hrus(max(precipitation_m_s, 0)) * dt_s * 1000;
swe_before = max(snow.swe_mm, 0);
[snow.swe_mm, snow.depth_mm, snow.melt_mm, P_snow, P_rain, snow.density_kg_m3] = ...
    Snow_Model_Function(snow.swe_mm, snow.depth_mm, ...
    snow.density_kg_m3, on_hrus(meteorology.temperature_c), ...
    on_hrus(meteorology.minimum_temperature_c), P_mm, ...
    on_hrus(meteorology.wind_speed_m_s), on_hrus(meteorology.latitude_deg), ...
    meteorology.day_of_year, snow.alpha, snow.epsilon, snow.C_e, snow.DDF, ...
    snow.T_snow_all, snow.T_rain_all, snow.rho_snow_init, snow.rho_max, ...
    snow.k_t, snow.k_swe, snow.k_D, dt_s);
rain_m_s = max(P_rain, 0) / 1000 / dt_s;
melt_volume_m3 = sum(max(snow.melt_mm, 0) / 1000 .* snow.area_m2, 2);
% Sublimation from the pack's own balance, so the ledger closes exactly.
snow.sublimation_mm = max(swe_before + P_snow - snow.melt_mm - snow.swe_mm, 0);
sublimation_volume_m3 = sum(snow.sublimation_mm / 1000 .* snow.area_m2, 'all');
snow.cumulative_sublimation_m3 = snow.cumulative_sublimation_m3 + sublimation_volume_m3;
end
