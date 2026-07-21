function [Snow_Properties, Audit] = hp2d_initialize_snow_state( ...
    SnowConfig, LULC, valid_mask, initial_swe_mm, initial_depth_mm, cell_area_m2)
%HP2D_INITIALIZE_SNOW_STATE Map LULC snow parameters and initialize snow.

[unmapped_mask, fallback, Audit] = hp2d_class_code_audit( ...
    LULC, valid_mask, SnowConfig.class_index, SnowConfig.parameter_values, ...
    SnowConfig.parameter_names, 'Snow', cell_area_m2);

map_size = size(LULC);
Snow_Properties = struct();
for j = 1:numel(SnowConfig.parameter_names)
    parameter_name = char(SnowConfig.parameter_names(j));
    Snow_Properties.(parameter_name) = nan(map_size);
end

for i = 1:numel(SnowConfig.class_index)
    class_mask = valid_mask & LULC == SnowConfig.class_index(i);
    for j = 1:numel(SnowConfig.parameter_names)
        parameter_name = char(SnowConfig.parameter_names(j));
        Snow_Properties.(parameter_name)(class_mask) = ...
            SnowConfig.parameter_values(i,j);
    end
end
for j = 1:numel(SnowConfig.parameter_names)
    parameter_name = char(SnowConfig.parameter_names(j));
    Snow_Properties.(parameter_name)(unmapped_mask) = fallback(j);
end

if any(~isfinite(Snow_Properties.rho_snow_init(valid_mask))) || ...
        any(~isfinite(Snow_Properties.rho_max(valid_mask)))
    error('Snow parameters are incomplete after class-code fallback.');
end

[swe_input, has_swe] = normalize_initial_raster(initial_swe_mm, map_size, valid_mask, 'Initial SWE');
[depth_input, has_depth] = normalize_initial_raster(initial_depth_mm, map_size, valid_mask, 'Initial snow depth');

rho_water = 1000;
SWE_t = zeros(map_size);
H_snow_t = zeros(map_size);
rho_snow = Snow_Properties.rho_snow_init;
given_swe = has_swe & isfinite(swe_input);
given_depth = has_depth & isfinite(depth_input);

both = given_swe & given_depth;
inconsistent = both & xor(swe_input > 0, depth_input > 0);
if any(inconsistent(valid_mask))
    error('Initial SWE and snow-depth rasters are inconsistent: a positive value requires both quantities to be positive.');
end

swe_only = given_swe & ~given_depth;
depth_only = given_depth & ~given_swe;
SWE_t(swe_only) = swe_input(swe_only);
H_snow_t(swe_only) = SWE_t(swe_only) .* rho_water ./ rho_snow(swe_only);

H_snow_t(depth_only) = depth_input(depth_only);
SWE_t(depth_only) = H_snow_t(depth_only) .* rho_snow(depth_only) ./ rho_water;

positive_both = both & swe_input > 0 & depth_input > 0;
SWE_t(positive_both) = swe_input(positive_both);
H_snow_t(positive_both) = depth_input(positive_both);
rho_from_rasters = rho_water .* SWE_t(positive_both) ./ H_snow_t(positive_both);
if any(rho_from_rasters <= 0 | rho_from_rasters > Snow_Properties.rho_max(positive_both))
    error('Initial SWE and snow-depth rasters imply a nonphysical snow density.');
end
rho_snow(positive_both) = rho_from_rasters;

SWE_t(~valid_mask) = NaN;
H_snow_t(~valid_mask) = NaN;
rho_snow(~valid_mask) = NaN;
Snow_Properties.SWE_t = SWE_t;
Snow_Properties.H_snow_t = H_snow_t;
Snow_Properties.rho_snow = rho_snow;
Snow_Properties.M_snow = zeros(map_size); Snow_Properties.M_snow(~valid_mask) = NaN;
Snow_Properties.P_snow = zeros(map_size); Snow_Properties.P_snow(~valid_mask) = NaN;
Snow_Properties.P_rain = zeros(map_size); Snow_Properties.P_rain(~valid_mask) = NaN;
Snow_Properties.E_s = zeros(map_size); Snow_Properties.E_s(~valid_mask) = NaN;
end

function [raster, has_input] = normalize_initial_raster(value, map_size, valid_mask, label)
has_input = ~isempty(value);
raster = nan(map_size);
if ~has_input
    return
end
if ~isequal(size(value), map_size)
    error('%s raster does not match the DEM grid after alignment.', label);
end
raster = double(value);
if any(raster(valid_mask) < 0)
    error('%s raster contains negative values.', label);
end
raster(~valid_mask) = NaN;
end
