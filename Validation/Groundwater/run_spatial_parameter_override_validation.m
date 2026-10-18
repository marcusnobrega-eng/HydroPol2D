function report = run_spatial_parameter_override_validation()
%RUN_SPATIAL_PARAMETER_OVERRIDE_VALIDATION Focused raster/table precedence checks.
root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(root, 'HydroPol2D_Functions'));

soil = struct();
soil.theta_r = [0.05, 0.06];
soil.theta_sat = [0.45, 0.46];
soil.alpha_vg = [2.0, 2.5];
soil.n_vg = [1.5, 1.6];
soil.Sy = [0.20, 0.22];
soil.ksat_gw = [20, 30];
soil.ksat = [10, 20];
soil.Ks_multiplier_near_surface = [2, 2];
soil.Ks_multiplier_root_zone = [0.5, 0.5];
soil.Ks_multiplier_transmission = [0.25, 0.25];
lulc.root_depth_m = [1.0, 1.5];
active = true(1, 2);

rasters = struct( ...
    'ThetaR', [0.08, nan], ...
    'ThetaS', [0.48, nan], ...
    'AlphaVG', [1.2, nan], ...
    'NVG', [1.8, nan], ...
    'SpecificYield', [0.12, nan], ...
    'GroundwaterKsat', [0.7, nan], ...
    'RootDepth', [0.9, nan], ...
    'KsatNearSurface', [7.0, nan], ...
    'KsatRootZone', [4.0, nan], ...
    'KsatTransmission', [2.0, nan]);

[mapped, mapped_lulc, audit] = hp2d_apply_spatial_parameter_overrides( ...
    soil, lulc, rasters, active);

assert(max(abs(mapped.theta_r - [0.08, 0.06])) < 1e-12);
assert(max(abs(mapped.theta_sat - [0.48, 0.46])) < 1e-12);
assert(max(abs(mapped_lulc.root_depth_m - [0.9, 1.5])) < 1e-12);
assert(max(abs(mapped.ksat - [7.0, 40.0])) < 1e-12);
assert(max(abs(mapped.ksat .* mapped.Ks_multiplier_root_zone - [4.0, 10.0])) < 1e-12);
assert(max(abs(mapped.ksat .* mapped.Ks_multiplier_transmission - [2.0, 5.0])) < 1e-12);
assert(audit.ThetaR == 1 && audit.KsatNearSurface == 1);

invalid = rasters;
invalid.SpecificYield = [0.0, nan];
rejected = false;
try
    hp2d_apply_spatial_parameter_overrides(soil, lulc, invalid, active);
catch error_info
    rejected = contains(error_info.message, 'specific yield');
end
assert(rejected, 'Invalid specific yield was not rejected.');

report = struct( ...
    'finite_raster_overrides_table', true, ...
    'nodata_uses_table_fallback', true, ...
    'layer_specific_ksat_preserved', true, ...
    'invalid_values_rejected', true);
end
