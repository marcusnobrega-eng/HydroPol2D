function [soil, lulc, audit] = hp2d_apply_spatial_parameter_overrides(soil, lulc, rasters, active)
%HP2D_APPLY_SPATIAL_PARAMETER_OVERRIDES Overlay finite maps on class values.
% Raster nodata remains non-finite and therefore retains the class-table
% value. Layer-specific Ksat maps are represented through the existing
% base-Ksat and multiplier fields without changing downstream interfaces.

audit = struct();
[soil.theta_r, audit.ThetaR] = local_override(soil.theta_r, rasters, 'ThetaR', active);
[soil.theta_sat, audit.ThetaS] = local_override(soil.theta_sat, rasters, 'ThetaS', active);
[soil.alpha_vg, audit.AlphaVG] = local_override(soil.alpha_vg, rasters, 'AlphaVG', active);
[soil.n_vg, audit.NVG] = local_override(soil.n_vg, rasters, 'NVG', active);
[soil.Sy, audit.SpecificYield] = local_override(soil.Sy, rasters, 'SpecificYield', active);
[soil.ksat_gw, audit.GroundwaterKsat] = local_override( ...
    soil.ksat_gw, rasters, 'GroundwaterKsat', active);
[lulc.root_depth_m, audit.RootDepth] = local_override( ...
    lulc.root_depth_m, rasters, 'RootDepth', active);

ksat_near = soil.ksat .* soil.Ks_multiplier_near_surface;
ksat_root = soil.ksat .* soil.Ks_multiplier_root_zone;
ksat_transmission = soil.ksat .* soil.Ks_multiplier_transmission;
[ksat_near, audit.KsatNearSurface] = local_override( ...
    ksat_near, rasters, 'KsatNearSurface', active);
[ksat_root, audit.KsatRootZone] = local_override( ...
    ksat_root, rasters, 'KsatRootZone', active);
[ksat_transmission, audit.KsatTransmission] = local_override( ...
    ksat_transmission, rasters, 'KsatTransmission', active);

soil.ksat = ksat_near;
positive_near = active & isfinite(ksat_near) & ksat_near > 0;
soil.Ks_multiplier_near_surface(active) = 1;
soil.Ks_multiplier_root_zone(positive_near) = ...
    ksat_root(positive_near) ./ ksat_near(positive_near);
soil.Ks_multiplier_transmission(positive_near) = ...
    ksat_transmission(positive_near) ./ ksat_near(positive_near);
soil.Ks_multiplier_root_zone(active & ~positive_near) = 0;
soil.Ks_multiplier_transmission(active & ~positive_near) = 0;
soil.m_vg = 1 - 1 ./ soil.n_vg;

local_validate(soil, lulc, active);
end

function [values, count] = local_override(values, rasters, name, active)
count = 0;
if ~isfield(rasters, name) || isempty(rasters.(name))
    return
end
candidate = rasters.(name);
if isobject(candidate) && isprop(candidate, 'Z')
    candidate = candidate.Z;
end
candidate = double(candidate);
if ~isequal(size(candidate), size(values))
    error('Spatial parameter raster %s does not match the model grid.', name);
end
use = active & isfinite(candidate);
values(use) = candidate(use);
count = nnz(use);
end

function local_validate(soil, lulc, active)
if any(~isfinite(soil.theta_r(active)) | soil.theta_r(active) < 0) || ...
        any(~isfinite(soil.theta_sat(active)) | soil.theta_sat(active) <= 0 | ...
        soil.theta_sat(active) > 1) || ...
        any(soil.theta_r(active) >= soil.theta_sat(active))
    error('Spatial theta_r/theta_s values are invalid. Require 0 <= theta_r < theta_s <= 1.');
end
if any(~isfinite(soil.alpha_vg(active)) | soil.alpha_vg(active) <= 0) || ...
        any(~isfinite(soil.n_vg(active)) | soil.n_vg(active) <= 1)
    error('Spatial van Genuchten values require alpha > 0 and n > 1.');
end
if any(~isfinite(soil.ksat(active)) | soil.ksat(active) < 0) || ...
        any(~isfinite(soil.Ks_multiplier_root_zone(active)) | ...
        soil.Ks_multiplier_root_zone(active) < 0) || ...
        any(~isfinite(soil.Ks_multiplier_transmission(active)) | ...
        soil.Ks_multiplier_transmission(active) < 0) || ...
        any(~isfinite(soil.ksat_gw(active)) | soil.ksat_gw(active) < 0)
    error('Spatial hydraulic conductivities must be finite and nonnegative.');
end
if any(~isfinite(soil.Sy(active)) | soil.Sy(active) <= 0 | soil.Sy(active) > 1)
    error('Spatial specific yield must be in (0, 1].');
end
if any(~isfinite(lulc.root_depth_m(active)) | lulc.root_depth_m(active) < 0)
    error('Spatial root depth must be finite and nonnegative.');
end
end
