function [Subsurface, BC_States] = hp2d_prepare_layered_subsurface(cfg, Soil_Properties, LULC_Properties, BC_States, elevation, idx_nan, ParameterRasters, dtb_raster_depth_m, layer_theta)
%HP2D_PREPARE_LAYERED_SUBSURFACE Build the layered column, bedrock and initial slot water.
%   Column bottom (m below ground): the DTB raster where it is finite (dtb_raster_depth_m,
%   [] if none), else cfg.subsurface_depth_to_bedrock_m, else Soil_Depth with a warning;
%   capped at cfg.subsurface_max_depth_m. Soil base: SoilBaseDepth raster, else the
%   scalar. Rock parameters: Regolith* rasters where finite, else the scalars; Kh falls
%   back to ksat_gw and Sy to the retention curve. Same precedence as HydroPol2D-Python.
%   layer_theta: {near, root, transmission} initial theta maps, or [] for theta_i.
[ny, nx] = size(elevation);
valid = ~idx_nan(:);
N = ny*nx;
    function v = cells(map)
        v = double(map(:)) + zeros(N, 1);
    end
    function v = raster_or(name, fallback)
        v = cells(fallback);
        if isfield(ParameterRasters, name) && ~isempty(ParameterRasters.(name))
            candidate = double(ParameterRasters.(name).Z(:));
            assert(numel(candidate) == N, 'HydroPol2D:Subsurface', ...
                'Raster %s does not match the model grid.', name);
            use = valid & isfinite(candidate);
            v(use) = candidate(use);
        end
    end

legacy_depth = cells(Soil_Properties.Soil_Depth);
if ~isempty(dtb_raster_depth_m)
    if isfield(ParameterRasters, 'SoilDepthParameter') && ~isempty(ParameterRasters.SoilDepthParameter)
        warning('HydroPol2D:Subsurface', ['SoilDepthParameter is not used by the layered ', ...
            'subsurface; the DTB raster sets the column bottom.']);
    end
    dtb = cells(dtb_raster_depth_m);
    fallback = legacy_depth;
    if ~isempty(cfg.subsurface_depth_to_bedrock_m)
        fallback(:) = cfg.subsurface_depth_to_bedrock_m;
    end
    dtb(~isfinite(dtb)) = fallback(~isfinite(dtb));
elseif ~isempty(cfg.subsurface_depth_to_bedrock_m)
    dtb = cfg.subsurface_depth_to_bedrock_m + zeros(N, 1);
else
    warning('HydroPol2D:Subsurface', ['Layered subsurface: no DTB raster or ', ...
        'subsurface_depth_to_bedrock_m, so Soil_Depth is used as the depth to bedrock.']);
    dtb = legacy_depth;
end
if isfield(ParameterRasters, 'SpecificYield') && ~isempty(ParameterRasters.SpecificYield)
    warning('HydroPol2D:Subsurface', ['SpecificYield is not used by the layered subsurface: ', ...
        'Sy comes from each slot''s retention curve, or regolith_specific_yield for rock.']);
end
base = min(dtb, cfg.subsurface_max_depth_m);
base(~valid) = 0;
soil_base = min(raster_or('SoilBaseDepth', cfg.subsurface_soil_base_depth_m), base);
soil_base(~valid) = 0;
root = cells(LULC_Properties.root_depth_m);
root(~isfinite(root)) = 0;
ltop = cells(Soil_Properties.Ltop);
ltop(~isfinite(ltop)) = 0.05;
[z_top, z_bot] = hp2d_slot_geometry(root, soil_base, base, ltop, cfg.deep_thicknesses_m);

L = Soil_Properties.Layers;
P.soil_theta_s = cells(Soil_Properties.theta_sat);
P.soil_theta_r = cells(Soil_Properties.theta_r);
P.soil_alpha_m_inverse = cells(Soil_Properties.alpha_vg);
P.soil_n_vg = cells(Soil_Properties.n_vg);
P.soil_kv_m_s = [cells(L.ksat_near_surface), cells(L.ksat_root_zone), cells(L.ksat_transmission)]/1000/3600;
P.soil_lateral_anisotropy = cfg.soil_lateral_anisotropy;
P.regolith_theta_s = raster_or('RegolithThetaS', cfg.regolith_theta_s);
P.regolith_theta_r = raster_or('RegolithThetaR', cfg.regolith_theta_r);
P.regolith_alpha_m_inverse = raster_or('RegolithAlphaVG', cfg.regolith_alpha_vg_m_inverse);
P.regolith_n_vg = raster_or('RegolithNVG', cfg.regolith_n_vg);
P.regolith_kv_m_s = raster_or('RegolithKsatV', cfg.regolith_ksat_v_mm_h)/1000/3600;
kh_fallback = cells(Soil_Properties.ksat_gw);
if ~isempty(cfg.regolith_ksat_h_mm_h), kh_fallback(:) = cfg.regolith_ksat_h_mm_h; end
P.regolith_kh_m_s = raster_or('RegolithKsatH', kh_fallback)/1000/3600;
sy_fallback = nan(N, 1);
if ~isempty(cfg.regolith_specific_yield), sy_fallback(:) = cfg.regolith_specific_yield; end
P.regolith_specific_yield = raster_or('RegolithSpecificYield', sy_fallback);
P.retention_head_m = cfg.capillary_field_capacity_head_m;
P.ground_m = cells(elevation);
P.column_base_depth_m = base;
P.valid = valid;
for name = {'soil_theta_s', 'soil_theta_r', 'soil_alpha_m_inverse', 'soil_n_vg', ...
        'regolith_theta_s', 'regolith_theta_r', 'regolith_alpha_m_inverse', 'regolith_n_vg', ...
        'regolith_kv_m_s', 'regolith_kh_m_s', 'ground_m'}
    values = P.(name{1});
    values(~valid) = 0;
    P.(name{1}) = values;
end
P.soil_kv_m_s(~valid, :) = 0;
bad = valid & ~(P.regolith_theta_r < P.regolith_theta_s & P.regolith_theta_s <= 1 & ...
    P.regolith_theta_r >= 0 & P.regolith_alpha_m_inverse > 0 & P.regolith_n_vg > 1 & ...
    P.regolith_kv_m_s >= 0 & P.regolith_kh_m_s >= 0);
if any(bad)
    error('HydroPol2D:Subsurface', 'Invalid regolith parameters in %d active cell(s).', nnz(bad));
end
sy = P.regolith_specific_yield;
bad = valid & isfinite(sy) & ~(sy > 0 & sy <= P.regolith_theta_s - P.regolith_theta_r + 1e-12);
if any(bad)
    error('HydroPol2D:Subsurface', ...
        'Regolith specific yield lies outside (0, theta_s - theta_r] in %d active cell(s).', nnz(bad));
end
sy(~valid) = NaN;
P.regolith_specific_yield = sy;
column = hp2d_slot_column(z_top, z_bot, P);

z0 = reshape(P.ground_m - base, ny, nx);
z0(idx_nan) = nan;
head = double(BC_States.h_t);
head = min(max(head, z0), elevation);
head(idx_nan) = nan;
BC_States.z0 = z0;
BC_States.h_0 = head;
BC_States.h_t = head;

if nargin > 8 && iscell(layer_theta) && numel(layer_theta) == 3
    theta3 = [cells(layer_theta{1}), cells(layer_theta{2}), cells(layer_theta{3})];
else
    theta3 = repmat(cells(Soil_Properties.theta_i), 1, 3);
end
theta3(~isfinite(theta3)) = 0;
head_cells = head(:);
head_cells(~valid) = 0;
water = hp2d_initial_slot_water(column, head_cells, theta3);

Subsurface = struct();
Subsurface.column = column;
Subsurface.water_mm = water;
Subsurface.ny = ny;
Subsurface.nx = nx;
Subsurface.n_slots = size(z_top, 2);
Subsurface.valid = valid;
Subsurface.config = cfg;
Subsurface.feddes = cfg.feddes;
Subsurface.drainage_substeps = cfg.subsurface_drainage_substeps;
Subsurface.column_base_depth_m = reshape(base, ny, nx);
Subsurface.soil_base_depth_m = reshape(soil_base, ny, nx);
Subsurface.column_base_depth_m(idx_nan) = nan;
Subsurface.soil_base_depth_m(idx_nan) = nan;
end
