function [cfg, flags] = hp2d_layered_subsurface_config(flags)
%HP2D_LAYERED_SUBSURFACE_CONFIG Settings of the layered subsurface and Feddes stress.
%   Defaults equal HydroPol2D-Python (vertical.subsurface, vertical.feddes). Values come
%   from flags fields (Excel 'Flags' sheet rows or InputData_Bypass.flags), and the
%   Builder's environment variables win: HYDROPOL_SUBSURFACE_SCHEME (legacy|layered),
%   HYDROPOL_ENABLE_FEDDES (0|1) and HYDROPOL_<NAME> in upper case for every numeric
%   setting below (for example HYDROPOL_REGOLITH_SPECIFIC_YIELD).
%   Empty depth_to_bedrock_m, regolith_ksat_h_mm_h and regolith_specific_yield keep
%   the fallbacks: DTB raster / Soil_Depth, ksat_gw, and the retention-curve Sy.
%   soil_ / regolith_lateral_conductivity_multiplier scale the lateral K of the soil slots
%   (Ksat x anisotropy) and of the deep slots (rock Kh); the recession calibration fits them.
defaults = { ...
    'subsurface_soil_base_depth_m', 2.0; 'subsurface_depth_to_bedrock_m', []; ...
    'subsurface_max_depth_m', 50.0; 'subsurface_deep_first_thickness_m', 0.25; ...
    'subsurface_deep_growth_ratio', 2.0; 'subsurface_deep_max_thickness_m', 4.0; ...
    'regolith_theta_s', 0.35; 'regolith_theta_r', 0.05; 'regolith_alpha_vg_m_inverse', 2.0; ...
    'regolith_n_vg', 1.5; 'regolith_ksat_v_mm_h', 5.0; 'regolith_ksat_h_mm_h', []; ...
    'regolith_specific_yield', []; 'soil_lateral_anisotropy', 1.0; ...
    'soil_lateral_conductivity_multiplier', 1.0; 'regolith_lateral_conductivity_multiplier', 1.0; ...
    'subsurface_drainage_substeps', 1; 'capillary_field_capacity_head_m', 3.3; ...
    'groundwater_target_dt_s', 3600; 'groundwater_min_dt_s', 60; 'groundwater_courant', 0.25; ...
    'feddes_h1_m', -0.1; 'feddes_h2_m', -0.25; 'feddes_h3_m', -4.0; 'feddes_h3_low_m', []; ...
    'feddes_h3_high_m', []; 'feddes_h4_m', -80.0};
    function value = setting(name, default)
        value = default;
        if isfield(flags, name) && ~isempty(flags.(name))
            candidate = double(flags.(name));
            if isscalar(candidate) && isfinite(candidate)
                value = candidate;
            end
        end
        text = strtrim(getenv(['HYDROPOL_' upper(name)]));
        if ~isempty(text)
            candidate = str2double(text);
            if isfinite(candidate)
                value = candidate;
            end
        end
    end
cfg = struct();
for i = 1:size(defaults, 1)
    cfg.(defaults{i, 1}) = setting(defaults{i, 1}, defaults{i, 2});
end

layered = isfield(flags, 'flag_layered_subsurface') && double(flags.flag_layered_subsurface) == 1;
scheme = lower(strtrim(getenv('HYDROPOL_SUBSURFACE_SCHEME')));
if ~isempty(scheme)
    assert(any(strcmp(scheme, {'legacy', 'layered'})), 'HydroPol2D:Subsurface', ...
        'HYDROPOL_SUBSURFACE_SCHEME must be legacy or layered.');
    layered = strcmp(scheme, 'layered');
end
feddes_on = isfield(flags, 'flag_feddes') && double(flags.flag_feddes) == 1;
feddes_text = strtrim(getenv('HYDROPOL_ENABLE_FEDDES'));
if ~isempty(feddes_text)
    feddes_on = strcmp(feddes_text, '1');
end
flags.flag_layered_subsurface = double(layered);
flags.flag_feddes = double(feddes_on);
cfg.layered = layered;

if cfg.regolith_theta_r >= cfg.regolith_theta_s
    error('HydroPol2D:Subsurface', 'regolith_theta_r must be below regolith_theta_s.');
end
if ~isempty(cfg.regolith_specific_yield) && ~(cfg.regolith_specific_yield > 0 && ...
        cfg.regolith_specific_yield <= cfg.regolith_theta_s - cfg.regolith_theta_r)
    error('HydroPol2D:Subsurface', ...
        'regolith_specific_yield must lie in (0, regolith_theta_s - regolith_theta_r].');
end
if cfg.subsurface_drainage_substeps < 1 || cfg.subsurface_drainage_substeps ~= round(cfg.subsurface_drainage_substeps)
    error('HydroPol2D:Subsurface', 'subsurface_drainage_substeps must be a positive integer.');
end
cfg.deep_thicknesses_m = hp2d_deep_layer_thicknesses(cfg.subsurface_max_depth_m, ...
    cfg.subsurface_deep_first_thickness_m, cfg.subsurface_deep_growth_ratio, ...
    cfg.subsurface_deep_max_thickness_m);

cfg.feddes = [];
if feddes_on
    h3 = cfg.feddes_h3_m;
    F = struct('h1_m', cfg.feddes_h1_m, 'h2_m', cfg.feddes_h2_m, ...
        'h3_high_m', h3, 'h3_low_m', h3, 'h4_m', cfg.feddes_h4_m);
    if ~isempty(cfg.feddes_h3_high_m), F.h3_high_m = cfg.feddes_h3_high_m; end
    if ~isempty(cfg.feddes_h3_low_m), F.h3_low_m = cfg.feddes_h3_low_m; end
    wettest = max(F.h3_high_m, F.h3_low_m);
    driest = min(F.h3_high_m, F.h3_low_m);
    if ~(0 >= F.h1_m && F.h1_m > F.h2_m && F.h2_m >= wettest && driest > F.h4_m)
        error('HydroPol2D:Feddes', 'Feddes heads must satisfy 0 >= h1 > h2 >= h3 > h4 (m).');
    end
    cfg.feddes = F;
end
end
