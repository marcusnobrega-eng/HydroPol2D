function L = hp2d_voronoi_layered_initialize(mesh, state, groundwater, cfg, hydrology)
%HP2D_VORONOI_LAYERED_INITIALIZE One layered column per mesh HRU (cell x HRU pairs).
%   Columns are the (cell, HRU) pairs of the HRU matrices in MATLAB linear order, with
%   host = cell index and area = HRU area. The HRUs of a cell share its ground
%   (surface_bed), column base (surface_bed - groundwater bottom) and head, as in
%   HydroPol2D-Python voronoi_layered.prepare_layered_mesh. Optional hydrology fields:
%   soil_base_depth_m and regolith_* ([n x 1] per cell, or scalars), ksat_root_m_s and
%   ksat_transmission_m_s ([n x h]); without them the cfg scalars and ksat_m_s apply.
if nargin < 5, hydrology = struct(); end
n = mesh.n_cells;
h = size(state.fraction, 2);
N = n*h;
host = repmat((1:n)', h, 1);
    function v = per_column(value)
        value = double(value);
        if isscalar(value)
            v = value + zeros(N, 1);
        elseif numel(value) == n && (h == 1 || iscolumn(value))
            v = value(host);
        elseif numel(value) == h && isrow(value)
            v = reshape(repmat(value, n, 1), [], 1);  % one value per HRU type
        else
            v = value(:);
        end
    end
    function v = cell_value(name, fallback)
        if isstruct(hydrology) && isfield(hydrology, name) && ~isempty(hydrology.(name))
            v = per_column(hydrology.(name));
            missing = ~isfinite(v);
            fb = per_column(fallback);
            v(missing) = fb(missing);
        else
            v = per_column(fallback);
        end
    end
surface_bed = mesh.surface_bed(:);
% Column bottom: the aquifer bottom, capped at the layered maximum depth (as on the grid).
base_cell = min(max(surface_bed - groundwater.bottom(:), 0), cfg.subsurface_max_depth_m);
base = base_cell(host);
soil_base = min(cell_value('soil_base_depth_m', cfg.subsurface_soil_base_depth_m), base);
root = state.root_depth_m(:);
near = state.surface_layer_depth_m(:);
[z_top, z_bot] = hp2d_slot_geometry(root, soil_base, base, near, cfg.deep_thicknesses_m);
ks_near = state.ksat_m_s(:);
ks_root = cell_value('ksat_root_m_s', ks_near);
ks_trans = cell_value('ksat_transmission_m_s', ks_near);
P.soil_theta_s = state.theta_sat(:);
P.soil_theta_r = state.theta_r(:);
P.soil_alpha_m_inverse = state.alpha_vg_1_m(:);
P.soil_n_vg = state.n_vg(:);
P.soil_kv_m_s = [ks_near, ks_root, ks_trans];
P.soil_lateral_anisotropy = cfg.soil_lateral_anisotropy;
P.regolith_theta_s = cell_value('regolith_theta_s', cfg.regolith_theta_s);
P.regolith_theta_r = cell_value('regolith_theta_r', cfg.regolith_theta_r);
P.regolith_alpha_m_inverse = cell_value('regolith_alpha_vg_m_inverse', cfg.regolith_alpha_vg_m_inverse);
P.regolith_n_vg = cell_value('regolith_n_vg', cfg.regolith_n_vg);
P.regolith_kv_m_s = cell_value('regolith_ksat_v_mm_h', cfg.regolith_ksat_v_mm_h)/1000/3600;
kh_fallback = per_column(groundwater.hydraulic_conductivity);
if ~isempty(cfg.regolith_ksat_h_mm_h), kh_fallback(:) = cfg.regolith_ksat_h_mm_h/1000/3600; end
if isfield(hydrology, 'regolith_ksat_h_mm_h') && ~isempty(hydrology.regolith_ksat_h_mm_h)
    kh = per_column(hydrology.regolith_ksat_h_mm_h)/1000/3600;
    kh(~isfinite(kh)) = kh_fallback(~isfinite(kh));
else
    kh = kh_fallback;
end
P.regolith_kh_m_s = kh;
sy_fallback = nan(N, 1);
if ~isempty(cfg.regolith_specific_yield), sy_fallback(:) = cfg.regolith_specific_yield; end
P.regolith_specific_yield = cell_value('regolith_specific_yield', sy_fallback);
P.retention_head_m = cfg.capillary_field_capacity_head_m;
P.ground_m = surface_bed(host);
P.column_base_depth_m = base;
P.valid = true(N, 1);
sy = P.regolith_specific_yield;
bad = isfinite(sy) & ~(sy > 0 & sy <= P.regolith_theta_s - P.regolith_theta_r + 1e-12);
if any(bad)
    error('HydroPol2D:Subsurface', ...
        'Regolith specific yield lies outside (0, theta_s - theta_r] in %d mesh column(s).', nnz(bad));
end
column = hp2d_slot_column(z_top, z_bot, P);
head = min(max(groundwater.head(:), surface_bed - base_cell), surface_bed);
saturation = 0.5;
if isfield(hydrology, 'initial_soil_saturation') && ~isempty(hydrology.initial_soil_saturation)
    saturation = per_column(hydrology.initial_soil_saturation);
end
theta0 = P.soil_theta_r + saturation.*(P.soil_theta_s - P.soil_theta_r);
L.column = column;
L.water_mm = hp2d_initial_slot_water(column, head(host), repmat(theta0, 1, 3));
L.head_m = head;
L.host = host;
L.area_m2 = state.area_m2(:);
L.ground_m = surface_bed;
L.base_depth_m = base_cell;
L.n_cells = n;
L.n_hru = h;
L.n_slots = size(z_top, 2);
L.drainage_substeps = cfg.subsurface_drainage_substeps;
L.feddes = cfg.feddes;
L.config = cfg;
end
