% Infiltration Routine
% Developer: Marcus Nobrega, Ph.D.
% Goal: Estimate effective precipitation and infiltration at cells domain
% Date: 2/22/2025
% Added recharge calculation. Added 2D boussinesq solution to the
% groundwater model

use_neal_subgrid_volume = flags.flag_subgrid == 1 && flags.flag_overbanks == 1;
use_lookup_subgrid_volume = flags.flag_subgrid == 1 && flags.flag_overbanks ~= 1 ...
    && exist('SubgridTables', 'var') && ~isempty(SubgridTables);

if use_neal_subgrid_volume
    neal_d_p_representative = depths.d_p;
    neal_volume_t = hp2d_neal_cell_volume( ...
        max(depths.d_t ./ 1000, 0), Wshed_Properties.River_Width, ...
        Wshed_Properties.River_Depth, Wshed_Properties.Resolution);
    neal_volume_p = hp2d_neal_cell_volume( ...
        max(depths.d_p ./ 1000, 0), Wshed_Properties.River_Width, ...
        Wshed_Properties.River_Depth, Wshed_Properties.Resolution);

    % Hydrologic fluxes are areal depths over the full coarse cell.
    depths.d_t = 1000 .* neal_volume_t ./ Wshed_Properties.cell_area;
    depths.d_p = 1000 .* neal_volume_p ./ Wshed_Properties.cell_area;
    C_a = Wshed_Properties.cell_area .* ones(size(depths.d_t), 'like', depths.d_t);
    C_a(idx_nan) = NaN;
elseif use_lookup_subgrid_volume
    d_t_before_hydrology = depths.d_t;
    C_a = Wshed_Properties.cell_area .* ones(size(depths.d_t), 'like', depths.d_t);
    C_a(idx_nan) = NaN;
end

% Interception Module
interception_module

% Snow Module
Snow_Module

% Evaporation / Evapotranspiration Module
Evaporation_Evapotranspiration_Module

% Infiltration Module
Infiltration_Module

% Recharge and Groundwater Module
Groundwater_Module

% Inflow
depths.d_t = depths.d_t + BC_States.inflow;

% Lookup-table subgrid stores water by volume, not by representative area.
% Hydrologic source/sink modules update depth in mm, so convert the net
% areal change to volume over the full coarse cell and invert the storage
% curve back to representative depth.
if use_lookup_subgrid_volume
    delta_depth_mm = depths.d_t - d_t_before_hydrology;
    delta_depth_mm(~isfinite(delta_depth_mm)) = 0;
    depths.d_t = 1000 .* hp2d_subgrid_apply_volume_change( ...
        max(d_t_before_hydrology ./ 1000, 0), delta_depth_mm, ...
        SubgridTables, Wshed_Properties.cell_area);
end

% Depths
if min(min(depths.d_t)) < -1e-8
    catch_index = catch_index + 1;
    warning('Negative depths. Please reduce the time-step.')
else
    depths.d_t(depths.d_t < 1e-6) = 0;
end

depths.d_t(depths.d_t < 0) = 0;

if use_neal_subgrid_volume
    neal_volume_t = max(depths.d_t, 0) ./ 1000 .* Wshed_Properties.cell_area;
    depths.d_t = 1000 .* hp2d_neal_depth_from_volume( ...
        neal_volume_t, Wshed_Properties.River_Width, ...
        Wshed_Properties.River_Depth, Wshed_Properties.Resolution);
    depths.d_p = neal_d_p_representative;
    C_a = hp2d_neal_cell_area( ...
        depths.d_t ./ 1000, Wshed_Properties.River_Width, ...
        Wshed_Properties.River_Depth, Wshed_Properties.Resolution);
    C_a(idx_nan) = NaN;
end

% Effective precipitation available to the routing model.
depths.d_tot = depths.d_t;
