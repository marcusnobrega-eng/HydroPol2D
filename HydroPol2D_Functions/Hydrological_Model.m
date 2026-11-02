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

% D8 local-inertial routing: the water in a river channel, up to its
% bankfull volume, neither infiltrates nor evaporates; the vertical modules
% see only the water above it (as in HydroPol2D-Python). It is set aside
% here and returned to the cell below.
use_d8_channel_storage = isfield(flags, 'flag_d8_local_inertial') && ...
    flags.flag_d8_local_inertial == 1 && exist('D8_Properties', 'var');
if use_d8_channel_storage
    d8_channel_mm = 1000 .* min(max(depths.d_t, 0) ./ 1000 .* Wshed_Properties.cell_area, ...
        D8_Properties.bankfull_volume_m3) ./ Wshed_Properties.cell_area;
    d8_channel_mm(~isfinite(d8_channel_mm)) = 0;
    d8_channel_mm = cast(d8_channel_mm, 'like', depths.d_t);
    depths.d_t = depths.d_t - d8_channel_mm;
end

% Green infrastructure covers a fraction f of some cells: the host soil,
% canopy and surface of those cells get the rain, ET and LAI over (1 - f), as
% in HydroPol2D-Python. The gross forcing is restored after the soil block.
use_green_infrastructure = isfield(flags, 'flag_green_infrastructure') && ...
    flags.flag_green_infrastructure == 1;
if use_green_infrastructure
    [BC_States, Hydro_States, LAI_raster, GI_Gross_Forcing] = ...
        hp2d_green_infrastructure_scale_forcing(BC_States, Hydro_States, LAI_raster, ...
        Green_Infrastructure.host_fraction_map, flags);
end

% Snow Module: splits the precipitation into rain and snowfall before the
% canopy. Snowfall goes to the snowpack and melt to the surface; only the
% rain reaches the canopy (as in HydroPol2D-Python).
Snow_Module

% Rain on a detention pond goes to its level pool (hp2d_detention_ponds_forcing_step
% before this script), not to the canopy and the soil of the pond cells.
use_detention_ponds = isfield(flags, 'flag_detention_ponds') && flags.flag_detention_ponds == 1;
if use_detention_ponds && ~(exist('detention_ponds_coupled', 'var') && detention_ponds_coupled)
    error('HydroPol2D:DetentionPonds', ...
        'Detention ponds are coupled in HydroPol2D_Main_While only; this time loop does not run them.');
end
if use_detention_ponds
    if isscalar(BC_States.Liquid_Precipitation)
        BC_States.Liquid_Precipitation = BC_States.Liquid_Precipitation .* ...
            ones(size(depths.d_t), 'like', depths.d_t);
    end
    BC_States.Liquid_Precipitation(Detention_Ponds.mask) = 0;
end

% Interception Module: intercepts the rain and adds the throughfall to the surface
interception_module

if isfield(flags, 'flag_layered_subsurface') && flags.flag_layered_subsurface == 1
    % Layered soil-over-regolith column: infiltration, ET (Feddes), drainage, the
    % water-table flux, the MMF table update, river exchange and lateral flow.
    Layered_Subsurface_Module
else
    % Infiltration Module
    Infiltration_Module

    % Evaporation / Evapotranspiration Module: after infiltration, so a cell
    % is wet only if water is left on the surface (as in HydroPol2D-Python and
    % the layered path). Open water evaporates; soil ET runs on dry cells.
    Evaporation_Evapotranspiration_Module

    % Recharge and Groundwater Module
    Groundwater_Module
end

% Green infrastructure: after the host cell's vertical step, before the
% inflow and the routing (HydroPol2D-Python _advance_green_infrastructure).
if use_green_infrastructure
    [BC_States, Hydro_States, LAI_raster, GI_Gross_Forcing] = ...
        hp2d_green_infrastructure_restore_forcing(BC_States, Hydro_States, LAI_raster, ...
        GI_Gross_Forcing);
    Green_Infrastructure_Module
end

% Inflow
depths.d_t = depths.d_t + BC_States.inflow;

% Water left on pond cells (groundwater exfiltration, an inflow cell) joins the pond.
if use_detention_ponds
    [Detention_Ponds, depths.d_t] = hp2d_detention_ponds_take_surface(Detention_Ponds, depths.d_t);
end

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

if use_d8_channel_storage
    depths.d_t = depths.d_t + d8_channel_mm;
end

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
