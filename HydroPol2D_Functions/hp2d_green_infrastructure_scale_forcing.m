function [BC_States, Hydro_States, LAI_raster, Gross] = hp2d_green_infrastructure_scale_forcing( ...
    BC_States, Hydro_States, LAI_raster, host_fraction, flags)
%HP2D_GREEN_INFRASTRUCTURE_SCALE_FORCING Give the host part of GI cells its share of the forcing.
%   A green-infrastructure facility covers a fraction f of some cells. As in
%   HydroPol2D-Python (_advance_vertical_hydrology, host_area_fraction), the
%   host soil, canopy and surface of such a cell get the rain, the canopy and
%   open-water evaporation (Ep), the soil ET demand (ETP, or the ET maps) and
%   the LAI over (1 - f); the GI compartments get the rest
%   (hp2d_green_infrastructure_step). Gross keeps the unscaled fields so
%   hp2d_green_infrastructure_restore_forcing can put them back after the
%   vertical step (ETP and the ET maps are held between forcing records, so
%   they must not stay scaled).
Gross = struct('active', true);
Gross.delta_p_agg = BC_States.delta_p_agg;
BC_States.delta_p_agg = BC_States.delta_p_agg .* host_fraction;
for name = {'ETP', 'Ep'}
    if isfield(Hydro_States, name{1}) && ~isempty(Hydro_States.(name{1}))
        Gross.(name{1}) = Hydro_States.(name{1});
        Hydro_States.(name{1}) = Hydro_States.(name{1}) .* host_fraction;
    end
end
for name = {'current_evaporation_map', 'current_transpiration_map'}
    if isfield(BC_States, name{1}) && ~isempty(BC_States.(name{1}))
        Gross.(name{1}) = BC_States.(name{1});
        BC_States.(name{1}) = BC_States.(name{1}) .* host_fraction;
    end
end
if isfield(flags, 'flag_abstraction') && flags.flag_abstraction == 1 && ~isempty(LAI_raster)
    Gross.LAI = LAI_raster.Z;
    LAI_raster.Z = LAI_raster.Z .* host_fraction;
end
end
