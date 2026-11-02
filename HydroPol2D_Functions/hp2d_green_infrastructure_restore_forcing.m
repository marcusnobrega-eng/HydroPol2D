function [BC_States, Hydro_States, LAI_raster, Gross] = hp2d_green_infrastructure_restore_forcing( ...
    BC_States, Hydro_States, LAI_raster, Gross)
%HP2D_GREEN_INFRASTRUCTURE_RESTORE_FORCING Put back the forcing that
%   hp2d_green_infrastructure_scale_forcing scaled to the host fraction.
if ~isstruct(Gross) || ~isfield(Gross, 'active') || ~Gross.active
    return
end
BC_States.delta_p_agg = Gross.delta_p_agg;
for name = {'ETP', 'Ep'}
    if isfield(Gross, name{1}), Hydro_States.(name{1}) = Gross.(name{1}); end
end
for name = {'current_evaporation_map', 'current_transpiration_map'}
    if isfield(Gross, name{1}), BC_States.(name{1}) = Gross.(name{1}); end
end
if isfield(Gross, 'LAI'), LAI_raster.Z = Gross.LAI; end
Gross.active = false;
end
