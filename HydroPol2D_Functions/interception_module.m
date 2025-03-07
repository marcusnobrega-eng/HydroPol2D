% Estimates interception, throughfall, and canopy evaporation
% Developer: Marcus Nobrega, Ph.D
% Date: 2/24/2025

if flags.flag_ETP == 1 && flags.flag_abstraction == 1
    [Hydro_States.S, Hydro_States.T, Hydro_States.E_int, Hydro_States.St_F] = interceptionModel(BC_States.delta_p_agg, Hydro_States.Ep*(time_step/60/24), LAI_raster.Z, Hydro_States.S, 0.2);
elseif flags.flag_ETP == 0 && ~isempty(LAI_raster.Z)
    [Hydro_States.S, Hydro_States.T, Hydro_States.E_int, Hydro_States.St_F] = interceptionModel(BC_States.delta_p_agg, zeros(size(elevation,1),size(elevation,2)), LAI_raster.Z, Hydro_States.S, 0.2);
else
    Hydro_States.T = BC_States.delta_p_agg; % [mm]
    Hydro_States.S = zeros(size(elevation,1),size(elevation,2));
    Hydro_States.E_int = zeros(size(elevation,1),size(elevation,2));
    Hydro_States.St_F = zeros(size(elevation,1),size(elevation,2)); 
end

% Actual Rain on the Grid that can be infiltrated
BC_States.Eff_Rainfall = Hydro_States.T + Hydro_States.St_F;
BC_States.Eff_Rainfall(idx_nan) = nan;

% Depth mass balance
depths.d_t = depths.d_p + BC_States.Eff_Rainfall; % mm
