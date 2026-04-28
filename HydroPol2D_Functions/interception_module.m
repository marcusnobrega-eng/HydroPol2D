%% ═══════════════════════════════════════════════════════════════════════
%  Module: interceptionAndThroughfall
%  🛠️ Developer: Marcus Nobrega, Ph.D.
%  📅 Date: 02/24/2025
% ─────────────────────────────────────────────────────────────────────────────
%  ➤ Purpose:
%      Estimate interception, throughfall, and canopy evaporation for the
%      hydrological model. The module adjusts rainfall based on subgrid and
%      overbank conditions, applies an interception model, and computes the 
%      effective rainfall available for infiltration.
%
%  ➤ Inputs:
%      • flags: Structure with boolean flags controlling processes such as
%               subgrid modeling (flag_subgrid), overbanks (flag_overbanks),
%               ETP (flag_ETP), and abstraction (flag_abstraction).
%      • BC_States: Structure containing boundary condition states including
%                   aggregated rainfall (delta_p_agg).
%      • Wshed_Properties: Watershed properties, including spatial resolution.
%      • C_a: Cell area used for conversion in subgrid calculations.
%      • Hydro_States: Structure holding hydrological states such as interception 
%                      storage (S), throughfall (T), intercepted evaporation (E_int),
%                      and storage fraction (St_F).
%      • LAI_raster: Raster containing Leaf Area Index data (Z).
%      • elevation: Matrix of elevation values.
%      • time_step: Duration of the current time step (in minutes).
%      • DEM_raster: Digital Elevation Model raster structure (used for cellsize).
%      • idx_nan: Logical mask indicating cells with invalid/missing data.
%
%  ➤ Outputs:
%      • Hydro_States.T: Throughfall (mm) – portion of rainfall not intercepted.
%      • Hydro_States.S: Updated canopy interception storage (mm).
%      • Hydro_States.E_int: Canopy evaporation (mm).
%      • Hydro_States.St_F: Fraction of storage released as throughfall (mm).
%      • BC_States.Eff_Rainfall: Effective rainfall on the grid that can be infiltrated (mm).
%      • errors(1): Interception mass balance error [m³] computed as error*DEM_raster.cellsize^2.
%
%  ➤ Notes:
%      • When subgrid and overbank conditions are active, the interception is
%        adjusted by scaling the aggregated rainfall to the grid cell area.
%      • If ETP and abstraction are enabled, the interception model is applied 
%        to partition rainfall between storage, throughfall, and evaporation.
%      • In the absence of these flags, the model defaults to assigning the full
%        rainfall as throughfall with zero interception and evaporation.
%      • The effective rainfall is masked by idx_nan to propagate missing data.
% ═══════════════════════════════════════════════════════════════════════


if flags.flag_subgrid == 1 && flags.flag_overbanks == 1
    P_interception = BC_States.delta_p_agg .* Wshed_Properties.Resolution^2 ./ C_a; % mm
else
    P_interception = BC_States.delta_p_agg; % mm
end


if flags.flag_ETP == 1 && flags.flag_abstraction == 1
    [Hydro_States.S, Hydro_States.T, Hydro_States.E_int, Hydro_States.St_F, error] = interceptionModel(P_interception, Hydro_States.Ep*(time_step/60/24), LAI_raster.Z, Hydro_States.S, 0.2, zero_matrix);
elseif flags.flag_ETP == 0 && flags.flag_abstraction == 1
    [Hydro_States.S, Hydro_States.T, Hydro_States.E_int, Hydro_States.St_F, error] = interceptionModel(P_interception, zeros(size(elevation,1),size(elevation,2)), LAI_raster.Z, Hydro_States.S, 0.2, zero_matrix);
else
    Hydro_States.T = BC_States.delta_p_agg; % [mm]
    if k == 1
        Hydro_States.S = zeros(size(elevation,1),size(elevation,2));
        Hydro_States.E_int = zeros(size(elevation,1),size(elevation,2));
        Hydro_States.St_F = zeros(size(elevation,1),size(elevation,2)); 
    end
    error = 0;
end

% Actual Rain on the Grid that can be infiltrated
BC_States.Eff_Rainfall = Hydro_States.T + Hydro_States.St_F;
BC_States.Eff_Rainfall(idx_nan) = nan;

errors(1) = error*DEM_raster.cellsize^2; % Interception mass balance error [m3]
% 
% % Depth mass balance
% depths.d_t = depths.d_p + BC_States.Eff_Rainfall; % mm
