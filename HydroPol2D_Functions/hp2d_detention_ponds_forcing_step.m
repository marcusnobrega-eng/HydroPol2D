function Detention_Ponds = hp2d_detention_ponds_forcing_step(Detention_Ponds, delta_p_agg, Ep_mm_day, flags, time_step_min)
%HP2D_DETENTION_PONDS_FORCING_STEP Rain, evaporation and seepage of the ponds for one step.
%   Runs before Hydrological_Model, as HydroPol2D-Python applies the pond forcing
%   before the vertical step. delta_p_agg: rain depth of this step [mm] (scalar or
%   grid; the pond cells' rain goes to the pond, not to the canopy or the soil);
%   Ep_mm_day: open-water evaporation (Hydro_States.Ep, used when flag_ETP = 1).
%   Seepage is internal recharge: with groundwater on, Detention_Ponds.recharge_rate_m_s
%   [m/s] goes to the aquifer under the pond (Groundwater_Module or
%   Layered_Subsurface_Module); with groundwater off it leaves the model and is booked
%   in last_seepage_out_m3 (Python's pond_seepage_m3).

P = Detention_Ponds;
dt_s = time_step_min * 60;
n_cells = numel(P.cells);
if isscalar(delta_p_agg)
    rain_mm = double(gather(delta_p_agg)) * ones(n_cells, 1);
else
    rain_mm = double(gather(delta_p_agg(P.cells)));
end
if flags.flag_ETP == 1 && ~isempty(Ep_mm_day)
    if isscalar(Ep_mm_day)
        evaporation = double(gather(Ep_mm_day)) * ones(n_cells, 1);
    else
        evaporation = double(gather(Ep_mm_day(P.cells)));
    end
    evaporation(~isfinite(evaporation)) = 0;
else
    evaporation = zeros(n_cells, 1);
end
rain_mm(~isfinite(rain_mm)) = 0;

F = hp2d_detention_pond_forcing(P.storage_m3, rain_mm * 3600 / dt_s, evaporation, P, dt_s);
P.storage_m3 = F.storage_m3;
P.cumulative_rainfall_m3 = P.cumulative_rainfall_m3 + F.rainfall_m3;
P.cumulative_evaporation_m3 = P.cumulative_evaporation_m3 + F.evaporation_m3;
P.cumulative_seepage_m3 = P.cumulative_seepage_m3 + F.seepage_m3;
P.last_evaporation_m3 = sum(F.evaporation_m3);
P.recharge_rate_m_s(:) = 0;
if flags.flag_groundwater_modeling == 1
    P.recharge_rate_m_s(P.cells) = F.recharge_rate_m_s;
    P.last_seepage_out_m3 = 0;
else
    P.cumulative_seepage_out_m3 = P.cumulative_seepage_out_m3 + F.seepage_m3;
    P.last_seepage_out_m3 = sum(F.seepage_m3);
end
Detention_Ponds = P;
end
