function Detention_Ponds = hp2d_detention_ponds_record(Detention_Ponds, record, time_min)
%HP2D_DETENTION_PONDS_RECORD Store the pond series at hydrograph record "record".
%   Stage [m], storage [m3], device discharge [m3/s] and cumulative volume [m3], and
%   cumulative critical-capture inflow [m3], at the outlet-hydrograph times (the
%   detention outputs of HydroPol2D-Python).

P = Detention_Ponds;
H = P.history;
H.time_min(record, 1) = double(gather(time_min));
H.stage_m(record, :) = hp2d_detention_pond_stage(P.storage_m3, P)';
H.storage_m3(record, :) = P.storage_m3';
H.device_discharge_m3_s(record, :) = P.device_discharge_m3_s';
H.device_cumulative_m3(record, :) = P.cumulative_device_m3';
H.lateral_inflow_m3(record, :) = P.cumulative_lateral_m3';
P.history = H;
Detention_Ponds = P;
end
