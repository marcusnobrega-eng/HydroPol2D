function Detention_Ponds = hp2d_detention_ponds_write_outputs(Detention_Ponds, folder)
%HP2D_DETENTION_PONDS_WRITE_OUTPUTS Final pond depth map and the pond series as CSV.
%   Writes Detention_Ponds_Time_Series.csv (one row per hydrograph record) and
%   Detention_Ponds_Water_Balance.csv (cumulative volumes of each pond) to folder,
%   and keeps the final depth map in Detention_Ponds.depth_map_m.

P = Detention_Ponds;
P.depth_map_m = hp2d_detention_pond_depth_map(P);
if nargin < 2 || isempty(folder) || ~isfolder(folder)
    Detention_Ponds = P;
    return
end
H = P.history;
keep = H.time_min > 0 | (1:numel(H.time_min))' == 1;
T = table(H.time_min(keep), 'VariableNames', {'time_min'});
for p = 1:P.n_ponds
    name = matlab.lang.makeValidName(P.identifier{p});
    T.([name '_stage_m']) = H.stage_m(keep, p);
    T.([name '_storage_m3']) = H.storage_m3(keep, p);
    T.([name '_lateral_inflow_m3']) = H.lateral_inflow_m3(keep, p);
end
for d = 1:P.n_devices
    name = matlab.lang.makeValidName(strrep(P.device_identifier{d}, '/', '__'));
    T.([name '_discharge_m3_s']) = H.device_discharge_m3_s(keep, d);
    T.([name '_cumulative_m3']) = H.device_cumulative_m3(keep, d);
end
writetable(T, fullfile(folder, 'Detention_Ponds_Time_Series.csv'));
device_net = accumarray(P.dev_pond, P.cumulative_device_m3, [P.n_ponds 1]);
residual = P.storage_m3 - P.initial_storage_m3 - P.cumulative_rainfall_m3 - ...
    P.cumulative_lateral_m3 - P.cumulative_surface_in_m3 + P.cumulative_evaporation_m3 + ...
    P.cumulative_seepage_m3 + device_net;
B = table(string(P.identifier(:)), P.initial_storage_m3, P.storage_m3, ...
    P.cumulative_rainfall_m3, P.cumulative_lateral_m3, P.cumulative_surface_in_m3, ...
    P.cumulative_evaporation_m3, P.cumulative_seepage_m3, P.cumulative_seepage_out_m3, ...
    device_net, P.cumulative_external_m3, residual, ...
    'VariableNames', {'pond','initial_storage_m3','final_storage_m3','rainfall_m3', ...
    'critical_inflow_m3','surface_and_exfiltration_inflow_m3','evaporation_m3', ...
    'seepage_m3','seepage_out_of_model_m3','device_net_outflow_m3', ...
    'device_outflow_out_of_model_m3','residual_m3'});
writetable(B, fullfile(folder, 'Detention_Ponds_Water_Balance.csv'));
Detention_Ponds = P;
end
