function S = hp2d_green_infrastructure_initial_state(G)
%HP2D_GREEN_INFRASTRUCTURE_INITIAL_STATE Starting GI state (Python engine initial state):
%   empty surface and gravel, media at the design initial water content, no
%   lateral flow, zero cumulative volumes.
n = G.compartment_count;
z = zeros(n, 1);
S = struct();
S.surface_m3 = z;
S.media_m3 = G.initial_media_water_volume_m3(:);
S.gravel_m3 = z;
S.edge_q_m3_s = zeros(G.connection_count, 1);
S.gravel_edge_q_m3_s = zeros(G.connection_count, 1);
for name = {'rain', 'capture', 'infiltration', 'overflow', 'et', 'underdrain', 'native_exfiltration'}
    S.(['step_' name{1} '_m3']) = z;
    S.(['cumulative_' name{1} '_m3']) = z;
end
S.step_external_m3_host = zeros(numel(G.host_cells), 1);
S.step_total_external_m3 = 0;
S.step_total_et_m3 = 0;
S.step_dt_s = 0;
S.cumulative_external_m3 = 0;
S.step_log = {};
end
