function F = hp2d_urban_drainage_node_fluxes(N, pipe_q_m3_s, inlet_q_m3_s)
%HP2D_URBAN_DRAINAGE_NODE_FLUXES Positive inflow and outflow components of
%   each node [m3/s] from the signed pipe and inlet flows
%   (urban_drainage_node_fluxes in HydroPol2D-Python).

n = numel(N.node_ids);
q = pipe_q_m3_s(:);
qi = inlet_q_m3_s(:);
pipe_in = accumarray(N.pipe_start, max(-q, 0), [n 1]) + accumarray(N.pipe_end, max(q, 0), [n 1]);
pipe_out = accumarray(N.pipe_start, max(q, 0), [n 1]) + accumarray(N.pipe_end, max(-q, 0), [n 1]);
if isempty(qi)
    capture = zeros(n, 1);
    surcharge = zeros(n, 1);
else
    capture = accumarray(N.inlet_node, max(qi, 0), [n 1]);
    surcharge = accumarray(N.inlet_node, max(-qi, 0), [n 1]);
end
internal_in = pipe_in + capture;
internal_out = pipe_out + surcharge;
boundary_net = internal_in - internal_out;
outfall_out = zeros(n, 1);
outfall_in = zeros(n, 1);
outfall_out(~N.dynamic) = max(boundary_net(~N.dynamic), 0);
outfall_in(~N.dynamic) = max(-boundary_net(~N.dynamic), 0);
F = struct('total_inflow_m3_s', internal_in + outfall_in, ...
    'total_outflow_m3_s', internal_out + outfall_out, ...
    'surface_capture_m3_s', capture, 'surface_surcharge_m3_s', surcharge, ...
    'outfall_inflow_m3_s', outfall_in, 'outfall_outflow_m3_s', outfall_out);
end
