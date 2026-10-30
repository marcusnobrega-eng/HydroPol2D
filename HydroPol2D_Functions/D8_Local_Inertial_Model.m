function [d_t, outlet_flow, outflow, Hf, qout_left, qout_right, qout_up, qout_down, D8_States] = ...
    D8_Local_Inertial_Model(d_tot, outflow_prev, D8, D8_States, time_step, idx_nan, ...
    flag_critical, flag_reservoir, Reservoir_Data)
%D8_LOCAL_INERTIAL_MODEL One routing step of the D8 local-inertial scheme.
%   The Wflow-style D8 scheme of HydroPol2D-Python
%   (routing/d8_river_local_inertial.py, flags.routing = d8_local_inertial
%   or full_d8_local_inertial). Water is held as cell volume; d_tot and d_t
%   are that volume over the cell area [mm], as in the D4 solvers.
%
%   Hybrid mode (D8.full_domain = 0), in this order:
%     1. D4 floodplain faces (hp2d_d8_floodplain_d4) with the D8 stage, the
%        faces of cardinal river links closed;
%     2. the D8 river links on the resulting storage (hp2d_d8_river_links);
%   Full-domain mode: step 2 only, on links between every valid cell and its
%   receiver. Then, in both modes,
%     3. reservoir controls (hp2d_reservoir_releases), the depth above the
%        channel bed from the cell tables;
%     4. the outlets at the terminal D8 nodes (hp2d_d8_outlet_flow). Their
%        outflow, and any reservoir release that leaves the model, is the
%        outlet flow [mm/h over the cell].
%
%   outflow_prev / outflow: pages 1-2 are the D4 floodplain face flows [mm/h]
%   (Local_Inertial_Model_D4 convention: east and north faces), page 3 the
%   outlet flow. Hf: face depths [m] (pages 1-2) and outlet depth (page 3).
%   D8_States.discharge_m3_s holds the link discharges between steps.

dt = time_step * 60; % [s]
A = D8.cell_area;
[ny, nx] = deal(D8.ny, D8.nx);
template = d_tot;

V = double(gather_if(d_tot(:))) / 1000 * A;
V(~isfinite(V)) = 0;
V = max(V, 0);
V(~D8.valid_grid(:)) = 0;

%% 1. Floodplain (hybrid mode)
qx = zeros(ny, nx);
qy = zeros(ny, nx);
hx = zeros(ny, nx);
hy = zeros(ny, nx);
if ~D8.full_domain
    qx_prev = double(gather_if(outflow_prev(:,:,1))) * A / 1000 / 3600; % [m3/s]
    qy_prev = double(gather_if(outflow_prev(:,:,2))) * A / 1000 / 3600;
    qx_prev(~isfinite(qx_prev)) = 0;
    qy_prev(~isfinite(qy_prev)) = 0;
    stage = hp2d_d8_cell_stage(D8.Cells, V);
    [V, qx, qy, hx, hy] = hp2d_d8_floodplain_d4(V, qx_prev, qy_prev, stage, D8, dt, flag_critical);
end

%% 2. D8 links
[V, D8_States.discharge_m3_s, cardinal_m3, diagonal_m3] = ...
    hp2d_d8_river_links(V, D8_States.discharge_m3_s, D8, dt);
D8_States.cumulative_cardinal_transfer_m3 = D8_States.cumulative_cardinal_transfer_m3 + cardinal_m3;
D8_States.cumulative_diagonal_transfer_m3 = D8_States.cumulative_diagonal_transfer_m3 + diagonal_m3;

%% 3. Reservoir controls
reservoir_external_m3 = zeros(ny, nx);
if flag_reservoir == 1
    geometry = struct('type', 'handles', ...
        'depth', @(volume, idx) cell_depth(D8, volume, idx), ...
        'volume', @(depth, idx) cell_volume(D8, depth, idx));
    [exchange_m3, reservoir_external_m3] = hp2d_reservoir_releases( ...
        reshape(V, ny, nx), Reservoir_Data.x_index, Reservoir_Data.y_index, ...
        Reservoir_Data.k1, Reservoir_Data.h1, Reservoir_Data.k2, Reservoir_Data.k3, ...
        Reservoir_Data.h2, Reservoir_Data.k4, Reservoir_Data.y_ds1_index, ...
        Reservoir_Data.x_ds1_index, Reservoir_Data.y_ds2_index, Reservoir_Data.x_ds2_index, ...
        dt, geometry);
    V = max(V + exchange_m3(:), 0);
end

%% 4. Outlets
[V, outlet_m3_s, outlet_depth_m] = hp2d_d8_outlet_flow(V, D8, dt);
outlet_mm_h = reshape(outlet_m3_s, ny, nx) / A * 1000 * 3600;
outlet_mm_h = hp2d_reservoir_book_external(outlet_mm_h, reservoir_external_m3, A, dt);
D8_States.last_outlet_m3_s = sum(outlet_m3_s);

%% Outputs in the conventions of the D4 solvers
d_t = reshape(V, ny, nx) / A * 1000;
d_t(idx_nan) = NaN;
outflow_mm_h = zeros(ny, nx, 3);
outflow_mm_h(:,:,1) = qx / A * 1000 * 3600;
outflow_mm_h(:,:,2) = qy / A * 1000 * 3600;
outflow_mm_h(:,:,3) = outlet_mm_h;
Hf_m = zeros(ny, nx, 3);
Hf_m(:,:,1) = hx;
Hf_m(:,:,2) = hy;
Hf_m(:,:,3) = reshape(outlet_depth_m, ny, nx);

d_t = cast_like(d_t, template);
outlet_flow = cast_like(outlet_mm_h, template);
outflow = cast_like(outflow_mm_h, template);
Hf = cast_like(Hf_m, template);
qout_right = outflow(:,:,1);
qout_left = -[zeros(ny, 1, 'like', outflow), outflow(:, 1:end-1, 1)];
qout_up = outflow(:,:,2);
qout_down = -[outflow(2:end, :, 2); zeros(1, nx, 'like', outflow)];

end

function h = cell_depth(D8, volume, idx)
% Depth above the cell invert: the channel bed on table cells.
if D8.Cells.is_table(idx)
    h = max(hp2d_d8_cell_stage(D8.Cells, volume, idx) - D8.Cells.datum(idx), 0);
else
    h = volume / D8.cell_area;
end
end

function volume = cell_volume(D8, depth, idx)
% Volume held up to a depth above the cell invert.
if D8.Cells.is_table(idx)
    volume = hp2d_d8_cell_volume(D8.Cells, D8.Cells.datum(idx) + depth, idx);
else
    volume = depth * D8.cell_area;
end
end

function x = gather_if(x)
if isa(x, 'gpuArray')
    x = gather(x);
end
end

function y = cast_like(x, template)
if isa(template, 'gpuArray')
    y = gpuArray(cast(x, classUnderlying(template)));
else
    y = cast(x, class(template));
end
end
