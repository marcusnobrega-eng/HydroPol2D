%@ -0,0 +1,242 @@
function [qout_left,qout_right,qout_up,qout_down,outlet_flow,qout_ne,qout_se,qout_sw,qout_nw,d_t,I_tot_end_cell,outflow,Hf] = Bates_Inertial_8D(reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2,flag_reservoir,z,d_tot,d_p,roughness_cell,cell_area,time_step,Resolution,outlet_index,outlet_type,slope_outlet,row_outlet,col_outlet,d_tolerance,outflow,idx_nan,flag_critical)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                 %
%                 Produced by Marcus Nobrega Gomes Junior         %
%                 e-mail:marcusnobrega.engcivil@gmail.com         %
%                           September 2021                        %
%                                                                 %
%                 Last Updated: 7 July, 2022                      %
%                                                                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% --------------------- ATTENTION ------------------------------- %

% Matrix matrix_store is actually used to store all arrays, not necessarily
% delta_h

%§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
%   -----  DESCRIPTION  -----
%   This function estimates the transfered volume from a cell to the
%   neighbour cells using a weigthing approach in terms of the available
%   volume in the neighbour cells
%§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
% ---------------% Finding Nans or Values Below a Threshold  % ---------------%
if isgpuarray(cell_area)
    d_t_min = gpuArray(d_tolerance/1000); % m
    ny = gpuArray(size(z,1));
    nx = gpuArray(size(z,2));
else
    d_t_min = d_tolerance/1000; % m
    ny = size(z,1);
    nx = size(z,2);  
end
% ---------------% Adding minimum slope to do calculations % ---------------%
% h_min = d_t_min;% m
h_min = 0;

% --------------- Notation  % ---------------%
%   <-I-> (left right up down) = [left ; right ; up; down]
% --------------- Cell Depth and Water surface Elevation  % ---------------%
depth_cell = max(d_p./1000,0); % meters (fixing zero values)

% Assuming a linearization for small depths
mask = depth_cell <= h_min;
depth_cell(mask) = 0;

% Chaning depth_cell of B.C cells to NaN
for i = 1:length(reservoir_y)
    depth_cell(reservoir_y(i),reservoir_x(i)) = NaN;
end

y = z + depth_cell;

%% --------------- Depth Differences (delta_h) for All Cells  % ---------------%
matrix_store(:,:,1) = [y(:,1), y(:,2:(nx)) - y(:,1:(nx-1))]; % left
matrix_store(:,:,2) = [y(:,1:(nx-1)) - y(:,2:(nx)), y(:,end)]; % right
matrix_store(:,:,3) = [y(1,:); y(2:(ny),:) - y(1:(ny-1),:)]; % up
matrix_store(:,:,4) = [y(1:(ny-1),:) - y(2:(ny),:); y(end,:)]; % down

% --------------- Depth Differences for All NE, SE, SW, NW  % ---------------%
% 6, 7, 8, 9, respectively
matrix_store(:,:,6) =  [y(1,:); [y(2:(ny),1:(nx-1)) - y(1:(ny-1),2:(nx)), y(1:(ny-1),end)]]; % NE
matrix_store(:,:,7) =  [[y(1:(ny-1),1:(nx-1)) - y(2:ny,2:(nx)), y(2:(ny),end)]; y(end,:)]; % SE
matrix_store(:,:,8) = [[y(2:(ny),1), y(1:(ny-1),2:(nx)) - y(2:(ny),1:(nx-1))]; y(end,:)]; % SW
matrix_store(:,:,9) = [y(1,:);[y(1:(ny-1),1), y(2:(ny),2:(nx)) - y(1:(ny-1),1:(nx-1))]]; % NW

%% ---------------- Hf (Effective Depth for Calculations) ----------- %
% max(wse,wsei) - max(z,zi)
Hf = 0*matrix_store;
Hf(:,:,1) = [y(:,1) - z(:,1), max(y(:,2:(nx)), y(:,1:(nx-1))) - max(z(:,2:(nx)), z(:,1:(nx-1)))]; % left
Hf(:,:,2) = [max(y(:,1:(nx-1)), y(:,2:(nx))) - max(z(:,1:(nx-1)), z(:,2:(nx))), y(:,end) - z(:,end)]; % right
Hf(:,:,3) = [y(1,:) - z(1,:); max(y(2:(ny),:), y(1:(ny-1),:)) - max(z(2:(ny),:), z(1:(ny-1),:))]; % up
Hf(:,:,4) = [max(y(1:(ny-1),:), y(2:(ny),:)) - max(z(1:(ny-1),:), z(2:(ny),:)); y(end,:) - z(end,:)]; % down
Hf(:,:,5) = y - z;
% --------------- Depth Differences for All NE, SE, SW, NW  % ---------------%
% 6, 7, 8, 9, respectively
Hf(:,:,6) =  [y(1,:) - z(1,:); [max(y(2:(ny),1:(nx-1)), y(1:(ny-1),2:(nx))) - max(z(2:(ny),1:(nx-1)), z(1:(ny-1),2:(nx))), y(1:(ny-1),end) - z(1:(ny-1),end)]]; % NE
Hf(:,:,7) =  [[max(y(1:(ny-1),1:(nx-1)), y(2:ny,2:(nx))) - max(z(1:(ny-1),1:(nx-1)), z(2:ny,2:(nx))), y(2:(ny),end) - z(2:(ny),end)]; y(end,:) - z(end,:)]; % SE
Hf(:,:,8) = [[y(2:(ny),1) - z(2:(ny),1), max(y(1:(ny-1),2:(nx)) , y(2:(ny),1:(nx-1))) - max(z(1:(ny-1),2:(nx)) , z(2:(ny),1:(nx-1)))]; y(end,:) - z(end,:)]; % SW
Hf(:,:,9) = [y(1,:) - z(1,:);[y(1:(ny-1),1) - z(1:(ny-1),1), max(y(2:(ny),2:(nx)), y(1:(ny-1),1:(nx-1))) - max(z(2:(ny),2:(nx)), z(1:(ny-1),1:(nx-1)))]]; % NW


%%%%% CHECKED 9/28/2021 %%%%%
% --------------- Outlet Calculations  % ---------------%
if outlet_type == 1
    S_0 = slope_outlet; % Normal slope
    % V = h * area
    matrix_store(:,:,5) = (S_0*Resolution*outlet_index); % Normal slope depth difference
    %%%%% CHECKED 9/28/2021 %%%%%
else
    S_0 = zeros(size(z));
    for i = 1:length(row_outlet)
        % Checking left, right up, and down
        row = row_outlet(i); % Row of outlet
        col = col_outlet(i); % Col of outlet
        S_0(row,col) = (depth_cell(row,col).^(-1/6)).*sqrt(9.81).*roughness_cell(row,col); % Critical Slope
    end
    matrix_store(:,:,5) = (S_0*Resolution);
end

% ---------------% Correcting delta_h  % ---------------%
% delta_h(logical(isnan(delta_h) + isinf(delta_h)) ) = 0; % TESTING
% delta_h(isnan(delta_h)) = 0; % TESTING

% ---------------% Boundary Conditions  % ---------------%
% These only applied when we are modeling a rectangular grid watershed
[matrix_store(:,1,1),matrix_store(:,end,2),matrix_store(1,:,3),matrix_store(end,:,4)] = deal(0);
[matrix_store(:,end,6),matrix_store(:,end,7),matrix_store(end,:,8),matrix_store(:,1,9)] = deal(0);
[matrix_store(1,:,6),matrix_store(end,:,7),matrix_store(:,1,8),matrix_store(1,:,9)] = deal(0);

if flag_reservoir == 1
    for ii = 1:length(reservoir_y)
        matrix_store(reservoir_y(ii),reservoir_x(ii),:) = 0; % Imposing no flow to all directions
    end
end

matrix_store(logical(matrix_store < h_min | isnan(matrix_store))) = 0;
% ---------------% Available Volumes  % ---------------%
%slope
matrix_store(:,:,1:4) = matrix_store(:,:,1:4)/Resolution; % Calculating the slopes
matrix_store(:,:,5) = matrix_store(:,:,5)/Resolution;
matrix_store(:,:,6:end) = matrix_store(:,:,6:end)/(sqrt(2)*Resolution);

% No negative gradients
matrix_store(matrix_store<0) = 0;

% Matrix_store is the wse gradient from cell i to cell i +1
% such that matrix_store = (wse(i) - wse(i+1) ) / dx
% by definition it should be the opposite in Bates Formulation, so we multiply by -1
matrix_store = -matrix_store;
%% Bates Formulation
outflow = outflow/1000/3600*Resolution^2/Resolution; % m2 per sec. It comes from mm/h
dt = time_step*60; % time-step in seconds
% q_t_dt = [q_t - g h dt (wse grad)] / (1 + g h dt n^2 q_t / h^(10/3))
% q_t is in m2/s, all other variables are in international units
g = 9.81;
outflow_prev = outflow;
outflow = (outflow_prev - g*Hf*dt.*matrix_store)./ ...
          (1 + g*Hf*dt.*roughness_cell.^2.*outflow_prev./(Hf.^(10/3))); % m2 per sec (Hf can be simplified)

%% Eliminating Surplus Velocities at Wet-Dry Interfaces
% --------------- All NE, SE, SW, NW  % ---------------%
% 6, 7, 8, 9, respectively
% Left 
% artificial_depth = 0;
% idx = Hf(:,:,1) > artificial_depth & [Hf(:,1:(end-1),1),zeros(size(depth_cell,1),1)] == artificial_depth;
% interface_flow = outflow(:,:,1);
% depth = Hf(:,:,1);
% if any(any(idx)) 
%     interface_flow(idx) = depth(idx).*sqrt(g*depth(idx)); % m2 per sec
%     outflow(:,:,1) = interface_flow;   
% end
% % Right
% idx = Hf(:,:,2) > artificial_depth & [zeros(size(depth_cell,1),1), Hf(:,2:(end),2)] == artificial_depth;
% interface_flow = outflow(:,:,2);
% depth = Hf(:,:,2);
% if any(any(idx))
%     interface_flow(idx) = depth(idx).*sqrt(g*depth(idx));
%     outflow(:,:,2) = interface_flow;
% end
% % Up
% idx = Hf(:,:,3) > artificial_depth & [zeros(1,size(depth_cell,2));Hf(1:(end-1),:,3)] == artificial_depth;
% interface_flow = outflow(:,:,3);
% depth = Hf(:,:,3);
% if any(any(idx))
%     interface_flow(idx) = depth(idx).*sqrt(g*depth(idx));
%     outflow(:,:,3) = interface_flow;
% end
% % Down
% idx = Hf(:,:,4) > artificial_depth & [Hf(2:(end),:,4);zeros(1,size(depth_cell,2))] == artificial_depth;
% interface_flow = outflow(:,:,4);
% depth = Hf(:,:,4);
% if any(any(idx))
%     interface_flow(idx) = depth(idx).*sqrt(g*depth(idx));
%     outflow(:,:,4) = interface_flow;
% end
% % NE
% idx = Hf(:,:,6) > artificial_depth & [zeros(1,size(depth_cell,2));Hf(1:(end-1),(2:end),6)] == artificial_depth;
% interface_flow = outflow(:,:,4);
% depth = Hf(:,:,6);
% if any(any(idx))
%     interface_flow(idx) = depth(idx).*sqrt(g*depth(idx));
%     outflow(:,:,6) = interface_flow;
% end
% % SE
% idx = Hf(:,:,7) > artificial_depth & [Hf(2:(end),(2:end),7),zeros(1,size(depth_cell,2))] == artificial_depth;
% interface_flow = outflow(:,:,4);
% depth = Hf(:,:,7);
% if any(any(idx))
%     interface_flow(idx) = depth(idx).*sqrt(g*depth(idx));
%     outflow(:,:,7) = interface_flow;
% end
% % SW
% idx = Hf(:,:,8) > artificial_depth & [Hf(2:(end),(1:end-1),8),zeros(1,size(depth_cell,2))] == artificial_depth;
% interface_flow = outflow(:,:,8);
% depth = Hf(:,:,8);
% if any(any(idx))
%     interface_flow(idx) = depth(idx).*sqrt(g*depth(idx));
%     outflow(:,:,8) = interface_flow;
% end
% % NW
% idx = Hf(:,:,9) > artificial_depth & [zeros(1,size(depth_cell,2));Hf(1:(end-1),(1:end-1),0)] == artificial_depth;
% interface_flow = outflow(:,:,9);
% depth = Hf(:,:,9);
% if any(any(idx))
%     interface_flow(idx) = depth(idx).*sqrt(g*depth(idx));
%     outflow(:,:,9) = interface_flow;
% end
%% Limiting outflow to critical velocity
if flag_critical == 1
    critical_velocity = Hf.*sqrt(g*Hf);
    outflow = min(outflow,critical_velocity);
end

% Taking out negative values
outflow = Resolution*outflow/(Resolution^2)*1000*3600; % mm per hour
outflow(outflow<0) = 0;
outflow(isnan(outflow)) = 0; outflow(isinf(outflow)) = 0;
outflow(mask) = 0;
% Taking out values outside of the domain
outflow(idx_nan) = nan;

%% Intercell Volume
I_tot_end_cell = dt*sum(1/1000*1/3600*outflow,3)*Resolution^2; % Total outflow in m3
% matrix_store now becomes outflow
matrix_store = outflow; % mm per hour

%% Reservoir Boundary Condition - We are assuming that all flow drains
% towards the spillway
if flag_reservoir == 1   
	for ii = 1:length(reservoir_y)
        if ~isnan(yds1(ii))
            dtsup = d_tot(reservoir_y(ii),reservoir_x(ii))./1000; % % Water depth in the cell that has the boundary condition (m)
            dt_h = (time_step)/60; % timestep in hours
            % ---- First Boundary Condition ----- %
            available_volume = 1000*(max(dtsup - h1(ii),0))/dt_h; %  mm/h
            dh = min(k1(ii)*(max(dtsup - h1(ii),0))^k2(ii)/cell_area*1000*3600,available_volume)*dt_h; % mm
            I_tot_end_cell(reservoir_y(ii),reservoir_x(ii)) = I_tot_end_cell(reservoir_y(ii),reservoir_x(ii)) + dh/1000*cell_area;
            dtsup = dtsup - dh/1000;
            % Refreshing downstream cell
            d_tot(yds1(ii),xds1(ii)) = d_tot(yds1(ii),xds1(ii)) + dh; 
        else
            dh = 0;
        end
        if ~isnan(yds2(ii))
            % Refreshing downstream cell
            d_tot(yds1(ii),xds1(ii)) = d_tot(yds1(ii),xds1(ii)) + dh;
            % ---- Second Boundary Condition ----- %
            available_volume = 1000*(max(dtsup - h2(ii),0))/dt_h; %  mm/h
            dh = min(k3(ii)*(max(dtsup - h2(ii),0))^k4(ii)/cell_area*1000*3600,available_volume)*dt_h; % mm
            I_tot_end_cell(reservoir_y(ii),reservoir_x(ii)) = I_tot_end_cell(reservoir_y(ii),reservoir_x(ii)) + dh/1000*cell_area;
            % Refreshing downstream cell
            d_tot(yds2(ii),xds2(ii)) = d_tot(yds2(ii),xds2(ii)) + dh;
        end
	end
end

qout_left = matrix_store(:,:,1);
qout_right = matrix_store(:,:,2);
qout_up = matrix_store(:,:,3);
qout_down = matrix_store(:,:,4);
outlet_flow = matrix_store(:,:,5);
% Inclined Directions
qout_ne = matrix_store(:,:,6);
qout_se = matrix_store(:,:,7);
qout_sw = matrix_store(:,:,8);
qout_nw = matrix_store(:,:,9);

%% ---------------% Final depth at the cell % ---------------%
d_t = d_tot - I_tot_end_cell/cell_area*1000; % final depth in mm;
% Chaning depth_cell of B.C cells to NaN
for i = 1:length(reservoir_y)
    d_t(reservoir_y(i),reservoir_x(i)) = NaN;
end
end
