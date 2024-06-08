%@ -0,0 +1,242 @@
function [qout_left,qout_right,qout_up,qout_down,outlet_flow,qout_ne,qout_se,qout_sw,qout_nw,d_t,I_tot_end_cell] = CA_Routing_8D(reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2,flag_reservoir,elevation_cell,d_tot,roughness_cell,cell_area,time_step,h_0_cell,Resolution,I_tot_end_cell,outlet_index,outlet_type,slope_outlet,row_outlet,col_outlet,idx_nan,flag_critical,d_tolerance)
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
    ny_max = gpuArray(size(elevation_cell,1));
    nx_max = gpuArray(size(elevation_cell,2));
else
    d_t_min = d_tolerance/1000; % m
    ny_max = size(elevation_cell,1);
    nx_max = size(elevation_cell,2);  
end
% ---------------% Adding minimum slope to do calculations % ---------------%
h_min = d_t_min;% m
I_tot_end_previous = I_tot_end_cell; 

% --------------- Notation  % ---------------%
%   <-I-> (left right up down) = [left ; right ; up; down]
% --------------- Cell Depth and Water surface Elevation  % ---------------%
depth_cell = max(d_tot./1000,0); % meters (fixing zero values)

wse_cell = elevation_cell + depth_cell;

% --------------- Depth Differences (delta_h) for All Cells  % ---------------%
matrix_store(:,:,1) = [wse_cell(:,1), wse_cell(:,2:(nx_max)) - wse_cell(:,1:(nx_max-1))]; % left
matrix_store(:,:,2) = [wse_cell(:,1:(nx_max-1)) - wse_cell(:,2:(nx_max)), wse_cell(:,end)]; % right
matrix_store(:,:,3) = [wse_cell(1,:); wse_cell(2:(ny_max),:) - wse_cell(1:(ny_max-1),:)]; % up
matrix_store(:,:,4) = [wse_cell(1:(ny_max-1),:) - wse_cell(2:(ny_max),:); wse_cell(end,:)]; % down

% --------------- Depth Differences for All NE, SE, SW, NW  % ---------------%
% 6, 7, 8, 9, respectively
matrix_store(:,:,6) =  [wse_cell(1,:); [wse_cell(2:(ny_max),1:(nx_max-1)) - wse_cell(1:(ny_max-1),2:(nx_max)), wse_cell(1:(ny_max-1),end)]]; % NE
matrix_store(:,:,7) =  [[wse_cell(1:(ny_max-1),1:(nx_max-1)) - wse_cell(2:ny_max,2:(nx_max)), wse_cell(2:(ny_max),end)]; wse_cell(end,:)]; % SE
matrix_store(:,:,8) = [[wse_cell(2:(ny_max),1), wse_cell(1:(ny_max-1),2:(nx_max)) - wse_cell(2:(ny_max),1:(nx_max-1))]; wse_cell(end,:)]; % SW
matrix_store(:,:,9) = [wse_cell(1,:);[wse_cell(1:(ny_max-1),1), wse_cell(2:(ny_max),2:(nx_max)) - wse_cell(1:(ny_max-1),1:(nx_max-1))]]; % NW

%%%%% CHECKED 9/28/2021 %%%%%
% --------------- Outlet Calculations  % ---------------%
if outlet_type == 1
    S_0 = slope_outlet; % Normal slope
    % V = h * area
    matrix_store(:,:,5) = (S_0*Resolution*outlet_index); % Normal slope depth difference
    %%%%% CHECKED 9/28/2021 %%%%%
else
    S_0 = zeros(size(elevation_cell));
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

matrix_store(logical(matrix_store < h_min | isnan(matrix_store))) = 0;
% ---------------% Available Volumes  % ---------------%
%slope
matrix_store(:,:,1:4) = matrix_store(:,:,1:4)/Resolution; % Calculating the slopes
matrix_store(:,:,5) = matrix_store(:,:,5)/Resolution;
matrix_store(:,:,6:end) = matrix_store(:,:,6:end)/(sqrt(2)*Resolution);
[sf,index] = max(matrix_store,[],3); % Maximum Slope

idx = index <= 5;
Distance = zeros(size(elevation_cell));
Distance(idx) = Resolution;
Distance(idx == 0) = sqrt(2)*Resolution;

matrix_store(:,:,1:4) = matrix_store(:,:,1:4)*Resolution; %Getting back the delta_h matrix
matrix_store(:,:,5) = matrix_store(:,:,5)*Resolution;
matrix_store(:,:,6:end) = matrix_store(:,:,6:end)*(sqrt(2)*Resolution);

% Vol 
matrix_store = cell_area*matrix_store; % Recycling variables; delta_h turns into Volume.

%%%%% CHECKED 9/28/2021 %%%%%
% ---------------% Minimum and Max Volumes  % ---------------%
% Minimum
c = Resolution^2*1000; % large number (100 m of water depth)
Vol_min_matrix = min(matrix_store+(matrix_store==0)*c,[],3);

% Total Volume
Vol_tot = sum(matrix_store,3);

% ---------------% Weights  % ---------------%
matrix_store = matrix_store./(Vol_tot + Vol_min_matrix);
matrix_store(isnan(matrix_store))=0;
weigth_max = max(matrix_store,[],3);

%% ---------------% Velocity Calculation %---------------%
% Velocity to the cell with the highest gradient
if flag_critical == 1
    v_m = min(sqrt(9.81.*depth_cell),1./roughness_cell.*((max(depth_cell - h_0_cell/1000,0))).^(2/3).*sqrt(sf));
else
    v_m = 1./roughness_cell.*((max(depth_cell - h_0_cell/1000,0))).^(2/3).*sqrt(sf);% Manning's equation only
end
%% ---------------% Intercell Volume Calculation %---------------%
% Estimated Volume that goes to the cell with the highest gradient
%         I_m = v_m.*depth_cell*Resolution*time_step*60; % Velocity x Area x Time_step (m3)

% ---------------% Total outflow volume % ---------------%
%%%% In case we won't want to input I_tot_end_previous %%%
if sum(sum(I_tot_end_previous)) == 0
    I_tot_end_previous = 0.9*cell_area*depth_cell;
end
% I_tot_end = min(depth_cell*cell_area,I_m/max(weight),Vol_min_matrix + I_tot_begin);
I_tot_end_cell = min(depth_cell.*cell_area,v_m.*depth_cell./weigth_max.*(Distance*time_step*60)); % m3
I_tot_end_cell = min(I_tot_end_cell,I_tot_end_previous + Vol_min_matrix);
I_tot_end_cell(logical(d_tot < d_t_min*1000) + idx_nan > 0) = 0;  % If the depth is smaller than dmin ATTENTION HERE

% ---------------% Outflow Volume to each direction % ---------------%
matrix_store = matrix_store.*I_tot_end_cell; % m3
I_tot_end_cell = sum(matrix_store,3); % Volume that leaves the cell in m3 
matrix_store = matrix_store./(time_step*60)/(cell_area)*1000*3600; % mm/hr

% Reservoir Boundary Condition - We are assuming that all flow drains
% towards the spillway
if flag_reservoir == 1
	for ii = 1:length(reservoir_y)
		dtsup = d_tot(reservoir_y(ii),reservoir_x(ii))./1000; % % Water depth in the cell that has the boundary condition (m)
        
        dt_h = (time_step)/60; % timestep in hours
        % ---- First Boundary Condition ----- %
        available_volume = 1000*(max(dtsup - h1(ii),0))/dt_h; %  mm/h
        dh = min(k1(ii)*(max(dtsup - h1(ii),0))^k2(ii)/cell_area*1000*3600,available_volume)*dt_h; % mm
        I_tot_end_cell(reservoir_y(ii),reservoir_x(ii)) = I_tot_end_cell(reservoir_y(ii),reservoir_x(ii)) + dh/1000*cell_area;
        dtsup = dtsup - dh/1000;
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


end