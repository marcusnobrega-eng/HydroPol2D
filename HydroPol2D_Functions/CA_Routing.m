function [qout_left,qout_right,qout_up,qout_down,outlet_flow,d_t,I_tot_end_cell,I] = CA_Routing(reservoir_dir,reservoir_x,reservoir_y,Kv,p,flag_reservoir,elevation_cell,d_tot,roughness_cell,cell_area,time_step,h_0_cell,Resolution,I_tot_end_cell,outlet_check,outlet_type,slope_outlet,row_outlet,col_outlet,idx_nan,flag_critical)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                 %
%                 Produced by Marcus Nobrega Gomes Junior         %
%                 e-mail:marcusnobrega.engcivil@gmail.com         %
%                           September 2021                        %
%                                                                 %
%                 Last Updated: 7 July, 2022                      %
%                                                                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
%   -----  DESCRIPTION  -----
%   This function estimates the transfered volume from a cell to the
%   neighbour cells using a weigthing approach in terms of the available
%   volume in the neighbour cells
%§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
% ---------------% Finding Nans or Values Below a Threshold  % ---------------%
if isgpuarray(cell_area)
    d_t_min = gpuArray(1e-6); % m
    ny_max = gpuArray(size(elevation_cell,1));
    nx_max = gpuArray(size(elevation_cell,2));
else
    d_t_min = 1e-6; % m
    ny_max = size(elevation_cell,1);
    nx_max = size(elevation_cell,2);    
end
% idx2 = logical(d_tot < d_t_min*1000) + idx_nan > 0;

% ---------------% Adding minimum slope to do calculations % ---------------%
h_min = d_t_min;% m
I_tot_end_previous = I_tot_end_cell; 

% --------------- Notation  % ---------------%
%   <-I-> (left right up down) = [left ; right ; up; down]

% --------------- Cell Depth and Water surface Elevation  % ---------------%
depth_cell = d_tot./1000; % meters
wse_cell = elevation_cell + depth_cell;

% --------------- Depth Differences for All Cells  % ---------------%
% 
% Elevation_Properties.elevation_cell = elevation;
% Elevation_Properties.elevation_left_t = [zeros(ny_max,1),elevation(:,1:(nx_max-1))];
% Elevation_Properties.elevation_right_t = [elevation(:,(2:(nx_max))) zeros(ny_max,1)];
% Elevation_Properties.elevation_up_t = [zeros(1,nx_max) ; elevation(1:(ny_max-1),:)];
% Elevation_Properties.elevation_down_t = [elevation(2:(ny_max),:) ; zeros(1,nx_max)];


delta_h(:,:,1) = (wse_cell - [zeros(ny_max,1),elevation_cell(:,1:(nx_max-1))] - [zeros(ny_max,1),depth_cell(:,1:(nx_max-1))]); % left
delta_h(:,:,2) = (wse_cell - [elevation_cell(:,(2:(nx_max))) zeros(ny_max,1)] - [depth_cell(:,(2:(nx_max))) zeros(ny_max,1)]); % right
delta_h(:,:,3) = (wse_cell - [zeros(1,nx_max) ; elevation_cell(1:(ny_max-1),:)] - [zeros(1,nx_max) ; depth_cell(1:(ny_max-1),:)]); % up
delta_h(:,:,4) = (wse_cell - [elevation_cell(2:(ny_max),:) ; zeros(1,nx_max)] - [depth_cell(2:(ny_max),:) ; zeros(1,nx_max)]); % down

%%%%% CHECKED 9/28/2021 %%%%%
% --------------- Outlet Calculations  % ---------------%
if outlet_type == 1
    S_0 = slope_outlet; % Normal slope
    % V = h * area
    delta_h(:,:,5) = (S_0*Resolution*outlet_check); % Normal slope depth difference
    %%%%% CHECKED 9/28/2021 %%%%%
else
    S_0 = zeros(size(elevation_cell));
    for i = 1:length(row_outlet)
        % Checking left, right up, and down
        row = row_outlet(i); % Row of outlet
        col = col_outlet(i); % Col of outlet
        S_0(row,col) = (depth_cell(row,col).^(-1/6)).*sqrt(9.81).*roughness_cell(row,col); % Critical Slope
    end
    delta_h(:,:,5) = (S_0*Resolution);
end

% ---------------% Correcting delta_h  % ---------------%
% delta_h(isnan(delta_h)) = 0; % TESTING
% delta_h(isinf(delta_h)) = 0; % TESTING

% ---------------% Boundary Conditions  % ---------------%
% These only applied when we are modeling a rectangular grid watershed
delta_h(:,1,1) = 0; % Left
%     direction = ones(size(delta_h(:,:,1)));
delta_h(:,end,2) = 0; % Right
delta_h(1,:,3) = 0; % Up
delta_h(end,:,4) = 0; % Down

if flag_reservoir == 1
    for ii = 1:length(reservoir_y)
        delta_h(reservoir_y(ii),reservoir_x(ii),:) = 0; % Imposing no flow to all directions       
    end
end

% idx = logical(delta_h < h_min | isnan(delta_h)); % A lot of attention here
delta_h(logical(delta_h < h_min | isnan(delta_h))) = 0;

% ---------------% Available Volumes  % ---------------%
Vol = cell_area*delta_h; % 3-D Array with pages following direction notation
% Max h
maxh = max(delta_h,[],3);
%%%%% CHECKED 9/28/2021 %%%%%
% ---------------% Minimum and Max Volumes  % ---------------%
% Minimum
c = Resolution^2*1000; % large number (100 m of water depth)
Vol_nan = Vol;
Vol_nan(logical(delta_h < h_min | isnan(delta_h))) = c;
Vol_min_matrix = min(Vol_nan,[],3);
% Total Volume
Vol_tot = sum(Vol,3);

% ---------------% Weights  % ---------------%
weight = Vol./(Vol_tot + Vol_min_matrix);
weight(isnan(weight)) = 0; % Attention here
weight_max = max(weight,[],3);

%% ---------------% Velocity Calculation %---------------%
% Velocity to the cell with the highest gradient
if flag_critical == 1
    v_m = min(sqrt(9.81.*depth_cell),1./roughness_cell.*((max(depth_cell - h_0_cell/1000,0))).^(2/3).*sqrt(maxh/Resolution));
else
    v_m = 1./roughness_cell.*((max(depth_cell - h_0_cell/1000,0))).^(2/3).*sqrt(maxh/Resolution);% Manning's equation only
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
I_tot_end_cell = min(depth_cell.*cell_area,v_m.*depth_cell./weight_max*(Resolution*time_step*60)); % m3
I_tot_end_cell = min(I_tot_end_cell,I_tot_end_previous + Vol_min_matrix);
% idx2 = logical(d_tot < d_t_min*1000) + idx_nan > 0
% I_tot_end_cell(idx2) = 0;
I_tot_end_cell(logical(d_tot < d_t_min*1000) + idx_nan > 0) = 0;  % If the depth is smaller than dmin ATTENTION HERE

% ---------------% Outflow Volume to each direction % ---------------%
I = weight.*I_tot_end_cell; % m3
I_tot_end_cell = sum(I,3);% Volume that leaves the cell in m3
qout = I./(time_step*60)/(cell_area)*1000*3600; % mm/hr

% Reservoir Boundary Condition - We are assuming that all flow drains
% towards the spillway
if flag_reservoir == 1
    for ii = 1:length(reservoir_y)
        % 1 left, 2 right , 3 up, 4 down
        if reservoir_dir(ii) == 1
            head = wse_cell(reservoir_y(ii),reservoir_x(ii) + 1); % Head of water in the neighbour cell
        elseif reservoir_dir(ii) == 2
            head = wse_cell(reservoir_y(ii),reservoir_x(ii) - 1);
        elseif reservoir_dir(ii) == 3
            head = wse_cell(reservoir_y(ii) + 1,reservoir_x(ii));
        elseif reservoir_dir(ii) == 4
            head = wse_cell(reservoir_y(ii) - 1,reservoir_x(ii) - 1);
        end
        qout(reservoir_y(ii),reservoir_x(ii),:) = 0;
        qout(reservoir_y(ii),reservoir_x(ii),reservoir_dir(ii)) = Kv(ii)*(max(head - p(ii),0))^(3/2)/cell_area*1000*3600; % m3/s to mm/h
        % Boundary Condition of Maximum Flow
        Available_Volume = 1/(time_step*60)*1000*depth_cell(reservoir_y(ii),reservoir_x(ii)); % mm/h
        qout(reservoir_y(ii),reservoir_x(ii),reservoir_dir(ii)) = min(qout(reservoir_y(ii),reservoir_x(ii),reservoir_dir(ii)),Available_Volume); % mm/h
        % Total Outflow from this cell
        I_tot_end_cell(reservoir_y(ii),reservoir_x(ii)) = (time_step*60)*1/(3600*1000)*cell_area*qout(reservoir_y(ii),reservoir_x(ii),reservoir_dir(ii)); % m3
    end    
end
    qout_left = qout(:,:,1);
    qout_right = qout(:,:,2);
    qout_up = qout(:,:,3);
    qout_down = qout(:,:,4);
    outlet_flow = qout(:,:,5);

%% ---------------% Final depth at the cell % ---------------%
d_t = d_tot - I_tot_end_cell/cell_area*1000; % final depth in mm;
end