function [qout_left,qout_right,qout_up,qout_down,outlet_flow,d_t,I_tot_end_cell,outflow,Hf,Qc,Qf,Qci,Qfi] = Local_Inertial_Model_D4(flag_numerical_scheme,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2,flag_reservoir,z,d_tot,d_p,roughness_cell,cell_area,time_step,Resolution,outlet_index,outlet_type,slope_outlet,row_outlet,col_outlet,d_tolerance,outflow,idx_nan,flag_critical,flag_subgrid,nc,nf,River_Width, River_Depth, idx_rivers,Qc_prev,Qf_prev,Qci_prev,Qfi_prev)

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                 %
%                 Produced by Marcus Nobrega Gomes Junior         %
%                 e-mail:marcusnobrega.engcivil@gmail.com         %
%                           September 2021                        %
%                                                                 %
%                 Last Updated: 5 August, 2024                    %
%                                                                 %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§



%%% Domain Dimensions
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
h_min = 0;  % In cases where inflow is being modeling, this value has to be 0

% --------------- Notation  % ---------------%
%   <-matrix3D-> (left right up down) = [left ; right ; up; down] going
%   outside of the cell

% --------------- Cell Depth and Water surface Elevation  % ---------------%
% depth_cell = 0.5*(d_tot + d_p); % This is important when inflow hydrograph is simulated
% depth_cell = d_tot;
depth_cell = d_p;
depth_cell = max(depth_cell/1000,0); % meters (fixing zero values)

% Water Surface Elevation [m]
y = z + depth_cell; % Total depth minus outlet flow

%% --------------- Slope for all cell boundaries % ---------------%
nan_col = nan*y(:,1); nan_row = nan*y(1,:);

matrix_store = 0*outflow;
% South to North and West to East Positive
% x-x
% matrix_store(:,:,1) = [nan_col, y(:,2:end) - y(:,1:(end-1))]/Resolution; % Left
% matrix_store(:,:,2) = [y(:,2:end) - y(:,1:(end-1)), nan_col]/Resolution; % Right
matrix_store(:,:,1) = [y(:,2:end) - y(:,1:(end-1)), nan_col]/Resolution; % Right

% y-y
% matrix_store(:,:,3) = [nan_row; y(1:(end-1),:) - y(2:end,:)]/Resolution; % Up
matrix_store(:,:,2) = [nan_row; y(1:(end-1),:) - y(2:end,:)]/Resolution; % Up
% matrix_store(:,:,4) = [y(1:(end-1),:) - y(2:end,:); nan_row]/Resolution; % Down

% if nansum(nansum(matrix_store(:,2,1) - matrix_store(:,1,2)))
%     ttt = 1;
% end

%% ---------------- Hf (Effective Water Depth for Flow) ----------- %
Hf = 0*outflow;
% x-x 
% Hf(:,:,1) = [nan_col, max(y(:,2:end),y(:,1:(end-1))) - max(z(:,2:end),z(:,1:(end-1)))]; % left
% Hf(:,:,2) = [max(y(:,2:nx), y(:,1:(nx-1))) - max(z(:,2:nx), z(:,1:(nx-1))), nan_col]; % right
Hf(:,:,1) = [max(y(:,2:nx), y(:,1:(nx-1))) - max(z(:,2:nx), z(:,1:(nx-1))), nan_col]; % right
% y-y
% Hf(:,:,3) = [nan_row; max(y(1:(end-1),:),y(2:end,:)) - max(z(1:(end-1),:),z(2:end,:))]; % up
Hf(:,:,2) = [nan_row; max(y(1:(end-1),:),y(2:end,:)) - max(z(1:(end-1),:),z(2:end,:))]; % up
% Hf(:,:,4) = [max(y(1:(end-1),:), y(2:end,:)) - max(z(1:(end-1),:), z(2:end,:)); nan_row]; % down

Hf(isnan(matrix_store)) = 0;

% if nansum(nansum(Hf(:,2,1) - Hf(:,1,2)))
%     ttt = 1;
% end

% Low depth cells are considered an artificial depth
mask = logical((Hf < h_min));
% Artificial Depth
artificial_depth = 0;
Hf(mask) = artificial_depth; % No outflow from cells with very low depth
% check_cells
%% Local-Inertial Formulation
outflow = outflow/1000/3600*Resolution^2/Resolution; % m2 per sec. It comes from mm/h
dt = time_step*60; % time-step in seconds
g = 9.81;
outflow_prev = outflow;

% Sub-grid channels
% Treatment using sub-grid approach for channels
if flag_subgrid == 1
    % All cells have sub-grid channels and floodplains following the same
    % roughness coefficients of the 2D domain
    % Only the original formulation is considered in this case
    [Q,Qc,Qf,Qci,Qfi,C_a] = subgrid_channel(depth_cell, River_Width, z + River_Depth, z , Resolution, nc, nf, Qc_prev, Qf_prev,Qci_prev,Qfi_prev,g,dt,idx_rivers);
    outflow = Q./(C_a/Resolution); % m2/s normalized by flow surface area
    % Treating Domain Issues
    outflow(isnan(outflow)) = 0; outflow(isinf(outflow)) = 0;    
else
    % We are modeling using the original coarse grid approach with
    % different numerical schemes
    Qc = 0;
    Qf = 0;
    Qci = 0;
    Qfi = 0;
    Q = 0;
    % -------------- Inertial Solver ------------- %
    [outflow] = Inertial_Solver(flag_numerical_scheme,outflow_prev,dt,Hf,matrix_store,roughness_cell,Resolution,idx_nan);
    % Treating Domain Issues
    outflow(isnan(outflow)) = 0; outflow(isinf(outflow)) = 0;  
    C_a = ones(size(roughness_cell))*Resolution^2; % cell area in m2;
end


%% Intercell Volume
% 1 - Right (Negative)
% 2 - Up (Negatove)
% 3 - Outlet (Negative)

% Vol_Flux = dt*Resolution*(outflow(:,:,1) - outflow(:,:,2) ...
%                          - outflow(:,:,3) + outflow(:,:,4));

Vol_Flux = dt*(C_a/Resolution).* ...
               ([zeros(ny,1) , outflow(:,1:(nx-1),1)] - outflow(:,:,1) ...
                - outflow(:,:,2) + [outflow(2:end,:,2); zeros(1,nx)]); % m3

if flag_subgrid == 1
    % % Adding NE and SE fluxes
    % outfluxes = Qci + Qfi;
    % Vol_Flux = Vol_flux + dt*(C_a/Resolution).* ...
    %                            ([zeros(ny,1) , outfluxes(:,1:(nx-1),1)] - outfluxes(:,:,1) ...
    %                             - outfluxes(:,:,2) + [outfluxes(2:end,:,2); zeros(1,nx)]); % m3
    % 
    % hf(2:end,1:(end-1),3) = max(y(1:(end-1),2:end), y(2:end,1:(end-1))) - max(z(1:(end-1),2:end), z(2:end,1:(end-1)));
    % hf(1:(end-1),1:(end-1),4) = max(y(2:end,2:(end)), y(1:(end-1),1:(end-1))) - max(z(2:end,2:(end)), z(1:(end-1),1:(end-1)));
end
%% Eliminating Surplus Velocities at Wet-Dry Interfaces
% Left 
% idx = Hf(:,:,1) > artificial_depth & [Hf(:,2:(end-1),1),zeros(size(depth_cell,1),1)] == artificial_depth;
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
%% Limiting outflow to critical velocity
if flag_critical == 1
    critical_velocity = Hf.*sqrt(g*Hf);
    outflow = min(outflow,critical_velocity);
end
outflow = Resolution*outflow/(Resolution^2)*1000*3600; % mm per hour
matrix_store = outflow; % mm per hour
%% Reservoir Boundary Condition - We are assuming that all flow drains
if flag_reservoir == 1   
	for ii = 1:length(reservoir_y)
        if ~isnan(yds1(ii))
            dtsup = d_tot(reservoir_y(ii),reservoir_x(ii))./1000; % % Water depth in the cell that has the boundary condition (m)
            dt_h = (time_step)/60; % timestep in hours
            % ---- First Boundary Condition ----- %
            available_volume = 1000*(max(dtsup - h1(ii),0))/dt_h; %  mm/h
            dh = min(k1(ii)*(max(dtsup - h1(ii),0))^k2(ii)/cell_area*1000*3600,available_volume)*dt_h; % mm
            Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh/1000*cell_area;
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
            Vol_Flux(reservoir_y(ii),reservoir_x(ii)) = Vol_Flux(reservoir_y(ii),reservoir_x(ii)) + dh/1000*cell_area;
            % Refreshing downstream cell
            d_tot(yds2(ii),xds2(ii)) = d_tot(yds2(ii),xds2(ii)) + dh;
        end
	end
end

%% Output Fluxes
qout_left = -[zeros(ny,1), matrix_store(:,:,1)]; % See the signals following the convention
qout_right = [matrix_store(:,:,1)];
qout_up = matrix_store(:,:,2);
qout_down = -[matrix_store(2:end,:,2); zeros(1,nx)];


%% ---------------% Final depth at the cell % ---------------%

% Cell mass balance
% d_tot = d_p + rainfall, inflow or anything else
d_t = d_tot + Vol_Flux./C_a*1000 ; % final depth in mm

if flag_subgrid == 1 % Maybe we have a change from inbank <-> overbank
    % Eq. 15 and 16 of A subgrid channel model for simulating river hydraulics andfloodplain inundation over large and data sparse areas
    % Inbank - Overbank
    idx = d_t/1000 > River_Depth & d_tot/1000 <= River_Depth & idx_rivers; % Cells in which there is a change from inbank to overbank
    if sum(sum(idx)) > 0
        ttt = 1;
    end
    d_t(idx) = 1000*(d_t(idx)/1000 - (d_t(idx)/1000 - z(idx) + (z(idx) - River_Depth(idx))).*(1 - River_Width(idx)/Resolution)); % mm 
    % Overbank - Inbank
    idx = d_t/1000 <= River_Depth & d_tot/1000 >= River_Depth & idx_rivers; % Cells in which there is a change from overbank to inbank
    if sum(sum(idx)) > 0
        ttt = 1;
    end 
    factor = d_t(idx)/1000 - z(idx) + (z(idx) - River_Depth(idx));
    d_t(idx) = 1000*(d_t(idx)/1000 + (Resolution*(factor))./River_Width(idx) ...
                        - (d_t(idx)/1000 - z(idx) + (z(idx) - River_Depth(idx)))); % mm
end

lost_mass = 1/1000*abs(sum(sum(d_t(d_t<0))))*cell_area;
if lost_mass > 1*cell_area
    warning('Lost mass too high.')
end
d_t(d_t<0) = 0;



% Outlet Flow
% --------------- Outlet Calculations  % ---------------%
if outlet_type == 1
    S_0 = slope_outlet*outlet_index; % Normal slope
else
    S_0 = zeros(size(z));
    for i = 1:length(row_outlet)
        % Checking left, right up, and down
        row = row_outlet(i); % Row of outlet
        col = col_outlet(i); % Col of outlet
        S_0(row,col) = (max((d_t(row,col)/1000),0).^(-1/6)).*sqrt(9.81).*roughness_cell(row,col); % Critical Slope
    end
    S_0(isinf(S_0)) = 0;
end
mask_outlet = S_0 ~= 0;

% Outlet Mass Balance
% Outlet Flow Depth
Hf_outlet = zeros(ny,nx);
Hf_outlet(mask_outlet) = max(d_t(mask_outlet),0)/1000; % meters
Hf(:,:,3) = Hf_outlet;

% Outlet Mass Balance
outlet_flow = 1./roughness_cell.*Resolution.*Hf(:,:,3).^(5/3).*abs(S_0).^(0.5)./(Resolution^2)*1000*3600; % mm/h
outlet_flow = min(outlet_flow,max((d_t-1000*h_min),0)/(time_step/60));
outlet_flow(~mask_outlet) = 0;

% Mass Balance Check
% Current Stored Volume
current_volume = nansum(nansum(d_tot/1000))*cell_area;
% Final Stored Volume
final_volume = nansum(nansum(d_t/1000))*cell_area;

error_vol = current_volume - final_volume;
if error_vol > 10
    ttt = 1;
end

% Final Depth
d_t = d_t  - outlet_flow*(time_step/60); % final depth in mm



%% ---------------% Total Flow that Leaves the Cell ----------- %
mask = outflow; mask(mask<0) = 0; % We want to get only outflow flux per cell
I_tot_end_cell = sum(mask,3)*dt/1000*1/3600*Resolution^2; % m3

end