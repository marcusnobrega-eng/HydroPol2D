% Diffive Wave Model

function [qout_left,qout_right,qout_up,qout_down,outlet_flow,d_t,I_tot_end_cell,outflow,Hf,Qc,Qf,Qci,Qfi,C_a] = Diffusive_Wave_Model(flag_numerical_scheme,reservoir_x,reservoir_y,k1,h1,k2,k3,h2,k4,yds1,xds1,yds2,xds2,flag_reservoir,z,d_tot,d_p,roughness_cell,cell_area,time_step,Resolution,outlet_index,outlet_type,slope_outlet,row_outlet,col_outlet,d_tolerance,outflow,idx_nan,flag_critical,flag_subgrid,nc,nf,River_Width, River_Depth,Qc_prev,Qf_prev,Qci_prev,Qfi_prev,C_a_prev)

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
% Rivers
% River_Width = River_Width*0;
% River_Depth = River_Depth*0;
% River_Width = ones(ny,nx)*100;
% River_Depth = ones(ny,nx)*10;
% River_Width(1,6) = 50;
% River_Depth(1,6) = 1;
% z(:,6) = 10;
% z(1,6) = 10;
% River_Depth(River_Depth > 0) = 1;
% River_Width(River_Width > 0) = 50;
% Adding a tributary
% River_Width(5,:) = 50;
% River_Depth(5,:) = 1;
idx_rivers = River_Width > 0; % Rivers are now cells with no zero width 
% ---------------% Adding minimum slope to do calculations % ---------------%
h_min = 5/1000;  % In cases where inflow is being modeling, this value has to be 0
dt = time_step*60;

% --------------- Notation  % ---------------%
%   <-matrix3D-> (left right up down) = [left ; right ; up; down] going
%   outside of the cell

% --------------- Cell Depth and Water surface Elevation  % ---------------%
% depth_cell = 0.5*(d_tot + d_p); % This is important when inflow hydrograph is simulated
depth_cell = d_tot;
% depth_cell = d_p;
depth_cell = max(depth_cell/1000,0); % meters (fixing zero values)

% Current Stored Volume
if flag_subgrid == 1
    current_volume = nansum(nansum((Resolution - River_Width).*Resolution.*max((d_tot/1000 - River_Depth),0))) + ...
                      nansum(nansum(Resolution.*River_Width.*d_tot/1000)); % m3
else
    current_volume = nansum(nansum(C_a_prev.*d_tot/1000)); % m3
end


% Water Surface Elevation [m]
y = z + depth_cell; % Total depth minus outlet flow

%% --------------- Slope for all cell boundaries % ---------------%
nan_col = nan*y(:,1); nan_row = nan*y(1,:);

matrix_store = 0*outflow;
% South to North and West to East Positive
% x-x
matrix_store(:,:,1) = [y(:,2:end) - y(:,1:(end-1)), nan_col]/Resolution; % Right
% y-y
matrix_store(:,:,2) = [nan_row; y(1:(end-1),:) - y(2:end,:)]/Resolution; % Up

% Limiter in water surface slope to avoid numerical instabilities
% dwse_dx_limiter = 0.10; % Limiter
% 
% matrix_store(matrix_store > dwse_dx_limiter) = dwse_dx_limiter;
% matrix_store(matrix_store < -dwse_dx_limiter) = -dwse_dx_limiter;

%% ---------------- Hf (Effective Water Depth for Flow) ----------- %
Hf = 0*outflow;
% x-x 
Hf(:,:,1) = [max(y(:,2:nx), y(:,1:(nx-1))) - max(z(:,2:nx), z(:,1:(nx-1))), nan_col]; % right
% y-y
Hf(:,:,2) = [nan_row; max(y(1:(end-1),:),y(2:end,:)) - max(z(1:(end-1),:),z(2:end,:))]; % up
Hf(isnan(matrix_store)) = 0;

% Low depth cells are considered an artificial depth
mask = logical((Hf < h_min));
mask_depth = depth_cell < h_min;

% New mask
mask = logical(mask + repmat(mask_depth,1,1,3)); % Fail in depth and fail in Hf

% Artificial Depth
% artificial_depth = 0;
% Hf(mask) = artificial_depth; % No outflow from cells with very low depth
% 
% if max(max(max(Hf))) > 0
%     ttt = 1;
% end
%% Diffusive Wave Formulation
   % We are modeling using the original coarse grid approach with
    % different numerical schemes
    Qc = 0;
    Qf = 0;
    Qci = 0;
    Qfi = 0;
    Q = 0;

    % -------------- Diffusive Wave ------------- %
    % x-x
    outflow(:,:,1) = Manning_Equation(roughness_cell, Hf(:,:,1), matrix_store(:,:,1),1); % m2/s or m3/s per unit width
    % y-y
    outflow(:,:,2) = Manning_Equation(roughness_cell, Hf(:,:,2), matrix_store(:,:,2),1); % m2/s or m3/s per unit width
    % Treating Domain Issues
    outflow(isnan(outflow)) = 0; outflow(isinf(outflow)) = 0;  
    C_a = ones(size(roughness_cell))*Resolution^2; % cell area in m2;
    cell_width = C_a/Resolution; % cell width in m


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
    critical_velocity = Hf(:,:,1:size(outflow,3)).*sqrt(9.81*Hf(:,:,1:size(outflow,3))); % m3/s per unit width
    outflow = min(outflow,critical_velocity); % m2/s (normalized by flow width)
end

%% Flow limiter
alpha = 0.05;
% x-x
Hf_ = Hf(:,:,1); 
flow = outflow(:,:,1); % m2/s or m3/s per unit width
flow(mask(:,:,1)) = min((alpha*Resolution/dt).*Hf_(mask(:,:,1)), outflow(mask(:,:,1))); % Flow Limiter
outflow(:,:,1) = flow;
% y-y
Hf_ = Hf(:,:,2); 
flow = outflow(:,:,2); % m2/s or m3/s per unit width
flow(mask(:,:,2)) = min((alpha*Resolution/dt).*Hf_(mask(:,:,2)), outflow(mask(:,:,2))); % Flow Limiter
outflow(:,:,2) = flow;
%% Intercell Volume

outflow_rate = outflow.*(cell_width); % m3/s
outflow = outflow_rate./(C_a)*1000*3600; % mm per hour
matrix_store = outflow; % mm per hour

% 1 - Right (Negative)
% 2 - Up (Negative)
% 3 - Outlet (Negative)

% Vol_Flux = dt*Resolution*(outflow(:,:,1) - outflow(:,:,2) ...
%                          - outflow(:,:,3) + outflow(:,:,4));

% Vol_Flux = dt*(C_a/Resolution).* ...
%                ([zeros(ny,1) , outflow(:,1:(nx-1),1)] - outflow(:,:,1) ...
%                 - outflow(:,:,2) + [outflow(2:end,:,2); zeros(1,nx)]); % m3

Vol_Flux = dt*([zeros(ny,1) , outflow_rate(:,1:(nx-1),1)] - outflow_rate(:,:,1) ...
                - outflow_rate(:,:,2) + [outflow_rate(2:end,:,2); zeros(1,nx)]); % m3


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
qout_left = -[zeros(ny,1), matrix_store(:,(1:end-1),1)]; % See the signals following the convention
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
    idx = d_t/1000 > River_Depth & d_p/1000 <= River_Depth & idx_rivers; % Cells in which there is a change from inbank to overbank
    if sum(sum(idx)) > 0
        d_t(idx) = 1000*(d_t(idx)/1000 - (d_t(idx)/1000 - z(idx) + (z(idx) - River_Depth(idx))).*(1 - River_Width(idx)/Resolution)); % mm 
        C_a(idx) = Resolution^2;
    end
    % Overbank - Inbank
    idx = d_t/1000 <= River_Depth & d_p/1000 >= River_Depth & idx_rivers; % Cells in which there is a change from overbank to inbank
    if sum(sum(idx)) > 0
        factor = d_t(idx)/1000 - z(idx) + (z(idx) - River_Depth(idx));
        d_t(idx) = 1000*(d_t(idx)/1000 + (Resolution*(factor))./River_Width(idx) ...
            - (d_t(idx)/1000 - z(idx) + (z(idx) - River_Depth(idx)))); % mm
        C_a(idx) = River_Width(idx)*Resolution;
    end 
end
lost_mass = 1/1000*abs(sum(sum(d_t(d_t<0).*C_a(d_t<0))));

% if lost_mass > 0.1*Resolution^2
%     error('Lost mass too high.')
% end

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
outlet_flow = 1./roughness_cell.*cell_width.*Hf(:,:,3).^(5/3).*abs(S_0).^(0.5)./(C_a)*1000*3600; % mm/h
outlet_flow = min(outlet_flow,max((d_t-1000*h_min),0)/(time_step/60));
outlet_flow(~mask_outlet) = 0;

% Final Depth
d_t = d_t  - outlet_flow*(time_step/60); % final depth in mm


% Mass Balance Check
% Final Stored Volume
if flag_subgrid == 1
    final_volume = nansum(nansum((Resolution - River_Width).*Resolution.*max((d_t/1000 - River_Depth),0))) + ...
                      nansum(nansum(Resolution.*River_Width.*d_t/1000)); % m3
else
    final_volume = nansum(nansum(C_a.*d_t/1000)); % m3
end

error_vol = (current_volume - final_volume) - nansum(nansum(C_a.*1/1000*time_step/60.*outlet_flow));
if abs(error_vol) > 10
    ttt = 1;
end




%% ---------------% Total Flow that Leaves the Cell ----------- %
mask = outflow; 
% mask(mask<0) = 0; % We want to get only outflow flux per cell (maybe we dont need it)
% I_tot_end_cell = sum(mask,3)*dt/1000*1/3600*Resolution^2; % m3
I_tot_end_cell = abs(sum(mask,3))*dt/1000*1/3600*Resolution^2; % m3

% surf(d_t); view(0,90);  colorbar; colormap('jet')
% pause(0.00001)

    function [Q] = Manning_Equation(n,h,Sf,dx)
        % Input:
        % - n> Manning's roughness coefficient [sm-1/3]
        % - h> Water depth [m]
        % - Sf> Friction slope [m/m]
        % - dx> Raster resolution
        % Output:
        % - Q> Flow rate [m3/s]
        Q = sign(Sf).*1./n.*dx.*h.^(5/3).*sqrt(abs(Sf));
    end
end

