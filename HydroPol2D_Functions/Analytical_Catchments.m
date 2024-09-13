function [h_analytical,flags,Resolution,elevation,Wshed_Properties,Inflow,idx_rivers,Elevation_Properties,d_p,outlet_index,LULC_Properties,Stage_Parameters,BC_States,Inflow_Parameters] = Analytical_Catchments(running_control,flags,Wshed_Properties,Elevation_Properties,LULC_Properties)

% %% - Non-Breaking Wave (Horizontal)
% Theoretical DEM
n = 0.005;
u = 0.635;
Resolution = 25; % m
dx = 240*Resolution; % m
dy = 32*Resolution; % m

nx = dx/Resolution;
ny = dy/Resolution;

% --- Outlet ---- %
DEM = ones(ny,nx);
slope = 0.00;
for i = 1:size(DEM,2)
    DEM(:,i) = 1 - (i-1)*slope*Resolution;
end
DEM(:,nx) = DEM(:,nx) - 0.001; % Outlet

% Adding some noise in the DEM
% variance = 0.01;
% for i = 1:size(DEM,1)
%     for j = 1:size(DEM,2)
%         DEM(i,j) = DEM(i,j) + sqrt(variance)*randn;
%     end
% end

% DEM(1,3) = nan; 
% DEM(ny,3) = nan;
% DEM(3,3) = nan;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Real World DEM
% load elevation_dam.mat
% DEM = elevation;
% Resolution = 9.98; % m
% dx = Resolution*size(DEM,2); % m
% dy = Resolution*size(DEM,1); % m
% 
% nx = dx/Resolution;
% ny = dy/Resolution;
% 
% n = 0.005; % not used
% u = 0.95; % not used


% Time for Stage Hydrograph or Inflow
Tf = 2.5*60; % min

%% Stage Hydrograph Properties
flags.flag_stage_hydrograph = 1;
Stage_Parameters.easting_inlet_cells = [1*ones(ny,1)]';
Stage_Parameters.northing_inlet_cells = [1:1:ny]';

Stage_Parameters.n_stage_gauges = 1;
Stage_Parameters.time_step_stage = 0.1/60; % minutes
Stage_Parameters.n_stage_obs = Tf/Stage_Parameters.time_step_stage;
Stage_Parameters.time_stage = Stage_Parameters.time_step_stage:Stage_Parameters.time_step_stage:(Stage_Parameters.n_stage_obs*Stage_Parameters.time_step_stage);
Stage_Parameters.time_stage = Stage_Parameters.time_stage';
t = Stage_Parameters.time_stage *60;

Stage_Parameters.stage = (7/3.*n.^2.*u^3.*t).^(3/7);

stage_cells = zeros(ny,nx);
for kk = 1:length(Stage_Parameters.easting_inlet_cells)
    stage_cells(Stage_Parameters.northing_inlet_cells(kk),Stage_Parameters.easting_inlet_cells(kk)) = 1;
end

Wshed_Properties.n_inlets_stage = sum(sum(stage_cells));


%% Inflow Hydrograph
flags.flag_inflow = 0;
% Theoretical
% Inflow_Parameters.easting_inlet_cells = [1*ones(ny,1)];
% Inflow_Parameters.northing_inlet_cells = [1:1:ny]';
% Real-World
Inflow_Parameters.easting_inlet_cells = [48,49,50,51]';
Inflow_Parameters.northing_inlet_cells = [205,203,202,200]';



Inflow_Parameters.n_stream_gauges = 1;
Inflow_Parameters.time_step_inflow = 0.1/60; % minutes
Inflow_Parameters.n_stream_obs = Tf/Inflow_Parameters.time_step_inflow;
Inflow_Parameters.time_inflow = Inflow_Parameters.time_step_inflow:Inflow_Parameters.time_step_inflow:(Inflow_Parameters.n_stream_gauges*Inflow_Parameters.time_step_inflow);
Inflow_Parameters.time_inflow = Inflow_Parameters.time_inflow';
running_control.steps = Tf/running_control.time_step_model;
Inflow_Parameters.inflow_discharge = 49000*ones(Inflow_Parameters.n_stream_obs,1); % m3/s

Wshed_Properties.n_inlets = 0;
Wshed_Properties.n_inlets = length(Inflow_Parameters.easting_inlet_cells);
Wshed_Properties.inflow_mask = zeros(ny,nx);
% Here we read where the stream gauges are located
if flags.flag_inflow == 1
    % Wshed_Properties.inflow_mask = zeros(size(DEM_raster.Z));
    Wshed_Properties.inflow_cells = zeros(ny,nx,Inflow_Parameters.n_stream_gauges);
    for i = 1:Inflow_Parameters.n_stream_gauges
        for z = 1:Wshed_Properties.n_inlets(i)
            x = Inflow_Parameters.easting_inlet_cells(z,i);
            y = Inflow_Parameters.northing_inlet_cells(z,i);
            Wshed_Properties.inflow_cells(y,x,i) = 1; % Coordinates of each stream gauge
            Wshed_Properties.inflow_mask(y,x) = 1;
        end
    end
    Wshed_Properties.inflow_mask = logical(Wshed_Properties.inflow_mask);
end
BC_States.inflow = zeros(ny,nx); % This is to solve spatially, don't delete

if flags.flag_inflow == 1
    for i = 1:Inflow_Parameters.n_stream_gauges
        BC_States.inflow = BC_States.inflow + 0*Wshed_Properties.inflow_cells(:,:,i); % mm
    end
end

Inflow_Parameters.inflow_hydrograph_rate = Inflow_Parameters.inflow_discharge./(sum(sum(Wshed_Properties.inflow_mask))*Resolution^2)*1000*3600; % mm/h

Inflow_Parameters.inflow_hydrograph_rate = Inflow_Parameters.inflow_hydrograph_rate*flags.flag_inflow;

if flags.flag_inflow == 1
    inflow_length = size(Inflow_Parameters.inflow_hydrograph_rate,1); % Number of interval
    inflow_discretized = zeros(size(Inflow_Parameters.inflow_hydrograph_rate,2),ceil(inflow_length*Inflow_Parameters.time_step_inflow/running_control.time_step_model)); % Preallocating
    for z = 1:Inflow_Parameters.n_stream_gauges
        for i =1:((inflow_length)*Inflow_Parameters.time_step_inflow/running_control.time_step_model)
            inflow_discretized(z,i) = Inflow_Parameters.inflow_hydrograph_rate(ceil((round(i*running_control.time_step_model/Inflow_Parameters.time_step_inflow,12))),z); % Discretized into moldel's time-step
        end
    end
end

if flags.flag_inflow == 1
    [~,BC_States.delta_inflow,inflow_intensity] = accumulated_incremental(running_control.steps,inflow_discretized,running_control.time_step_model);
    BC_States.time_deltainflow = cumsum(ones(size(BC_States.delta_inflow,2),1)*running_control.time_step_model); % Vector with indices
    BC_States.delta_inflow = BC_States.delta_inflow;
end


%% Extra Parameters
Wshed_Properties.cell_area = Resolution^2;
Wshed_Properties.stage_mask = logical(stage_cells);
Wshed_Properties.roughness = n*ones(ny,nx);
LULC_Properties.roughness = n*ones(ny,nx);
flags.flag_rainfall = 0;
elevation = DEM;
Wshed_Properties.Resolution = Resolution;
Wshed_Properties.rainfall_matrix = ones(ny,nx);
Inflow = ones(ny,nx);
idx_rivers = zeros(ny,nx);
Elevation_Properties.elevation_cell = DEM;
BC_States.delta_p_agg = 0;
BC_States.outflow_volume = 0;
BC_States.inflow_volume = 0;
d_p = zeros(ny,nx);
outlet_index = DEM == min(min(DEM));
[Wshed_Properties.row_outlet,Wshed_Properties.col_outlet] = find(outlet_index == 1); % finding rows and cols of these cells

% Solving Analytical Values
x = Resolution:Resolution:Resolution*size(DEM,2);
time = [2700 5400 9000]; % sec;
for i = 1:length(time)
    h_analytical(:,i) = (-7/3*(n^2*u^2.*(min(x - u.*time(i),0)))).^(3/7);
end
%% - Non-Breaking Wave (Vertical)
% n = 0.005;
% u = 0.635;
% Resolution = 25; % m
% dx = 32*Resolution; % m
% dy = 240*Resolution; % m
%
% nx = dx/Resolution;
% ny = dy/Resolution;
%
% % --- Outlet ---- %
% DEM = ones(ny,nx);
% slope = 0.00;
% for i = 1:size(DEM,2)
%     DEM(i,:) = 1 - (i-1)*slope*Resolution;
% end
% DEM(ny,:) = DEM(ny,:) - 0.001; % Outlet
%
% % DEM(round(ny/3):round(2*ny/3),round(nx/3):round(2*nx/3)) = nan;
%
% % DEM(25,2:3) = nan;
%
% % Time for Stage Hydrograph or Inflow
% Tf = 2.5*60; % min
%
% %% Stage Hydrograph Properties
% flags.flag_stage_hydrograph = 1;
%
% flags.flag_waterbalance = 0;
%
% Stage_Parameters.easting_inlet_cells = [1:1:nx]';
% Stage_Parameters.northing_inlet_cells = [1*ones(nx,1)];
% Stage_Parameters.n_stage_gauges = 1;
% Stage_Parameters.time_step_stage = 0.1/60; % minutes
% Stage_Parameters.n_stage_obs = Tf/Stage_Parameters.time_step_stage;
% Stage_Parameters.time_stage = Stage_Parameters.time_step_stage:Stage_Parameters.time_step_stage:(Stage_Parameters.n_stage_obs*Stage_Parameters.time_step_stage);
% Stage_Parameters.time_stage = Stage_Parameters.time_stage';
% t = Stage_Parameters.time_stage *60;
%
% Stage_Parameters.stage = (7/3.*n.^2.*u^3.*t).^(3/7);
%
% stage_cells = zeros(ny,nx);
% for kk = 1:length(Stage_Parameters.easting_inlet_cells)
%     stage_cells(Stage_Parameters.northing_inlet_cells(kk),Stage_Parameters.easting_inlet_cells(kk)) = 1;
% end
%
% Wshed_Properties.n_inlets_stage = sum(sum(stage_cells));
%
%
% %% Inflow Hydrograph
% flags.flag_inflow = 0;
% Inflow_Parameters.easting_inlet_cells = [1:1:nx]';
% Inflow_Parameters.northing_inlet_cells = [1*ones(nx,1)];
%
% Inflow_Parameters.n_stream_gauges = 1;
% Inflow_Parameters.time_step_inflow = 0.1/60; % minutes
% Inflow_Parameters.n_stream_obs = Tf/Inflow_Parameters.time_step_inflow;
% Inflow_Parameters.time_inflow = Inflow_Parameters.time_step_inflow:Inflow_Parameters.time_step_inflow:(Inflow_Parameters.n_stream_gauges*Inflow_Parameters.time_step_inflow);
% Inflow_Parameters.time_inflow = Inflow_Parameters.time_inflow';
% running_control.steps = Tf/running_control.time_step_model;
%
% Inflow_Parameters.inflow_discharge = 40*ones(Inflow_Parameters.n_stream_obs,1); % m3/s
%
% Wshed_Properties.n_inlets = 0;
% Wshed_Properties.n_inlets = length(Inflow_Parameters.easting_inlet_cells);
% % Here we read where the stream gauges are located
% Wshed_Properties.inflow_mask = zeros(ny,nx);
% if flags.flag_inflow == 1
%     % Wshed_Properties.inflow_mask = zeros(size(DEM_raster.Z));
%     Wshed_Properties.inflow_cells = zeros(ny,nx,Inflow_Parameters.n_stream_gauges);
%     for i = 1:Inflow_Parameters.n_stream_gauges
%         for z = 1:Wshed_Properties.n_inlets(i)
%             x = Inflow_Parameters.easting_inlet_cells(z,i);
%             y = Inflow_Parameters.northing_inlet_cells(z,i);
%             Wshed_Properties.inflow_cells(y,x,i) = 1; % Coordinates of each stream gauge
%             Wshed_Properties.inflow_mask(y,x) = 1;
%         end
%     end
%     Wshed_Properties.inflow_mask = logical(Wshed_Properties.inflow_mask);
% end
% BC_States.inflow = zeros(ny,nx); % This is to solve spatially, don't delete
% if flags.flag_inflow == 1
%     BC_States.inflow = zeros(ny,nx); % This is to solve spatially, don't delete
%     for i = 1:Inflow_Parameters.n_stream_gauges
%         BC_States.inflow = BC_States.inflow + 0*Wshed_Properties.inflow_cells(:,:,i); % mm
%     end
% end
%
% Inflow_Parameters.inflow_hydrograph_rate = Inflow_Parameters.inflow_discharge./(sum(sum(Wshed_Properties.inflow_mask))*Resolution^2)*1000*3600; % mm/h
%
% Inflow_Parameters.inflow_hydrograph_rate = Inflow_Parameters.inflow_hydrograph_rate*flags.flag_inflow;
%
% if flags.flag_inflow == 1
%     inflow_length = size(Inflow_Parameters.inflow_hydrograph_rate,1); % Number of interval
%     inflow_discretized = zeros(size(Inflow_Parameters.inflow_hydrograph_rate,2),ceil(inflow_length*Inflow_Parameters.time_step_inflow/running_control.time_step_model)); % Preallocating
%     for z = 1:Inflow_Parameters.n_stream_gauges
%         for i =1:((inflow_length)*Inflow_Parameters.time_step_inflow/running_control.time_step_model)
%             inflow_discretized(z,i) = Inflow_Parameters.inflow_hydrograph_rate(ceil((round(i*running_control.time_step_model/Inflow_Parameters.time_step_inflow,12))),z); % Discretized into moldel's time-step
%         end
%     end
% end
%
% if flags.flag_inflow == 1
%     [~,BC_States.delta_inflow,inflow_intensity] = accumulated_incremental(running_control.steps,inflow_discretized,running_control.time_step_model);
%     BC_States.time_deltainflow = cumsum(ones(size(BC_States.delta_inflow,2),1)*running_control.time_step_model); % Vector with indices
%     BC_States.delta_inflow = BC_States.delta_inflow;
% end
%
% %% Extra Parameters
% Wshed_Properties.cell_area = Resolution^2;
% Wshed_Properties.stage_mask = logical(stage_cells);
% Wshed_Properties.roughness = n*ones(ny,nx);
% LULC_Properties.roughness = n*ones(ny,nx);
% flags.flag_rainfall = 0;
% elevation = DEM;
% Wshed_Properties.Resolution = Resolution;
% Wshed_Properties.rainfall_matrix = ones(ny,nx);
% Inflow = ones(ny,nx);
% idx_rivers = zeros(ny,nx);
% Elevation_Properties.elevation_cell = DEM;
% BC_States.delta_p_agg = 0;
% BC_States.outflow_volume = 0;
% BC_States.inflow_volume = 0;
% d_p = zeros(ny,nx);
% outlet_index = DEM == min(min(DEM));
% [Wshed_Properties.row_outlet,Wshed_Properties.col_outlet] = find(outlet_index == 1); % finding rows and cols of these cells
%
%
% % Solving Analytical Values
% x = Resolution:Resolution:Resolution*size(DEM,2);
% time = [2700 5400 9000]; % sec;
% for i = 1:length(time)
%     h_analytical(:,i) = (-7/3*(n^2*u^2.*(min(x - u.*time(i),0)))).^(3/7);
% end
end

