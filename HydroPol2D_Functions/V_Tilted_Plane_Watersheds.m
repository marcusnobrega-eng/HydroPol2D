function [dem] = V_Tilted_Plane_Watersheds(Delta_x,length_x,channel_x,length_y,datum_zero,slope_x,slope_y)
%% V-Tilted Watershed
% Goal - Assign parameters for cells according to the land use and land
% cover given by the imperviousness map

%% variables to determine the DEM
% Delta_x = 20;
a_grid = Delta_x; % size of the grid cell
% length_x = 800; % length of the hillslope
% channel_x = 20; % width of the channel
% length_y = 1000; % length of the channel
% datum_zero = 0; %m - min level
resolution = a_grid; %m
% slope_y = 0.02; %i %slope in y direction
% slope_x = 0.05; %j %slope in x direction
channel_grid = channel_x;

%% Variables to assign to each cell: Soil Properties, Roughness Coefficient and other properties
%%%% Testing Overland Flow %%%%
% n_per = 0.015;
% n_imp = 0.15;
% ksat_per = 0.0001;
% ksat_imp = 0.0001;
% d_0_per = 0;
% d_0_imp = 0;
% h_0_per = 0;
% h_0_imp = 0;
% psi_per = 0.0001;
% psi_imp = 0.0001;
% I_0_per = 1000;
% I_0_imp = 1000;
% teta_sat_per = 0.09;
% teta_sat_imp = 0.09;
% teta_i_per = 0.08;
% teta_i_imp = 0.08;

%%%% Testing infiltration %%%%
% n_per = 0.3;
% n_imp = 0.018;
% ksat_per = 10.922;
% ksat_imp = 2;
% d_0_per = 0;
% d_0_imp = 0;
% h_0_per = 10;
% h_0_imp = 0.5;
% psi_per = 110;
% psi_imp = 5;
% I_0_per = 10;
% I_0_imp = 10;
% teta_sat_per = 0.454;
% teta_sat_imp = 0.2;
% teta_i_per = 0.01;
% teta_i_imp = 0.1;



%% Calculations
%%% Preallocating arrays

y_1 = (length_y)/resolution;
x_1 = (length_x + channel_x)/resolution;
elevation_ = zeros(y_1,x_1); roughness_ = zeros(y_1,x_1);
desloc_j = (length_x+channel_x)/resolution;

% For loops for each plane
for i=1:(length_y)/resolution
    for j =1:(length_x + channel_x)/resolution
        if i==1
            if (j)*resolution <= channel_grid
                elevation_(i,j) = datum_zero;
             
            else                
                elevation_(i,j) = datum_zero + slope_x*(j-1)*resolution;    
            end           
        else
            if (j)*resolution <= channel_grid
                elevation_(i,j) = elevation_(i-1,j) + resolution*slope_y;
            else
                elevation_(i,j) = elevation_(i-1,j) + resolution*slope_y;
            end
        end
    end
end
for i = 1:(length_y)/resolution
    for j = 1:(2*desloc_j - 1)
        if j<=desloc_j
            elevation(i,j +desloc_j-1) = elevation_(i,j);
        else
            number_j = j-desloc_j;
            elevation(i,desloc_j - number_j) = elevation(i,j);
        end
    end
end
% Chaning names to match with the main file
dem = elevation;

% GRIDobj2geotiff(dem,'DEM_V_Tilted.tif')
% GRIDobj2geotiff(dem*0,'Initial_Depths_V_Tilted.tif');
% GRIDobj2geotiff(dem./dem*5,'Initial_Soil_Moisture_V_Tilted.tif');
% GRIDobj2geotiff(dem./dem,'LULC_V_Tilted.tif');
% GRIDobj2geotiff(dem./dem,'SOIL_Tilted.tif');




end