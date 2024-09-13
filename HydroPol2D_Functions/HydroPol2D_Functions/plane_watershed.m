%% Watershed, LULC, and Soil Maps for a Plane
% Developer: Marcus Nobregs
% This is the code to create a 1-D plan with slope s0

function [DEM,SOIL,LULC] = plane_watershed(pixelsize,s0,width,length,elevation_end)
% s0 in degrees

% Input Data
% pixelsize = 0.02; % Meters
% s0 = 0.01; % deg
% width = 1.48; % width of the plane
% length = 2.96; % length of the plane
% elevation_end = 0; % elevation of the outlet cells (m)
nx = width/pixelsize;
ny = length/pixelsize;

% Converting From Degrees to Percentage
s0 = tan(deg2rad(s0));


% Elevation Data
DEM = zeros(ny,nx);
for i = 1:nx
    for j = 1:ny
        if j == 1
            DEM(j,i) = elevation_end;
        else
            DEM(j,i) = DEM(j-1,i) + s0*pixelsize;
        end
    end
end

SOIL = ones(size(DEM));

LULC = ones(size(DEM));

% grid
x_grid = [1:1:nx]*pixelsize;
y_grid = [1:1:ny]*pixelsize;
[X,Y] = meshgrid(x_grid,y_grid);

%% SurfPlot
% sf = 12;
% h = surf(X,Y,DEM);
% view(0,90);
% colormap(linspecer); view(0,90);
% c = colorbar;
% xlabel('x (m)','Interpreter','latex','FontSize',sf)
% ylabel('y (m)','Interpreter','latex','FontSize',sf)
% ylabel(c,'Elevation (m)','Interpreter','Latex','FontSize',sf)
% axis tight
% axis equal


%% Export Rasters
xllcorner = pixelsize/2;
yllcorner = 0;

SaveAsciiRaster(DEM,pixelsize,xllcorner,yllcorner,'DEM_Plane',-9999);
SaveAsciiRaster(SOIL,pixelsize,xllcorner,yllcorner,'SOIL_Plane',-9999);
SaveAsciiRaster(LULC,pixelsize,xllcorner,yllcorner,'LULC_Plane',-9999);
end

