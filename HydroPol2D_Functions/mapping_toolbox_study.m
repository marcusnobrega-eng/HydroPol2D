% Studying Geo Maps
close all

land = readgeotable("landareas.shp");
rivers = readgeotable("worldrivers.shp");
cities = readgeotable("worldcities.shp");

%% Example 1 
figure
worldmap Brazil

load coastlines
plotm(coastlat,coastlon)

figure
worldmap europe

ax = gca;
getm(ax,"MapProjection")

%% Display data
geoshow(land,"FaceColor",[0.88 0.95 0.81])
geoshow(rivers,"Color",[0 0.4470 0.7410])
geoshow(cities,"Marker",".","MarkerEdgeColor",[0.6350 0.0780 0.1840])

% Label
textm(35,14,"Mediterranean Sea")
