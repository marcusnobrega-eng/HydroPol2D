function [rainfall_raster, register, register_data] = NWP_rainfall_processing(t_d,ax_system_output,register,cycler,DEM_raster,dates,first_register)

    % Add the nctoolbox package to the MATLAB path
    if register == 0
        addpath('D:\Google_drive\Drives compartilhados\2DCAWQ_Model\Model_Functions_Directory\nctoolbox');
        setup_nctoolbox; % Packages for reading GRIB2 files (standar for NWPs according to the WMO)
        % Create a temporal folder to storage rainfall data
        mkdir('NWP_rainfall_data');
        register_data = datetime(dates(end));
    end
    
    % Donwloading the filtered data for precipitation rate surface data
    cycle = ["00","06","12","18"];
    % Extractring the date from the last update from the system
    hour = sprintf('%03d',register);
    [year,month,day] = ymd(datetime(first_register));
    year = num2str(year);month = num2str(month);day = num2str(day);
    if strlength(num2str(month))==1
        month=strcat('0',month);
    end
    if strlength(day)==1
        day=strcat('0',day);
    end
    % Converting hour time into intro for the url
    url = strcat('https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?dir=%2Fgfs.', ...
        year,month,day,'%2F',cycle(cycler),'%2Fatmos&file=gfs.t',cycle(cycler),'z.pgrb2.0p25.f',hour,'&var_PRATE=on&lev_surface=on');
    % Fetch the data from the URL
    trier = 1;
    while trier==1
        try
            % pause between images download, server restriction
            pause(15);
            file = websave('NWP_rainfall_data/file_NWP',url);
            t_d.TimerFcn{29} = strcat("GFS image No. ~",num2str(register)," donwloaded");
            trier=0;
        catch
            fprintf('Lost connection with NOAA server, retrying...');
            t_d.TimerFcn{29} = "Lost connection with NOAA server or newest file is not ready, retrying...";
            pause(300);
        end
    end
    
    % Open the GRIB2 file
    nc = ncgeodataset(file);
    % Get a list of variables in the GRIB2 file
    base_info = nc.data('Precipitation_rate_surface');
    data = reshape(base_info,721,1440).*3600; %Kg/m2/s to mm/hr we multiply by 60x60
    lats = nc.data('lat');
    lons = nc.data('lon');
    resolution = lats(1)-lats(2); % for regular meshes

    % Due to some bug about of fread() function, half of the data is
    % inverted in the horizontal axis, the invert them again.
    data= horzcat(data(:,721:1440),data(:,1:720));
    
    % Determing the timezone to correct local time
    [lat,lon] = projinv(DEM_raster.georef.SpatialRef.ProjectedCRS,DEM_raster.georef.SpatialRef.XWorldLimits,DEM_raster.georef.SpatialRef.YWorldLimits);
    zone = timezone(mean(lon));

    %Message for cycle computing
    if cycler == 1 && zone > 0
        t_d.TimerFcn{29} = strcat("Main forecast computing from ~", num2str(24 - zone), " hrs update");
    elseif cycler == 1 && zone < 0
        t_d.TimerFcn{29} = strcat("Main forecast computing from ~", num2str(0 - zone), " hrs update");
    else
        t_d.TimerFcn{29} = strcat("Main forecast computing from ~", num2str((cycler-1)*6 - zone), " hrs update");
    end
    % Create a spatial referencing object
    % Making the raster
    % EPGS 4326 for WGS 84
    crs = geocrs(4326);
    % Extent and resolution definition in Degrees
    % ref = georasterref(worldFile, rasterSize, rasterInterpretation);
    ref = georasterref('LatitudeLimits',[-90 90], ...
        'LongitudeLimits',[-180 180], ...
        'RasterSize',[721 1440], ...
        'RasterInterpretation','cells');
    ref.GeographicCRS = crs;
    ref.ColumnsStartFrom = 'north';
    % Croping the raster to the country area extent, this is to reduce
    % computational effort in future processing
    % Convert the extracted components into numbers
    register_data = datetime(dates(end)) + hours(1);

    % Creating the extent to crop the data
    [rr2, rR] = geocrop(data,ref,[lat(1)-1 lat(2)+1],[lon(1)-1 lon(2)+1]);
    % rr2 = rr2(end:-1:1,:);
    % Reproject the coordinates from EPSG:4326 to EPSG:3857
    crs = projcrs(3857);
    % Reprojecting the latitude and longitude
    [x, y] = projfwd(crs, rR.LatitudeLimits, rR.LongitudeLimits);
    % Average dimention of the pixel size acordding the ecuator
    % Accuracy could be affected according the latitude
    cs = deg2km(rR.CellExtentInLatitude)*1000;
    nx = rR.RasterSize(2); ny = rR.RasterSize(1);
    rasterSize = [ny nx];
    % Creating the reference for the cropped area
    ref2 = maprefcells(x, y, rasterSize);
    ref2.ProjectedCRS = crs;
    % Making a copy of the DEM_raster to generate the NWP raster
    rr = DEM_raster; rr.Z = rr2; rr.name = 'NWP_raster'; rr.cellsize = cs;
    rr.size = size(rr2); rr.georef.SpatialRef = ref2; rr.georef.Height = size(rr2,1);
    rr.georef.Width = size(rr2,2);
    rr.refmat = [0 -cs; cs 0; ref2.XWorldLimits(1) ref2.YWorldLimits(2)];
    rr.georef.RefMatrix = rr.refmat;

    % Resampling with the DEM resolution and nearest method
    rainfall_raster = resample(rr,DEM_raster.cellsize,'nearest');
    clear rr;
    rainfall_raster = crop(rainfall_raster,DEM_raster.georef.SpatialRef.XWorldLimits,DEM_raster.georef.SpatialRef.YWorldLimits);
    % to match with the reference DEM array size
    rainfall_raster = resample(rainfall_raster,DEM_raster);
    temp = DEM_raster; temp.Z = ~isnan(DEM_raster.Z);
    rainfall_raster = clip(rainfall_raster,temp);
    register = register + 1; % for the next hour of NWP
    clear temp file nc;
    delete('NWP_rainfall_data/file_NWP.gbx9');
    delete('NWP_rainfall_data/file_NWP.ncx');
    delete('NWP_rainfall_data/file_NWP');
end