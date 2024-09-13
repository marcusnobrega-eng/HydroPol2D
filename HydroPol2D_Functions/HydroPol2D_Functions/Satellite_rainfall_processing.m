function [rainfall_raster, register, register_data, register_data_2] = Satellite_rainfall_processing(t_d,ax_system_output,register, product,date_begin,date_end,offline_flag,real_time_flag,DEM_raster)
%% Satellite data processing

% If the model is set up for off-line processing,
% register is the variable that control with witch file we are working
if register==0
    % Create a temporal folder to storage rainfall data
    mkdir('Spatial_rainfall_data');
end

if offline_flag == 1 && real_time_flag == 1
    error('There is no way to run both cases together. Please choose either offline or real-time modeling with satellite data.');
end

% satellite data is downloaded once
if offline_flag == 1 % We are running with satellite data - offline
    if register == 0
        trier_e = 1;
        while trier_e == 1 
            try
                % Stablish the File Tranfer Protocol
                FTP = ftp('ftp://persiann.eng.uci.edu');
                productID = product;
                % Reading the strar and end date
                [s_y,s_m,s_d] = ymd(date_begin);
                [e_y,e_m,e_d] = ymd(date_end);
                s_h = date_begin.Hour;s_mn = minute(date_begin);
                e_h = date_end.Hour;e_mn = minute(date_end);
                % extranting year month and day, then converting into strings
                % s_m
                s_y=num2str(s_y);s_m=num2str(s_m);s_d=num2str(s_d);
                e_y=num2str(e_y);e_m=num2str(e_m);e_d=num2str(e_d);
                s_h = num2str(s_h);e_h = num2str(e_h);
                s_mn = num2str(s_mn);e_mn = num2str(e_mn);
                if strlength(s_m)==1
                    s_m=strcat('0',s_m);
                end
                if strlength(e_m)==1
                    e_m=strcat('0',e_m);
                end
                if strlength(s_d)==1
                    s_d=strcat('0',s_d);
                end
                if strlength(e_d)==1
                    e_d=strcat('0',e_d);
                end
                if strlength(s_h)==1
                    s_h =strcat('0',s_h);
                end
                if strlength(e_h)==1
                    e_h =strcat('0',e_h);
                end            
                if strlength(s_mn)==1
                    s_mn=strcat('0',s_mn);
                end
                if strlength(e_mn)==1
                    e_mn=strcat('0',e_mn);
                end
                
                % Changing the directoy of the FTP for the desired product required
                directory = append('CHRSdata/PDIRNow/',productID,'/',num2str(s_y)); cd(FTP,directory);
                % Extract the directory with all data available
                list = dir(FTP); temp_list = {list.name};
                % Finding the index from the list of the star and end date
                star_name_founder = strcat('pdirnow1h',s_y(3:4),s_m,s_d,s_h,'.bin.gz');
                end_name_founder = strcat('pdirnow1h',e_y(3:4),e_m,e_d,e_h,'.bin.gz');
                start_name = find(ismember(temp_list, star_name_founder))+1; % the +1 is to gather the forward rainfall
                end_name = find(ismember(temp_list, end_name_founder))+1;
                
                % Downloading the data
                finder = temp_list(start_name); 
                downloaded_dir = dir('Spatial_rainfall_data'); temp_downloaded = {downloaded_dir.name};
                if register == 0
                    register_data = 1; % Please check it, Luis.
                else
                    register_data = register; % Please check it, Luis.
                end
                if sum(ismember(temp_downloaded,list(register_data,1).name))==0
                    mget(FTP,list(start_name,1).name,'Spatial_rainfall_data');
                end
                trier_e = 0;
                register = start_name +1;
            catch
                fprintf('FTP lost connection');
                fprintf('trying again...');
                pause(5)
            end
        end
    else
        trier_e = 1;
        while trier_e == 1 
            try
                % Stablish the File Tranfer Protocol
                FTP = ftp('ftp://persiann.eng.uci.edu');
                productID = product;
                % Reading the strar and end date
                [s_y,s_m,s_d] = ymd(date_begin);
                [e_y,e_m,e_d] = ymd(date_end);
                s_h = date_begin.Hour;
                e_h = date_end.Hour;
                s_mn = minute(date_begin);
                e_mn = minute(date_end);
                % extranting year month and day, then converting into strings
                % s_m
                s_y=num2str(s_y);s_m=num2str(s_m);s_d=num2str(s_d);
                e_y=num2str(e_y);e_m=num2str(e_m);e_d=num2str(e_d);
                s_h = num2str(s_h);e_h = num2str(e_h);
                s_mn = num2str(s_mn);e_mn = num2str(e_mn);
                if strlength(s_m)==1
                    s_m=strcat('0',s_m);
                end
                if strlength(e_m)==1
                    e_m=strcat('0',e_m);
                end
                if strlength(s_d)==1
                    s_d=strcat('0',s_d);
                end
                if strlength(e_d)==1
                    e_d=strcat('0',e_d);
                end
                if strlength(s_h)==1
                    s_h =strcat('0',s_h);
                end
                if strlength(e_h)==1
                    e_h =strcat('0',e_h);
                end            
                if strlength(s_mn)==1
                    s_mn=strcat('0',s_mn);
                end
                if strlength(e_mn)==1
                    e_mn=strcat('0',e_mn);
                end
                
                % Changing the directoy of the FTP for the desired product required
                directory = append('CHRSdata/PDIRNow/',productID,'/',num2str(s_y)); cd(FTP,directory);
                % Extract the directory with all data available
                list = dir(FTP); temp_list = {list.name};
                finder = temp_list(register);
                downloaded_dir = dir('Spatial_rainfall_data'); temp_downloaded = {downloaded_dir.name};
                if sum(ismember(temp_downloaded,list(register,1).name))==0
                    mget(FTP,list(register,1).name,'Spatial_rainfall_data');
                end
                trier_e = 0;
                register = register+1;
             catch
                fprintf('FTP lost connection');
                fprintf('trying again...');
                pause(5)
             end
        end
    end
end

if real_time_flag == 1
    if register == 0
        try
            % Stablish the File Tranfer Protocol
            FTP = ftp('ftp://persiann.eng.uci.edu');
            % Reading the strar and end date
            [s_y,s_m,s_d] = ymd(date_begin);
            s_h = date_begin.Hour;
            % extranting year month and day, then converting into strings
            s_y=num2str(s_y);s_m=num2str(s_m);s_d=num2str(s_d);s_h=num2str(s_h);
            if strlength(s_m)==1
                s_m=strcat('0',s_m);
            end
            if strlength(s_d)==1
                s_d=strcat('0',s_d);
            end
            if strlength(s_h)==1
                s_h =strcat('0',s_h);
            end
            productID = product;
            % Changing the directoy of the FTP for the desired product required
            directory = append('CHRSdata/PDIRNow/',productID,'/',num2str(s_y)); cd(FTP,directory);
            % Extract the directory with all data available
            list = dir(FTP); temp_list = {list.name};
            % Finding the index from the list of the star and end date
            star_name_founder = strcat('pdirnow1h',s_y(3:4),s_m,s_d,s_h,'.bin.gz');
            % Reading the last file in the directory
            end_name_founder = temp_list(end);
            start_name = find(ismember(temp_list, star_name_founder));
            end_name = find(ismember(temp_list, end_name_founder));  
            % Downloading the data
            % register controls wich file will be donwloaded, its set up with the
            % index of the start_date, then it will increase ++1 and could inform
            % how long the model is running.
            register = start_name;
            % downloading the data
            downloaded_dir = dir('Spatial_rainfall_data'); temp_downloaded = {downloaded_dir.name};
            if sum(ismember(temp_downloaded,list(register,1).name))==0
                 mget(FTP,list(register,1).name,'Spatial_rainfall_data');
            end
            finder = temp_list(register);
            register = register+1;
        catch
            fprintf('FTP lost connection');
            fprintf('trying again...');
            pause(5)
        end
    else
        trier_e = 1;
        while trier_e == 1 
            try
                % Stablish the File Tranfer Protocol
                FTP = ftp('ftp://persiann.eng.uci.edu');
                % Reading the strar and end date
                [s_y,s_m,s_d] = ymd(date_begin);
                % extranting year month and day, then converting into strings
                s_y=num2str(s_y);s_m=num2str(s_m);s_d=num2str(s_d);
                if strlength(s_m)==1
                    s_m=strcat('0',s_m);
                end
                if strlength(s_d)==1
                    s_d=strcat('0',s_d);
                end
                productID = product;
                % Changing the directoy of the FTP for the desired product required
                directory = append('CHRSdata/PDIRNow/',productID,'/',num2str(s_y)); cd(FTP,directory);
                % Extract the directory with all data available
                list = dir(FTP); temp_list = {list.name};
                % Reading the last file in the directory
                end_name_founder = temp_list(end);
                end_name = find(ismember(temp_list, end_name_founder));  
                trier_e = 0;
                if register == end_name+1
                    trier = 1;
                    while trier == 1
                        try
                            % Updating the directory
                            FTP = ftp('ftp://persiann.eng.uci.edu');
                            directory = append('CHRSdata/PDIRNow/',productID,'/',num2str(s_y)); cd(FTP,directory);
                            list = dir(FTP); temp_list = {list.name};
                            % Finding the next file with register
                            finder = temp_list(register);
                            % downloading the data
                            mget(FTP,list(register,1).name,'Spatial_rainfall_data');
                            register = register+1;
                            trier = 0;
                        catch
                            fprintf('Waiting for the next satellite update...');
                            system_ouput="Waiting for the next satellite update...";
                            flag_message=1;
                            HydroPol2D_real_time_dashboard([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],t_d,system_ouput,ax_system_output,flag_message,[],[],0,0);
                            % Save variables from the main workspace
                            % evalin('caller', 'save(''HydroPol2D_Main_While_luis.mat'')');
                            % % Save variables from this function's workspace
                            pause(600)
                        end
                    end
                else
                    % downloading the data
                    finder = temp_list(register);
                    downloaded_dir = dir('Spatial_rainfall_data'); temp_downloaded = {downloaded_dir.name};
                    if sum(ismember(temp_downloaded,list(register,1).name))==0
                        mget(FTP,list(register,1).name,'Spatial_rainfall_data');
                    end
                    register = register+1;
                    system_ouput="Main model routing computing...";
                    flag_message=1;
                    HydroPol2D_real_time_dashboard([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],t_d,system_ouput,ax_system_output,flag_message,[],[],0,0);
                end
            catch
                fprintf('FTP lost connection');
                fprintf('trying again...');
                system_ouput="FTP lost connection. Trying again...";
                flag_message=1;
                HydroPol2D_real_time_dashboard([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],t_d,system_ouput,ax_system_output,flag_message,[],[],0,0);
                pause(10)
            end
        end
    end
end

    %%%%% Reading PDIR-now data for the study area %%%%%%

    register_data = finder;
    
    % Unzip the data
    gunzip(strcat('Spatial_rainfall_data\',finder));
    % Open the binary file
    fid = fopen(strcat('Spatial_rainfall_data\',finder{1}(1:end-3)),'r'); 
    % Read the binary data with the desired format
    data = fread(fid,[9000 3000],'int16','l')';
    % close the binary file
    fclose(fid);
    % Delete the unziped file
    delete(strcat('Spatial_rainfall_data\',finder{1}(1:end-3)))
    % -99 data is NoData
    data(data<0) = NaN;
    % Scaled data according data provider
    data = data / 100;
    % Due to some bug about of fread() function, half of the data is
    % inverted in the horizontal axis, the invert them again.
    data= horzcat(data(:,4501:9000),data(:,1:4500));

    % Making the raster
    % EPGS 4326 for WGS 84 
    crs = geocrs(4326);
    % Extent and resolution definition in Degrees
    % Making the base matrix for the raster
    ref2 = georasterref('LatitudeLimits',[-60 60], ...
        'LongitudeLimits',[-180 180], ...
        'RasterSize',[3000 9000], ...
        'RasterInterpretation','cells');
    ref2.GeographicCRS = crs; ref2.ColumnsStartFrom = 'north';
    % Croping the raster to the country area extent, this is to reduce
    % computational effort in future processing
    [lat,lon] = projinv(DEM_raster.georef.SpatialRef.ProjectedCRS,DEM_raster.georef.SpatialRef.XWorldLimits,DEM_raster.georef.SpatialRef.YWorldLimits);
    % Determing the timezone to correct local time
    zone = timezone(mean(lon));
    % Convert the extracted components into numbers
    finder = finder{1};
    year = str2double(finder(10:11))+2000;
    month = str2double(finder(12:13));
    day = str2double(finder(14:15));
    hour = str2double(finder(16:17));
    % Create a datetime object
    register_data_2 = datetime(year, month, day, hour, 0, 0) - hours(zone);
    % Creating the extent to crop the data
    [rr2, rR] = geocrop(data,ref2,[lat(1)-1 lat(2)+1],[lon(1)-1 lon(2)+1]);
    % rr2 = rr2(end:-1:1,:);
    % Reproject the coordinates from EPSG:4326 to EPSG:3857
    crs = projcrs(3857);
    
    %Cutting the raster according the DEM_raster
    rainfall_raster = raster_cutter(DEM_raster,rR,rr2,1);
    clear rr;
    clear temp FTP;

end

