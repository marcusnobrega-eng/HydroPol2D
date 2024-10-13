%% Normalized Discharge gauges%%

zero_matrix = zeros(size(elevation,1),size(elevation,2));
if flags.flag_obs_gauges == 1 && flags.flag_rainfall == 1
    % Rainfall Std Deviation
    if flags.flag_spatial_rainfall == 1
        store=1;
        flag_loader=1;
        for i = 1:length(running_control.time_records)
            if i > saver_memory_maps*store
                store = store + 1;
                load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
            else
                if flag_loader == 1
                  load(strcat('Temporary_Files\save_map_hydro_',num2str(store)),'Maps');
                  flag_loader=0;
                end
            end
            for j = 1:length(gauges.easting_obs_gauges)
                zgauges = Maps.Hydro.spatial_rainfall_maps(:,:,i - ((store-1)*saver_memory_maps)).*subcatchments{1,j};
                Rainfall_Parameters.std_dev_gauges{1,j}(i,1) = nanstd(zgauges(:));
            end
        end    
    end
    % Catchment area of each gauge and graphics
    for i = 1:length(gauges.easting_obs_gauges)
        gauges.catchment_area(i,1) = Wshed_Properties.fac_area(gauges.northing_obs_gauges(i,1),gauges.easting_obs_gauges(i,1)); % km2
    end

    % Creating the folder
    mkdir(strcat(folderName,'\Normalized_Discharge_gauges'));        
    % Specify the folder where the files live.
    myFolder_wd = strcat(pwd,'\',folderName,'\Normalized_Discharge_gauges'); % Current folder
    % Check to make sure that folder actually exists.  Warn user if it doesn't.
    if ~isfolder(myFolder_wd)
        errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder_wd);
        uiwait(warndlg(errorMessage));
        return;
    end    
    
    % Plot the graphics
    labels_depth = gauges.labels_observed_string; labels_depth{gauges.num_obs_gauges+1} = 'Outlet';
    labels_depth{gauges.num_obs_gauges + 2} = 'Rainfall Intensity';

    for i = 1:length(gauges.easting_obs_gauges)
        figure;
        subplot(2, 1, 1);
        imageschs(DEM_raster, DEM_raster.*subcatchments{1,i});
        title(labels_depth(i))
        xlabel('Longitude');
        ylabel('Latitude');
        axis equal; % Keep the aspect ratio    

        subplot(2, 1, 2);
        set(gcf,'units','inches','position',[3,0,9,5]);
        set(gca,'ycolor','black');
        specific_discharge = gauges.hydrograph_cell(:,i)/gauges.catchment_area(i,1) ; % m3/s/km2
        plot(gather(running_control.time_hydrograph),specific_discharge,'LineWidth',1.5,'linestyle',ls);
        hold on
    
        % Outlet
        %plot(gather(running_control.time_hydrograph),gather(outlet_states.outlet_hydrograph)/(Wshed_Properties.drainage_area/1000/1000),'LineWidth',1.5,'linestyle',ls,'marker','.','color','red');
        xlabel('Elapsed Time (min)','interpreter','latex'); ylabel('Specific Discharge $(\mathrm{m^3/s/km^2})$','interpreter','latex'); set(gca,'FontSize',12);
        hold on
        yyaxis right; set(gca,'ydir','reverse','ycolor','black');
        if flags.flag_rainfall == 1
            if flags.flag_spatial_rainfall ~=1
                bar(gather(Rainfall_Parameters.time_rainfall),gather(Rainfall_Parameters.intensity_rainfall),'FaceColor',[0 .55 .55],'EdgeColor',[0 .5 .5],'LineWidth',1.5)
                ylabel('Rainfall Intensity (mm/h)','interpreter','latex');
                ylim([0 max(gather(Rainfall_Parameters.intensity_rainfall))*6])
            else
                bar(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)),gather(BC_States.average_spatial_rainfall_gauges{1,i}),'FaceColor',[0 .5 .5],'EdgeColor',[0 .55 .55],'LineWidth',1.5);
                ylabel('Aerial Mean Rainfall Intensity (mm/h)','interpreter','latex');
                hold on
                try
                    er = errorbar(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)),BC_States.average_spatial_rainfall_gauges{1,i},Rainfall_Parameters.std_dev_gauges{1,i}(1:(dim),1),Rainfall_Parameters.std_dev_gauges{1,i}(1:(dim),1));
                    er.Color = [0 0 0];
                    er.LineStyle = 'none';
                end
                plot((gather(Spatial_Rainfall_Parameters.rainfall_spatial_duration(1,1:(dim)))),gather(BC_States.average_spatial_rainfall_gauges{1,i}),'LineWidth',1.5,'color','blue')
                ylim([0 max(max(gather(BC_States.average_spatial_rainfall_gauges{1,i})))*6])
            end
        end
        
        legend([labels_depth{i},labels_depth(end),'Standar Dev'],'FontName','Garamond','FontSize',8,'location','bestoutside')

        title('Specific Discharge','interpreter','latex','fontsize',12)
        set(gca, 'TickLength', [0.02 0.01]);
        set(gca,'Tickdir','out')
        set(gca,'FontName','Garamond','FontSize',12)
        try
            exportgraphics(gcf,fullfile(myFolder_wd,'Specific_Discharge_Gauges_'+labels_depth{i}+'.pdf'),'ContentType','vector')
        catch
            fprintf('Specific discharge gauges no exported, PDF export error')
        end
        saveas(gcf,fullfile(myFolder_wd,'Specific_Discharge_Gauges_'+labels_depth{i}+'.fig'))
        close all
    end
end    
   

