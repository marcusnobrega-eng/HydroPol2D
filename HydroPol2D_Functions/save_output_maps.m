% Save Output Maps
% Goal: Save model output maps
% Developer: Marcus Nobrega, Ph.D.

% Maps with generally coarser resolution
recording_parameters.actual_record_state = find(running_control.time_records < t_save,1,'last');
recording_parameters.delta_record = recording_parameters.actual_record_state - recording_parameters.last_record_maps;
recording_parameters.last_record_maps = recording_parameters.actual_record_state;

if flags.flag_automatic_calibration ~= 1
    if k == 1 % First time-step
        Maps.Hydro.d(:,:,1) = depths.d_t; 
        if flags.flag_human_instability == 1
            Maps.Hydro.risk(:,:,1) = Human_Instability.risk_t;
        elseif flags.flag_human_instability == 2
        elseif flags.flag_human_instability == 3
            Maps.Hydro.risk_cm(:,:,1) = Human_Instability.risk_t_cm;
            Maps.Hydro.risk_tm(:,:,1) = Human_Instability.risk_t_tm;
            Maps.Hydro.risk_am(:,:,1) = Human_Instability.risk_t_am;
            Maps.Hydro.risk_om(:,:,1) = Human_Instability.risk_t_om;
            Maps.Hydro.risk_cf(:,:,1) = Human_Instability.risk_t_cf;
            Maps.Hydro.risk_tf(:,:,1) = Human_Instability.risk_t_tf;
            Maps.Hydro.risk_af(:,:,1) = Human_Instability.risk_t_af;
            Maps.Hydro.risk_of(:,:,1) = Human_Instability.risk_t_of;
        end
        Maps.Hydro.I_t(:,:,1) = Soil_Properties.I_t;
    elseif recording_parameters.delta_record > 0 % Saving maps
        t_store = recording_parameters.actual_record_state;
        % Maps.Hydro.d(:,:,t_store) = depths.d_t;
        Maps.Hydro.d(:,:,saver_count) = depths.d_t;
        if flags.flag_human_instability == 1
            % Maps.Hydro.risk(:,:,t_store) = Human_Instability.risk_t;
            Maps.Hydro.risk(:,:,saver_count) = Human_Instability.risk_t;
        elseif flags.flag_human_instability == 2
        elseif flags.flag_human_instability == 3
            Maps.Hydro.risk_cm(:,:,saver_count) = Human_Instability.risk_t_cm;
            Maps.Hydro.risk_tm(:,:,saver_count) = Human_Instability.risk_t_tm;
            Maps.Hydro.risk_am(:,:,saver_count) = Human_Instability.risk_t_am;
            Maps.Hydro.risk_om(:,:,saver_count) = Human_Instability.risk_t_om;
            Maps.Hydro.risk_cf(:,:,saver_count) = Human_Instability.risk_t_cf;
            Maps.Hydro.risk_tf(:,:,saver_count) = Human_Instability.risk_t_tf;
            Maps.Hydro.risk_af(:,:,saver_count) = Human_Instability.risk_t_af;
            Maps.Hydro.risk_of(:,:,saver_count) = Human_Instability.risk_t_of;
        end
        % Maps.Hydro.I_t(:,:,t_store) = Soil_Properties.I_t;
        Maps.Hydro.I_t(:,:,saver_count) = Soil_Properties.I_t;
        if flags.flag_waterquality ~= 1
            saver_count = saver_count+1;
            if saver_count > 12
                saver_count = 1;
                save(strcat('Temporary_Files\save_map_hydro_',num2str(store),'.mat'),'Maps');
                store = store + 1;
            end
        end
    end % Calls the sub
end

% Hydrographs and Pollutographs with same time resolution
recording_parameters.actual_record_hydrograph = find(running_control.time_record_hydrograph <= t_save,1,'last');
recording_parameters.delta_record_hydrograph = recording_parameters.actual_record_hydrograph - recording_parameters.last_record_hydrograph;
recording_parameters.last_record_hydrograph = recording_parameters.actual_record_hydrograph;

if flags.flag_automatic_calibration ~= 1
    if k == 1
        outlet_states.outlet_hydrograph(1,1) = nansum(nansum(outlet_states.outlet_flow))/1000*Wshed_Properties.cell_area/3600; % m3/s
        running_control.time_hydrograph(1,1) = running_control.time_calculation_routing(k,1)/60;
        outlet_states.depth_outlet(1,1) = mean(depths.d_t(idx_outlet));
        % Maximum Flodded Areas
        Flooded_Area = max(sum(sum((depths.d_t > 150)))*Wshed_Properties.Resolution^2,Flooded_Area); % Areas larger than 15 cm
        % Maximum Risk Areas
        if flags.flag_human_instability == 1
            % Maximum Risk Areas
            Risk_Area = max(sum(sum((Human_Instability.risk_t == 1)))*Wshed_Properties.Resolution^2,Risk_Area);
        end
        if flags.flag_waterquality(1,1) == 1
            Maps.WQ_States.Pol_Conc_Map(:,:,1) = WQ_States.P_conc;
            Maps.WQ_States.Pol_mass_map(:,:,1) = WQ_States.B_t;
            WQ_States.EMC_outlet(1,1) = WQ_States.mass_outlet/WQ_States.vol_outlet; % mg/L
            WQ_States.mass_outlet_save(1,1) = WQ_States.mass_outlet;
            WQ_States.vol_outlet_save(1,1) = WQ_States.vol_outlet;
        end
        % Saving Data of Input Gauges
        if flags.flag_obs_gauges == 1
            for i = 1:gauges.num_obs_gauges
                gauges.x_cell = gauges.easting_obs_gauges(i); gauges.y_cell = gauges.northing_obs_gauges(i);
                gauges.wse_cell(1,i) = depths.d_t(gauges.y_cell,gauges.x_cell)/1000 + Elevation_Properties.elevation_cell(gauges.y_cell,gauges.x_cell); % m
                gauges.depth_cell(1,i) = depths.d_t(gauges.y_cell,gauges.x_cell)/1000; % m
                gauges.hydrograph_cell(1,i) = CA_States.I_tot_end_cell(gauges.y_cell,gauges.x_cell)/(running_control.time_calculation_routing(k)); % m3/s
            end
        end
        if flags.flag_waterquality == 1
            WQ_States.WQ_States.outet_pollutograph(1,1) = Out_Conc; % Already averaged for all outlet cells
        end
    elseif recording_parameters.delta_record_hydrograph > 0
        t_store = recording_parameters.actual_record_hydrograph;
        outlet_states.outlet_hydrograph(t_store,1) = nansum(nansum(outlet_states.outlet_flow))/1000*Wshed_Properties.cell_area/3600; % m3/s
        running_control.time_hydrograph(t_store,1) = t;
        outlet_states.depth_outlet(t_store,1) = mean(depths.d_t(idx_outlet));
        % Maximum Flodded Areas
        Flooded_Area = max(sum(sum((depths.d_t > 150)))*Wshed_Properties.Resolution^2,Flooded_Area); % Areas larger than 15 cm
        if flags.flag_human_instability == 1
            % Maximum Risk Areas
            Risk_Area = max(sum(sum((Human_Instability.risk_t == 1)))*Wshed_Properties.Resolution^2,Risk_Area);
        end
        if flags.flag_waterquality == 1
            % Maps.WQ_States.Pol_Conc_Map(:,:,t_store) = WQ_States.P_conc;
            % Maps.WQ_States.Pol_mass_map(:,:,t_store) = WQ_States.B_t;
            Maps.WQ_States.Pol_Conc_Map(:,:,saver_count) = WQ_States.P_conc;
            Maps.WQ_States.Pol_mass_map(:,:,saver_count) = WQ_States.B_t;
            saver_count = saver_count+1;
            if saver_count > 12
                saver_count = 1;
                save(strcat('Temporary_Files\save_map_hydro_',num2str(store),'.mat'),'Maps');
                store = store + 1;
            end

            WQ_States.EMC_outlet(t_store,1) = WQ_States.mass_outlet/WQ_States.vol_outlet; % mg/L
            WQ_States.mass_outlet_save(t_store,1) = WQ_States.mass_outlet;
            WQ_States.vol_outlet_save(t_store,1) = WQ_States.vol_outlet;
        end
        if flags.flag_obs_gauges == 1
            % Saving Data of Input Gauges
            for i = 1:gauges.num_obs_gauges
                gauges.x_cell = gauges.easting_obs_gauges(i); gauges.y_cell = gauges.northing_obs_gauges(i);
                gauges.wse_cell(t_store,i)= depths.d_t(gauges.y_cell,gauges.x_cell)/1000 + Elevation_Properties.elevation_cell(gauges.y_cell,gauges.x_cell); % m
                gauges.depth_cell(t_store,i) = depths.d_t(gauges.y_cell,gauges.x_cell)/1000; % m
                gauges.hydrograph_cell(t_store,i) = CA_States.I_tot_end_cell(gauges.y_cell,gauges.x_cell)/(running_control.time_calculation_routing(k)); % m3/s
            end
        end

        if flags.flag_waterquality == 1
            WQ_States.outet_pollutograph(t_store,1) = Out_Conc;
            % Saving Data of Input Gauges
            if flags.flag_obs_gauges == 1
                for i = 1:gauges.num_obs_gauges
                    gauges.x_cell = gauges.northing_obs_gauges(i); gauges.y_cell = gauges.northing_obs_gauges(i);
                    WQ_States.pollutograph_cell(t_store,i)= WQ_States.P_conc(gauges.y_cell,gauges.x_cell); % mg/L
                end
            end
        end
    end % Calls the sub
end


% Saving Maximum Depths and Outlet Flows
if flags.flag_automatic_calibration ~= 1
    if k == 1
        depths.dmax_final = depths.d_t;
        velocities.vmax_final = velocities.velocity_raster;
        if flags.flag_waterquality == 1
            WQ_States.max_Conc = zeros(size(Elevation_Properties.elevation_cell));
        end
    else
        depths.dmax_final = max(depths.d_t,depths.dmax_final);
        velocities.vmax_final = max(velocities.velocity_raster,velocities.vmax_final);
        if flags.flag_waterquality == 1
            WQ_States.max_Conc = max(WQ_States.max_Conc,WQ_States.P_conc);
        end
    end
end

outlet_runoff_volume = nansum(nansum(outlet_states.outlet_flow))*time_step/60*Wshed_Properties.cell_area/Wshed_Properties.drainage_area + outlet_runoff_volume; % mm
