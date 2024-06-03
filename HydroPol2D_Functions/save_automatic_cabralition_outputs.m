% Saving Automatic Calibration Outputs
% Goal: Saving outputs from automatic calibration
% Developer: Marcus Nobrega, Ph.D.    

if flags.flag_automatic_calibration == 1
        z1_save = find(time_obs <= t,1,'last');
        if ~isempty(z1_save)
            if z1_save > z2_save
                Qmod(z1_save,1) = nansum(nansum(outlet_states.outlet_flow))/1000*Wshed_Properties.cell_area/3600; % m3/s
                if flags.flag_waterquality == 1
                    Cmod(z1_save,1) = Out_Conc; % mg/L
                else
                    Cmod(z1_save,1) = 0; % mg/L
                end
                % Observed Gauges
                for jj = 1:length(northing_obs_cell)
                    Qmod_gauges(z1_save,jj) = CA_States.I_tot_end_cell(northing_obs_cell(jj),easting_obs_cell(jj))/(running_control.time_calculation_routing(k));
                    if flags.flag_waterquality == 1
                        Cmod_gauges(z1_save,jj) = Maps.WQ_States.Pol_Conc(northing_obs_cell(jj),easting_obs_cell(jj));
                    else
                        Cmod_gauges(z1_save,jj) = 0;
                    end
                end
            end
            z2_save = z1_save;
        else
            z1_save = 0;
            z2_save = z1_save;
        end
    end