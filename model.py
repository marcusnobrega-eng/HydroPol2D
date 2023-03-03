###################################################################
#                                                                 #
#                        Produced by:                             #
#                 Marcus Nobrega Junior gomes                     #
#              (marcusnobrega.engcivil@gmail.com)                 #
#                               &                                 #
#                 Luis Miguel Castillo Rápalo                     #
#                   (luis.castillo@unah.hn)                       #
#                    September 2021                               #
#                                                                 #
#       Last update : 7 July, 2021                                #
###################################################################

# --- Import Python Libraries --- #
def model(path):
    import numpy as np
    import pandas as pd
    import math
    import time
    import matplotlib.pylab as plt
    from matplotlib import rc
    import matplotlib as mpl
    import os
    import warnings
    import glob
    import csv
    from osgeo import gdal
    # --- Import created libraries/functions --- #
    import Raster_exporter
    import resize_matrix
    import CA_Routing
    import build_up_wash_off
    import input_data
    import accumulated_incremental
    import Raster_exporter
    os.chdir(path)
    # --------------- Initial Data - General -------------------- #
    general_data = input_data.GeneralData('general_data.ini')
    lulc_data = input_data.LULCData('LULC_data.ini')
    lulc_parameters = lulc_data.variable_to_matrix()
    lulc_index = lulc_data.lulc_index_array()
    rainfall = input_data.PrecipitationData('Rainfall_Intensity_Data.csv')
    inflow = input_data.InflowData('Inflow_Hydrograph.csv', general_data.resolution)

    # -- loading DEM and LULC data into in two matrix's
    dem, imp = np.loadtxt('Elevation_Data.txt'), np.loadtxt('Land_Cover_Data.txt')
    # -- stating the matrix size
    cell_area = np.power(general_data.resolution, 2)
    ny_max, nx_max = dem.shape[0], dem.shape[1]

    # ---------------- Inflow Cells -------------------------
    inflow_cells = np.zeros((ny_max, nx_max))
    for i in range(inflow.n_inlets):
        inflow_cells[int(inflow.inlet_cells.iloc[i, 1])-1, int(inflow.inlet_cells.iloc[i, 0])-1] = 1  # from dataframe to integer
        #we subtract 1 to the inlets coordinates, Python starts with zero.

    # -- Rectangle at the beginning of the Catchment --
    # ---------------- Rainfall Matrix ----------------------
    # flag_rainfall = 0; if 1, all watershed receive rainfall, else 0
    if general_data.flag_rainfall == 0:
        rainfall_matrix = np.zeros((ny_max, nx_max))
    else:
        rainfall_matrix = np.ones((ny_max, nx_max))

    # ------------------ Outlet ---------------------------
    outlet_index_fulldomain = np.zeros((ny_max, nx_max))
    # the following for-loop is used when we enter the boundaries of the domain as the outlet
    # for i = range(int(general_data.x_outlet_begin), int(general_data.x_outlet_end)):
    #       for j = range(int(general_data.y_outlet_begin), int(general_data.y_outlet_end)):
    #           outlet_index_fulldomain[j,i]=1

    # %%%% Finding Outlets with Filled DEMs %%%%
    dem = np.where(dem < 0, np.NaN, dem)  # replacing -9999 with NaN in dem matrix
    row_outlet, col_outlet = np.zeros((1)), np.zeros((1))
    idx_outlet = list(zip(np.where(dem == np.nanmin(dem))[0], np.where(dem == np.nanmin(dem))[1]))
    row_outlet[0], col_outlet[0] = int(idx_outlet[0][0]), int(idx_outlet[0][1])  # outlet coordinates
    outlet_index = outlet_index_fulldomain  # unnecessary i think, two matrixes with the same data !!??

    # --- Checking if it is in the boundary of the domain (Old way) -- #
    # for i in range(len(idx_outlet)):
    #       if (col_outlet[i] == 1) or (col_outlet[i] == np.shape(dem)[1]) or (row_outlet[i] == 1) or (row_outlet[i] == np.shape(dem)[0]):

    # -- New way: Assuming that the outlets should be at the border of the domain, we can do:
    # for i in range(len([row_outlet])):
    #     # Checking left, right, up, and down
    #     row, col = coors[i,0], coors[i, 1]  # Row and Col  for multiple outlets
    #     # Left
    #     if (row - 1 == 0) or math.isnan(dem[row - 1, col]):
    #         outlet_index[coors[i, 0], coors[i, 1]] = 1  # This is an outlet
    #     # Right
    #     if (row + 1 > np.shape(dem)[0]) or math.isnan(dem[row + 1, col]):
    #         outlet_index[coors[i, 0], coors[i, 1]] = 1  # This is an outlet
    #     # Up
    #     if (col - 1 == 0) or math.isnan(dem[row, col - 1]):
    #         outlet_index[coors[i, 0], coors[i, 1]] = 1  # This is an outlet
    #     # Down
    #     if (col + 1 == np.shape(dem)[1]) or math.isnan(dem[row, col + 1]):  # why is not >??
    #         outlet_index[coors[i, 0], coors[i, 1]] = 1  # This is an outlet
    # row_outlet, col_outlet = np.where(outlet_index == 1)[0], np.where(outlet_index == 1)[1]  # refreshing outlet
    #
    # # -- Extending the outlet domain by a matrix of nXn size
    # extent = 2  # this can be changed according the extent from the main oultlet
    # u_row_outlet, d_row_outlet = row_outlet - extent, row_outlet + extent
    # l_col_outlet, r_col_outlet = col_outlet - extent, col_outlet + extent

    # for i in range(int(l_col_outlet), int(r_col_outlet)):
    #     for j in range(int(u_row_outlet), int(d_row_outlet)):
    #         if math.isnan(dem[j, i]):
    #             continue
    #         if (math.isnan(dem[j-1, i-1])) or (math.isnan(dem[j-1, i])) or (math.isnan(dem[j-1, i+1])) or (math.isnan(dem[j, i-1])) or (math.isnan(dem[j,i+1])) or (math.isnan(dem[j+1,i-1])) or (math.isnan(dem[j+1,i])) or (math.isnan(dem[j+1,i+1])):
    #             outlet_index[j,i] = 1
    # row_outlet_2, col_outlet_2 = np.where(outlet_index > 0)[0], np.where(outlet_index > 0)[1] #alocating the rows and columns of the "possible" boundary of the outlet.

    # %%%% Recording time of outputs (i.e. flows, concentrations....%%%%#
    # -- Calculations
    steps = general_data.routing_time / general_data.time_step_model  # Number of calculation steps
    number_of_records = math.floor(steps * general_data.time_step_model / general_data.record_time_maps)  # Number of stored data (size of the vector)
    # -- Checking if the recording time is working
    if number_of_records == 0:
        raise Exception('The recording time is larger than the routing time, please change it')

    time_records = np.arange(0, general_data.routing_time + general_data.record_time_maps, general_data.record_time_maps)  # time in minutes
    time_record_hydrograph = np.arange(0, general_data.routing_time + general_data.record_time_hydrographs, general_data.record_time_hydrographs)  # time in minutes
    time_change_records = np.arange(0, general_data.routing_time + ((general_data.time_step_change / 60)*(general_data.time_step_change / 60)), general_data.time_step_change / 60)  # time in minutes
    time_change_matrices = np.arange(0, general_data.routing_time, general_data.time_step_matrices / 60)  # time in minutes
    # -- vector to store data
    time_store = time_records / general_data.time_step_model  # number of steps necessary to reach the record vector
    time_store[0] = 1  # the zero is the firt time step

    # %%% Inflow and Precipitation data %%%#
    inflow.inflow_hydrograph_rate = inflow.inflow_hydrograph_rate * general_data.flag_inflow
    rainfall.intensity_rainfall = rainfall.intensity_rainfall * general_data.flag_rainfall  # If flag_rainfall is zero, no rainfall is considered

    # %%% Grid domain %%%#
    # --------------- Original Grid -------------------
    # In case flag_abstraction == 1, we cut domain assuming the cells entered in the input_data file
    if general_data.flag_abstraction == 0:
        xmin = 0  # initial position x in the grid (columns)
        ymin = 0  # lines (up to down)
        xmax = nx_max
        ymax = ny_max

    # -------------- Cutting Cells -----------------
    inflow_cells = inflow_cells[ymin:ymax, xmin:xmax]  # check indices due to that python don't consider the last row/column
    outlet_index = outlet_index[ymin:ymax, xmin:xmax]
    elevation = dem[ymin:ymax, xmin:xmax]  # Using only specified grid
    spatial_domain = np.zeros((ymax, xmax))
    dem = dem[ymin:ymax, xmin:xmax]
    lulc_matrix = imp[ymin:ymax, xmin:xmax]  # Using only specified grid
    rainfall_matrix = np.double(elevation >= 0)
    rainfall_matrix = np.where(rainfall_matrix == 0, math.inf, rainfall_matrix)

    # %%% Determining soil properties for each cell, according to the imperviousness and Cutting all full-domain matrices %%%#
    # --- testing if Imp and DEM have the same -9999 values
    idx_elevation = np.where(elevation < 0, True, False)  # logical
    idx_LULC = np.where(lulc_matrix < 0, True, False)  # logical
    # ----------- Checking if Data are Equivalent -------------

    # %%% Pre-allocating Arrays and Checking no-data Values
    # ------------ Old-way ------------: DEM and IMP must have the same no-data matching each other. In this case, the following if-statement must be true:
    # if max(max((idx_elevation - idx_LULC))) ~= 0
    #     error('Please, check your input data. Extent of DEM and LULC are to be the same.')
    # end
    # idx = idx_elevation;
    # idx_nan = idx_elevation; %  saving nan_matrix
    # ------------ New-way ------------:  Most of the time these data have different resolutions and no-data values. In this case, if at least one of the bothr rasters (i.e., DEM and IMP have a missing data, we assume both don't have)
    idx = (idx_elevation + idx_LULC)
    idx = idx > 0
    idx_nan = idx.copy()
    idx_nan_5 = np.zeros((5, ymax, xmax))  # 5D array of zeros / preallocate

    for i in range(1, 5):
        idx_nan_5[i, :, :] = idx_nan  # this is used in CA model / here could be upgrade for 9 "moore" model
    idx_nan_5 = idx_nan_5 > 0  # turning the arrays into logicals
    elevation[idx] = math.inf  # change those values for inf
    lulc_matrix[idx] = math.inf  # change those values for inf
    spatial_domain[idx] = math.inf  # change those values for inf

    # %%% Preallocate arrays to avoid large computational efforts
    y_1 = math.dist((ymin,), (ymax,))
    x_1 = math.dist((xmin,), (xmax,))
    roughness, d_0, h_0, psi = spatial_domain.copy(), spatial_domain.copy(), spatial_domain.copy(), spatial_domain.copy()
    ksat, teta_i, teta_sat, I_0 = spatial_domain.copy(), spatial_domain.copy(), spatial_domain.copy(), spatial_domain.copy()  # not necessary for now, they are prelocated again
    if general_data.flag_waterquality == 1: #Build-up and Wash-off matrices
        C_1, C_2, C_3, C_4, B_0 = spatial_domain.copy(), spatial_domain.copy(), spatial_domain.copy(), spatial_domain.copy(), spatial_domain.copy()
    # if general_data.flag_initial_buildup == 1: # Only if initial buildup is being modeled to prealcoate variables
    #     B_0 = spatial_domain.copy()

    # ----------------- Fill the variables -----------------
    # --- Correcting LULC_index ---
    lulc_matrix[idx_nan] = -9999
    idx_1 = np.where(lulc_matrix == lulc_index[0], 1, 0)  # find imp1 surfaces
    idx_2 = np.where(lulc_matrix == lulc_index[1], 1, 0)  # find per2 surfaces
    idx_3 = np.where(lulc_matrix == lulc_index[2], 1, 0)  # find per3 surfaces
    idx_4 = np.where(lulc_matrix == lulc_index[3], 1, 0)  # find per4 surfaces
    idx_5 = np.where(lulc_matrix == lulc_index[4], 1, 0)  # find per5 surfaces
    idx_6 = np.where(lulc_matrix == lulc_index[5], 1, 0)  # find per6 surfaces
    idx = np.stack([idx_1, idx_2, idx_3, idx_4, idx_5, idx_6])
    impervious_cells = np.sum(idx_1)
    pervious_cells = np.sum(idx[1:5, :, :])

    for i in range(len(lulc_parameters[:, 0])): #number of lulc classes, it could be upgrade for n-classes
        roughness = np.where(idx[i, :, :] == True, lulc_parameters[i, 0], roughness)
        ksat = np.where(idx[i, :, :] == True, lulc_parameters[i, 1], ksat)
        h_0 = np.where(idx[i, :, :] == True, lulc_parameters[i, 3], h_0)
        psi = np.where(idx[i, :, :] == True, lulc_parameters[i, 4], psi)
        teta_i = np.where(idx[i, :, :] == True, lulc_parameters[i, 7], teta_i)
        teta_sat = np.where(idx[i, :, :] == True, lulc_parameters[i, 6], teta_sat)
        I_0 = np.where(idx[i, :, :] == True, lulc_parameters[i, 5], I_0)
        # ------------ Warm-up Data ------------
        if general_data.flag_warmup == 1:
            d_0 = np.loadtxt('Warmup_Depths.txt')*1000 # Depths in meter converted to mm
            # Treat inf values
            d_0[idx_nan] = math.inf #carefull here!
        else:
            d_0 = np.where(idx[i, :, :] == True, lulc_parameters[i, 2], d_0)
        # --------- Water-Quality Data ---------
        if  general_data.flag_waterquality == 1: # Only if water quality is being modeled
            C_1 = np.where(idx[i, :, :] == True, lulc_parameters[i, 8], C_1)
            C_2 = np.where(idx[i, :, :] == True, lulc_parameters[i, 9], C_2)
            C_3 = np.where(idx[i, :, :] == True, lulc_parameters[i, 10], C_3)
            C_4 = np.where(idx[i, :, :] == True, lulc_parameters[i, 11], C_4)
            if i + 1 == len(lulc_parameters[:, 0]):
                # Initial Build-up Calculation
                if general_data.flag_initial_buildup == 1:
                    B_0 = np.loadtxt('Initial_Buildup.txt')  # Kg/cell
                    # Treat inf values
                    B_0[idx_nan] = math.inf  # carefull here
                    B_t = B_0;
                else:  # Let's calculate it
                    B_0 = np.array(C_1) * np.array(1 - np.exp(-1 * C_2 * general_data.ADD))  # Kg/ha, look variables
                    B_t = B_0 * cell_area / math.pow(10, 4)  # Kg

    # d_0 = warm_up_depth //// to delete it
    # clear wipe space in memory to free memory
    del idx, idx_1, idx_2, idx_3, idx_4, idx_5, idx_6

    # --------- Full Domain Cells -----------------
    rainfall_matrix_full_domain, teta_sat_fulldomain, teta_i_fulldomain, psi_fulldomain = rainfall_matrix.copy(), teta_sat.copy(), teta_i.copy(), psi.copy()
    roughness_fulldomain, ksat_fulldomain, h_0_fulldomain, I_p_fulldomain, I_t_fulldomain= roughness.copy(), ksat.copy(), h_0.copy(), I_0.copy(), I_0.copy()
    outlet_index_fulldomain = outlet_index_fulldomain[ymin:ymax, xmin:xmax]
    # ---- New xmax and ymax ----
    xmin, ymin = 1, 1
    zzz = np.shape(elevation)
    xmax, ymax = zzz[1], zzz[0]

    #%%%% Conversion of inflow into the time-step of calculations %%%%
    inflow_length = len(inflow.inflow_hydrograph_rate)-1 #Number of intervals
    inflow_discretized = np.zeros((1, math.ceil(inflow_length*inflow.time_step_inflow/general_data.time_step_model)))  # preallocating
    for i in np.arange((inflow_length-1)*inflow.time_step_inflow/general_data.time_step_model):
        inflow_discretized[0,int(i)] = inflow.inflow_hydrograph_rate[math.ceil((round(i*general_data.time_step_model/inflow.time_step_inflow,12)))]

    #%%% Conversion of rainfall into the time-step of calculation %%%
    intensity_rainfall_lenght = len(rainfall.intensity_rainfall)-1 #Number of intervals
    intensity_discretized = np.zeros((1, math.ceil(intensity_rainfall_lenght*rainfall.time_step_rainfall/general_data.time_step_model))) #preallocating
    for i in np.arange(intensity_rainfall_lenght*rainfall.time_step_rainfall/general_data.time_step_model):
        intensity_discretized[0, int(i)] = rainfall.intensity_rainfall[math.ceil((round((i + 1) * general_data.time_step_model / rainfall.time_step_rainfall, 12))) - 1]  # Discretized // check round problem

    # %%% Determination of grid parameters and outlet coordinates (Whole Domain)
    nx_max = len(elevation[0, :])  # number of x cells
    ny_max = len(elevation[:, 0])  # number of y cells
    # --- calculates the number of non-inf cells to determine the watershed area ---
    number_inf = np.sum(np.isinf(elevation))
    drainage_area = (nx_max * ny_max - number_inf) * cell_area
    impervious_area = np.sum(impervious_cells) * cell_area
    pervious_area = np.sum(pervious_cells) * cell_area
    impervious_rate = impervious_area / drainage_area
    pervious_rate = pervious_area / drainage_area
    # --- Pollutant Mass ---
    if general_data.flag_waterquality == 1:
        initial_mass = np.ma.masked_invalid(B_t).sum()  # Kg of pollutant / check sum ingnoring inf cells / not tested

    #%%% Calculation of accumulated and incremental inflow on the inflow cells %%%
    delta_inflow, inflow_intensity = accumulated_incremental.accumulated_incremental(steps, inflow_discretized, general_data.time_step_model)
    time_deltainflow = np.cumsum(np.ones((len(delta_inflow)))*general_data.time_step_model)
    #%%% Calculation of accumulated and incremental precipitation in each cells %%%
    delta_p, rainfall_intensity = accumulated_incremental.accumulated_incremental(steps, intensity_discretized, general_data.time_step_model)
    time_deltap = np.cumsum(np.ones((len(delta_p)))*general_data.time_step_model)
    outlet_runoff_volume = np.array([])

    #--- Pre allocating more arrays ---#
    Total_Inflow = 0
    d_t = spatial_domain.copy()
    risk_t = spatial_domain.copy()
    risk_t_bti = spatial_domain.copy()
    risk_t_fti = spatial_domain.copy()
    outflow_cms_t = spatial_domain.copy()
    inf_m3s_t = spatial_domain.copy()
    q_exut_t = spatial_domain.copy()
    #--- Space dependent arrays ---#
    qout_left_t = spatial_domain.copy()
    qout_right_t = spatial_domain.copy()
    qout_up_t = spatial_domain.copy()
    qout_down_t = spatial_domain.copy()
    d_p = spatial_domain.copy()
    hydraulic_radius_t = spatial_domain.copy()
    d_avg = spatial_domain.copy()
    qin_t = spatial_domain.copy()
    vel_left = spatial_domain.copy()
    vel_right = spatial_domain.copy()
    vel_up = spatial_domain.copy()
    vel_down = spatial_domain.copy()
    I_tot_end = np.ones(np.shape(spatial_domain)) #ATENTION HERE!
    I_tot_end_fulldomain = spatial_domain.copy()
    wse_outlet_calculated = 0
    #inundated_cels = spatial_domain

    #---------- Time dependent arrays ---------------
    time_size = len(time_store)
    d = np.zeros((time_size, ny_max, nx_max))
    risk = np.zeros((time_size, ny_max, nx_max))
    risk_bti = np.zeros((time_size, ny_max, nx_max))
    risk_fti = np.zeros((time_size, ny_max, nx_max))
    wse_outlet = np.zeros((time_size, 1)) #number of identified outlets
    wse_outlet_2 = np.zeros((time_size, 1)) #for a especific outlet
    outlet_hydrograph = np.zeros((time_size, 1))
    outlet_pollutograph = np.zeros((time_size, 1))
    time_hydrograph = np.zeros((time_size, 1))
    time_step_save_1 = np.array([])
    time_step_save_2 = np.array([])
    flow_velocity = np.zeros((int(steps), 1))
    flow_acceleration = np.zeros((int(steps), 1))
    alfa_save = np.zeros((int(steps), 1))
    alfa_min = np.asarray(0.05)
    slope_alfa = np.asarray(0.225)  # 25% of alfa
    v_threshold = np.asarray(3)  # % m/s
    if general_data.flag_waterquality == 1:
        Pol_Conc_Map = np.zeros((time_size, ny_max, nx_max))
        Pol_Mass_Map = np.zeros((time_size, ny_max, nx_max))
        outet_pollutograph = np.zeros((time_size, 1))

    # --- Clearing a few variables
    del intensity_discretized, inflow_discretized

    # %%% Main Routing %%%#
    tic = time.time()  # Star counting time
    # ------------ Find Initial Domain ------------- / to exclude
    yx_depth_stored = np.argwhere(d_0 > 0)
    rowcol_check = np.argwhere(inflow_cells > 0)
    row_col_check_rainfall = np.argwhere(rainfall_matrix > 0)
    cells_check = np.unique(np.vstack((yx_depth_stored, rowcol_check, row_col_check_rainfall)), axis=0)
    # flows_cells = gpuArray(zeros(ny_max + 2, nx_max + 2))

    # ---------- Initialize Variables --------------
    time_step_factor = general_data.alfa
    flows_cells = np.zeros((ny_max, nx_max))  # Attention here
    I_t, I_p = I_0.copy(), I_0.copy()
    time_calculation_routing = np.zeros((1, 1))  # Undefined due to the adaptively time-step
    k, actual_record_state, last_record_maps, last_record_hydrograph = 0, 0, 0, 0
    delta_inflow_agg, delta_p_agg = delta_inflow[0], delta_p[0] * general_data.flag_rainfall
    time_step = general_data.time_step_model
    t = general_data.time_step_model  # in minutes
    I = np.zeros((4, 1))  # Volumes in each direction
    time_save_previous, Previous_Volume, t_previous, coordinate_x, coordinate_y, time_matrices_previous = 0, 0, 0, 1, 1, 0
    # delta_h_matrix = np.zeros((1,5))
    step_error, time_step_cell, relative_vol_error = np.zeros((1, int(steps))), general_data.max_time_step * np.ones(
        (1, int(steps))), np.zeros((1, int(steps)))
    I_cell = np.zeros((5, ny_max, nx_max))
    slope = np.zeros((4, ny_max, nx_max))
    C, P_conc, mass_lost, Tot_Washed = 0, 0, 0, 0
    tmin_wq = general_data.max_time_step
    dmax_final = np.zeros((np.shape(elevation)))
    risk_final = np.zeros((np.shape(elevation)))
    risk_bti_final = np.zeros((np.shape(elevation)))
    risk_fti_final = np.zeros((np.shape(elevation)))

    # %%%% ------ Main while ------ %%%%#

    while t <= (general_data.routing_time + general_data.min_time_step / 60):
        # ----------- Infiltration and Available Depth -------------
        if general_data.flag_waterquality == 1:
            print("Perc(%) = " + str(np.round_((t / general_data.routing_time) * 100, 2)) + " || t(sec) = " + str(
                np.round(time_step * 60, 2)) + " || dt(mm) = " + str(
                np.round(np.nanmax(d_t), 2)) + " || inf(mm/hr) = " + str(
                np.round(np.nanmax(C), 2)) + "|| C(mg/l) = " + str(np.round_(np.max(P_conc), 2)) + " || dtmWQ = " + str(
                np.round(tmin_wq, 2)))
        else:
            print("Perc(%) = " + str(np.round_((t / general_data.routing_time) * 100, 2)) + " || t(sec) = " + str(
                np.round(time_step * 60, 2)) + " || dt(mm) = " + str(
                np.round(np.nanmax(d_t), 2)) + " || inf(mm/hr) = " + str(
                np.round(np.nanmax(C), 2)))
        if tmin_wq < 0 or np.isnan(tmin_wq) or np.isinf(tmin_wq):
            ttt = 1
        if k == 0:
            if general_data.flag_infiltration == 1:
                C = np.array(ksat) * np.array(
                    1 + np.array(np.array(d_0 + psi) * np.array(teta_sat - teta_i - 0.001)) / np.array(I_0))  # Matrix form. The 0.001 avoids teta_s - teta_i = 0
                Inflow_Rate = np.array(delta_p_agg * rainfall_matrix + delta_inflow_agg * inflow_cells + d_0) / np.array(
                    time_step / 60)
                tx = np.minimum(C, Inflow_Rate)
                I_t = I_0 + tx * time_step / 60
                pef = delta_p_agg * rainfall_matrix + delta_inflow_agg * inflow_cells - tx * time_step / 60
                d_t = d_0 + pef
                inf_m3s_t = np.array(tx / 1000) * (cell_area / 3600)
            else:
                Inflow_Rate = np.array(delta_p_agg * rainfall_matrix + delta_inflow_agg * inflow_cells + d_0) / np.array(
                    time_step / 60)
                pef = delta_p_agg * rainfall_matrix + delta_inflow_agg * inflow_cells
                d_t = d_0 + pef
        else:
            if general_data.flag_infiltration == 1:
                ## Effective precipitation - Green-Ampt (1911)
                C = np.array(ksat) * np.array(
                    1 + np.array(np.array(d_p + psi) * np.array(teta_sat - teta_i - 0.001)) / np.array(I_p))  # Matrix form. The 0.001 avoids 0/0 values
                Inflow_Rate = np.array(delta_p_agg * rainfall_matrix + delta_inflow_agg * inflow_cells + d_p) / np.array(
                    time_step / 60)
                tx = np.minimum(C, Inflow_Rate)
                I_t = I_p + tx * time_step / 60
                pef = delta_p_agg * rainfall_matrix + delta_inflow_agg * inflow_cells - tx * time_step / 60
                d_t = d_p + pef  # ATTENTION HERE
                inf_m3s_t = np.array(tx / 1000) * (cell_area / 3600)
            else:
                Inflow_Rate = np.array(delta_p_agg * rainfall_matrix + delta_inflow_agg * inflow_cells + d_p) / np.array(
                    time_step / 60)
                pef = delta_p_agg * rainfall_matrix + delta_inflow_agg * inflow_cells
                d_t = d_p + pef
        total_available_depth = d_t

        #  ---- ELEVATIONS ----  #
        elevation_left_t = np.hstack([np.zeros((ny_max, 1)), elevation[:, 0:(nx_max-1)]])
        elevation_right_t = np.hstack([elevation[:, 1:nx_max], np.zeros((ny_max, 1))])
        elevation_up_t = np.vstack([np.zeros((1, nx_max)), elevation[0:(ny_max-1), :]])
        elevation_down_t = np.vstack([elevation[1: ny_max], np.zeros((1, nx_max))])
        #  ---- WATER DEPTHS ----  #
        d_left_cell = np.hstack([np.zeros((ny_max, 1)), total_available_depth[:, 0:(nx_max-1)]])
        d_right_cell = np.hstack([total_available_depth[:, 1:nx_max], np.zeros((ny_max, 1))])
        d_up_cell = np.vstack([np.zeros((1, nx_max)), total_available_depth[0:(ny_max-1), :]])
        d_down_cell = np.vstack([total_available_depth[1: ny_max], np.zeros((1, nx_max))])
        #  ---- ASSIGNING VALUES ----  #
        elevation_cell = elevation.copy()
        d_t_cell = total_available_depth.copy()
        hydraulic_radius_t_cell = hydraulic_radius_t.copy()
        roughness_cell = roughness.copy()
        h_0_cell = h_0.copy()
        I_tot_end_cell = I_tot_end.copy()  # ATTENTION HERE!

        # %%% flood routing for each cell %%% #
        if np.max(np.where(np.isinf(d_t) == True, 1, 0)) == 1:  # check if Nan or Inf occurred in d_t
            ttt = 1
        flag_ca = 0
        if flag_ca == 0:  # Different CA's model
            qout_left_t, qout_right_t, qout_up_t, qout_down_t, outlet_flow, d_t, I_tot_end, I_cell = CA_Routing.CA_Routing(
                elevation_cell, elevation_left_t, elevation_right_t, elevation_up_t, elevation_down_t, d_t_cell,
                d_left_cell, d_right_cell, d_up_cell, d_down_cell, roughness_cell, cell_area, time_step, h_0_cell,
                general_data.resolution, I_tot_end_cell, outlet_index, general_data.outlet_type, general_data.slope_outlet,
                row_outlet, col_outlet, idx_nan, general_data.flag_critical)

        else:
            qout_left, qout_right, qout_up, qout_down, outlet_flow, d_t_cell, I_tot_end_cell, I = CA_Routing.CA_Routing_2(
                elevation_cell, elevation_left_t, elevation_right_t, elevation_up_t, elevation_down_t, d_t_cell,
                d_left_cell, d_right_cell, d_up_cell, d_down_cell, roughness_cell, cell_area, time_step, h_0_cell,
                general_data.resolution, I_tot_end_cell, outlet_index, general_data.outlet_type, general_data.slope_outlet,
                row_outlet, col_outlet, idx_nan, general_data.flag_critical)

        # --- Calculating Qin --- #
        qin_left_t = np.hstack([np.zeros((ny_max, 1)), qout_right_t[:, 0:(nx_max-1)]])
        qin_right_t = np.hstack([qout_left_t[:, 1:nx_max], np.zeros((ny_max, 1))])
        qin_up_t = np.vstack([np.zeros((1, nx_max)), qout_down_t[0:(ny_max-1), :]])
        qin_down_t = np.vstack([qout_up_t[1: ny_max], np.zeros((1, nx_max))])
        qin_t = qin_left_t + qin_right_t + qin_up_t + qin_down_t
        idx = qin_t > general_data.flow_tolerance
        idx2 = qin_t <= general_data.flow_tolerance
        flows_cells[idx] = 1
        flows_cells[idx2] = 0
        qin_t[np.isinf(qin_t)] = 0  # Adding a constraint at q_in such that it does not change d_t for inf
        qin_t[idx2] = 0  # Adding a constraint at q_in such that it does not change d_t for nans

        if general_data.flag_waterquality == 1:
            B_t, P_conc, Out_Conc, tmin_wq, tot_W_out, mass_lost, Tot_Washed = build_up_wash_off.build_up_wash_off(C_3, C_4, qout_left_t,
                                                                                            qout_right_t, qout_up_t,
                                                                                            qout_down_t, outlet_flow, B_t,
                                                                                            time_step, nx_max, ny_max,
                                                                                            cell_area, outlet_index,
                                                                                            idx_nan_5,
                                                                                            general_data.flag_wq_model,
                                                                                            mass_lost,
                                                                                            Tot_Washed)

        # %%% Checking Mass Balance %%% #
        if general_data.flag_waterquality == 1:
            if np.ma.masked_invalid(B_t).sum() > initial_mass:
                ttt = 1
        # %%% New Time-step Calculation %%% #
        pos_save = np.argwhere(time_change_records < t)[-1]
        time_save = time_change_records[pos_save]  # min
        delta_time_save = time_save - time_save_previous
        time_save_previous = time_save
        actual_record_timestep = np.argwhere(time_change_records < t)[-1]
        if delta_time_save > 0 or k == 0:  # First time-step
            if general_data.flag_timestep == 0:
                # %%% SOLUTION FOR COURANT METHOD %%% #
                d_left_cell = np.hstack([np.zeros((ny_max, 1)), d_t[:, 0:(nx_max - 1)]])
                d_right_cell = np.hstack([d_t[:, 1:nx_max], np.zeros((ny_max, 1))])
                d_up_cell = np.vstack([np.zeros((1, nx_max)), d_t[0:(ny_max - 1), :]])
                d_down_cell = np.vstack([d_t[1: ny_max], np.zeros((1, nx_max))])
                # Considering Velocity as Celerity
                vel_left = (np.array(I_cell[0, :, :]) / np.array(
                    (0.5 / 1000) * np.array(d_t + d_left_cell) * general_data.resolution)) / (time_step * 60) + np.sqrt(
                    9.81 * (d_t + d_left_cell) / 2 / 1000)  # m/s
                vel_right = (np.array(I_cell[1, :, :]) / np.array(
                    (0.5 / 1000) * np.array(d_t + d_right_cell) * general_data.resolution)) / (time_step * 60) + np.sqrt(
                    9.81 * (d_t + d_right_cell) / 2 / 1000)
                vel_up = (np.array(I_cell[2, :, :]) / np.array(
                    (0.5 / 1000) * np.array(d_t + d_up_cell) * general_data.resolution)) / (time_step * 60) + np.sqrt(
                    9.81 * (d_t + d_up_cell) / 2 / 1000)
                vel_down = (np.array(I_cell[3, :, :]) / np.array(
                    (0.5 / 1000) * np.array(d_t + d_down_cell) * general_data.resolution)) / (time_step * 60) + np.sqrt(
                    9.81 * (d_t + d_down_cell) / 2 / 1000)

                # ----------- Find the Maximum Velocity -------------- #
                max_velocity_left = np.nanmax(vel_left)
                max_velocity_right = np.nanmax(vel_right)
                max_velocity_up = np.nanmax(vel_up)
                max_velocity_down = np.nanmax(vel_down)

                # --- Maximum of all of them --- #
                max_velocity = np.max([max_velocity_left, max_velocity_right, max_velocity_up, max_velocity_down])

                # --- Velocity Raster --- #
                velocity_raster = np.maximum(vel_left, vel_right, vel_up)
                velocity_raster = np.maximum(velocity_raster, vel_down)

                # --- Calculation of varying Alpha --- #
                alfa_max = time_step_factor
                if k == 0:
                    dvdt = 0
                    max_velocity_previous = max_velocity
                else:
                    dvdt = np.maximum((max_velocity-max_velocity_previous)/(general_data.time_step_change), 0)  #m/s2
                flow_velocity[pos_save, 0] = max_velocity  # m/s
                flow_acceleration[pos_save, 0] = time_step  # m/s2
                #  if flag_D8 ==1:
                #     factor_grid = cp.sqrt(1/2)
                # else:
                #     factor_grid = 1
                factor_grid =1
                new_timestep = (factor_grid*general_data.resolution / max_velocity)  # seconds
                time_step_factor = np.maximum(alfa_max - slope_alfa*(np.maximum(max_velocity - v_threshold, 0)), alfa_min)
                alfa_save[pos_save, 0] = time_step_factor

                # --- Calculation of Stability --- #
                mu = 0.5  # Friction coefficient
                Cd = 1.1
                ro_person = 1000  # kg/m3
                weight_person = 75  # kg
                height_person = 1.75  # meters
                width1_person = 0.3  # meters
                width2_person = 0.3  # meters
                ro_water = 1000  # kg/m3
                gravity = 9.81  # m/s2
                F_person = weight_person*gravity  # N
                # Buyoance
                F_buoy = (width1_person*width2_person)*d_t/1000*ro_water*gravity  # N
                # Available Friction
                available_friction = mu*np.maximum(F_person - F_buoy, 0)
                # Hydrodynamic Force
                hydro_force = np.maximum((0.5*(Cd*width1_person*d_t /1000 * np.power(velocity_raster,2))), 0)
                # Risk Factor
                risk_t = np.minimum(hydro_force/available_friction,1)
                risk_t[idx_nan] = -9999


                # --- Calculation of BTI and FTI Stability --- #
                vel_t_zzz = velocity_raster.copy()
                d_t_zzz = d_t.copy()
                d_t_zzz[np.isnan(d_t_zzz)] = -9999
                vel_t_zzz[np.isnan(vel_t_zzz)] = -9999

                risk_t_bti = np.logical_or(np.logical_or(np.logical_and(np.logical_and((0.37 < vel_t_zzz), (vel_t_zzz <= 3)), (d_t_zzz/1000 >= 1.63*np.exp(-1.37*vel_t_zzz) + 0.22)),
                                            np.logical_and((vel_t_zzz < 0.37), (d_t_zzz/ 1000 >= 1.2))), vel_t_zzz >= 3).astype(int)

                risk_t_fti = np.logical_or(np.logical_or(np.logical_and(np.logical_and((0.45 < vel_t_zzz), (vel_t_zzz <= 3)), (d_t_zzz/1000 >= 1.43*np.exp(-0.88*vel_t_zzz) + 0.24)),
                                            np.logical_and((vel_t_zzz < 0.45), (d_t_zzz/ 1000 >= 1.2))), vel_t_zzz >= 3).astype(int)

                risk_t_bti[idx_nan] = -9999
                risk_t_fti[idx_nan] = -9999  # no data

                # ------- Adding water Quality minim Time-step ------- #
                # new_timestep = np.minimum(0.9 * tmin_wq, time_step_factor * new_timestep)
                t_previous = t
                previous_time_step = time_step
            else:  # Solution for the Stable Method - Time-step update
                # Water slopes Calculation (THERE IS A MISTAKE HERE)
                wse_cell = elevation_cell + d_t_cell / 1000
                slope[0, :, :] = (1 / general_data.resolution) * (wse_cell - elevation_left_t - d_left_cell / 1000)
                slope[1, :, :] = (1 / general_data.resolution) * (wse_cell - elevation_right_t - d_right_cell / 1000)
                slope[2, :, :] = (1 / general_data.resolution) * (wse_cell - elevation_up_t - d_up_cell / 1000)
                slope[3, :, :] = (1 / general_data.resolution) * (wse_cell - elevation_down_t - d_down_cell / 1000)
                slope_cell = np.maximum(np.nanmin(slope, axis=0), general_data.slope_min)  # Minimum assumed slope
                time_step_cell = np.maximum(np.power(general_data.resolution, 2) / 4 * (
                            2 * np.array(roughness) / np.array(np.power(np.array(d_t_cell / 1000), 5 / 3)) * np.array(
                        np.sqrt(slope_cell))), general_data.min_time_step)
                time_step_cell[idx_nan_5[0, :, :]] = math.inf
                new_timestep = np.nanmin(time_step_cell)
                # Adding water Quality Minimum Time-step
                new_timestep = np.minimum(tmin_wq, new_timestep)
                t_previous = t
                previous_time_step = time_step
            # Rounding time-step to min_time-step or max_time-step with the required precision
            if general_data.flag_wq_model==1:
                # --- Adding Water Quality factor of water quality
                alfa_wq = 1  # Redcution factor of water quality
                new_timestep = np.minimum(alfa_wq*tmin_wq, new_timestep)
            else:
                new_timestep = time_step_factor*new_timestep
            if k == 0:
                time_calculation_routing[k] = new_timestep
            else:
                time_calculation_routing = np.vstack((time_calculation_routing, new_timestep))
            time_calculation_routing[k] = np.maximum(general_data.time_step_increments * np.floor(time_calculation_routing[k] / general_data.time_step_increments), general_data.min_time_step)
            time_calculation_routing[k] = np.minimum(time_calculation_routing[k], general_data.max_time_step)
            time_step = time_calculation_routing[k] / 60  # New time-step for the next run
            if time_calculation_routing[k] == general_data.min_time_step:  # Check if we reached the minimum time-step
                # print("the model has got instability. Please, decrease the minimum time-step or decrease the duration where the new time-step is calculated.")
                # break
                unstable = 1
        elif k == 0:
            time_calculation_routing[k] = np.vstack((time_calculation_routing, time_step * 60))  # sec
        else:
            time_calculation_routing = np.vstack((time_calculation_routing, time_calculation_routing[k - 1]))

        # -----  Aggregating Inflows to the New Time-step ----- #
        if general_data.flag_inflow > 0:
            z1 = int(np.argwhere(time_deltainflow > t_previous)[0])  # begin of the time-step
            z2 = int(np.argwhere(time_deltainflow <= t)[-1])  # end of the time-step
            if z2 < z1:
                z2 = z1
            if time_step >= general_data.time_step_model:
                delta_inflow_agg = np.mean(delta_inflow[z1:z2]) / (general_data.time_step_model * 60) * time_step * 60
            else:
                delta_inflow_agg = delta_inflow[z1] / (general_data.time_step_model * 60) * time_step * 60
            if k == 0:
                t_previous = time_calculation_routing[k] / 60
            else:
                t_previous = t
            if np.isnan(delta_inflow_agg):
                ttt = 1

        # -----  Aggregating Precipitation to the New Time-step ----- #
        if general_data.flag_rainfall > 0:
            z1 = int(np.argwhere(time_deltap > t_previous)[0])  # begin of the time-step
            z2 = int(np.argwhere(time_deltap <= t)[-1])  # end of the time-step
            if z2 < z1:
                z2 = z1
            if time_step >= general_data.time_step_model:
                if z1 == z2:  # Mean function does not work if z1 = z2
                    delta_p_agg = np.mean(delta_p[z1]) / (general_data.time_step_model * 60) * time_step * 60
                else:
                    delta_p_agg = np.mean(delta_p[z1:z2]) / (general_data.time_step_model * 60) * time_step * 60
            else:
                delta_p_agg = delta_p[z1] / (general_data.time_step_model * 60) * time_step * 60
            if k == 0:
                t_previous = time_calculation_routing[k] / 60
            else:
                t_previous = t

        # %%% New Domain - Matrices %%% #
        flag_deactivate = 1
        if (general_data.flag_rainfall == 0) and (general_data.flag_inflow == 1) and (
                flag_deactivate == 0):  # We are simulating only inflow hydrograph
            pos_matrices = int(np.argwhere(time_change_matrices < t)[-1])
            time_matrices = time_change_matrices[pos_matrices]  # min
            if k == 0:
                time_matrices = 1
            delta_time_matrices = time_matrices - time_matrices_previous
            time_matrices_previous = time_matrices

            if delta_time_matrices > 0:
                # --- Finding cells receiving water --- #
                # --- Local Coordinate System from the beginning of the Domain --- #
                yx_receiving = np.argwhere(flows_cells > 0)
                # --- Finding cells with stored water --- #
                yx_depth_stored = np.argwhere(d_t > general_data.depth_tolerance)
                row_col_check = np.argwhere(inflow_cells > 0)
                cells_check = np.unique(np.vstack((yx_receiving, yx_depth_stored, row_col_check)),
                                        axis=0)  # Domain Fronteirs - Global Coordinate System

                xmin_cells = np.min(cells_check[:, 1] - 1) + coordinate_x
                xmin_cells_t = xmin_cells  # Update
                xmax_cells = np.max(cells_check[:, 1] - 1) + coordinate_x
                ymin_cells = np.min(cells_check[:, 0] - 1) + coordinate_y
                ymin_cells_t = ymin_cells  # Update
                ymax_cells = np.max(cells_check[:, 0] - 1) + coordinate_y

                # %%% Min Values %%% # (only works in the absolute coordinate reference system
                xmin_domain_new = np.maximum(xmin_cells - general_data.factor_cells, 1)
                ymin_domain_new = np.maximum(ymin_cells - general_data.factor_cells, 1)

                # %%% Max Values %%% #
                xmax_domain_new = np.minimum(xmax_cells + general_data.factor_cells, xmax)  # xmax is wrong here
                ymax_domain_new = np.minimum(ymax_cells + general_data.factor_cells, ymax)
                domain = np.zeros((int(ymax_domain_new - ymin_domain_new + 1), int(xmax_domain_new - xmin_domain_new + 1)))
                # %%% Domain is in the local coordinate system %%% #
                # --- New Coordinates - Local of the 1st cells in the Global Coordinate System --- #
                coordinate_x_previous = coordinate_x
                coordinate_y_previous = coordinate_y
                if coordinate_x == 1:  # First cut in the domain (outside border)
                    coordinate_x = xmin_domain_new  # Using the outside reference (new reference)
                    coordinate_y = ymin_domain_new
                    begin = 1  # Initial Abstraction of the Matrices
                else:
                    begin = 0
                    coordinate_x = xmin_domain_new
                    coordinate_y = ymin_domain_new

                # Matrices Sizes - Deslocation (local Coordinate system) #
                dymax_matrix = ymax_cells - coordinate_y + 1  # Top of the matrix
                dymin_matrix = ymin_cells - coordinate_y + 1  # Bottom
                dxmax_matrix = xmax_cells - coordinate_x + 1  # Right side
                dxmin_matrix = xmin_cells - coordinate_x + 1  # Left side
                ny_max, nx_max = np.shape(domain)[0], np.shape(domain)[1]

                # --- Abstracted Matrices --- #
                if begin == 1:
                    ymin_domain_abstracted = ymin_cells - 0  # Already in the right coordinate system
                    ymax_domain_abstracted = ymax_cells - 0
                    xmin_domain_abstracted = xmin_cells - 0
                    xmax_domain_abstracted = xmax_cells - 0
                else:
                    ymin_domain_abstracted = ymin_cells - coordinate_y_previous + 1
                    ymax_domain_abstracted = ymax_cells - coordinate_y_previous + 1
                    xmin_domain_abstracted = xmin_cells - coordinate_x_previous + 1
                    xmax_domain_abstracted = xmax_cells - coordinate_x_previous + 1

                # New Matrices
                # --- Depths --- #
                d_t = resize_matrix.resize_matrix(d_t, domain, ymin_domain_abstracted, ymax_domain_abstracted,
                                                  xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix,
                                                  dymax_matrix, dxmin_matrix, dxmax_matrix, 1)
                # --- inflow Cells --- #
                inflow_cells = resize_matrix.resize_matrix(inflow_cells, domain, ymin_domain_abstracted,
                                                           ymax_domain_abstracted, xmin_domain_abstracted,
                                                           xmax_domain_abstracted, dymin_matrix, dymax_matrix, dxmin_matrix,
                                                           dxmax_matrix, 1)
                # --- Rainfall Matrix --- #
                rainfall_matrix = resize_matrix.resize_matrix(rainfall_matrix, domain, ymin_domain_abstracted,
                                                              ymax_domain_abstracted, xmin_domain_abstracted,
                                                              xmax_domain_abstracted, dymin_matrix, dymax_matrix,
                                                              dxmin_matrix, dxmax_matrix, 1)
                # --- Flow cells --- #
                flows_cells = resize_matrix.resize_matrix(flows_cells, domain, ymin_domain_abstracted,
                                                          ymax_domain_abstracted, xmin_domain_abstracted,
                                                          xmax_domain_abstracted, dymin_matrix, dymax_matrix, dxmin_matrix,
                                                          dxmax_matrix, 1)
                # --- Inundation cells --- #
                # inundated_cells = resize_matrix.resize_matrix(inundated_cells, domain, ymin_domain_abstracted, ymax_domain_abstracted, xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix, dymax_matrix, dxmin_matrix, dxmax_matrix, 1)
                # --- Qins --- #
                qin_t = resize_matrix.resize_matrix(qin_t, domain, ymin_domain_abstracted, ymax_domain_abstracted,
                                                    xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix,
                                                    dymax_matrix, dxmin_matrix, dxmax_matrix, 1)
                # qout_left_t = resize_matrix(qout_left_t, domain, ymin_domain_abstracted, ymax_domain_abstracted, xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix, dymax_matrix, dxmin_matrix, dxmax_matrix, 1)
                # qout_right_t = resize_matrix(qout_right_t, domain, ymin_domain_abstracted, ymax_domain_abstracted, xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix, dymax_matrix, dxmin_matrix, dxmax_matrix, 1)
                # qout_up_t = resize_matrix(qout_up_t, domain, ymin_domain_abstracted, ymax_domain_abstracted, xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix, dymax_matrix, dxmin_matrix, dxmax_matrix, 1)
                # qout_down_t = resize_matrix(qout_down_t, domain, ymin_domain_abstracted, ymax_domain_abstracted, xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix, dymax_matrix, dxmin_matrix, dxmax_matrix, 1)
                # --- Velocities --- #
                vel_left = resize_matrix.resize_matrix(vel_left, domain, ymin_domain_abstracted, ymax_domain_abstracted,
                                                       xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix,
                                                       dymax_matrix, dxmin_matrix, dxmax_matrix, 1)
                vel_rigt = resize_matrix.resize_matrix(vel_right, domain, ymin_domain_abstracted, ymax_domain_abstracted,
                                                       xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix,
                                                       dymax_matrix, dxmin_matrix, dxmax_matrix, 1)
                vel_up = resize_matrix.resize_matrix(vel_up, domain, ymin_domain_abstracted, ymax_domain_abstracted,
                                                     xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix,
                                                     dymax_matrix, dxmin_matrix, dxmax_matrix, 1)
                vel_down = resize_matrix.resize_matrix(vel_down, domain, ymin_domain_abstracted, ymax_domain_abstracted,
                                                       xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix,
                                                       dymax_matrix, dxmin_matrix, dxmax_matrix, 1)
                # --- I_tot_end_cell --- #
                I_tot_end = resize_matrix.resize_matrix(I_tot_end, domain, ymin_domain_abstracted, ymax_domain_abstracted,
                                                        xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix,
                                                        dymax_matrix, dxmin_matrix, dxmax_matrix, 1)
                I_tot_end_cell = I_tot_end
                # --- I_cell --- #
                I_cell = resize_matrix.resize_matrix(I_cell, domain, ymin_domain_abstracted, ymax_domain_abstracted,
                                                     xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix,
                                                     dymax_matrix, dxmin_matrix, dxmax_matrix, 1)

                # %%% All Domain Cells %%% #
                # --- Elevation --- #
                elevation = dem[int(ymin_domain_new):int(ymax_domain_new), int(xmin_domain_new):int(xmax_domain_new)]
                # --- Infiltration Parameters --- #
                teta_sat = teta_sat_fulldomain[int(ymin_domain_new):int(ymax_domain_new),
                           int(xmin_domain_new):int(xmax_domain_new)]
                teta_i = teta_i_fulldomain[int(ymin_domain_new):int(ymax_domain_new),
                         int(xmin_domain_new):int(xmax_domain_new)]
                psi = psi_fulldomain[int(ymin_domain_new):int(ymax_domain_new), int(xmin_domain_new):int(xmax_domain_new)]
                ksat = ksat_fulldomain[int(ymin_domain_new):int(ymax_domain_new), int(xmin_domain_new):int(xmax_domain_new)]
                # --- Roughness --- #
                roughness = roughness_fulldomain[int(ymin_domain_new):int(ymax_domain_new),
                            int(xmin_domain_new):int(xmax_domain_new)]
                # --- Initial Abstration --- #
                h_0 = h_0_fulldomain[int(ymin_domain_new):int(ymax_domain_new), int(xmin_domain_new):int(xmax_domain_new)]
                # --- Accumulated Infiltration --- #
                domain_infiltration = I_0[int(ymin_domain_new):int(ymax_domain_new),
                                      int(xmin_domain_new):int(xmax_domain_new)]
                I_t = resize_matrix.resize_matrix(I_t, domain, ymin_domain_abstracted, ymax_domain_abstracted,
                                                  xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix,
                                                  dymax_matrix, dxmin_matrix, dxmax_matrix, 1)
                I_p = I_t
                # --- Outlet Cells --- #
                outlet_index = outlet_index_fulldomain[int(ymin_domain_new):int(ymax_domain_new),
                               int(xmin_domain_new):int(xmax_domain_new)]
                # --- Rainfall cells --- #
                rainfall_matrix = rainfall_matrix_full_domain[int(ymin_domain_new):int(ymax_domain_new),
                                  int(xmin_domain_new):int(xmax_domain_new)]

                # --- Cells Check --- #
                yx_receiving = np.argwhere(flows_cells > 0)
                yx_depth_stored = np.argwhere(d_t > general_data.depth_tolerance)
                row_col_check = np.argwhere(inflow_cells > 0)
                cells_check = np.unique(np.vstack((yx_receiving, yx_depth_stored, row_col_check)),
                                        axis=0)  # Domain Fronteirs - Global Coordinate System

        # %%% Inflows and Depth Update %%% #
        d_t = d_t + qin_t * time_step / 60
        # Inf Values
        d_t[np.isinf(d_t)] = np.nan  # Replacing inf to nan
        # --- Clearing stored values --- #
        # qin_t = np.zeros(ny_max,nx_max)
        qout_left_t = np.zeros((ny_max, nx_max))
        qout_right_t = np.zeros((ny_max, nx_max))
        qout_up_t = np.zeros((ny_max, nx_max))
        qout_down_t = np.zeros((ny_max, nx_max))

        # --------- Checking Water Balance and Saving variables --------- #
        if general_data.flag_waterbalance == 1:
            Delta_Total_inflow = delta_inflow_agg / 1000 * inflow.n_inlets * cell_area + delta_p_agg / 1000 * (
                        ny_max * nx_max) * cell_area - np.sum(outlet_flow) / 1000 * cell_area * time_step / 60  # m3
            if k == 0:  # add initial depths
                Total_Inflow = Delta_Total_inflow + np.sum(d_0 / 1000 * cell_area)  # Cumulative in m3
            else:
                Total_Inflow = Total_Inflow + Delta_Total_inflow  # Cumulative in m3
            # --- stored Water --- #
            Stored_Volume = np.sum(np.array(cell_area / 1000) * np.array(d_t)) + np.sum(
                (I_t - I_p) / 1000 * cell_area)  # m3 of water depth + infiltration
            Delta_Stored = Stored_Volume - Previous_Volume  # m3
            # jump the previous volume and calculate at the end #
            step_error[k] = Delta_Stored / Delta_Total_inflow
            if step_error[k] < 0.9:
                tttt = 1
            # --- Water Balance Error --- #
            relative_vol_error[k] = (Stored_Volume - Total_Inflow) / Total_Inflow * 100
            extra_depth = (Delta_Total_inflow - Delta_Stored) / (inflow.n_inlets * cell_area) * 1000
            d_t = d_t + inflow_cells * extra_depth
            Stored_Volume = np.sum(cell_area / 1000 * d_t)  # m3 (d_t already updated)
            Delta_Stored = Stored_Volume - Previous_Volume  # m3
            # --- updating Previous Volume --- #
            Previous_Volume = Stored_Volume  # m3
            step_error[k] = Delta_Stored / Delta_Total_inflow

            # --- Water Balance Error --- #
            relative_vol_error[k] = (Stored_Volume - Total_Inflow) / Total_Inflow * 100

        # %%% Saving Plotting Values - Recording Time %%% #
        # Maps of Flood Depth, WSE and Pollutant Concentrations #
        t_save = t + time_calculation_routing[k] / 60
        actual_record_state = np.argwhere(time_records < t_save)[-1]
        delta_record = actual_record_state - last_record_maps
        last_record_maps = actual_record_state

        if k == 0:
            d[0, :, :] = d_t
            risk[0, :, :] = risk_t
            risk_bti[0, :, :] = risk_t_bti
            risk_fti[0, :, :] = risk_t_fti
            # Water surface elevation outlet
            acumulator = []  # np.zeros()
            wse_outlet[0, :] = d_t[idx_outlet[0][0], idx_outlet[0][1]] / 1000 + elevation[
                idx_outlet[0][0], idx_outlet[0][1]]  # extracting the water surface elevation
            wse_outlet_2[0, :] = d_t[36, 22] / 1000 + elevation[36, 22]  # Modify this!!!!!!
            if general_data.flag_waterquality == 1:
                Pol_Conc_Map[0, :, :] = P_conc
                Pol_Mass_Map[0, :, :] = B_t
        elif delta_record > 0:
            t_store = actual_record_state
            d[t_store, :, :] = d_t
            risk[t_store, :, :] = risk_t
            risk_bti[t_store, :, :] = risk_t_bti
            risk_fti[t_store, :, :] = risk_t_fti
            # Water surface elevation outlet #
            wse_outlet[t_store, :] = d_t[idx_outlet[0][0], idx_outlet[0][1]] / 1000 + elevation[
                idx_outlet[0][0], idx_outlet[0][1]]
            wse_outlet_2[t_store, :] = d_t[36, 22] / 1000 + elevation[36, 22]
            if general_data.flag_waterquality == 1:
                Pol_Conc_Map[t_store, :, :] = P_conc
                Pol_Mass_Map[t_store, :, :] = B_t
        # Calls the sub
        # --- Hydrographs --- #
        actual_record_hydrograph = np.argwhere(time_record_hydrograph < t_save)[-1]
        delta_record_hydrograph = actual_record_hydrograph - last_record_hydrograph
        last_record_hydrograph = actual_record_hydrograph

        if k == 0:
            outlet_hydrograph[0, 0] = np.nansum(outlet_flow) / 1000 * cell_area / 3600  # m^3/s
            time_hydrograph[0, 0] = time_calculation_routing[k] / 60
            if general_data.flag_waterquality == 1:
                outlet_pollutograph[0, 0] = Out_Conc
        elif delta_record_hydrograph > 0:
            t_store = actual_record_hydrograph
            outlet_hydrograph[t_store, 0] = np.nansum(outlet_flow) / 1000 * cell_area / 3600  # m3/s
            time_hydrograph[t_store, 0] = t
            if general_data.flag_waterquality == 1:
                outlet_pollutograph[t_store, 0] = Out_Conc
            # --- Waitbar (t/routing_time)
        # Calls the sub
        # Saving Maximum Depths and Outlet flows
        dmax_final = np.maximum(d_t, dmax_final)
        risk_final = np.maximum(risk_t, risk_final)
        risk_bti_final = np.maximum(risk_t_bti, risk_bti_final)
        risk_fti_final = np.maximum(risk_t_fti, risk_fti_final)


        outlet_runoff_volume = np.append(outlet_runoff_volume,
                                         np.nansum(outlet_flow) * time_step / 60 * cell_area / drainage_area)  # mm
        # Previous Depths
        d_p = d_t
        I_p = I_t

        # Check if nan or Inf occured in d_t
        if np.max(np.where(np.isinf(d_t) == True, 1, 0)) == 1:
            ttt = 1
        # increase the counter
        t = time_calculation_routing[k] / 60 + t
        time_step_save_1 = np.append(time_step_save_1, time_calculation_routing[k])
        time_step_save_2 = np.append(time_step_save_2, t)
        k = k + 1

    # %%% Post-Processing Results %%% #
    toc = time.time()
    # --- Outlet Hydrograph --- #
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)
    mpl.rcParams.update(mpl.rcParamsDefault)
    # latex font is missing #""""""
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax2 = ax1.twinx()
    fig_1 = ax1.plot(time_hydrograph, outlet_hydrograph, 'black', linewidth=1.5, label='Outflow')
    ax1.set_xlabel('Time $(min)$')
    ax1.set_ylabel('Flow Discharge $(m^3/s)$')
    fig_2 = ax2.plot(rainfall.time_rainfall, rainfall.intensity_rainfall, 'blue', linewidth=1.5, label='Rainfall')
    ax2.set_ylim([0, 200])
    ax2.invert_yaxis()
    ax2.set_ylabel('Rainfall intensity $(mm/hr)$')
    if general_data.flag_inflow == 1:  # If inflow is included this will be in the plot
        fig_3 = ax1.plot(inflow.time_inflow, inflow.inflow_discharge, 'black', linestyle='dashed', linewidth=1.5,
                         label='Inflow')
        fig.legend(prop={'size': 10}, loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)
    else:
        fig.legend(prop={'size': 10}, loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)
    # plt.show()
    plt.savefig('Hydrograph.pdf', dpi=300)
    # plt.close()

    Outlet_Hydrograph_Data = np.hstack((time_hydrograph, outlet_hydrograph))
    Outlet_Hydrolevels_Data = np.hstack((time_hydrograph, wse_outlet_2))
    alfa_vel_acel_data = np.hstack((alfa_save, flow_velocity, flow_acceleration))
    with open('Outlet_hydrograph_Data.txt', 'w', encoding='UTF8') as f:  # Exporting Discharge data to .txt file
        writer = csv.writer(f, delimiter=',')
        writer.writerow(['Time (min)', 'Discharge (m3/s)'])
        writer.writerows(Outlet_Hydrograph_Data)
    with open('Outlet_hydrolevels_Data.txt', 'w', newline='', encoding='UTF8') as f:  # Exporting leves data to .txt file
        writer = csv.writer(f, delimiter=',')
        writer.writerow(['Time (min)','levels (m)'])
        writer.writerows(Outlet_Hydrolevels_Data)
    with open('alpha_max_vel_time_step.txt', 'w', newline='', encoding='UTF8') as f:  # Exporting Discharge data to .txt file
        writer = csv.writer(f, delimiter=',')
        writer.writerow(['alpha(-)', 'max_velocity (m/s)', 'time_step (s)'])
        writer.writerows(alfa_vel_acel_data)
    if general_data.flag_waterquality == 1:  # Plotting and exporting Pollutograph
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        fig_1 = ax1.plot(time_hydrograph, outlet_pollutograph, 'red', linewidth=1.5)
        ax1.set_xlabel('Time $(min)$')
        ax1.set_ylabel('Pol. Concentration $(mg/L)$')
        fig.legend(prop={'size': 10}, loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=ax1.transAxes)
        plt.show()
        plt.savefig('Pollutograph.pdf', dpi=300)
        plt.close()
        Outlet_Pollutograph_Data = np.hstack((time_hydrograph, outlet_pollutograph))
        with open('Outlet_Pollutograph_Data.txt', 'w', encoding='UTF8') as f:  # Exporting Discharge data to .txt file
            writer = csv.writer(f, delimiter=',')
            writer.writerow(['Time (min)', 'Pol. Concentration (mg/L)'])
            writer.writerows(Outlet_Pollutograph_Data)

    # --- Exporting Raters --- #
    # Delete  previous raster in the folder
    no_data_value = -9999
    # Specify the path where the files are.
    # path = os.getcwd()
    # Check to make sure that folder actually exists. Warn user if it doesn't.
    if os.path.isdir(path) == False:
        warnings.warn('Error: The main path does not exist')
    # Get a list of all files in the folder with the desired file name pattern.
    filePattern = [f for f in glob.glob("*.asc")]
    for i in range(len(filePattern)):  # Deleting previous maps in ASCII format
        os.remove(str(path) + "/" + str(filePattern[i]))
    filePattern = [f for f in glob.glob("*.gif")]
    for i in range(len(filePattern)):  # Deleting previous .gifs
        os.remove(str(path) + "/" + str(filePattern[i]))

    # Changing Nan Values
    d[np.isnan(d)] = 0
    d[np.isinf(d)] = 0
    # Raster - Depths, WSE, and Pol. Conc
    for i in range(len(time_records)):
        idx = d[i, :, :] < general_data.depth_wse * 1000
        time_map = time_records[i]
        if general_data.flagwse == 0:
            FileName = "Flood_depths_t_" + str(int(time_map)) + "_min"
            raster_exportion = d[i, :, :] / 1000
            raster_exportion[idx] = no_data_value
        else:
            FileName = "Water_Surface_Elevation_t_" + str(int(time_map)) + "_min"
            raster_exportion = d[i, :, :] / 1000
            raster_exportion[idx] = no_data_value
            raster_exportion = pd.DataFrame(raster_exportion + dem)
        # Risk rasters #
        FileName_risk = 'Risk_t_' + str(int(time_map)) + "_min"
        FileName_risk_bti = 'Risk_BTI_t_' + str(int(time_map)) + "_min"
        FileName_risk_fti = 'Risk_FTI_t_' + str(int(time_map)) + "_min"
        raster_exportion_risk = risk[i, :, :]
        raster_exportion_risk_bti = risk_bti[i, :, :]
        raster_exportion_risk_fti = risk_fti[i, :, :]
        # Save Map in ASCII format and TIFF
        Raster_exporter.raster_exporter(FileName, np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                        general_data.yllupcorner, general_data.resolution, no_data_value, raster_exportion)
        Raster_exporter.raster_exporter(FileName_risk_bti, np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                        general_data.yllupcorner, general_data.resolution, no_data_value,
                                        raster_exportion_risk_bti)
        Raster_exporter.raster_exporter(FileName_risk_fti, np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                        general_data.yllupcorner, general_data.resolution, no_data_value,
                                        raster_exportion_risk_fti)
        # Water Quality
        if general_data.flag_waterquality == 1:
            FileName = "Pollutant_Concentration" + str(int(time_map)) + "_min"
            raster_exportion = Pol_Conc_Map[i, :, :] / 1000
            idx = raster_exportion < general_data.pol_min  # Finding values below Pol_min
            raster_exportion[idx] = no_data_value
            # Save Map in ASCII format and TIFF
            Raster_exporter.raster_exporter(FileName, np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                            general_data.yllupcorner, general_data.resolution, no_data_value,
                                            raster_exportion)

    # --- Point of accumulation --- #
    if general_data.flag_waterquality == 1:
        FileName = 'Accumulation_Areas_1m'
        depth_accumulation = 1
        idx_depth = d_t > depth_accumulation * 1000
        raster_exportion = no_data_value * np.ones((np.shape(d_t)))
        raster_exportion[idx_depth] = 1
        Raster_exporter.raster_exporter(FileName, np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                        general_data.yllupcorner, general_data.resolution, no_data_value, raster_exportion)
    # --- Point of accumulation --- #
    if general_data.flag_waterquality == 1:
        FileName = "Accumulation_Areas"
        pol_accumulation = 0.01  # Kg/m^2
        zzz = B_t / cell_area
        zzz[np.isinf(zzz)] = no_data_value
        idx_bt = zzz > pol_accumulation  # Kg/m^2
        raster_exportion = no_data_value * np.ones((np.shape(B_t)))
        final_mass = Pol_Mass_Map[-1, :, :] / cell_area  # map in Kg divided by area
        raster_exportion[idx_bt] = final_mass[idx_bt]
        Raster_exporter.raster_exporter(FileName, np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                        general_data.yllupcorner, general_data.resolution, no_data_value, raster_exportion)
    # --- Maximum --- #
    zzz = dmax_final / 1000  # m
    idx = zzz < general_data.depth_wse  # Finding values below the threshold
    if general_data.flagwse == 0:  # Save only max depth
        # mMaximum_depths
        FileName = "Maximum_Depths"
        zzz[np.isinf(zzz)] = no_data_value
        zzz[np.isnan(zzz)] = no_data_value
        zzz[idx] = no_data_value
        raster_exportion = zzz
        Raster_exporter.raster_exporter(FileName, np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                        general_data.yllupcorner, general_data.resolution, no_data_value, raster_exportion)
    else:
        raster_exportion = zzz + dem  # Surface water elevation
        raster_exportion[idx] = no_data_value
        FileName = "Max_Water_Surface_Elevation"
        Raster_exporter.raster_exporter(FileName, np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                        general_data.yllupcorner, general_data.resolution, no_data_value, raster_exportion)
    # Risk rasters #
    # Save Map in ASCII format and TIFF
    Raster_exporter.raster_exporter('Risk_final', np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                    general_data.yllupcorner, general_data.resolution, no_data_value,
                                    risk_final)
    Raster_exporter.raster_exporter('Risk_BTI_final', np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                    general_data.yllupcorner, general_data.resolution, no_data_value,
                                    risk_bti_final)
    Raster_exporter.raster_exporter('Risk_FTI_final', np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                    general_data.yllupcorner, general_data.resolution, no_data_value,
                                    risk_fti_final)
    # --- DEM --- #
    FileName = "DEM"
    zzz = elevation
    zzz[np.isinf(zzz)] = no_data_value
    raster_exportion = zzz
    Raster_exporter.raster_exporter(FileName, np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                    general_data.yllupcorner, general_data.resolution, no_data_value, raster_exportion)
    # --- LULC --- #
    FileName = "Land_Use_&_Land_Cover_Data"
    zzz = imp
    zzz[np.isinf(zzz)] = no_data_value
    zzz[np.isnan(zzz)] = no_data_value
    raster_exportion = zzz
    Raster_exporter.raster_exporter(FileName, np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                    general_data.yllupcorner, general_data.resolution, no_data_value, raster_exportion)

    if general_data.flag_waterquality == 1:
        for i in range(len(time_records)):
            if i == len(time_records):  # Water Quality - Mass of Pollutant
                time_map = time_records[i]
                FileName = "Mass_of_pollutant" + str(time_map) + "_min"
                raster_exportion = Pol_Mass_Map[i, :, :]
                idx = raster_exportion < 0.01  # Finding values below Pol_min
                raster_exportion[idx] = no_data_value
                raster_exportion = raster_exportion / cell_area  # Kg/m^2
                Raster_exporter.raster_exporter(FileName, np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                                general_data.yllupcorner, general_data.resolution, no_data_value,
                                                raster_exportion)
        # Maximum
        FileName = "Maximum_Pol_Conc" + str(time_map) + "_min"
        raster_exportion = np.max(Pol_Conc_Map, axis=0)
        idx = raster_exportion < general_data.pol_min  # Finding values below Pol_min
        raster_exportion[idx] = no_data_value
        Raster_exporter.raster_exporter(FileName, np.size(dem[0, :]), np.size(dem[:, 0]), general_data.xllupcorner,
                                        general_data.yllupcorner, general_data.resolution, no_data_value, raster_exportion)
        # Final Pollutant Mass

    # Generate GIF Files of the Simulation
    # Inundtion_maps !!!!!!!!!!!!!!!!!!!!

    # Show Summary Results
    print("Drainage Area = " + str(round(drainage_area / 1000 / 1000, 3)) + " Km2")
    print("Impervious Area = " + str(round(impervious_area / 1000 / 1000, 3)) + " Km2")
    print("pervious Area = " + str(round(pervious_area / 1000 / 1000, 3)) + " Km2")
    print("Simulation time = " + str(round((toc - tic) / 60, 3)) + " Minutes")
    print("Rainfall volume = " + str(round(np.sum(delta_p), 3)) + " mm")
    print("Runoff volume = " + str(round(np.sum(outlet_runoff_volume), 3)) + " mm")
    print("Runoff Coefficient = " + str(round((np.sum(outlet_runoff_volume) / np.sum(delta_p)), 3)) + " ")
    # Summary
    print("Maximum Flood Depth = " + str(round(1 / 1000 * np.max(d), 3)) + " m")
    print("Maximum Outflow = " + str(round(np.max(outlet_hydrograph), 4)) + " m^3/s")

    with open('Summary_table.txt', 'w', newline='', encoding='UTF8') as f:  # Exporting leves data to .txt file
        writer = csv.writer(f, delimiter=',')
        writer.writerow(['Variable', 'value'])
        writer.writerow(['Drainage_area (Km2)', round(60000000 / 1000 / 1000, 3)])
        writer.writerow(['Impervious_area (km2)', round(impervious_area / 1000 / 1000, 3)])
        writer.writerow(['Pervious_area (km2)', round(pervious_area / 1000 / 1000, 3)])
        writer.writerow(['Simulation_time (min)', round((toc - tic) / 60, 3)])
        writer.writerow(['Rainfall volume (mm)', round(np.sum(delta_p), 3)])
        writer.writerow(['Runoff Coefficient', round((np.sum(outlet_runoff_volume) / np.sum(delta_p)), 3)])
    if general_data.flag_waterquality == 1:
        print("Maximum Concentration = " + str(round(np.max(Pol_Conc_Map), 3)) + " mg/L")
        print("Maximum Stored Mass of Pollutant = " + str(round(np.max(B_t), 4)) + " Kg")
        print("Initial Pollutant Mass = " + str(round(initial_mass / 1000, 4)) + " ton")
        print("Final Pollutant Mass = " + str(round(np.ma.masked_invalid(B_t).sum(), 3)) + " ton")
