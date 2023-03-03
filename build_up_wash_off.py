###################################################################
#                                                                 #
#                        Produced by:                             #
#                 Marcus Nobrega Junior gomes                     #
#              (marcusnobrega.engcivil@gmail.com)                 #
#                               &                                 #
#                 Luis Miguel Castillo Rápalo                     #
#                   (luis.castillo@unah.hn)                       #
#                        September 2021                           #
#                                                                 #
#       Last update : 1 July, 2021                                #
#    Goal: Solve the Build-up and Wash-off Model                  #
###################################################################
import numpy as np
import math

def build_up_wash_off(C_3, C_4, qout_left_t, qout_right_t, qout_up_t, qout_down_t, outlet_flow, B_t, time_step, nx_max, ny_max, cell_area, outlet_index, idx_nan_5, flag_wq_model, mass_lost, Tot_Washed):
    # Pollutan Rate Flux
    # ----------- Vector of Flows ------------- #
    # 3-D array with layers 0, 1, 2, 3, 4, and 5 being left, right, up, down, and
    # outlet flow rates in mm/hr
    q_out_t = np.stack([qout_left_t, qout_right_t, qout_up_t, qout_down_t, outlet_flow])
    # Minimum Values
    min_Bt = 1  # g/m2
    idx_bt_1 = np.where((B_t * 1000 / cell_area < min_Bt) + np.isnan(B_t) + np.isinf(B_t) == True, 1,
                        0)  # Nans, Infs, or values below min_Bt
    idx_bt = np.tile(idx_bt_1, [5, 1, 1])  # Creating a 3-D array
    B_begin = B_t
    Bmin = 10  # g / m2
    Bmax = 100  # g / m2
    Bmin_kg = Bmin / 1000 * cell_area  # kg
    Bmax_kg = Bmax / 1000 * cell_area  # kg
    # ----------- Choosing Which Wash-off Equation ---------- #
    if flag_wq_model == 1:  # Rating Curve
        B_t_extra = np.where(B_t * 1000 / cell_area <= Bmin, 0, B_t)  # 100g / m2
        B_t_extra = np.where(B_t * 1000 / cell_area >= Bmax, Bmax_kg, B_t_extra)  # 100g / m2
        W_out_t = np.array(C_3) * np.array(np.power(np.array(q_out_t), np.array(C_4)))  # Kg/hr (rating curve)
        f_bt = (1 + (np.maximum(B_t_extra - Bmin_kg, 0)))  # f(B(t)) - water quality rating curve selection
        W_out_t = np.array(C_3) * np.array(
            np.power(np.array(q_out_t / 1000 / 3600 * cell_area), np.array(C_4))) * np.array(
            f_bt)  # kg/hr (rating curve)
    else:
        W_out_t = np.array(C_3) * np.array(np.power(np.array(q_out_t / 1000 / 3600), np.array(C_4))) * np.array(
            B_t)  # kg/hr (mass based curve)
        W_out_t[idx_bt] = 0  # No flux if mass is too low (below 1 g/m2)

        # ----------- Total Outflows per Cell ---------- #
    tot_W_out = np.sum(W_out_t, axis=0)  # kg/hr
    tot_q_out = np.sum(q_out_t, axis=0)  # mm/hr
    small_number = 1e-16

    # ----------- Replacing Nans to 0 ---------- #
    tot_W_out[idx_nan_5[0, :, :]] = 0  # Values with nan are replaced by 0
    tot_q_out[idx_nan_5[0, :, :]] = 0  # Values with nan are replaced by 0
    W_out_t[idx_nan_5] = 0  # Values with nan are replaced by 0

    # ----------- Assigning Wash-off Rates for Each Direction ---------- #
    W_out_left_t = W_out_t[0, :, :]
    W_out_right_t = W_out_t[1, :, :]
    W_out_up_t = W_out_t[2, :, :]
    W_out_down_t = W_out_t[3, :, :]
    W_out_outlet_t = W_out_t[4, :, :]
    # ----------- Creating Pollutant Inflow Matrices ---------- #
    W_in_left_t = np.hstack([np.zeros((ny_max, 1)), W_out_right_t[:, 0:(nx_max - 1)]])  # Left
    W_in_right_t = np.hstack([W_out_left_t[:, 1:nx_max], np.zeros((ny_max, 1))])  # Right
    W_in_up_t = np.vstack([np.zeros((1, nx_max)), W_out_down_t[0:(ny_max - 1), :]])  # Up
    W_in_down_t = np.vstack([W_out_up_t[1: ny_max], np.zeros((1, nx_max))])  # Down
    W_in_t = np.stack([W_in_left_t, W_in_right_t, W_in_up_t, W_in_down_t])  # outlet
    # ----------- Net Pollutant Rate ---------- #
    # At the same time, pollutants are entering and leaving the cells. the
    # matrix dW measures the difference between inflows and outflows
    dW = np.array(np.sum(W_in_t, axis=0)) - tot_W_out  # inflow pol - Outflow pol (kg/hr)
    B_t = np.where(np.isnan(B_t) == True, 0, B_t)
    # %%% Minimum Time-Step (adaptative time-step) %%% #
    if np.max(dW) <= 0:
        tmin_wq = 10
    elif np.min(B_t) <= 0:
        # Find values equal 0
        zzzz = B_t.copy()  # Copying the values of B_t
        zzzz[zzzz * 1000 / cell_area < min_Bt] = math.inf  # Replacing negative values to inf
        dW_zzz = dW.copy()  # Flux of new pollutant discharge (kg/h)
        dW_zzz[dW_zzz > 0] = 0  # Replacing positive values to 0
        tmin_wq = 3600 * np.nanmin(np.array(zzzz) / np.array(np.abs(dW_zzz) + small_number))  # time-step in seconds
    else:  # Maybe we should delete since both conditions are very similar if not equal
        zzzz = B_t.copy()  # Copying the values of B_t
        zzzz[zzzz * 1000 / cell_area < min_Bt] = math.inf  # Replacing negative values to inf
        dW_zzz = dW.copy()  # Flux of new pollutant discharge (kg/h)
        dW_zzz[dW_zzz > 0] = 0  # Replacing positive values to 0
        tmin_wq = 3600 * np.nanmin(np.array(zzzz) / np.array(np.abs(dW_zzz) + small_number))  # time-step in seconds

    # %%%% Overall Mass Balance at Cells %%%% #
    B_t_mid = B_t + dW * time_step / 60  # Refreshing Pollutant Mass (mass balance)
    if np.min(B_t_mid) < 0 and tmin_wq < time_step * 60:  # Break time-step internally
        n_steps_decimal = time_step * 60 / (tmin_wq)  # decimal time - step
        n_steps_floor = np.floor(n_steps_decimal)  # integer lowest value
        steps = (n_steps_floor + 1)  # number of steps within a model time-step
        outflow_mass = 0  # Starting to measure outflow pollutant mass (kg)
        outlet_mass = 0  # Starting to measure outlet outflow pollutant mass (kg)
        for i in range(int(steps)):
            if i == (steps - 1):
                # Break time step internally. Flow rates won't change but B_t will
                dt = (n_steps_decimal - n_steps_floor) * tmin_wq
            else:
                dt = time_step / (n_steps_floor)
            # ----------- Choosing Which Wash-off Equation ---------- #
        if flag_wq_model == 1:  # Rating Curve
            B_t_extra = np.where(B_t * 1000 / cell_area <= Bmin, 0, B_t)  # 100g / m2
            B_t_extra = np.where(B_t * 1000 / cell_area >= Bmax, Bmax_kg, B_t_extra)  # 100g / m2
            W_out_t = np.array(C_3) * np.array(np.power(np.array(q_out_t), np.array(C_4)))  # Kg/hr (rating curve)
            f_bt = (1 + (np.maximum(B_t_extra - Bmin_kg, 0)))  # f(B(t)) - water quality rating curve selection
            W_out_t = np.array(C_3) * np.array(
                np.power(np.array(q_out_t / 1000 / 3600 * cell_area), np.array(C_4))) * np.array(
                f_bt)  # kg/hr (rating curve)
        else:
            W_out_t = np.array(C_3) * np.array(np.power(np.array(q_out_t / 1000 / 3600), np.array(C_4))) * np.array(
                B_t)  # kg/hr (mass based curve)
            W_out_t[idx_bt] = 0  # No flux if mass is too low (below 1 g/m2)

            # ----------- Total Outflows per Cell ---------- #
        tot_W_out = np.sum(W_out_t, axis=0)  # kg/hr
        tot_q_out = np.sum(q_out_t, axis=0)  # mm/hr
        small_number = 1e-16

        # ----------- Replacing Nans to 0 ---------- #
        tot_W_out[idx_nan_5[0, :, :]] = 0  # Values with nan are replaced by 0
        tot_q_out[idx_nan_5[0, :, :]] = 0  # Values with nan are replaced by 0
        W_out_t[idx_nan_5] = 0  # Values with nan are replaced by 0

        # ----------- Assigning Wash-off Rates for Each Direction ---------- #
        W_out_left_t = W_out_t[0, :, :]
        W_out_right_t = W_out_t[1, :, :]
        W_out_up_t = W_out_t[2, :, :]
        W_out_down_t = W_out_t[3, :, :]
        W_out_outlet_t = W_out_t[4, :, :]
        # ----------- Creating Pollutant Inflow Matrices ---------- #
        W_in_left_t = np.hstack([np.zeros((ny_max, 1)), W_out_right_t[:, 0:(nx_max - 1)]])  # Left
        W_in_right_t = np.hstack([W_out_left_t[:, 1:nx_max], np.zeros((ny_max, 1))])  # Right
        W_in_up_t = np.vstack([np.zeros((1, nx_max)), W_out_down_t[0:(ny_max - 1), :]])  # Up
        W_in_down_t = np.vstack([W_out_up_t[1: ny_max], np.zeros((1, nx_max))])  # Down
        W_in_t = np.stack([W_in_left_t, W_in_right_t, W_in_up_t, W_in_down_t])  # outlet
        # ----------- Net Pollutant Rate ---------- #
        # At the same time, pollutants are entering and leaving the cells. the
        # matrix dW measures the difference between inflows and outflows
        dW = np.array(np.sum(W_in_t, axis=0)) - tot_W_out  # inflow pol - Outflow pol (kg/hr)
        outflow_mass = outflow_mass + tot_W_out * dt / 60  # kg of pollutant that left the cell
        outlet_mass = outlet_mass + W_out_outlet_t * dt / 60  # kg of pollutant that left the outlet cells
        B_t = B_t + dW * dt / 60  # Mass Balance for the incremental time-step

        # ---------- Average Pollutant Flux in the Time-Step ----------- #
        tot_W_out = outflow_mass / (time_step / 60)  # kg / hr
        # ---------- Average Pollutant Flux a the outlet in the Time-step ----------- #
        W_out_outlet_t = outlet_mass / (time_step / 60)  # kg / hr
        # ---------- Rounding Negative Values ----------- #
        mass_lost = np.sum(B_t[B_t < 0]) + mass_lost
        B_t[B_t < 0] = 0  # Replacing negative values to 0
    else:
        B_t = B_t_mid  # We can use the same time-step from the hydrodynamic model
        # ---------- Rounding Negative Values ----------- #
        mass_lost = np.sum(B_t[B_t < 0]) + mass_lost
        B_t[B_t < 0] = 0  # Replacing negative values to 0

    # %%% Minimum Time-Step (adaptative time-step) %%% #
    if np.max(dW) <= 0:
        tmin_wq = 10
    elif np.min(B_t) <= 0:
        # Find values equal 0
        zzzz = B_t.copy()  # Copying the values of B_t
        zzzz[zzzz * 1000 / cell_area < min_Bt] = math.inf  # Replacing negative values to inf
        dW_zzz = dW.copy()  # Flux of new pollutant discharge (kg/h)
        dW_zzz[dW_zzz > 0] = 0  # Replacing positive values to 0
        tmin_wq = 3600 * np.nanmin(np.array(zzzz) / np.array(np.abs(dW_zzz) + small_number))  # time-step in seconds
    else:  # Maybe we should delete since both conditions are very similar if not equal
        zzzz = B_t.copy()  # Copying the values of B_t
        zzzz[zzzz * 1000 / cell_area < min_Bt] = math.inf  # Replacing negative values to inf
        dW_zzz = dW.copy()  # Flux of new pollutant discharge (kg/h)
        dW_zzz[dW_zzz > 0] = 0  # Replacing positive values to 0
        tmin_wq = 3600 * np.nanmin(np.array(zzzz) / np.array(np.abs(dW_zzz) + small_number))  # time-step in seconds

    B_t = np.round(B_t, 6)
    # ----------- Adding a constraint in q_out to avoid huge numbers ----------- #
    tot_q_out = np.maximum(tot_q_out, 10)  # Constraint at the minimum flow
    # ----------- Pollutant Concentration ----------- #
    P_conc = np.maximum(np.power(10, 6) * np.array(tot_W_out) / np.array(tot_q_out * cell_area), 0)  # mg / L
    # ----------- Outlet Concentration ----------- #
    Out_Conc = np.maximum((1000 * np.sum(W_out_outlet_t)) / ((np.sum(outlet_flow) / 1000) * cell_area), 0)  # mg / L
    Tot_Washed = (np.max(B_begin - B_t, 0)) + Tot_Washed

    return B_t, P_conc, Out_Conc, tmin_wq, tot_W_out, mass_lost, Tot_Washed


def build_up_wash_off_GPU(C_3, C_4, qout_left_t, qout_right_t, qout_up_t, qout_down_t, outlet_flow, B_t, time_step, nx_max, ny_max, cell_area, outlet_index, idx_nan_5, flag_wq_model):
    import cupy as cp
    # Pollutan Rate Flux
    # ----------- Vector of Flows ------------- #
    # 3-D array with layers 0, 1, 2, 3, 4, and 5 being left, right, up, down, and
    # outlet flow rates in mm/hr
    q_out_t = cp.stack([qout_left_t, qout_right_t, qout_up_t, qout_down_t, outlet_flow])
    # ----------- Choosing Which Wash-off Equation ---------- #
    if flag_wq_model == 1:  # Rating Curve
        W_out_t = cp.multiply((C_3), cp.power(q_out_t, C_4))  # Kg/hr (rating curve)
    else:
        W_out_t = cp.multiply(C_3, cp.multiply(cp.power(q_out_t, C_4), B_t))  # kg/hr (mass based curve)
    # ----------- Total Outflows per Cell ---------- #
    tot_W_out = cp.sum(W_out_t, axis=0)  # kg/hr
    tot_q_out = cp.sum(q_out_t, axis=0)  # mm/hr
    small_number = 0.000001
    small_number = cp.asarray(small_number)

    # ----------- Replacing Nans to 0 ---------- #
    tot_W_out = cp.where(cp.isnan(tot_W_out) == True, 0, tot_W_out)
    W_out_t = cp.where(cp.isnan(W_out_t) == True, 0, W_out_t)
    tot_q_out = cp.where(cp.isnan(tot_q_out) == True, 0, tot_q_out)
    # ----------- Assigning Wash-off Rates for Each Direction ---------- #
    W_out_left_t = W_out_t[0, :, :]
    W_out_right_t = W_out_t[1, :, :]
    W_out_up_t = W_out_t[2, :, :]
    W_out_down_t = W_out_t[3, :, :]
    W_out_outlet_t = W_out_t[4, :, :]
    # ----------- Creating Pollutant Inflow Matrices ---------- #
    W_in_left_t = cp.hstack([cp.zeros((int(ny_max), 1)), W_out_right_t[:, 0:(nx_max-1)]])  # Left
    W_in_right_t = cp.hstack([W_out_left_t[:, 1:nx_max], cp.zeros((int(ny_max), 1))])  # Right
    W_in_up_t = cp.vstack([cp.zeros((1, int(nx_max))), W_out_down_t[0:(ny_max-1), :]])  # Up
    W_in_down_t = cp.vstack([W_out_up_t[1: ny_max], cp.zeros((1, int(nx_max)))])  # Down
    W_in_t = cp.stack([W_in_left_t, W_in_right_t, W_in_up_t, W_in_down_t])  # outlet
    # ----------- Net Pollutant Rate ---------- #
    # At the same time, pollutants are entering and leaving the cells. the
    # matrix dW measures the difference between inflows and outflows
    dW = (cp.sum(W_in_t, axis=0)) - tot_W_out  # inflow pol - Outflow pol (kg/hr)
    dW = cp.where(cp.isnan(dW) == True, 0, dW)
    B_t = cp.where(cp.isnan(B_t) == True, 0, B_t)
    # %%% Minimum Time-Step (adaptative time-step) %%% #
    if cp.max(dW) <= 0:
        tmin_wq = 10
    elif cp.min(B_t) == 0:
        # Find values equal 0
        zzzz = B_t.copy()
        zzzz[zzzz <= 0] = math.inf
        dW_zzz = dW.copy()
        dW_zzz[dW_zzz < 0] = 0
        tmin_wq = 3600 * cp.nanmin(cp.divide(zzzz,(dW_zzz + small_number)))   # seconds
    else:
        dW_zzz = dW.copy()
        dW_zzz[dW_zzz < 0] = 0
        tmin_wq = 3600 * cp.nanmin(cp.divide(B_t, (dW_zzz + small_number)))   # seconds

    # %%%% Overall Mass Balance at Cells %%%% #
    B_t = B_t + dW*time_step/60  # Refreshing Pollutant Mass
    # --------- Imposing Constraint at B_t ------------- #
    # this constraint rounds any inf or Nan at B_t to 0
    # Since we are changing time-steps each 60 seconds (usually), we could
    # experience instability within this period. To this end, we force any
    # weird number to 0
    # Another important consideration would be round very tiny pollutant
    # massess to zero. So let's assume we have a minimum value in kg/m2
    # If we have a minimum of 0.01 kg/ha, we would have 10^-6 kg/m2 as a
    # minimum bound. Therefore, if we round B_t to the 6th decimal place, we
    # sole this problem.
    B_t[B_t < 0.00001] = 0
    # Adding a constraint in q_out to avoid huge numbers
    # ---------- Imposing Constraint at Flow Rate ------------ #
    tot_q_out = cp.maximum(tot_q_out, 10)  # 10 mm/hr
    # ---------- Pollutant Concentration --------------------- #
    P_conc = cp.maximum(cp.divide((cp.power(10, 6) * tot_W_out), (tot_q_out*cell_area)), 0)  # mg/L
    idx = outlet_index == 1
    # ---------- Outlet Concentration ------------------------ #
    Out_Conc = cp.maximum(cp.divide(1000*cp.sum(W_out_outlet_t[idx]), cp.divide(cp.sum(outlet_flow), 1000) * cell_area), 0)  # mg/L

    return B_t, P_conc, Out_Conc, tmin_wq, tot_W_out