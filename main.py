###################################################################
#                                                                 #
#                        Produced by:                             #
#                 Marcus Nobrega Junior gomes                     #
#              (marcusnobrega.engcivil@gmail.com)                 #
#                               &                                 #
#                 Luis Miguel Castillo Rápalo                     #
#                   (luis.castillo@unah.hn)                       #
#                         February 2023                           #
#                                                                 #
#                  Last Updated: 2 February23                     #
#                                                                 #
###################################################################

#§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§
#   -----  DESCRIPTION  -----
#   Name : HydroPol2D
#   Objetive: Hydrodynamic modeling for discrete rainfall events
#   version: 1
#§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§§

# --- requeried libraries --- ##
import math
import os
import pandas as pd
import numpy as np
import model_GPU
import model
import hydroeval as he
import matplotlib.pylab as plt
from matplotlib import rc
import matplotlib as mpl
from sklearn.metrics import r2_score
import cProfile
import main_plotter



# --- Nash-Sutcliffe Efficiency (NSE) --- #
def nse(targets,predictions):
    return 1-(np.sum((targets-predictions)**2)/np.sum((targets-np.mean(targets))**2))

# --- Data Path --- #
# It should be provided the directory of the events to be simulated by the model
path ='I:/Meu Drive/Papers/Paper - 2Dmodeling + loss curves/runs'  # Provide the path for the descrete events
# runs = sorted(os.listdir(path))
runs = ['HydroPol2D']

# --- Non-Parallel Running --- #  # Use if only the CPU is available (run slow)
# for run in runs:
#     path_run = path + '/' + run
#     model.model(path_run)
#     print(run + ' = ok')

# --- Parallel Running --- #  # Use if the GPU is available (run faster)
# with cProfile.Profile() as profile:
for run in runs:
    path_run = path + '/' + run
    model_GPU.model_gpu(path_run)
    print(run + ' = ok')
# results = pstats.Stats(profile)
# results.sort_stats(pstats.SortKey.TIME)
# results.print_stats()



# --- Plotting results --- of Calibrated and validated runs --- #
# summary = np.empty((len(runs),4), dtype=object)
# fig = plt.figure(figsize=(8, 8))
# spec = fig.add_gridspec(4, 2)
# row, col, i = 0, 0, 0
# letter= ['$A$', '$B$', '$C$', '$D$', '$E$', '$F$', '$G$']
# for run in runs:
#     obs = pd.read_csv(path + '/' + run + '/level_rain.csv', delimiter=',')
#     cal = pd.read_csv(path + '/' + run + '/Outlet_hydrolevels_Data.txt', delimiter=',')
#     slope = 0.00397
#     manning = 0.015
#     base = 5  # m
#     height = 4.1  # m
#     streamflow = (1 / manning) * (base * (obs.iloc[:, 1] - np.min(obs['Water Height (masl)']))) * np.power(
#          (base * height) / (base + 2 * height), 2 / 3) * np.power(slope, 1 / 2)
#     streamflow_estimated = pd.read_csv(path + '/' + run + '/Outlet_hydrograph_Data.txt', delimiter=',')
#     if run == runs[0]:
#          axs = fig.add_subplot(spec[row, :])
#          axs.plot(obs['Date/time'].astype(str).str[5:-3], streamflow, 'black', linewidth=1.5, label='$Observed$')
#          axs.plot(obs['Date/time'].astype(str).str[5:-3], streamflow_estimated['Discharge (m3/s)'], 'red', linewidth=1, label='$Simulated$', color = 'tab:red')
#          axs_r = axs.twinx()
#          axs_r.bar(obs['Date/time'].astype(str).str[5:-3], obs['Rainfall_Intensity (mm/h)'], color='dodgerblue', width=1, alpha=0.7, label="$Rainfall$")
#          fig.legend(prop={'size': 10}, loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=axs_r.transAxes)
#          axs.xaxis.set_major_locator(plt.MaxNLocator(7))
#          axs_r.xaxis.set_major_locator(plt.MaxNLocator(7))
#          axs_r.set(ylabel='$Rainfall\ Intensity\ (mm/h)$')
#     else:
#          axs = fig.add_subplot(spec[row, col])
#          axs.plot(obs['Date/time'].astype(str).str[5:-3], streamflow, 'black', linewidth=1.5)
#          axs.plot(obs['Date/time'].astype(str).str[5:-3], streamflow_estimated['Discharge (m3/s)'], 'red', linewidth=1, color = 'tab:red')
#          axs_r = axs.twinx()
#          axs_r.bar(obs['Date/time'].astype(str).str[5:-3], obs['Rainfall_Intensity (mm/h)'], color='dodgerblue', width=1, alpha=0.7)
#          axs.xaxis.set_major_locator(plt.MaxNLocator(3))
#          axs_r.xaxis.set_major_locator(plt.MaxNLocator(3))
#          axs.tick_params(axis='x', labelrotation=0)
#          plt.tight_layout()
#     axs_r.set_axisbelow(True)
#     axs_r.set_ylim([0, obs['Rainfall_Intensity (mm/h)'].max() * 3])
#     axs_r.invert_yaxis()
#     axs.annotate(letter[i], (0.05, 0.85), xycoords='axes fraction', bbox=dict(facecolor='none', boxstyle='round,pad=0.2'))
#     axs.set_title(run[3:5]+'-'+run[0:2]+'-'+ run[6:10])
#     summary[i,0] = run
#     summary[i,1] = np.round(nse(streamflow, streamflow_estimated['Discharge (m3/s)']), 2)  # NSE
#     summary[i,2] = np.round(r2_score(streamflow, streamflow_estimated['Discharge (m3/s)']),2)
#     summary[i,3] = np.round(he.pbias(np.asarray(streamflow),np.asarray(streamflow_estimated['Discharge (m3/s)'])),2)
#     i += 1
#     if row == 0:
#          row += 1
#          axs.set(ylabel='$Flow\ Discharge\ (m^3/s)$')
#          continue
#     elif col == 1:
#          axs_r.set(ylabel='$Rainfall\ Intensity\ (mm/h)$')
#          row += 1
#          col = 0
#          continue
#     else:
#          col = 1
#          axs.set(ylabel='$Flow\ Discharge\ (m^3/s)$')
#
#
# fig.savefig('I:/Meu Drive/Papers/Paper - 2Dmodeling + level calibration/figures/'+ 'Calibration_100dpi.jpg', format='jpg', dpi=100)
# fig.savefig('I:/Meu Drive/Papers/Paper - 2Dmodeling + level calibration/figures/'+ 'calibration_600dpi.jpg', format='jpg', dpi=600)

# --- Plotting results --- of sensitivity runs --- #

# path ='I:/Meu Drive/Papers/Paper - 2Dmodeling + level calibration/591_runs/sensitivity/v_treshold'
# runs = sorted(os.listdir(path))
# #
# # # # # --- Plotting results --- of sensitive runs V_treshold #
#
# rc_fonts = {
#     "text.usetex": False,
# }
# plt.rcParams.update(rc_fonts)
#
# fig = plt.figure(figsize=(8.27, 4))
# spec = fig.add_gridspec(2, 3)
# i = 0
# slope = 0.00397
# manning = 0.035
# base = 10  # m
# height = 4.1  # m
# letter = ['A', 'D']
# # labels = ['Slope 90%',
# #          'Slope 15%',
# #          'Slope 25%',
# #          'Slope 50%',
# #          'Slope 75%']
# labels = [r'$Vel_t\ >\ 0.5\ m/s$',
#          r'$Vel_t\ >\ 1\ m/s$',
#          r'$Vel_t\ >\ 3\ m/s$',
#          r'$Vel_t\ >\ 5\ m/s$',
#          r'$Vel_t\ >\ 10\ m/s$']
# # labels = ['Alpha min 0.05',
# #          'Alpha min 0.1',
# #          'Alpha min 0.2',
# #          'Alpha min 0.3',
# #          'Alpha min 0.5']
# color = ['tab:red', 'tab:blue', 'tab:green', 'tab:orange', 'tab:grey']
# style = ['solid','solid','solid','solid','solid']
# marker = ['p', 'o', 'v', 's', '*']
# markersize = 3
# times = []
# metrics = []
# for run in runs:
#      obs = pd.read_csv(path + '/' + run + '/level_rain.csv', delimiter=',')
#      cal = pd.read_csv(path + '/' + run + '/Outlet_hydrolevels_Data.txt', delimiter=',')
#      streamflow = (1 / manning) * (base * (obs.iloc[:, 1] - np.min(obs['Water Height (masl)']))) * np.power((base*height)/(base + 2*height), 2/3) * np.power(slope, 1 / 2)
#      streamflow_estimated = np.asarray(pd.read_csv(path + '/' + run + '/Outlet_hydrograph_Data.txt', delimiter=','))
#      velocities = np.asarray(pd.read_csv(path + '/' + run + '/alpha_max_vel_time_step.txt', delimiter=','))
#      alfa = velocities[1:299, 0]
#      velocities = velocities[1:299, 2]
#      t = np.linspace(1, 298, num=298)
#      time = pd.read_csv(path + '/' + run + '/Summary_table.txt', delimiter=',')
#      time = round(time.iloc[3, 1]/60, 2)  # How long was the simulation time
#      label = '\nNSE: {}'.format(np.round(nse(streamflow, streamflow_estimated[:, 1]), 2))
#      times.append((time))
#      metrics.append((np.round(nse(obs['Water Height (masl)'] - obs['Water Height (masl)'].min(), cal['levels (m)'] - cal['levels (m)'].min()), 2)))
#      if run == runs[0]:
#           axs3 = fig.add_subplot(spec[1, 0])
#           axs2 = fig.add_subplot(spec[0, 0])
#
#           axs3.plot(obs['Date/time'].astype(str).str[5:-3], streamflow, 'black', linewidth=1.5, label= '$Observed$')
#           axs3.plot(obs['Date/time'].astype(str).str[5:-3], streamflow_estimated[:,1], color=color[i], linestyle= style[i],linewidth=1)
#           axs_r = axs3.twinx()
#           axs_r.bar(obs['Date/time'].astype(str).str[5:-3], obs['Rainfall_Intensity (mm/h)'], color='dodgerblue', width=1, alpha=0.7, label="Rainfall")
#           axs3.xaxis.set_major_locator(plt.MaxNLocator(2))
#           axs_r.xaxis.set_major_locator(plt.MaxNLocator(2))
#           # axs_r.set(ylabel='Rainfall Intensity $(mm/h)$')
#           axs_r.set_axisbelow(False)
#           axs_r.axis('off')
#           axs_r.set_ylim([0, obs['Rainfall_Intensity (mm/h)'].max() * 3])
#           axs_r.invert_yaxis()
#           axs3.set(ylabel=r'$Flow\ Discharge$ $(m^3/s)$')
#
#           axs3.axes.get_xaxis().set_ticks([])
#           axs3.set(xlabel=r'$Time$ $(hrs)$')
#
#           axs2.plot(t, velocities*60, color=color[i], linestyle= style[i], linewidth=1.5, label=labels[i])
#           axs2.set(xlabel=r'$Maximum\ Velocity$ $(m/s)$')
#
#           axs2.set(ylabel=r'$Time-Step\ (s)$')
#           axs2.axes.get_xaxis().set_visible(False)
#
#           i = i + 1
#      else:
#           axs3.plot(obs['Date/time'].astype(str).str[5:-3], streamflow_estimated[:,1], color=color[i], linestyle= style[i],  linewidth=1.5)
#           axs2.plot(t, velocities*60, color=color[i], linestyle= style[i], linewidth=1.5, label=labels[i])
#           # axs2.plot(obs['Date/time'].astype(str).str[5:-3], streamflow_estimated['Discharge (m3/s)'], color=color[i], linestyle=style[i], linewidth=1)
#           i = i + 1
#
#     plt.tight_layout()
# axs3.annotate("$D$", (0.05, 0.85), xycoords='axes fraction',
#               bbox=dict(facecolor='none', boxstyle='round,pad=0.2'))
# axs2.annotate("$A$", (0.3, 0.85), xycoords='axes fraction',
#               bbox=dict(facecolor='none', boxstyle='round,pad=0.2'))
# axs3.legend(prop={'size': 7},  loc='upper right', bbox_transform=axs3.transAxes)
# axs2.legend(prop={'size': 7},  loc='upper right', bbox_transform=axs2.transAxes)
# # axs3.legend(prop={'size': 7},  loc='upper right', bbox_transform=axs2.transAxes)
# #
# #
# path ='I:/Meu Drive/Papers/Paper - 2Dmodeling + level calibration/591_runs/sensitivity/slope'
# runs = sorted(os.listdir(path))
# # --- Plotting results --- of sensitive runs slope #
# # fig = plt.figure(figsize=(8, 5))
# # spec = fig.add_gridspec(2, 1)
# i = 0
# slope = 0.00397
# manning = 0.035
# base = 10  # m
# height = 4.1  # m
# letter = ['B', 'E']
# labels = [r'$\alpha\ 15\%$',
#          r'$\alpha\ 25\%$',
#          r'$\alpha\ 50\%$',
#          r'$\alpha\ 75\%$',
#          r'$\alpha\ 90\%$']
# # labels = ['Vel > 0.5 $m/s$',
# #          'Vel > 1 $m/s$',
# #          'Vel > 3 $m/s$',
# #          'Vel > 5 $m/s$',
# #          'Vel > 10 $m/s$']
# # labels = ['Alpha min 0.05',
# #          'Alpha min 0.1',
# #          'Alpha min 0.2',
# #          'Alpha min 0.3',
# #          'Alpha min 0.5']
# color = ['tab:red', 'tab:blue', 'tab:green', 'tab:orange', 'tab:grey']
# style = ['solid','solid','solid','solid','solid']
# marker = ['p', 'o', 'v', 's', '*']
# markersize = 3
# times = []
# metrics = []
# for run in runs:
#      obs = pd.read_csv(path + '/' + run + '/level_rain.csv', delimiter=',')
#      cal = pd.read_csv(path + '/' + run + '/Outlet_hydrolevels_Data.txt', delimiter=',')
#      streamflow = (1 / manning) * (base * (obs.iloc[:, 1] - np.min(obs['Water Height (masl)']))) * np.power((base*height)/(base + 2*height), 2/3) * np.power(slope, 1 / 2)
#      streamflow_estimated = np.asarray(pd.read_csv(path + '/' + run + '/Outlet_hydrograph_Data.txt', delimiter=','))
#      velocities = np.asarray(pd.read_csv(path + '/' + run + '/alpha_max_vel_time_step.txt', delimiter=','))
#      alfa = velocities[1:299, 0]
#      velocities = velocities[1:299, 2]
#      t = np.linspace(1, 298, num=298)
#      time = pd.read_csv(path + '/' + run + '/Summary_table.txt', delimiter=',')
#      time = round(time.iloc[3, 1]/60, 2)  # How long was the simulation time
#      label = '\nNSE: {}'.format(np.round(nse(streamflow, streamflow_estimated[:, 1]), 2))
#      times.append((time))
#      metrics.append((np.round(nse(obs['Water Height (masl)'] - obs['Water Height (masl)'].min(), cal['levels (m)'] - cal['levels (m)'].min()), 2)))
#      if run == runs[0]:
#           axs8 = fig.add_subplot(spec[1, 1])
#           axs4 = fig.add_subplot(spec[0, 1])
#
#           axs8.plot(obs['Date/time'].astype(str).str[5:-3], streamflow, 'black', linewidth=1.5)
#           axs8.plot(obs['Date/time'].astype(str).str[5:-3], streamflow_estimated[:,1], color=color[i], linestyle= style[i],linewidth=1.5, label=labels[i])
#           axs_r = axs8.twinx()
#           axs_r.bar(obs['Date/time'].astype(str).str[5:-3], obs['Rainfall_Intensity (mm/h)'], color='dodgerblue', width=1, alpha=0.7, label="Rainfall")
#           axs8.xaxis.set_major_locator(plt.MaxNLocator(2))
#           axs_r.xaxis.set_major_locator(plt.MaxNLocator(2))
#           # axs_r.set(ylabel='Rainfall Intensity $(mm/h)$')
#           axs_r.set_axisbelow(False)
#           # axs_r.axis('off')
#           axs_r.set_ylim([0, obs['Rainfall_Intensity (mm/h)'].max() * 3])
#           axs_r.invert_yaxis()
#           # axs.set(ylabel='Flow Discharge $(m^3/s)$')
#           axs8.get_yaxis().set_visible(False)
#           axs_r.get_yaxis().set_visible(False)
#           axs8.axes.get_xaxis().set_ticks([])
#           axs8.set(xlabel='$Time$ $(hrs)$')
#
#           axs4.plot(t, velocities*60, color=color[i], linestyle= style[i], linewidth=1.5, label=labels[i])
#           # axs4.set(xlabel='$Maximum-Velocity$ $(m/s)$')
#           axs4.axes.get_xaxis().set_visible(False)
#           # axs2.set(ylabel='Alpha Value')
#           axs4.get_yaxis().set_visible(False)
#
#           # axs2.xaxis.set_major_locator(plt.MaxNLocator(5))
#
#           i = i + 1
#      else:
#           axs8.plot(obs['Date/time'].astype(str).str[5:-3], streamflow_estimated[:,1], color=color[i], linestyle= style[i],  linewidth=1.5)
#           axs4.plot(t, velocities*60, color=color[i], linestyle= style[i], linewidth=1.5, label=labels[i])
#           # axs2.plot(obs['Date/time'].astype(str).str[5:-3], streamflow_estimated['Discharge (m3/s)'], color=color[i], linestyle=style[i], linewidth=1)
#           i = i + 1
#
#      plt.tight_layout()
# axs8.annotate("$E$", (0.05, 0.85), xycoords='axes fraction',
#               bbox=dict(facecolor='none', boxstyle='round,pad=0.2'))
# axs4.annotate("$B$", (0.3, 0.85), xycoords='axes fraction',
#               bbox=dict(facecolor='none', boxstyle='round,pad=0.2'))
# axs4.legend(prop={'size': 7}, loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=axs4.transAxes)
#
#
# path ='I:/Meu Drive/Papers/Paper - 2Dmodeling + level calibration/591_runs/sensitivity/alfa_min'
# runs = sorted(os.listdir(path))
#
# # # --- Plotting results --- of sensitive runs alfa_min #
# # fig = plt.figure(figsize=(8, 5))
# # spec = fig.add_gridspec(1, 2)
# i = 0
# slope = 0.00397
# manning = 0.035
# base = 10  # m
# height = 4.1  # m
# letter = ['C', 'F']
# # labels = ['Slope 15%',
# #          'Slope 25%',
# #          'Slope 50%',
# #          'Slope 75%',
# #          'Slope 90%']
# # labels = ['Vel > 0.5 $m/s$',
# #          'Vel > 1 $m/s$',
# #          'Vel > 3 $m/s$',
# #          'Vel > 5 $m/s$',
# #          'Vel > 10 $m/s$']
# labels = ['$CFL_\mathrm{min}\ =\ 0.05$',
#          '$CFL_\mathrm{min}\ =\ 0.1$',
#          '$CFL_\mathrm{min}\ =\ 0.2$',
#          '$CFL_\mathrm{min}\ =\ 0.3$',
#          '$CFL_\mathrm{min}\ =\ 0.5$']
# color = ['tab:red', 'tab:blue', 'tab:green', 'tab:orange', 'tab:grey']
# style = ['solid','solid','solid','solid','solid']
# marker = ['p', 'o', 'v', 's', '*']
# markersize = 3
# times = []
# metrics = []
# for run in runs:
#      obs = pd.read_csv(path + '/' + run + '/level_rain.csv', delimiter=',')
#      cal = pd.read_csv(path + '/' + run + '/Outlet_hydrolevels_Data.txt', delimiter=',')
#      streamflow = (1 / manning) * (base * (obs.iloc[:, 1] - np.min(obs['Water Height (masl)']))) * np.power((base*height)/(base + 2*height), 2/3) * np.power(slope, 1 / 2)
#      streamflow_estimated = np.asarray(pd.read_csv(path + '/' + run + '/Outlet_hydrograph_Data.txt', delimiter=','))
#      velocities = np.asarray(pd.read_csv(path + '/' + run + '/alpha_max_vel_time_step.txt', delimiter=','))
#      alfa = velocities[1:299, 0]
#      velocities = velocities[1:299, 2]
#      t = np.linspace(1, 298, num=298)
#      time = pd.read_csv(path + '/' + run + '/Summary_table.txt', delimiter=',')
#      time = round(time.iloc[3, 1]/60, 2)  # How long was the simulation time
#      label = '\nNSE: {}'.format(np.round(nse(streamflow, streamflow_estimated[:, 1]), 2))
#      times.append((time))
#      metrics.append((np.round(nse(obs['Water Height (masl)'] - obs['Water Height (masl)'].min(), cal['levels (m)'] - cal['levels (m)'].min()), 2)))
#      if run == runs[0]:
#           axs9 = fig.add_subplot(spec[1, 2])
#           axs5 = fig.add_subplot(spec[0, 2])
#
#           axs9.plot(obs['Date/time'].astype(str).str[5:-3], streamflow, 'black', linewidth=1.5)
#           axs9.plot(obs['Date/time'].astype(str).str[5:-3], streamflow_estimated[:,1], color=color[i], linestyle= style[i],linewidth=1.5, label=labels[i])
#           axs_r = axs9.twinx()
#           axs_r.bar(obs['Date/time'].astype(str).str[5:-3], obs['Rainfall_Intensity (mm/h)'], color='dodgerblue', width=1, alpha=0.7, label="Rainfall")
#           axs9.xaxis.set_major_locator(plt.MaxNLocator(2))
#           axs_r.xaxis.set_major_locator(plt.MaxNLocator(2))
#           axs_r.set(ylabel='$Rainfall\ Intensity\ (mm/h)$')
#           # axs_r.set_axisbelow(False)
#           axs9.get_yaxis().set_visible(False)
#           axs_r.set_ylim([0, obs['Rainfall_Intensity (mm/h)'].max() * 3])
#           axs_r.invert_yaxis()
#           # axs.set(ylabel='Flow Discharge $(m^3/s)$')
#           axs9.axes.get_xaxis().set_ticks([])
#           axs9.set(xlabel='$Time$ $(hrs)$')
#           axs5.axes.get_xaxis().set_visible(False)
#
#           axs5.plot(t, velocities*60, color=color[i], linestyle= style[i], linewidth=1.5, label=labels[i])
#           # axs5.set(xlabel='$Maximumz Velocity$ $(m/s)$')
#           # axs2.set(ylabel='Alpha Value')
#           axs5.get_yaxis().set_visible(False)
#           # axs2.xaxis.set_major_locator(plt.MaxNLocator(5))
#
#           i = i + 1
#      else:
#           axs9.plot(obs['Date/time'].astype(str).str[5:-3], streamflow_estimated[:,1], color=color[i], linestyle= style[i],  linewidth=1.5)
#           axs5.plot(t, velocities*60, color=color[i], linestyle= style[i], linewidth=1.5, label=labels[i])
#           # axs2.plot(obs['Date/time'].astype(str).str[5:-3], streamflow_estimated['Discharge (m3/s)'], color=color[i], linestyle=style[i], linewidth=1)
#           i = i + 1
#
#
# axs9.annotate("$F$", (0.05, 0.85), xycoords='axes fraction',
#               bbox=dict(facecolor='none', boxstyle='round,pad=0.2'))
# axs5.annotate("$C$", (0.3, 0.85), xycoords='axes fraction',
#               bbox=dict(facecolor='none', boxstyle='round,pad=0.2'))
# axs5.legend(prop={'size': 7}, loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=axs5.transAxes)
# axs_r.legend(prop={'size': 7}, loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=axs_r.transAxes)
# plt.tight_layout()
# #
# fig.savefig('I:/Meu Drive/Papers/Paper - 2Dmodeling + level calibration/figures/'+ 'sensitivity_100dpi.jpg', format='jpg', dpi=100)
# fig.savefig('I:/Meu Drive/Papers/Paper - 2Dmodeling + level calibration/figures/'+ 'sensitivity_600dpi.jpg', format='jpg', dpi=600)