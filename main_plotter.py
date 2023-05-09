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

import pandas as pd
import numpy as np
import matplotlib.pylab as plt

# --- Data Path --- #
# It should be provided the directory of the events to be simulated by the model
path ='I:/Meu Drive/Papers/Paper - 2Dmodeling + loss curves/runs'  # Provide the path for the descrete events
# runs = sorted(os.listdir(path))
runs = ['HydroPol2D']

# --- Plotting results --- of Calibrated and validated runs --- #

from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Garamond']}, size = '13')
rc('text', usetex=False)

summary = np.empty((len(runs),4), dtype=object)
fig = plt.figure(figsize=(8, 8))
spec = fig.add_gridspec(5, 2)
row, col, i = 0, 0, 0
letter= ['a)', 'b)', 'c)', 'd)', 'e)', 'f)', 'g)','h)','i)','j)']
for run in runs:
    obs = pd.read_csv('I:/Meu Drive/Papers/Paper - 2Dmodeling + loss curves/runs/HydroPol2D/observed_data/level_rain2010-10-16-20-40-00.csv', delimiter=',', skiprows=2)
    cal = pd.read_csv('I:/Meu Drive/Papers/Paper - 2Dmodeling + loss curves/runs/HydroPol2D/Outlet_hydrograph_Data.txt', delimiter=',')
    # height = 4.1  # m
    streamflow = obs['Unnamed: 3']
    streamflow_estimated = cal['Discharge (m3/s)']
    if run == runs[0]:
         axs = fig.add_subplot(spec[row, :])
         axs.plot(obs['DATA'].astype(str).str[5:-3], streamflow, 'black', linewidth=1.5, label='$Observed$')
         axs.plot(obs['DATA'].astype(str).str[5:-3], streamflow_estimated, 'red', linewidth=1, label='$Simulated$', color = 'tab:red')
         axs_r = axs.twinx()
         axs_r.bar(obs['DATA'].astype(str).str[5:-3], obs['Unnamed: 2'], color='dodgerblue', width=1, alpha=0.7, label="$Rainfall$")
         fig.legend(prop={'size': 10}, loc='upper right', bbox_to_anchor=(1, 1), bbox_transform=axs_r.transAxes)
         axs.xaxis.set_major_locator(plt.MaxNLocator(7))
         axs_r.xaxis.set_major_locator(plt.MaxNLocator(7))
         axs_r.set(ylabel='Rainfall Intensity (mm/h)')
    else:
         axs = fig.add_subplot(spec[row, col])
         axs.plot(obs['DATA'].astype(str).str[5:-3], streamflow, 'black', linewidth=1.5)
         axs.plot(obs['DATA'].astype(str).str[5:-3], streamflow_estimated, 'red', linewidth=1, color = 'tab:red')
         axs_r = axs.twinx()
         axs_r.bar(obs['DATA'].astype(str).str[5:-3], obs['Unnamed: 2'], color='dodgerblue', width=1, alpha=0.7)
         axs.xaxis.set_major_locator(plt.MaxNLocator(3))
         axs_r.xaxis.set_major_locator(plt.MaxNLocator(3))
         axs.tick_params(axis='x', labelrotation=0)
         plt.tight_layout()
    axs_r.set_axisbelow(True)
    axs_r.set_ylim([0, obs['Unnamed: 2'].max() * 3])
    axs_r.invert_yaxis()
    axs.annotate(letter[i], (-0.12, 1.05), xycoords='axes fraction', fontsize=16, weight='bold')
    axs.set_title(run[3:5]+'-'+run[0:2]+'-'+ run[6:10])
    # summary[i,0] = run
    # summary[i,1] = np.round(nse(streamflow, streamflow_estimated['Discharge (m3/s)']), 2)  # NSE
    # summary[i,2] = np.round(r2_score(streamflow, streamflow_estimated['Discharge (m3/s)']),2)
    # summary[i,3] = np.round(he.pbias(np.asarray(streamflow),np.asarray(streamflow_estimated['Discharge (m3/s)'])),2)
    i += 1
    if row == 0:
         row += 1
         axs.set(ylabel='Flow Discharge $\mathregular{(m^3/s)}$')
         continue
    elif col == 1:
         axs_r.set(ylabel='Rainfall Intensity (mm/h)')
         row += 1
         col = 0
         continue
    else:
         col = 1
         axs.set(ylabel='Flow Discharge $\mathregular{(m^3/s)}$')


fig.savefig('I:/Meu Drive/Papers/Paper - 2Dmodeling + level calibration/figures/'+ 'Calibration_100dpi.jpg', format='jpg', dpi=100)
fig.savefig('I:/Meu Drive/Papers/Paper - 2Dmodeling + level calibration/figures/'+ 'calibration_600dpi.jpg', format='jpg', dpi=600)
