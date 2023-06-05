###################################################################
#                                                                 #
#                        Produced by:                             #
#                 Marcus Nobrega Junior gomes                     #
#              (marcusnobrega.engcivil@gmail.com)                 #
#                               &                                 #
#                 Luis Miguel Castillo Rápalo                     #
#                   (luis.castillo@unah.hn)                       #
#                        March 2023                               #
#                                                                 #
#                 Last update : 1 March, 2023                     #
###################################################################

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# -------- DESCRIPTION ----------------
# This function calculates human risk results from all the runs for C.C.
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# Import required libraries
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

path ='D:/Google_drive/Meu Drive/Papers/Paper - 2Dmodeling + human risk/runs/HydroPol2D_ssp245/Output/'  # Provide the path for the discrete events'  # Provide the path for the discrete events
huff_path = 'D:/Google_drive/Meu Drive/Papers/Congreso SBRH XXV/huff.txt'
cc_path = 'D:/Google_drive/Meu Drive/Papers/Paper - 2Dmodeling + level calibration/CLIMBra/Data_ssp245/mpi_esm1_ssp245_daily.txt'
rainfalls_path = 'D:/Google_drive/Meu Drive/Papers/Paper - 2Dmodeling + human risk/runs/HydroPol2D_ssp245/Rainfall_Intensity_Data.csv'
runs = sorted(os.listdir(path))

# Determining the rainfalls Return Periods
huff = np.array(pd.read_csv(huff_path,sep=','))
cc = pd.read_csv(cc_path, sep='\t', header=None)
cc = np.asarray(cc.groupby(cc[2]).max())
hours = 2 # horus
time_interval = 10 # mins
t_p = np.zeros((int((hours*60)/time_interval)+1,2))

for i in range(len(t_p)):
    t_p[i,0] = i*time_interval/(hours*60)
    if (i == 0):
        t_p[i,1] = 0
        continue
    elif (i == len(t_p)-1):
        t_p[i,1] = 1
        continue
    j1 = np.argwhere(t_p[i,0]>huff[:,0])[-1]
    j2 = np.argwhere(t_p[i,0]<huff[:,0])[0]
    t_p[i,1] = huff[j1,1] - ((huff[j1,1]-huff[j2,1])/(huff[j1,0]-huff[j2,0]))*(huff[j1,0]-t_p[i,0])

# Making the rainfalls
rainfalls_huff = RPs = np.zeros((int((hours*60)/time_interval)+1, int(cc.shape[0])*2))
RPs = np.zeros((int(cc.shape[0]),2))
j = 1
for i in range(int(cc.shape[0])):
    rainfalls_huff[:,j] = t_p[:,1]*cc[i,2]
    rainfalls_huff[1:,j] = (np.diff(rainfalls_huff[:,j]).T)*6
    rainfalls_huff[:,j-1] = t_p[:,0]*hours*60
    RPs[i,0] = np.sum(rainfalls_huff[:,j])/6/hours
    RPs[i,1] = ((RPs[i,0]*(9.227+hours*60)**0.707)/(511.124))**(1/0.158)
    j += 2

rainfalls_huff = pd.DataFrame(rainfalls_huff)
rainfalls_huff.to_csv('D:/Google_drive/Meu Drive/Papers/Congreso SBRH XXV/rainfall_intensity_ssp245', sep=',', index=False)






rainfalls = np.array(pd.read_csv(rainfalls_path, sep=','))
RPs = np.zeros((3,int(rainfalls.shape[1]/2)))

j=1
for i in range(int(rainfalls.shape[1]/2)):
    arg = np.argwhere(rainfalls[:,j] !=0)
    RPs[0, i], RPs[1, i] = arg[0], arg[-1]
    time =  (arg[-1]*10 - arg[0]*10)
    intensity = np.sum(rainfalls[arg,j]/6)/(time/60)
    RPs[2, i] = ((intensity*(9.227+time)**0.707)/(511.124))**(1/0.158)
    j += 2

















