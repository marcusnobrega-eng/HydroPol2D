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
# This function calculates human risk of being dragged by water forces
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

import math
import os
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
from labellines import labelLines
from matplotlib import rc
import matplotlib as mpl
from matplotlib import rc

# define some functions
def submerged(h, Y, d, D, slope):
    if (h / math.cos(slope) < Y / 2):
        V_s = 2 * (h / math.cos(slope)) * math.pi * ((d ** 2) / 4)
        A_s = 2 * (h / math.cos(slope)) * d
        X_gs = d / 2
        Y_gs = h / (2 * math.cos(slope))
        return V_s, A_s, X_gs, Y_gs
    elif (h / math.cos(slope) > Y / 2):
        V_s = 2 * (Y / 2) * math.pi * ((d ** 2) / 4) + ((h / math.cos(slope)) - (Y / 2)) * math.pi * ((D ** 2) / 4)
        A_s = 2 * (Y / 2) * d + ((h / math.cos(slope)) - (Y / 2)) * D
        X_gs = (2 * ((math.pi * d ** 2) / 4) * ((Y / 2) * (d / 2)) + ((math.pi * D ** 2) / 4) * (
                    (h / math.cos(slope)) - (Y / 2)) * (D / 2)) / (
                           2 * ((math.pi * d ** 2) / 4) * (Y / 2) + ((math.pi * D ** 2) / 4) * (
                               (h / math.cos(slope)) - (Y / 2)))
        Y_gs = (2 * ((math.pi * d ** 2) / 4) * ((Y ** 2) / 8) + ((math.pi * D ** 2) / 4) * (
                    (h / math.cos(slope)) - (Y / 2)) * ((0.5 * ((h / math.cos(slope)) - (Y / 2))) + (Y / 2))) / (
                           2 * ((math.pi * d ** 2) / 4) * (Y / 2) + ((math.pi * D ** 2) / 4) * (
                               (h / math.cos(slope)) - (Y / 2)))
        return V_s, A_s, X_gs, Y_gs

plt.style.use('seaborn-darkgrid')
rc('font', **{'family': 'serif', 'serif': ['Garamond']}, size = '12')
rc('text', usetex=False)

summary = pd.read_excel('I:/Meu Drive/Papers/Paper - 2Dmodeling + human risk/SISVAN/summary.xlsx')
names = ['Boy', 'Teenage Boy', 'Man', 'Elderly Man', 'Girl', 'Teenage Girl', 'Women', 'Elderly Women']


# state some variables #
for k in np.arange(0, 1, 1):
    for item in range(8):
        m_p = summary['NU_PESO'].loc[item]
        y_p = summary['NU_ALTURA'].loc[item]/100 #  To meters
        D_p = summary['D'].loc[item]
        d_p = summary['d'].loc[item]

        slope_ = k  # %
        slope = math.atan(slope_/100)
        density = 1000  # kg/m3
        gravity = 9.81  # m/s2
        Cc = 1  # Drag coefficient
        Fc = 0.46  # Friction coefficient

        person = np.empty((len(np.arange(0.01, y_p+0.3, 0.01)), 5))*np.nan
        row = 0
        slide = 1
        toppling = 2
        drowning = 3

        person[:, drowning] = (13/16)*y_p
        for i in np.arange(0.01, y_p+0.3, 0.01):
            person[row, 0] = i
            flag_s = 0
            flag_t = 0
            for j in np.arange(0.0, 5, 0.001):
                Buoy_person, area_person, xgs_person, ygs_person = submerged(i, y_p, d_p, D_p, slope)
                L = 0.5*density*Cc*(math.sin((math.pi/2)-slope)**2)*math.cos((math.pi/2)-slope)*(j**2)*area_person
                D = 0.5*density*Cc*(math.sin((math.pi/2)-slope)**3)*(j**2)*area_person
                Wp = m_p*gravity*math.sin(slope)
                Wn = m_p*gravity*math.cos(slope)
                Bn = density*gravity*Buoy_person*math.cos(slope)
                T = Fc*(Wn-Bn-L)
                if (D + Wp > T) and (flag_s == 0):  # no longer stability by slipping
                    person[row, slide] = j
                    flag_s = 1
                if (D*(i/2)+(Wp*((7/12)*y_p*math.cos(slope)))+(Bn*((xgs_person/math.cos(slope))+(ygs_person*math.sin(slope))))+(L*((xgs_person/math.cos(slope))+((i/2)*(math.tan(slope))))) > Wn*(((5/6)*d_p)/math.cos(slope)+((7/12)*y_p)*math.sin(slope))) and (flag_t==0): # no longer stability by toppling
                    person[row, toppling] = j
                    flag_t = 1
                if flag_s == 1 and flag_t == 1:
                    break
            row += 1

        person[:, 4] = np.minimum(np.nan_to_num(person[:, slide], nan=np.inf),np.nan_to_num(person[:, toppling], nan=np.inf))
        person = np.delete(person, np.argwhere(person[:, 4] == 0)[1:], axis=0)
        person = np.delete(person, np.argwhere(person[:, 4] == np.inf), axis=0)

        # plt.plot(person[:,slide], person[:,0])
        # plt.plot(person[:,toppling], person[:,0])
        # plt.plot(person[:,0], person[:,drowning])
        if item < 4:
            plt.plot(person[:, 4], np.minimum(person[:, 0], person[:, drowning]), label=names[item])
        else:
            plt.plot(person[:, 4], np.minimum(person[:, 0], person[:, drowning]), label=names[item], linestyle='dashed')
    plt.xlim(0,4)
    plt.xlabel('Flow velocity (m/s)')
    plt.ylabel('Flow depth (m)')
    plt.legend()
    plt.savefig('I:/Meu Drive/Papers/Paper - 2Dmodeling + human risk/Figures/risk_levels.png', dpi=600)
# lines = plt.gca().get_lines()
# labelLines(lines, xvals=np.linspace(0.5, 3.5, len(lines)), fontsize=7)











