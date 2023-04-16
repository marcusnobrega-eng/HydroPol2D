###################################################################
#                                                                 #
#                        Produced by:                             #
#                 Marcus Nobrega Junior gomes                     #
#              (marcusnobrega.engcivil@gmail.com)                 #
#                               &                                 #
#                 Luis Miguel Castillo Rápalo                     #
#                   (luis.castillo@unah.hn)                       #
#                       September 2021                            #
#                                                                 #
#       Last update : 11 september, 2021                          #
###################################################################

# Calculates accumulated and incremental variable, assuming the model's
# time-step and the timeseries time-step
# 09/30/2020

import numpy as np

# delta_p is the variable discretized into model's time-step

def accumulated_incremental(steps, variable_discretized, time_step_variable):
    delta_variable = np.zeros((int(steps)))
    accum_variable = np.cumsum(variable_discretized[0,:]*time_step_variable/60)

    for t in np.arange(steps):
        if t == 0:
            delta_variable[int(t)] = accum_variable[int(t)]
        elif t >= len(variable_discretized[0,:]):
            delta_variable[int(t)] = 0
        else:
            delta_variable[int(t)] = accum_variable[int(t)]-accum_variable[int(t)-1] #mm

    variable_intensity = delta_variable/(time_step_variable/60) # mm/h trhoughall time-steps
    delta_variable = np.round(delta_variable, 6) #ROUNDED VALUE, CAREFULL HERE!

    return delta_variable, variable_intensity