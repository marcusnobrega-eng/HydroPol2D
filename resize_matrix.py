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

# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# -------- DESCRIPTION ----------------
# This function created a new expanded output_matrix from a input_matrix
# using a new coordinate system
# domain (mandatory) : expanded domain (zero matrix)
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

import numpy as np

def resize_matrix(input_matrix, domain, ymin_domain_abstracted, ymax_domain_abstracted, xmin_domain_abstracted, xmax_domain_abstracted, dymin_matrix, dymax_matrix, dxmin_matrix, dxmax_matrix, dimension):
    scratch_matrix = input_matrix[ymin_domain_abstracted:ymax_domain_abstracted, xmin_domain_abstracted:xmax_domain_abstracted]
    output_matrix = np.tile(domain, [dimension, 1, 1])
    output_matrix[:, int(dymin_matrix):int(dymax_matrix), int(dxmin_matrix):int(dxmax_matrix)] = scratch_matrix
    return output_matrix