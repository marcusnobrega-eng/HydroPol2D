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
# This function calculates the average of weigth, height and age according SISVAN data.
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
import pandas as pd
import numpy as np

# Read the CSV file into a DataFrame
path = 'I:/Meu Drive/Papers/Paper - 2Dmodeling + human risk/SISVAN/'
encodings = ['utf-8', 'latin-1', 'iso-8859-1']
df = pd.read_csv(path+'sisvan_estado_nutricional_2021.csv', delimiter=';', encoding='latin-1')

# Filter the data based on the specified conditions
filtered_df = df[(df['CO_MUNICIPIO_IBGE'] == 'SP') & (df['NU_FASE_VIDA'].isin([8, 7, 6, 4]))]

# Group the filtered data by 'CO_MUNICIPIO_IBGE'
grouped_df = filtered_df.groupby('CO_MUNICIPIO_IBGE')

# Iterate over the groups and perform desired operations
for group_name, group_data in grouped_df:
    print(f"Group: {group_name}")
    print(group_data)
    print()

# If you want to perform further operations on the grouped data, you can access it using 'grouped_df'
# For example, you can calculate the sum of a specific column in each group:
sum_column = grouped_df['column_name'].sum()

