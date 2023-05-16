import pandas as pd

#  File directory of the huge csv
file_path = 'E:/Google_drive/Meu Drive/Papers/Paper - 2Dmodeling + human risk/SISVAN/sisvan_estado_nutricional_2021.csv'
# The desired columns to extract and not to load all the data
columns = ['NO_MUNICIPIO', 'NU_IDADE_ANO', 'NU_FASE_VIDA', 'SG_SEXO', 'NU_PESO',  'NU_ALTURA',  'DS_IMC']
# Create a new dataframe to store the filtered data
database = pd.DataFrame()
summary = pd.DataFrame(columns=['Person'])
# Municipal code of São Paulo according the IBGE
SP_IBGE = 'SAO PAULO'
# Main loop to load, filter and store the data from the csv
nrows = 0
for chunk in pd.read_csv(file_path, chunksize=100000, sep=';', encoding='latin1'):
    nrows += len(chunk)
    temp = pd.concat((database, chunk.loc[:, columns]), axis=0)
    database = temp[temp['NO_MUNICIPIO']==SP_IBGE]
    print('Proceessed {0}'.format(nrows))

# Fix commna separator for decimals
database['NU_PESO'] = database['NU_PESO'].str.replace(',', '.').astype(float)
database['NU_FASE_VIDA'] = database['NU_FASE_VIDA'].str.replace(',', '.').astype(float)
database['NU_ALTURA'] = database['NU_ALTURA'].str.replace(',', '.').astype(float)
database['DS_IMC'] = database['DS_IMC'].str.replace(',', '.').astype(float)

# Mean of height, weight and IMC
summary['Person'] = ['mas_child', 'mas_teen', 'mas_adult', 'mas_older', 'fem_child', 'fem_teen', 'fem_adult', 'fem_older']
metrics = ['NU_PESO','NU_ALTURA','DS_IMC']

for i in metrics:
    summary.loc[0, i] = database[(database['SG_SEXO'] == 'M') & (database['NU_FASE_VIDA'] == 4)][i].mean()
    summary.loc[1, i] = database[(database['SG_SEXO'] == 'M') & (database['NU_FASE_VIDA'] == 6)][i].mean()
    summary.loc[2, i] = database[(database['SG_SEXO'] == 'M') & (database['NU_FASE_VIDA'] == 7)][i].mean()
    summary.loc[3, i] = database[(database['SG_SEXO'] == 'M') & (database['NU_FASE_VIDA'] == 8)][i].mean()
    summary.loc[4, i] = database[(database['SG_SEXO'] == 'F') & (database['NU_FASE_VIDA'] == 4)][i].mean()
    summary.loc[5, i] = database[(database['SG_SEXO'] == 'F') & (database['NU_FASE_VIDA'] == 6)][i].mean()
    summary.loc[6, i] = database[(database['SG_SEXO'] == 'F') & (database['NU_FASE_VIDA'] == 7)][i].mean()
    summary.loc[7, i] = database[(database['SG_SEXO'] == 'F') & (database['NU_FASE_VIDA'] == 8)][i].mean()

summary.loc[0, 'Count'] = len(database[(database['SG_SEXO'] == 'M') & (database['NU_FASE_VIDA'] == 4)][i])
summary.loc[1, 'Count'] = len(database[(database['SG_SEXO'] == 'M') & (database['NU_FASE_VIDA'] == 6)][i])
summary.loc[2, 'Count'] = len(database[(database['SG_SEXO'] == 'M') & (database['NU_FASE_VIDA'] == 7)][i])
summary.loc[3, 'Count'] = len(database[(database['SG_SEXO'] == 'M') & (database['NU_FASE_VIDA'] == 8)][i])
summary.loc[4, 'Count'] = len(database[(database['SG_SEXO'] == 'F') & (database['NU_FASE_VIDA'] == 4)][i])
summary.loc[5, 'Count'] = len(database[(database['SG_SEXO'] == 'F') & (database['NU_FASE_VIDA'] == 6)][i])
summary.loc[6, 'Count'] = len(database[(database['SG_SEXO'] == 'F') & (database['NU_FASE_VIDA'] == 7)][i])
summary.loc[7, 'Count'] = len(database[(database['SG_SEXO'] == 'F') & (database['NU_FASE_VIDA'] == 8)][i])

summary.to_excel('E:/Google_drive/Meu Drive/Papers/Paper - 2Dmodeling + human risk/SISVAN/summary.xlsx')
