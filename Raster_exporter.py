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
# This function export data to ASCI and GeoTiff format
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
import csv
import glob
import os
from osgeo import gdal

def raster_exporter(FileName ,ncols, nrows, xllcorner, yllcorner, cellsize, no_data_value, raster_exportion):
    with open(str(FileName)+".asc", 'w', newline='') as f:  #
        writer = csv.writer(f, delimiter=" ")
        writer.writerow(["ncols","","","","","","","", str(ncols)])
        writer.writerow(["nrows","","","","","","","", str(nrows)])
        writer.writerow(["xllcorner","","","", str(xllcorner)])
        writer.writerow(["yllcorner","","","", str(yllcorner)])
        writer.writerow(["cellsize","","","","", str(cellsize)])
        writer.writerow(["NODATA_value","", str(no_data_value)])
        writer.writerows(raster_exportion)
    options = gdal.TranslateOptions(noData="-9999")
    tif_1 = gdal.Translate(FileName+".tif", FileName+".asc", outputSRS="EPSG:31983", options=options)
    os.remove(FileName+".asc")
    return