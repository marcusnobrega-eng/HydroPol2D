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
# This function export data to GIF format
# $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

#libraries
import numpy as np
import xarray as xr
import pyproj
import imageio
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as mcolors
from matplotlib.ticker import FixedLocator
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

def Spatial_data_exporter(FileName ,ncols, nrows, xllcorner, yllcorner, cellsize, no_data_value, raster_exportion):
    # Create a figure and axis object
    # Define the function to update the animation
    if FileName == 'Water Flood Depth':
        fig, ax = plt.subplots()
        def update(frame):
            ax.clear()
            im = ax.imshow(raster_exportion[frame], cmap='jet')
            ax.set_title(str(FileName).format(frame))
            ax.set_xlabel("x (m)")
            ax.set_ylabel("y (m)")
            fig.colorbar(im, cax=cax).set_label('Depths (m)')
        cax = fig.add_axes([0.85,0.15, 0.05, 0.7])
        # Create the animation object
        ani = FuncAnimation(fig, update, frames=raster_exportion.shape[0], interval=300)
        # Save the animation as a GIF
        ani.save(str(FileName)+'.gif', writer='pillow', dpi=600)
        plt.close()
    elif FileName == 'Water Surface Elevation':
        fig, ax = plt.subplots()

        def update(frame):
            ax.clear()
            im = ax.imshow(raster_exportion[frame], cmap='jet')
            ax.set_title(str(FileName).format(frame))
            ax.set_xlabel("x (m)")
            ax.set_ylabel("y (m)")
            fig.colorbar(im, cax=cax).set_label('Depths (m)')

        cax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        # Create the animation object
        ani = FuncAnimation(fig, update, frames=raster_exportion.shape[0], interval=300)
        # Save the animation as a GIF
        ani.save(str(FileName) + '.gif', writer='pillow', dpi=600)
        plt.close()
    elif (FileName == 'Digital Elevation Model'):
        raster = raster_exportion  # Example random raster data
        # Create a figure
        fig = plt.figure()

        def find_max_ignore_inf(arr):
            arr = np.array(arr)
            arr[np.isinf(arr)] = np.nan
            return np.nanmax(arr)

        # Create a 3D plot
        ax = fig.add_subplot(111, projection='3d')
        ax.set_title(str(FileName))
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        # ax.set_zlabel("Elevation\n (masl)")

        # Set the extent of the plot
        x, y = np.meshgrid(np.arange(0, raster.shape[1]), np.arange(0, raster.shape[0]))
        ax.plot_surface(x, y, raster, cmap='gist_earth', edgecolor='none')
        ax.set_zlim(find_max_ignore_inf(raster) * 0.8, find_max_ignore_inf(raster))

        # Define the animation update function
        def update(frame):
            # Rotate the plot around the y-axis
            ax.view_init(elev=30, azim=frame)  # Increase elev to move camera up

        # Create the animation
        ani = FuncAnimation(fig, update, frames=np.linspace(0, 360, 36), interval=200)

        # Save the animation as a GIF
        ani.save('raster_rotation.gif', writer='pillow', dpi=600)
        plt.close()
    elif (FileName == 'Flow Velocity'):
        fig, ax = plt.subplots()
        def update(frame):
            ax.clear()
            im = ax.imshow(raster_exportion[frame], cmap='jet')
            ax.set_title(str(FileName).format(frame))
            ax.set_xlabel("x (m)")
            ax.set_ylabel("y (m)")
            fig.colorbar(im, cax=cax).set_label('Velocity (m/s)')

        cax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        # Create the animation object
        ani = FuncAnimation(fig, update, frames=raster_exportion.shape[0], interval=300)
        # Save the animation as a GIF
        ani.save(str(FileName) + '.gif', writer='pillow', dpi=600)
        plt.close()
    else: # for Human instability flooding risk
        fig, ax = plt.subplots()
        # Define custom colormap
        cmap = ListedColormap(['#00FFFFFF', '#52BE80','#229954','#145A32','#F4D03F','#D4AC0D','#9A7D0A','#E67E22','#CA6F1E','#935116','#C0392B','#A93226','#641E16',])
        boundaries = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5]  # Boundaries between discrete values
        norm = BoundaryNorm(boundaries, ncolors=cmap.N, clip=True)

        def update(frame):
            ax.clear()
            im = ax.imshow(raster_exportion[frame], cmap=cmap, norm=norm)  # Use custom colormap and norm
            ax.set_title(str(FileName).format(frame))
            ax.set_xlabel("x (m)")
            ax.set_ylabel("y (m)")
            fig.colorbar(im, cax=cax, ticks=[0, 1, 2, 3])
            cax.get_yaxis().set_major_locator(FixedLocator([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]))
            cax.set_yticklabels(['1', '2', '3', '4','5','6','7','8','9','10','11','12'])

        cax = fig.add_axes([0.83, 0.15, 0.02, 0.7])
        # Create the animation object
        ani = FuncAnimation(fig, update, frames=raster_exportion.shape[0], interval=300)
        # Save the animation as a GIF
        ani.save(str(FileName) + '.gif', writer='pillow', dpi=600)
        plt.close()
    # wgs84 = pyproj.Proj('EPSG:4326')
    # # Define the source and target coordinate systems
    # utm = pyproj.Proj('EPSG:31983')
    # # Define the dimensions and coordinates of the dataset
    # lat = [yllcorner + cellsize * i for i in range(ncols)]
    # lon = [xllcorner + cellsize * i for i in range(nrows)]
    # #
    # # # to wgs84
    # n,lat=pyproj.transform(utm, wgs84, np.zeros_like(lat), lat)
    # lon,n=pyproj.transform(utm, wgs84, np.zeros_like(lon), lon)
    # #
    # time = [i for i in range(len(raster_exportion))]
    # #
    # data = xr.Dataset({'data':(['time','lat',"lon"],raster_exportion)},
    #                   coords={'time':time,
    #                           'lon':lon,
    #                           'lat':lat}
    # )
    # data = data.sortby('lat',ascending=False)
    #
    # data.to_netcdf('test.nc')

    # ds = nc.Dataset('test.nc', 'w', format='NETCDF4')
    #
    # lat_dim = ds.createDimension('lat', ncols)  # latitude axis
    # lon_dim = ds.createDimension('lon', nrows)  # longitude axis
    # time_dim = ds.createDimension('time', len(raster_exportion))
    #
    # times = ds.createVariable('time', 'f4', ('time',))
    # lats = ds.createVariable('lat', 'f4', ('lat',))
    # lons = ds.createVariable('lon', 'f4', ('lon',))
    # value = ds.createVariable('value', 'f4', ('time', 'lat', 'lon',))
    # value.units = 'Unknown'
    #
    # lats[:] = lat
    # lons[:] = lon
    # times[:] = time
    # value[:] = raster_exportion

    # ds.close()
    return