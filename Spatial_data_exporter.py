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
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
# from matplotlib import rc
# rc('font', **{'family': 'serif', 'serif': ['Garamond']}, size = '13')
# rc('text', usetex=False)
import matplotlib.font_manager as fm
from matplotlib.ticker import FixedLocator
from matplotlib.colors import ListedColormap, BoundaryNorm
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation
import geopandas as gpd
import contextily as ctx
import time



def Spatial_data_exporter(FileName ,ncols, nrows, xllcorner, yllcorner, cellsize, no_data_value, raster_exportion, rain, rain_time):
    # Create a figure and axis object
    # Define the function to update the animation
    if FileName == 'Water Flood Depth':
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(6,8), gridspec_kw={'height_ratios': [0.10, 0.90]})
        shape = gpd.read_file('E:/Google_drive/Meu Drive/Papers/Paper - 2Dmodeling + human risk/runs/Example/591_basin.shp').to_crs(epsg=3857)
        bounds = shape.total_bounds
        extent = [bounds[0], bounds[2], bounds[1], bounds[3]]
        plt.subplots_adjust(hspace=-0.20)
        def update(frame):
            ax1.clear()
            ax2.clear()
            #  rainfall plot with
            ax1.plot(range(0, len(rain)*rain_time,rain_time), rain, zorder=1)
            ax1.scatter(frame*rain_time, rain[frame], s=50, color='tab:green', edgecolor='black', zorder=2)
            ax1.invert_yaxis()
            ax1.xaxis.set_ticks_position('top')
            ax1.xaxis.set_label_position('top')
            ax1.set_xlabel('Time (min)')
            ax1.set_ylabel('Rainfall \n Intensity (mm/h)')

            shape.boundary.plot(ax=ax2, color='black', linewidth=1)
            while True:
                try:
                    ctx.add_basemap(ax=ax2, zoom=14, alpha=0.3)
                    break
                except ConnectionError:
                    print("Connection timed out. Retrying...")
                    time.sleep(120)  # wait for 5 seconds before trying again
                except Exception as e:
                    print(f"An exception ocurred: {str(e)}")
                    break
            im = ax2.imshow(raster_exportion[frame], extent=extent, cmap='jet')
            ax2.set_xlabel("Longitude (m)")
            ax2.set_ylabel("Latitude (m)")
            ax2.ticklabel_format(useOffset=False, style='plain')
            ax2.tick_params(axis='y', labelrotation=90)

            fig.colorbar(im, cax=cax, orientation='horizontal').set_label('Flow Depth (m)')

        cax = fig.add_axes([0.15, 0.12, 0.70, 0.03])
        # Create the animation object

        ani = FuncAnimation(fig, update, frames=raster_exportion.shape[0], interval=1000)
        # Save the animation as a GIF
        ani.save(str(FileName)+'.gif', writer='pillow', dpi=600)
        plt.close()
    elif FileName == 'Water Surface Elevation':
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(6,8), gridspec_kw={'height_ratios': [0.10, 0.90]})
        shape = gpd.read_file('E:/Google_drive/Meu Drive/Papers/Paper - 2Dmodeling + human risk/runs/Example/591_basin.shp').to_crs(epsg=3857)
        bounds = shape.total_bounds
        extent = [bounds[0], bounds[2], bounds[1], bounds[3]]
        plt.subplots_adjust(hspace=-0.20)

        def update(frame):
            ax1.clear()
            ax2.clear()
            #  rainfall plot with
            ax1.plot(range(0, len(rain)*rain_time,rain_time), rain, zorder=1)
            ax1.scatter(frame*rain_time, rain[frame], s=50, color='tab:green', edgecolor='black', zorder=2)
            ax1.invert_yaxis()
            ax1.xaxis.set_ticks_position('top')
            ax1.xaxis.set_label_position('top')
            ax1.set_xlabel('Time (min)')
            ax1.set_ylabel('Rainfall \n Intensity (mm/h)')

            shape.boundary.plot(ax=ax2, color='black', linewidth=1)
            while True:
                try:
                    ctx.add_basemap(ax=ax2, zoom=14, alpha=0.3)
                    break
                except ConnectionError:
                    print("Connection timed out. Retrying...")
                    time.sleep(120)  # wait for 5 seconds before trying again
                except Exception as e:
                    print(f"An exception ocurred: {str(e)}")
                    break
            im = ax2.imshow(raster_exportion[frame], extent=extent, cmap='jet')
            fig.colorbar(im, cax=cax, orientation='horizontal').set_label('Flow elevations (masl)')
            ax2.set_xlabel("Longitude (m)")
            ax2.set_ylabel("Latitude (m)")
            ax2.ticklabel_format(useOffset=False, style='plain')
            ax2.tick_params(axis='y', labelrotation=90)

        cax = fig.add_axes([0.15, 0.12, 0.70, 0.03])
        # Create the animation object
        ani = FuncAnimation(fig, update, frames=raster_exportion.shape[0], interval=1000)
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
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(6,8), gridspec_kw={'height_ratios': [0.10, 0.90]})
        shape = gpd.read_file('E:/Google_drive/Meu Drive/Papers/Paper - 2Dmodeling + human risk/runs/Example/591_basin.shp').to_crs(epsg=3857)
        bounds = shape.total_bounds
        extent = [bounds[0], bounds[2], bounds[1], bounds[3]]
        plt.subplots_adjust(hspace=-0.20)
        def update(frame):
            ax1.clear()
            ax2.clear()
            #  rainfall plot with
            ax1.plot(range(0, len(rain)*rain_time,rain_time), rain, zorder=1)
            ax1.scatter(frame*rain_time, rain[frame], s=50, color='tab:green', edgecolor='black', zorder=2)
            ax1.invert_yaxis()
            ax1.xaxis.set_ticks_position('top')
            ax1.xaxis.set_label_position('top')
            ax1.set_xlabel('Time (min)')
            ax1.set_ylabel('Rainfall \n Intensity (mm/h)')

            shape.boundary.plot(ax=ax2, color='black', linewidth=1)
            while True:
                try:
                    ctx.add_basemap(ax=ax2, zoom=14, alpha=0.3)
                    break
                except ConnectionError:
                    print("Connection timed out. Retrying...")
                    time.sleep(120)  # wait for 5 seconds before trying again
                except Exception as e:
                    print(f"An exception ocurred: {str(e)}")
                    break
            im = ax2.imshow(raster_exportion[frame], extent=extent, cmap='jet')
            ax2.set_xlabel("Longitude (m)")
            ax2.set_ylabel("Latitude (m)")
            ax2.ticklabel_format(useOffset=False, style='plain')
            ax2.tick_params(axis='y', labelrotation=90)

            fig.colorbar(im, cax=cax, orientation='horizontal').set_label('Velocity (m/s)')

        cax = fig.add_axes([0.15, 0.12, 0.70, 0.03])
        # Create the animation object
        ani = FuncAnimation(fig, update, frames=raster_exportion.shape[0], interval=1000)
        # Save the animation as a GIF
        ani.save(str(FileName) + '.gif', writer='pillow', dpi=600)
        plt.close()
    elif (FileName == 'FTI Risk') or (FileName == 'BTI Risk'):
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(6,8), gridspec_kw={'height_ratios': [0.10, 0.90]})
        # Define custom colormap
        shape = gpd.read_file('E:/Google_drive/Meu Drive/Papers/Paper - 2Dmodeling + human risk/runs/Example/591_basin.shp').to_crs(epsg=3857)
        bounds = shape.total_bounds
        extent = [bounds[0], bounds[2], bounds[1], bounds[3]]
        # Remove space between subplots
        plt.subplots_adjust(hspace=-0.20)
        # cmap = ListedColormap(['#ffffff00', '#52BE80','#229954','#145A32','#F4D03F','#D4AC0D','#9A7D0A','#E67E22','#CA6F1E','#935116','#C0392B','#A93226','#641E16',])
        cmap = ListedColormap(
            ['#ffffff00', '#E67E22', '#F4D03F', '#C0392B'])
        boundaries = [-0.5, 0.5, 1.5, 2.5, 3.5]  # Boundaries between discrete values
        norm = BoundaryNorm(boundaries, ncolors=cmap.N, clip=True)

        def update(frame):
            ax1.clear()
            ax2.clear()
            #  rainfall plot with
            ax1.plot(range(0, len(rain)*rain_time,rain_time), rain, zorder=1)
            ax1.scatter(frame*rain_time, rain[frame], s=50, color='tab:green', edgecolor='black', zorder=2)
            ax1.invert_yaxis()
            ax1.xaxis.set_ticks_position('top')
            ax1.xaxis.set_label_position('top')
            ax1.set_xlabel('Time (min)')
            ax1.set_ylabel('Rainfall \n Intensity (mm/h)')
            #  risk plot
            shape.boundary.plot(ax=ax2, color='black', linewidth=1)
            # try to add basemap, if it fails due to a ConnectTimeout exception, wait for 5 seconds and try again
            while True:
                try:
                    ctx.add_basemap(ax=ax2, zoom=14, alpha=0.3)
                    break
                except ConnectionError:
                    print("Connection timed out. Retrying...")
                    time.sleep(120)  # wait for 5 seconds before trying again
                except Exception as e:
                    print(f"An exception ocurred: {str(e)}")
                    break

            im = ax2.imshow(raster_exportion[frame], extent=extent, cmap=cmap, norm=norm)  # Use custom colormap and norm
            ax2.set_xlabel("Longitude (m)")
            ax2.set_ylabel("Latitude (m)")
            ax2.ticklabel_format(useOffset=False, style='plain')
            ax2.tick_params(axis='y', labelrotation=90)

            fig.colorbar(im, cax=cax, orientation='horizontal')
            cax.get_xaxis().set_major_locator(FixedLocator([0, 1, 2, 3]))
            cax.yaxis.set_major_locator(plt.NullLocator())
            cax.set_xticklabels(['No \n Riks', 'Sliding \n Risk', 'Sliding + Toppling \n Risk', 'Toppling \n Risk'])

        cax = fig.add_axes([0.15, 0.12, 0.70, 0.03])
        # fig.tight_layout()
        # Create the animation object
        ani = FuncAnimation(fig, update, frames=raster_exportion.shape[0], interval=1000)
        # Save the animation as a GIF
        ani.save(str(FileName) + '.gif', writer='pillow', dpi=600)
        plt.close()

    else: # for Human instability flooding risk
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(6,8), gridspec_kw={'height_ratios': [0.10, 0.90]})
        # Define custom colormap
        shape = gpd.read_file('E:/Google_drive/Meu Drive/Papers/Paper - 2Dmodeling + human risk/runs/Example/591_basin.shp').to_crs(epsg=3857)
        bounds = shape.total_bounds
        extent = [bounds[0], bounds[2], bounds[1], bounds[3]]
        # Remove space between subplots
        plt.subplots_adjust(hspace=-0.20)
        # cmap = ListedColormap(['#ffffff00', '#52BE80','#229954','#145A32','#F4D03F','#D4AC0D','#9A7D0A','#E67E22','#CA6F1E','#935116','#C0392B','#A93226','#641E16',])
        cmap = ListedColormap(
            ['#ffffff00', '#52BE80', '#52BE80', '#52BE80', '#F4D03F', '#F4D03F', '#F4D03F', '#E67E22', '#E67E22',
             '#E67E22', '#C0392B', '#C0392B', '#C0392B'])
        boundaries = [-0.5, 0.5, 3.5, 6.5, 9.5, 12.5]  # Boundaries between discrete values
        norm = BoundaryNorm(boundaries, ncolors=cmap.N, clip=True)

        def update(frame):
            ax1.clear()
            ax2.clear()
            #  rainfall plot with
            ax1.plot(range(0, len(rain)*rain_time,rain_time), rain, zorder=1)
            ax1.scatter(frame*rain_time, rain[frame], s=50, color='tab:green', edgecolor='black', zorder=2)
            ax1.invert_yaxis()
            ax1.xaxis.set_ticks_position('top')
            ax1.xaxis.set_label_position('top')
            ax1.set_xlabel('Time (min)')
            ax1.set_ylabel('Rainfall \n Intensity (mm/h)')
            #  risk plot
            shape.boundary.plot(ax=ax2, color='black', linewidth=1)
            # try to add basemap, if it fails due to a ConnectTimeout exception, wait for 5 seconds and try again
            while True:
                try:
                    ctx.add_basemap(ax=ax2, zoom=14, alpha=0.3)
                    break
                except ConnectionError:
                    print("Connection timed out. Retrying...")
                    time.sleep(120)  # wait for 5 seconds before trying again
                except Exception as e:
                    print(f"An exception ocurred: {str(e)}")
                    break

            im = ax2.imshow(raster_exportion[frame], extent=extent, cmap=cmap, norm=norm)  # Use custom colormap and norm
            ax2.set_xlabel("Longitude (m)")
            ax2.set_ylabel("Latitude (m)")
            ax2.ticklabel_format(useOffset=False, style='plain')
            ax2.tick_params(axis='y', labelrotation=90)

            fig.colorbar(im, cax=cax, orientation='horizontal')
            cax.get_xaxis().set_major_locator(FixedLocator([0, 2, 5, 8, 11]))
            cax.yaxis.set_major_locator(plt.NullLocator())
            cax.set_xticklabels(['No \n Riks', 'Child \n Risk', 'Teen \n Risk', 'Elderly \n Risk', 'Adult \n Risk'])

        cax = fig.add_axes([0.15, 0.12, 0.70, 0.03])
        # fig.tight_layout()
        # Create the animation object
        ani = FuncAnimation(fig, update, frames=raster_exportion.shape[0], interval=1000)
        # Save the animation as a GIF
        ani.save(str(FileName) + '.gif', writer='pillow', dpi=600)
        plt.close()


    return