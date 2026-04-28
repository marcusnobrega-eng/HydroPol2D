import os
import rasterio
from rasterio.mask import mask
import geopandas as gpd

# === INPUTS ===
raster_folder = r"/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/India_5000m_cpu/India_Full_5000m/India_Full/Static"
shapefile_path = r"/oak/stanford/groups/gorelick/Marcus/India/Shapefile/India_Country_Boundary.shp"

# === OUTPUT FOLDER ===
output_folder = os.path.join(raster_folder, "Clipped")
os.makedirs(output_folder, exist_ok=True)

# === LOAD SHAPEFILE ===
gdf = gpd.read_file(shapefile_path)

# Ensure geometry is in same CRS as rasters later
shapes = gdf.geometry.values

# === LOOP THROUGH RASTERS ===
for filename in os.listdir(raster_folder):
    if filename.lower().endswith((".tif", ".tiff")):
        raster_path = os.path.join(raster_folder, filename)
        
        with rasterio.open(raster_path) as src:
            # Reproject shapefile if CRS mismatch
            if gdf.crs != src.crs:
                gdf_proj = gdf.to_crs(src.crs)
                shapes = gdf_proj.geometry.values
            
            # Clip raster
            out_image, out_transform = mask(src, shapes, crop=True)
            
            # Update metadata
            out_meta = src.meta.copy()
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform
            })
            
            # Output path
            out_path = os.path.join(output_folder, f"clipped_{filename}")
            
            # Save clipped raster
            with rasterio.open(out_path, "w", **out_meta) as dest:
                dest.write(out_image)

        print(f"Clipped: {filename}")

print("✅ All rasters clipped successfully!")