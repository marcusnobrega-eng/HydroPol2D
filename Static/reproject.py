from pathlib import Path
import rasterio
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.crs import CRS


def reproject_rasters_in_folder(
    input_folder,
    output_folder_name="reprojected_7755",
    target_epsg=7755,
    raster_extensions=(".tif", ".tiff", ".img", ".vrt")
):

    input_folder = Path(input_folder)
    output_folder = input_folder / output_folder_name
    output_folder.mkdir(exist_ok=True)

    target_crs = CRS.from_epsg(target_epsg)

    raster_files = [
        f for f in input_folder.iterdir()
        if f.is_file() and f.suffix.lower() in raster_extensions
    ]

    print(f"Found {len(raster_files)} raster(s)\n")

    for raster_path in raster_files:
        output_path = output_folder / raster_path.name
        print(f"Processing: {raster_path.name}")

        try:
            with rasterio.open(raster_path) as src:
                if src.crs is None:
                    print(f"  Skipped (no CRS): {raster_path.name}")
                    continue

                transform, width, height = calculate_default_transform(
                    src.crs,
                    target_crs,
                    src.width,
                    src.height,
                    *src.bounds
                )

                kwargs = src.meta.copy()
                kwargs.update({
                    "crs": target_crs,
                    "transform": transform,
                    "width": width,
                    "height": height,
                    "compress": "lzw"  # keeps output size reasonable
                })

                with rasterio.open(output_path, "w", **kwargs) as dst:
                    for band in range(1, src.count + 1):
                        reproject(
                            source=rasterio.band(src, band),
                            destination=rasterio.band(dst, band),
                            src_transform=src.transform,
                            src_crs=src.crs,
                            dst_transform=transform,
                            dst_crs=target_crs,
                            resampling=Resampling.bilinear  # better for continuous data
                        )

                print(f"  Saved → {output_path}\n")

        except Exception as e:
            print(f"  Error: {e}\n")


if __name__ == "__main__":
    input_folder = r"/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/India_Full/Static"

    reproject_rasters_in_folder(
        input_folder=input_folder,
        output_folder_name="reprojected_7755",
        target_epsg=7755
    )