"""Compare cell-mean and terrain-aware inundation for one Pune output time."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import rasterio
import xarray as xr
from matplotlib.colors import ListedColormap


ROOT = Path(__file__).resolve().parents[2]
HYDROBATHY = ROOT.parent / "HydroBathyDEM"
CASE = ROOT / "Applications" / "Pune_Voronoi_Storm"
RESULTS = CASE / "Outputs" / "100mm_1h_12h"
VISIBLE = CASE / "Outputs" / "postprocess_acceptance"
MESH_DIR = (
    HYDROBATHY
    / "examples/pune_catchment/outputs/pune_design_corridor_mesh_complete/mesh"
)
DEM = (
    HYDROBATHY
    / "examples/pune_catchment/outputs/hydrobathy_20km2_corrected_dem/dem"
    / "DEM_hydraulic_conditioned.tif"
)
TIME_S = 25_200.0
THRESHOLD_M = 0.01


def main() -> None:
    archive = RESULTS / "pune-100mm-1h-12h-results.nc"
    mesh_file = MESH_DIR / "hydropol_hybrid_mesh.nc"
    overlap_file = MESH_DIR / "hydropol_mesh_overlap.nc"
    current_file = (
        VISIBLE
        / "Rasters_Water_Depths"
        / "Flood_Depth_0000025200s.tif"
    )
    output_file = (
        VISIBLE
        / "Rasters_Water_Depths"
        / "Flood_Depth_Terrain_VolumeConservative_0000025200s.tif"
    )
    wse_output_file = (
        VISIBLE
        / "Rasters_Water_Depths"
        / "Flood_Depth_Terrain_HorizontalWSE_0000025200s.tif"
    )

    with xr.open_dataset(archive, decode_cf=False) as dataset:
        times = np.asarray(dataset["time_s"], dtype=np.float64)
        time_index = int(np.argmin(np.abs(times - TIME_S)))
        assert abs(times[time_index] - TIME_S) < 1e-6
        face_depth = np.asarray(
            dataset["surface_depth_m"].isel(time=time_index), dtype=np.float64
        )
    with xr.open_dataset(mesh_file, decode_cf=False) as dataset:
        face_area = np.asarray(dataset["cell_area_m2"], dtype=np.float64)
        face_bed = np.asarray(dataset["cell_bed_elevation_m"], dtype=np.float64)
    with xr.open_dataset(overlap_file, decode_cf=False) as dataset:
        face_index = np.asarray(dataset["overlap_mesh_index"], dtype=np.int64)
        pixel_index = np.asarray(dataset["overlap_raster_index"], dtype=np.int64)
        overlap_area = np.asarray(dataset["overlap_area_m2"], dtype=np.float64)
        rows = int(dataset.attrs["raster_rows"])
        columns = int(dataset.attrs["raster_cols"])

    with rasterio.open(DEM) as source:
        terrain_north = source.read(1).astype(np.float64)
        profile = source.profile.copy()
        nodata = source.nodata
    assert terrain_north.shape == (rows, columns)
    terrain = terrain_north[::-1].reshape(-1)
    if nodata is not None:
        assert not np.any(terrain[pixel_index] == nodata)
    assert np.isfinite(terrain[pixel_index]).all()

    volume = face_area * face_depth
    minimum = np.full(face_area.size, np.inf)
    maximum = np.full(face_area.size, -np.inf)
    np.minimum.at(minimum, face_index, terrain[pixel_index])
    np.maximum.at(maximum, face_index, terrain[pixel_index])
    assert np.isfinite(minimum).all() and np.isfinite(maximum).all()

    lower = minimum.copy()
    upper = maximum + np.divide(
        volume, face_area, out=np.zeros_like(volume), where=face_area > 0
    )
    for _ in range(36):
        stage = 0.5 * (lower + upper)
        trial = np.bincount(
            face_index,
            weights=overlap_area * np.maximum(stage[face_index] - terrain[pixel_index], 0),
            minlength=face_area.size,
        )
        low = trial < volume
        lower[low] = stage[low]
        upper[~low] = stage[~low]
    stage = 0.5 * (lower + upper)

    contribution = overlap_area * np.maximum(
        stage[face_index] - terrain[pixel_index], 0
    )
    pixel_volume = np.bincount(
        pixel_index, weights=contribution, minlength=rows * columns
    )
    covered_area = np.bincount(
        pixel_index, weights=overlap_area, minlength=rows * columns
    )
    pixel_area = abs(profile["transform"].a * profile["transform"].e)
    terrain_depth = (pixel_volume / pixel_area).reshape(rows, columns)[::-1]
    covered = (covered_area > 0).reshape(rows, columns)[::-1]

    model_wse = face_bed + face_depth
    wse_contribution = overlap_area * np.maximum(
        model_wse[face_index] - terrain[pixel_index], 0
    )
    wse_pixel_volume = np.bincount(
        pixel_index, weights=wse_contribution, minlength=rows * columns
    )
    wse_depth = (wse_pixel_volume / pixel_area).reshape(rows, columns)[::-1]

    profile.update(dtype="float32", count=1, nodata=-9999.0, compress="lzw")
    written = np.where(covered, terrain_depth, profile["nodata"]).astype(np.float32)
    with rasterio.open(output_file, "w", **profile) as target:
        target.write(written, 1)
    wse_written = np.where(covered, wse_depth, profile["nodata"]).astype(np.float32)
    with rasterio.open(wse_output_file, "w", **profile) as target:
        target.write(wse_written, 1)

    with rasterio.open(current_file) as source:
        current = source.read(1).astype(np.float64)
    native_volume = float(volume.sum())
    terrain_volume = float(terrain_depth.sum() * pixel_area)
    wse_volume = float(wse_depth.sum() * pixel_area)
    current_volume = float(current.sum() * pixel_area)
    statistics = {
        "native_volume_m3": native_volume,
        "cell_mean_raster_volume_m3": current_volume,
        "terrain_raster_volume_m3": terrain_volume,
        "terrain_relative_volume_error": (terrain_volume - native_volume) / native_volume,
        "horizontal_wse_raster_volume_m3": wse_volume,
        "horizontal_wse_volume_fraction": wse_volume / native_volume,
        "cell_mean_area_gt_1cm_km2": float(np.sum(current > THRESHOLD_M) * pixel_area / 1e6),
        "terrain_area_gt_1cm_km2": float(np.sum(terrain_depth > THRESHOLD_M) * pixel_area / 1e6),
        "horizontal_wse_area_gt_1cm_km2": float(np.sum(wse_depth > THRESHOLD_M) * pixel_area / 1e6),
        "cell_mean_max_depth_m": float(current.max()),
        "terrain_max_depth_m": float(terrain_depth.max()),
        "horizontal_wse_max_depth_m": float(wse_depth.max()),
    }
    _write_statistics(VISIBLE / "Tables_CSV" / "Terrain_Aware_Comparison_25200s.csv", statistics)
    _plot_comparison(current, terrain_depth, covered, VISIBLE)
    _plot_three_methods(current, wse_depth, terrain_depth, VISIBLE)
    print(output_file)
    print(wse_output_file)
    for key, value in statistics.items():
        print(f"{key}={value:.12g}")


def _write_statistics(path: Path, values: dict[str, float]) -> None:
    path.write_text(
        "Metric,Value\n" + "".join(f"{key},{value:.15g}\n" for key, value in values.items()),
        encoding="utf-8",
    )


def _plot_comparison(
    current: np.ndarray, terrain: np.ndarray, covered: np.ndarray, root: Path
) -> None:
    wet_current = current > THRESHOLD_M
    wet_terrain = terrain > THRESHOLD_M
    positive = np.concatenate((current[wet_current], terrain[wet_terrain]))
    vmax = float(np.quantile(positive, 0.995))
    difference = terrain - current
    dmax = float(np.quantile(np.abs(difference[covered]), 0.995))
    classes = np.zeros(current.shape, dtype=np.uint8)
    classes[wet_current & ~wet_terrain] = 1
    classes[wet_current & wet_terrain] = 2
    classes[~wet_current & wet_terrain] = 3

    plt.rcParams.update({"font.family": "Helvetica", "font.size": 10})
    figure, axes = plt.subplots(2, 2, figsize=(13, 11), constrained_layout=True)
    maps = [
        (np.ma.masked_where(~wet_current, current), "Cell-mean depth", "viridis", THRESHOLD_M, vmax),
        (np.ma.masked_where(~wet_terrain, terrain), "Terrain-aware depth", "viridis", THRESHOLD_M, vmax),
        (np.ma.masked_where(~covered, difference), "Terrain-aware minus cell mean", "RdBu", -dmax, dmax),
    ]
    for axis, (values, title, cmap, vmin, upper) in zip(axes.flat[:3], maps, strict=True):
        image = axis.imshow(values, cmap=cmap, vmin=vmin, vmax=upper)
        axis.set_title(title, fontweight="normal")
        axis.set_axis_off()
        figure.colorbar(image, ax=axis, shrink=0.82, label="Depth (m)")
    axes[1, 1].imshow(
        np.ma.masked_where(classes == 0, classes),
        cmap=ListedColormap(["#E25759", "#C8C8C8", "#0B81A2"]),
        vmin=1,
        vmax=3,
    )
    axes[1, 1].set_title("Wet extent at 1 cm", fontweight="normal")
    axes[1, 1].set_axis_off()
    axes[1, 1].text(
        0.02,
        0.02,
        "Red: cell mean only   Gray: both   Blue: terrain-aware only",
        transform=axes[1, 1].transAxes,
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "none"},
    )
    for folder, suffix in (("Figures_PNG", "png"), ("Figures_PDF", "pdf"), ("Figures_SVG", "svg")):
        figure.savefig(root / folder / f"Terrain_Aware_Comparison_25200s.{suffix}", dpi=250)
    plt.close(figure)


def _plot_three_methods(
    current: np.ndarray, horizontal_wse: np.ndarray, conservative: np.ndarray, root: Path
) -> None:
    arrays = (current, horizontal_wse, conservative)
    positive = np.concatenate(tuple(values[values > THRESHOLD_M] for values in arrays))
    vmax = float(np.quantile(positive, 0.995))
    figure, axes = plt.subplots(1, 3, figsize=(16, 6), constrained_layout=True)
    titles = (
        "Cell-mean depth",
        "Terrain-aware: model WSE",
        "Terrain-aware: volume conservative",
    )
    for axis, values, title in zip(axes, arrays, titles, strict=True):
        image = axis.imshow(
            np.ma.masked_less_equal(values, THRESHOLD_M),
            cmap="viridis",
            vmin=THRESHOLD_M,
            vmax=vmax,
        )
        axis.set_title(title, fontweight="normal")
        axis.set_axis_off()
        figure.colorbar(image, ax=axis, shrink=0.75, label="Depth (m)")
    for folder, suffix in (
        ("Figures_PNG", "png"),
        ("Figures_PDF", "pdf"),
        ("Figures_SVG", "svg"),
    ):
        figure.savefig(root / folder / f"Terrain_Aware_Three_Methods_25200s.{suffix}", dpi=250)
    plt.close(figure)


if __name__ == "__main__":
    main()
