#!/usr/bin/env python3
"""Create QGIS rasters and compact diagnostics for the Pune Voronoi event."""

from pathlib import Path

import geopandas as gpd
import matplotlib as mpl
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import pandas as pd
import rasterio
from matplotlib.colors import LinearSegmentedColormap
from rasterio.features import rasterize
from rasterio.transform import from_origin


HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "Outputs" / "100mm_1h_12h"
RESULT = OUTPUT / "pune-100mm-1h-12h-results.nc"
SUMMARY = OUTPUT / "pune-100mm-1h-12h-summary.csv"
DIAGNOSTICS = OUTPUT / "pune-100mm-1h-12h-diagnostics.csv"
MESH = (
    HERE.parents[2]
    / "HydroBathyDEM/examples/pune_catchment/outputs/"
    "pune_design_corridor_mesh_complete/mesh/hydropol_hybrid_mesh.gpkg"
)
RASTER_RESOLUTION_M = 30.0
NODATA = -9999.0


def variable_cell_time(dataset: netCDF4.Dataset, name: str) -> np.ndarray:
    variable = dataset.variables[name]
    values = np.asarray(variable[:], dtype=np.float32)
    dimensions = tuple(variable.dimensions)
    if dimensions == ("cell", "time"):
        return values
    if dimensions == ("time", "cell"):
        return values.T
    raise ValueError(f"Unexpected {name} dimensions: {dimensions}")


def write_raster(path: Path, values: np.ndarray, mesh: gpd.GeoDataFrame) -> np.ndarray:
    xmin, ymin, xmax, ymax = mesh.total_bounds
    width = int(np.ceil((xmax - xmin) / RASTER_RESOLUTION_M))
    height = int(np.ceil((ymax - ymin) / RASTER_RESOLUTION_M))
    transform = from_origin(xmin, ymax, RASTER_RESOLUTION_M, RASTER_RESOLUTION_M)
    image = rasterize(
        zip(mesh.geometry, values, strict=True),
        out_shape=(height, width),
        transform=transform,
        fill=NODATA,
        dtype="float32",
    )
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=height,
        width=width,
        count=1,
        dtype="float32",
        crs=mesh.crs,
        transform=transform,
        nodata=NODATA,
        compress="deflate",
        predictor=3,
        tiled=True,
    ) as destination:
        destination.write(image, 1)
    return np.ma.masked_equal(image, NODATA)


def positive_limit(values: np.ndarray, percentile: float = 99.5) -> float:
    finite = np.asarray(values)[np.isfinite(values) & (np.asarray(values) > 0)]
    return max(float(np.percentile(finite, percentile)), 1e-6) if finite.size else 1.0


def write_video(
    path: Path,
    values: np.ndarray,
    time_hours: np.ndarray,
    mesh: gpd.GeoDataFrame,
    title: str,
    units: str,
    cmap: mpl.colors.Colormap,
) -> None:
    xmin, ymin, xmax, ymax = mesh.total_bounds
    shape = (
        int(np.ceil((ymax - ymin) / RASTER_RESOLUTION_M)),
        int(np.ceil((xmax - xmin) / RASTER_RESOLUTION_M)),
    )
    transform = from_origin(xmin, ymax, RASTER_RESOLUTION_M, RASTER_RESOLUTION_M)
    vmax = positive_limit(values)
    figure, axis = plt.subplots(figsize=(9, 7), constrained_layout=True)
    first = rasterize(zip(mesh.geometry, values[:, 0], strict=True), out_shape=shape, transform=transform,
                      fill=NODATA, dtype="float32")
    plotted = axis.imshow(np.ma.masked_equal(first, NODATA), cmap=cmap, vmin=0, vmax=vmax)
    axis.set_axis_off()
    heading = axis.set_title(f"{title} — {time_hours[0]:.1f} h", fontweight="normal")
    colorbar = figure.colorbar(plotted, ax=axis, shrink=0.82, extend="max")
    colorbar.set_label(units)
    colorbar.outline.set_linewidth(1.5)

    def update(frame: int):
        image = rasterize(zip(mesh.geometry, values[:, frame], strict=True), out_shape=shape,
                          transform=transform, fill=NODATA, dtype="float32")
        plotted.set_data(np.ma.masked_equal(image, NODATA))
        heading.set_text(f"{title} — {time_hours[frame]:.1f} h")
        return plotted, heading

    movie = animation.FuncAnimation(figure, update, frames=values.shape[1], blit=False)
    movie.save(path, writer=animation.FFMpegWriter(fps=2, bitrate=4000))
    plt.close(figure)


def main() -> None:
    if not RESULT.exists():
        raise FileNotFoundError(f"Solver output is not ready: {RESULT}")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    mesh = gpd.read_file(MESH, layer="mesh", columns=["face_id", "geometry"])
    mesh = mesh.sort_values("face_id").reset_index(drop=True)

    with netCDF4.Dataset(RESULT) as dataset:
        depth = variable_cell_time(dataset, "surface_depth_m")
        velocity = variable_cell_time(dataset, "surface_velocity_m_s")
        time_hours = np.asarray(dataset.variables["time_s"][:], dtype=float) / 3600.0
    if len(mesh) != depth.shape[0]:
        raise ValueError(f"Mesh/result mismatch: {len(mesh)} polygons vs {depth.shape[0]} cells")

    fields = {
        "maximum_depth_m": np.nanmax(depth, axis=1),
        "final_depth_m": depth[:, -1],
        "maximum_velocity_m_s": np.nanmax(velocity, axis=1),
        "final_velocity_m_s": velocity[:, -1],
    }
    images = {
        name: write_raster(OUTPUT / f"pune-100mm-1h-12h-{name}.tif", values, mesh)
        for name, values in fields.items()
    }

    mpl.rcParams.update(
        {
            "font.family": ["Helvetica Neue", "Helvetica", "Arial", "sans-serif"],
            "font.size": 9,
            "axes.linewidth": 1.5,
            "xtick.major.width": 1.5,
            "ytick.major.width": 1.5,
        }
    )
    depth_cmap = LinearSegmentedColormap.from_list(
        "depth", ["#F7FBFF", "#8CC5E3", "#1A80BB", "#082A54"]
    )
    velocity_cmap = LinearSegmentedColormap.from_list(
        "velocity", ["#FFF7EC", "#F0C571", "#EA801C", "#9D2C00"]
    )
    write_video(OUTPUT / "pune-100mm-1h-12h-depth.mp4", depth, time_hours, mesh,
                "Surface-water depth", "Depth (m)", depth_cmap)
    write_video(OUTPUT / "pune-100mm-1h-12h-velocity.mp4", velocity, time_hours, mesh,
                "Surface-water velocity", "Velocity (m s$^{-1}$)", velocity_cmap)
    panels = [
        ("maximum_depth_m", "Maximum depth", "Depth (m)", depth_cmap),
        ("final_depth_m", "Depth after 12 h", "Depth (m)", depth_cmap),
        ("maximum_velocity_m_s", "Maximum velocity", "Velocity (m s$^{-1}$)", velocity_cmap),
        ("final_velocity_m_s", "Velocity after 12 h", "Velocity (m s$^{-1}$)", velocity_cmap),
    ]
    figure, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for axis, (name, title, label, cmap) in zip(axes.flat, panels, strict=True):
        limit = positive_limit(fields[name])
        plotted = axis.imshow(images[name], cmap=cmap, vmin=0, vmax=limit)
        axis.set_title(title, fontweight="normal")
        axis.set_axis_off()
        colorbar = figure.colorbar(plotted, ax=axis, shrink=0.82, extend="max")
        colorbar.set_label(label)
        colorbar.outline.set_linewidth(1.5)
        axis.text(
            0.01,
            0.01,
            f"True maximum: {np.nanmax(fields[name]):.3g}",
            transform=axis.transAxes,
            ha="left",
            va="bottom",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 2},
        )
    figure.savefig(OUTPUT / "pune-100mm-1h-12h-state-maps.png", dpi=300)
    figure.savefig(OUTPUT / "pune-100mm-1h-12h-state-maps.pdf")
    plt.close(figure)

    diagnostics = pd.read_csv(DIAGNOSTICS)
    diagnostic_time_h = diagnostics["time_s"].to_numpy() / 3600.0
    dt = diagnostics["dt_s"].to_numpy()
    boundary_volume = diagnostics["boundary_net_inflow_volume_m3"].to_numpy()
    outlet_discharge = np.maximum(-boundary_volume, 0) / dt
    source_volume = diagnostics["surface_source_volume_m3"].to_numpy()
    cumulative_source = np.cumsum(source_volume)
    cumulative_outlet = np.cumsum(np.maximum(-boundary_volume, 0))
    residual = diagnostics["step_mass_residual_m3"].to_numpy()
    residual_fraction = residual / max(cumulative_source[-1], 1.0)

    figure, axes = plt.subplots(2, 2, figsize=(10, 6.5), constrained_layout=True)
    axes[0, 0].plot(diagnostic_time_h, diagnostics["max_surface_depth_m"], color="#1A80BB", lw=2)
    axes[0, 0].set(ylabel="Maximum depth (m)")
    axes[0, 1].plot(diagnostic_time_h, diagnostics["max_surface_velocity_m_s"], color="#EA801C", lw=2)
    axes[0, 1].set(ylabel="Maximum velocity (m s$^{-1}$)")
    axes[1, 0].plot(diagnostic_time_h, dt, color="#298C8C", lw=1.5)
    axes[1, 0].set(xlabel="Simulation time (h)", ylabel="Adaptive timestep (s)")
    axes[1, 1].plot(diagnostic_time_h, outlet_discharge, color="#7E4794", lw=1.5, label="Outlet discharge")
    axes[1, 1].set(xlabel="Simulation time (h)", ylabel="Discharge (m$^3$ s$^{-1}$)")
    for axis in axes.flat:
        axis.grid(color="#D4D4D4", linewidth=0.6)
        axis.set_xlim(0, max(12, diagnostic_time_h[-1]))
    figure.savefig(OUTPUT / "pune-100mm-1h-12h-diagnostics.png", dpi=300)
    figure.savefig(OUTPUT / "pune-100mm-1h-12h-diagnostics.pdf")
    plt.close(figure)

    summary = pd.read_csv(SUMMARY).iloc[0].to_dict()
    summary.update(
        {
            "final_stored_surface_volume_m3": float(np.sum(depth[:, -1] * mesh.geometry.area.to_numpy())),
            "cumulative_outlet_volume_m3": float(cumulative_outlet[-1]),
            "maximum_absolute_mass_residual_fraction": float(np.max(np.abs(residual_fraction))),
            "map_display_percentile": 99.5,
            "map_rasterization_note": "30 m visualization raster; native polygon NetCDF is authoritative",
            "number_of_saved_times": int(depth.shape[1]),
            "last_saved_time_h": float(time_hours[-1]),
        }
    )
    pd.DataFrame([summary]).to_csv(OUTPUT / "pune-100mm-1h-12h-postprocessing-summary.csv", index=False)
    print(f"Postprocessing complete: {OUTPUT}")


if __name__ == "__main__":
    main()
