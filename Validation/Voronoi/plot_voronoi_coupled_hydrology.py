#!/usr/bin/env python3
"""Plot and export the variable-mesh V-tilted coupled hydrology states."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
import rasterio
from rasterio.transform import from_origin
from scipy import sparse


def remapper(path: Path):
    with netCDF4.Dataset(path) as ds:
        mesh = np.asarray(ds["overlap_mesh_index"][:], dtype=int)
        raster = np.asarray(ds["overlap_raster_index"][:], dtype=int)
        area = np.asarray(ds["overlap_area_m2"][:], dtype=float)
        raster_area = np.asarray(ds["raster_area_m2"][:], dtype=float)
        x = np.asarray(ds["x_edges"][:]); y = np.asarray(ds["y_edges"][:])
        shape = int(ds.raster_rows), int(ds.raster_cols)
    matrix = sparse.csr_matrix((area / raster_area[raster], (raster, mesh)), shape=(len(raster_area), mesh.max() + 1))
    return lambda values: np.asarray(matrix @ values).reshape(shape), x, y


def write_tiff(path: Path, values: np.ndarray, x: np.ndarray, y: np.ndarray) -> None:
    transform = from_origin(x[0], y[-1], x[1] - x[0], y[1] - y[0])
    with rasterio.open(path, "w", driver="GTiff", height=values.shape[0], width=values.shape[1],
                       count=1, dtype="float64", transform=transform, compress="deflate", predictor=3) as dst:
        dst.write(np.flipud(values), 1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("results", type=Path)
    parser.add_argument("overlap", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    remap, x, y = remapper(args.overlap)
    with netCDF4.Dataset(args.results) as ds:
        native = {
            "peak-surface-depth-m": np.max(ds["surface_depth_m"][:], axis=0),
            "cumulative-infiltration-m": ds["cumulative_infiltration_m"][-1, :],
            "soil-water-m": ds["soil_water_m"][-1, :],
            "cumulative-actual-et-m": ds["cumulative_actual_et_m"][-1, :],
            "cumulative-recharge-m": ds["cumulative_recharge_m"][-1, :],
            "depth-to-groundwater-m": ds["depth_to_groundwater_m"][-1, :],
        }
    maps = {name: remap(np.asarray(value)) for name, value in native.items()}
    for name, values in maps.items(): write_tiff(args.output / f"vtilted-variable-{name}.tif", values, x, y)

    plt.rcParams.update({"font.family": ["Helvetica Neue", "Helvetica", "Arial", "sans-serif"],
                         "font.size": 9, "axes.linewidth": 1.5})
    panels = [
        ("peak-surface-depth-m", "Peak surface depth", "m", "Blues"),
        ("cumulative-infiltration-m", "Cumulative infiltration", "m", "YlGnBu"),
        ("soil-water-m", "Soil-water storage", "m", "YlGn"),
        ("cumulative-actual-et-m", "Cumulative actual ET", "m", "YlOrBr"),
        ("cumulative-recharge-m", "Cumulative recharge", "m", "PuBuGn"),
        ("depth-to-groundwater-m", "Depth to groundwater", "m", "cividis"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(11.2, 6.6), constrained_layout=True)
    for label, ax, (name, title, unit, cmap) in zip("abcdef", axes.flat, panels):
        image = ax.imshow(maps[name], origin="lower", extent=[x[0], x[-1], y[0], y[-1]], cmap=cmap, aspect="equal")
        ax.set_title(f"({label})  {title}", fontweight="normal", loc="left")
        ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
        colorbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.025)
        colorbar.set_label(unit); colorbar.outline.set_linewidth(1.5)
    figure = args.output / "vtilted-variable-coupled-hydrology"
    fig.savefig(figure.with_suffix(".png"), dpi=300)
    fig.savefig(figure.with_suffix(".pdf"))
    fig.savefig(figure.with_suffix(".svg"))
    plt.close(fig)


if __name__ == "__main__":
    main()
