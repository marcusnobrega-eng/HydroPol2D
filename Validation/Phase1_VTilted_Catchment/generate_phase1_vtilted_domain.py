#!/usr/bin/env python3
"""Generate the Phase 1 v-tilted synthetic catchment rasters.

This is the Python equivalent of `generate_phase1_vtilted_domain.m`, provided
so the domain can be generated in environments without MATLAB.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin


def read_config(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["parameter"]: row["value"] for row in csv.DictReader(handle)}


def f(config: dict[str, str], name: str) -> float:
    return float(config[name])


def write_raster(path: Path, data: np.ndarray, transform, epsg: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(
        path,
        "w",
        driver="GTiff",
        height=data.shape[0],
        width=data.shape[1],
        count=1,
        dtype="float32",
        crs=f"EPSG:{epsg}",
        transform=transform,
        nodata=None,
        compress="lzw",
    ) as dst:
        dst.write(data.astype("float32"), 1)


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    case_dir = Path(__file__).resolve().parent
    static_dir = case_dir / "Static"
    output_dir = case_dir / "Outputs" / "Validation"
    config = read_config(case_dir / "Config" / "Domain_Config.csv")

    dx = f(config, "dx")
    length_y = f(config, "length_y")
    left_width = f(config, "left_width")
    channel_width = f(config, "channel_width")
    right_width = f(config, "right_width")
    left_slope = f(config, "left_slope")
    right_slope = f(config, "right_slope")
    long_slope = f(config, "longitudinal_slope")
    z_outlet = f(config, "z_outlet")
    channel_incision = f(config, "channel_incision")
    epsg = int(f(config, "epsg"))

    width_total = left_width + channel_width + right_width
    nx = int(round(width_total / dx))
    ny = int(round(length_y / dx))

    col = np.arange(nx, dtype=float)
    row = np.arange(ny, dtype=float)
    x = col[None, :] * dx
    y = row[:, None] * dx

    x_left_end = left_width
    x_channel_end = left_width + channel_width

    left_mask = np.broadcast_to(x < x_left_end, (ny, nx))
    channel_mask = np.broadcast_to((x >= x_left_end) & (x <= x_channel_end), (ny, nx))
    right_mask = np.broadcast_to(x > x_channel_end, (ny, nx))

    dem = z_outlet + long_slope * np.broadcast_to(y, (ny, nx))
    channel_base = dem - channel_incision
    dem[left_mask] = channel_base[left_mask] + left_slope * (x_left_end - np.broadcast_to(x, (ny, nx))[left_mask])
    dem[channel_mask] = channel_base[channel_mask]
    dem[right_mask] = channel_base[right_mask] + right_slope * (np.broadcast_to(x, (ny, nx))[right_mask] - x_channel_end)

    lulc = np.zeros((ny, nx), dtype=float)
    soil = np.zeros((ny, nx), dtype=float)
    zone_id = np.zeros((ny, nx), dtype=float)

    for array in (lulc, soil, zone_id):
        array[left_mask] = 1
        array[channel_mask] = 2
        array[right_mask] = 3

    def zoned(left_name: str, channel_name: str, right_name: str) -> np.ndarray:
        arr = np.zeros((ny, nx), dtype=float)
        arr[left_mask] = f(config, left_name)
        arr[channel_mask] = f(config, channel_name)
        arr[right_mask] = f(config, right_name)
        return arr

    dtb = zoned("dtb_left", "dtb_channel", "dtb_right")
    albedo = zoned("albedo_left", "albedo_channel", "albedo_right")
    lai = zoned("lai_left", "lai_channel", "lai_right")
    initial_sm = zoned("initial_sm_left", "initial_sm_channel", "initial_sm_right")
    gw_depth = zoned("gw_depth_left", "gw_depth_channel", "gw_depth_right")
    gw_table = dem - dtb + gw_depth

    rasters = {
        "DEM.tif": dem,
        "LULC.tif": lulc,
        "SOIL.tif": soil,
        "DTB.tif": dtb,
        "Albedo.tif": albedo,
        "LAI.tif": lai,
        "Initial_SM.tif": initial_sm,
        "GW_depth.tif": gw_depth,
        "GW_table.tif": gw_table,
        "Zone_ID.tif": zone_id,
    }

    for name, arr in rasters.items():
        if not np.isfinite(arr).all():
            raise ValueError(f"{name} contains NaN or Inf")

    transform = from_origin(0.0, length_y, dx, dx)
    for name, arr in rasters.items():
        write_raster(static_dir / name, arr, transform, epsg)

    summary_rows = [
        {
            "zone": "left_hillslope",
            "cell_count": int(left_mask.sum()),
            "lai": f(config, "lai_left"),
            "initial_sm_mm": f(config, "initial_sm_left"),
            "dtb_m": f(config, "dtb_left"),
        },
        {
            "zone": "channel_strip",
            "cell_count": int(channel_mask.sum()),
            "lai": f(config, "lai_channel"),
            "initial_sm_mm": f(config, "initial_sm_channel"),
            "dtb_m": f(config, "dtb_channel"),
        },
        {
            "zone": "right_hillslope",
            "cell_count": int(right_mask.sum()),
            "lai": f(config, "lai_right"),
            "initial_sm_mm": f(config, "initial_sm_right"),
            "dtb_m": f(config, "dtb_right"),
        },
    ]
    write_csv(output_dir / "VTilted_Domain_Summary.csv", summary_rows)

    print(f"Generated {len(rasters)} rasters in {static_dir}")
    print(f"Wrote summary to {output_dir / 'VTilted_Domain_Summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
