# ============================================================
# HYDROPOL2D TEMP ANIMATION FRAME GENERATOR
# Uses contextily for automatic basemap handling
# Python 3.6 compatible
# ============================================================

from pathlib import Path
from datetime import datetime, timedelta
import re
import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import BoundaryNorm
from matplotlib import font_manager

import rasterio
from rasterio.warp import calculate_default_transform, reproject, transform_geom
from rasterio.enums import Resampling
from rasterio.transform import array_bounds
from rasterio import features

import contextily as cx

try:
    from scipy.io import loadmat as scipy_loadmat
except Exception:
    scipy_loadmat = None

try:
    import h5py
except Exception:
    h5py = None


# ============================================================
# USER INPUTS - PATHS
# ============================================================
TEMP_DIR = r"/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/Stanford/Stanford_Big_Flood_10m/Outputs/Temporary_Files"
DEM_TIF = r"/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/Stanford/Stanford_Big_Flood_10m/Static/DEM.tif"
OUTPUT_FRAMES_DIR = r"/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/Stanford/Stanford_Big_Flood_10m/Outputs/Png_Animations"
TEST_OUTPUT_PNG = r"/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/Stanford/Stanford_Big_Flood_10m/Outputs/Png_Animations/test_frame.png"


# ============================================================
# USER INPUTS - RUN MODE
# ============================================================
TEST_SINGLE_FRAME = False
TEST_FRAME_INDEX = 300


# ============================================================
# USER INPUTS - TIME SETTINGS
# ============================================================
DATE_BEGIN = datetime(2022,12,31,0,0,0);
MAP_SAVE_TIMESTEP_HOURS = 0.25
FLAG_ELAPSED_TIME = False

# ============================================================
# USER INPUTS - DEPTH / MASK SETTINGS
# ============================================================
DEPTH_THRESHOLD_M = 0.01
DISPLAY_THRESHOLD_M = 0.02


# ============================================================
# USER INPUTS - SCALE MODE
# ============================================================
SCALE_MODE = "manual"   # "manual" or "percentile"
MANUAL_DEPTH_MIN = 0.0
MANUAL_DEPTH_MAX = 2.0
PERCENTILE_FOR_MAX = 98.0
PERCENTILE_SAMPLE_COUNT = 12


# ============================================================
# USER INPUTS - FONT CUSTOMIZATION
# ============================================================
FONT_PATH = "/oak/stanford/groups/gorelick/Marcus/misc/Avenir/AvenirNextCyr-Regular.ttf"
FALLBACK_FONT_FAMILY = "DejaVu Sans"


# ============================================================
# USER INPUTS - TEXT CUSTOMIZATION
# ============================================================
TITLE_TEMPLATE = "Flood Depth\n{timestamp}"
CBAR_LABEL = "Flood depth (m)"
X_LABEL = "X [m]"
Y_LABEL = "Y [m]"
FRAME_LABEL_TEMPLATE = "Frame {current}/{total}"

TITLE_FONTSIZE = 22
AXIS_LABEL_FONTSIZE = 16
AXIS_TICK_FONTSIZE = 12
CBAR_LABEL_FONTSIZE = 16
CBAR_TICK_FONTSIZE = 12
FRAME_TEXT_FONTSIZE = 12


# ============================================================
# USER INPUTS - FIGURE / LAYOUT CUSTOMIZATION
# ============================================================
BASE_FIGSIZE = (10, 10)
FIG_WIDTH_SCALE = 2.0
FIG_HEIGHT_SCALE = 1.0
FIGSIZE = (BASE_FIGSIZE[0] * FIG_WIDTH_SCALE, BASE_FIGSIZE[1] * FIG_HEIGHT_SCALE)
DPI = 150

MAP_PAD_FRACTION_X = 0.02
MAP_PAD_FRACTION_Y = 0.02


# ============================================================
# USER INPUTS - FLOOD VISUALIZATION CUSTOMIZATION
# ============================================================
CMAP = "cool"
ALPHA = 0.75


# ============================================================
# USER INPUTS - COLORBAR CUSTOMIZATION
# ============================================================
COLORBAR_EXTEND = "max"
COLORBAR_FRACTION = 0.08
COLORBAR_SHRINK = 0.50
COLORBAR_PAD = 0.02
COLORBAR_EXTEND_FRAC = 0.08
COLORBAR_ASPECT = 20
N_COLOR_LEVELS = 10


# ============================================================
# USER INPUTS - CATCHMENT / OUTLINE CUSTOMIZATION
# ============================================================
OUTLINE_COLOR = "black"
OUTLINE_LINEWIDTH = 1.5
OUTLINE_ZORDER = 4


# ============================================================
# USER INPUTS - BASEMAP CUSTOMIZATION
# ============================================================
BASEMAP_SOURCE = cx.providers.CartoDB.DarkMatter
BASEMAP_ZOOM = "auto"
BASEMAP_ATTRIBUTION = False


# ============================================================
# USER INPUTS - EXPORT CUSTOMIZATION
# ============================================================
SAVE_BBOX_INCHES = "tight"
PNG_FACE_COLOR = "white"
PROGRESS_EVERY_N_FRAMES = 10


# ============================================================
# FONT SETUP
# ============================================================
def setup_font(font_path, fallback_family):
    if font_path is not None:
        font_path = Path(font_path)
        if not font_path.exists():
            raise FileNotFoundError("FONT_PATH does not exist: {}".format(font_path))
        font_manager.fontManager.addfont(str(font_path))
        font_name = font_manager.FontProperties(fname=str(font_path)).get_name()
        mpl.rcParams["font.family"] = font_name
        print("Using custom font: {}".format(font_name))
    else:
        mpl.rcParams["font.family"] = fallback_family
        print("Using fallback font: {}".format(fallback_family))


# ============================================================
# HELPERS - FILES / TIME
# ============================================================
def list_chunk_files(temp_dir):
    temp_dir = Path(temp_dir)
    if not temp_dir.exists():
        raise FileNotFoundError("TEMP_DIR does not exist: {}".format(temp_dir))

    files = list(temp_dir.glob("save_map_hydro_*.mat"))
    if not files:
        raise FileNotFoundError("No save_map_hydro_*.mat files found in {}".format(temp_dir))

    def extract_num(p):
        m = re.search(r"save_map_hydro_(\d+)\.mat$", p.name)
        if m is None:
            return math.inf
        return int(m.group(1))

    return sorted(files, key=extract_num)


def frame_time_from_index(frame_index):
    return DATE_BEGIN + timedelta(hours=frame_index * MAP_SAVE_TIMESTEP_HOURS)


# ============================================================
# HELPERS - MATLAB READING
# ============================================================
def _is_hdf5_mat(path):
    path = Path(path)
    with open(str(path), "rb") as f:
        header = f.read(1024)

    hdf_sig = b"\x89HDF\r\n\x1a\n"
    return (hdf_sig in header) or (b"MATLAB 7.3 MAT-file" in header)


def _safe_squeeze(arr):
    arr = np.asarray(arr)
    arr = np.squeeze(arr)
    return arr


def _convert_hdf5_node(node, root):
    if isinstance(node, h5py.Dataset):
        data = node[()]

        if isinstance(data, bytes):
            return data.decode("utf-8", errors="ignore")

        if getattr(data, "dtype", None) is not None and data.dtype.kind == "O":
            out = []
            for ref in data.flat:
                out.append(_convert_hdf5_node(root[ref], root))
            return np.array(out, dtype=object).reshape(data.shape)

        return np.array(data)

    if isinstance(node, h5py.Group):
        out = {}
        for k in node.keys():
            if k == "#refs#":
                continue
            out[k] = _convert_hdf5_node(node[k], root)
        return out

    return node


def load_mat_file(path):
    path = Path(path)

    if _is_hdf5_mat(path):
        if h5py is None:
            raise ImportError(
                "This MAT file is MATLAB v7.3/HDF5, but h5py is not installed. "
                "Install it with: pip install h5py"
            )
        with h5py.File(str(path), "r") as f:
            return _convert_hdf5_node(f, f)

    if scipy_loadmat is None:
        raise ImportError("scipy is required to read non-v7.3 .mat files.")

    try:
        data = scipy_loadmat(str(path), squeeze_me=True, struct_as_record=False)
        out = {}
        for k, v in data.items():
            if k.startswith("__"):
                continue
            out[k] = v
        return out

    except NotImplementedError as e:
        msg = str(e)
        if "Please use HDF reader for matlab v7.3 files" in msg:
            if h5py is None:
                raise ImportError(
                    "This MAT file is MATLAB v7.3/HDF5, but h5py is not installed. "
                    "Install it with: pip install h5py"
                )
            with h5py.File(str(path), "r") as f:
                return _convert_hdf5_node(f, f)
        raise


def maybe_getattr(obj, name, default=None):
    if isinstance(obj, dict):
        return obj.get(name, default)
    return getattr(obj, name, default)


def get_maps_struct(loaded):
    if "Maps" in loaded:
        return loaded["Maps"]
    for v in loaded.values():
        hydro = maybe_getattr(v, "Hydro", None)
        if hydro is not None:
            return v
    return None


def extract_depth_3d_from_maps(maps):
    hydro = maybe_getattr(maps, "Hydro", None)
    if hydro is None:
        raise KeyError("Could not find Maps.Hydro in MAT file")

    d = maybe_getattr(hydro, "d", None)
    if d is None:
        raise KeyError("Could not find Maps.Hydro.d in MAT file")

    d = _safe_squeeze(d)
    d = np.asarray(d)

    if d.ndim == 2:
        d = d[:, :, np.newaxis]
    elif d.ndim != 3:
        raise ValueError("Maps.Hydro.d must be 2D or 3D. Got shape {}".format(d.shape))

    d = d.astype(np.float32)

    # HydroPol recovery code indicates mm -> m
    depth_m = d / 1000.0
    return depth_m


def standardize_depth_cube_to_dem(depth_3d, dem_shape):
    """
    Expected output:
        (dem_rows, dem_cols, n_frames)

    Handles common MATLAB/HDF5 axis-order issues and the common HydroPol
    case where the saved map is missing a 1-cell border relative to the DEM.
    """
    dem_rows, dem_cols = dem_shape

    if depth_3d.ndim != 3:
        raise ValueError("depth_3d must be 3D. Got shape {}".format(depth_3d.shape))

    s = depth_3d.shape

    # Case 1: already correct (rows, cols, t)
    if s[0] == dem_rows and s[1] == dem_cols:
        return depth_3d

    # Case 2: row/col swapped (cols, rows, t)
    if s[0] == dem_cols and s[1] == dem_rows:
        print("Depth cube appears transposed in row/col. Fixing with transpose(1,0,2).")
        return np.transpose(depth_3d, (1, 0, 2))

    # Case 3: time first (t, rows, cols)
    if s[1] == dem_rows and s[2] == dem_cols:
        print("Depth cube appears to be (t, rows, cols). Fixing with transpose(1,2,0).")
        return np.transpose(depth_3d, (1, 2, 0))

    # Case 4: time first, row/col swapped (t, cols, rows)
    if s[1] == dem_cols and s[2] == dem_rows:
        print("Depth cube appears to be (t, cols, rows). Fixing with transpose(2,1,0).")
        return np.transpose(depth_3d, (2, 1, 0))

    # Case 5: (t, cols-2, rows-2)
    if s[1] == dem_cols - 2 and s[2] == dem_rows - 2:
        print("Depth cube appears to be (t, cols-2, rows-2).")
        print("Fixing with transpose(2,1,0) and 1-cell edge padding.")
        arr = np.transpose(depth_3d, (2, 1, 0))
        arr = np.pad(
            arr,
            pad_width=((1, 1), (1, 1), (0, 0)),
            mode="constant",
            constant_values=np.nan
        )
        return arr

    # Case 6: (t, rows-2, cols-2)
    if s[1] == dem_rows - 2 and s[2] == dem_cols - 2:
        print("Depth cube appears to be (t, rows-2, cols-2).")
        print("Fixing with transpose(1,2,0) and 1-cell edge padding.")
        arr = np.transpose(depth_3d, (1, 2, 0))
        arr = np.pad(
            arr,
            pad_width=((1, 1), (1, 1), (0, 0)),
            mode="constant",
            constant_values=np.nan
        )
        return arr

    # Case 7: (cols-2, rows-2, t)
    if s[0] == dem_cols - 2 and s[1] == dem_rows - 2:
        print("Depth cube appears to be (cols-2, rows-2, t).")
        print("Fixing with transpose(1,0,2) and 1-cell edge padding.")
        arr = np.transpose(depth_3d, (1, 0, 2))
        arr = np.pad(
            arr,
            pad_width=((1, 1), (1, 1), (0, 0)),
            mode="constant",
            constant_values=np.nan
        )
        return arr

    raise ValueError(
        "Could not align depth cube shape {} with DEM shape {}. "
        "Need to inspect MAT layout.".format(s, dem_shape)
    )


# ============================================================
# HELPERS - DEM / REPROJECTION
# ============================================================
def load_dem_single_band(path):
    with rasterio.open(str(path)) as src:
        dem = src.read(1).astype(np.float32)
        meta = {
            "crs": src.crs,
            "transform": src.transform,
            "width": src.width,
            "height": src.height,
            "bounds": src.bounds,
            "nodata": src.nodata,
        }
    return dem, meta


def extent_from_transform(transform, width, height):
    west, south, east, north = array_bounds(height, width, transform)
    return [west, east, south, north]


def build_target_grid_from_dem(dem_meta, target_crs):
    transform, width, height = calculate_default_transform(
        dem_meta["crs"],
        target_crs,
        dem_meta["width"],
        dem_meta["height"],
        *dem_meta["bounds"]
    )

    return {
        "crs": target_crs,
        "transform": transform,
        "width": width,
        "height": height,
    }


def reproject_array_to_target(src_array, src_crs, src_transform, target_grid,
                              src_nodata=np.nan,
                              resampling=Resampling.bilinear,
                              dst_dtype=np.float32):
    dst = np.full(
        (target_grid["height"], target_grid["width"]),
        np.nan,
        dtype=dst_dtype
    )

    reproject(
        source=src_array,
        destination=dst,
        src_transform=src_transform,
        src_crs=src_crs,
        src_nodata=src_nodata,
        dst_transform=target_grid["transform"],
        dst_crs=target_grid["crs"],
        dst_nodata=np.nan,
        resampling=resampling
    )

    return dst


def compute_plot_window_from_extent(extent, pad_frac_x, pad_frac_y):
    xmin, xmax, ymin, ymax = extent

    if xmax < xmin:
        xmin, xmax = xmax, xmin
    if ymax < ymin:
        ymin, ymax = ymax, ymin

    dx = xmax - xmin
    dy = ymax - ymin

    padx = pad_frac_x * dx
    pady = pad_frac_y * dy

    return [xmin - padx, xmax + padx, ymin - pady, ymax + pady]


def build_dem_valid_mask(dem, dem_meta):
    if dem_meta["nodata"] is not None:
        valid_mask = np.isfinite(dem) & (dem != dem_meta["nodata"])
    else:
        valid_mask = np.isfinite(dem)

    return valid_mask.astype(np.uint8)


def extract_outline_geometries_from_mask(mask_uint8, transform):
    """
    Convert a raster valid-data mask into polygon geometries.
    Returns a list of GeoJSON-like geometries in the DEM CRS.
    """
    geoms = []

    for geom, value in features.shapes(mask_uint8, mask=mask_uint8.astype(bool), transform=transform):
        if int(value) == 1:
            geoms.append(geom)

    if not geoms:
        raise RuntimeError("No valid polygons could be extracted from the DEM mask.")

    return geoms


def reproject_geometries(geoms, src_crs, dst_crs):
    out = []
    for geom in geoms:
        out.append(transform_geom(src_crs, dst_crs, geom))
    return out


def _append_polygon_rings_to_lines(geom, lines_xy):
    gtype = geom.get("type", None)
    coords = geom.get("coordinates", None)

    if gtype == "Polygon":
        # exterior ring only
        exterior = coords[0]
        xs = [pt[0] for pt in exterior]
        ys = [pt[1] for pt in exterior]
        lines_xy.append((xs, ys))
        return

    if gtype == "MultiPolygon":
        for poly in coords:
            exterior = poly[0]
            xs = [pt[0] for pt in exterior]
            ys = [pt[1] for pt in exterior]
            lines_xy.append((xs, ys))
        return

    raise ValueError("Unsupported geometry type for plotting: {}".format(gtype))


def geoms_to_plot_lines(geoms):
    lines_xy = []
    for geom in geoms:
        _append_polygon_rings_to_lines(geom, lines_xy)
    return lines_xy


# ============================================================
# HELPERS - CHUNK / FRAME PROCESSING
# ============================================================
def scan_chunks(temp_dir, dem_shape):
    chunk_files = list_chunk_files(temp_dir)
    valid_mask = None
    total_frames = 0
    max_depth = None
    min_positive = np.inf

    print("Found {} chunk files".format(len(chunk_files)))

    for i, fp in enumerate(chunk_files, start=1):
        loaded = load_mat_file(fp)
        maps = get_maps_struct(loaded)
        if maps is None:
            raise RuntimeError("Could not find Maps structure in {}".format(fp))

        depth_m_raw = extract_depth_3d_from_maps(maps)
        print("Raw depth cube shape from MAT during scan: {}".format(depth_m_raw.shape))

        depth_m = standardize_depth_cube_to_dem(depth_m_raw, dem_shape)
        print("Standardized depth cube shape during scan: {}".format(depth_m.shape))

        depth_m[~np.isfinite(depth_m)] = np.nan

        local_valid = np.any(np.isfinite(depth_m) & (depth_m > DEPTH_THRESHOLD_M), axis=2)
        if valid_mask is None:
            valid_mask = local_valid.copy()
        else:
            valid_mask |= local_valid

        masked = depth_m.copy()
        masked[masked <= DEPTH_THRESHOLD_M] = np.nan

        if np.any(np.isfinite(masked)):
            local_max = np.nanmax(masked)
            if np.isfinite(local_max):
                if max_depth is None:
                    max_depth = local_max
                else:
                    max_depth = max(max_depth, local_max)

            finite_pos = masked[np.isfinite(masked)]
            if finite_pos.size:
                min_positive = min(min_positive, float(np.nanmin(finite_pos)))

        total_frames += depth_m.shape[2]
        print("Scanned chunk {}/{}: {} -> {} frames".format(
            i, len(chunk_files), fp.name, depth_m.shape[2]
        ))

    if valid_mask is None:
        raise RuntimeError("No valid depth data found in temp files")

    if max_depth is None or not np.isfinite(max_depth):
        max_depth = 1.0

    if not np.isfinite(min_positive):
        min_positive = DEPTH_THRESHOLD_M

    return {
        "chunk_files": chunk_files,
        "valid_mask": valid_mask,
        "total_frames": total_frames,
        "max_depth": float(max_depth),
        "min_positive": float(min_positive),
    }


def sample_chunk_files(chunk_files, sample_count):
    if sample_count is None or sample_count >= len(chunk_files):
        return chunk_files
    if sample_count <= 1:
        return [chunk_files[0]]
    idx = np.linspace(0, len(chunk_files) - 1, sample_count).round().astype(int)
    idx = np.unique(idx)
    return [chunk_files[i] for i in idx]


def compute_percentile_vmax_from_temp(chunk_files, percentile, sample_count, dem_shape):
    sampled = sample_chunk_files(chunk_files, sample_count)
    values = []

    for fp in sampled:
        loaded = load_mat_file(fp)
        maps = get_maps_struct(loaded)

        depth_m_raw = extract_depth_3d_from_maps(maps)
        depth_m = standardize_depth_cube_to_dem(depth_m_raw, dem_shape)

        depth_m[~np.isfinite(depth_m)] = np.nan
        depth_m[depth_m <= DEPTH_THRESHOLD_M] = np.nan

        finite = depth_m[np.isfinite(depth_m)]
        if finite.size:
            values.append(float(np.nanpercentile(finite, percentile)))

    if not values:
        return 1.0

    return max(values)


def get_depth_scale(scan_info, dem_shape):
    if SCALE_MODE == "manual":
        depth_min = float(MANUAL_DEPTH_MIN)
        depth_max = float(MANUAL_DEPTH_MAX)

        if depth_max <= depth_min:
            raise ValueError("MANUAL_DEPTH_MAX must be greater than MANUAL_DEPTH_MIN")

        print("Using manual depth scale from {:.2f} to {:.3f} m".format(depth_min, depth_max))
        return depth_min, depth_max

    if SCALE_MODE == "percentile":
        vmax = compute_percentile_vmax_from_temp(
            scan_info["chunk_files"],
            percentile=PERCENTILE_FOR_MAX,
            sample_count=PERCENTILE_SAMPLE_COUNT,
            dem_shape=dem_shape
        )

        if vmax <= 0:
            vmax = 1.0

        print("Using percentile depth scale from 0.00 to {:.3f} m".format(vmax))
        return 0.0, float(vmax)

    raise ValueError("SCALE_MODE must be 'manual' or 'percentile'")


# ============================================================
# PLOTTING
# ============================================================
def add_plot_layers(ax, flood_arr, flood_extent, norm, cmap, colorbar_ticks,
                    plot_window, catchment_outline_lines):
    flood_im = ax.imshow(
        flood_arr,
        extent=flood_extent,
        origin="upper",
        cmap=cmap,
        norm=norm,
        alpha=ALPHA,
        zorder=3
    )

    ax.set_xlim(plot_window[0], plot_window[1])
    ax.set_ylim(plot_window[2], plot_window[3])

    cx.add_basemap(
        ax,
        crs="EPSG:3857",
        source=BASEMAP_SOURCE,
        zoom=BASEMAP_ZOOM,
        attribution=BASEMAP_ATTRIBUTION,
        reset_extent=False
    )

    for xs, ys in catchment_outline_lines:
        ax.plot(
            xs,
            ys,
            color=OUTLINE_COLOR,
            linewidth=OUTLINE_LINEWIDTH,
            zorder=OUTLINE_ZORDER
        )

    cbar = plt.colorbar(
        flood_im,
        ax=ax,
        fraction=COLORBAR_FRACTION,
        pad=COLORBAR_PAD,
        shrink=COLORBAR_SHRINK,
        aspect=COLORBAR_ASPECT,
        extend=COLORBAR_EXTEND,
        extendfrac=COLORBAR_EXTEND_FRAC
    )
    cbar.set_ticks(colorbar_ticks)
    cbar.set_label(CBAR_LABEL, fontsize=CBAR_LABEL_FONTSIZE)
    cbar.ax.tick_params(labelsize=CBAR_TICK_FONTSIZE)

    ax.set_xlabel(X_LABEL, fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(Y_LABEL, fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(axis="both", labelsize=AXIS_TICK_FONTSIZE)


def save_frame_png(output_png, flood_arr, timestamp, frame_index, total_frames,
                   flood_extent, norm, cmap, colorbar_ticks, plot_window,
                   catchment_outline_lines):
    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)

    add_plot_layers(
        ax=ax,
        flood_arr=flood_arr,
        flood_extent=flood_extent,
        norm=norm,
        cmap=cmap,
        colorbar_ticks=colorbar_ticks,
        plot_window=plot_window,
        catchment_outline_lines=catchment_outline_lines
    )

    if FLAG_ELAPSED_TIME:
        time_str = "t = {} h".format(frame_index * MAP_SAVE_TIMESTEP_HOURS)
    else:
        time_str = timestamp.strftime("%Y-%m-%d %H:%M:%S")

    ax.set_title(TITLE_TEMPLATE.format(timestamp=time_str), fontsize=TITLE_FONTSIZE)

    ax.text(
        0.02,
        0.02,
        FRAME_LABEL_TEMPLATE.format(current=frame_index + 1, total=total_frames),
        transform=ax.transAxes,
        fontsize=FRAME_TEXT_FONTSIZE,
        bbox=dict(facecolor="white", alpha=0.8, edgecolor="none")
    )

    output_png = Path(output_png)
    output_png.parent.mkdir(parents=True, exist_ok=True)

    fig.savefig(
        str(output_png),
        dpi=DPI,
        bbox_inches=SAVE_BBOX_INCHES,
        facecolor=PNG_FACE_COLOR
    )
    plt.close(fig)

    print("Saved frame: {}".format(output_png))


# ============================================================
# MAIN
# ============================================================
def main():
    print("Checking inputs...")
    print("TEMP_DIR: {}".format(TEMP_DIR))
    print("DEM_TIF: {}".format(DEM_TIF))
    print("OUTPUT_FRAMES_DIR: {}".format(OUTPUT_FRAMES_DIR))
    print("TEST_OUTPUT_PNG: {}".format(TEST_OUTPUT_PNG))
    print("SCALE_MODE: {}".format(SCALE_MODE))

    if not Path(TEMP_DIR).exists():
        raise FileNotFoundError("TEMP_DIR does not exist: {}".format(TEMP_DIR))
    if not Path(DEM_TIF).exists():
        raise FileNotFoundError("DEM_TIF does not exist: {}".format(DEM_TIF))

    setup_font(FONT_PATH, FALLBACK_FONT_FAMILY)

    dem, dem_meta = load_dem_single_band(DEM_TIF)

    if dem_meta["crs"] is None:
        raise ValueError("DEM has no CRS defined")

    dem_shape = dem.shape
    print("DEM shape: {}".format(dem_shape))
    print("DEM native bounds: {}".format(dem_meta["bounds"]))
    print("DEM CRS: {}".format(dem_meta["crs"]))

    # Build catchment-like outline from valid DEM mask
    dem_valid_mask = build_dem_valid_mask(dem, dem_meta)
    dem_mask_geoms = extract_outline_geometries_from_mask(
        dem_valid_mask,
        dem_meta["transform"]
    )
    dem_mask_geoms_3857 = reproject_geometries(
        dem_mask_geoms,
        dem_meta["crs"],
        "EPSG:3857"
    )
    catchment_outline_lines = geoms_to_plot_lines(dem_mask_geoms_3857)
    print("Extracted {} outline polygon(s) from DEM valid-data mask.".format(
        len(catchment_outline_lines)
    ))

    scan_info = scan_chunks(TEMP_DIR, dem_shape)
    print("Total recovered frames: {}".format(scan_info["total_frames"]))

    target_grid = build_target_grid_from_dem(dem_meta, "EPSG:3857")

    flood_extent = extent_from_transform(
        target_grid["transform"],
        target_grid["width"],
        target_grid["height"]
    )

    plot_window = compute_plot_window_from_extent(
        flood_extent,
        MAP_PAD_FRACTION_X,
        MAP_PAD_FRACTION_Y
    )

    print("Projected plotting extent (EPSG:3857): {}".format(flood_extent))
    print("Plot window: {}".format(plot_window))

    depth_min, depth_max = get_depth_scale(scan_info, dem_shape)

    levels = np.linspace(depth_min, depth_max, N_COLOR_LEVELS + 1)
    cmap = plt.get_cmap(CMAP, N_COLOR_LEVELS)
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
    colorbar_ticks = levels

    output_dir = Path(OUTPUT_FRAMES_DIR)
    output_dir.mkdir(parents=True, exist_ok=True)

    valid_mask = scan_info["valid_mask"]
    total_frames = scan_info["total_frames"]

    global_frame_idx = 0

    for chunk_idx, fp in enumerate(scan_info["chunk_files"], start=1):
        loaded = load_mat_file(fp)
        maps = get_maps_struct(loaded)

        depth_3d_raw = extract_depth_3d_from_maps(maps)
        print("Raw depth cube shape from MAT: {}".format(depth_3d_raw.shape))

        depth_3d = standardize_depth_cube_to_dem(depth_3d_raw, dem_shape)
        print("Standardized depth cube shape: {}".format(depth_3d.shape))

        n_local = depth_3d.shape[2]
        print("Rendering chunk {}/{}: {} ({} frames)".format(
            chunk_idx, len(scan_info["chunk_files"]), fp.name, n_local
        ))

        for j in range(n_local):
            depth_m = depth_3d[:, :, j].astype(np.float32)
            depth_m[~np.isfinite(depth_m)] = np.nan
            depth_m[depth_m <= DEPTH_THRESHOLD_M] = np.nan
            depth_m[~valid_mask] = np.nan

            if depth_m.shape != dem_shape:
                raise ValueError(
                    "Depth slice shape {} does not match DEM shape {} after standardization.".format(
                        depth_m.shape, dem_shape
                    )
                )

            flood_proj = reproject_array_to_target(
                src_array=depth_m,
                src_crs=dem_meta["crs"],
                src_transform=dem_meta["transform"],
                target_grid=target_grid,
                src_nodata=np.nan,
                resampling=Resampling.bilinear
            )

            flood_proj = np.where(
                np.isfinite(flood_proj) & (flood_proj >= DISPLAY_THRESHOLD_M),
                flood_proj,
                np.nan
            )

            timestamp = frame_time_from_index(global_frame_idx)

            if TEST_SINGLE_FRAME:
                if global_frame_idx != TEST_FRAME_INDEX:
                    global_frame_idx += 1
                    continue

                save_frame_png(
                    output_png=TEST_OUTPUT_PNG,
                    flood_arr=flood_proj,
                    timestamp=timestamp,
                    frame_index=global_frame_idx,
                    total_frames=total_frames,
                    flood_extent=flood_extent,
                    norm=norm,
                    cmap=cmap,
                    colorbar_ticks=colorbar_ticks,
                    plot_window=plot_window,
                    catchment_outline_lines=catchment_outline_lines
                )
                print("Saved single test frame to {}".format(TEST_OUTPUT_PNG))
                return

            output_png = output_dir / "frame_{:04d}.png".format(global_frame_idx)
            save_frame_png(
                output_png=output_png,
                flood_arr=flood_proj,
                timestamp=timestamp,
                frame_index=global_frame_idx,
                total_frames=total_frames,
                flood_extent=flood_extent,
                norm=norm,
                cmap=cmap,
                colorbar_ticks=colorbar_ticks,
                plot_window=plot_window,
                catchment_outline_lines=catchment_outline_lines
            )

            global_frame_idx += 1

            if (global_frame_idx % PROGRESS_EVERY_N_FRAMES == 0) or (global_frame_idx == total_frames):
                print("Saved {}/{} frames".format(global_frame_idx, total_frames))

    print("Done. Frames saved to {}".format(output_dir))


if __name__ == "__main__":
    main()