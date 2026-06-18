
'''
source /oak/stanford/groups/gorelick/Marcus/python_env/my_env/bin/activate
cd /oak/stanford/groups/gorelick/HydroPol2D/HydroPol2D_Functions
python /oak/stanford/groups/gorelick/HydroPol2D/HydroPol2D_Functions/generate_animations.py
'''

# ============================================================
# HYDROPOL2D MAP EXPORT + TEMP ANIMATION FRAME GENERATOR
# Flood depth PNGs + optional GeoTIFFs + maximum flood depth GeoTIFF
# I_t PNGs + optional GeoTIFFs from save_input_maps_*.mat
# Hillshade background + optional basemap + manual/auto scales
# Python 3.6 compatible
# ============================================================

from pathlib import Path
from datetime import datetime, timedelta
import argparse
import re
import math

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.animation as mpl_animation
from matplotlib.colors import BoundaryNorm, Normalize
from matplotlib import font_manager

import rasterio
from rasterio.warp import calculate_default_transform, reproject, transform_geom
from rasterio.enums import Resampling
from rasterio.transform import array_bounds, from_bounds
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
TEMP_DIR = r"/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/India/2000m/Outputs/Temporary_Files"
DEM_TIF = r"/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/India/2000m/Static/DEM_fabdem.tif"

# Existing PNG frame output folder for flood depth.
OUTPUT_FRAMES_DIR = r"/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/India/2000m/Outputs/Png_Animations"

# NEW: choose the folder where all GeoTIFFs and new map products are exported.
# Change only this path when you want outputs somewhere else.
EXPORT_DATA_DIR = r"/oak/stanford/groups/gorelick/HydroPol2D/Case_Studies/India/2000m/Outputs/Exported_GeoTIFFs"

# Folder containing save_input_maps_*.mat. Usually it is the same Temporary_Files folder.
INPUT_MAPS_DIR = TEMP_DIR

# PNG folder for I_t frames.
OUTPUT_IT_FRAMES_DIR = str(Path(EXPORT_DATA_DIR) / "I_t_Png_Animations")

TEST_OUTPUT_PNG = str(Path(OUTPUT_FRAMES_DIR) / "test_image.png")
TEST_OUTPUT_IT_PNG = str(Path(OUTPUT_IT_FRAMES_DIR) / "test_I_t_image.png")


# ============================================================
# USER INPUTS - MAIN PROCESS SWITCHES
# ============================================================
# Keep True if you still want to generate flood-depth PNG frames.
PROCESS_FLOOD_DEPTH_PNGS = True

# Save maximum flood depth over all saved depth maps as one GeoTIFF.
EXPORT_MAX_FLOOD_DEPTH_TIF = True

# Save maximum flood depth as one PNG using the same visual style.
EXPORT_MAX_FLOOD_DEPTH_PNG = True

# Save every saved flood depth map as an individual GeoTIFF.
# WARNING: this can create many files.
EXPORT_EACH_FLOOD_DEPTH_TIF = False

# Process I_t from save_input_maps_*.mat.
PROCESS_IT_MAPS = True

# Plot every I_t map as PNG using the same style.
PLOT_EACH_IT_PNG = True

# Save every I_t map as individual GeoTIFF.
# WARNING: this can create many files.
EXPORT_EACH_IT_TIF = False

# Save maximum I_t over all saved input maps as GeoTIFF and PNG.
EXPORT_MAX_IT_TIF = True
EXPORT_MAX_IT_PNG = True

# NEW PERFORMANCE SWITCHES
# False = do not read all .mat files once just to count frames for progress reporting.
# This avoids reading save_map_hydro_*.mat twice.
READ_MAT_FILES_FOR_FRAME_COUNTS = False

# False = do not scan save_input_maps_*.mat before processing just to auto-detect
# the I_t colorbar maximum. This avoids reading save_input_maps_*.mat twice.
# If this is False and IT_MAX_MM is None, per-frame I_t PNGs use IT_DEFAULT_MAX_MM.
SCAN_IT_FOR_AUTO_COLORBAR = False

# Used only when IT_MAX_MM is None and SCAN_IT_FOR_AUTO_COLORBAR = False.
# Change this based on the expected I_t range for your event.
IT_DEFAULT_MAX_MM = 100.0


# ============================================================
# USER INPUTS - RUN MODE
# ============================================================
# TEST_SINGLE_FRAME applies only to flood-depth PNG generation.
# GeoTIFF maximum export needs the full loop, so it is skipped in test mode.
TEST_SINGLE_FRAME = False
TEST_FRAME_INDEX = 1500

# I_t test mode. If True, only one I_t PNG is generated.
TEST_SINGLE_IT_FRAME = False
TEST_IT_FRAME_INDEX = 0


# ============================================================
# USER INPUTS - TIME SETTINGS
# ============================================================
DATE_BEGIN = datetime(2024, 7, 17, 0, 0, 0)
MAP_SAVE_TIMESTEP_HOURS = 1
FLAG_ELAPSED_TIME = False


# ============================================================
# USER INPUTS - DEPTH / MASK SETTINGS
# ============================================================
# The model depth d is commonly saved in mm. The script converts d/1000 to meters.
DEPTH_SAVED_IN_MM = True

DEPTH_THRESHOLD_M = 0.00
DISPLAY_THRESHOLD_M = 0.00

MIN_FLOOD_DEPTH_M = 0.0
MAX_FLOOD_DEPTH_M = 2.0

FLOOD_CMAP = "cool"
N_COLOR_INTERVALS = 10
FLOOD_ALPHA = 0.72

COLOR_MODE = "discrete"  # "continuous" OR "discrete"


# ============================================================
# USER INPUTS - I_t SETTINGS
# ============================================================
# I_t is saved in mm according to your note. No conversion is applied.
IT_VARIABLE_CANDIDATES = [
    "I_t", "It", "I", "I_mm", "Infiltration", "infiltration", "infiltration_mm"
]

IT_DISPLAY_THRESHOLD_MM = 0.00
IT_MIN_MM = 0.0

# Use None for automatic maximum from all save_input_maps files.
# Or set a number, e.g., 100.0, for a fixed colorbar.
IT_MAX_MM = None

IT_CMAP = "turbo"
IT_ALPHA = 0.72
IT_TITLE_TEMPLATE = "Infiltration $I_t$\n{timestamp}"
IT_MAX_TITLE_TEMPLATE = "Maximum Infiltration $I_t$"
IT_CBAR_LABEL = "$I_t$ (mm)"


# ============================================================
# USER INPUTS - GEOTIFF EXPORT SETTINGS
# ============================================================
# "dem" = write GeoTIFFs in the DEM native CRS/resolution.
# "web_mercator" = write GeoTIFFs in EPSG:3857, matching the basemap plot grid.
# Recommendation for scientific analysis: keep "dem".
GEOTIFF_EXPORT_CRS_MODE = "dem"  # "dem" OR "web_mercator"
GEOTIFF_NODATA_VALUE = -9999.0
GEOTIFF_COMPRESS = "deflate"


# ============================================================
# USER INPUTS - HILLSHADE CUSTOMIZATION
# ============================================================
USE_HILLSHADE = True
HILLSHADE_ALPHA = 0.25
HILLSHADE_AZIMUTH = 315.0
HILLSHADE_ALTITUDE = 45.0
HILLSHADE_ZORDER = 1


# ============================================================
# USER INPUTS - BASEMAP CUSTOMIZATION
# ============================================================
USE_CONTEXTILY_BASEMAP = True

BASEMAP_SOURCE = cx.providers.CartoDB.DarkMatter
BASEMAP_ZOOM = "auto"
BASEMAP_ATTRIBUTION = False
BASEMAP_ZORDER = 0


# ============================================================
# USER INPUTS - FONT CUSTOMIZATION
# ============================================================
FONT_PATH = "/oak/stanford/groups/gorelick/Marcus/misc/Avenir/AvenirNextCyr-Regular.ttf"
FALLBACK_FONT_FAMILY = "DejaVu Sans"


# ============================================================
# USER INPUTS - TEXT CUSTOMIZATION
# ============================================================
TITLE_TEMPLATE = "Flood Depth\n{timestamp}"
MAX_DEPTH_TITLE_TEMPLATE = "Maximum Flood Depth"
CBAR_LABEL = "Flood depth (m)"
X_LABEL = "X [m]"
Y_LABEL = "Y [m]"
FRAME_LABEL_TEMPLATE = "Frame {current}"

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
# USER INPUTS - COLORBAR CUSTOMIZATION
# ============================================================
COLORBAR_EXTEND = "max"
COLORBAR_FRACTION = 0.08
COLORBAR_SHRINK = 0.50
COLORBAR_PAD = 0.02
COLORBAR_EXTEND_FRAC = 0.08
COLORBAR_ASPECT = 20


# ============================================================
# USER INPUTS - CATCHMENT / OUTLINE CUSTOMIZATION
# ============================================================
OUTLINE_COLOR = "white"
OUTLINE_LINEWIDTH = 2
OUTLINE_ZORDER = 4


# ============================================================
# USER INPUTS - EXPORT CUSTOMIZATION
# ============================================================
SAVE_BBOX_INCHES = "tight"
PNG_FACE_COLOR = "white"
PROGRESS_EVERY_N_FRAMES = 10


# ============================================================
# GENERAL UTILITIES
# ============================================================
def setup_font(font_path, fallback_family):
    if font_path is not None:
        font_path = Path(font_path)
        if font_path.exists():
            font_manager.fontManager.addfont(str(font_path))
            font_name = font_manager.FontProperties(fname=str(font_path)).get_name()
            mpl.rcParams["font.family"] = font_name
            print("Using custom font: {}".format(font_name))
        else:
            mpl.rcParams["font.family"] = fallback_family
            print("FONT_PATH not found. Using fallback font: {}".format(fallback_family))
    else:
        mpl.rcParams["font.family"] = fallback_family
        print("Using fallback font: {}".format(fallback_family))


def _extract_number_from_filename(path, pattern):
    m = re.search(pattern, path.name)
    if m is None:
        return math.inf
    return int(m.group(1))


def list_hydro_chunk_files(temp_dir):
    temp_dir = Path(temp_dir)
    if not temp_dir.exists():
        raise FileNotFoundError("TEMP_DIR does not exist: {}".format(temp_dir))

    files = list(temp_dir.glob("save_map_hydro_*.mat"))
    if not files:
        raise FileNotFoundError("No save_map_hydro_*.mat files found in {}".format(temp_dir))

    return sorted(files, key=lambda p: _extract_number_from_filename(p, r"save_map_hydro_(\d+)\.mat$"))


def list_input_map_files(input_dir):
    input_dir = Path(input_dir)
    if not input_dir.exists():
        raise FileNotFoundError("INPUT_MAPS_DIR does not exist: {}".format(input_dir))

    files = list(input_dir.glob("save_input_maps_*.mat"))
    if not files:
        raise FileNotFoundError("No save_input_maps_*.mat files found in {}".format(input_dir))

    return sorted(files, key=lambda p: _extract_number_from_filename(p, r"save_input_maps_(\d+)\.mat$"))


def frame_time_from_index(frame_index):
    return DATE_BEGIN + timedelta(hours=frame_index * MAP_SAVE_TIMESTEP_HOURS)


def format_time_label(frame_index):
    if FLAG_ELAPSED_TIME:
        return "t = {} h".format(frame_index * MAP_SAVE_TIMESTEP_HOURS)
    return frame_time_from_index(frame_index).strftime("%Y-%m-%d %H:%M:%S")


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
    if isinstance(loaded, dict) and "Maps" in loaded:
        return loaded["Maps"]

    if isinstance(loaded, dict):
        for v in loaded.values():
            hydro = maybe_getattr(v, "Hydro", None)
            if hydro is not None:
                return v

    return None


# ============================================================
# MATLAB STRUCT / VARIABLE SEARCH UTILITIES
# ============================================================
def _normalize_key(name):
    return str(name).lower().replace("_", "")


def _is_numeric_array(obj):
    try:
        arr = np.asarray(obj)
    except Exception:
        return False

    if arr.dtype.kind in ("f", "i", "u", "b"):
        return True
    return False


def _iter_named_children(obj):
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield str(k), v
        return

    # scipy.io.loadmat MATLAB struct objects often have _fieldnames.
    fieldnames = getattr(obj, "_fieldnames", None)
    if fieldnames is not None:
        for k in fieldnames:
            yield str(k), getattr(obj, k)
        return

    # Some loaded structs can be object arrays.
    if isinstance(obj, np.ndarray) and obj.dtype.kind == "O" and obj.size <= 100:
        flat = obj.flat
        for idx, v in enumerate(flat):
            yield "item{}".format(idx), v
        return


def find_first_named_variable(obj, candidate_names, root_name="root", max_depth=12):
    """
    Recursively search a loaded MATLAB structure for the first numeric array whose
    field/key name matches one of candidate_names.
    """
    exact = set([str(x).lower() for x in candidate_names])
    normalized = set([_normalize_key(x) for x in candidate_names])

    def _search(current, path, depth):
        if depth > max_depth:
            return None, None

        # First pass: prefer exact field/key name matches at this level.
        children = list(_iter_named_children(current) or [])
        for child_name, child_value in children:
            child_name_lower = child_name.lower()
            child_name_norm = _normalize_key(child_name)
            if (child_name_lower in exact) or (child_name_norm in normalized):
                if _is_numeric_array(child_value):
                    return child_value, path + "." + child_name

        # Second pass: recurse.
        for child_name, child_value in children:
            if _is_numeric_array(child_value):
                continue
            found_value, found_path = _search(child_value, path + "." + child_name, depth + 1)
            if found_value is not None:
                return found_value, found_path

        return None, None

    return _search(obj, root_name, 0)


# ============================================================
# ARRAY STANDARDIZATION
# ============================================================
def _numeric_to_3d_cube(arr, variable_name):
    arr = _safe_squeeze(arr)
    arr = np.asarray(arr)

    if arr.dtype.kind not in ("f", "i", "u", "b"):
        raise ValueError("{} is not numeric. dtype={}".format(variable_name, arr.dtype))

    if arr.ndim == 2:
        arr = arr[:, :, np.newaxis]
    elif arr.ndim != 3:
        raise ValueError("{} must be 2D or 3D. Got shape {}".format(variable_name, arr.shape))

    return arr.astype(np.float32)


def _cube_orientation_options(cube_3d):
    """
    Return possible interpretations of a MATLAB map cube as (rows, cols, time).

    HydroPol2D/MATLAB outputs are sometimes saved as:
      - (rows, cols, time)
      - (cols, rows, time)
      - (time, rows, cols)
      - (time, cols, rows)

    Some model arrays may also omit the outer ghost/boundary cells, so padded
    variants are included in standardize_cube_to_grid().
    """
    s = cube_3d.shape
    opts = []

    # Native MATLAB-like order already interpreted as rows, cols, time.
    opts.append(("native rows,cols,time", (s[0], s[1]), lambda a: a))

    # Common row/column transpose.
    opts.append(("transposed cols,rows,time", (s[1], s[0]), lambda a: np.transpose(a, (1, 0, 2))))

    # Time-first cases.
    opts.append(("time,rows,cols", (s[1], s[2]), lambda a: np.transpose(a, (1, 2, 0))))
    opts.append(("time,cols,rows", (s[2], s[1]), lambda a: np.transpose(a, (2, 1, 0))))

    return opts


def infer_model_grid_shape_from_cube(cube_3d, reference_shape, variable_name="map"):
    """
    Infer the model map grid shape from a raw output cube, using the DEM only as
    a reference for orientation. This is intentionally tolerant of small DEM vs.
    output-map dimension differences.
    """
    if cube_3d.ndim != 3:
        raise ValueError("{} cube must be 3D. Got shape {}".format(variable_name, cube_3d.shape))

    ref_rows, ref_cols = reference_shape
    candidates = []

    for desc, shape, _ in _cube_orientation_options(cube_3d):
        rows, cols = shape
        if rows <= 0 or cols <= 0:
            continue

        # Penalize implausible maps where one spatial dimension is tiny compared
        # with the DEM. This helps avoid interpreting the time dimension as space.
        size_score = abs(rows - ref_rows) + abs(cols - ref_cols)
        ratio_score = abs(math.log(max(rows, 1) / float(max(ref_rows, 1)))) + abs(math.log(max(cols, 1) / float(max(ref_cols, 1))))
        score = size_score + 1000.0 * ratio_score
        candidates.append((score, desc, (int(rows), int(cols))))

    if not candidates:
        raise ValueError("Could not infer model grid shape from {} cube shape {}".format(
            variable_name, cube_3d.shape
        ))

    candidates.sort(key=lambda x: x[0])
    best_score, best_desc, best_shape = candidates[0]

    print(
        "Inferred {} model grid shape {} from raw cube shape {} using orientation '{}' "
        "and reference DEM shape {}.".format(
            variable_name, best_shape, cube_3d.shape, best_desc, reference_shape
        )
    )

    return best_shape


def standardize_cube_to_grid(cube_3d, target_shape, variable_name="map"):
    """
    Standardize any supported HydroPol2D/MATLAB map cube to (rows, cols, time)
    using the target model grid shape, not necessarily the original DEM shape.
    """
    target_rows, target_cols = target_shape

    if cube_3d.ndim != 3:
        raise ValueError("{} cube must be 3D. Got shape {}".format(variable_name, cube_3d.shape))

    # Exact-orientation matches.
    for desc, shape, op in _cube_orientation_options(cube_3d):
        if shape == (target_rows, target_cols):
            if desc != "native rows,cols,time":
                print("{} cube appears to be {}. Fixing orientation.".format(variable_name, desc))
            return op(cube_3d).astype(np.float32)

    # Interior maps that need a 1-cell edge pad.
    for desc, shape, op in _cube_orientation_options(cube_3d):
        if shape == (target_rows - 2, target_cols - 2):
            print("{} cube appears to be {} without outer 1-cell edge.".format(variable_name, desc))
            print("Padding 1 cell on all sides to match target grid {}.".format(target_shape))
            arr = op(cube_3d)
            return np.pad(arr, ((1, 1), (1, 1), (0, 0)), mode="constant", constant_values=np.nan).astype(np.float32)

    # One-column/one-row missing cases sometimes occur after clipping.
    for desc, shape, op in _cube_orientation_options(cube_3d):
        rows, cols = shape
        if rows == target_rows and cols == target_cols - 1:
            print("{} cube appears to be {} with one missing column. Padding right edge.".format(variable_name, desc))
            arr = op(cube_3d)
            return np.pad(arr, ((0, 0), (0, 1), (0, 0)), mode="constant", constant_values=np.nan).astype(np.float32)
        if rows == target_rows - 1 and cols == target_cols:
            print("{} cube appears to be {} with one missing row. Padding bottom edge.".format(variable_name, desc))
            arr = op(cube_3d)
            return np.pad(arr, ((0, 1), (0, 0), (0, 0)), mode="constant", constant_values=np.nan).astype(np.float32)

    raise ValueError(
        "Could not align {} cube shape {} with target model grid shape {}. "
        "Need to inspect MAT layout.".format(variable_name, cube_3d.shape, target_shape)
    )


# Backward-compatible name used in older parts of this script.
def standardize_cube_to_dem(cube_3d, dem_shape, variable_name="map"):
    return standardize_cube_to_grid(cube_3d, dem_shape, variable_name=variable_name)

def extract_depth_3d_from_maps(maps):
    hydro = maybe_getattr(maps, "Hydro", None)
    if hydro is None:
        raise KeyError("Could not find Maps.Hydro in MAT file")

    d = maybe_getattr(hydro, "d", None)
    if d is None:
        raise KeyError("Could not find Maps.Hydro.d in MAT file")

    d = _numeric_to_3d_cube(d, "Maps.Hydro.d")

    if DEPTH_SAVED_IN_MM:
        depth_m = d / 1000.0
    else:
        depth_m = d

    return depth_m.astype(np.float32)


def extract_it_3d_from_loaded_input(loaded, target_shape):
    it_raw, found_path = find_first_named_variable(
        loaded,
        IT_VARIABLE_CANDIDATES,
        root_name="input_mat",
        max_depth=12
    )

    if it_raw is None:
        raise KeyError(
            "Could not find I_t in save_input_maps file. Tried candidate names: {}".format(
                ", ".join(IT_VARIABLE_CANDIDATES)
            )
        )

    it_cube_raw = _numeric_to_3d_cube(it_raw, found_path)
    it_cube = standardize_cube_to_grid(it_cube_raw, target_shape, variable_name="I_t")
    print("Detected I_t variable at: {} | standardized shape: {}".format(found_path, it_cube.shape))

    return it_cube.astype(np.float32), found_path


# ============================================================
# RASTER / PROJECTION UTILITIES
# ============================================================
def load_dem_single_band(path):
    with rasterio.open(str(path)) as src:
        dem = src.read(1).astype(np.float32)
        nodata = src.nodata

        if nodata is not None:
            dem = np.where(dem == nodata, np.nan, dem)

        meta = {
            "crs": src.crs,
            "transform": src.transform,
            "width": src.width,
            "height": src.height,
            "bounds": src.bounds,
            "nodata": nodata,
        }

    return dem, meta


def infer_model_grid_shape(reference_dem_shape):
    """
    Read only the first available output-map MAT file to determine the model map
    grid. This grid becomes the working grid for the DEM, hillshade, masks, PNGs,
    and native GeoTIFF exports.
    """
    # Prefer flood-depth outputs when they are part of the requested processing,
    # because save_output_maps/save_map_hydro maps define the hydraulic output grid.
    use_depth = (
        PROCESS_FLOOD_DEPTH_PNGS
        or EXPORT_MAX_FLOOD_DEPTH_TIF
        or EXPORT_EACH_FLOOD_DEPTH_TIF
        or EXPORT_MAX_FLOOD_DEPTH_PNG
    )

    if use_depth:
        chunk_files = list_hydro_chunk_files(TEMP_DIR)
        fp = chunk_files[0]
        loaded = load_mat_file(fp)
        maps = get_maps_struct(loaded)
        if maps is None:
            raise RuntimeError("Could not find Maps structure in {}".format(fp))
        depth_raw = extract_depth_3d_from_maps(maps)
        return infer_model_grid_shape_from_cube(
            depth_raw,
            reference_shape=reference_dem_shape,
            variable_name="depth"
        )

    if PROCESS_IT_MAPS:
        input_files = list_input_map_files(INPUT_MAPS_DIR)
        fp = input_files[0]
        loaded = load_mat_file(fp)
        it_raw, found_path = find_first_named_variable(
            loaded,
            IT_VARIABLE_CANDIDATES,
            root_name="input_mat",
            max_depth=12
        )
        if it_raw is None:
            raise KeyError(
                "Could not find I_t in first save_input_maps file. Tried candidate names: {}".format(
                    ", ".join(IT_VARIABLE_CANDIDATES)
                )
            )
        it_cube_raw = _numeric_to_3d_cube(it_raw, found_path)
        return infer_model_grid_shape_from_cube(
            it_cube_raw,
            reference_shape=reference_dem_shape,
            variable_name="I_t"
        )

    return reference_dem_shape


def resample_dem_to_model_grid_if_needed(dem, dem_meta, model_shape):
    """
    If the DEM dimensions do not match the saved model maps, resample the DEM to
    the model-map matrix size while preserving the DEM CRS and spatial bounds.

    IMPORTANT:
      - This does NOT reproject the DEM to a new CRS.
      - It only changes the DEM raster dimensions/resolution so its matrix size
        matches the HydroPol2D saved output maps.
      - The resulting working grid remains in the DEM CRS and uses the same
        spatial bounds as the original DEM.
    """
    model_rows, model_cols = model_shape

    if dem.shape == (model_rows, model_cols):
        print("DEM shape already matches model output grid: {}".format(dem.shape))
        return dem, dem_meta

    print(
        "DEM shape {} does not match model output grid {}. Resampling DEM to model grid.".format(
            dem.shape,
            (model_rows, model_cols)
        )
    )

    left, bottom, right, top = dem_meta["bounds"]
    dst_transform = from_bounds(left, bottom, right, top, model_cols, model_rows)

    # Same CRS, same bounds, new row/column count. This aligns DEM-derived
    # hillshade/mask/outline with the HydroPol2D map matrices.
    target_grid_native = {
        "crs": dem_meta["crs"],
        "transform": dst_transform,
        "width": int(model_cols),
        "height": int(model_rows),
    }

    dst = reproject_array_to_target(
        src_array=dem.astype(np.float32),
        src_crs=dem_meta["crs"],
        src_transform=dem_meta["transform"],
        target_grid=target_grid_native,
        src_nodata=np.nan,
        resampling=Resampling.bilinear,
        dst_dtype=np.float32
    )

    new_meta = dict(dem_meta)
    new_meta["transform"] = dst_transform
    new_meta["width"] = int(model_cols)
    new_meta["height"] = int(model_rows)
    new_meta["bounds"] = rasterio.coords.BoundingBox(left=left, bottom=bottom, right=right, top=top)

    print("Resampled DEM shape: {}".format(dst.shape))
    print("Working DEM/model CRS: {}".format(new_meta["crs"]))
    print("Updated DEM/model transform: {}".format(dst_transform))

    return dst.astype(np.float32), new_meta

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
    """
    Reproject/resample a 2D array to target_grid.

    This helper is used for two different purposes:
      1) same-CRS resampling of the DEM to the model output matrix size; and
      2) map plotting reprojection to EPSG:3857 when using a contextily basemap.

    It handles NaN nodata robustly by temporarily replacing NaNs with a numeric
    sentinel before calling rasterio.warp.reproject. This avoids cases where
    hillshade or data layers become fully NaN after reprojection.
    """
    src = np.asarray(src_array, dtype=np.float32)

    internal_nodata = -9999999.0  # safe numeric float32 nodata sentinel for rasterio/GDAL

    if src_nodata is None:
        src_for_warp = src
        actual_src_nodata = None
        actual_dst_nodata = internal_nodata
    else:
        try:
            src_nodata_is_nan = bool(np.isnan(src_nodata))
        except Exception:
            src_nodata_is_nan = False

        if src_nodata_is_nan:
            src_for_warp = np.where(np.isfinite(src), src, internal_nodata).astype(np.float32)
            actual_src_nodata = internal_nodata
            actual_dst_nodata = internal_nodata
        else:
            src_for_warp = src
            actual_src_nodata = src_nodata
            actual_dst_nodata = src_nodata

    dst = np.full(
        (target_grid["height"], target_grid["width"]),
        actual_dst_nodata,
        dtype=dst_dtype
    )

    reproject(
        source=src_for_warp,
        destination=dst,
        src_transform=src_transform,
        src_crs=src_crs,
        src_nodata=actual_src_nodata,
        dst_transform=target_grid["transform"],
        dst_crs=target_grid["crs"],
        dst_nodata=actual_dst_nodata,
        resampling=resampling
    )

    if actual_dst_nodata is not None:
        dst = np.asarray(dst, dtype=np.float32)
        dst[dst == actual_dst_nodata] = np.nan

    return dst.astype(dst_dtype)

def compute_hillshade(dem, transform, azimuth=315.0, altitude=45.0):
    """
    Compute a DEM hillshade on the working DEM/model grid.

    The hillshade is first computed in the DEM/model CRS. If a web basemap is
    used, the hillshade is later reprojected to EPSG:3857 together with the
    flood-depth/I_t layers and the outline.
    """
    dem = dem.astype(np.float32).copy()
    valid = np.isfinite(dem)

    if not np.any(valid):
        raise ValueError("DEM has no valid finite cells for hillshade.")

    fill_value = float(np.nanmedian(dem))
    dem_filled = np.where(valid, dem, fill_value)

    dx = abs(transform.a)
    dy = abs(transform.e)

    if dx <= 0 or dy <= 0:
        raise ValueError("Invalid DEM transform for hillshade. dx={}, dy={}".format(dx, dy))

    grad_y, grad_x = np.gradient(dem_filled, dy, dx)

    slope = np.pi / 2.0 - np.arctan(np.sqrt(grad_x * grad_x + grad_y * grad_y))
    aspect = np.arctan2(-grad_x, grad_y)

    azimuth_rad = np.deg2rad(azimuth)
    altitude_rad = np.deg2rad(altitude)

    shaded = (
        np.sin(altitude_rad) * np.sin(slope)
        + np.cos(altitude_rad) * np.cos(slope) * np.cos(azimuth_rad - aspect)
    )

    hillshade = 255.0 * (shaded + 1.0) / 2.0
    hillshade = np.clip(hillshade, 0.0, 255.0)
    hillshade[~valid] = np.nan

    return hillshade.astype(np.float32)

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
    return np.isfinite(dem).astype(np.uint8)


def extract_outline_geometries_from_mask(mask_uint8, transform):
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
# GEOTIFF EXPORT UTILITIES
# ============================================================
def prepare_array_for_geotiff(native_arr, dem_meta, target_grid, resampling=Resampling.bilinear):
    mode = str(GEOTIFF_EXPORT_CRS_MODE).lower()

    if mode == "dem":
        return native_arr, dem_meta["crs"], dem_meta["transform"]

    if mode == "web_mercator":
        arr_proj = reproject_array_to_target(
            src_array=native_arr,
            src_crs=dem_meta["crs"],
            src_transform=dem_meta["transform"],
            target_grid=target_grid,
            src_nodata=np.nan,
            resampling=resampling
        )
        return arr_proj, target_grid["crs"], target_grid["transform"]

    raise ValueError("Invalid GEOTIFF_EXPORT_CRS_MODE: {}. Use 'dem' or 'web_mercator'.".format(
        GEOTIFF_EXPORT_CRS_MODE
    ))


def write_geotiff(output_tif, arr, crs, transform, units="", description=""):
    output_tif = Path(output_tif)
    output_tif.parent.mkdir(parents=True, exist_ok=True)

    arr = np.asarray(arr, dtype=np.float32)
    out = np.where(np.isfinite(arr), arr, GEOTIFF_NODATA_VALUE).astype(np.float32)

    profile = {
        "driver": "GTiff",
        "height": out.shape[0],
        "width": out.shape[1],
        "count": 1,
        "dtype": "float32",
        "crs": crs,
        "transform": transform,
        "nodata": GEOTIFF_NODATA_VALUE,
        "compress": GEOTIFF_COMPRESS,
    }

    with rasterio.open(str(output_tif), "w", **profile) as dst:
        dst.write(out, 1)
        tags = {}
        if units:
            tags["units"] = units
        if description:
            tags["description"] = description
        if tags:
            dst.update_tags(**tags)

    print("Saved GeoTIFF: {}".format(output_tif))


def update_nanmax(max_arr, new_arr):
    if max_arr is None:
        return new_arr.copy().astype(np.float32)
    return np.fmax(max_arr, new_arr).astype(np.float32)


# ============================================================
# PLOTTING UTILITIES
# ============================================================
def build_norm_and_ticks(min_value, max_value, cmap_name, n_intervals, color_mode):
    if max_value is None:
        raise ValueError("max_value cannot be None when building a norm.")

    if max_value <= min_value:
        max_value = min_value + 1.0

    levels = np.linspace(float(min_value), float(max_value), int(n_intervals) + 1)

    if color_mode == "discrete":
        cmap = plt.get_cmap(cmap_name, int(n_intervals))
        norm = BoundaryNorm(levels, ncolors=cmap.N, clip=False)
        colorbar_ticks = levels

    elif color_mode == "continuous":
        cmap = plt.get_cmap(cmap_name)
        norm = Normalize(vmin=float(min_value), vmax=float(max_value), clip=False)
        colorbar_ticks = levels

    else:
        raise ValueError("Invalid COLOR_MODE: {}. Use 'continuous' or 'discrete'.".format(color_mode))

    return norm, cmap, colorbar_ticks


def add_plot_layers(ax, data_arr, data_extent, norm, cmap, colorbar_ticks,
                    plot_window, catchment_outline_lines, cbar_label,
                    hillshade_arr=None, hillshade_extent=None,
                    data_alpha=0.72, interpolation_mode="bilinear"):

    ax.set_xlim(plot_window[0], plot_window[1])
    ax.set_ylim(plot_window[2], plot_window[3])

    if USE_CONTEXTILY_BASEMAP:
        cx.add_basemap(
            ax,
            crs="EPSG:3857",
            source=BASEMAP_SOURCE,
            zoom=BASEMAP_ZOOM,
            attribution=BASEMAP_ATTRIBUTION,
            reset_extent=False,
            zorder=BASEMAP_ZORDER
        )

    if USE_HILLSHADE and hillshade_arr is not None:
        ax.imshow(
            hillshade_arr,
            extent=hillshade_extent,
            origin="upper",
            cmap="gray",
            vmin=0,
            vmax=255,
            alpha=HILLSHADE_ALPHA,
            zorder=HILLSHADE_ZORDER
        )

    data_im = ax.imshow(
        data_arr,
        extent=data_extent,
        origin="upper",
        cmap=cmap,
        norm=norm,
        alpha=data_alpha,
        interpolation=interpolation_mode,
        zorder=3
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
        data_im,
        ax=ax,
        fraction=COLORBAR_FRACTION,
        pad=COLORBAR_PAD,
        shrink=COLORBAR_SHRINK,
        aspect=COLORBAR_ASPECT,
        extend=COLORBAR_EXTEND,
        extendfrac=COLORBAR_EXTEND_FRAC
    )

    cbar.set_ticks(colorbar_ticks)
    cbar.set_label(cbar_label, fontsize=CBAR_LABEL_FONTSIZE)
    cbar.ax.tick_params(labelsize=CBAR_TICK_FONTSIZE)

    ax.set_xlabel(X_LABEL, fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel(Y_LABEL, fontsize=AXIS_LABEL_FONTSIZE)
    ax.tick_params(axis="both", labelsize=AXIS_TICK_FONTSIZE)


def save_map_png(output_png, data_arr, title_text, frame_text,
                 data_extent, norm, cmap, colorbar_ticks, plot_window,
                 catchment_outline_lines, cbar_label,
                 hillshade_arr=None, hillshade_extent=None,
                 data_alpha=0.72):

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)

    interpolation_mode = "bilinear" if COLOR_MODE == "continuous" else "nearest"

    add_plot_layers(
        ax=ax,
        data_arr=data_arr,
        data_extent=data_extent,
        norm=norm,
        cmap=cmap,
        colorbar_ticks=colorbar_ticks,
        plot_window=plot_window,
        catchment_outline_lines=catchment_outline_lines,
        cbar_label=cbar_label,
        hillshade_arr=hillshade_arr,
        hillshade_extent=hillshade_extent,
        data_alpha=data_alpha,
        interpolation_mode=interpolation_mode
    )

    ax.set_title(title_text, fontsize=TITLE_FONTSIZE)

    if frame_text is not None:
        ax.text(
            0.02,
            0.02,
            frame_text,
            transform=ax.transAxes,
            fontsize=FRAME_TEXT_FONTSIZE,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"),
            zorder=10
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
    print("Saved PNG: {}".format(output_png))


# ============================================================
# FRAME COUNT / I_t RANGE UTILITIES
# ============================================================
def count_total_depth_frames(chunk_files, dem_shape):
    total_frames = 0

    for i, fp in enumerate(chunk_files, start=1):
        loaded = load_mat_file(fp)
        maps = get_maps_struct(loaded)
        if maps is None:
            raise RuntimeError("Could not find Maps structure in {}".format(fp))

        depth_raw = extract_depth_3d_from_maps(maps)
        depth_std = standardize_cube_to_dem(depth_raw, dem_shape, variable_name="depth")
        total_frames += depth_std.shape[2]

        print("Counted depth chunk {}/{}: {} -> {} frames".format(
            i, len(chunk_files), fp.name, depth_std.shape[2]
        ))

    return total_frames


def scan_it_files_for_count_and_range(input_files, dem_shape):
    total_frames = 0
    global_min = np.inf
    global_max = -np.inf
    detected_path = None

    for i, fp in enumerate(input_files, start=1):
        loaded = load_mat_file(fp)
        it_cube, found_path = extract_it_3d_from_loaded_input(loaded, dem_shape)
        detected_path = found_path

        it_cube = np.where(
            np.isfinite(it_cube) & (it_cube > IT_DISPLAY_THRESHOLD_MM),
            it_cube,
            np.nan
        )

        total_frames += it_cube.shape[2]

        if np.any(np.isfinite(it_cube)):
            local_min = float(np.nanmin(it_cube))
            local_max = float(np.nanmax(it_cube))
            global_min = min(global_min, local_min)
            global_max = max(global_max, local_max)

        print("Scanned I_t chunk {}/{}: {} -> {} frames".format(
            i, len(input_files), fp.name, it_cube.shape[2]
        ))

    if not np.isfinite(global_min):
        global_min = IT_MIN_MM
    if not np.isfinite(global_max):
        global_max = IT_MIN_MM + 1.0

    return total_frames, global_min, global_max, detected_path


# ============================================================
# DEPTH PROCESSING
# ============================================================
def process_depth_maps(dem, dem_meta, target_grid, flood_extent, plot_window,
                       catchment_outline_lines, hillshade_proj):

    chunk_files = list_hydro_chunk_files(TEMP_DIR)
    print("Found {} hydro chunk files.".format(len(chunk_files)))

    depth_norm, depth_cmap, depth_ticks = build_norm_and_ticks(
        MIN_FLOOD_DEPTH_M,
        MAX_FLOOD_DEPTH_M,
        FLOOD_CMAP,
        N_COLOR_INTERVALS,
        COLOR_MODE
    )

    output_png_dir = Path(OUTPUT_FRAMES_DIR)
    output_png_dir.mkdir(parents=True, exist_ok=True)

    export_dir = Path(EXPORT_DATA_DIR)
    depth_tif_dir = export_dir / "Flood_Depth_TimeStep_TIFs"
    max_depth_dir = export_dir / "Maximum_Flood_Depth"

    dem_shape = dem.shape

    if TEST_SINGLE_FRAME:
        if EXPORT_MAX_FLOOD_DEPTH_TIF or EXPORT_EACH_FLOOD_DEPTH_TIF or EXPORT_MAX_FLOOD_DEPTH_PNG:
            print("WARNING: TEST_SINGLE_FRAME=True. Full depth GeoTIFF exports are skipped in test mode.")

        fp0 = chunk_files[0]
        loaded0 = load_mat_file(fp0)
        maps0 = get_maps_struct(loaded0)
        depth_3d_raw0 = extract_depth_3d_from_maps(maps0)
        depth_3d0 = standardize_cube_to_dem(depth_3d_raw0, dem_shape, variable_name="depth")
        frames_per_chunk = depth_3d0.shape[2]

        target_chunk_zero_based = TEST_FRAME_INDEX // frames_per_chunk
        local_frame_idx = TEST_FRAME_INDEX % frames_per_chunk

        if target_chunk_zero_based >= len(chunk_files):
            raise ValueError(
                "TEST_FRAME_INDEX={} is too large. With {} chunks and {} frames/chunk, "
                "maximum valid frame is {}.".format(
                    TEST_FRAME_INDEX,
                    len(chunk_files),
                    frames_per_chunk,
                    len(chunk_files) * frames_per_chunk - 1
                )
            )

        fp = chunk_files[target_chunk_zero_based]
        print("TEST depth mode: chunk={}, local_frame={}".format(fp.name, local_frame_idx))

        loaded = load_mat_file(fp)
        maps = get_maps_struct(loaded)
        depth_3d_raw = extract_depth_3d_from_maps(maps)
        depth_3d = standardize_cube_to_dem(depth_3d_raw, dem_shape, variable_name="depth")

        depth_m = depth_3d[:, :, local_frame_idx].astype(np.float32)
        depth_m[~np.isfinite(depth_m)] = np.nan
        depth_m[depth_m <= DEPTH_THRESHOLD_M] = np.nan

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

        time_str = format_time_label(TEST_FRAME_INDEX)
        save_map_png(
            output_png=TEST_OUTPUT_PNG,
            data_arr=flood_proj,
            title_text=TITLE_TEMPLATE.format(timestamp=time_str),
            frame_text=FRAME_LABEL_TEMPLATE.format(current=TEST_FRAME_INDEX + 1),
            data_extent=flood_extent,
            norm=depth_norm,
            cmap=depth_cmap,
            colorbar_ticks=depth_ticks,
            plot_window=plot_window,
            catchment_outline_lines=catchment_outline_lines,
            cbar_label=CBAR_LABEL,
            hillshade_arr=hillshade_proj,
            hillshade_extent=flood_extent,
            data_alpha=FLOOD_ALPHA
        )
        return

    if READ_MAT_FILES_FOR_FRAME_COUNTS:
        total_frames = count_total_depth_frames(chunk_files, dem_shape)
        print("Total depth frames: {}".format(total_frames))
    else:
        total_frames = None
        print("Skipping depth pre-count to avoid reading save_map_hydro_*.mat twice.")

    max_depth_m = None
    global_frame_idx = 0

    for chunk_idx, fp in enumerate(chunk_files, start=1):
        loaded = load_mat_file(fp)
        maps = get_maps_struct(loaded)

        if maps is None:
            raise RuntimeError("Could not find Maps structure in {}".format(fp))

        depth_3d_raw = extract_depth_3d_from_maps(maps)
        print("Raw depth cube shape from MAT: {}".format(depth_3d_raw.shape))

        depth_3d = standardize_cube_to_dem(depth_3d_raw, dem_shape, variable_name="depth")
        print("Standardized depth cube shape: {}".format(depth_3d.shape))

        n_local = depth_3d.shape[2]
        print("Rendering/exporting depth chunk {}/{}: {} ({} frames)".format(
            chunk_idx, len(chunk_files), fp.name, n_local
        ))

        for j in range(n_local):
            depth_m = depth_3d[:, :, j].astype(np.float32)
            depth_m[~np.isfinite(depth_m)] = np.nan
            depth_m[depth_m <= DEPTH_THRESHOLD_M] = np.nan

            if depth_m.shape != dem_shape:
                raise ValueError(
                    "Depth slice shape {} does not match DEM shape {} after standardization.".format(
                        depth_m.shape, dem_shape
                    )
                )

            max_depth_m = update_nanmax(max_depth_m, depth_m)

            if EXPORT_EACH_FLOOD_DEPTH_TIF:
                tif_arr, tif_crs, tif_transform = prepare_array_for_geotiff(
                    depth_m,
                    dem_meta,
                    target_grid,
                    resampling=Resampling.bilinear
                )
                output_tif = depth_tif_dir / "flood_depth_m_{:04d}.tif".format(global_frame_idx)
                write_geotiff(
                    output_tif,
                    tif_arr,
                    tif_crs,
                    tif_transform,
                    units="m",
                    description="Flood depth at saved model frame {}".format(global_frame_idx)
                )

            if PROCESS_FLOOD_DEPTH_PNGS:
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

                time_str = format_time_label(global_frame_idx)
                output_png = output_png_dir / "frame_{:04d}.png".format(global_frame_idx)

                save_map_png(
                    output_png=output_png,
                    data_arr=flood_proj,
                    title_text=TITLE_TEMPLATE.format(timestamp=time_str),
                    frame_text=FRAME_LABEL_TEMPLATE.format(current=global_frame_idx + 1),
                    data_extent=flood_extent,
                    norm=depth_norm,
                    cmap=depth_cmap,
                    colorbar_ticks=depth_ticks,
                    plot_window=plot_window,
                    catchment_outline_lines=catchment_outline_lines,
                    cbar_label=CBAR_LABEL,
                    hillshade_arr=hillshade_proj,
                    hillshade_extent=flood_extent,
                    data_alpha=FLOOD_ALPHA
                )

            global_frame_idx += 1

            if total_frames is None:
                if (global_frame_idx % PROGRESS_EVERY_N_FRAMES == 0):
                    print("Processed depth {} frames".format(global_frame_idx))
            else:
                if (global_frame_idx % PROGRESS_EVERY_N_FRAMES == 0) or (global_frame_idx == total_frames):
                    print("Processed depth {}/{} frames".format(global_frame_idx, total_frames))

    if max_depth_m is not None:
        if EXPORT_MAX_FLOOD_DEPTH_TIF:
            tif_arr, tif_crs, tif_transform = prepare_array_for_geotiff(
                max_depth_m,
                dem_meta,
                target_grid,
                resampling=Resampling.bilinear
            )
            output_tif = max_depth_dir / "maximum_flood_depth_m.tif"
            write_geotiff(
                output_tif,
                tif_arr,
                tif_crs,
                tif_transform,
                units="m",
                description="Maximum flood depth over all saved model maps"
            )

        if EXPORT_MAX_FLOOD_DEPTH_PNG:
            max_depth_proj = reproject_array_to_target(
                src_array=max_depth_m,
                src_crs=dem_meta["crs"],
                src_transform=dem_meta["transform"],
                target_grid=target_grid,
                src_nodata=np.nan,
                resampling=Resampling.bilinear
            )
            max_depth_proj = np.where(
                np.isfinite(max_depth_proj) & (max_depth_proj >= DISPLAY_THRESHOLD_M),
                max_depth_proj,
                np.nan
            )
            save_map_png(
                output_png=max_depth_dir / "maximum_flood_depth_m.png",
                data_arr=max_depth_proj,
                title_text=MAX_DEPTH_TITLE_TEMPLATE,
                frame_text=None,
                data_extent=flood_extent,
                norm=depth_norm,
                cmap=depth_cmap,
                colorbar_ticks=depth_ticks,
                plot_window=plot_window,
                catchment_outline_lines=catchment_outline_lines,
                cbar_label=CBAR_LABEL,
                hillshade_arr=hillshade_proj,
                hillshade_extent=flood_extent,
                data_alpha=FLOOD_ALPHA
            )

    print("Depth processing complete.")


# ============================================================
# I_t PROCESSING
# ============================================================
def process_it_maps(dem, dem_meta, target_grid, flood_extent, plot_window,
                    catchment_outline_lines, hillshade_proj):

    input_files = list_input_map_files(INPUT_MAPS_DIR)
    print("Found {} input map files.".format(len(input_files)))

    dem_shape = dem.shape

    if SCAN_IT_FOR_AUTO_COLORBAR:
        total_it_frames, detected_min, detected_max, detected_path = scan_it_files_for_count_and_range(
            input_files,
            dem_shape
        )
    else:
        total_it_frames = None
        detected_min = np.nan
        detected_max = np.nan
        detected_path = None
        print("Skipping I_t pre-scan to avoid reading save_input_maps_*.mat twice.")

    if IT_MAX_MM is None:
        if SCAN_IT_FOR_AUTO_COLORBAR and np.isfinite(detected_max):
            it_max_for_plot = detected_max
        else:
            it_max_for_plot = float(IT_DEFAULT_MAX_MM)
            if PLOT_EACH_IT_PNG:
                print("Using IT_DEFAULT_MAX_MM={} mm for per-frame I_t PNG colorbar.".format(IT_DEFAULT_MAX_MM))
                print("Set IT_MAX_MM manually or SCAN_IT_FOR_AUTO_COLORBAR=True if you need an auto colorbar before plotting every frame.")
    else:
        it_max_for_plot = float(IT_MAX_MM)

    if detected_path is not None:
        print("I_t detected path example: {}".format(detected_path))
    if np.isfinite(detected_min) and np.isfinite(detected_max):
        print("I_t detected finite range after threshold: min={:.4f} mm, max={:.4f} mm".format(
            float(detected_min), float(detected_max)
        ))
    print("I_t colorbar range for per-frame PNGs: {:.4f} to {:.4f} mm".format(
        float(IT_MIN_MM), float(it_max_for_plot)
    ))
    if total_it_frames is not None:
        print("Total I_t frames: {}".format(total_it_frames))

    it_norm, it_cmap, it_ticks = build_norm_and_ticks(
        IT_MIN_MM,
        it_max_for_plot,
        IT_CMAP,
        N_COLOR_INTERVALS,
        COLOR_MODE
    )

    output_png_dir = Path(OUTPUT_IT_FRAMES_DIR)
    output_png_dir.mkdir(parents=True, exist_ok=True)

    export_dir = Path(EXPORT_DATA_DIR)
    it_tif_dir = export_dir / "I_t_TimeStep_TIFs"
    max_it_dir = export_dir / "Maximum_I_t"

    if TEST_SINGLE_IT_FRAME:
        if EXPORT_MAX_IT_TIF or EXPORT_EACH_IT_TIF or EXPORT_MAX_IT_PNG:
            print("WARNING: TEST_SINGLE_IT_FRAME=True. Full I_t GeoTIFF exports are skipped in test mode.")

        fp0 = input_files[0]
        loaded0 = load_mat_file(fp0)
        it_cube0, _ = extract_it_3d_from_loaded_input(loaded0, dem_shape)
        frames_per_chunk = it_cube0.shape[2]

        target_chunk_zero_based = TEST_IT_FRAME_INDEX // frames_per_chunk
        local_frame_idx = TEST_IT_FRAME_INDEX % frames_per_chunk

        if target_chunk_zero_based >= len(input_files):
            raise ValueError(
                "TEST_IT_FRAME_INDEX={} is too large. With {} chunks and {} frames/chunk, "
                "maximum valid frame is {}.".format(
                    TEST_IT_FRAME_INDEX,
                    len(input_files),
                    frames_per_chunk,
                    len(input_files) * frames_per_chunk - 1
                )
            )

        fp = input_files[target_chunk_zero_based]
        print("TEST I_t mode: chunk={}, local_frame={}".format(fp.name, local_frame_idx))
        loaded = load_mat_file(fp)
        it_cube, _ = extract_it_3d_from_loaded_input(loaded, dem_shape)

        it_mm = it_cube[:, :, local_frame_idx].astype(np.float32)
        it_mm[~np.isfinite(it_mm)] = np.nan
        it_mm[it_mm <= IT_DISPLAY_THRESHOLD_MM] = np.nan

        it_proj = reproject_array_to_target(
            src_array=it_mm,
            src_crs=dem_meta["crs"],
            src_transform=dem_meta["transform"],
            target_grid=target_grid,
            src_nodata=np.nan,
            resampling=Resampling.bilinear
        )

        time_str = format_time_label(TEST_IT_FRAME_INDEX)
        save_map_png(
            output_png=TEST_OUTPUT_IT_PNG,
            data_arr=it_proj,
            title_text=IT_TITLE_TEMPLATE.format(timestamp=time_str),
            frame_text=FRAME_LABEL_TEMPLATE.format(current=TEST_IT_FRAME_INDEX + 1),
            data_extent=flood_extent,
            norm=it_norm,
            cmap=it_cmap,
            colorbar_ticks=it_ticks,
            plot_window=plot_window,
            catchment_outline_lines=catchment_outline_lines,
            cbar_label=IT_CBAR_LABEL,
            hillshade_arr=hillshade_proj,
            hillshade_extent=flood_extent,
            data_alpha=IT_ALPHA
        )
        return

    max_it_mm = None
    global_frame_idx = 0
    detected_path_runtime = None

    for chunk_idx, fp in enumerate(input_files, start=1):
        loaded = load_mat_file(fp)
        it_cube, found_path_runtime = extract_it_3d_from_loaded_input(loaded, dem_shape)
        if detected_path_runtime is None:
            detected_path_runtime = found_path_runtime
            print("I_t detected path example: {}".format(detected_path_runtime))

        n_local = it_cube.shape[2]
        print("Rendering/exporting I_t chunk {}/{}: {} ({} frames)".format(
            chunk_idx, len(input_files), fp.name, n_local
        ))

        for j in range(n_local):
            it_mm = it_cube[:, :, j].astype(np.float32)
            it_mm[~np.isfinite(it_mm)] = np.nan
            it_mm[it_mm <= IT_DISPLAY_THRESHOLD_MM] = np.nan

            if it_mm.shape != dem_shape:
                raise ValueError(
                    "I_t slice shape {} does not match DEM shape {} after standardization.".format(
                        it_mm.shape, dem_shape
                    )
                )

            max_it_mm = update_nanmax(max_it_mm, it_mm)

            if EXPORT_EACH_IT_TIF:
                tif_arr, tif_crs, tif_transform = prepare_array_for_geotiff(
                    it_mm,
                    dem_meta,
                    target_grid,
                    resampling=Resampling.bilinear
                )
                output_tif = it_tif_dir / "I_t_mm_{:04d}.tif".format(global_frame_idx)
                write_geotiff(
                    output_tif,
                    tif_arr,
                    tif_crs,
                    tif_transform,
                    units="mm",
                    description="I_t at saved input frame {}".format(global_frame_idx)
                )

            if PLOT_EACH_IT_PNG:
                it_proj = reproject_array_to_target(
                    src_array=it_mm,
                    src_crs=dem_meta["crs"],
                    src_transform=dem_meta["transform"],
                    target_grid=target_grid,
                    src_nodata=np.nan,
                    resampling=Resampling.bilinear
                )

                it_proj = np.where(
                    np.isfinite(it_proj) & (it_proj >= IT_DISPLAY_THRESHOLD_MM),
                    it_proj,
                    np.nan
                )

                time_str = format_time_label(global_frame_idx)
                output_png = output_png_dir / "I_t_frame_{:04d}.png".format(global_frame_idx)

                save_map_png(
                    output_png=output_png,
                    data_arr=it_proj,
                    title_text=IT_TITLE_TEMPLATE.format(timestamp=time_str),
                    frame_text=FRAME_LABEL_TEMPLATE.format(current=global_frame_idx + 1),
                    data_extent=flood_extent,
                    norm=it_norm,
                    cmap=it_cmap,
                    colorbar_ticks=it_ticks,
                    plot_window=plot_window,
                    catchment_outline_lines=catchment_outline_lines,
                    cbar_label=IT_CBAR_LABEL,
                    hillshade_arr=hillshade_proj,
                    hillshade_extent=flood_extent,
                    data_alpha=IT_ALPHA
                )

            global_frame_idx += 1

            if total_it_frames is None:
                if (global_frame_idx % PROGRESS_EVERY_N_FRAMES == 0):
                    print("Processed I_t {} frames".format(global_frame_idx))
            else:
                if (global_frame_idx % PROGRESS_EVERY_N_FRAMES == 0) or (global_frame_idx == total_it_frames):
                    print("Processed I_t {}/{} frames".format(global_frame_idx, total_it_frames))

    if max_it_mm is not None:
        if EXPORT_MAX_IT_TIF:
            tif_arr, tif_crs, tif_transform = prepare_array_for_geotiff(
                max_it_mm,
                dem_meta,
                target_grid,
                resampling=Resampling.bilinear
            )
            output_tif = max_it_dir / "maximum_I_t_mm.tif"
            write_geotiff(
                output_tif,
                tif_arr,
                tif_crs,
                tif_transform,
                units="mm",
                description="Maximum I_t over all saved input maps"
            )

        if EXPORT_MAX_IT_PNG:
            if IT_MAX_MM is None:
                finite_max_it = max_it_mm[np.isfinite(max_it_mm)]
                if finite_max_it.size > 0:
                    max_it_png_vmax = float(np.nanmax(finite_max_it))
                else:
                    max_it_png_vmax = float(IT_DEFAULT_MAX_MM)
            else:
                max_it_png_vmax = float(IT_MAX_MM)

            max_it_norm, max_it_cmap, max_it_ticks = build_norm_and_ticks(
                IT_MIN_MM,
                max_it_png_vmax,
                IT_CMAP,
                N_COLOR_INTERVALS,
                COLOR_MODE
            )

            max_it_proj = reproject_array_to_target(
                src_array=max_it_mm,
                src_crs=dem_meta["crs"],
                src_transform=dem_meta["transform"],
                target_grid=target_grid,
                src_nodata=np.nan,
                resampling=Resampling.bilinear
            )
            max_it_proj = np.where(
                np.isfinite(max_it_proj) & (max_it_proj >= IT_DISPLAY_THRESHOLD_MM),
                max_it_proj,
                np.nan
            )
            save_map_png(
                output_png=max_it_dir / "maximum_I_t_mm.png",
                data_arr=max_it_proj,
                title_text=IT_MAX_TITLE_TEMPLATE,
                frame_text=None,
                data_extent=flood_extent,
                norm=max_it_norm,
                cmap=max_it_cmap,
                colorbar_ticks=max_it_ticks,
                plot_window=plot_window,
                catchment_outline_lines=catchment_outline_lines,
                cbar_label=IT_CBAR_LABEL,
                hillshade_arr=hillshade_proj,
                hillshade_extent=flood_extent,
                data_alpha=IT_ALPHA
            )

    print("I_t processing complete.")


# ============================================================
# MAIN
# ============================================================
def main():
    print("Checking inputs...")
    print("TEMP_DIR: {}".format(TEMP_DIR))
    print("INPUT_MAPS_DIR: {}".format(INPUT_MAPS_DIR))
    print("DEM_TIF: {}".format(DEM_TIF))
    print("OUTPUT_FRAMES_DIR: {}".format(OUTPUT_FRAMES_DIR))
    print("OUTPUT_IT_FRAMES_DIR: {}".format(OUTPUT_IT_FRAMES_DIR))
    print("EXPORT_DATA_DIR: {}".format(EXPORT_DATA_DIR))
    print("GeoTIFF CRS mode: {}".format(GEOTIFF_EXPORT_CRS_MODE))
    print("Manual max flood depth: {} m".format(MAX_FLOOD_DEPTH_M))
    print("Flood colormap: {}".format(FLOOD_CMAP))
    print("I_t colormap: {}".format(IT_CMAP))
    print("Number of color intervals: {}".format(N_COLOR_INTERVALS))

    if not Path(TEMP_DIR).exists():
        raise FileNotFoundError("TEMP_DIR does not exist: {}".format(TEMP_DIR))
    if not Path(DEM_TIF).exists():
        raise FileNotFoundError("DEM_TIF does not exist: {}".format(DEM_TIF))
    if PROCESS_IT_MAPS and not Path(INPUT_MAPS_DIR).exists():
        raise FileNotFoundError("INPUT_MAPS_DIR does not exist: {}".format(INPUT_MAPS_DIR))

    Path(EXPORT_DATA_DIR).mkdir(parents=True, exist_ok=True)

    setup_font(FONT_PATH, FALLBACK_FONT_FAMILY)

    dem, dem_meta = load_dem_single_band(DEM_TIF)

    if dem_meta["crs"] is None:
        raise ValueError("DEM has no CRS defined")

    print("Original DEM shape: {}".format(dem.shape))
    print("DEM native bounds: {}".format(dem_meta["bounds"]))
    print("DEM CRS: {}".format(dem_meta["crs"]))

    # ------------------------------------------------------------------
    # IMPORTANT GRID ALIGNMENT FIX
    # ------------------------------------------------------------------
    # The saved HydroPol2D output maps define the authoritative matrix size.
    # If the DEM raster has a slightly different resolution/number of rows or
    # columns, we resample the DEM to match the output-map grid BEFORE building
    # the hillshade, mask, outline, plotting grid, and GeoTIFF metadata.
    model_shape = infer_model_grid_shape(reference_dem_shape=dem.shape)
    dem, dem_meta = resample_dem_to_model_grid_if_needed(dem, dem_meta, model_shape)

    dem_shape = dem.shape
    print("Working DEM/model grid shape: {}".format(dem_shape))

    # ------------------------------------------------------------------
    # CRS LOGIC
    # ------------------------------------------------------------------
    # Native/model grid: DEM CRS. Flood-depth and I_t arrays are assumed to
    # be on this grid after the DEM has been resampled to the model matrix size.
    #
    # Plot grid: EPSG:3857. This is required for the contextily web basemap.
    # Every visual layer shown in the PNGs is therefore reprojected to EPSG:3857:
    #   - flood depth / I_t
    #   - hillshade
    #   - DEM valid-data outline
    #   - contextily basemap
    #
    # GeoTIFF exports: controlled by GEOTIFF_EXPORT_CRS_MODE.
    # Default is "dem", meaning GeoTIFFs are written in the working DEM/model CRS.
    target_grid = build_target_grid_from_dem(dem_meta, "EPSG:3857")

    print("Working native/model CRS for arrays and default GeoTIFFs: {}".format(dem_meta["crs"]))
    print("PNG plotting CRS: {}".format(target_grid["crs"]))
    if USE_CONTEXTILY_BASEMAP and str(target_grid["crs"]).upper() != "EPSG:3857":
        raise ValueError("Contextily basemap requires the plotting grid to be EPSG:3857.")

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

    print("Projected plotting extent EPSG:3857: {}".format(flood_extent))
    print("Plot window: {}".format(plot_window))

    dem_valid_mask = build_dem_valid_mask(dem, dem_meta)

    dem_mask_geoms = extract_outline_geometries_from_mask(
        dem_valid_mask,
        dem_meta["transform"]
    )

    dem_mask_geoms_plot_crs = reproject_geometries(
        dem_mask_geoms,
        dem_meta["crs"],
        target_grid["crs"]
    )

    catchment_outline_lines = geoms_to_plot_lines(dem_mask_geoms_plot_crs)

    print("Extracted {} outline polygon(s) from DEM valid-data mask.".format(
        len(catchment_outline_lines)
    ))

    hillshade_proj = None

    if USE_HILLSHADE:
        print("Computing DEM hillshade...")

        hillshade_native = compute_hillshade(
            dem,
            dem_meta["transform"],
            azimuth=HILLSHADE_AZIMUTH,
            altitude=HILLSHADE_ALTITUDE
        )

        hillshade_proj = reproject_array_to_target(
            src_array=hillshade_native,
            src_crs=dem_meta["crs"],
            src_transform=dem_meta["transform"],
            target_grid=target_grid,
            src_nodata=np.nan,
            resampling=Resampling.bilinear
        )

        if np.any(np.isfinite(hillshade_proj)):
            print(
                "Hillshade projected to plot CRS. min={:.2f}, max={:.2f}".format(
                    float(np.nanmin(hillshade_proj)),
                    float(np.nanmax(hillshade_proj))
                )
            )
        else:
            print("WARNING: projected hillshade is fully NaN. Disabling hillshade layer.")
            hillshade_proj = None

    if MAX_FLOOD_DEPTH_M <= MIN_FLOOD_DEPTH_M:
        raise ValueError("MAX_FLOOD_DEPTH_M must be greater than MIN_FLOOD_DEPTH_M.")

    if PROCESS_FLOOD_DEPTH_PNGS or EXPORT_MAX_FLOOD_DEPTH_TIF or EXPORT_EACH_FLOOD_DEPTH_TIF or EXPORT_MAX_FLOOD_DEPTH_PNG:
        process_depth_maps(
            dem=dem,
            dem_meta=dem_meta,
            target_grid=target_grid,
            flood_extent=flood_extent,
            plot_window=plot_window,
            catchment_outline_lines=catchment_outline_lines,
            hillshade_proj=hillshade_proj
        )

    if PROCESS_IT_MAPS:
        process_it_maps(
            dem=dem,
            dem_meta=dem_meta,
            target_grid=target_grid,
            flood_extent=flood_extent,
            plot_window=plot_window,
            catchment_outline_lines=catchment_outline_lines,
            hillshade_proj=hillshade_proj
        )

    print("Done.")
    print("PNG flood frames: {}".format(OUTPUT_FRAMES_DIR))
    print("PNG I_t frames: {}".format(OUTPUT_IT_FRAMES_DIR))
    print("GeoTIFF/data exports: {}".format(EXPORT_DATA_DIR))

# ============================================================
# DYNAMIC MULTI-PANEL ANIMATION CLI
# ============================================================
def _resolve_modeling_results(results_dir):
    results_dir = Path(results_dir)
    if (results_dir / "Modeling_Results").exists():
        return results_dir / "Modeling_Results"
    return results_dir


def _parse_time_hours_from_raster_name(path):
    m = re.search(r"_t_([0-9.+\-eE]+)_h\.tif$", path.name)
    if m is None:
        return None
    try:
        return float(m.group(1))
    except Exception:
        return None


def _read_single_band(path):
    with rasterio.open(str(path)) as src:
        arr = src.read(1).astype(np.float32)
        nodata = src.nodata
        if nodata is not None:
            arr = np.where(arr == nodata, np.nan, arr)
    arr[~np.isfinite(arr)] = np.nan
    return arr


def _read_single_band_with_grid(path):
    with rasterio.open(str(path)) as src:
        arr = src.read(1).astype(np.float32)
        nodata = src.nodata
        if nodata is not None:
            arr = np.where(arr == nodata, np.nan, arr)
        grid = {
            "crs": src.crs,
            "transform": src.transform,
            "width": src.width,
            "height": src.height,
        }
    arr[~np.isfinite(arr)] = np.nan
    return arr, grid


def _safe_percentile(arrays, percentile, fallback):
    values = []
    for arr in arrays:
        data = np.asarray(arr, dtype=np.float32)
        data = data[np.isfinite(data) & (data > 0)]
        if data.size:
            if data.size > 250000:
                step = int(math.ceil(data.size / 250000.0))
                data = data[::step]
            values.append(data)
    if not values:
        return fallback
    merged = np.concatenate(values)
    vmax = float(np.nanpercentile(merged, percentile))
    if not np.isfinite(vmax) or vmax <= 0:
        return fallback
    return max(vmax, fallback)


def _discover_raster_series(modeling_results):
    specs = [
        {
            "key": "depth",
            "title": "Flood Depth",
            "units": "m",
            "folder": "Rasters_Water_Depths",
            "glob": "Flood_Depths_t_*_h.tif",
            "cmap": "Blues",
            "fallback_vmax": 0.25,
        },
        {
            "key": "velocity",
            "title": "Velocity",
            "units": "m/s",
            "folder": "Rasters_Velocity",
            "glob": "Velocity_t_*_h.tif",
            "cmap": "viridis",
            "fallback_vmax": 0.5,
        },
        {
            "key": "hazard",
            "title": "Hazard Metric",
            "units": "m2/s",
            "folder": "Rasters_Hazard",
            "glob": "Flood_Hazard_t_*_h.tif",
            "cmap": "inferno",
            "fallback_vmax": 0.25,
        },
        {
            "key": "infiltration",
            "title": "Infiltration",
            "units": "mm/h",
            "folder": "Rasters_Infiltration",
            "glob": "Infiltration_t_*_h.tif",
            "cmap": "YlGnBu",
            "fallback_vmax": 25.0,
        },
    ]

    panels = []
    for spec in specs:
        folder = modeling_results / spec["folder"]
        if not folder.exists():
            continue

        files = []
        for fp in folder.glob(spec["glob"]):
            t_h = _parse_time_hours_from_raster_name(fp)
            if t_h is not None:
                files.append((t_h, fp))

        if not files:
            continue

        files = sorted(files, key=lambda x: x[0])
        times = [x[0] for x in files]
        paths = [x[1] for x in files]

        panels.append({
            "source": "raster",
            "key": spec["key"],
            "title": spec["title"],
            "units": spec["units"],
            "times_h": times,
            "paths": paths,
            "cmap": spec["cmap"],
            "fallback_vmax": spec["fallback_vmax"],
        })

    return panels


def _extract_temp_panel_arrays(run_root, dem_path, map_step_hours):
    temp_dir = Path(run_root) / "Temporary_Files"
    if not temp_dir.exists():
        temp_dir = Path(run_root)
    chunk_files = sorted(temp_dir.glob("save_map_hydro_*.mat"),
                         key=lambda p: _extract_number_from_filename(p, r"save_map_hydro_(\d+)\.mat$"))
    if not chunk_files:
        return []

    dem, dem_meta = load_dem_single_band(dem_path)
    ref_shape = dem.shape

    specs = [
        ("depth", "Flood Depth", "m", "d", "Blues", 0.25, 1.0 / 1000.0),
        ("velocity", "Velocity", "m/s", "velocity", "viridis", 0.5, 1.0),
        ("hazard", "Hazard Metric", "m2/s", "hazard_dv", "inferno", 0.25, 1.0),
        ("infiltration", "Infiltration", "mm/h", "f", "YlGnBu", 25.0, 1.0),
    ]

    loaded_chunks = []
    target_shape = None
    for fp in chunk_files:
        loaded = load_mat_file(fp)
        maps = get_maps_struct(loaded)
        if maps is None:
            continue
        hydro = maybe_getattr(maps, "Hydro", None)
        if hydro is None:
            continue
        loaded_chunks.append(hydro)
        if target_shape is None:
            raw_d = maybe_getattr(hydro, "d", None)
            if raw_d is not None:
                raw_d = _numeric_to_3d_cube(raw_d, "Maps.Hydro.d")
                target_shape = infer_model_grid_shape_from_cube(raw_d, ref_shape, "depth")

    if not loaded_chunks or target_shape is None:
        return []

    left, bottom, right, top = dem_meta["bounds"]
    model_rows, model_cols = target_shape
    model_transform = from_bounds(left, bottom, right, top, model_cols, model_rows)

    panels = []
    for key, title, units, field, cmap, fallback_vmax, scale in specs:
        cubes = []
        for hydro in loaded_chunks:
            raw = maybe_getattr(hydro, field, None)
            if raw is None:
                continue
            try:
                raw_cube = _numeric_to_3d_cube(raw, "Maps.Hydro.{}".format(field))
                cube = standardize_cube_to_grid(raw_cube, target_shape, variable_name=field)
                cube = cube.astype(np.float32) * np.float32(scale)
                cube[~np.isfinite(cube)] = np.nan
                cubes.append(cube)
            except Exception as exc:
                print("Skipping temporary field {}: {}".format(field, exc))
                cubes = []
                break
        if not cubes:
            continue

        cube = np.concatenate(cubes, axis=2)
        n_frames = cube.shape[2]
        times_h = [i * float(map_step_hours) for i in range(n_frames)]
        panels.append({
            "source": "array",
            "key": key,
            "title": title,
            "units": units,
            "times_h": times_h,
            "cube": cube,
            "crs": dem_meta["crs"],
            "transform": model_transform,
            "cmap": cmap,
            "fallback_vmax": fallback_vmax,
        })

    return panels


def _nearest_panel_index(panel, time_h):
    times = np.asarray(panel["times_h"], dtype=float)
    if times.size == 0:
        return None
    return int(np.nanargmin(np.abs(times - float(time_h))))


def _panel_array(panel, idx):
    if idx is None:
        return None
    if panel["source"] == "raster":
        return _read_single_band(panel["paths"][idx])
    return panel["cube"][:, :, idx].astype(np.float32)


def _panel_frame_native(panel, idx):
    if idx is None:
        return None, None
    if panel["source"] == "raster":
        return _read_single_band_with_grid(panel["paths"][idx])
    arr = panel["cube"][:, :, idx].astype(np.float32)
    grid = {
        "crs": panel["crs"],
        "transform": panel["transform"],
        "width": arr.shape[1],
        "height": arr.shape[0],
    }
    return arr, grid


def _compute_panel_scales(panels, scale_percentile):
    for panel in panels:
        samples = []
        n = len(panel["times_h"])
        if n <= 0:
            panel["vmax"] = panel["fallback_vmax"]
            continue
        sample_indices = sorted(set([0, n // 4, n // 2, (3 * n) // 4, n - 1]))
        for idx in sample_indices:
            arr = _panel_array(panel, idx)
            if arr is not None:
                samples.append(arr)
        panel["vmax"] = _safe_percentile(samples, scale_percentile, panel["fallback_vmax"])
        print("{} scale vmax ({:.1f}th percentile sample): {:.4g} {}".format(
            panel["title"], scale_percentile, panel["vmax"], panel["units"]
        ))


def _set_panel_image(ax, arr, panel):
    data = np.asarray(arr, dtype=np.float32)
    data[~np.isfinite(data)] = np.nan
    data[data <= 0] = np.nan
    im = ax.imshow(
        data,
        origin="upper",
        cmap=panel["cmap"],
        vmin=0,
        vmax=panel["vmax"],
        interpolation="nearest",
    )
    ax.set_title("{} ({})".format(panel["title"], panel["units"]), fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor("#f2f2f2")
    return im


def _default_dem_for_results(modeling_results, run_root):
    candidates = [
        modeling_results / "Rasters_Static" / "DEM_resampled.tif",
        run_root / "Static" / "DEM.tif",
        run_root / "DEM.tif",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _project_panel_frame(panel, idx, target_grid):
    arr, grid = _panel_frame_native(panel, idx)
    if arr is None or grid is None:
        return None
    arr = np.asarray(arr, dtype=np.float32)
    arr[~np.isfinite(arr)] = np.nan
    projected = reproject_array_to_target(
        src_array=arr,
        src_crs=grid["crs"],
        src_transform=grid["transform"],
        target_grid=target_grid,
        src_nodata=np.nan,
        resampling=Resampling.bilinear,
        dst_dtype=np.float32,
    )
    projected = np.asarray(projected, dtype=np.float32)
    projected[~np.isfinite(projected)] = np.nan
    projected[projected <= 0] = np.nan
    return projected


def _apply_wet_display_mask(arr, panel, depth_arr, wet_threshold_m=0.01):
    """Mask dry-cell velocity and hazard displays using the paired depth map."""
    if arr is None:
        return None
    if depth_arr is None:
        return arr
    if panel.get("key") not in ("velocity", "hazard"):
        return arr
    masked = np.asarray(arr, dtype=np.float32).copy()
    wet = np.asarray(depth_arr, dtype=np.float32) >= wet_threshold_m
    masked[~wet] = np.nan
    return masked


def _prepare_website_base_layers(dem_path):
    dem, dem_meta = load_dem_single_band(dem_path)
    target_grid = build_target_grid_from_dem(dem_meta, "EPSG:3857")
    extent = extent_from_transform(
        target_grid["transform"],
        target_grid["width"],
        target_grid["height"]
    )
    plot_window = compute_plot_window_from_extent(
        extent,
        MAP_PAD_FRACTION_X,
        MAP_PAD_FRACTION_Y
    )

    hillshade_proj = None
    try:
        hillshade_native = compute_hillshade(
            dem,
            dem_meta["transform"],
            azimuth=HILLSHADE_AZIMUTH,
            altitude=HILLSHADE_ALTITUDE
        )
        hillshade_proj = reproject_array_to_target(
            src_array=hillshade_native,
            src_crs=dem_meta["crs"],
            src_transform=dem_meta["transform"],
            target_grid=target_grid,
            src_nodata=np.nan,
            resampling=Resampling.bilinear,
            dst_dtype=np.float32,
        )
    except Exception as exc:
        print("WARNING: could not prepare hillshade: {}".format(exc))

    outline_lines = []
    try:
        dem_valid_mask = build_dem_valid_mask(dem, dem_meta)
        dem_mask_geoms = extract_outline_geometries_from_mask(
            dem_valid_mask,
            dem_meta["transform"]
        )
        dem_mask_geoms_plot_crs = reproject_geometries(
            dem_mask_geoms,
            dem_meta["crs"],
            target_grid["crs"]
        )
        outline_lines = geoms_to_plot_lines(dem_mask_geoms_plot_crs)
    except Exception as exc:
        print("WARNING: could not prepare catchment outline: {}".format(exc))

    return {
        "dem": dem,
        "dem_meta": dem_meta,
        "target_grid": target_grid,
        "extent": extent,
        "plot_window": plot_window,
        "hillshade": hillshade_proj,
        "outline_lines": outline_lines,
    }


def generate_dynamic_panel_animation(results_dir, output_path=None, dem_path=None,
                                     map_step_hours=0.25, fps=3, dpi=140,
                                     scale_percentile=99.0, max_frames=None,
                                     title=None):
    results_dir = Path(results_dir)
    modeling_results = _resolve_modeling_results(results_dir)
    run_root = modeling_results.parent if modeling_results.name == "Modeling_Results" else results_dir

    panels = _discover_raster_series(modeling_results)

    if not panels:
        if dem_path is None:
            dem_path = _default_dem_for_results(modeling_results, run_root)
        if dem_path is None:
            raise FileNotFoundError(
                "No dynamic raster folders found and no DEM was provided for temporary MAT fallback."
            )
        panels = _extract_temp_panel_arrays(run_root, dem_path, map_step_hours)

    if not panels:
        raise FileNotFoundError("No dynamic depth/velocity/hazard/infiltration data found in {}".format(results_dir))

    if dem_path is None:
        dem_path = _default_dem_for_results(modeling_results, run_root)
    if dem_path is None:
        raise FileNotFoundError("Could not find DEM for website-style animation base layers.")
    dem_path = Path(dem_path)

    _compute_panel_scales(panels, scale_percentile)

    # Prefer depth times as the animation clock; otherwise use the first panel.
    clock_panel = next((p for p in panels if p["key"] == "depth"), panels[0])
    frame_times = list(clock_panel["times_h"])
    if max_frames is not None and max_frames > 0:
        frame_times = frame_times[:int(max_frames)]

    if not frame_times:
        raise ValueError("No frame times available for animation.")

    if output_path is None:
        output_dir = modeling_results / "GIFs_MP4"
        output_dir.mkdir(parents=True, exist_ok=True)
        output_path = output_dir / "Dynamic_Flood_Panels.mp4"
    else:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

    poster_path = output_path.with_suffix(".png")
    frames_dir = output_path.parent / (output_path.stem + "_frames")
    frames_dir.mkdir(parents=True, exist_ok=True)

    base = _prepare_website_base_layers(dem_path)
    target_grid = base["target_grid"]
    extent = base["extent"]
    plot_window = base["plot_window"]
    hillshade_proj = base["hillshade"]
    outline_lines = base["outline_lines"]

    ncols = len(panels)
    fig_width = max(5.2 * ncols, 6.5)
    fig, axes = plt.subplots(1, ncols, figsize=(fig_width, 5.8), dpi=dpi, squeeze=False)
    axes = axes.ravel()
    artists = []

    first_time_h = frame_times[0]
    depth_panel = next((p for p in panels if p["key"] == "depth"), None)
    for ax, panel in zip(axes, panels):
        ax.set_xlim(plot_window[0], plot_window[1])
        ax.set_ylim(plot_window[2], plot_window[3])
        ax.set_aspect("equal", adjustable="box")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_facecolor("#111111")

        if USE_CONTEXTILY_BASEMAP:
            try:
                cx.add_basemap(
                    ax,
                    crs="EPSG:3857",
                    source=BASEMAP_SOURCE,
                    zoom=BASEMAP_ZOOM,
                    attribution=BASEMAP_ATTRIBUTION,
                    reset_extent=False,
                    zorder=BASEMAP_ZORDER
                )
            except Exception as exc:
                print("WARNING: basemap unavailable for panel '{}': {}".format(panel["title"], exc))

        if USE_HILLSHADE and hillshade_proj is not None:
            ax.imshow(
                hillshade_proj,
                extent=extent,
                origin="upper",
                cmap="gray",
                vmin=0,
                vmax=255,
                alpha=HILLSHADE_ALPHA,
                zorder=HILLSHADE_ZORDER,
            )

        idx0 = _nearest_panel_index(panel, first_time_h)
        arr0 = _project_panel_frame(panel, idx0, target_grid)
        if depth_panel is not None:
            depth_idx0 = _nearest_panel_index(depth_panel, first_time_h)
            depth0 = _project_panel_frame(depth_panel, depth_idx0, target_grid)
            arr0 = _apply_wet_display_mask(arr0, panel, depth0)
        cmap = plt.get_cmap(panel["cmap"]).copy()
        cmap.set_bad((0, 0, 0, 0))
        im = ax.imshow(
            arr0,
            extent=extent,
            origin="upper",
            cmap=cmap,
            vmin=0,
            vmax=panel["vmax"],
            alpha=0.78,
            interpolation="bilinear",
            zorder=3,
        )

        for xs, ys in outline_lines:
            ax.plot(xs, ys, color=OUTLINE_COLOR, linewidth=1.1, zorder=4)

        ax.set_title("{} ({})".format(panel["title"], panel["units"]), fontsize=11, color="#222222")
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.018)
        cbar.set_label(panel["units"], fontsize=9)
        cbar.ax.tick_params(labelsize=8)
        artists.append(im)

    def draw_frame(frame_number, save_png=False):
        time_h = frame_times[frame_number]
        depth_arr = None
        if depth_panel is not None:
            depth_idx = _nearest_panel_index(depth_panel, time_h)
            depth_arr = _project_panel_frame(depth_panel, depth_idx, target_grid)
        for im, panel in zip(artists, panels):
            idx = _nearest_panel_index(panel, time_h)
            arr = _project_panel_frame(panel, idx, target_grid)
            arr = _apply_wet_display_mask(arr, panel, depth_arr)
            im.set_data(arr)

        heading = title or "HydroPol2D Dynamic Flood Maps"
        fig.suptitle("{} | t = {:.2f} h".format(heading, time_h), fontsize=15, y=0.98)
        fig.tight_layout(rect=[0, 0, 1, 0.93])

        if save_png:
            frame_path = frames_dir / "frame_{:04d}.png".format(frame_number)
            fig.savefig(str(frame_path), dpi=dpi, facecolor="white")
            return frame_path
        return []

    draw_frame(0, save_png=False)
    fig.savefig(str(poster_path), dpi=dpi, facecolor="white")
    print("Saved poster PNG: {}".format(poster_path))

    if output_path.suffix.lower() == ".gif":
        writer_name = "pillow"
    else:
        writer_name = "ffmpeg"

    try:
        writer = mpl_animation.writers[writer_name](fps=fps)
        with writer.saving(fig, str(output_path), dpi=dpi):
            for frame_number in range(len(frame_times)):
                draw_frame(frame_number, save_png=False)
                writer.grab_frame(facecolor="white")
                if (frame_number + 1) % 5 == 0 or frame_number == len(frame_times) - 1:
                    print("Rendered frame {}/{}".format(frame_number + 1, len(frame_times)))
        print("Saved animation: {}".format(output_path))
    except Exception as exc:
        print("Animation writer failed: {}".format(exc))
        print("Falling back to PNG frames in {}".format(frames_dir))
        for frame_number in range(len(frame_times)):
            draw_frame(frame_number, save_png=True)
        print("Saved PNG frame sequence: {}".format(frames_dir))
    finally:
        plt.close(fig)

    return output_path


def dynamic_cli(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate a robust multi-panel HydroPol2D dynamic-map animation."
    )
    parser.add_argument("--results-dir", required=True,
                        help="Run root or Modeling_Results folder.")
    parser.add_argument("--output", default=None,
                        help="Output MP4/GIF path. Defaults to Modeling_Results/GIFs_MP4/Dynamic_Flood_Panels.mp4.")
    parser.add_argument("--dem", default=None,
                        help="DEM path used only when falling back to Temporary_Files/save_map_hydro_*.mat.")
    parser.add_argument("--map-step-hours", type=float, default=0.25,
                        help="Fallback MAT map timestep in hours.")
    parser.add_argument("--fps", type=float, default=3.0)
    parser.add_argument("--dpi", type=int, default=140)
    parser.add_argument("--scale-percentile", type=float, default=99.0)
    parser.add_argument("--max-frames", type=int, default=None)
    parser.add_argument("--title", default=None)
    args = parser.parse_args(argv)

    generate_dynamic_panel_animation(
        results_dir=args.results_dir,
        output_path=args.output,
        dem_path=args.dem,
        map_step_hours=args.map_step_hours,
        fps=args.fps,
        dpi=args.dpi,
        scale_percentile=args.scale_percentile,
        max_frames=args.max_frames,
        title=args.title,
    )


if __name__ == "__main__":
    import sys
    if "--results-dir" in sys.argv:
        dynamic_cli()
    else:
        main()
