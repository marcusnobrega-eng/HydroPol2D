"""Build Examples/Stanford/Static_10m_LiDAR/DEM.tif (done 2026-11-03).

Starting from the earlier 10 m 3DEP file of the example, which has no data in a 5.9 km2 western
headwater patch:
1. the gap is filled from today's USGS 3DEPElevation service at 10 m, plus the difference
   (old file - today's 3DEP) interpolated harmonically in from the gap edge, so the fill meets
   the old file without a step;
2. two 3DEP tile edges inside the gap (10-15 m steps) are inpainted harmonically;
3. 26 single-cell pits and spikes (more than 10 m from the median of their 8 neighbours, up to
   82 m deep; also in today's 3DEP) take that median;
4. the DEM is clipped to the 10 m catchment mask of Static_10m_LiDAR_ModelDEM.

Today's 3DEP for the same box and grid (EPSG:3310, 10 m, 1923 x 1602 cells):
  https://elevation.nationalmap.gov/arcgis/rest/services/3DEPElevation/ImageServer/exportImage?
  bbox=-204750,-75540,-185520,-59520&bboxSR=3310&imageSR=3310&size=1923,1602&format=tiff&
  pixelType=F32&noData=-9999&interpolation=RSP_BilinearInterpolation&f=image

  python fill_10m_lidar_gap.py OLD_DEM TODAYS_3DEP MODEL_DEM OUT
"""
import sys

import numpy as np
import rasterio
from scipy import ndimage, sparse
from scipy.sparse.linalg import spsolve

raw_path, dep_path, mask_path, out_path = sys.argv[1:5]
ND = -9999.0
SPIKE_M = 10.0

with rasterio.open(raw_path) as d:
    z = d.read(1).astype(np.float64)
    profile = d.profile
    T = d.transform
with rasterio.open(dep_path) as d:
    g = d.read(1).astype(np.float64)
    assert d.transform == T and d.shape == z.shape
with rasterio.open(mask_path) as d:  # same 10 m lattice, shifted by whole cells
    m = d.read(1) != d.nodata
    r_off = int(round((T.f - d.transform.f) / T.a))
    c_off = int(round((d.transform.c - T.c) / T.a))
    inside = np.zeros(z.shape, bool)
    inside[r_off:r_off + m.shape[0], c_off:c_off + m.shape[1]] = m

valid = z != ND
gap = inside & ~valid
print(f"catchment {inside.sum() / 1e4:.2f} km2, gap {gap.sum() / 1e4:.2f} km2")

# 1. Offset D = raw - 3DEP, known on valid cells; harmonic (Laplace) interpolation into the gap,
#    so the fill equals the raw DEM at the gap edge and follows 3DEP's shape inside.
D = np.where(valid, z - g, 0.0)
idx = -np.ones(z.shape, np.int64)
idx[gap] = np.arange(gap.sum())
rows, cols = np.nonzero(gap)
A = sparse.lil_matrix((len(rows), len(rows)))
b = np.zeros(len(rows))
for k, (r, c) in enumerate(zip(rows, cols)):
    n = 0
    for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        rr, cc = r + dr, c + dc
        if not (0 <= rr < z.shape[0] and 0 <= cc < z.shape[1]):
            continue
        if gap[rr, cc]:
            A[k, idx[rr, cc]] = -1.0
            n += 1
        elif valid[rr, cc]:
            b[k] += D[rr, cc]
            n += 1
        # neighbours that are no data outside the catchment: zero gradient (left out)
    A[k, k] = max(n, 1)
D_gap = spsolve(A.tocsr(), b)
filled = z.copy()
filled[gap] = g[gap] + D_gap
print(f"gap offset raw-3DEP interpolated: mean {D_gap.mean():+.2f} m, range {D_gap.min():+.2f} to {D_gap.max():+.2f} m")

# 1b. Today's 3DEP has two source (tile) edges inside the gap: an east-west one across rows
#     391-402 and a north-south one across columns 95-104 of this grid, with 10-15 m steps.
#     Those strips of the gap are inpainted harmonically from the terrain on both sides.
def harmonic_fill(values, region, known):
    """Solve Laplace's equation for values in region; known cells around it are fixed."""
    ids = -np.ones(values.shape, np.int64)
    ids[region] = np.arange(region.sum())
    rr, cc = np.nonzero(region)
    M = sparse.lil_matrix((len(rr), len(rr)))
    rhs = np.zeros(len(rr))
    for k, (r, c) in enumerate(zip(rr, cc)):
        n = 0
        for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            r2, c2 = r + dr, c + dc
            if not (0 <= r2 < values.shape[0] and 0 <= c2 < values.shape[1]):
                continue
            if region[r2, c2]:
                M[k, ids[r2, c2]] = -1.0
                n += 1
            elif known[r2, c2]:
                rhs[k] += values[r2, c2]
                n += 1
        M[k, k] = max(n, 1)
    out = values.copy()
    out[region] = spsolve(M.tocsr(), rhs)
    return out


row_index, col_index = np.indices(z.shape)
seam = gap & (((row_index >= 391) & (row_index <= 402)) | ((col_index >= 95) & (col_index <= 104)))
filled = harmonic_fill(filled, seam, inside & (filled != ND) & ~seam)
print(f"3DEP tile edges inside the gap: {seam.sum() / 1e4:.3f} km2 inpainted")

# 2. Single-cell pits and spikes: more than SPIKE_M from the 3x3 median of the other cells.
work = np.where(filled != ND, filled, np.nan)
pad = np.pad(work, 1, constant_values=np.nan)
stack = np.stack([pad[1 + dr:1 + dr + z.shape[0], 1 + dc:1 + dc + z.shape[1]]
                  for dr in (-1, 0, 1) for dc in (-1, 0, 1) if (dr, dc) != (0, 0)])
med = np.nanmedian(stack, axis=0)
count = np.sum(np.isfinite(stack), axis=0)
# Single-cell holes and spikes (also in today's 3DEP); they would trap or block flow.
spike = inside & valid & np.isfinite(work) & (count >= 6) & (np.abs(work - med) > SPIKE_M)
dev = (work - med)[spike]
print(f"single-cell pits/spikes > {SPIKE_M:g} m: {int(spike.sum())} ({int((dev < 0).sum())} pits, "
      f"{int((dev > 0).sum())} spikes), largest {np.abs(dev).max():.1f} m")
filled[spike] = med[spike]

# 3. Clip to the catchment.
filled[~inside] = ND
print(f"valid {np.sum(filled != ND) / 1e4:.2f} km2, range {filled[filled != ND].min():.2f} "
      f"to {filled[filled != ND].max():.2f} m")

# Seam check: steps between the gap's edge cells and their outside neighbours vs. elsewhere.
edge = gap & ndimage.binary_dilation(valid & inside)
steps, normal, interior = [], [], []
for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
    nb = np.roll(np.roll(filled, -dr, 0), -dc, 1)
    sel = edge & np.roll(np.roll(valid & inside, -dr, 0), -dc, 1)
    steps.append(np.abs(filled - nb)[sel])
    sel2 = valid & inside & np.roll(np.roll(valid & inside, -dr, 0), -dc, 1)
    normal.append(np.abs(filled - nb)[sel2])
    sel3 = gap & ~edge & np.roll(np.roll(gap, -dr, 0), -dc, 1)
    interior.append(np.abs(filled - nb)[sel3])
steps, normal, interior = np.concatenate(steps), np.concatenate(normal), np.concatenate(interior)
dy = np.abs(np.diff(filled, axis=0))
both = gap[1:] & gap[:-1]
per_row = [np.mean(dy[r][both[r]]) for r in range(dy.shape[0]) if both[r].sum() > 20]
print(f"largest mean step between two rows inside the gap: {max(per_row):.2f} m (median {np.median(per_row):.2f} m)")
print(f"step between neighbours inside the gap: p50 {np.median(interior):.2f} m, p99 {np.percentile(interior, 99):.2f} m")
print(f"step across the gap edge: p50 {np.median(steps):.2f} m, p99 {np.percentile(steps, 99):.2f} m; "
      f"between valid neighbours: p50 {np.median(normal):.2f} m, p99 {np.percentile(normal, 99):.2f} m")

profile.update(dtype="float32", nodata=ND, compress="deflate", predictor=3, tiled=False)
with rasterio.open(out_path, "w", **profile) as o:
    o.write(filled.astype(np.float32), 1)
