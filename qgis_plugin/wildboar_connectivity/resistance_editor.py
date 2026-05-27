# -*- coding: utf-8 -*-
"""
Apply user-drawn fences and overpasses to a working copy of the
resistance raster.

A fence:    a polyline that gets rasterised as a ribbon of very high
            resistance, so the LCP and circuit solver treat it like a
            wall.
An overpass: a point that gets stamped as a small disc of very low
            resistance, so the solver can cross an otherwise impassable
            barrier (motorway, river) at exactly that spot.

These modifications happen on a numpy COPY of the resistance window so
the underlying raster file is never touched.
"""

import json

import numpy as np
from rasterio.features import rasterize


# Resistance values used for the modifications.
#
# A plugin-drawn fence is treated as a TRUE WALL: the cells in the
# fence ribbon become NaN, which propagates through Dijkstra (NaN
# becomes np.inf at the MCP step), Circuitscape, and the output
# masking, so the LCP and current-flow code cannot cross. NaN
# OVERWRITES any underlying resistance value, including pre-existing
# NaN walls from the notebook (highways, large lakes, ASP fences).
#
# A plugin-placed overpass is treated as the EASIEST POSSIBLE WALK:
# its disc gets a very low finite resistance, again OVERWRITING any
# underlying value. This lets the user punch a passable corridor
# through a highway / lake / fence by simply placing an overpass on
# it - exactly what a real wildlife crossing does in the landscape.
import numpy as _np

FENCE_RESISTANCE    = _np.nan
OVERPASS_RESISTANCE = 1e-3    # equals the lower clip floor used downstream


def apply_modifications(
        resistance: np.ndarray,
        transform,
        fences_geojson: list,
        overpasses_xy: list,
        fence_width_cells: int = 1,
        overpass_radius_cells: int = 2,
):
    """Return a new resistance grid with fences and overpasses applied.

    Parameters
    ----------
    resistance : 2D numpy array
        Working copy of the resistance window (will not be mutated).
    transform : affine.Affine
        The window's geo-transform; used by rasterio.features.rasterize.
    fences_geojson : list of dict
        GeoJSON-like geometry dicts (LineString) in the raster CRS.
    overpasses_xy : list of (x, y) tuples
        Overpass centres in the raster CRS.
    fence_width_cells : int
        Half-width of the fence ribbon in pixels. 1 means a one-pixel
        line (use rasterize all_touched=True); >1 means a thicker
        ribbon via morphological dilation.
    overpass_radius_cells : int
        Radius of the low-resistance disc, in pixels.
    """
    out = resistance.astype(np.float64, copy=True)
    rows, cols = out.shape

    # ---- Fences --------------------------------------------------
    if fences_geojson:
        # Burn every fence onto a single 0/1 mask, then dilate if asked.
        mask = rasterize(
            [(g, 1) for g in fences_geojson],
            out_shape=(rows, cols),
            transform=transform,
            fill=0,
            dtype="uint8",
            all_touched=True,
        ).astype(bool)
        if fence_width_cells > 1:
            mask = _dilate_disc(mask, fence_width_cells)
        out[mask] = FENCE_RESISTANCE

    # ---- Overpasses ----------------------------------------------
    if overpasses_xy:
        for (x, y) in overpasses_xy:
            r, c = _xy_to_rc(transform, x, y)
            if r is None:
                continue
            r1, r2 = max(0, r - overpass_radius_cells), \
                     min(rows, r + overpass_radius_cells + 1)
            c1, c2 = max(0, c - overpass_radius_cells), \
                     min(cols, c + overpass_radius_cells + 1)
            # Disc mask
            rr_idx, cc_idx = np.ogrid[r1:r2, c1:c2]
            disc = ((rr_idx - r) ** 2 + (cc_idx - c) ** 2
                    <= overpass_radius_cells ** 2)
            sub = out[r1:r2, c1:c2]
            sub[disc] = OVERPASS_RESISTANCE
            out[r1:r2, c1:c2] = sub

    # Re-clip to the same range used in _load_window so the solver stays
    # numerically stable.
    np.clip(out, 1e-3, 1e6, out=out)
    return out


def _xy_to_rc(transform, x, y):
    """Map a world coordinate to (row, col) using an inverse affine."""
    try:
        inv = ~transform
    except Exception:
        return None, None
    c, r = inv * (x, y)
    return int(round(r)), int(round(c))


def _dilate_disc(mask: np.ndarray, radius: int) -> np.ndarray:
    """Cheap morphological dilation with a disc structuring element.

    Used when fences are drawn as 1-pixel lines but should be displayed
    or computed as thicker ribbons. Pure numpy, no scipy.ndimage.
    """
    out = mask.copy()
    rows, cols = mask.shape
    for dr in range(-radius, radius + 1):
        for dc in range(-radius, radius + 1):
            if dr * dr + dc * dc > radius * radius:
                continue
            # Shift the mask by (dr, dc) and OR into out.
            r0, r1 = max(0, dr), rows + min(0, dr)
            c0, c1 = max(0, dc), cols + min(0, dc)
            src_r0, src_r1 = max(0, -dr), rows + min(0, -dr)
            src_c0, src_c1 = max(0, -dc), cols + min(0, -dc)
            out[r0:r1, c0:c1] |= mask[src_r0:src_r1, src_c0:src_c1]
    return out
