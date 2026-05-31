# -*- coding: utf-8 -*-
"""
Apply user-drawn fences and overpasses to a working copy of the resistance
raster (the source GeoTIFF is never modified).

Fence   — polyline rasterised as a ribbon of NaN (impassable wall).
Overpass — point stamped as a small disc of minimum resistance (1e-3),
           overriding any barrier including highways and large lakes.
"""

import numpy as np
from rasterio.features import rasterize

# NaN turns fence cells into true walls recognised by MCP, Circuitscape,
# and the AOI masking step. 1e-3 matches the lower clip floor used in
# _load_window so an overpass is always the easiest passable cell.
FENCE_RESISTANCE    = np.nan
OVERPASS_RESISTANCE = 1e-3


def apply_modifications(
        resistance: np.ndarray,
        transform,
        fences_geojson: list,
        overpasses_xy: list,
        fence_width_cells: int = 1,
        overpass_radius_cells: int = 2,
) -> np.ndarray:
    """Return a new resistance grid with fences and overpasses applied.

    Parameters
    ----------
    resistance : 2D float64 array
        Working copy of the resistance window (not mutated).
    transform : affine.Affine
        Geo-transform of the window.
    fences_geojson : list of GeoJSON-like dicts (LineString, raster CRS)
    overpasses_xy : list of (x, y) tuples (raster CRS)
    fence_width_cells : int
        Half-width of the fence ribbon in pixels.
    overpass_radius_cells : int
        Radius of the low-resistance disc, in pixels.
    """
    out = resistance.astype(np.float64, copy=True)
    rows, cols = out.shape

    # ---- Fences ---------------------------------------------------------
    if fences_geojson:
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

    # ---- Overpasses -----------------------------------------------------
    if overpasses_xy:
        for (x, y) in overpasses_xy:
            r, c = _xy_to_rc(transform, x, y)
            if r is None:
                continue
            r1 = max(0, r - overpass_radius_cells)
            r2 = min(rows, r + overpass_radius_cells + 1)
            c1 = max(0, c - overpass_radius_cells)
            c2 = min(cols, c + overpass_radius_cells + 1)
            rr_idx, cc_idx = np.ogrid[r1:r2, c1:c2]
            disc = (rr_idx - r) ** 2 + (cc_idx - c) ** 2 \
                   <= overpass_radius_cells ** 2
            out[r1:r2, c1:c2][disc] = OVERPASS_RESISTANCE

    # Clip only finite cells so fence NaN walls are preserved exactly.
    finite = np.isfinite(out)
    out[finite] = np.clip(out[finite], 1e-3, 1e6)
    return out


def _xy_to_rc(transform, x, y):
    """Map a world coordinate to (row, col) via the inverse affine."""
    try:
        inv = ~transform
    except Exception:
        return None, None
    c, r = inv * (x, y)
    return int(round(r)), int(round(c))


def _dilate_disc(mask: np.ndarray, radius: int) -> np.ndarray:
    """Morphological dilation with a disc structuring element (pure numpy)."""
    out = mask.copy()
    rows, cols = mask.shape
    for dr in range(-radius, radius + 1):
        for dc in range(-radius, radius + 1):
            if dr * dr + dc * dc > radius * radius:
                continue
            r0, r1 = max(0, dr), rows + min(0, dr)
            c0, c1 = max(0, dc), cols + min(0, dc)
            src_r0 = max(0, -dr);  src_r1 = rows + min(0, -dr)
            src_c0 = max(0, -dc);  src_c1 = cols + min(0, -dc)
            out[r0:r1, c0:c1] |= mask[src_r0:src_r1, src_c0:src_c1]
    return out
