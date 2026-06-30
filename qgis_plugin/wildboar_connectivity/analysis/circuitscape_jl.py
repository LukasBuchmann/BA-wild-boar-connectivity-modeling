# -*- coding: utf-8 -*-
"""
Subprocess wrapper around the official Circuitscape.jl (Julia).

Used as the primary pinchpoint solver. Falls back to the scipy.sparse
implementation in connectivity_task._single_source_current_scipy()
if Julia or Circuitscape.jl is not available, or if the subprocess
fails.

Setup the user needs (one-time):
    1. Install Julia from https://julialang.org (>=1.6)
    2. In a Julia REPL:  using Pkg; Pkg.add("Circuitscape")

Invoke Julia via subprocess and pass it an INI file. Inputs and
outputs are ESRI ASCII grids written in a tempdir that is cleaned up
afterwards.
"""

import os
import shutil
import subprocess
import tempfile

import numpy as np


class CircuitscapeJlError(RuntimeError):
    pass


# ---------------------------------------------------------------------
# Availability checks
# ---------------------------------------------------------------------
def is_julia_available() -> bool:
    """True if a `julia` executable is reachable via PATH."""
    return shutil.which("julia") is not None


_CIRCUITSCAPE_INSTALLED_CACHE = None


def is_circuitscape_installed(timeout: int = 90) -> bool:
    """True if `using Circuitscape` succeeds in Julia. Cached."""
    global _CIRCUITSCAPE_INSTALLED_CACHE
    if _CIRCUITSCAPE_INSTALLED_CACHE is not None:
        return _CIRCUITSCAPE_INSTALLED_CACHE
    if not is_julia_available():
        _CIRCUITSCAPE_INSTALLED_CACHE = False
        return False
    try:
        res = subprocess.run(
            ["julia", "--startup-file=no", "-e",
             "using Circuitscape"],
            capture_output=True, text=True, timeout=timeout,
        )
        ok = res.returncode == 0
    except Exception:
        ok = False
    _CIRCUITSCAPE_INSTALLED_CACHE = ok
    return ok


# ---------------------------------------------------------------------
# ESRI ASCII grid I/O (Circuitscape's preferred format)
# ---------------------------------------------------------------------
def _write_ascii_grid(arr: np.ndarray, transform, path: str,
                      nodata: float = -9999.0) -> None:
    """Write a 2D float array to ESRI .asc.

    Assumes a north-up, square-pixel affine transform - which is what
    rasterio.window_transform always produces for us.
    """
    rows, cols = arr.shape
    cellsize = float(transform.a)
    if abs(abs(transform.e) - cellsize) > 1e-6:
        raise CircuitscapeJlError(
            "Non-square pixels are not supported by ESRI ASCII grid.")
    xll = float(transform.c)
    yll = float(transform.f) + float(transform.e) * rows

    out = np.where(np.isfinite(arr), arr, nodata).astype(np.float64)
    with open(path, "w") as f:
        f.write(f"ncols         {cols}\n")
        f.write(f"nrows         {rows}\n")
        f.write(f"xllcorner     {xll}\n")
        f.write(f"yllcorner     {yll}\n")
        f.write(f"cellsize      {cellsize}\n")
        f.write(f"NODATA_value  {nodata}\n")
        np.savetxt(f, out, fmt="%.6g")


def _read_ascii_grid(path: str):
    """Read an ESRI .asc into a numpy array, mapping NoData -> NaN."""
    header = {}
    with open(path) as f:
        for _ in range(6):
            line = f.readline()
            if not line:
                break
            key, val = line.split()
            header[key.lower()] = val
        data = np.loadtxt(f, dtype=np.float64)
    nodata = float(header.get("nodata_value", -9999.0))
    data[data == nodata] = np.nan
    return data, header


# ---------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------
def run_one_to_all(resistance: np.ndarray,
                   transform,
                   source_mask: np.ndarray,
                   sink_mask: np.ndarray,
                   timeout: int = 900,
                   log: callable = None) -> np.ndarray:
    """Run Circuitscape.jl in pairwise (regions) mode with 2 polygons.

    Parameters
    ----------
    resistance : 2D array (float), shape (rows, cols)
        Resistance grid for the analysis window. NaN cells are treated
        as NoData by Circuitscape (region in landscape is excluded).
    transform : affine.Affine
        Geo-transform of the resistance window.
    source_mask, sink_mask : 2D bool arrays
        True where the origin habitat / the union of all other habitats
        sits. Each is treated as one Circuitscape "region polygon" so
        the solve becomes a single source-to-sink problem.
    timeout : seconds
    log : callable(str) | None
        Optional logger.

    Returns
    -------
    np.ndarray
        The cumulative current density grid (rows, cols), same shape
        and transform as the input resistance, NaN outside the AOI.
    """
    if not is_julia_available():
        raise CircuitscapeJlError("Julia not found on PATH.")
    if not is_circuitscape_installed():
        raise CircuitscapeJlError(
            "Circuitscape.jl not installed. Run in Julia: "
            "using Pkg; Pkg.add(\"Circuitscape\")")

    tmpdir = tempfile.mkdtemp(prefix="cs_jl_")
    try:
        rast_path  = os.path.join(tmpdir, "resistance.asc")
        focal_path = os.path.join(tmpdir, "focal_regions.asc")
        ini_path   = os.path.join(tmpdir, "config.ini")
        out_prefix = os.path.join(tmpdir, "out")

        _write_ascii_grid(resistance, transform, rast_path)

        # Focal regions: 1 = source (origin habitat), 2 = sink (all
        # other habitats unioned). Cells outside are NoData.
        focal = np.full(resistance.shape, -9999.0, dtype=np.float64)
        focal[source_mask] = 1.0
        focal[sink_mask & ~source_mask] = 2.0
        _write_ascii_grid(focal, transform, focal_path)

        # Circuitscape is configured entirely through this INI file rather
        # than function arguments, since it is invoked as an external Julia
        # process. "pairwise" mode with two focal regions (source/sink, see
        # the focal raster above) reduces the general all-pairs solve to the
        # single source-to-sink current map this plugin needs.
        ini = f"""[Circuitscape Mode]
data_type = raster
scenario = pairwise

[Habitat raster or graph]
habitat_file = {rast_path}
habitat_map_is_resistances = True

[Connection scheme for raster habitat data]
connect_four_neighbors_only = False
connect_using_avg_resistances = True

[Short circuit regions (aka polygons)]
use_polygons = False

[Options for pairwise and one-to-all and all-to-one modes]
point_file = {focal_path}
point_file_contains_polygons = True

[Mask file]
use_mask = False

[Calculation options]
solver = cg+amg
preemptive_memory_release = False
print_timings_in_logs = False
parallelize = False

[Output options]
output_file = {out_prefix}
write_cur_maps = True
write_cum_cur_map_only = True
log_transform_maps = False
write_volt_maps = False
"""
        with open(ini_path, "w") as f:
            f.write(ini)

        if log:
            log(f"Circuitscape.jl: launching Julia ({resistance.shape[0]}x"
                f"{resistance.shape[1]} window)...")
        script = (
            "using Circuitscape; "
            f'compute(raw"{ini_path}")'
        )
        result = subprocess.run(
            ["julia", "--startup-file=no", "-e", script],
            capture_output=True, text=True, timeout=timeout,
        )
        if result.returncode != 0:
            raise CircuitscapeJlError(
                f"Circuitscape.jl failed: {result.stderr.strip() or result.stdout[-500:]}")

        # Circuitscape names its output file differently depending on the
        # version and on whether it ran in "cumulative" or "pairwise" mode,
        # so try the two known patterns first ...
        candidates = [
            f"{out_prefix}_cum_curmap.asc",
            f"{out_prefix}_curmap_1_2.asc",
        ]
        # ... and fall back to scanning the tempdir for anything that looks
        # like a current map, so an unexpected naming scheme in a future
        # Circuitscape release does not silently break the pipeline.
        for f in os.listdir(tmpdir):
            if f.endswith(".asc") and "curmap" in f.lower():
                candidates.append(os.path.join(tmpdir, f))
        out_path = next((p for p in candidates if os.path.exists(p)), None)
        if out_path is None:
            files = ", ".join(os.listdir(tmpdir))
            raise CircuitscapeJlError(
                f"Circuitscape produced no current map. Files: {files}")

        if log:
            log(f"Circuitscape.jl: reading {os.path.basename(out_path)}")
        current, _ = _read_ascii_grid(out_path)
        return current

    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)
