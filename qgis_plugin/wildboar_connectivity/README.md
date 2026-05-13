# Wildboar Connectivity v2 (QGIS plugin)

Interactive point-to-point connectivity across a resistance raster.

## What it does
Pick a resistance raster from the layer panel, click a **Start** and **End**
point on the map canvas, choose one or more methods, and hit **Run**:

| Method | Output | Library |
|---|---|---|
| Least-Cost Path (Dijkstra) | red corridor line (memory vector) | `scikit-image` |
| Circuit theory (current flow) | "Current flow / pinchpoints" raster (Magma ramp) | `scipy.sparse` |
| Graph-based (k-shortest) | k alternative corridor lines (cyan) | `networkx` |

All math runs in a `QgsTask` worker thread, so the QGIS UI stays responsive.

## Install (dev)
Symlink the plugin folder into the QGIS plugins directory:

```powershell
# from this repo root
cmd /c mklink /D `
  "$env:APPDATA\QGIS\QGIS3\profiles\default\python\plugins\wildboar_connectivity" `
  "$pwd\qgis_plugin\wildboar_connectivity"
```

Then enable **Wildboar Connectivity v2** in *Plugins -> Manage and Install*.

## Python dependencies
Beyond what ships with QGIS 3:

```text
scikit-image
networkx
scipy
rasterio
```

Install into the QGIS Python environment (OSGeo4W shell):

```cmd
python -m pip install scikit-image networkx scipy rasterio
```

## Files
* `__init__.py` - `classFactory` (QGIS entry point)
* `main_plugin.py` - plugin class + `ConnectivityTask` (background math)
* `point_tool.py` - `PointSelectionTool` (QgsMapToolEmitPoint subclass)
* `connectivity_dialog.py` / `connectivity_dialog_base.ui` - dialog
* `metadata.txt` - QGIS plugin metadata

## Math note: resistance to conductance
For the circuit step, each pair of 4-connected cells `(i, j)` becomes a
resistor with edge resistance `R_edge = 0.5 * (R_i + R_j)`. The Laplacian
uses conductance `g = 1 / R_edge`. Solving `L v = b` with `b = +1` at the
source and `-1` at the sink (sink grounded to keep `L` non-singular) yields
voltages `v`; node currents `0.5 * Sum_j |g_ij (v_i - v_j)|` highlight
pinchpoints. Full derivation lives in the docstring of `main_plugin.py`.
