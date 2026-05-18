# Wildboar Connectivity v2 (QGIS plugin)

Local-scale habitat connectivity for wild boar: pick **one point on the
map**, the plugin defines a circular movement neighbourhood around it,
finds every core habitat in that circle, and computes pairwise LCPs, a
cumulative Circuitscape pinchpoint map, and a network graph with
centrality scores.

## Workflow

1. Load `Potenzial_Wildschwein_Kerneinstaende_ZH_2024_WILMA.shp` and your
   resistance raster into QGIS.
2. Open the plugin (toolbar icon or *Plugins -> &Wildboar Connectivity v2*).
3. Pick the resistance raster and the habitats layer from the dropdowns.
4. Click **Pick center point on map**, then click anywhere — this point
   is the centre of the analysis circle (no need to hit a habitat).
5. Set the **wild boar range radius** (default 7 km — slightly larger
   than a typical home range, smaller than dispersal).
6. Enable the methods you want, press **Run**.

Outputs (added as layers, top to bottom):

| Layer | Type | Content |
|---|---|---|
| `Pinchpoints (log10 cum. current, N pairs)` | raster (Magma ramp) | sum of \|edge currents\| over every habitat pair, log10 stretched, habitats masked |
| `LCP corridors (N pairs)` | line vector | one feature per pair, attrs `from_id`, `to_id`, `cost` |
| `Habitat graph - edges (E)` | line vector | straight edges between habitat centroids, attrs `cost`, `weight = 1/cost` |
| `Habitat graph - nodes (N)` | point vector | habitat centroids with `degree`, `betweenness`, `area_km2` |
| `Selected habitats (N)` | polygon vector | the habitats inside the circle |
| `Analysis range (buffer)` | polygon vector | the wild-boar range circle |
| `Analysis centre` | point vector | the click point |

## Why a small scale?

Wild boar daily movement is ~1-5 km; mean dispersal ~10 km. Modeling
connectivity across an entire canton mixes ecologically incompatible
scales. Restricting each analysis to a local neighborhood around one
core habitat is both faster and biologically meaningful — and it
mirrors how an ASF response zone would be defined.

## What the pinchpoint map actually shows

Each pixel value = sum of |edge current| flowing through that cell,
summed across every habitat pair. Bright streaks = corridors used by
**multiple** pairs. The brightest narrow streaks are **pinchpoints**:
cells where many corridors converge, whose loss would disconnect the
local network. These are the natural targets for fencing decisions or
landscape restoration.

Values are log10(amps) because raw current spans 3+ orders of magnitude.

## Mathematical core

For 4-connected cells (i, j):

    R_edge = 0.5 * (R_i + R_j)         # cells in series, half each
    g_ij   = 1 / R_edge                 # conductance

Weighted graph Laplacian:

    L = D - A,   D_ii = sum_j g_ij,   A_ij = g_ij

For each pair of habitats (A, B):

    b[A] = +1/|A|,  b[B] = -1/|B|        # area sources (no point singularities)
    solve  L v = b   with one corner grounded
    I_e = g_ij (v_i - v_j)                # edge currents

Cumulative across pairs:

    Cum_e = sum_pairs |I_e|
    node  = 0.5 * sum_{j ~ i} Cum(i,j)

The Laplacian is factored **once** (`scipy.sparse.linalg.splu`); every
pairwise solve is a cheap back-substitution.

## Files

* `main_plugin.py` - plugin class, UI wiring, habitat picking, buffering
* `connectivity_task.py` - `QgsTask` with all the math (rasterise habitats, build L, factor, solve, accumulate)
* `point_tool.py` - `QgsMapToolEmitPoint` subclass that emits a signal carrying a role tag
* `connectivity_dialog.py` / `connectivity_dialog_base.ui` - dialog
* `__init__.py`, `metadata.txt`, `icon.png`

## Install (dev)

Copy to the QGIS plugins folder:

```powershell
Copy-Item -Recurse -Force `
  "$pwd\qgis_plugin\wildboar_connectivity" `
  "$env:APPDATA\QGIS\QGIS3\profiles\default\python\plugins\"
```

Or, for live editing, create a directory symlink (Admin PowerShell):

```powershell
Remove-Item "$env:APPDATA\QGIS\QGIS3\profiles\default\python\plugins\wildboar_connectivity" -Recurse -Force
New-Item -ItemType SymbolicLink `
  -Path  "$env:APPDATA\QGIS\QGIS3\profiles\default\python\plugins\wildboar_connectivity" `
  -Target "$pwd\qgis_plugin\wildboar_connectivity"
```

Then enable **Wildboar Connectivity v2** in *Plugins -> Manage and
Install*.

## Python dependencies

Beyond what ships with QGIS 3 (rasterio is usually bundled):

```cmd
python -m pip install scikit-image scipy
```

Run that in the OSGeo4W Shell so the packages land in QGIS's Python.

## Known limits / future work

- **4-connected grid.** Mild Manhattan bias in the current field. An
  8-connected variant is a 10-line change.
- **One starting habitat per run.** True all-to-all across the entire
  canton is intentionally not done — see "Why a small scale" above.
- **Single Laplacian factorization** assumes the resistance window fits
  in memory comfortably. Up to ~250 k cells (e.g. 500x500) is fine on
  a laptop; beyond that, switch to CG.
