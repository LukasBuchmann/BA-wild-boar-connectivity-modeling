# Wildboar Connectivity (ASF) — QGIS Plugin

Decision-support tool for African Swine Fever (ASF) spread modelling. No
habitat polygons required — just a resistance raster and a single click.

## Workflow

1. Load your Keeley-transformed iSSF resistance raster into QGIS.
2. Open the plugin (toolbar icon or *Plugins → Wildboar Connectivity (ASF)*).
3. Select the resistance raster and (optionally) a core-habitat polygon layer.
4. Click **Pick outbreak point** and click the infection location on the map.
5. Enable the desired outputs and press **Run**.
6. Optionally draw fences or place wildlife overpasses, then press **Run** again
   to compare.

## Outputs

| Layer | Content |
|---|---|
| Infection risk (resistant kernel) | `exp(−cost/D)` decay from origin, green/yellow/red |
| LCP corridors | Least-cost paths to nearest low-resistance clusters, weighted by cluster size |
| Pinchpoints (Circuit theory) | Circuitscape.jl current density; bright = movement bottleneck |
| BCRW random-walk density | Biased correlated random-walk visit counts (optional) |
| iSSF-IBMM contamination probability | log₁₀(p) contamination per cell from agent simulation (optional) |

## Analysis method

**Source disc** — a small disc (radius = 2 pixels, ≈ 50 m) at the clicked
pixel. Snaps to the nearest passable cell if the click lands on a wall.

**Dijkstra** — one sweep from the source disc gives cost-distances to all
reachable cells. The AOI is every cell with finite cost.

**Infection risk** — `risk = exp(−cost / D_cost)` where `D_cost` is the cost
of walking 4 km (published mean wild-boar natal dispersal) through
median-resistance terrain. Maps directly to EFSA (2018) operational zones:
risk > 0.5 → protection zone, 0.15–0.5 → surveillance zone.

**LCP corridors** — for each low-resistance cluster in the AOI, the cheapest
MCP traceback from its entry cell back to the source disc is found. Segments
shared by many LCPs accumulate strength; the backbone corridor stands out.

**Circuit theory** — Circuitscape.jl (Julia) if available, otherwise a scipy
sparse Laplacian solve. Source = small disc at origin; sink = AOI boundary
cells. Source cells are masked to NaN before display to remove the injection
artefact.

**BCRW** — biased correlated random walk. Each step selects a neighbour
proportional to `1/(R × step_dist) × exp(κ cos θ)` (cost bias × von Mises
directional persistence, κ = 2).

**iSSF-IBMM** — per-agent stochastic infectious period from
`Gamma(shape=3, scale=3.5 d)` × 20 steps/day (EFSA 2018). Selection surface
recovered from the resistance raster via the inverse Keeley transform
`S = 1 − log(R)/4`. Output is `log₁₀(contamination probability)`.

## Landscape modifications

Draw a **fence** (left-click vertices, right-click to finish) to burn NaN
walls into a working copy of the resistance grid. Place an **overpass** to
punch a minimum-resistance disc through any barrier. Press **Reset** to clear
all modifications.

## Installation (development)

```powershell
# Symlink (Admin PowerShell — live-reload on save)
New-Item -ItemType SymbolicLink `
  -Path   "$env:APPDATA\QGIS\QGIS3\profiles\default\python\plugins\wildboar_connectivity" `
  -Target "$pwd\qgis_plugin\wildboar_connectivity"
```

Then enable **Wildboar Connectivity (ASF)** in *Plugins → Manage and Install*.

## Python dependencies

```cmd
# Run in the OSGeo4W Shell
python -m pip install scikit-image scipy
```

Julia + Circuitscape.jl (optional, for the pinchpoint solver):
```julia
using Pkg; Pkg.add("Circuitscape")
```

## File structure

| File | Role |
|---|---|
| `main_plugin.py` | Plugin entry, UI wiring, tool activation |
| `connectivity_task.py` | QgsTask: all analysis (Dijkstra, LCPs, Circuit, BCRW, IBMM) |
| `connectivity_dialog_base.ui` | Qt Designer dialog layout |
| `connectivity_dialog.py` | uic loader for the dialog |
| `circuitscape_jl.py` | Subprocess wrapper for Circuitscape.jl |
| `resistance_editor.py` | Applies fences/overpasses to the resistance array |
| `fence_tool.py` | Polyline-drawing map tool |
| `point_tool.py` | Single-point capture map tool |
| `layer_utils.py` | Layer z-ordering and project management helpers |
| `__init__.py`, `metadata.txt`, `icon.png` | QGIS plugin boilerplate |
