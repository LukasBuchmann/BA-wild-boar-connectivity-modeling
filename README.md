# Wild Boar Connectivity and ASF Containment

This repository contains the code and report for a Bachelor's thesis in Applied Digital Life Science (ZHAW). The project builds an interactive QGIS planning tool that helps wildlife administrations and authorities prepare for African Swine Fever (ASF) outbreaks among wild boar (*Sus scrofa*) in the Swiss canton of Zurich.

The work has two parts that build on each other: an R analysis pipeline that turns GPS telemetry into a wild boar resistance surface and a QGIS plugin that uses these results to run interactive, scenario-based ASF spread simulations directly on a map.

## Repository layout

```
.
├── data/
│   ├── raw/         Source data: telemetry, SwissTLM3D, DHM25, boundaries, fences, passages
│   ├── processed/   Pipeline outputs: resistance and suitability rasters, model results, figures
│   └── temp/        Scratch files created while the notebooks run
├── docs/            Quarto book project containing the thesis report (chapters as .qmd files)
├── qgis_plugin/
│   └── wildboar_connectivity/   The "Wildboar Connectivity (ASF)" QGIS plugin
└── scripts/
    ├── 01_create_resistance_surface.ipynb   Telemetry to resistance surface (R / iSSF)
    └── 02_connectivity.ipynb                unused testing file
```

All paths inside the notebooks and the plugin are written relative to their own location, so the project can be moved or cloned anywhere on disk without editing any code.

## What the pipeline does

**`scripts/01_create_resistance_surface.ipynb`** loads cleaned GPS telemetry for Canton Aargau, computes movement metrics, and classifies each step as either "in-patch" (resting in habitat) or "in-matrix" (travelling between patches) with a rule-based classifier. It then fits an integrated step-selection function (iSSF, conditional logistic regression via `survival::clogit`) on the in-matrix steps to learn which landscape covariates wild boar prefer or avoid while moving. The notebook keeps both an exploratory 13-covariate model and the final, parsimonious 7-covariate model side by side, with comments explaining how the variable selection was done (VIF checks and plausibility considerations). The resulting selection coefficients are converted into a 25 m resistance surface for Aargau via the inverse Keeley exponential transformation, and the same statistics are then used to transfer the model and produce a resistance surface for Canton Zurich.

**`qgis_plugin/wildboar_connectivity/`** is the interactive planning tool itself. Loaded into QGIS together with a Keeley-transformed resistance raster, it lets a user click a single point as an outbreak origin and computes, on demand: an infection-risk raster, least-cost path corridors to the nearest low-resistance habitat clusters, a Circuitscape/circuit-theory pinchpoint raster, an optional biased correlated random walk (BCRW) density surface, and an optional iSSF-based individual-based movement model (iSSF-IBMM) that estimates contamination probability over time. Fences and wildlife overpasses can be drawn directly on the map, and re-running the analysis immediately shows how they change the predicted spread.

**`docs/`** is a Quarto book project containing the full thesis report, built from the chapter files listed in `docs/_quarto.yml`.

## Setting up the project on your own machine

### 1. Get the code and data

Clone or copy this repository. The `data/raw/` folder must contain the source datasets referenced by the notebooks (telemetry CSV, SwissTLM3D and swissBOUNDARIES3D geopackages, the DHM25 elevation raster, the wildlife passages geodatabase, and the fence shapefiles). These are large geodata files and are not tracked in git; obtain them from the data sources cited in the report and place them under `data/raw/` using the same file names that appear in the notebook configuration cells, so the relative paths resolve correctly.

### 2. R environment (for the analysis notebooks)

Both notebooks run on an R kernel (`ir`, i.e. the IRkernel for Jupyter) and install their own dependencies on first run via `pacman::p_load(...)`. You need:

- R (version 4.x) and Jupyter with the R kernel (`IRkernel::installspec()`)
- The `pacman` package, which the notebooks install automatically if missing

The notebooks then pull in the packages they need themselves: `amt`, `terra`, `sf`, `dplyr`, `tidyr`, `purrr`, `lubridate`, `survival`, `broom`, `ggplot2` for the resistance surface pipeline, and additionally `igraph`, `gdistance`, `car` for the connectivity analysis and VIF checks.

Open the notebooks in Jupyter (or in RStudio via `jupytext`/the Jupyter extension) and run the cells in order, starting with `01_create_resistance_surface.ipynb`. Its outputs (resistance and suitability rasters, model coefficients, diagnostic plots) are written to `data/processed/`.

### 3. Julia and Circuitscape (for circuit-theory analysis)

Both the connectivity notebook and the QGIS plugin's pinchpoint analysis use Circuitscape.jl as their primary solver, invoked as a Julia subprocess. To enable it:

1. Install Julia (>= 1.6) from <https://julialang.org> and make sure the `julia` executable is on your PATH.
2. In a Julia REPL, run `using Pkg; Pkg.add("Circuitscape")`.

If Julia or Circuitscape.jl is not available, the QGIS plugin automatically falls back to a pure Python solver based on `scipy.sparse` (see `qgis_plugin/wildboar_connectivity/analysis/connectivity_task.py`), so the plugin remains usable without a Julia installation, just slower on large rasters.

### 4. QGIS and the plugin's Python dependencies

The plugin targets QGIS 3.16 or newer (see `qgis_plugin/wildboar_connectivity/metadata.txt`) and uses QGIS's bundled Python environment, so no separate virtual environment is needed for the plugin itself. It does rely on a few third-party Python packages that are not always bundled with QGIS by default: `numpy`, `rasterio`, `scipy`, and `scikit-image`. If QGIS reports an import error when loading the plugin, install the missing packages into QGIS's Python environment, for example through the OSGeo4W shell on Windows:

```
python -m pip install rasterio scipy scikit-image
```

### 5. Installing the plugin in QGIS

1. Locate your QGIS profile's plugin folder. On Windows this is typically `%APPDATA%\QGIS\QGIS3\profiles\default\python\plugins\`; on Linux/macOS it is `~/.local/share/QGIS/QGIS3/profiles/default/python/plugins/`.
2. Copy (or symlink) the `qgis_plugin/wildboar_connectivity/` folder into that plugin directory, keeping the folder name `wildboar_connectivity`.
3. Restart QGIS, then enable the plugin via *Plugins → Manage and Install Plugins → Installed* and ticking "Wildboar Connectivity (ASF)".
4. Load a Keeley-transformed resistance raster (produced by `01_create_resistance_surface.ipynb` and found in `data/processed/`) into your QGIS project, open the plugin, click an outbreak origin on the map, and press Run.
