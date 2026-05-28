# -*- coding: utf-8 -*-
"""
iSSF-based Individual-Based Movement Model (IBMM) for ASF spread in wild boar.

Each agent represents one infected wild boar starting at the outbreak point.
The agent moves through the landscape using a step-selection rule derived from
the fitted iSSF: at each step, IBMM_N_CANDIDATES tentative steps are generated
(Gamma step-length distribution, von Mises turning-angle distribution) and the
agent moves to the candidate with probability proportional to 1/resistance at
its endpoint. The agent stops when it reaches the end of its stochastic
infectious period (time from infection to death, sampled per agent from an
ASF-calibrated Gamma distribution).

The output raster is the contamination-probability surface: the fraction of
agents that visit each cell, i.e. the probability that a single outbreak at the
selected point leads to an infected animal reaching that cell before dying.

iSSF coefficients
-----------------
From issf_coefficients.csv (scaled covariates, all positive):
    forest_density_sc      0.5614   boars prefer denser forest
    dist_railway_sc        0.4244   prefer distance from railway
    dist_settlement_sc     0.2633   prefer distance from settlements
    dist_main_road_sc      0.2624   prefer distance from main roads
    forest_density_g_sc    0.2057   prefer forest-density gradients
    dist_secondary_road_sc 0.1967   prefer distance from secondary roads
    dist_water_sc          0.0778   prefer distance from water

Movement kernel
---------------
Fitted to WildBoar_InMatrix_Steps.csv (in-matrix / dispersal state, 15-min GPS
interval): Gamma(shape=2.5, scale=160 m) → mean ≈ 400 m per active step;
von Mises(kappa=0.5) → turning-angle SD ≈ 85° (moderate persistence).

ASF infectious period
---------------------
Time from infection to death: Gamma(shape=3, scale=3.5 d) → mean ≈ 10.5 d,
SD ≈ 6 d. Consistent with EFSA (2018), Probst et al. (2017), and Chenais et
al. (2019). During the infectious period the animal moves normally (no
behavioural change before clinical signs). Daily active movement steps derived
from GPS telemetry: ~20 in-matrix steps per day (wild boar mean daily path
≈ 4–8 km; at 400 m/step → 10–20 active dispersal steps/day).

References
----------
Avgar T et al. (2016). Integrated step selection analysis. MEE 7: 619-630.
EFSA (2018). Scientific Opinion on ASF. EFSA J 16(7): 5344.
Probst C et al. (2017). Combining hunting mortality and wild boar movement.
    Transbound Emerg Dis 64: 1820-1833.
Chenais E et al. (2019). Epidemiological considerations on ASF in Europe.
    Porcine Health Manag 5: 10.
"""

import numpy as np

# --------------------------------------------------------------------------
# iSSF coefficients — stored for documentation/reference. The simulation
# uses 1/resistance as the selection surface, which captures the same
# information because the resistance raster was built by inverting the
# iSSF suitability surface.
# --------------------------------------------------------------------------
ISSF_COEFFICIENTS = {
    "forest_density_sc":      0.5614,
    "dist_railway_sc":        0.4244,
    "dist_settlement_sc":     0.2633,
    "dist_main_road_sc":      0.2624,
    "forest_density_g_sc":    0.2057,
    "dist_secondary_road_sc": 0.1967,
    "dist_water_sc":          0.0778,
}

# --------------------------------------------------------------------------
# Movement kernel — Gamma step-length distribution fitted to
# WildBoar_InMatrix_Steps.csv (15-min GPS interval, in-matrix state).
# --------------------------------------------------------------------------
IBMM_SL_GAMMA_SHAPE   = 2.5    # dimensionless
IBMM_SL_GAMMA_SCALE_M = 160.0  # metres  →  mean = 2.5 × 160 = 400 m

# Von Mises concentration for turning angles.
# kappa = 0.5 → SD ≈ 85°, moderate directional persistence.
IBMM_KAPPA = 0.5

# Candidate steps evaluated at each move ("available steps" in iSSF parlance).
IBMM_N_CANDIDATES = 25

# --------------------------------------------------------------------------
# ASF infectious period — Gamma distribution for time from infection to death.
# Mean ≈ 10.5 days, SD ≈ 6 days (EFSA 2018; Probst et al. 2017).
# --------------------------------------------------------------------------
ASF_IP_GAMMA_SHAPE      = 3.0   # dimensionless
ASF_IP_GAMMA_SCALE_DAYS = 3.5   # days  →  mean = 3 × 3.5 = 10.5 d

# Active in-matrix movement steps per day (GPS telemetry).
# Wild boar mean daily path ≈ 4–8 km; at mean step length 400 m →
# 10–20 dispersal steps/day. Use 20 as the upper bound for dispersing
# animals (which are the epidemiologically relevant ones).
ASF_ACTIVE_STEPS_PER_DAY = 20

# Hard cap: absolute maximum steps regardless of the sampled lifespan.
# Prevents runaway agents in edge cases (e.g. very long lifespan samples).
IBMM_MAX_STEPS_CAP = 3000


def simulate_ibmm(
    hab_suitability: np.ndarray,
    nodata_mask,
    origin_cells: np.ndarray,
    win_tf,
    aoi_mask,
    n_agents: int = 2000,
    rng=None,
    cancel_check=None,
) -> np.ndarray:
    """Simulate n_agents infected wild boar and return a contamination-
    probability raster.

    Each agent:
      1. Starts at a random cell within the outbreak zone.
      2. Lives for a stochastic number of active movement steps sampled from
         Gamma(ASF_IP_GAMMA_SHAPE, ASF_IP_GAMMA_SCALE_DAYS)
         × ASF_ACTIVE_STEPS_PER_DAY (= infectious period in steps).
      3. At each step, selects its next cell from IBMM_N_CANDIDATES
         tentative endpoints proportionally to hab_suitability there
         (the iSSF selection rule).
      4. Records a visit at every cell it passes through.

    The output is contamination probability: visits[cell] / n_agents, i.e.
    the probability that at least one agent visits the cell (for large
    n_agents this converges to the true visitation probability under the
    model). Values range from 0 (no agent ever reached this cell) to 1
    (every agent passed through — only the origin cells under most
    parameterisations).

    Parameters
    ----------
    hab_suitability : np.ndarray (rows, cols), float64
        Habitat suitability / selection surface (≥ 0). Cells with value 0
        are impassable; in practice pass 1/resistance from the caller.
    nodata_mask : np.ndarray (rows, cols bool) or None
        True = impassable cell (wall, fence, NoData).
    origin_cells : np.ndarray (int64, flat indices)
        Pool of cells from which agents are started.
    win_tf : rasterio.transform.Affine
        Affine transform — used only to read cell size.
    aoi_mask : np.ndarray (rows, cols bool) or None
        Agents that step outside the AOI are stopped immediately.
    n_agents : int
        Number of infected-boar trajectories to simulate.
    rng : np.random.Generator or None
    cancel_check : callable() → bool or None
        Called every 100 agents; returns True if the QGIS task was
        cancelled.

    Returns
    -------
    np.ndarray (float32)
        Contamination-probability grid. Zero cells are left as 0;
        the caller should replace 0 with NaN before writing to QGIS.
    """
    if rng is None:
        rng = np.random.default_rng()

    rows, cols = hab_suitability.shape
    cellsize_m = float(abs(win_tf.a))
    sl_scale_px = IBMM_SL_GAMMA_SCALE_M / cellsize_m

    # Prepare selection surface.
    hs = hab_suitability.astype(np.float64, copy=True)
    hs[~np.isfinite(hs)] = 0.0
    hs[hs < 0.0] = 0.0
    if nodata_mask is not None:
        hs[nodata_mask] = 0.0

    # visit[r, c] = number of agents that stepped on cell (r, c).
    visit = np.zeros((rows, cols), dtype=np.float64)

    orr, occ = np.unravel_index(origin_cells, (rows, cols))
    n_starts = len(orr)

    for agent_idx in range(n_agents):
        if cancel_check is not None and agent_idx % 100 == 0:
            if cancel_check():
                break

        # ---- Sample this agent's infectious lifespan -----------------
        life_days = rng.gamma(ASF_IP_GAMMA_SHAPE, ASF_IP_GAMMA_SCALE_DAYS)
        agent_steps = min(
            int(round(life_days * ASF_ACTIVE_STEPS_PER_DAY)),
            IBMM_MAX_STEPS_CAP,
        )
        if agent_steps < 1:
            agent_steps = 1

        # ---- Initialise position and heading -------------------------
        k0 = int(rng.integers(0, n_starts))
        r, c = int(orr[k0]), int(occ[k0])
        heading = float(rng.uniform(0.0, 2.0 * np.pi))

        # ---- Walk for agent_steps steps ------------------------------
        for _ in range(agent_steps):
            if aoi_mask is not None and not aoi_mask[r, c]:
                break

            visit[r, c] += 1.0

            # Generate IBMM_N_CANDIDATES tentative steps.
            sls = rng.gamma(IBMM_SL_GAMMA_SHAPE, sl_scale_px,
                            size=IBMM_N_CANDIDATES)
            tas = rng.vonmises(0.0, IBMM_KAPPA, size=IBMM_N_CANDIDATES)
            dirs = heading + tas

            # Pixel coordinates of candidate endpoints.
            # Column increases East; row increases South.
            # Heading 0 = East, π/2 = North (row decreases).
            nc = c + sls * np.cos(dirs)
            nr = r - sls * np.sin(dirs)
            nc_i = np.round(nc).astype(np.intp)
            nr_i = np.round(nr).astype(np.intp)

            valid = ((nr_i >= 0) & (nr_i < rows) &
                     (nc_i >= 0) & (nc_i < cols))
            if not valid.any():
                break

            nr_safe = np.where(valid, nr_i, 0).astype(np.intp)
            nc_safe = np.where(valid, nc_i, 0).astype(np.intp)
            w = np.where(valid, hs[nr_safe, nc_safe], 0.0)

            wsum = float(w.sum())
            if wsum <= 0.0:
                break

            w /= wsum
            k = int(rng.choice(IBMM_N_CANDIDATES, p=w))

            if not valid[k]:
                break

            heading = float(np.arctan2(
                -(float(nr_i[k]) - r),
                float(nc_i[k]) - c,
            ))
            r, c = int(nr_i[k]), int(nc_i[k])

    # ---- Contamination probability -----------------------------------
    # Divide by n_agents (not by total visits) so each cell's value is
    # the fraction of agents that visited it = P(contamination | outbreak
    # at this origin). Interpretable on a 0–1 scale without further
    # normalisation.
    prob = visit / max(n_agents, 1)
    return prob.astype(np.float32)
