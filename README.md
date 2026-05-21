# SWOTCalVal_Yukon24

Calibration/validation (CalVal) workflow for SWOT water surface elevation (WSE), slope, and river width against field data collected during the Yukon River field campaign (July–September 2024). The study area spans the Yukon Flats domain, covering the Yukon River (YR), Porcupine River (PR), Chandalar River (CD), Coleen River (CL), and Sheenjek River (SJ).

SWOT products evaluated: **RiverSP** (versions PIC0/SWORD v16 and PGD0/SWORD v17b) (and **RiverTile**).

---

## Workflow Overview

Scripts are numbered by phase. Phase 0 prepares input data; Phases 1–3 compute in situ vs. SWOT residuals by in situ measurement type; Phase 4 aggregates results and produces figures and tables for the manuscript.

---

### Phase 0 — Data Preparation

| Script | Description |
|--------|-------------|
| `0.0_SWOT_OpenNC_SaveCSV.R` | Batch-reads SWOT L2 PIXC and PIXCVec NetCDF files, extracts key variables (including derived WSE and uncertainty fields), optionally applies quality filters, and writes one CSV per file. |
| `0.1_RiverTile_timeseries_creation.ipynb` | Builds node- and reach-level SWOT RiverTile time series CSVs from individual overpass files for use in downstream comparisons. |
| `0.2_ortho_resampling.ipynb` | Downsamples RGB/NIR orthomosaic rasters from 1 m to 3 m resolution using rasterio (average resampling). |
| `0.3_PT_GNSS_reach_slope_bias_corr.R` | Applies per-sensor (PT) and per-drift (GNSS) WSE bias corrections derived from node-level comparisons to reach boundary WSEs, then recomputes bias-corrected reach slopes. |

---

### Phase 1 — Pressure Transducer (PT) vs. SWOT

| Script | Description |
|--------|-------------|
| `1.1_PT_consistency_checks.R` | Verifies PT WSE records are physically consistent (monotonically decreasing downstream, no unexplained jumps between adjacent sensors). |
| `1.2_PT_node_wse_diff.R` | Time- (±7.5 min) and space-matches PT WSE to SWOT node WSE; computes absolute residuals, removes per-PT median bias, and exports per-cluster CSVs and a merged CalVal CSV. |
| `1.3_PT_reach_wse_slope_diff.R` | Matches PT reach WSE and bias-corrected slope to SWOT reach observations; computes absolute and relative residuals and exports to CSV. |

---

### Phase 2 — GNSS Drift Survey vs. SWOT

| Script | Description |
|--------|-------------|
| `2.1_GNSS_node_wse_diff.R` | Time- (±5 hr) and space-matches GNSS drift survey WSE to SWOT node WSE; computes absolute residuals, removes per-drift median bias, and exports to CSV. |
| `2.2_GNSS_reach_wse_slope_diff.R` | Matches GNSS reach WSE and slope to SWOT reach observations; computes absolute and bias-corrected residuals for both variables and exports to CSV. |

---

### Phase 3 — Orthomosaic Width vs. SWOT

| Script | Description |
|--------|-------------|
| `3.1_PIXCVec_concatenate_node_polygons.R` | Reads PIXCVec tiles, groups by overpass timestamp, concatenates same-overpass tiles into sf point datasets, and generates Voronoi node polygons for the AOI. |
| `3.2_ortho_PIXCVec_node_polygon_widths.ipynb` | Computes water area within each SWORD node polygon from PIXCVec and orthomosaic data; outputs per-node water area CSVs. |
| `3.3_ortho_node_width_diff.R` | Joins orthomosaic water area to SWOT node observations on the overpass date; derives ortho width from area ÷ node length; computes absolute residuals and percent differences, removes per-reach-code median bias, and exports per-cluster and merged CSVs. |
| `3.4_ortho_reach_width_diff.R` | Aggregates node-level ortho widths to reach means, matches to SWOT reach observations, computes residuals and percent differences, and produces summary tables by river. |

---

### Phase 4 — Aggregated Validation & Figures

| Script | Description |
|--------|-------------|
| `4.1_wse_comparisons.R` | Consolidates PT and GNSS node- and reach-level WSE residuals across SWOT versions; produces manuscript Tables 3–4, S1–S2 and Figures 4, 6. |
| `4.2_slope_comparisons.R` | Consolidates PT and GNSS reach-level slope residuals (in cm/km) across SWOT versions; produces manuscript Tables 5–6, S3–S4 and Figure 5. |
| `4.3_width_comparisons.R` | Consolidates orthomosaic node-level width residuals across SWOT versions (outliers |residual| ≥ 1500 m excluded); produces manuscript Tables 7–8, S5 and Figure 7. |
| `4.4_wse_slope_domain_inclusion.R` | Determines which nodes/reaches were observed by Version C only, Version D only, both, or neither; writes shapefiles and produces domain summary tables and in situ summary figures (Figures 2, 3; Table 1). |
| `4.5_river_kernel_densityoverviews.R` | Generates per-river kernel density distributions of SWORD prior width and slope, and SWOT functional repeat time statistics at node and reach scales. |

---

## Data Dependencies

- **SWOT products**: PIXC, PIXCVec, RiverSP (PIC0/v16, PGD0/v17b), RiverTile (v16, v17b)
- **SWORD**: SWORD v16 and v17b (North America, HB81)
- **In situ**: Pressure transducer (PT) WSE time series (processed with CalVal toolboxes); GNSS boat-drift WSE surveys (processed with CalVal toolboxes & JPL GipsyX software); RGB/NIR orthomosaics water masks (manually digitized)

## Software Requirements

- **R** (tidyverse, sf, lubridate, dplyr, ncdf4, deldir, ggtext)
- **Python** (rasterio, numpy, and Jupyter for `.ipynb` scripts)
