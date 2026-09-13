# Validation of SWOT water surface elevation, slope, and inundation extent in the Yukon River Basin, Alaska, USA

This repository contains the code for the article:

> Kluetmeier, C.L., Pavelsky, T.M., Rowley, T., Bennitt F.B., Gleason, C.J.,
> Minear, J.T., Smith, L.C., Thellman, A., Tooley, W., Olson, A.J., Wald, M.,
> Cooley S.W. Validation of SWOT water surface elevation, slope, and
> inundation extent in the Yukon River Basin, Alaska, USA. (Submitted).

Calibration/validation (CalVal) workflow for SWOT water surface elevation
(WSE), slope, and river width against field data collected in the Yukon
River Basin (July–September 2024). The domain covers the Yukon River (YR), 
Porcupine River (PR), Chandalar River (CD), Coleen River (CL), Sheenjek River 
(SJ), and Black River (BL).

SWOT products evaluated: **RiverSP** version C0 (PIC0, SWORD v16) and version
D0 (PGD0, SWORD v17b).

---

## Repository structure

```
SWOTCalVal_Yukon24/
├── 0_fetch_SWOT_data/
│   ├── 0.0_hydrocron_SWOT_fetch.ipynb
│   └── 0.1_PT_GNSS_reach_slope_bias_corr.R
├── 1_PT_SWOT_harmonization/
│   ├── 1.1_PT_consistency_checks.R
│   ├── 1.2_PT_node_wse_diff.R
│   └── 1.3_PT_reach_wse_slope_diff.R
├── 2_GNSS_SWOT_harmonization/
│   ├── 2.1_GNSS_node_wse_diff.R
│   └── 2.2_GNSS_reach_wse_slope_diff.R
├── 3_Orthomosaic_SWOT_harmonization/
│   ├── 3.1_RiverObs/
│   ├── 3.2_ortho_PIXCVec_node_polygon_widths.ipynb
│   ├── 3.3_ortho_node_width_diff.R
│   └── 3.4_ortho_reach_width_diff.R
├── 4_SWOT_validation/
│   ├── 4.0_comparison_helpers.R
│   ├── 4.1_wse_comparisons.R
│   ├── 4.2_slope_comparisons.R
│   ├── 4.3_width_comparisons.R
│   ├── 4.4_wse_slope_domain_inclusion.R
│   ├── 4.5_river_kernel_densityoverviews.R
│   └── 4.6_domain_flags.R
├── .gitignore
└── README.md
```

Scripts are numbered by phase. Phase 0 prepares input data; Phases 1–3
compute in situ vs. SWOT residuals by in situ measurement type; Phase 4
aggregates results and produces the figures and tables in the manuscript.

---

## Data availability

SWOT L2_HR_RiverSP versions C0 and D0 were accessed via the Hydrocron API
(Greguska et al., 2024). The Level 2 Hydrology data are also available through
NASA EarthData (https://search.earthdata.nasa.gov/search) and CNES
HydrowebNext (https://hydroweb.next.theia-land.fr). The SWOT River Database v16 
and v17b are available to download from the SWORD website 
(https://www.swordexplorer.com/; Altenau et al., 2021). Compiled in situ
observations, software, and corresponding metadata for this analysis are
available in the following data release (Kluetmeier et al., 2026). Software
for RiverObs and upstream processing of PT and GNSS observations to SWORD
products are available in the following GitHub repositories
(https://github.com/SWOTAlgorithms/RiverObs,
https://github.com/cjgleason/calval_toolbox).

---

## Data dependencies

- **SWOT products**: PIXCVec, RiverSP (PIC0/v16, PGD0/v17b)
- **SWORD**: SWORD v16 and v17b (North America, HB81), both the shapefile and netCDF distributions
- **SWORD ID translators**: `NA_NodeIDs_v17b_vs_v16.csv`, `NA_ReachIDs_v17b_vs_v16.csv`
- **In situ**: pressure transducer (PT) WSE time series and GNSS boat-drift WSE surveys (both processed with the CalVal toolboxes); RGB/NIR orthomosaic water masks (manually digitized)

## Software requirements

- **R**: tidyverse, lubridate, sf, scales
- **Python**: numpy, pandas, geopandas, rasterio, shapely, pyproj, netCDF4,
  requests, matplotlib, folium, Jupyter
- **RiverObs**: installed separately from its own repository; see
  `3.1_RiverObs/README.md`
