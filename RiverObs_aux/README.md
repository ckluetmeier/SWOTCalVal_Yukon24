# riverobs-watermask-width

**Deriving SWOT-comparable river width from digitized water masks, using the SWOT RiverObs processor.**

This repository provides a reproducible pipeline that takes a manually digitized
water mask — from aerial orthomosaics, drone imagery, or any high-resolution
optical source — and processes it through
[RiverObs](https://github.com/SWOTAlgorithms/RiverObs), the same software that
generates the operational SWOT River Single-Pass (RiverSP) vector product.

The result is a reference width dataset on the SWORD node and reach grid that
differs from the SWOT product **only in the input pixel cloud** — not in how
water is attributed to nodes, not in how area is aggregated, and not in how area
is converted to width.

- **Pipeline version:** 1.1.0
- **Tested against RiverObs commit:** `dc0c7fd` (2026-07-29)
- **Tested prior databases:** SWORD v16, v17b (North America)

---

## Table of contents

1. [Why this exists](#1-why-this-exists)
2. [Method](#2-method)
3. [Repository contents](#3-repository-contents)
4. [Installation](#4-installation)
5. [Patching RiverObs](#5-patching-riverobs)
6. [Input requirements](#6-input-requirements)
7. [Quick start](#7-quick-start)
8. [Step reference](#8-step-reference)
9. [Configuration reference](#9-configuration-reference)
10. [Output data dictionary](#10-output-data-dictionary)
11. [Verification](#11-verification)
12. [Design decisions](#12-design-decisions)
13. [Known limitations](#13-known-limitations)
14. [Adapting to a new study area](#14-adapting-to-a-new-study-area)
15. [Citation and license](#15-citation-and-license)

---

## 1. Why this exists

Validating SWOT river width requires a reference measurement defined the same
way SWOT's is. The obvious approach — draw polygons around each node, count
water pixels inside — introduces a subtle but consequential problem if those
polygons are derived from the SWOT product itself (for example, Thiessen
polygons built from PIXCVec point locations):

1. **The reference becomes conditional on the observation.** PIXCVec points
   exist only where SWOT detected water. A tessellation of those points is
   clipped to the envelope of SWOT detections, so reference water outside that
   envelope is never counted. The design then captures SWOT's *commission*
   errors at full strength while truncating its *omission* errors.

2. **The reference varies with viewing geometry.** PIXCVec point density depends
   on cross-track distance, so tessellation cell size does too. Any analysis of
   width residuals against cross-track distance has the reference itself as a
   function of the regressor.

3. **No cross-channel cutoff.** Dissolved tessellation cells extend laterally as
   far as the point envelope allows, whereas RiverObs applies an explicit
   cross-channel limit plus a segmentation-based extension. Reference and
   observation search different footprints.

Running the reference mask through RiverObs removes all three. The node grid,
the centerline, the node length, the cross-channel search corridor, and the
area-to-width conversion all come from the prior database and from the
processor, not from the observation being validated.

---

## 2. Method

### 2.1 The width definition is unchanged

RiverObs computes node width in `RiverObs/RiverNode.py`:

```python
def width_area(self, area_var='pixel_area'):
    area = self.sum(area_var)
    width_area = area / self.ds     # self.ds is the SWORD prior node_length
    return width_area
```

That is `width = summed water area / prior node length` — the same definition
most manual approaches already use. What this pipeline changes is not the
formula but **which water is summed for which node**.

Reach width is `sum(node area) / sum(node p_length)` over observed nodes: a
length-weighted mean, matching the operational reach product, not an unweighted
mean of node widths.

### 2.2 How RiverObs assigns water to nodes

RiverObs does not use polygons. Each water cell is assigned to the node nearest
in along-track distance, in curvilinear (along-reach, cross-reach) coordinates
on a spline through the prior-database nodes, subject to a cross-channel cutoff
(`max_width / 3`) and a segmentation-based extension governed by
`ext_dist_coef`. The region a node's area is summed over is a curvilinear cell
whose lateral edges follow the water mask — not a rectangle, and not expressible
as a simple polygon. [Step 3b](#step-3b--node-footprints) recovers the actual
assigned regions from the processor's own output.

### 2.3 The search corridor

RiverObs only sums water inside a cross-channel corridor sized from the **prior
database's** channel width. Two gates apply, both scaled by the per-node
`wth_coef` and `ext_dist_coef`:

| Gate | Where | Half-corridor |
|---|---|---|
| 1 | `assign_reaches` → `RiverObs.get_ext_dist_threshold` | `(2 × prior_width × wth_coef) / 3` |
| 2 | `assign_reaches_ext_dist_coef` | `ext_dist_coef × max(node_spacing, max(prior_max_width, prior_width) × wth_coef)` |

The tighter of the two binds. For a SWOT pixel cloud this is a sensible filter:
it keeps distant lakes and unrelated water out. For a **digitized reference mask
in which every polygon is known-good river water, it is a loss, not a filter** —
anabranches and braidplain threads that sit outside the prior channel are
excluded, and the node width comes out too narrow.

Two environment variables widen both gates together:

```bash
export RIVEROBS_WTH_COEF_FACTOR=3.0
export RIVEROBS_EXT_DIST_COEF_FACTOR=3.0
```

They multiply the prior coefficients, which are corridor knobs only and are not
written to any output product. `1.0` is stock RiverObs. The shim logs the
resulting half-corridor in metres per reach, so the effect is explicit:

```
shim: corridor reach 0 -- prior width 120 m, half-corridor gate1 480 m,
      gate2 21600 m (binding: 480 m)
```

**Tune, don't guess.** `util_assignment_audit.py` reports the fraction of each
mask that reached a node. Raise the factors until that plateaus at ~100%, then
stop — an over-wide corridor lets adjacent reaches claim the same water, which
the audit reports as double counting. Worked example on a synthetic braidplain
(120 m main channel plus anabranches at ±320 m and ±380 m):

| Factor | Assigned | Median node width | Double counting |
|---|---|---|---|
| 1.0 (stock) | 41.38% | 119.9 m | none |
| 2.0 | 62.27% | — | none |
| **3.0** | **100.00%** | **289.7 m** | none |
| 4.0 | 100.00% | — | none |

The true mask area over the 8 km reach corresponds to a mean width of 290.0 m,
so the relaxed run recovers it to 0.1%. The stock run recovers only the main
channel — which is exactly the visual signature of an over-tight corridor: node
polygons that trace the main thread and stop at the edge of the digitized water.

### 2.4 Pipeline

```
   digitized water mask (.shp)
              │
              ▼
   ┌─────────────────────┐
   │ step1               │  burn every feature to a regular grid
   │ rasterize           │  → <survey>_water_<res>m.tif
   └─────────────────────┘
              │
              ▼
   ┌─────────────────────┐
   │ step2               │  SimplePixelCloud → CalValToRiverTile
   │ RiverObs            │  → RiverTile .nc, PIXCVec .nc, combined CSVs
   └─────────────────────┘        ▲
              │                   └── SWORD prior database (.nc)
              ▼
   ┌─────────────────────┐
   │ step3b              │  dissolve water cells by assigned node_id
   │ node footprints     │  → true_node_polygons_<survey>_<version>.shp
   └─────────────────────┘
              │
              ▼
      ╔═══════════════════╗
      ║  MANUAL: QGIS     ║  digitize cloud / ragged survey edges into ONE
      ╚═══════════════════╝  exclusion layer covering every survey
              │
              ▼
   ┌─────────────────────┐
   │ step3c --batch      │  one call, every survey x every version
   │ exclusion list      │  → node_qc_all.csv
   └─────────────────────┘     (survey, sword_version, node_id, keep, excl_frac)
```

Everything before the manual step is scripted and deterministic. No node is
excluded automatically — see [§12.2](#122-why-there-is-no-automated-coverage-gate).

---

## 3. Repository contents

| File | Purpose |
|---|---|
| `apply_riverobs_patches.py` | Applies the four source edits to a RiverObs checkout. Idempotent, with `--status` and `--revert`. |
| `riverobs_shim.py` | Five runtime patches applied by import (four guards plus the opt-in corridor relaxation). Touches no source files. |
| `run_calval2rivertile.py` | Wrapper around RiverObs's `calval2rivertile.py` that loads the shim first. |
| `step0_check_sword_nc.py` | Verifies a prior-database netCDF is usable, and that width-critical variables are present. |
| `step1_rasterize_watermask.py` | Burns a water-mask shapefile to a regular grid. |
| `step2_run_riverobs.py` | Runs every survey × prior-database version; flattens outputs to CSV. |
| `step3b_true_node_polygons.py` | Recovers the actual per-node water footprints from the PIXCVec assignment. |
| `step3c_manual_node_qc.py` | Converts a QGIS review into a node exclusion list. `--batch` applies one exclusion layer to every survey and version in a single call. |
| `util_assignment_audit.py` | Reports what fraction of each mask reached a node, exports the unassigned water, and detects cross-reach double counting. |
| `util_watermask_feature_audit.py` | Audits water-mask shapefiles for the multi-feature undercount. |
| `run_all.sh` | Batch driver for steps 1, 2, 3b and 3c. Sources `config.sh` if present. |
| `config.example.sh` | Site paths and the survey manifest. Copy to `config.sh`, which is gitignored. |
| `riverobs_ortho.template.rdf` | Annotated RiverObs configuration template. Copy to `riverobs_ortho_<version>.rdf` per prior-database version; the copies are gitignored. |
| `environment.yml` | Verified minimal conda environment. |

Every script prints its pipeline version as the first line of output.

---

## 4. Installation

### 4.1 Environment

Use `environment.yml` in this repository rather than RiverObs's own. The latter
requires `pysal`, GDAL/`osgeo`, `Rtree` and `scikit-image`, none of which are
imported anywhere on the code path this pipeline uses, and pinning them makes
the solve slow or unsatisfiable on some platforms.

```bash
conda env create -f environment.yml
conda activate RiverObs
```

Two dependencies are worth calling out because they are imported at the top of
`SWOTRiverEstimator.py` but are absent from RiverObs's `Install.md` package
list: **`bottleneck`** and **`piecewise-regression`**.

### 4.2 RiverObs

RiverObs is used from source, not installed:

```bash
git clone https://github.com/SWOTAlgorithms/RiverObs.git
export PYTHONPATH="/path/to/RiverObs/src:$PYTHONPATH"
```

Add the `export` to your shell profile to make it persistent. Edits under
`src/` take effect on the next run with no rebuild.

Verify:

```bash
python -c "import SWOTRiver.Estimate; from RiverObs.ReachDatabase import ReachExtractor; print('ok')"
```

A `SyntaxWarning` from `calval.py` is expected and harmless — it is present
upstream.

### 4.3 Prior database

RiverObs reads the SWORD **netCDF** distribution, not the shapefiles. Point
`reach_db_path` at a continent-level `.nc` file or at a directory of them.

---

## 5. Patching RiverObs

Two source files, four edits, plus four runtime patches. All are guards and
plumbing: **none touch node assignment, area aggregation, or the area-to-width
conversion.**

```bash
python apply_riverobs_patches.py --riverobs-root /path/to/RiverObs
```

### 5.1 Source edits (`apply_riverobs_patches.py`)

| File | Edit | Rationale |
|---|---|---|
| `src/SWOTRiver/products/calval.py` | Implement the `from_airborne_imagery` stub | Upstream reserves this classmethod for exactly this input type but leaves it `NotImplementedError`. |
| `src/SWOTRiver/products/calval.py` | Declare `geolocation_qual`, `classification_qual`, `sig0_qual` on `SimplePixelCloud` | Two are commented out upstream and one is misspelled `classificaiton_qual`; `sig0_qual` is absent. Without all three, `process_node()` raises `TypeError: 'NoneType' object is not subscriptable`. |
| `src/bin/calval2rivertile.py` | Register the `airborne_watermask` format | Makes the converter reachable from the CLI. |
| `src/bin/calval2rivertile.py` | Add the dispatch branch | As above. |

The converter itself is adapted from `from_cnes_watermask()` with two
corrections: the CRS is read from the file rather than hard-coded to EPSG:32610,
and `rasterio.transform.xy()` output is reshaped before boolean masking (it
returns flat sequences for 2-D row/column input).

### 5.2 Runtime patches (`riverobs_shim.py`)

Applied on import; no files are modified. These live in modules whose source
text varies between RiverObs revisions, so they key on class and function names
instead.

| Patch | Symptom without it |
|---|---|
| `SWOTRiverEstimator.__init__` — default `reach_pct_good_sus_thresh` to 0 | `TypeError: '<=' not supported between instances of 'float' and 'NoneType'`. `CalValToRiverTile` never passes this argument although `L2PixcToRiverTile` does. |
| `RiverReach.__init__` — broadcast scalar attributes to node count | `ValueError: zero-dimensional arrays cannot be concatenated`. Under `height_agg_method='orig'`, `wse_s_u` is a scalar constant. |
| `ReachExtractor.__init__` — pad `centerline_lon`/`centerline_lat` to equal length | `ValueError: setting an array element with a sequence ... inhomogeneous shape`. Prior-database reaches have different centerline point counts; since numpy 1.24 a ragged list raises. The product declares these as 2-D `(reaches, centerlines)`, so a rectangular array is required anyway, and `rivertile.py` already filters the padding via `np.abs(lats) < 90`. |
| `RiverObs.get_node_stat` — substitute for two absent variables | `AttributeError: 'RiverNode' object has no attribute 'time_from_prev_xover'`. Per-line crossover timing does not exist in a water mask. Only these two are substituted; every other missing variable re-raises so RiverObs's own handlers run unchanged. |

Because of the shim, **run jobs through `run_calval2rivertile.py`**, not
`src/bin/calval2rivertile.py` directly. `step2_run_riverobs.py` does this
automatically.

---

## 6. Input requirements

### 6.1 Water mask

- Polygon shapefile (or any OGR-readable vector format)
- Projected CRS with metre units
- Water only — the pipeline treats every polygon as water
- Multipart or multi-feature both fine; **every** feature is rasterized

Run `util_watermask_feature_audit.py` once over your masks to confirm they
contain what you expect.

### 6.2 Orthomosaic (optional but recommended)

Used as the QGIS backdrop for the manual review, and optionally by step 1 for a
coverage cross-check. Not a processing input.

### 6.3 Prior database

SWORD netCDF for the region of interest, one file per version you intend to
compare against.

### 6.4 Survey manifest

One line per survey in `config.sh`, plus a matching entry in the `SURVEYS` list
at the top of `step2_run_riverobs.py` carrying the observation date.

---

## 7. Quick start

```bash
conda activate RiverObs
export PYTHONPATH="/path/to/RiverObs/src:$PYTHONPATH"

# one-time
python apply_riverobs_patches.py --riverobs-root /path/to/RiverObs
python step0_check_sword_nc.py /path/to/na_sword_v17b.nc --verify
cp config.example.sh config.sh && $EDITOR config.sh
cp riverobs_ortho.template.rdf riverobs_ortho_v17b.rdf && $EDITOR riverobs_ortho_v17b.rdf

# batch: steps 1, 2 and 3b produce the node polygons; nothing is excluded yet
bash run_all.sh 1 2 3b
bash run_all.sh 1               # a single step

# manual QGIS review: digitize cloud / ragged edges into ONE exclusion layer
# covering every survey, at the EXCLUDE_SHP path set in config.sh

# step 3c applies it to every survey and version at once
bash run_all.sh 3c              # re-run this alone after any edit to the layer

# or the whole thing, once the exclusion layer exists
bash run_all.sh
```

---

## 8. Step reference

### Step 0 — prior-database check

```bash
python step0_check_sword_nc.py /path/na_sword_v17b.nc --verify
```

Read-only. Reports how many of RiverObs's declared variables the file carries,
how many of the width-critical subset are present, and whether the required
`/reaches` subgroups exist. `--verify` loads the database through RiverObs and
calls one reach, which exercises the full read path.

**Missing variables are normal and harmless.** RiverObs declares far more than
the public releases carry; `Product.__getattr__` returns a fully-masked array
for anything absent. `/reaches/type` is missing from public SWORD and RiverObs
derives it from `reach_id % 10`, which is exact.

The check matters precisely *because* absence is silent rather than fatal: a
database missing `node_length` would not crash, it would produce plausible and
wrong widths. `node_length` is the divisor in `width = area / node_length`;
`max_width` and `ext_dist_coef` size the search corridor.

**Exit codes:** 0 usable, 1 not usable.

---

### Step 1 — rasterize the water mask

```bash
python step1_rasterize_watermask.py \
    /path/water_masks/<survey>_watermask.shp \
    /path/output/rasterized_water_masks \
    --name <survey> --res 3 \
    [--ortho-tif /path/orthomosaics/<survey>_ortho.tif]
```

Burns every polygon feature to a `uint8` grid, snapped so cell area is exactly
`res²`. `all_touched=False` — a cell must have its centre inside the polygon —
which keeps the area estimate unbiased.

**Output:** `<name>_water_<res>m.tif` (the only file step 2 needs).

**Grid resolution.** Node area is a sum of cell areas, so spacing sets only the
quantization of that sum. At 3 m a 200 m node in a 50 m channel holds ~3,300
cells and quantization error is a fraction of a percent — far below digitization
and centerline error. Going finer multiplies pixel-cloud size quadratically and
forces `SWOTRiverEstimator`'s connected-component segmentation to allocate an
image the size of the full survey grid. **3–5 m is the recommended range.**
Measured accuracy at 3 m against a channel of known area: **0.025%**.

`--pad` (default 2000 m) must exceed the largest prior `max_width` in the area,
or the search corridor runs off the edge of the grid.

`--ortho-tif` additionally writes a footprint raster and reports water cells
falling outside the imaged area — a useful one-off check that digitizing stayed
within the survey. Nothing downstream reads it.

---

### Step 2 — run RiverObs

```bash
python step2_run_riverobs.py \
    --riverobs-root /path/to/RiverObs \
    --raster-dir  /path/output/rasterized_water_masks \
    --rdf-dir     . \
    --out-dir     /path/output/riverobs \
    --res 3 --versions v17b v16 --keep-going
```

Loops over the `SURVEYS` list × `--versions`, invoking `run_calval2rivertile.py`
for each pair, then flattens the RiverTile node and reach groups into two
combined CSVs.

| Option | Effect |
|---|---|
| `--survey NAME [NAME ...]` | Restrict to named surveys |
| `--keep-going` | Continue past a failed survey, still write CSVs for the rest, list failures and exit non-zero |
| `--log-level` | `debug` / `info` / `warning` / `error` |

**Outputs**

```
<out-dir>/<version>/ortho_<survey>_<version>.nc            RiverTile
<out-dir>/<version>/ortho_<survey>_<version>_pixcvec.nc    per-pixel node assignment
<out-dir>/<version>/ortho_<survey>_<version>_pixc.nc       reformatted pixel cloud
<out-dir>/ortho_riverobs_nodes_all.csv
<out-dir>/ortho_riverobs_reaches_all.csv
```

**Expected log noise.** These are normal and can be ignored: `PIXC Quality flags
not found in input file`; `converting a masked element to nan`; `divide by zero
encountered in divide` at `ww = 1/(wse_r_u**2)`; `Mean of empty slice`;
`assign_reaches: no observations mapped to nodes in this reach` for reaches
outside the mask; and the shim's `time_from_prev_xover` / `time_to_next_xover`
substitution notice. All originate in the deliberately disabled WSE machinery
(see [§9](#9-configuration-reference)).

---

### Step 3b — node footprints

```bash
python step3b_true_node_polygons.py \
    /path/rasterized_water_masks/<survey>_water_3m.tif \
    /path/riverobs/<version>/ortho_<survey>_<version>_pixcvec.nc \
    /path/node_qc/true_node_polygons_<survey>_<version>.shp \
    --rivertile /path/riverobs/<version>/ortho_<survey>_<version>.nc \
    [--reaches]
```

The PIXCVec output carries a `node_id` for every water cell. Dissolving by it
recovers the exact region each node's area and width were summed over. Verified
against the RiverTile: polygon areas match `area_total` to **0.000 m²**, and
`polygon area / p_length` reproduces `width` exactly.

`--rivertile` joins `width`, `p_length`, `p_dist_out` and `n_good_pix` so the
layer can be styled by width. The script also prints nodes narrower than 60% of
the median as a **hint** for the review — nothing is excluded.

**Attributes written:** `node_id`, `node_str` (string copy, authoritative if the
numeric field is ever mangled by a DBF), `n_cells`, `area_m2`, `keep`
(initialized to 1), plus the RiverTile join if requested.

---

### Manual review

Load into QGIS:

1. the orthomosaic
2. the water mask, for context
3. `true_node_polygons_<survey>_<version>.shp`

Style the node layer graduated on `width`. Nodes truncated by the survey edge or
with part of the channel under cloud read as anomalously narrow against their
neighbours.

---

### Step 3c — exclusion list

**Route A — exclusion polygon (recommended).** Digitize a polygon layer covering
the areas you do not trust. **One layer holding polygons for every survey is
fine** — the test is spatial, so a polygon only affects nodes it overlaps.

Batch, every survey and every version in one call (this is what `run_all.sh 3c`
runs):

```bash
python step3c_manual_node_qc.py --batch \
    --nodes-dir   node_qc \
    --exclude-shp digitizing/bad_ends_watermasks.shp \
    --versions v17b v16 \
    --out-dir node_qc \
    --out node_qc/node_qc_all.csv
```

It discovers every `true_node_polygons_<survey>_<version>.shp` under
`--nodes-dir`, writes one CSV per survey/version plus the merged `--out`, and
prints a kept/dropped table. A survey that fails is reported and skipped rather
than taking the rest of the run down.

Single survey/version:

```bash
python step3c_manual_node_qc.py \
    --nodes  node_qc/true_node_polygons_<survey>_<version>.shp \
    --exclude-shp digitizing/bad_ends_watermasks.shp \
    --survey <survey> --sword-version <version> \
    --out node_qc/node_qc_<survey>_<version>.csv
```

Nodes overlapping the exclusion area by more than `--max-overlap` (default 5%)
are dropped; near-misses are listed. Because the exclusion is spatial, one act
of judgement filters **every** prior-database version consistently — different
versions place nodes differently, but the cloud is in the same place on the
ground. It is also the reproducible artefact: a reader sees which ground area
was excluded rather than an unexplained list of identifiers.

> **If two surveys overlap on the ground.** The whole layer is applied to every
> survey. That is right when the polygons mark ground that is bad in all of
> them, and wrong when two surveys cover the same reach on different dates with
> different bad areas — a polygon drawn for one date would also clip the other.
> In that case add a text field naming the survey each polygon belongs to and
> pass `--exclude-field <name>`. Polygons whose field is empty or `all` still
> apply everywhere. The script prints the layer's field names on startup.

**Route B — `keep` field.** Set `keep = 0` in QGIS and omit `--exclude-shp`.
Deleting rows works too, with `--original` pointing at an untouched copy.
`--exclude-ids` accepts a comma-separated list. This route is per-layer by
nature and has no batch form.

**Merge existing CSVs:**

```bash
python step3c_manual_node_qc.py --merge 'node_qc/node_qc_*.csv' \
    --out node_qc/node_qc_all.csv
```

Step 3c reads nothing but the step 3b shapefiles and the exclusion layer, so
after any edit to that layer in QGIS, `bash run_all.sh 3c` alone refreshes
`node_qc_all.csv` in seconds — nothing upstream is repeated.

---

## 9. Configuration reference

The RiverObs `.rdf` passed to `calval2rivertile.py`. Values not listed take
RiverObs defaults. Keep every parameter except `reach_db_path` identical across
prior-database versions, or a cross-version comparison measures your settings
rather than the products.

| Key | Value | Why |
|---|---|---|
| `reach_db_path` | path to SWORD `.nc` | The only line that changes between versions. |
| `class_list` | `[4]` | Must equal `class_value` in `from_airborne_imagery`. |
| `use_fractional_inundation` | `[False]` | A digitized mask is binary. |
| `use_segmentation` | `[True]` | **Do not set False.** `SWOTRiverEstimator` always passes a `max_width`, and `flag_out_channel_and_label()` calls `.copy()` on the labels unconditionally, so `False` raises `AttributeError`. Because the converter sets `range_index`/`azimuth_index` to raster column/row, segmentation is a proper connected-component labelling on the mask grid. |
| `area_agg_method` | `orig` | Node area = `sum(pixel_area)`, width = area / `p_length`. The `composite` method is built around SWOT's water fraction, sig0 and edge classes, which a digitized mask does not have. |
| `height_agg_method` | `orig` | **Must also be `orig`.** RiverObs skips `get_node_agg()` only when height *and* area are both `orig`; any other value enters the interferometric height aggregator, which needs `phase_noise_std` / `dh_dphi`. |
| `ds` | `None` | Keeps the prior node spacing, so output `node_id` matches the SWOT product one-to-one. |
| `minobs` | `10` | Permissive by design — report thin coverage and filter it explicitly in the review. |
| `min_fit_points` | `1000000` | Set above any real node count so every WSE/slope fit is skipped. Heights are constant zero, so `wse_r_u` is 0 and the WLS weights `1/wse_r_u²` are infinite (`LinAlgError: SVD did not converge`). No fit touches area or width — reach width uses `mask_area`, not `mask_wse`. |
| `outlier_method` | `None` | Nothing for outlier rejection to work on. |
| `scalar_max_width` | `600.0` | Fallback only; RiverObs scales the corridor per node from prior `max_width` and `ext_dist_coef`. |
| `reach_pct_good_sus_thresh` | `0` | Every mask cell is good, so the tested percentage is always 100 and any lower value is neutral. |
| quality thresholds | `0` | A digitized mask carries no per-pixel quality information. |

> **`wse`, `wse_u` and `slope` in the output are meaningless.** Pixel-cloud
> heights are zero and all WSE fitting is disabled. Read only `width`,
> `area_total`, `area_detct`, `p_length`, `n_good_pix`, and the reach-level
> `obs_frac_n` / `partial_f`.

---

## 10. Output data dictionary

### `ortho_riverobs_nodes_all.csv`

| Column | Units | Description |
|---|---|---|
| `survey` | — | Survey name from the manifest |
| `SWOT_date` | ISO date | Matching observation date (UTC) |
| `sword_version` | — | Prior-database version tag |
| `reach_id`, `node_id` | — | Prior-database identifiers |
| `lat`, `lon` | deg | Node location |
| `ortho_width_m` | m | **Reference width** = assigned water area / `p_length` |
| `ortho_width_u_m` | m | Width uncertainty (`NaN` under `area_agg_method = orig`) |
| `ortho_area_total_m2` | m² | Summed area of assigned water cells |
| `ortho_area_detct_m2` | m² | Detected area; equals `area_total` for a binary mask |
| `n_good_pix` | count | Cells contributing |
| `p_length` | m | Prior node length — the width divisor |
| `p_width` | m | Prior width |
| `p_dist_out` | m | Distance to outlet |

### `ortho_riverobs_reaches_all.csv`

| Column | Units | Description |
|---|---|---|
| `reach_id` | — | Prior-database identifier |
| `ortho_width_m` | m | `sum(node area) / sum(node p_length)`, length-weighted |
| `ortho_area_total_m2` | m² | Reach water area |
| `obs_frac_n` | 0–1 | Fraction of the reach's prior nodes that received data |
| `partial_f` | 0/1 | 1 when `obs_frac_n < 0.5` |
| `n_good_nod` | count | Nodes contributing |

`obs_frac_n` is computed **before** the manual node exclusions. If reach-scale
completeness matters, recompute it from the `keep`-filtered node set.

### `node_qc_all.csv`

| Column | Description |
|---|---|
| `survey`, `sword_version`, `node_id` | Join keys |
| `keep` | 1 = retain, 0 = exclude |
| `excl_frac` | Fraction inside the exclusion polygon (route A only) |

**Join:** `nodes_all` ⋈ `node_qc_all` on (`survey`, `sword_version`, `node_id`),
then filter `keep == 1`.

---

## 11. Verification

The pipeline was validated against synthetic inputs with analytically known
answers, run end to end through the patched RiverObs.

| Check | Result |
|---|---|
| Rasterized area vs. true polygon area (3 m grid) | 479,880 m² vs. 480,000 m² — **0.025% error** |
| Pixel-cloud conversion | Correct lat/lon, `pixel_area` = grid cell area, netCDF round-trip clean |
| Segmentation | Exactly 1 connected component for a single channel |
| **Node width vs. known channel width** | **119.925 m against a true 120.000 m** |
| `area_total / p_length == width` | Exact |
| Area conservation | `sum(area_total)` = **100.0%** of rasterized mask area |
| Step 3b polygons vs. RiverTile `area_total` | Match to **0.000 m²** |
| Prior database with 35 node / 63 reach variables absent | Identical widths to a complete database |
| Prior database with root-level, non-standard dimension names | Identical widths |
| Prior database with ragged centerline counts (15 / 25 / 30) | Correct after the shim; `centerline_lat` rectangular, real point counts recovered |
| Edge behaviour | Nodes where the mask stops mid-node return ~51% low — the artefact the manual review exists to catch |

To reproduce the width verification, rasterize a rectangle of known dimensions
along a synthetic reach and confirm `width ≈ known width` and
`sum(area_total) ≈ polygon area`.

---

## 12. Design decisions

### 12.1 Why `calval2rivertile.py` rather than a bespoke RiverObs harness

RiverObs ships a cal/val entry point that runs the same `SWOTRiverEstimator`,
the same `ReachExtractor`, and the same `process_reaches()` as the operational
processor. Using it means the reference and the SWOT product share the pixel-to-
node path by construction, rather than by careful re-implementation. It also
keeps the modification surface to one filled-in stub and one variable list.

### 12.2 Why there is no automated coverage gate

A footprint derived from an orthomosaic's nodata mask **cannot see cloud** —
cloud is valid data. An automated gate built on it would pass exactly the nodes
most in need of removal, which is worse than no gate because it looks like the
problem was handled. Node exclusion is therefore manual and recorded explicitly.

Two alternatives that do not work, for the record:

- **The water mask as footprint.** The coverage fraction reduces to
  `width / (2 × cross_frac × max_width)` — every `node_length` cancels — so the
  gate becomes a filter on the very width being validated, preferentially
  dropping narrow nodes. Measured on a test channel: 0 of 21 nodes pass, median
  coverage fraction 0.289, exactly `120 m ÷ 400 m`.
- **A study-area or floodplain polygon.** It describes where the floodplain is,
  not where an individual survey flew, so it passes every node in the domain. In
  a test where a survey covered 6 km of a 12 km reach chain, a per-survey
  footprint correctly rejected the 2 half-covered nodes while the domain-wide
  polygon passed all 61.

Nodes a survey never touched need no gate — they receive no water cells and
never reach the output. Only partial cases at the edge and under cloud are in
question, and both are judgement calls.

### 12.3 Why the shim rather than more source edits

`SWOTRiverEstimator.py`, `Estimate.py` and `ReachDatabase.py` differ between
RiverObs revisions, so exact-text patches are fragile. Runtime patches keyed on
class and function names are not. The shim route was verified to produce a
bit-identical RiverTile to an equivalent set of source edits.

---

## 13. Known limitations

- **Heights are not meaningful.** This pipeline is for area and width only.
- **Reference quality is the digitizer's.** Nothing here validates the mask
  itself.
- **Prior-database centerline misalignment is preserved** in both the reference
  and the SWOT product. This is intentional when the goal is to evaluate the
  instrument rather than the hydrography — but it means neither dataset is free
  of it, and the comparison is like-for-like rather than absolute.
- **Braided and anabranching channels** exercise the segmentation and
  cross-channel corridor logic hardest. Stock RiverObs will exclude threads that
  lie outside the prior channel. Always run `util_assignment_audit.py` in these
  settings and tune the corridor factors ([§2.3](#23-the-search-corridor)) before
  trusting a width.
- **The corridor factors are a global multiplier**, applied to every reach in a
  run. A survey containing both a single-thread reach and a wide braidplain is
  tuned to the braidplain; the single-thread reach then carries a wider corridor
  than it needs. This is only a risk where non-river water is present, which by
  assumption it is not in a hand-digitized river mask — but it is the reason to
  use the smallest factor that reaches ~100% rather than a large one.
- **Not tested at sub-metre grid spacing.** See [§8, step 1](#step-1--rasterize-the-water-mask).
- **Reaches dropped mid-run** (ghost reaches, single-node reaches) are lightly
  exercised.

---

## 14. Adapting to a new study area

1. **`config.sh`** — paths, `RES`, `VERSIONS`, and the survey manifest.
2. **`step2_run_riverobs.py`** — the `SURVEYS` list (names + observation dates)
   and `SWORD_VERSIONS` (version tag → `.rdf` filename).
3. **`.rdf`** — `reach_db_path` for your region.
4. **CRS** — `--crs` on steps 1 and 3c, default `EPSG:32606`. Any projected CRS
   with metre units works; the converter reads the CRS from the raster.
5. **`--pad`** — must exceed the largest prior `max_width` in the area.
6. **`--res`** — 3–5 m for channels of tens to hundreds of metres. Scale up for
   very wide rivers; verify against a polygon of known area either way.

Nothing in the pipeline is specific to orthomosaics. Any binary water mask in a
projected CRS works: drone imagery, classified satellite scenes, lidar-derived
water extent.

---

## 15. Citation and license

If you use this pipeline, please cite RiverObs and the SWORD prior database in
addition to this repository:

- RiverObs — https://github.com/SWOTAlgorithms/RiverObs
- Altenau, E. H., et al. (2021), The Surface Water and Ocean Topography (SWOT)
  Mission River Database (SWORD): A global river network for satellite data
  products, *Water Resources Research*, 57, e2021WR030054.

> Verify the SWORD citation against the version you use before publishing.

Patches in `apply_riverobs_patches.py` and `riverobs_shim.py` address defects in
the RiverObs cal/val code path that are not specific to this application; they
are suitable for upstream contribution.

**License:** add a `LICENSE` file. If you intend to redistribute the RiverObs
patches, check compatibility with the RiverObs license (Caltech / JPL).

---

## Appendix — troubleshooting

| Symptom | Cause and fix |
|---|---|
| `ModuleNotFoundError: No module named 'RiverObs'` | `PYTHONPATH` not set in this shell. |
| `No module named 'bottleneck'` / `'piecewise_regression'` | Missing from RiverObs's own package list. Install from `environment.yml`. |
| `TypeError: 'NoneType' object is not subscriptable` in `process_node` | Source patches not applied. Run `apply_riverobs_patches.py --status`. |
| `AttributeError: no attribute 'phase_noise_std'` | `height_agg_method` is not `orig`. |
| `LinAlgError: SVD did not converge` | `min_fit_points` is not set high enough to skip the WSE fits. |
| `ValueError: ... inhomogeneous shape after 1 dimensions` | Shim not loaded, or an older shim. Run jobs through `run_calval2rivertile.py`. |
| `No valid reaches in PRD for this PIXC data` | `reach_db_path` wrong, or step 1's `--pad` too small. |
| `MemoryError` | Grid too large. Increase `--res` or reduce `--pad`. |
| `IndexError` in the discharge model | Prior database `fit_coeffs` has fewer than 3 regions — a database problem, not a pipeline one. |
| Step 3b: `PIXCVec indices exceed the raster shape` | The water `.tif` and PIXCVec came from different runs. |
| Node polygons clip water that is plainly in the mask; `util_assignment_audit.py` reports well under 100% assigned | The cross-channel search corridor is sized from the *prior* channel width, so water outside the prior channel is never offered to a node. Widen it with `RIVEROBS_WTH_COEF_FACTOR` / `RIVEROBS_EXT_DIST_COEF_FACTOR` and re-run step 2. See [§2.3](#23-the-search-corridor). |
| Audit reports `DOUBLE COUNTING` | The corridor factors are wide enough that neighbouring reaches claim the same water. Each reach is processed against the full pixel cloud independently. Lower the factors to the smallest value that still plateaus at ~100% assigned. |
