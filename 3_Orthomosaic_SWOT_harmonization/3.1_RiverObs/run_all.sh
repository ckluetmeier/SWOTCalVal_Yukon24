#!/usr/bin/env bash
# =============================================================================
# run_all.sh -- steps 1, 2, 3b and 3c for all six orthomosaics
# -----------------------------------------------------------------------------
#   conda activate RiverObs
#   cd .../SWOTCalVal_Yukon24/3_Orthomosaic_SWOT_harmonization/3.1_RiverObs
#   bash run_all.sh                  # everything
#   bash run_all.sh 1                # step 1 only
#   bash run_all.sh 2 3b             # steps 2 and 3b only
#   bash run_all.sh 3c               # just re-apply the exclusion layer
#
# =============================================================================

set -uo pipefail

# Paths
BASE=/Users/camryn/Documents/UNC/_Tier1_sites/expanded_Yukon_Flats
SCRIPTS=/Users/camryn/Documents/UNC/_Tier1_sites/_data_management/YR2024_scripts/SWOTCalVal_Yukon24/3_Orthomosaic_SWOT_harmonization/3.1_RiverObs
RIVEROBS=/Users/camryn/Documents/UNC/_Tier1_sites/_data_management/YR2024_scripts/RiverObs

MASK_DIR=$BASE/Orthomosaics/width_validation
ORTHO_DIR=$BASE/Orthomosaics/orthos_resampled_3m/gap_filled
RASTER_DIR=$BASE/Orthomosaics/width_validation/rasterized_water_masks
OUT_DIR=$BASE/Orthomosaics/width_validation/RiverObs_output
QC_DIR=$BASE/Orthomosaics/width_validation/RiverObs_output/node_qc

# Hand-digitized polygons over cloud, haze and ragged survey edges. One layer
# covering all six orthos. The test is spatial, so it is valid for v16 and v17b
# alike -- node ids differ between SWORD versions, the ground does not.
EXCLUDE_SHP=$BASE/Orthomosaics/width_validation/digitizing/bad_ends_watermasks.shp

# If a polygon in EXCLUDE_SHP should apply to only ONE survey (two surveys
# covering the same reach on different dates, with different bad areas), add a
# text field to the layer naming the survey and set EXCLUDE_FIELD to its name.
# Leave empty to apply every polygon to every survey.
EXCLUDE_FIELD=""

# Fraction of a node polygon inside the exclusion area that is tolerated
# before the node is dropped.
MAX_OVERLAP=0.05

RES=3
VERSIONS="v17b v16"
LOG_LEVEL=info

# Search-corridor relaxation
# RiverObs sizes its cross-channel search corridor from the PRIOR channel
# width, so water lying well outside the prior channel -- anabranches, braidplain
# threads, wide side channels -- is excluded even when it is unambiguously river.
# These factors widen the corridor. 1.0 is stock RiverObs.
#
# Raise the factors until the fraction of the water mask assigned to a node
# plateaus near 100%, then stop; step 3b reports that fraction per survey. A
# wider corridor than you need lets adjacent reaches claim the same water.
export RIVEROBS_WTH_COEF_FACTOR=8.0
export RIVEROBS_EXT_DIST_COEF_FACTOR=8.0

# Factors used per survey: CD = 8, upperYR = 3, upperPR_CL = 3

# Survey definitions
# name | water mask basename | orthomosaic basename
# Names must match the SURVEYS list at the top of 3.1.2_run_RiverObs.py.
SURVEYS="
CD_071024|chandalar_240710_ortho_25cm_watermask|chandalar_240710_ortho_3m_epsg32606
upperPR_CL_071024|coleen_240710_ortho_25cm_watermask|coleen_240710_ortho_3m_epsg32606
upperPR_CL_071624|coleen_240716_ortho_25cm_watermask|coleen_240716_ortho_3m_epsg32606
lowerPR_SJ_072624|sheenjek_240725_ortho_25cm_watermask|sheenjek_240725_ortho_3m_epsg32606
lowerYR_071624|yukonDS_240716_ortho_25cm_watermask|yukonDS_240716_ortho_3m_epsg32606
upperYR_071024|yukonUS_240710_ortho_25cm_watermask|yukonUS_240710_ortho_3m_epsg32606
"

STEPS="${*:-1 2 3b 3c}"
want () { case " $STEPS " in *" $1 "*) return 0;; *) return 1;; esac; }

mkdir -p "$RASTER_DIR" "$OUT_DIR" "$QC_DIR"
cd "$SCRIPTS" || exit 1

FAILED=""

# =============================================================================
# STEP 1 -- rasterize each water mask to a 3 m grid
# =============================================================================
if want 1; then
echo
echo "################ STEP 1: rasterize water masks ################"
echo "$SURVEYS" | while IFS='|' read -r NAME MASK ORTHO; do
  [ -z "$NAME" ] && continue
  echo
  echo "--- $NAME ---"
  python 3.1.1_rasterize_watermask.py \
      "$MASK_DIR/$MASK.shp" \
      "$RASTER_DIR" \
      --name "$NAME" \
      --res "$RES" \
      --ortho-tif "$ORTHO_DIR/$ORTHO.tif" \
    || echo "!!! step1 FAILED for $NAME"
done
fi

# =============================================================================
# STEP 2 -- run RiverObs
# =============================================================================
if want 2; then
echo
echo "################ STEP 2: RiverObs ################"
python 3.1.2_run_RiverObs.py \
    --riverobs-root "$RIVEROBS" \
    --raster-dir "$RASTER_DIR" \
    --rdf-dir . \
    --out-dir "$OUT_DIR" \
    --res "$RES" \
    --versions $VERSIONS \
    --keep-going \
    --wth-coef-factor "$RIVEROBS_WTH_COEF_FACTOR" \
    --ext-dist-coef-factor "$RIVEROBS_EXT_DIST_COEF_FACTOR" \
    --log-level "$LOG_LEVEL"
STEP2_RC=$?
[ $STEP2_RC -ne 0 ] && FAILED="$FAILED step2(rc=$STEP2_RC)"
fi

# =============================================================================
# STEP 3b -- node polygons for review
# =============================================================================
if want 3b; then
echo
echo "################ STEP 3b: node polygons ################"
for V in $VERSIONS; do
  echo "$SURVEYS" | while IFS='|' read -r NAME MASK ORTHO; do
    [ -z "$NAME" ] && continue
    WATER="$RASTER_DIR/${NAME}_water_${RES}m.tif"
    PIXCVEC="$OUT_DIR/$V/ortho_${NAME}_${V}_pixcvec.nc"
    RIVERTILE="$OUT_DIR/$V/ortho_${NAME}_${V}.nc"
    OUT_SHP="$QC_DIR/true_node_polygons_${NAME}_${V}.shp"

    if [ ! -f "$RIVERTILE" ]; then
      echo "--- $NAME / $V : SKIPPED, no RiverTile (step 2 failed for this one)"
      continue
    fi
    echo
    echo "--- $NAME / $V ---"
    python 3.1.3_true_node_polygons.py \
        "$WATER" "$PIXCVEC" "$OUT_SHP" \
        --rivertile "$RIVERTILE" \
      || echo "!!! step3b FAILED for $NAME / $V"
  done
done
fi

# =============================================================================
# STEP 3c -- apply the hand-drawn exclusion layer to every survey and version
# =============================================================================
if want 3c; then
echo
echo "################ STEP 3c: apply exclusion layer ################"
if [ ! -f "$EXCLUDE_SHP" ]; then
  echo "SKIPPED: no exclusion layer at"
  echo "  $EXCLUDE_SHP"
  echo "Digitize one in QGIS over the orthomosaics, then re-run: bash run_all.sh 3c"
  FAILED="$FAILED step3c(no-exclusion-layer)"
else
  FIELD_ARG=""
  [ -n "$EXCLUDE_FIELD" ] && FIELD_ARG="--exclude-field $EXCLUDE_FIELD"
  python 3.1.4_manual_node_qc.py --batch \
      --nodes-dir "$QC_DIR" \
      --exclude-shp "$EXCLUDE_SHP" \
      --versions $VERSIONS \
      --max-overlap "$MAX_OVERLAP" \
      --out-dir "$QC_DIR" \
      --out "$QC_DIR/node_qc_all.csv" \
      $FIELD_ARG
  STEP3C_RC=$?
  [ $STEP3C_RC -ne 0 ] && FAILED="$FAILED step3c(rc=$STEP3C_RC)"
fi
fi

# =============================================================================
echo
echo "############################################################"
if [ -n "$FAILED" ]; then
  echo "Finished with problems:$FAILED"
  echo "Scroll up for '!!!' and 'FAILED' lines."
else
  echo "Finished."
fi
echo
echo "Outputs:"
echo "  water rasters   $RASTER_DIR"
echo "  RiverObs        $OUT_DIR"
echo "     combined CSVs: ortho_riverobs_nodes_all.csv, ortho_riverobs_reaches_all.csv"
echo "  node polygons   $QC_DIR/true_node_polygons_<survey>_<version>.shp"
echo "  node QC         $QC_DIR/node_qc_all.csv   <- join downstream analysis on this"
echo
echo "Corridor factors in effect: wth_coef x$RIVEROBS_WTH_COEF_FACTOR,"
echo "                            ext_dist_coef x$RIVEROBS_EXT_DIST_COEF_FACTOR"
echo "Step 3b reports, per survey, how much of each water mask reached a node."
echo
echo "Edited the exclusion layer in QGIS? Re-run step 3c alone:"
echo "  bash run_all.sh 3c"
echo "############################################################"
