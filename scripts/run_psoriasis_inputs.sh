#!/usr/bin/env bash
# Build every guidebook input for the psoriasis dataset (GEO GSE173706), in order.
# Needs ftp.ncbi.nlm.nih.gov reachable and `pip install -e ".[scrna]"`. Data stays out of Git.
set -euo pipefail
RAW=data/raw/GSE173706
OUT=data/psoriasis

# A1: download and build one raw-count AnnData (check $OUT/GSE173706_raw.samples.tsv afterwards)
python scripts/download_gse173706.py download "$RAW"
python scripts/download_gse173706.py build "$RAW" "$OUT/GSE173706_raw.h5ad" ${SAMPLES:+--samples "$SAMPLES"}

# A2: QC, integration by donor, two-level marker annotation, frequency tables
python scripts/sc_annotate.py "$OUT/GSE173706_raw.h5ad" "$OUT/annotated" \
  --markers resources/markers/skin_v1.yaml --sample-key sample_id --batch-key donor_id --condition-key condition

# A3 and A4 are independent of each other; both read the annotated object
python scripts/sc_trajectory.py "$OUT/annotated/annotated.h5ad" "$OUT/trajectory" \
  --lineage keratinocyte:keratinocyte_basal --lineage fibroblast:fibroblast_SFRP2
python scripts/sc_communication.py "$OUT/annotated/annotated.h5ad" "$OUT/communication"
