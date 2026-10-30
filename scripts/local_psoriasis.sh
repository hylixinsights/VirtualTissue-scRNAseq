#!/usr/bin/env bash
# Run the psoriasis input pipeline on your own computer (macOS, Linux, or Windows WSL/Git Bash).
#
#   bash scripts/local_psoriasis.sh            # set up, download GSE173706 (~250 MB), build, annotate,
#                                              # trajectory, communication
#   bash scripts/local_psoriasis.sh download   # only download + build data/psoriasis/GSE173706_raw.h5ad
#
# Needs Python 3.12 or 3.13 and internet access to ftp.ncbi.nlm.nih.gov. Safe to rerun:
# downloads resume, finished steps are reused. Everything lands under data/ (git-ignored).
# Set PYTHON=/path/to/python3.12 to pick an interpreter.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PYTHON:-python3}
VENV=.venv-scrna

if [ ! -x "$VENV/bin/python" ] && [ ! -x "$VENV/Scripts/python.exe" ]; then
  "$PY" -m venv "$VENV"
fi
VPY="$VENV/bin/python"; [ -x "$VPY" ] || VPY="$VENV/Scripts/python.exe"
"$VPY" -m pip install -q --upgrade pip
"$VPY" -m pip install -q -e ".[scrna]"

RAW=data/raw/GSE173706
OUT=data/psoriasis
"$VPY" scripts/download_gse173706.py download "$RAW"
"$VPY" scripts/download_gse173706.py build "$RAW" "$OUT/GSE173706_raw.h5ad" ${SAMPLES:+--samples "$SAMPLES"}
echo "Check the sample table: $OUT/GSE173706_raw.samples.tsv"
[ "${1:-all}" = download ] && exit 0

"$VPY" scripts/sc_annotate.py "$OUT/GSE173706_raw.h5ad" "$OUT/annotated" \
  --markers resources/markers/skin_v1.yaml --sample-key sample_id --batch-key donor_id --condition-key condition
"$VPY" scripts/sc_trajectory.py "$OUT/annotated/annotated.h5ad" "$OUT/trajectory" \
  --lineage keratinocyte:keratinocyte_basal --lineage fibroblast:fibroblast_SFRP2
"$VPY" scripts/sc_communication.py "$OUT/annotated/annotated.h5ad" "$OUT/communication"

# small bundle to share back (tables and reports, no expression matrices)
"$VPY" - <<'PYEOF'
import zipfile, pathlib
out = pathlib.Path("data/psoriasis")
with zipfile.ZipFile(out / "psoriasis_results_small.zip", "w", zipfile.ZIP_DEFLATED) as z:
    for p in sorted(out.rglob("*")):
        if p.is_file() and p.suffix in {".tsv", ".md", ".json"}:
            z.write(p, p.relative_to(out))
    for p in ("data/raw/GSE173706/SHA256SUMS", "data/raw/GSE173706/SOURCE.json"):
        z.write(p, pathlib.Path(p).name)
print(f"wrote {out / 'psoriasis_results_small.zip'}")
PYEOF
