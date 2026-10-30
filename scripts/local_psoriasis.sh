#!/usr/bin/env bash
# Run the psoriasis input pipeline on your own computer (macOS, Linux, or Windows WSL/Git Bash).
#
#   bash scripts/local_psoriasis.sh            # set up, download GSE173706 (~250 MB), build, annotate,
#                                              # trajectory, communication
#   bash scripts/local_psoriasis.sh download   # only download + build data/psoriasis/GSE173706_raw.h5ad
#
# Needs Python 3.12 or 3.13 and internet access to ftp.ncbi.nlm.nih.gov. Safe to rerun:
# downloads resume, finished steps are reused. Everything lands under data/ (git-ignored).
# Finds python3.13/python3.12 (or lets uv install 3.12); set PYTHON=/path/to/python3.12 to choose.
set -euo pipefail
cd "$(dirname "$0")/.."
VENV=.venv-scrna

ok_version() {  # true if interpreter $1 is Python 3.12 or 3.13
  "$1" -c 'import sys; sys.exit(0 if (3, 12) <= sys.version_info[:2] < (3, 14) else 1)' 2>/dev/null
}

find_python() {
  local c
  for c in ${PYTHON:-} python3.13 python3.12 python3 python; do
    if command -v "$c" >/dev/null 2>&1 && ok_version "$(command -v "$c")"; then
      command -v "$c"; return 0
    fi
  done
  if command -v uv >/dev/null 2>&1; then  # uv can fetch a standalone Python without admin rights
    uv python install 3.12 >&2 && uv python find 3.12 && return 0
  fi
  return 1
}

VPY="$VENV/bin/python"; [ -x "$VPY" ] || VPY="$VENV/Scripts/python.exe"
if [ -x "$VPY" ] && ! ok_version "$VPY"; then
  echo "Removing $VENV: it was created with $("$VPY" --version 2>&1), which is too old."
  rm -rf "$VENV"
fi
if [ ! -x "$VENV/bin/python" ] && [ ! -x "$VENV/Scripts/python.exe" ]; then
  if ! PY=$(find_python); then
    cat >&2 <<'MSG'
This pipeline needs Python 3.12 or 3.13, and none was found (macOS ships 3.9).
Install one of these, then rerun this script:
  - uv (recommended, no admin rights):  curl -LsSf https://astral.sh/uv/install.sh | sh
    then open a new terminal; the script will let uv download Python 3.12 itself
  - Homebrew:                           brew install python@3.12
  - python.org installer:               https://www.python.org/downloads/
Or point to an interpreter explicitly:   PYTHON=/path/to/python3.12 bash scripts/local_psoriasis.sh
MSG
    exit 1
  fi
  echo "Using $PY ($("$PY" --version 2>&1))"
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
