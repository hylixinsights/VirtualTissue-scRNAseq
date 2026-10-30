"""Live Jev (TypeSafe) smoke test for the psoriasis context.

Reads the key from TYPESAFE_API_KEY (or JEV_API_KEY) only. Sends a small, fixed
set of Choice requests and writes the parsed answers to an output directory.
The key, cookies and response headers are never written.

    TYPESAFE_API_KEY=... python scripts/jev_smoke_test.py reports/jev-smoke/2026-10-30

Eight requests, five billed (~3.5k input tokens). Error probes (bad key,
unknown model, missing field) are rejected before inference.
"""
import json
import os
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ENDPOINT = "https://api.typesafe.ai/v1/systemone"
MODEL = "jev-1.14.0"
OPTIONS = {
    "proliferate": "Enter the cell cycle and divide if space allows.",
    "differentiate": "Progress toward terminal cornified differentiation and move suprabasally.",
    "secrete_inflammatory": "Secrete antimicrobial peptides and chemokines (e.g. S100A7/8/9, CXCL1/8, CCL20) that recruit immune cells.",
    "rest": "Maintain current state until the next checkpoint.",
}
SEMANTICS = ("Observations are inferred RNA summaries, not protein activity. Choose a permitted "
             "behavior; probabilities rank choices, not kinetic rates.")
LESIONAL = "Psoriatic lesional skin (plaque), whole-tissue scRNA-seq reconstruction; positions assumed."
NONLESIONAL = "Non-lesional skin from a psoriasis patient, whole-tissue scRNA-seq reconstruction; positions assumed."

BASAL_INFLAMED = {"cell_type": "basal keratinocyte", "neighbors": ["Th17 cell", "dendritic cell"],
                  "exposure": {"IL17A": "high", "IL22": "high", "TNF": "moderate"},
                  "markers_up": ["KRT16", "KRT6A", "MKI67"]}
SPINOUS_QUIET = {"cell_type": "spinous keratinocyte", "neighbors": ["keratinocyte"],
                 "exposure": {"IL17A": "none", "IL22": "low"}, "markers_up": ["KRT10", "KRT1"]}
BASAL_QUIET = {"cell_type": "basal keratinocyte", "neighbors": ["keratinocyte"],
               "exposure": {"IL17A": "none"}, "markers_up": ["KRT5", "KRT14"]}
BASAL_WITH_PRIOR = dict(BASAL_QUIET, guidebook_prior={
    "source": "empirical frequencies from lesional scRNA-seq",
    "proliferate": 0.2, "differentiate": 0.2, "secrete_inflammatory": 0.55, "rest": 0.05})


def question(obs):
    return {"type": "choice",
            "instructions": {"task": "Choose one permitted local behavior for this cell.", "local_observation": obs},
            "criteria": OPTIONS}


def body(setting, observations, model=MODEL, state=True):
    b = {"model": model, "questions": {f"c{i + 1}": question(o) for i, o in enumerate(observations)}}
    if state:
        b["state"] = {"setting": setting, "prompt_version": "vtsc-probe-psoriasis-1", "semantics": SEMANTICS}
    return b


CASES = {
    "lesional_two_cells": (body(LESIONAL, [BASAL_INFLAMED, SPINOUS_QUIET]), "key"),
    "nonlesional_two_cells": (body(NONLESIONAL, [BASAL_INFLAMED, SPINOUS_QUIET]), "key"),
    "lesional_quiet_basal": (body(LESIONAL, [BASAL_QUIET]), "key"),
    "lesional_quiet_basal_repeat": (body(LESIONAL, [BASAL_QUIET]), "key"),
    "lesional_quiet_basal_with_prior": (body(LESIONAL, [BASAL_WITH_PRIOR]), "key"),
    "error_unknown_model": (body(LESIONAL, [BASAL_QUIET], model="jev-0.0.1"), "key"),
    "error_missing_state": (body(LESIONAL, [BASAL_QUIET], state=False), "key"),
    "error_bad_key": (body(LESIONAL, [BASAL_QUIET]), "invalid"),
}


def call(payload, key):
    req = urllib.request.Request(ENDPOINT, data=json.dumps(payload).encode(), method="POST",
                                 headers={"Authorization": "Bearer " + key, "Content-Type": "application/json"})
    start = time.time()
    try:
        with urllib.request.urlopen(req, timeout=180) as r:
            status, raw = r.status, r.read()
    except urllib.error.HTTPError as e:
        status, raw = e.code, e.read()[:4000]
    latency = round(time.time() - start, 2)
    try:
        data = json.loads(raw)
    except ValueError:
        data = raw.decode(errors="replace")[:500]
    if status == 422 and isinstance(data, dict):
        # Validation errors echo the request body; keep only type and location.
        data = {"detail": [{k: d.get(k) for k in ("type", "loc", "msg")} for d in data.get("detail", [])]}
    return status, latency, data


def main(out):
    key = os.environ.get("TYPESAFE_API_KEY") or os.environ.get("JEV_API_KEY")
    if not key:
        sys.exit("Set TYPESAFE_API_KEY (or JEV_API_KEY).")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    results = {}
    for name, (payload, which) in CASES.items():
        status, latency, data = call(payload, key if which == "key" else "apikey_invalid")
        results[name] = {"http_status": status, "latency_s": latency, "request": payload, "response": data}
        print(f"{name}: HTTP {status} in {latency}s")
    text = json.dumps({"endpoint": ENDPOINT, "date_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                       "results": results}, indent=1)
    if key in text:
        sys.exit("Refusing to write: the key appeared in the output.")
    (out / "smoke_results.json").write_text(text + "\n", encoding="utf8")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "reports/jev-smoke/latest")
