"""Build the three reviewable v0.5 notebooks (English)."""
from pathlib import Path
import nbformat as n
root=Path(__file__).resolve().parents[1]
def md(s):return n.v4.new_markdown_cell(s)
def code(s):return n.v4.new_code_cell(s)
def save(name,cells):
 nb=n.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}})
 n.write(nb,root/'notebooks'/name)
setup='''from pathlib import Path
import os, sys, subprocess, json
INSTALL_PACKAGE = True
# For a local checkout: pip install -e ".[test,notebooks]" first.
# No untested hosted clone URL is executed by this notebook.
if INSTALL_PACKAGE:
    subprocess.check_call([sys.executable, "-m", "pip", "install", "-e", ".[notebooks]"])
WORK = Path(os.environ.get("VTSC_NOTEBOOK_WORK", "outputs/notebook-v05"))
WORK.mkdir(parents=True, exist_ok=True)
'''
save('Part1_Prepare_scRNAseq.ipynb',[
 md('# Part 1 · One melanoma biopsy → versioned input\n\nExploratory software, not validated patient prediction. Raw RNA stays at t0; inferred CD8 labels and assumed kinetic profiles are separate. Use a new WORK directory for each accepted preparation.'),
 code(setup),
 code('''DATA_SOURCE = "melanoma"  # "synthetic" or "own"
OWN_INPUT = Path("my_biopsy.h5ad")
OWN_METADATA = None  # optional TSV with unique cell_id, cell_type, cd8_state, sample_id
OWN_LAYER = "counts"  # use "raw" for .raw.X, or None if X already contains raw counts
OWN_SAMPLE = None  # required when the source contains more than one biopsy
OWN_STATE_RECIPE = None  # explicitly choose "rank-markers-v1" if no reviewed states
SEED = 21
MAX_CELLS = 120
from vtsc.tissue_inputs import download_melanoma, synthetic_melanoma, prepare_tissue
bundle = WORK / "input.zip"
if DATA_SOURCE == "melanoma":
    record = download_melanoma(WORK / "source")
    provenance = prepare_tissue(record["input"], bundle, recipe="melanoma-markers", max_cells=MAX_CELLS, seed=SEED)
elif DATA_SOURCE == "synthetic":
    source = WORK / "synthetic_melanoma.h5ad"
    synthetic_melanoma(source, cells=MAX_CELLS, seed=SEED)
    provenance = prepare_tissue(source, bundle, max_cells=MAX_CELLS, seed=SEED, source_license="MIT; generated fixture")
elif DATA_SOURCE == "own":
    provenance = prepare_tissue(OWN_INPUT, bundle, metadata=OWN_METADATA, layer=OWN_LAYER, sample_id=OWN_SAMPLE,
                                state_recipe=OWN_STATE_RECIPE, max_cells=MAX_CELLS, seed=SEED)
else:
    raise ValueError("Choose melanoma, synthetic or own")
print(json.dumps({k:provenance[k] for k in ["sampled_composition", "sampling", "state_source", "source_license"]}, indent=2))
'''),
 code('''audit = provenance["program_audit"]
print(json.dumps({"CD8_states":audit["cd8_state_counts"], "held_out_markers":audit["heldout_markers"], "coverage":audit["coverage"]}, indent=2))
(WORK / "preparation-report.json").write_text(json.dumps(provenance, indent=2))
'''),
 md('Review Unknown and coverage. Marker concordance is a diagnostic, not validated accuracy. Missing Cytotoxic or Dysfunctional cells are not added. The 10x source is CC BY 4.0; its metadata discrepancy is recorded. SPECIAL exports remain separate and restricted; no paired outcome analysis is performed.')])
save('Part2-3_Run_Simulation.ipynb',[
 md('# Parts 2–3 · Dynamic CD8 and explicit anti-PD-1\n\nMemory, Activated, Cytotoxic and Dysfunctional share a hazard kernel with the nonspatial solver. Every number in the default model is an assumption. Exposure is normalized, not a clinical dose.'),
 code('''from pathlib import Path
from copy import deepcopy
import os, json, yaml
from vtsc.melanoma import DEFAULT_CONFIG
from vtsc.melanoma_results import run_melanoma, compare_melanoma, preflight_melanoma
WORK = Path(os.environ.get("VTSC_NOTEBOOK_WORK", "outputs/notebook-v05"))
SEED = 21
HOURS = 12
DECISION_MODE = "mock-jev"  # free contract test; "local" or "jev" for an explicitly authorized live run
ALLOW_PAID = False         # live Jev requires True, an API key in environment, and a request ceiling
MAX_HTTP_ATTEMPTS = 500    # cumulative HTTP-attempt cap, not a currency budget
JEV_MODEL = "jev-1.14.0"
config = deepcopy(DEFAULT_CONFIG)
config["model"]["schedule"] = [{"start_hour":2, "end_hour":96, "exposure":0.8}]
config_path = WORK / "assumptions.yaml"
if config_path.exists():
    raise FileExistsError("Use a new WORK directory; do not overwrite an accepted configuration")
config_path.write_text(yaml.safe_dump(config, sort_keys=False))
print(config["model"]["parameter_provenance"])
'''),
 code('''preflight = preflight_melanoma(WORK / "input.zip", WORK / "preflight", config_path=config_path, seed=SEED)
print(preflight)
summaries = {}
for arm in ["anti-pd1", "exposure-zero", "recognition-off", "transitions-off"]:
    summaries[arm] = run_melanoma(WORK / "input.zip", WORK / arm, config_path=config_path, arm=arm, hours=HOURS, seed=SEED, decision_mode=DECISION_MODE, allow_paid=ALLOW_PAID, max_http_attempts=MAX_HTTP_ATTEMPTS, jev_model=JEV_MODEL)
print(json.dumps(summaries, indent=2))
comparison = compare_melanoma(WORK / "anti-pd1", WORK / "exposure-zero", WORK / "comparison")
'''),
 code('''from vtsc.dynamics import solve, STATES
from vtsc.tissue_inputs import read_tissue_bundle
p = read_tissue_bundle(WORK / "input.zip")
initial = [sum(c["cd8_state"] == s for c in p["cells"]) for s in STATES]
context = {"contact":1.,"target":1.,"reactive":1.,"tumor_fraction":.5,"myeloid_fraction":.1}
expected = solve(initial, [0, 6, 12], config["model"], context)
print("Homogeneous hypothetical context, not a prediction for the spatial run:")
print(dict(zip(STATES, expected[-1])))
'''),
 md('Compare multiple seeds with `vtsc ensemble-melanoma`. Recognition-off disables recognition, while exposure-zero preserves recognition and checkpoint inhibition. Checkpoint-off decouples the drug. No seed is selected by apparent benefit. Group-level variation is simulation uncertainty, not a patient confidence interval.')])
save('Part4_Explore_Results.ipynb',[
 md('# Part 4 · Verify, inspect and recover\n\nThe viewer replays recorded frames, behavioral choices and separate stochastic CD8 fate decisions. It does not modify the engine or synthesize new expression. Open the comparison HTML in a current browser; it is self-contained and requires no network.'),
 code('''from pathlib import Path
import os, json
from IPython.display import FileLink, display
from vtsc.results import read_run
from vtsc.melanoma_results import verify_dynamic, write_melanoma_html, recover_melanoma
WORK = Path(os.environ.get("VTSC_NOTEBOOK_WORK", "outputs/notebook-v05"))
result = verify_dynamic(read_run(WORK / "anti-pd1.zip"))
comparison = verify_dynamic(read_run(WORK / "comparison.zip"))
reopened = WORK / "comparison-reopened.html"
write_melanoma_html(comparison, reopened)
display(FileLink(str(reopened)))
print(json.dumps(result["summary"], indent=2))
'''),
 code('''transitions = [d for d in result["decisions"] if d["state_decision"] and d["state_decision"]["choice"].startswith("transition:")]
if transitions:
    print(json.dumps(transitions[0], indent=2))
else:
    print("No state transition was sampled in this finite run. Inspect probabilities; do not pick a seed to force a desired event.")
'''),
 code('''# Deterministic recovery recomputes into a fresh directory and verifies the saved prefix.
recovery = recover_melanoma(WORK / "anti-pd1.zip", WORK / "recovered")
assert read_run(WORK / "recovered")["frames"] == result["frames"]
print(recovery["stop_reason"])
'''),
 md('The inspector exposes state age, capacities, allowed options, probabilities, draw, protein assumptions, contacts and hash-chained journal. The original SPECIAL study still needs eligibility, training-only reference, fit/identifiability, LOPO, lockbox and external/spatial validation. These software runs do not pass those gates.')])
