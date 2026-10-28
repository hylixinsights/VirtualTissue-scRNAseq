# Portable formats and compatibility

The packaged `resources/*_v2.schema.json` files describe structural contracts.
`read_tissue_bundle`, `validate_config` and `verify_dynamic` enforce additional
semantic constraints: hashes, finite values, permitted lineages/states, one
biopsy, rate bounds, schedule order and legal transitions. JSON Schema alone
cannot verify biology, identity alignment or the journal hash chain.

## Input ZIP

Exactly `input.json` and `manifest.json`; total uncompressed size <=50 MB.
The manifest contains SHA-256 of input.json and schema_version 2. Each cell has
an immutable cell_id, lineage/cell_type, sample_id, nullable patient_id/clone_id,
cd8_state, molecular_source_id, initial observed RNA marker counts and observed
rank axes. The full RNA matrix remains at the source, identified by SHA-256.
`observed_rna_t0` contains selected **observed** counts; it is never generated
from a future state. `reference_id` names an assumed simulated profile.
Source uncertainty and licences remain in provenance. Unknown is explicit.

## Scenario YAML

Generate with `vtsc melanoma-template scenario.yaml`. All template keys are
required; no silent deep merge. Mechanics use hours and micrometers, hazards
1/hour, influx cells/hour and exposure a normalized fraction. Every run carries
config (base contract), effective_model (control modifications) and
parameter_ledger. No kinetic rate is a clinical estimate by default.
The product library freeze is `src/vtsc/resources/FREEZE.json`; it is independent
of the protected clinical preregistration freeze.

## Recording and results ZIP

`run.json` stores schema_version 2, kind `vtsc_melanoma_run`, source hashes,
initial input, runtime environment, catalog, frames, decisions, events and summary.
`frame_columns` is authoritative: index, x, y, radius, alive, action, target,
damage, fatigue, cycle, CD8 state, state age, dysfunction, proliferation,
cytotoxicity and persistence. Cleared corpses leave frames but remain in catalog
and event history. Index is a stable run-local address, distinct from source ID.

Decisions have actual options, chosen action/outcome, context, seed and nested
state_decision. This records hazards/probabilities, draw, sampled_choice,
provider choice, prior/new state and before/after age. State choice happens
before mechanical action, so prior-state and post-transition checkpoint values
can differ. Hashes chain all decisions; ZIP manifest hashes every exported file.
`events.jsonl`, `decisions.jsonl` and `trajectory.tsv` are convenient projections.

A comparison recording has kind `vtsc_melanoma_comparison`, two embedded runs,
shared time horizon and explicit effective-config differences. Input, seed,
base config and layout must match; initial geometry is checked. The renderer
uses the latest frame at or before the shared clock. HTML embeds a gzip/base64
copy and uses the browser DecompressionStream API; no network/CDN is required.
Use a current Chromium, Firefox or Safari browser.

Recovery verifies the source/decision chain and recomputes from recorded input,
seed and config into a new destination, checking both frame and decision prefixes.
It is deterministic replay, not a serialized hidden process checkpoint. Extending
a partial final time step can change discretization and is rejected if the prefix
differs. Source-version changes are rejected when hashes are present.

## Limits

Initial population 2–2000, created population <=5000, horizon <=96 hours,
dt 0.05–0.5 hours, default budget 100000 cell frames (configurable to 200000).
The runtime stops explicitly at the population/recording cap. Check complete,
stop_reason and observed_hours before interpreting results. Group examples use
120 synthetic or 400 real cells for 24 h. Input preparation may use substantially
more memory than the small simulation. Resource measurements are in VALIDATION.

Legacy v1 tissue bundles remain readable for `run-tissue`; they are rejected by
`run-melanoma`. Reprepare their raw data with explicit CD8 annotation to migrate.
Legacy v0.2 nonspatial and v0.3 tissue runs are preserved, not relabeled as v0.4.
