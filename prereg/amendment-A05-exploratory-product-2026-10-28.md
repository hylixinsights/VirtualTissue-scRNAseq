# A05 — Prospective exploratory melanoma product contract

Date: 2026-10-28. Authorization: NEW_CONVERSATION_MELANOMA_ANTI_PD1.md,
explicitly adopted by the user. Parent source commit: d85b252.

Build and test an independent melanoma/anti-PD-1 product now. This authorizes
implementation of shared state hazards, estimation infrastructure, spatial agents,
controls, user inputs, offline inspection and synthetic recovery. It does not
waive G1/G2/G3/G4/G5/G6, assign a lockbox, or validate prediction in SPECIAL.
The original FREEZE, frozen plans, study.yaml and T/G board remain unchanged.

Development inputs: independent 10x melanoma, synthetic fixtures and explicitly
supplied single-biopsy inputs. No SPECIAL pre/on state differences or clinical
outcomes will inform rules. No author contact, external publication, push, PR,
paid provider calls or changes to the Xenium engine are authorized.

## Prospective implementation experiments

Before executing new analyses, fix: four states, a labelled assumed graph;
rank-based positive-gene program scoring; separate dysfunction, proliferation,
cytotoxicity and persistence axes; normalized assumed drug exposure. Use seeds
21–25 for example/control ensembles, never choose a seed by treatment effect.
Report 10x full/subsample lineage and CD8-state counts, Unknown, program coverage,
overlap, held-out marker concordance and association with depth/dissociation.
These diagnostics assess a provisional recipe, not a validated classifier.

Synthetic tests cover analytic cases, mass balance, stochastic recovery,
identifiability and common-kernel spatial consistency. The consistency test
replaces local inputs with constants and disables space-limited births and
extra mechanics, reports Monte Carlo uncertainty and time-step convergence.
Sensitivity varies reactivity, checkpoint strength and transition timescale
by fixed factors (0.5, 1, 2); no tuning to improve tumor killing.

## Public TME audit (outcome blind)

Acquire only Entire_TME.rds and Malignant_cells.rds; require published size and
MD5, record SHA-256. Inspect structural metadata names, not clinical values or
state compositions. Export RNA counts aligned with barcodes to new directories;
allow only sample identifiers, total counts/genes/mitochondrial fraction and
broad lineage labels. Archive other metadata separately with access restrictions.
Use CD3D/CD3E/TRAC plus CD8A/CD8B as a fixed provisional CD8 screen (>=2 detected
T markers, >=1 CD8 marker); no state-level or paired analysis. Report ambiguous
lineage; no relabelling to reach G1. Exact sparse counts and ordering must
round-trip. If resources prevent completion, record the failure and continue
independent product work. Lin stays metadata-only until actual matrix alignment
is validated; do not repeat the failed 29 GB conversion.

## Evidence classes

RNA: measured with assay bias. Labels/programs: inferred or user supplied.
State graph, rates, protein checkpoint activities, recognition, geometry and
normalized exposure: assumed. Fitted synthetic parameters: estimated from
synthetic observations. Simulated profiles are capacities, never future RNA.
