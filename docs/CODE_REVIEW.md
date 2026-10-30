# Code review: VirtualTissue-scRNAseq v0.5.0

Reviewed commit: `79d94c6` (2026-10-30). Scope: `src/vtsc/` (about 2,900 lines),
tests, CI, packaging and documentation. The goal used as the yardstick is the one in
[PLAN.md](PLAN.md): a **tissue-agnostic** pipeline that turns an annotated
single-cell object plus a declared context into a per-cell-type guidebook, and
simulates a virtual tissue in which Jev decides each cell's behavior.

## Summary

The engineering discipline is strong: provenance hashes, a durable and hash-chained
Jev request journal, explicit paid opt-in, zero-cost preflight, deterministic replay,
a clean split between *behavior* (Jev) and *fate* (kernel hazards), and 60 passing
tests (`python -m pytest -q`: 60 passed, 1 skipped, Python 3.13, 2026-10-30).

The main problem is **scope, not quality**. Almost every layer is hard-coded to one
scenario (melanoma + CD8 T cells + anti-PD-1). Psoriasis, nasal mucosa, lymph node
or brain cannot be run without editing source code in at least seven modules. The
repository also mixes a pre-registered research audit (SPECIAL/Lin cohorts) with the
product, which will confuse outside collaborators.

| Priority | Theme | Effort |
|---|---|---|
| High | Remove melanoma/CD8 hard-coding; introduce a tissue/guidebook abstraction | Large |
| High | Send declared context and guidebook priors to Jev; use Jev probabilities, not only the argmax | Medium |
| High | One source of truth for action menus (today there are three) | Medium |
| High | Accept trajectory and CellChat results as inputs | Medium |
| Medium | Split `simulate_melanoma` into small, testable steps | Medium |
| Medium | Separate the research audit from the product package | Small |
| Medium | Loosen exact dependency pins in `pyproject.toml` | Small |
| Medium | Spatial queries are O(N²); use the KD-tree everywhere | Small |
| Low | Code style: dense one-line statements; add a formatter and linter | Small |
| Low | Generalize the viewer and its labels | Medium |

## What works well (keep it)

1. **Behavior / fate separation** (`behavior.py`, `dynamics.py`). Jev ranks only
   kernel-permitted actions; death, division and state transitions stay with the
   kernel. This is the right design for a guidebook: the guidebook defines what a cell
   *can* do, Jev chooses among those options, and physics decides what actually happens.
2. **Safe provider** (`jev_provider.py`). No retries or redirects, write-ahead journal
   with `fsync`, response whitelisting, key redaction, request ceilings, and replay
   from cache. Keep all of this as the generalization proceeds.
3. **Credentials** are read from `TYPESAFE_API_KEY` / `JEV_API_KEY`
   (`melanoma_results.py:98`), never from config. Good.
4. **Input validation** rejects normalized counts, duplicate IDs and mismatched
   metadata (`inputs.py`, `tissue_inputs.py`).
5. **Reproducibility**: keyed RNG per cell/step (`dynamics.keyed_rng`), deterministic
   recovery, checksummed outputs and an offline viewer.
6. **Honest scoping** in docs: assumed parameters are labelled as assumed.

## High priority

### H1. Melanoma, CD8 and anti-PD-1 are hard-coded across the stack

Evidence:

- `tissue_inputs.py:44` fixes the cell-type vocabulary to
  `['Melanoma', 'CD8_T', 'T_cell', 'NK', 'Myeloid', 'Fibroblast', 'Endothelial', 'B_cell', 'Unknown']`;
  `prepare_tissue` rejects any other label (`tissue_inputs.py:169`). A keratinocyte cannot be loaded.
- `tissue_inputs.py:32` `PANELS` and `PROGRAMS` contain melanoma markers only.
- `tissue.py:26-29` `RADII` and `MOBILE` are keyed by those labels.
- `dynamics.py` `STATES`, `PROFILES`, `DEFAULT_MODEL` describe four CD8 states and a PD-1 checkpoint.
  `checked_model` requires exactly these four states.
- `jev_payload.py:21` sends the fixed setting "exploratory melanoma tissue" to Jev;
  `DESCRIPTIONS` has six melanoma actions.
- `jev_payload.py:30` batches by `CD8_T` / `NK` / other.
- `behavior.py:24` `options_for` and `melanoma.py:165` encode "only CD8/NK attack Melanoma",
  "only Myeloid clears", "only Melanoma grows".
- `melanoma.py` `ARMS` are anti-PD-1 specific; CLI help text says "melanoma biopsy".

**Recommendation.** Move all biology out of code into versioned, schema-validated data:

- `context.yaml` (declared by the user: tissue, disease, perturbation, library type; see PLAN.md).
- `guidebook.json` (generated): cell types, states, permitted actions with preconditions,
  context-conditioned action priors, interaction partners, radii, mobility, division ability.
- An **action registry** in code: each action (`migrate`, `proliferate`, `differentiate`,
  `secrete`, `attack`, `phagocytose`, `rest`, `die`, ...) is a small class with
  `is_legal(cell, neighborhood, guidebook)` and `apply(cell, world)`. The guidebook only
  references actions by name, so new tissues do not need new code unless they need a
  genuinely new mechanic.
- Keep the melanoma/CD8 model as the first **scenario pack** (`scenarios/melanoma_anti_pd1/`)
  that the new engine must reproduce bit-for-bit (regression test), then add psoriasis as
  the second pack.

### H2. Jev never sees the user's context, and its probabilities are discarded

- The prompt (`jev_payload.payload`) has a fixed melanoma setting. The tissue, disease and
  perturbation, which are the "prompt" in the psoriasis experiment, are never sent.
- `BehaviorController.choose` (`behavior.py:71-72`) executes `answer['choice']`, the
  argmax. The returned `probabilities` are only logged. With identical observations every
  cell then makes the same choice, so the tissue loses the cell-to-cell heterogeneity
  that a probability guidebook is meant to express.

**Recommendation.**

1. Payload v2 (`PROMPT_VERSION = 'vtsc-behavior-2'`): add a `context` block (tissue,
   disease, perturbation, library type) and, per question, the guidebook entry for that
   cell type and state: allowed actions, prior probabilities, evidence notes and
   the interaction partners currently in contact. Never send raw counts, barcodes or
   patient IDs (keep the existing privacy rule).
2. Add a `sampling` option: `argmax` (current behavior, keep for regression) or
   `sample` (draw from Jev's probabilities with `keyed_rng(seed, cell_id, step, 'behavior')`),
   with an optional mixing weight with the guidebook prior
   (`p = (1-w)·p_jev + w·p_prior`). Record prior, Jev output, mixed distribution and draw
   in `decisions.jsonl`.
3. Add a "guidebook-only" decision mode (sample from priors, no Jev) as the free
   baseline and control arm. This also makes the pipeline usable without an API key.

### H3. Three copies of the action-menu logic

The legal-action rules exist in `tissue.decision_menu` (`tissue.py:88`),
`behavior.options_for` (`behavior.py:24`) and inline in `simulate_melanoma`
(`melanoma.py:165-170`). They already differ: the inline and `tissue.py` menus offer
`divide`, while `options_for` offers `approach` and never `divide`. `local_context` in
`behavior.py` re-implements the neighbor/contact computation from `melanoma.py`.
A change in one place will silently desynchronize the Jev menu and the executed menu.

**Recommendation.** One `legal_actions(cell, neighborhood, guidebook)` function (or the
action registry from H1) used by the preflight, the Jev request builder and the kernel,
plus a test asserting that the offered menu equals the executed menu at the same
tick.

### H4. No path for trajectory or CellChat results

The goal is to use trajectory and cell-cell communication analyses when available.
Neither is accepted today; `prepare_tissue` reads only counts plus `cell_type` /
`cd8_state`.

**Recommendation.** Optional inputs with documented, tool-neutral formats:

- Trajectory: `obs` columns `pseudotime`, `lineage` / `trajectory_branch`, optional
  `terminal_state` probabilities (CellRank) or a TSV of state-to-state transition weights.
  Use them to derive state-transition priors (who can become whom, and roughly how fast
  in pseudotime units), labelled `evidence_class: inferred`.
- CellChat: the long table from `subsetCommunication(cellchat)` (`source, target, ligand,
  receptor, prob, pval, pathway_name`), one per condition. Also accept
  LIANA/CellPhoneDB tables with the same columns. Use them to set interaction partners,
  secreted fields and the strength of neighbor effects in the guidebook.

## Medium priority

### M1. `simulate_melanoma` is a 200-line monolith

`melanoma.py:48-250` handles validation, initialization, field diffusion, fate draws,
behavior execution, influx, recording and summary in one function with nested closures.
It is hard to test any mechanic alone, and a collaborator adding a tissue will have to
edit it. Split it into `World` (cells, fields, KD-tree), `step_fields`, `step_fates`,
`step_behaviors`, `step_demography` and `Recorder`, each unit-tested.

### M2. Research audit and product share one repository root

172 of 320 tracked files are under `reports/`; there are also `prereg/`, cohort-audit
scripts (`scripts/audit_*`, `extract_lin_*`, `decode_qs2_blocks.py`), and four
lock files at the root. A new collaborator cannot tell which commands matter.
Move the pre-registered SPECIAL/Lin work to `research/` (or a separate repository),
keep `scripts/` for product utilities only, and put lock files under `env/`.
`docs/NEW_CONVERSATION_MELANOMA_ANTI_PD1.md` is an AI hand-off prompt and should not
be in user documentation.

### M3. Exact dependency pins in `pyproject.toml`

`pyproject.toml:11` pins `numpy==2.5.4`, `pandas==3.1.0`, `scipy==1.18.1`,
`anndata==0.13.4`. Installing into an existing scanpy/R-reticulate environment will
fail or downgrade packages. Use compatible ranges (`numpy>=2.1,<3`) in
`pyproject.toml` and keep exact pins only in the lock files used by CI and
reproducibility runs.

### M4. Quadratic spatial checks

`tissue.free_disc` (`tissue.py:82`) loops over every cell, and it is called up to 24
times per division attempt and up to five times per move. With 2,000 cells this is
millions of distance computations per step. Query the step's `cKDTree` (rebuilt after
moves, or a uniform grid) for candidates within `2 * max_radius`.

### M5. Input coverage for the target audience

- Seurat users must export to H5AD themselves (documented in `docs/02_your_own_data.md`).
  Provide `r/export_for_vtsc.R` that writes counts, `meta.data`, reductions, and
  optional CellChat tables in one call (via `anndataR` or `zellkonverter`), and validates
  required columns.
- `max_cells` defaults to 400 and caps initial cells at 2,000. That is fine for the
  agent engine, but the **guidebook** must be computed on all cells; only the virtual
  tissue should be subsampled.
- Single-nucleus data (brain) needs a flag: nuclear transcripts under-represent some
  programs, and neuron mobility/division must default to off.
- Enriched libraries (for example CD45+ sorted) must not be used to set tissue
  frequencies without a declared correction (see PLAN.md, context schema).

### M6. Multi-sample and condition structure

Inputs are "one biopsy". Psoriasis needs paired lesional / non-lesional samples from
the same donors, plus healthy controls. Add `sample_id`, `donor_id` and `condition`
as required `obs` columns, compute frequencies per sample, and let the user choose which
sample (or condition mean) initializes the virtual tissue.

## Low priority

- **Style.** `behavior.py`, `jev_payload.py`, `melanoma.py` and `dynamics.py` use
  semicolon-chained statements and very long lines, which makes review hard for
  collaborators. Add `ruff` (lint + format) and a pre-commit hook; enforce in CI.
- **Type hints and docstrings** on public functions (`prepare_tissue`, `simulate_*`,
  `JevProvider`).
- **Viewer** (`web/melanoma.html`, `web/tissue.html`): title, metrics and colors are
  melanoma-specific. Read labels and palette from the guidebook.
- **CI** runs Ubuntu and Windows but no macOS, although the README says local checks were
  on macOS; add `macos-latest` or remove the claim.
- **Endpoint** `ENDPOINT` is a module constant (`jev_provider.py:21`); allow an
  environment override (`JEV_ENDPOINT`) for staging or self-hosted gateways, still with
  no redirects.
- **Warnings**: 215 test warnings, mostly from anndata on pandas 3; filter or fix so that
  real warnings stay visible.

## Verification performed for this review

- `pip install -e ".[test]"` and `python -m pytest -q` on Python 3.13: 60 passed, 1 skipped.
- A live Jev smoke test (one question, read from `TYPESAFE_API_KEY`) could not be run from
  the review environment because outbound access to `api.typesafe.ai` was denied by its
  network policy. The provider code path up to the HTTP call was exercised. This test should
  be repeated by task T0.3 in [PLAN.md](PLAN.md).
