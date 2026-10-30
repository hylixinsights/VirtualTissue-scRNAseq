# Plan: from an annotated single-cell object to a virtual tissue driven by Jev

Status: proposal, 2026-10-30. Companion: [CODE_REVIEW.md](CODE_REVIEW.md).

## 1. Goal in one paragraph

A user brings a dissociated single-cell (or single-nucleus) RNA-seq object with
annotated cell types, optionally with trajectory and CellChat results, and declares
the **context**: tissue of origin, disease, perturbation (treatment, mutation, lesion),
and whether the library is enriched or represents the whole tissue. The pipeline
learns a **guidebook**: for each cell type (and state), which actions it can and
cannot perform, with context-conditioned prior probabilities and the partners it
interacts with. It then builds a **virtual tissue** whose cell frequencies match a
chosen sample, and lets every cell ask **Jev** what to do next, given its
guidebook entry and its neighborhood. A perturbation is a change of the context
"prompt" (for example *non-lesional skin* → *psoriatic lesion*), and we watch the
cells' decisions and the tissue's drift.

First use case: psoriasis skin (lesional vs non-lesional). Later: melanoma, nasal
mucosa, lymph node, brain (single-nucleus).

## 2. Architecture

```
 user object (.h5ad / Seurat .rds)        context.yaml (declared)
 + optional trajectory, CellChat                 │
            │                                    │
            ▼                                    ▼
   [1] ingest & validate  ──────────────►  [2] context resolution
            │                                    │
            ▼                                    │
   [3] cell-type mapping (canonical vocabulary)  │
            │                                    │
            ▼                                    ▼
   [4] guidebook builder  ◄──────────────────────┘
       per cell type × state × condition:
       can / cannot, action priors, partners, fields, transitions
            │                 guidebook.json + guidebook.md (human review)
            ▼
   [5] virtual tissue initializer (frequencies of a chosen sample, layout template)
            │
            ▼
   [6] simulation kernel ── each decision checkpoint ──► Jev (context + guidebook + neighborhood)
            │                  ◄── probabilities over permitted actions
            │  kernel samples action, applies physics, fates, fields
            ▼
   [7] outputs: journals, trajectories, composition over time, viewer
            │
            ▼
   [8] evaluation: simulated drift vs observed perturbed sample
```

Design rules carried over from v0.5 (see CODE_REVIEW.md, "What works well"):
Jev ranks only kernel-permitted actions; the kernel owns death, division and physics;
every call is journaled and replayable; the key is read only from `TYPESAFE_API_KEY`.

## 3. The context declaration

Every run starts from a `context.yaml`. Example for this project:
[`examples/contexts/psoriasis_skin.yaml`](../examples/contexts/psoriasis_skin.yaml).

| Field | Meaning | Example |
|---|---|---|
| `organism` | Species | `human` |
| `tissue` | Tissue of origin, free text plus optional UBERON ID | `skin (trunk/limb), UBERON:0002097` |
| `assay` | `scRNA-seq` or `snRNA-seq`, chemistry | `scRNA-seq, 10x 3' v3` |
| `library` | `whole_tissue` or `enriched` (+ how: sort markers, depletion) | `whole_tissue` |
| `dissociation_bias` | Known under/over-represented types | `neutrophils under-represented` |
| `conditions` | Named contexts, each mapped to `obs` values | `non_lesional`, `lesional`, `healthy` |
| `baseline_condition` | Context used to initialize the tissue | `non_lesional` |
| `perturbations` | Named prompts: what changes, mechanism text, optional schedule | `psoriatic_lesion` |
| `disease` | Disease name, optional MONDO ID | `psoriasis, MONDO:0005083` |
| `layout_template` | Geometry template | `skin_layered` |

If `library: enriched`, the guidebook still learns per-type behavior, but tissue
frequencies must come from a declared correction (`frequency_source:
reference_table.tsv` or `user_declared`), never from the enriched counts.

## 4. The guidebook

`guidebook.json` (schema-validated, versioned) plus a human-readable
`guidebook.md`. One entry per canonical cell type, with optional states:

```yaml
Keratinocyte_basal:
  ontology: CL:0002187
  mobile: false
  can_divide: true
  radius_um: 6.0
  states: [quiescent, cycling]                 # from clustering or trajectory
  transitions:                                 # from trajectory, else assumed
    - {to: Keratinocyte_spinous, prior_rate_per_h: 0.01, evidence: inferred_pseudotime}
  actions:
    proliferate:  {allowed: true,  prior: {non_lesional: 0.05, lesional: 0.20}, evidence: "MKI67/TOP2A program score"}
    differentiate:{allowed: true,  prior: {non_lesional: 0.10, lesional: 0.06}, evidence: "pseudotime toward spinous"}
    secrete:      {allowed: true,  signals: [IL36G, S100A8, CXCL1], evidence: "CellChat lesional"}
    rest:         {allowed: true,  prior: {non_lesional: 0.85, lesional: 0.74}}
    migrate:      {allowed: false, reason: "anchored to basement membrane"}
  partners:                                    # from CellChat (prob, pathway)
    - {cell_type: T_CD8_Tc17, ligand: IL17A, receptor: IL17RA, condition: lesional, weight: 0.42}
  markers_used: [KRT5, KRT14, KRT15]
```

How each part is derived (all values carry an `evidence` tag: `observed`,
`inferred`, `assumed`, or `literature`):

| Guidebook part | Source | Method |
|---|---|---|
| Frequencies per condition | `obs` counts per `sample_id` | Proportions per sample; condition mean with donor-level CI |
| Program scores | Expression | Gene-set scores (proliferation, cytokine, differentiation, stress, death) per type × condition |
| Action priors | Program scores, contrasts | Monotone map from score to prior, normalized over allowed actions; lesional vs non-lesional shift from paired donors |
| Can / cannot | Canonical vocabulary + defaults per lineage | Rule table (e.g. neurons: no division, no migration), user-overridable |
| Transitions | Trajectory (pseudotime, CellRank) | Allowed edges and relative speed; absent trajectory → `assumed` |
| Partners and signals | CellChat / LIANA / CellPhoneDB tables | Top ligand-receptor pairs per sender-receiver per condition |
| Fields | Secreted ligands from CellChat | Named diffusing fields (IL17A, IL23, IL36G, TNF, IFNG, CXCL12...) replacing the single generic cue |

The guidebook is reviewed by a human before simulation (task C3). It is the main
scientific artifact and can be shared on its own.

## 5. How Jev is used

At each decision checkpoint the kernel builds, per cell, the menu of legal actions
(from the guidebook plus local physics) and sends:

- the **run context** (tissue, disease, active perturbation prompt);
- the cell's **guidebook entry** for its type and state, including priors for the active condition;
- the **local neighborhood** (counts of neighbor types, contacts, local field levels);
- the **simulated state** (fatigue, damage, cycle, time in state).

Jev returns a probability over the permitted actions. The kernel **samples** from
`(1-w)·p_jev + w·p_prior` with a keyed RNG (`w` is configurable; `w = 1` is the free
"guidebook only" mode) and records prior, Jev answer, mixture and draw. Perturbation
= changing the context block from `non_lesional` to `lesional` (the "prompt"), while
the initial tissue stays the same.

## 6. First dataset: psoriasis skin, lesional vs non-lesional

**Primary: GEO [GSE173706](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173706)**
— Ma F, Plazyo O, Billi AC, ... Gudjonsson JE. *Single cell and spatial sequencing
define processes by which keratinocytes and fibroblasts amplify inflammatory
responses in psoriasis.* Nat Commun 14, 3455 (2023).
[doi:10.1038/s41467-023-39020-4](https://doi.org/10.1038/s41467-023-39020-4).

- 33 samples: 14 psoriatic lesional (`PP`), 11 non-lesional skin from the same patients (`PN`),
  8 healthy donors (`NS`). Paired lesional / non-lesional donors are ideal for the
  "lesion as a prompt" design.
- Illumina NovaSeq 6000; BioProject PRJNA726808; SRA SRP318151.
- Supplementary file: `GSE173706_RAW.tar` (about 252 MB, CSV count matrices).
  The series combines single-cell and spatial samples; the sample list must be
  filtered to the scRNA-seq samples (task A1).
- Per-cell annotations are **not** on GEO. Re-annotate with the published marker
  sets (keratinocyte layers, SFRP2+ fibroblasts, CD8 Tc17, LAMP3+ DC, etc.) and
  check counts against the paper's figures; or request the annotated object from
  the authors.

Download:

```bash
mkdir -p data/raw/GSE173706 && cd data/raw/GSE173706
curl -L -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173706/suppl/GSE173706_RAW.tar"
curl -L -O "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE173nnn/GSE173706/matrix/GSE173706_series_matrix.txt.gz"
tar -xf GSE173706_RAW.tar
sha256sum * > SHA256SUMS   # record, then commit only the checksums, never the data
```

**Secondary (annotated, cross-check): Reynolds G et al.** *Developmental cell programs are
co-opted in inflammatory skin disease.* Science 371, eaba6500 (2021).
[doi:10.1126/science.aba6500](https://doi.org/10.1126/science.aba6500). Healthy adult skin,
atopic dermatitis and psoriasis, lesional and non-lesional, with author cell
annotations; raw data in ArrayExpress E-MTAB-8142 and processed objects distributed
through the Human Cell Atlas skin portal / CELLxGENE. The exact download link for
the annotated `.h5ad` was not verified in this review (task A1b).

**Not suitable as primary:** GSE151177 (Kim J, Krueger JG; emigrating cells after 48 h
culture, lesional vs healthy only, no non-lesional group).

## 7. Task board

Legend: **SEQ** = must wait for the listed dependencies; **PAR** = can run in
parallel with other tasks in the same phase. **Review gate** = code review (PR
reviewed and merged) before dependants start. Owners are suggestions.

### Phase 0. Foundations (sequential, do first)

| ID | Task | Depends on | Mode | Output | Done when |
|---|---|---|---|---|---|
| T0.1 | Merge this plan and review | none | SEQ | `docs/PLAN.md`, `docs/CODE_REVIEW.md` | PR merged |
| T0.2 | Restructure repo: move SPECIAL/Lin audit to `research/`, lock files to `env/`, loosen pins (review M2, M3) | T0.1 | SEQ | Clean root; `pip install -e .` works in a scanpy env | CI green; review gate |
| T0.3 | Jev smoke test from a network that allows `api.typesafe.ai`: one question, key from `TYPESAFE_API_KEY`, journal saved | none | PAR with T0.2 | `reports/jev-smoke/` (journal, no key) | Valid Choice response parsed |

### Phase 1. Inputs and engine generalization (parallel tracks)

| ID | Task | Depends on | Mode | Output |
|---|---|---|---|---|
| A1 | Download GSE173706, select scRNA-seq samples, build one AnnData with raw counts and `obs`: `sample_id`, `donor_id`, `condition` | T0.1 | PAR | `data/psoriasis/GSE173706_raw.h5ad` (not in Git) + checksums |
| A1b | Locate and download Reynolds 2021 annotated object; record link and checksum | T0.1 | PAR | Download script + `SOURCE.json` |
| B1 | Context schema (`context.schema.json`), loader, `vtsc context-template` | T0.2 | PAR | Schema, tests, example YAMLs |
| B2 | Canonical cell-type vocabulary + mapping tool (`vtsc map-celltypes`), with Cell Ontology IDs | T0.2 | PAR | `resources/celltypes_v1.json`, mapping TSV format |
| B3 | Action registry + single `legal_actions` (review H1, H3); port melanoma as scenario pack; bit-identical regression test vs v0.5 | T0.2 | PAR, review gate | `vtsc/actions/`, `scenarios/melanoma_anti_pd1/` |
| B4 | Named diffusing fields replacing the generic cue; KD-tree spatial queries (review M4) | B3 | SEQ | `vtsc/fields.py` |
| D1 | Jev payload v2: context block, guidebook entry, neighborhood; `sampling: argmax | sample`, prior mixing, `guidebook-only` mode (review H2) | T0.2 | PAR, review gate | `PROMPT_VERSION = vtsc-behavior-2` |
| R1 | R helper `r/export_for_vtsc.R`: Seurat → H5AD with counts, metadata, reductions; CellChat `subsetCommunication` → TSV per condition | T0.2 | PAR | Script + example |

### Phase 2. Psoriasis analysis (after A1)

| ID | Task | Depends on | Mode | Output |
|---|---|---|---|---|
| A2 | QC, integration, annotation of GSE173706 with published markers; compare type counts with the paper | A1 | SEQ | Annotated `.h5ad`, `annotation_report.md` |
| A3 | Trajectory: keratinocyte basal → spinous → granular; fibroblast → inflammatory SFRP2+ state (e.g. scanpy DPT / Slingshot / CellRank) | A2 | PAR | `obs.pseudotime`, transitions TSV |
| A4 | CellChat (or LIANA) separately for lesional, non-lesional, healthy | A2 | PAR | 3 communication TSVs |
| A5 | Cross-check annotations against Reynolds 2021 (label transfer) | A2, A1b | PAR | Agreement table |

### Phase 3. Guidebook (sequential core)

| ID | Task | Depends on | Mode | Output |
|---|---|---|---|---|
| C1 | Guidebook schema (`guidebook.schema.json`) and markdown renderer | B1, B2 | SEQ | Schema, renderer, tests |
| C2 | Guidebook builder `vtsc build-guidebook object.h5ad --context context.yaml [--trajectory t.tsv] [--cellchat dir/]` | C1, A3, A4 | SEQ | `guidebook.json`, `guidebook.md` |
| C3 | Human review of the psoriasis guidebook (dermatology / immunology reader) | C2 | SEQ, review gate | Signed-off guidebook v1 |

### Phase 4. Virtual tissue and simulation

| ID | Task | Depends on | Mode | Output |
|---|---|---|---|---|
| E1 | Initializer `vtsc init-tissue`: sample N cells matching frequencies of a chosen sample/condition; layout templates `skin_layered`, `mixed`, `tumor_nests` | B3, C1 | PAR with C2 | `tissue_init.zip` |
| E2 | Generic run command `vtsc simulate --guidebook g.json --tissue t.zip --context-prompt lesional --decision-mode guidebook-only|mock-jev|jev` | B4, D1, E1, C3 | SEQ | Run directory, journals |
| E3 | Generalize viewer labels, palettes and metrics from the guidebook | E2 | PAR | `explorer.html` |

### Phase 5. Psoriasis experiments and evaluation

| ID | Task | Depends on | Mode | Output |
|---|---|---|---|---|
| X1 | Free runs (`guidebook-only`): non-lesional tissue under `non_lesional` vs `lesional` prompt; 10 seeds | E2 | SEQ | Paired runs |
| X2 | Same arms with `mock-jev` (pipeline test), then live `jev` with a request ceiling | X1, T0.3 | SEQ | Paired runs + usage |
| X3 | Controls: frozen, shuffled context, prompt-off, interaction-off | E2 | PAR with X2 | Control runs |
| X4 | Evaluation: does the simulated drift (composition, keratinocyte proliferation, T-cell expansion, inflammatory fibroblast fraction) move toward the observed lesional sample? Leave-donors-out | X1–X3 | SEQ | `reports/psoriasis/EVALUATION.md` |

### Phase 6. Generalization (each tissue is a parallel track after X4)

| ID | Tissue | Specific issues |
|---|---|---|
| G1 | Melanoma | Port v0.5 anti-PD-1 scenario onto the generic engine (already regression-tested in B3) |
| G2 | Nasal mucosa | Epithelial layer template; ciliated / goblet / basal cells; airway perturbations (allergen, virus) |
| G3 | Lymph node | Follicle / T-zone template; enriched libraries common (declare `library: enriched`) |
| G4 | Brain (snRNA-seq) | `assay: snRNA-seq`; neurons non-dividing and immobile; glia actions; nuclear transcript bias noted in evidence |

## 8. What to run, in order (operator view)

```
Step 1  T0.1  merge plan                         (now)
Step 2  T0.2 + T0.3                              (in parallel)
Step 3  A1, A1b, B1, B2, B3, D1, R1              (all in parallel)
Step 4  A2 → (A3 ‖ A4 ‖ A5)       B4             (parallel tracks)
Step 5  C1 → C2 → C3                             (sequential)   ‖ E1, E3 prep
Step 6  E2                                       (after C3, B4, D1, E1)
Step 7  X1 → X2 ‖ X3 → X4                        
Step 8  G1 ‖ G2 ‖ G3 ‖ G4
```

Every code task ends in a pull request with tests and a review before merge. Data
never enters Git; only download scripts, checksums and `SOURCE.json` records do.

## 9. Proposed user-facing commands (target for v0.6)

```bash
vtsc context-template my_context.yaml
vtsc ingest my_object.h5ad outputs/ingested.zip --context my_context.yaml \
     --cell-type-column cell_type --sample-column sample_id --condition-column condition
vtsc build-guidebook outputs/ingested.zip outputs/guidebook \
     --trajectory pseudotime.tsv --cellchat cellchat_tables/
vtsc init-tissue outputs/ingested.zip outputs/tissue.zip --from-condition non_lesional --cells 800
vtsc simulate outputs/tissue.zip outputs/run_lesional --guidebook outputs/guidebook/guidebook.json \
     --context-prompt lesional --decision-mode guidebook-only --hours 48 --seed 1
vtsc compare outputs/run_nonlesional outputs/run_lesional outputs/comparison
vtsc explore outputs/comparison --open
```

Seurat users first run `Rscript r/export_for_vtsc.R my_seurat.rds my_object.h5ad`.

## 10. Risks and open questions

- **No ground-truth dynamics.** Dissociated data are snapshots; transitions and rates
  remain priors. Evaluation (X4) tests only whether the drift points toward the observed
  perturbed state, not kinetics.
- **Geometry is assumed.** Layout templates encode tissue architecture by hand; spatial
  data (the GSE173706 spatial samples, or the VirtualTissue-Xenium project) can later
  calibrate them.
- **Jev variability and cost.** Sampling from Jev probabilities adds variance; use seeds,
  ceilings and journals. The guidebook-only mode is the free baseline.
- **Annotation quality** drives everything; C3 is a hard gate.
- **Network access.** Jev calls need outbound access to `api.typesafe.ai`.
