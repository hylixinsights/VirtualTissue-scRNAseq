# VirtualTissue-SC × SPECIAL: predicting CD8 T-cell state dynamics under anti-PD-1, then adding space and agent decisions

**Execution plan and supplementary-material blueprint**

| Field | Value |
|---|---|
| Version | 1.0 (plan only; nothing executed) |
| Date | 2026-10-27 |
| Owner | Helder Nakaya (CSBL) |
| Supersedes | `claude/virtualtissue-sc-trajectory-plan-v0-2026-10-27.md` (exploratory v0) |
| Related | `claude/virtualtissue-xenium-0.1.1-review-2026-10-26.md`; VirtualTissue-Xenium v0.1.2 (commit `5482cbd`) |
| Language | English (code, docs, outputs) |

---

## How to use this document in a new conversation

1. Read §0–§4 once. They define the question, the boundaries of what can be claimed, and the data.
2. Open the **Task board (§11)**. Pick the first unticked task whose dependencies are complete.
3. For every task:
   - follow its steps;
   - write outputs to the stated paths;
   - check its acceptance criteria;
   - tick the box;
   - add an entry to the **Decision log** or **Deviation log** (§13) whenever a choice departs from this plan.
4. Save each task report to the Project as `claude/vtsc-<TASK-ID>-<slug>.md`, then update the task status here (re-upload this file with the ticks).
5. **Hard rules for every session:**
   - Never inspect or compute state changes for **lockbox** patients before task T4.1 (§7.1).
   - Never use the future response label, or on-treatment data of a test patient, as model input.
   - Spend no Jev credits without explicit approval and a request cap.
   - Do not edit the validated VirtualTissue-Xenium engine. New biology goes into a new, versioned runtime (Phase 5).
   - Every number placed in a report must come from a saved output file with a recorded hash.

**Size key:** S = under half a day; M = 1–2 days; L = 3 or more days of focused work.

---

## 0. Summary

We will test whether the change in CD8 T-cell state composition between a pre-treatment and an early on-treatment melanoma biopsy can be **predicted in patients not used for training**, from the pre-treatment biopsy alone (CD8 states plus tumour-microenvironment features). Data come from the SPECIAL cohort (Poźniak, Marine and colleagues; KU Leuven).

The project proceeds as a ladder. Each rung must beat the simpler rungs on held-out patients before the next is built:

1. A **feasibility audit** of the public SPECIAL objects.
2. A **CD8 state reference** (states, programs, transition graph) that does not impose a single naive → exhausted sequence.
3. A **non-spatial, few-parameter dynamics model** (transition, net proliferation, influx), compared with persistence and mean-shift baselines under patient-level cross-validation and a lockbox.
4. **External validation** and complementary datasets.
5. **Coupling to VirtualTissue** (spatial agents), kept only if space improves prediction or explains a pre-specified spatial pattern.
6. **Jev as an alternative transition policy**, kept only if it matches or improves on estimated probabilities in held-out patients.

Treatment contexts beyond anti-PD-1 monotherapy (combination therapy, other sites, healthy skin, untreated counterfactuals) are deferred.

---

## 1. Question, hypotheses and estimand

### 1.1 Primary question

> From the pre-treatment biopsy and its microenvironment, can we predict the change in the distribution of CD8 T-cell states after the first cycle of anti-PD-1 monotherapy in patients who were not used for training?

### 1.2 Estimand

For each eligible patient *i* on nivolumab monotherapy, the estimand is the **CD8 state composition of the early on-treatment biopsy**, $\pi_i^{on}$, predicted from information available at baseline: $\pi_i^{pre}$, environment features $z_i$, and the inter-biopsy interval $\Delta t_i$.

- Unit of analysis: **the patient**.
- Population: SPECIAL anti-PD-1 monotherapy patients with paired, eligible biopsies.

### 1.3 Hypotheses (pre-specified, tested in fixed sequence)

| ID | Hypothesis | Test |
|---|---|---|
| H1 | A state-transition model (M2) predicts $\pi^{on}$ better than persistence (M0) and the cohort mean shift (M1) | Paired, one-sided comparison of per-patient Aitchison error in held-out patients |
| H2 | Adding pre-specified baseline environment features (M3) improves on M2 | Same test, M3 vs M2; only if H1 holds (fixed-sequence testing controls the family-wise error) |
| H3 (secondary) | A spatial agent-based version (Phase 5) is non-inferior to the best non-spatial model and reproduces at least one pre-specified spatial pattern that the non-spatial model cannot express | Non-inferiority margin fixed in T5.2 |
| H4 (secondary) | Jev as transition policy is non-inferior to estimated transition probabilities | Same metrics; prompt and label controls |

**Exploratory analyses** (no confirmatory claims):
- whether prediction error differs between later responders and non-responders, with response used for stratification only, never as input;
- combination-therapy patients;
- mechanism attribution (transition vs proliferation vs influx).

### 1.4 What this project does **not** claim

- **No untreated counterfactual.** A pre-treatment biopsy is **not** a longitudinal untreated control. Any simulated "no treatment" arm is a hypothesis-dependent counterfactual and is labelled as such. Non-responders are not used as surrogates for untreated patients.
- **No response prediction.** Clinical response is not a model output.
- **No mother–daughter lineage.** Shared TCR clonotypes, where available, are clonal evidence, not a mother–daughter tree. Pseudotime is not time.
- **No assumed equilibrium.** A single biopsy does not demonstrate a stable state distribution. Models are not fitted to reproduce stability.
- **No absolute dynamics.** With proportions only, expansion of one state and loss of others can be indistinguishable (§6.3).

---

## 2. Design principles and claim boundaries

1. **Simple before complex.** Each added layer (environment, space, Jev) must earn its place on held-out patients.
2. **The patient is the unit** for splitting, validation and inference. Cells from one patient are never split between training and test.
3. **Separate the functional dimensions.** Dysfunction, proliferation, persistence and cytotoxicity are distinct axes. Dysfunctional CD8 cells in human melanoma can be highly proliferative (Li et al., 2019). No state is assumed to "do little of everything".
4. **Start from observed precursor states.** We make no assumption that every TCF7+ cell is a progenitor-exhausted (Tpex) cell, or that priming happens only in lymph nodes. Intratumoural naive T-cell activation has been shown experimentally (Thompson et al., 2010).
5. **The observed manifold constrains states, not histories.** Avatars, if used, represent a state; the sequence of jumps is simulated. An agent's patient, identity and clone are kept separate from the reference profile used to represent its state.
6. **Pre-specify, then log deviations.** This plan is frozen by hash at T0.4. Every later change goes to the Deviation log.
7. **Traceability.** Every reported number links to a task output with a SHA-256 hash.

### 2.1 Evidence class of every quantity

| Quantity | Class | Source |
|---|---|---|
| Cell counts per sample, gene counts | Measured (dissociation- and chemistry-biased) | SPECIAL objects |
| Cell type and state labels | Inferred (author annotation, audited) | SPECIAL objects; T2.1 |
| Program scores (AUCell or alternatives) | Inferred from counts | T2.3 |
| State graph, pseudotime | Inferred; robustness-gated | T2.4 |
| Inter-biopsy interval | Measured if per-patient dates exist; otherwise nominal | Clinical table; T1.2 |
| Transition, net-growth and influx rates | Estimated (model-dependent; partially identifiable) | Phase 3 |
| Environment features $z$ | Inferred from the pre-treatment biopsy | T2.6 |
| Geometry, contacts (Phase 5) | Assigned (synthetic) or measured (Xenium, other samples or patients) | Phase 5 |
| Agent choices (Phase 6) | Provider output (Jev) | Phase 6 |
| Response labels | Clinical; used for stratification and exploratory analysis only | Clinical table |

---

## 3. Data sources

### 3.1 SPECIAL cohort: facts confirmed from public sources

| Item | Value | Source |
|---|---|---|
| Study | SPECIAL trial (UZ/KU Leuven S62275): "Dissecting the melanoma ecosystem one cell at the time during immunotherapy" | EGA study EGAS00001006488; bioRxiv 2023.12.14.571631 |
| Patients / biopsies | 23 treatment-naive patients, stage IIIB–IV melanoma; 46 biopsies | EGA EGAD00001009291; bioRxiv |
| Pairing | Matched pre- and on-treatment biopsies in 20 of 23 patients | EGA; bioRxiv |
| On-treatment timing | Early, before the second ICB cycle (preprint: "within ±2 weeks"). Per-patient interval to be confirmed (T1.2) | bioRxiv |
| Treatment arms | Nivolumab monotherapy, n = 17; ipilimumab + nivolumab, n = 6 | bioRxiv |
| Biopsy sites | 12 cutaneous, 8 subcutaneous, 26 lymph node | bioRxiv |
| scRNA-seq chemistry | 10x Chromium: 16 libraries Single Cell 3′ v3; the rest 5′. Illumina NextSeq / HiSeq 4000 / NovaSeq 6000; GRCh38 | EGA dataset description |
| TCR / VDJ | Not mentioned in the EGA dataset description or the preprint. **Treat as unavailable** until shown otherwise | EGA; bioRxiv |
| Immune cells after QC | 25,401 (preprint; the deposited object may differ) | bioRxiv |
| CD8 states (preprint) | Memory, Activated, Cytotoxic, Dysfunctional | bioRxiv |
| Raw data | Controlled access, EGA (VIB-UZL DAC, EGAC00001002806) | EGA |
| Processed objects | Public, KU Leuven RDR, DOI 10.48804/GSAXBN; licence **CC BY-NC-SA 4.0** | researchportal.be record |
| Clinical metadata | Table S1 of Poźniak et al., *Cell* 2024 (per the deposit description); the 2025 NK paper also has supplements. Reconcile both (T1.2) | Deposit description |
| Related publications | Poźniak J et al., *Cell* 2024 (TCF4 network; PMID 38181739). Roda N et al., 2025, "Cytotoxic NK cells impede response to checkpoint immunotherapy in melanoma with an immune-excluded phenotype" (PMID 40530504; PMC12409281) | PubMed / PMC |

**Facts with direct consequences for the design:**
- **26 of 46 biopsies are lymph node metastases.** Biopsy site is a major covariate (T1.6).
- **Chemistry is mixed (3′ v3 vs 5′).** Chemistry must be audited for confounding with timepoint and patient (T1.6).
- **Pre- and on-treatment biopsies may come from different lesions or sites.** This must be audited per patient (T1.6).
- **Published NK findings come from the same cohort.** Cytotoxic NK cells were enriched in on-treatment non-responder lesions and excluded from tumour parenchyma. Any NK-based environment feature is therefore informed by prior analysis of these data (circularity flag; exploratory only).

### 3.2 Files to obtain (Phase 1)

Sizes are as reported on the RDR page (dataset version 8). Verify them at download (T0.2).

| Priority | File | Reported size | Use |
|---|---|---:|---|
| 1 | `GC_all_immune_sharing.rds` | 945.2 MB | Seurat object, annotated immune cells with response information |
| 1 | Clinical Table S1 (Poźniak 2024 *Cell*; plus the 2025 supplements) | — | Patient, timepoint, treatment, site, response, dates |
| 2 | `Entire_TME.rds` | 2.5 GB | Microenvironment features (T2.6) |
| 2 | `Malignant_cells.rds` | 824.4 MB | Malignant states (column `Malignant_clusters`) for environment features |
| 3 (Phase 5) | Xenium bundles (e.g., `GC19_GC16.zip`, 6.6 GB) | 6.6 GB and larger | Spatial validation. Panel: skin panel plus custom add-on, to verify |

**Do not download the full ~69 GB deposit.** The total for priorities 1–2 is about 4.3 GB.

### 3.3 External and complementary datasets (Phase 4)

| Dataset | Role | Caveat |
|---|---|---|
| GSE120575 (Sade-Feldman et al., 2018) | External check of the **direction** of CD8 state change pre vs post ICB in melanoma | Smart-seq2, CD45+ sorted: no tumour or stroma; states must be mapped by label transfer |
| GSE272993 (melanoma PBMC; RNA + TCR + surface protein; baseline and every ~3 weeks, as reported) | Temporal and clonal dynamics of circulating CD8; plausibility of influx and timescale priors | Blood, not tumour; details to verify on GEO |
| GSE123813 (Yost et al., 2019) | 10x + TCR pre/post anti-PD-1; clonal replacement; tests of the influx component | Basal and squamous cell carcinoma, not melanoma |

### 3.4 Licence, ethics and data use

- **CC BY-NC-SA 4.0.** Derived data (exported matrices, bundles, avatars, simulation inputs) must carry the same licence, keep attribution and **must not be used commercially**.
  - Keep this work separate from Hylix commercial services.
  - The MIT licence of VirtualTissue code is unaffected, but **released derived data** cannot be MIT.
- No raw (EGA) data are required. No attempt at re-identification. Patient IDs remain the study pseudonyms.
- Draft attribution: *"Processed single-cell data from the SPECIAL cohort (Poźniak et al., Cell 2024; Roda et al., 2025), KU Leuven RDR, doi:10.48804/GSAXBN, CC BY-NC-SA 4.0."*

### 3.5 Data logistics

- Store raw downloads outside any repository, e.g. `<DATA_ROOT>/SPECIAL/raw/`, and derived files in `<DATA_ROOT>/SPECIAL/derived/`.
- The `.rds` files need **R with a Seurat/SeuratObject version able to read them**. Record the versions.
- Convert once to a neutral format (T1.1): per-sample MatrixMarket counts, `obs.parquet`, `var.tsv`, embeddings. Then work in Python.
- Large single files may not transfer between environments (e.g., a ~400 MB per-file limit into a cloud session). Do the conversion where the `.rds` lives (local machine or Colab), then move the smaller per-sample files.
- Expect several GB of RAM to load `Entire_TME.rds`. Use a local workstation or high-RAM Colab.

---

## 4. Cohort definition and eligibility

### 4.1 Primary analysis set (pre-specified)

A patient is **eligible** if all of the following hold:
1. Nivolumab monotherapy.
2. Both pre- and on-treatment biopsies are present in the immune object.
3. Each biopsy has at least **30 CD8 T cells** passing QC. Sensitivity thresholds: 15 and 50.
4. The interval between biopsies is known or can be bounded.

Non-exclusion **flags**, used for sensitivity analyses:
- different lesion or anatomical site between timepoints;
- different chemistry between timepoints;
- lymph-node vs non-lymph-node biopsy.

### 4.2 Secondary sets

- **Combination-therapy pairs** (n ≤ 6): a separate exploratory analysis, never pooled with monotherapy.
- **Unpaired biopsies:** may inform the state reference (T2.x) but never the dynamics evaluation.

### 4.3 Flow diagram

A CONSORT-style diagram (Supplementary Figure S1) runs: patients → biopsies → cells → CD8 → eligible pairs → development / lockbox. It is produced in T1.3.

---

## 5. CD8 state reference

### 5.1 Starting point

- Use the authors' annotation as baseline (preprint: Memory, Activated, Cytotoxic, Dysfunctional). Audit it (T2.1) and do not replace it unless a documented failure is found.
- **Primary state granularity K = 4** (the published states). Finer substates (e.g., a precursor-like or cycling substate) are a sensitivity analysis.
- **Cycling** is treated as a **functional axis** (proliferation), not necessarily a separate state.

### 5.2 Program library (T2.2)

Each program has: name, genes, source (database or publication, with version), category, and coverage in the data. The library is frozen as a versioned JSON (Supplementary Table S4).

| Category | Programs (example genes; to be finalised from cited sources) |
|---|---|
| State identity | naive/stem-like memory (TCF7, SELL, CCR7, LEF1, IL7R); effector memory (GZMK, CXCR3); cytotoxic effector (GZMB, PRF1, GNLY, NKG7, FGFBP2, CX3CR1); dysfunctional (PDCD1, HAVCR2, LAG3, TIGIT, ENTPD1, TOX, CXCL13, LAYN); tissue residency (ITGAE, ZNF683, CXCR6) |
| Function (separate axes) | proliferation (G2M and S lists; MKI67, TOP2A); cytotoxic capacity; effector cytokines (IFNG, TNF); persistence/survival (BCL2, IL7R) |
| Response / exposure | recent TCR signalling (NR4A1–3, EGR2); IFN-γ response (Hallmark); TGF-β response; hypoxia (Hallmark); IL-2/STAT5 |
| Technical controls | dissociation-induced immediate-early genes (FOS, JUN, EGR1, HSP genes); mitochondrial fraction; ribosomal fraction |

**Overlap policy** (replaces the "fully disjoint" rule of v0):
- Report the Jaccard overlap between every pair of programs.
- Use disjoint versions only as a **sensitivity control**, since full disjointness removes shared biology and does not remove confounding.
- Validate state assignments with **markers not used to define the states**.
- Keep test patients out of every fitted step.

**Known confounder:** NR4A and EGR genes respond both to TCR signalling and to dissociation stress. Always report the TCR-signalling program jointly with the dissociation program.

### 5.3 Scoring and calibration (T2.3)

- Default method: AUCell (decoupler; pinned version).
  - `n_up` is chosen from the distribution of detected genes per CD8 cell. Cells with fewer detected genes than `n_up` are excluded and counted.
  - Alternatives for sensitivity: UCell, ULM on raw programs.
- Percentile calibration:
  - **within the CD8 lineage** (not within state), conditioned on detection-depth decile;
  - reference distribution fitted **only on training patients of each fold** and then applied to the held-out patient (fold-aware calibration).
- Chemistry check: a mixed model of score ~ chemistry + state + (1 | patient). If the chemistry effect exceeds a pre-set fraction (e.g., 25%) of the between-state difference, score-based analyses become chemistry-stratified or chemistry-adjusted.

### 5.4 State graph and pseudotime (T2.4)

- **Allowed transitions** (edge set *E*):
  - from PAGA connectivity on CD8 cells of development patients;
  - repeated over patient bootstrap resamples (≥ 200);
  - an edge is retained if present in ≥ 80% of resamples.
- **Direction:** pseudotime (diffusion pseudotime, Palantir, Slingshot) is computed with alternative roots and representations (HVG-PCA vs program space, with and without cycling genes). Direction is **"defensible"** only if:
  - orderings agree across methods (Spearman ρ ≥ 0.8);
  - the order is stable under patient bootstrap;
  - marker trends are monotonic as expected.
- **If direction is not defensible,** *E* is used **undirected**: transitions are allowed both ways and the direction is estimated from paired pre/on data in Phase 3. The paired design informs direction in a way a single snapshot cannot.
- Outputs: graph with edge confidence (Supplementary Figure S4); robustness table (Supplementary Table S5).

### 5.5 Functional axes (T2.5)

- For each state, compute the distributions of four axes: dysfunction, proliferation (cycling fraction), cytotoxicity and persistence.
- These later define capabilities in Phase 5 and priors on net growth in Phase 3.
- No axis is derived from another.

### 5.6 Environment features (T2.6)

**Pre-specified primary set** (at most 3 features, given n ≈ 15):

| Feature | Definition | Source object |
|---|---|---|
| $z_1$ tumour content | Malignant fraction among all QC cells of the biopsy | `Entire_TME.rds` |
| $z_2$ inflammation | Mean IFN-γ-response program (or CXCL9/CXCL10 program) in myeloid cells | `Entire_TME.rds` |
| $z_3$ site | Lymph node vs non-lymph-node | Clinical table |

**Exploratory features:**
- cytotoxic NK fraction (circularity flag, §3.1);
- malignant state composition (Poźniak 2024 states);
- MHC-I program in malignant cells;
- CD8 fraction of immune cells;
- Treg fraction.

**Rules:**
- Standardise within each training fold.
- Match sample IDs across objects before any feature is computed.
- Never double-count cells present in more than one object.

---

## 6. Non-spatial dynamics model

### 6.1 Notation

- $K$ states; patient $i$; observed counts $n_i^{pre}, n_i^{on} \in \mathbb{N}^K$.
- Totals: $N_i^{pre}, N_i^{on}$.
- Interval: $\Delta t_i$ (days).
- Environment: $z_i \in \mathbb{R}^p$.
- Composition: $\pi = n/N$ (with uncertainty).
- $\mathrm{clr}(p) = \log p - \overline{\log p}$.
- Aitchison distance: $d_A(p,q) = \lVert \mathrm{clr}(p) - \mathrm{clr}(q) \rVert_2$, with pseudocount 0.5 per state.

### 6.2 Model ladder

| Model | Prediction of $\hat\pi_i^{on}$ | Parameters |
|---|---|---|
| **M0 persistence** | $\hat\pi_i^{on} = \pi_i^{pre}$ | none (plus dispersion) |
| **M1 mean shift** | $\mathrm{clr}(\hat\pi_i^{on}) = \mathrm{clr}(\pi_i^{pre}) + \bar\delta$, where $\bar\delta$ is the mean training CLR change. M1b scales $\bar\delta$ by $\Delta t_i / \overline{\Delta t}$ | $K-1$ |
| **M2 state-only dynamics** | Linear compartment model (below), shared parameters | $\lvert E \rvert$ + growth + influx |
| **M3 environment-conditioned** | M2 with $q_{jk,i} = q^0_{jk} \exp(\beta_{c(jk)}^\top z_i)$; β shared within edge classes $c$ (e.g., "toward dysfunction", "toward cytotoxic", "other") | M2 + $p \times$ (number of classes) |

**Compartment model (M2/M3).** Let $x_i(t)$ be unnormalised state abundances:

$$
\frac{dx_i}{dt} = A_i x_i + \lambda_i e, \qquad A_i = Q_i^\top + \operatorname{diag}(\rho), \qquad \pi_i(t) = \frac{x_i(t)}{\mathbf{1}^\top x_i(t)}
$$

where:
- $Q_i$ is a generator with $q_{jk} \ge 0$ only for $(j \to k) \in E$;
- $\rho_k$ is the net growth of state $k$ (proliferation minus loss);
- $\lambda_i$ is influx relative to the initial CD8 abundance, with entry distribution $e$ on the simplex.

Initial condition: $x_i(0) \propto \pi_i^{pre}$, with $\pi_i^{pre}$ treated as latent, $\pi_i^{pre} \sim \mathrm{Dirichlet}(n_i^{pre} + \tfrac{1}{2})$.

**Likelihood:** $n_i^{on} \sim \mathrm{DirichletMultinomial}(N_i^{on}, \alpha\, \pi_i(\Delta t_i))$, with a shared dispersion $\alpha$.

**Within-state heterogeneity (score distributions).** The predicted on-treatment program-score distribution is a mixture over states, weighted by $\hat\pi_i^{on}$, of the empirical state-specific distributions from **training** patients' on-treatment cells:
- variant P-on (primary);
- variant P-pre: training pre-treatment cells, as a sensitivity analysis.

Avatars (real training cells) are drawn from these distributions. They never come from the test patient.

### 6.3 Mechanism variants and identifiability (T3.2–T3.3)

With two timepoints and proportions only, different mechanisms can produce the same $\pi^{on}$. We therefore fit and compare five variants:

| Variant | Transitions $Q$ | Net growth $\rho$ | Influx $\lambda$ |
|---|---|---|---|
| V-T | ✔ | 0 | 0 |
| V-P | 0 | ✔ | 0 |
| V-I | 0 | 0 | ✔ |
| V-TP | ✔ | ✔ | 0 |
| V-TPI | ✔ | ✔ (prior from cycling fraction) | ✔ (prior from external data) |

**Identifiability procedure:**
- **Simulation-based parameter recovery:** simulate cohorts with the observed sample sizes, $\Delta t$ and dispersion under known parameters; fit; report bias, coverage and correlation of estimates.
- **Profile likelihood / posterior correlations** on the real data.
- Variants whose held-out predictions are indistinguishable are reported as an **equivalence class**, not as mechanistic evidence.
- **Priors:**
  - $\rho_k$ centred on the state's cycling fraction divided by a cycle-duration prior, minus a loss prior;
  - transition timescales weakly informative (days to weeks);
  - influx prior from complementary data (T4.3) or left wide.
  - Each prior gets a prior-sensitivity analysis.

### 6.4 Inference

- Bayesian (NumPyro or Stan) or penalised maximum likelihood. Pin the choice in T3.1.
- Matrix exponential for the linear system; ODE solver if influx makes it inhomogeneous.
- Diagnostics: R-hat, effective sample size, divergences, posterior predictive checks.

---

## 7. Evaluation protocol

### 7.1 Splits

1. After eligibility (T1.3) and **before any state-composition change is computed**, randomly assign a **lockbox** of 4 eligible monotherapy patients:
   - stratified by later response and by site pairing;
   - recorded seed (master seed `20261027`).
2. The **development set** is everything else. All model development happens here by **leave-one-patient-out (LOPO)** cross-validation:
   - state calibration, scoring references, feature standardisation and every fitted quantity are refit inside each fold;
   - hyperparameters are pre-specified or chosen by nested LOPO.
3. Freeze the pipeline (T3.6), then evaluate:
   - (a) on the lockbox: unbiased, small n;
   - (b) by full LOPO on all eligible patients with the frozen pipeline.
4. Report both. Label (b) as development-informed.

### 7.2 Metrics

| Rank | Metric | Notes |
|---|---|---|
| **Primary** | Per-patient Aitchison distance $d_A(\hat\pi_i^{on}, \pi_i^{on})$ | Pseudocount 0.5 |
| Key secondary | Held-out log predictive density of $n_i^{on}$ under the Dirichlet-multinomial | Proper scoring rule; uses uncertainty |
| Secondary | Per-state CLR error; sign accuracy of change per state | Interpretability |
| Secondary | Coverage of 80% predictive intervals per state | Calibration |
| Secondary | Energy distance between predicted and observed on-treatment program-score distributions | Within-state heterogeneity |
| Exploratory | Error stratified by later response; by site; by chemistry | No confirmatory claim |

### 7.3 Statistics

- Paired one-sided Wilcoxon signed-rank tests on per-patient errors, in fixed sequence: H1 (M2 vs M1, then M2 vs M0), then H2 (M3 vs M2). α = 0.05.
- Effect sizes: median paired difference with patient-level bootstrap 95% CI (≥ 2,000 resamples).
- Small n is expected. **Estimation with CIs** is the main output; p-values are secondary.

### 7.4 Controls

| Control | Construction | Expected if the model is real |
|---|---|---|
| Permuted environment | Shuffle $z$ across training patients within each fold | M3's gain over M2 disappears |
| Mismatched pairs | Predict patient *i*'s on-treatment composition from patient *j*'s pre-treatment composition | Error rises; quantifies patient-specific persistence |
| Disjoint programs | Re-score with disjoint program versions | Conclusions stable |
| Marker hold-out | Validate state labels with markers excluded from state definition | Agreement ≥ pre-set κ |
| Synthetic truth | Parameter recovery (T3.2) | Bias near zero; nominal coverage |

### 7.5 Sensitivity analyses (pre-specified)

1. State granularity: K = 4 vs finer substates.
2. Minimum CD8 count: 15 / 30 / 50.
3. Scoring method: AUCell / UCell / ULM. Calibration: within lineage vs within state.
4. Program versions: full vs disjoint.
5. $\Delta t$: per-patient vs nominal.
6. Site: exclude pairs with site mismatch; lymph-node-only vs non-lymph-node-only.
7. Chemistry: exclude chemistry-mismatched pairs.
8. Mechanism variants V-T … V-TPI; prior widths.
9. Environment feature set: primary vs exploratory additions.
10. Combination-therapy pairs analysed separately (exploratory).

---

## 8. Spatial coupling (Phase 5, after Gate G4)

### 8.1 Rationale and gate

Space is added **only** to test whether local context (contacts, gradients, geometry) improves prediction or explains patterns that the non-spatial model cannot express. Passing Phase 1 bundle validation (a technical import test) is **not** evidence of biological adequacy.

### 8.2 New runtime (not a patch of the validated engine)

Following VirtualTissue-Xenium's maintenance rules, dynamic states require a **new, versioned runtime**. Existing bundles become incompatible by design.

**Agent state:**
- CD8 state $k$ and time-in-state;
- the four functional axes (dysfunction, proliferation propensity, cytotoxic capacity, persistence), mapped from state and avatar by a declared table;
- `clone_id` (simulated), patient of origin, and avatar ID (reference profile), kept as separate fields.

**Transition hazards:** derived from the frozen non-spatial model, $h_{jk}(t) = q^0_{jk} \exp(\beta^\top z_{local}(t))$, where $z_{local}$ are the local counterparts of the sample-level features:
- local tumour density;
- local IFN-like field produced by myeloid agents;
- site template.

**Consistency requirement:** when local features are replaced by their sample means, the ABM must reproduce M3 predictions within tolerance (mean-field check, T5.4). Differences beyond that are attributable to spatial heterogeneity (Jensen gap) and are reported.

**Mechanisms kept minimal at first:**
- antigen recognition (declared fraction of reactive CD8, per agent);
- contact;
- checkpoint inhibition (PD-1/PD-L1 term modulated by anti-PD-1 with a schedule);
- one or two environmental conditions.

### 8.3 Initialisation: virtual biopsies

- Agents are sampled from each held-out patient's **pre-treatment** cells by stratified proportional sampling within the patient. R virtual biopsies × S seeds per patient.
- The molecular reference (avatars, calibrations) comes from **training** patients only.
- **Geometry options:**
  - G-a: synthetic templates (inflamed / excluded / desert; lymph node vs skin);
  - G-b: Xenium geometry from the same patient, if available;
  - G-c: Xenium geometry from other patients (counterfactual arrangement).

### 8.4 SPECIAL Xenium audit and spatial validation targets (T5.1)

- Verify the panel (reported as a skin panel plus custom add-on), segmentation, sample–patient–timepoint mapping and coverage of state markers.
- With a few-hundred-gene panel, AUCell is fragile. Assign spatial CD8 states by marker rules or label transfer from scRNA-seq, with uncertainty.
- The VirtualTissue-Xenium preparation is written for the Prime 5K schema, so an adapter is required.
- **Pre-specify ≤ 2 spatial targets** before simulation, for example:
  - (i) enrichment of dysfunctional-state CD8 near malignant cells relative to memory-state CD8;
  - (ii) NK exclusion from tumour parenchyma in non-responders, as reported by Roda et al. (2025).
  - Only targets supported by the cohort's own spatial data are used.

### 8.5 Gate G5 (pre-specified)

The spatial model is retained if **both** hold:
- (a) it is non-inferior to the frozen non-spatial model on held-out patients (Aitchison error margin fixed in T5.2, e.g., +10%);
- (b) it reproduces at least one pre-specified spatial target better than a spatially shuffled control.

### 8.6 Documentation

The ABM is described with the **ODD protocol** (Overview, Design concepts, Details; Grimm et al., 2020) as Supplementary Note S10.

---

## 9. Jev as an alternative transition policy (Phase 6, after Gate G5 or G4)

- **Where Jev enters:** at fate checkpoints, Jev chooses among the allowed moves {stay, move to neighbour state *k* ∈ E} given:
  - the agent's state programs;
  - time-in-state;
  - local exposures.
  The kernel executes the move. Jev therefore **participates in determining fate**. It is compared with, not substituted for, the estimated probabilities.
- **Comparators:**
  - uniform among allowed moves;
  - frozen estimated hazards (M3 / ABM);
  - Jev with named states;
  - Jev with **blind labels** (S1 … SK);
  - prompt-only control (identical menus, one irrelevant field changed).
- **Determinism and cost:**
  - repeat identical requests (5 × 20 questions) to measure nondeterminism;
  - pilot ≤ 100 calls;
  - hard request caps per arm and a common stop time;
  - hash-chained journal, as in VirtualTissue-Xenium.
- **Gate G6:** Jev is retained only if non-inferior on held-out metrics. Any claimed advantage must survive blind labels; otherwise it reflects literature priors.

---

## 10. Risk register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Too few eligible monotherapy pairs or CD8 cells | Medium | High | Gate G1; broader states; report estimation, not significance |
| Site or lesion mismatch between timepoints | High (26/46 lymph node biopsies) | High | Flag; sensitivity analyses; site covariate |
| Chemistry confounded with timepoint | Medium | Medium | Audit (T1.6); stratify or adjust |
| States not separable or unstable | Medium | High | Coarser states; undirected graph |
| Mechanisms unidentifiable | High | Medium | Report equivalence classes; priors with sensitivity |
| Leakage through calibration or feature choice | Medium | High | Fold-aware fitting; lockbox; access log |
| NK feature circularity | High | Low | Exploratory only |
| Licence misuse (NC) | Low | High | Separate from commercial work; licence notices in derived files |
| `.rds` unreadable with current Seurat | Medium | Medium | Pin compatible R/Seurat; record versions |
| Jev reproduces textbook priors | High | Medium | Blind labels; prompt-only controls |

---

## 11. Task board

Legend: `[ ]` to do · `[~]` in progress · `[x]` done · **G** = gate (stop and decide).

### Phase 0: Setup and governance

- [ ] **T0.1 — Repository and environment** · S · Depends: —
  - **Goal:** an analysis repository with pinned environments.
  - **Steps:**
    - Decide the location (Decision D9): recommended a new analysis repository `VirtualTissue-SC` (Phases 1–4), with Phase 5 implemented as a new runtime inside the VirtualTissue platform.
    - Create the layout in Appendix A.
    - Pin Python (3.12; scanpy, anndata, decoupler, numpy/scipy, numpyro or cmdstanpy, scikit-learn, pyarrow) and R (version, Seurat, SeuratObject, Matrix).
    - Add a `config/study.yaml` (Appendix B).
  - **Outputs:** repository skeleton; `environment.yml`; `renv.lock` or `sessionInfo.txt`; `config/study.yaml`.
  - **Acceptance:** a clean install reproduces the environment; tests run (even if empty).
  - **Supplement:** S12 (software).

- [ ] **T0.2 — Data acquisition and integrity** · S · Depends: T0.1
  - **Goal:** priority-1 files obtained and verified.
  - **Steps:**
    - Download manually from the RDR page (browser) into `<DATA_ROOT>/SPECIAL/raw/`.
    - Record the URL, dataset version, date, byte size and SHA-256 in `data/MANIFEST.json`.
    - Obtain clinical Table S1 (2024 *Cell*) and the 2025 supplements.
  - **Outputs:** `data/MANIFEST.json`; `data/clinical/` (original files, unmodified).
  - **Acceptance:** sizes match the RDR listing; hashes recorded.
  - **Supplement:** S1, Supplementary Table S1.

- [ ] **T0.3 — Licence and data-use note** · S · Depends: T0.2
  - **Goal:** written terms for derived data.
  - **Outputs:** `DATA_LICENSE.md` (CC BY-NC-SA 4.0 notice, attribution text, NC restriction).
  - **Supplement:** S1, Data availability statement.

- [ ] **T0.4 — Pre-registration freeze** · S · Depends: T0.1
  - **Goal:** freeze this plan before any outcome-relevant analysis.
  - **Steps:**
    - Commit this file and `config/study.yaml`; record their SHA-256 in `prereg/FREEZE.json` with the date.
    - Optional: deposit on OSF or Zenodo for a timestamp.
  - **Acceptance:** freeze hash exists before T1.4 runs.
  - **Supplement:** S7 (pre-specification statement).

### Phase 1: Feasibility audit → first deliverable: Feasibility report

- [ ] **T1.1 — Export Seurat objects to a neutral format** · M · Depends: T0.2
  - **Steps:**
    - In R: read `GC_all_immune_sharing.rds`; record object class, Seurat version, assays, layers, `meta.data` columns and reductions.
    - Export raw counts (RNA assay, counts layer) per sample as MatrixMarket; write `obs.parquet` (all metadata), `var.tsv` and embeddings.
    - Verify round-trip: total counts per cell equal across formats.
  - **Outputs:** `derived/immune/<sample>/…`; `reports/T1.1_object_inventory.md` (assays, metadata dictionary with a description per column).
  - **Acceptance:** 100% of cells exported; count totals identical; metadata dictionary complete.
  - **Supplement:** S2, Supplementary Table S2 (metadata dictionary).

- [ ] **T1.2 — Patient and sample table** · M · Depends: T1.1
  - **Steps:** join object metadata with the clinical table(s). One row per biopsy with:
    - patient, timepoint, treatment arm;
    - biopsy site and lesion ID (if available);
    - date or day relative to treatment start;
    - chemistry;
    - total cells, immune cells, CD8 cells;
    - response (kept in a separate column, flagged "not an input").
  - **Reconciliation:** document every conflict between sources.
  - **Outputs:** `tables/samples.parquet`; `reports/T1.2_sample_table.md`.
  - **Acceptance:** 46 biopsies or an explained difference; arms total 17 + 6; pairing 20/23 or explained.
  - **Supplement:** Supplementary Table S3.

- [ ] **T1.3 — Eligibility and flow diagram** · S · Depends: T1.2
  - **Steps:** apply §4.1. Compute **total** CD8 per biopsy only (no state composition yet). Draw the flow diagram.
  - **Outputs:** `tables/eligibility.parquet`; `figures/S1_flow.svg`.
  - **Supplement:** Supplementary Figure S1.

- [ ] **T1.4 — Lockbox assignment** · S · Depends: T1.3, T0.4
  - **Steps:** stratified random assignment of 4 eligible monotherapy patients (§7.1) with the recorded seed. Write IDs to `prereg/LOCKBOX.json` with a hash. All later code must refuse to load lockbox on-treatment cells unless `--final-evaluation` is set, and must log any such access.
  - **Acceptance:** assignment hash recorded before T1.5.
  - **Supplement:** S7.

- [ ] **T1.5 — CD8 state availability (development patients only)** · S · Depends: T1.4
  - **Steps:**
    - Counts per state × biopsy for development patients.
    - Number of patients contributing ≥ 5 cells of each state at each timepoint.
    - For lockbox patients: totals only.
  - **Outputs:** `tables/state_counts_dev.parquet`; `figures/S2_state_counts.svg`.
  - **Supplement:** Supplementary Figure S2.

- [ ] **T1.6 — Technical confounders** · M · Depends: T1.2
  - **Steps:**
    - Cross-tabulate chemistry × timepoint × patient.
    - Site and lesion match per pair.
    - QC distributions (genes, UMIs, mitochondrial fraction) by chemistry and site.
    - Dissociation-program levels by sample.
  - **Outputs:** `reports/T1.6_confounders.md`; `figures/S3_confounders.svg`.
  - **Acceptance:** each confounder classified as negligible, adjust, stratify or exclude, with justification.
  - **Supplement:** S3; Supplementary Figure S3.

- [ ] **T1.7 — Feasibility report · Gate G1** · S · Depends: T1.3–T1.6
  - **Contents:**
    - eligible pairs;
    - cells per state and patient;
    - program coverage (preview from T2.2 if ready);
    - train/test separation;
    - evidence that an estimable change exists (development patients only: distribution of per-patient CLR changes vs mismatched-pair noise).
  - **G1 decision:**
    - **Go** if ≥ 10 eligible monotherapy pairs in development and each primary state is present in ≥ 70% of development biopsies.
    - **Otherwise:** coarsen states, relax the threshold (documented), or stop.
  - **Outputs:** `reports/T1.7_feasibility.md`. Save it to the Project as well.

### Phase 2: CD8 state reference

- [ ] **T2.1 — CD8 subset QC and annotation audit** · M · Depends: G1
  - **Steps:**
    - Subset CD8 cells (author labels); check CD8A/CD8B vs CD4.
    - Detect doublets and contamination (e.g., NK/CD8 mixing; the published NK compartment is relevant).
    - Audit state labels against held-out markers.
  - **Outputs:** `derived/cd8/cd8.h5ad`; `reports/T2.1_annotation_audit.md`.
  - **Acceptance:** documented label agreement (κ) with marker-based assignment.
  - **Supplement:** S4.

- [ ] **T2.2 — Program library and coverage** · M · Depends: T2.1
  - **Steps:** assemble programs (§5.2) from cited sources with versions; compute coverage; compute pairwise overlap (Jaccard); create disjoint variants.
  - **Outputs:** `resources/programs_v1.json` with provenance; `tables/program_coverage.tsv`.
  - **Supplement:** Supplementary Table S4.

- [ ] **T2.3 — Scoring and fold-aware calibration** · M · Depends: T2.2
  - **Steps:**
    - Choose `n_up` from the depth distribution.
    - Implement fold-aware percentile calibration (within lineage × depth).
    - Run the chemistry mixed-model check.
    - Run alternative scoring methods.
  - **Outputs:** `src/vtsc/scoring.py` with tests; `reports/T2.3_scoring.md`.
  - **Acceptance:** a unit test confirms that test-patient cells never enter the reference.
  - **Supplement:** S4.

- [ ] **T2.4 — State graph, pseudotime and robustness** · L · Depends: T2.3
  - **Steps:** apply §5.4: PAGA with patient bootstrap; three pseudotime methods; multiple roots and representations; apply the defensible-direction criteria.
  - **Outputs:** `resources/state_graph_v1.json` (edges, confidence, directed flag); `figures/S4_graph.svg`; `tables/S5_robustness.tsv`.
  - **Gate G2a:** graph frozen (directed or undirected).
  - **Supplement:** S5; Supplementary Figure S4; Supplementary Table S5.

- [ ] **T2.5 — Functional axes per state** · S · Depends: T2.3
  - **Outputs:** `tables/functional_axes.tsv`; `figures/S5_axes.svg`. In particular, the proliferation of dysfunctional cells.
  - **Supplement:** S4.

- [ ] **T2.6 — Environment features** · M · Depends: T1.2; requires `Entire_TME.rds` and `Malignant_cells.rds` exported as in T1.1
  - **Steps:** match IDs across objects; compute the primary $z_1$–$z_3$ and the exploratory features; fold-aware standardisation code.
  - **Outputs:** `tables/environment_features.parquet`; `reports/T2.6_features.md`.
  - **Supplement:** S6; Supplementary Table S6.

- [ ] **G2 — State reference frozen** · Depends: T2.1–T2.6
  - Freeze `programs_v1.json`, `state_graph_v1.json` and the feature definitions (hashes in `prereg/FREEZE.json`).

### Phase 3: Non-spatial dynamics model

- [ ] **T3.1 — Implement the model ladder** · L · Depends: G2
  - **Steps:** M0, M1 (and M1b), M2, M3 with the Dirichlet-multinomial likelihood; latent $\pi^{pre}$; matrix exponential; inference backend chosen and pinned.
  - **Outputs:** `src/vtsc/models/`; tests including analytic cases (e.g., K = 2).
  - **Supplement:** S6 (equations exactly as implemented).

- [ ] **T3.2 — Simulation-based parameter recovery** · M · Depends: T3.1
  - **Steps:** simulate cohorts matching the observed n, $N_i$, $\Delta t_i$ and dispersion for variants V-T … V-TPI; fit; report bias, 80%/95% coverage and correlations.
  - **Outputs:** `reports/T3.2_recovery.md`; `figures/S6_recovery.svg`.
  - **Supplement:** S8.

- [ ] **T3.3 — Identifiability on real data** · M · Depends: T3.2
  - **Steps:** profile likelihoods or posterior correlations; define the equivalence classes of mechanisms.
  - **Outputs:** `reports/T3.3_identifiability.md`.
  - **Supplement:** S8.

- [ ] **T3.4 — Development LOPO evaluation** · M · Depends: T3.1, T2.6
  - **Steps:** run §7.2–§7.4 on development patients; all controls.
  - **Outputs:** `results/dev_lopo/*.parquet`; `reports/T3.4_dev_results.md`.
  - **Supplement:** S7; Supplementary Table S7.

- [ ] **T3.5 — Sensitivity analyses** · L · Depends: T3.4
  - **Steps:** items 1–10 of §7.5.
  - **Outputs:** `results/sensitivity/`; `tables/S8_sensitivity.tsv`.
  - **Supplement:** S9; Supplementary Table S8.

- [ ] **T3.6 — Freeze the model · Gate G3** · S · Depends: T3.4, T3.5
  - **Steps:** choose the final model family and variant by the pre-specified rule: best development LOPO log predictive density, with ties going to the simpler model. Record the hash.
  - **G3 decision:** proceed to the lockbox regardless of the result (no outcome-dependent stopping); record the choice.

### Phase 4: Lockbox and external validation

- [ ] **T4.1 — Lockbox evaluation** · S · Depends: G3
  - **Steps:** run the frozen pipeline with `--final-evaluation`; also run full LOPO on all eligible patients.
  - **Outputs:** `results/final/`; `reports/T4.1_final.md` (H1/H2 results with CIs).
  - **Supplement:** S7; main-text Figure 3.

- [ ] **T4.2 — External direction check (GSE120575)** · M · Depends: G3
  - **Steps:** map CD8 states by label transfer (training reference only); compare the predicted direction of change with observed pre/post changes; report mapping confidence.
  - **Supplement:** S11 (external validation).

- [ ] **T4.3 — Complementary temporal and clonal data (optional)** · M · Depends: G3
  - **Steps:**
    - GSE272993: plausibility of timescales and influx priors in blood.
    - GSE123813: clonal replacement as a check on the influx component.
  - **Supplement:** S11.

- [ ] **G4 — Go/no-go for spatial work**
  - Proceed if H1 holds, or if the non-spatial model is at least informative (CIs reported). Otherwise, write up the non-spatial results and stop.

### Phase 5: Spatial coupling

- [ ] **T5.1 — SPECIAL Xenium audit** · M · Depends: G4
  - **Steps:** verify the panel, segmentation, sample mapping and marker coverage; define ≤ 2 spatial targets (§8.4) before any simulation.
  - **Outputs:** `reports/T5.1_xenium_audit.md`; `prereg/SPATIAL_TARGETS.json`.

- [ ] **T5.2 — Runtime design document (ODD) and margins** · M · Depends: T5.1
  - **Steps:** write the ODD description (§8.2–§8.3); fix the non-inferiority margin and gate criteria; plan the version bump.
  - **Outputs:** `docs/ODD_runtime_v3.md`.
  - **Supplement:** S10.

- [ ] **T5.3 — Implement the dynamic-state runtime** · L · Depends: T5.2
  - **Steps:** new runtime version; agents with state, time-in-state, functional axes, clone and avatar fields; hazards from the frozen model; local exposures; anti-PD-1 schedule; event recording for replay.
  - **Acceptance:** engine tests; integrity hashes; no change to the validated 0.1.2 runtime.

- [ ] **T5.4 — Mean-field consistency test** · M · Depends: T5.3
  - **Steps:** with local features set to sample means, the aggregated ABM output must match M3 predictions within tolerance. Report the Jensen gap when local features vary.

- [ ] **T5.5 — Geometry experiments** · L · Depends: T5.4
  - **Steps:** G-a/G-b/G-c × R virtual biopsies × S seeds per held-out patient; spatially shuffled control.
  - **Outputs:** `results/spatial/`.

- [ ] **T5.6 — Spatial evaluation · Gate G5** · M · Depends: T5.5
  - **Steps:** apply §8.5. Write `reports/T5.6_spatial.md`.

### Phase 6: Jev policy

- [ ] **T6.1 — Pilot (≤ 100 calls, explicit approval)** · S · Depends: G5 (or G4 if spatial is skipped)
  - **Steps:** prompt `vtsc-transition-1` (no outcome targets); blind-label variant; nondeterminism test; cost log.

- [ ] **T6.2 — Policy comparison · Gate G6** · M · Depends: T6.1
  - **Steps:** comparators of §9 on held-out patients with matched request caps and a common stop time.
  - **Supplement:** S12 (Jev section).

### Phase 7: Deferred extensions (not started before G4)

- [ ] T7.1 Combination-therapy pairs (separate model or test of transfer).
- [ ] T7.2 Site-specific templates (lymph node vs skin) beyond the covariate.
- [ ] T7.3 Healthy-skin context (needs a tissue profile, new data, new menus).
- [ ] T7.4 Explicit untreated counterfactual arm (labelled hypothetical).
- [ ] T7.5 Coupling to VirtualTissue-LN (priming compartment → influx).

### Phase 8: Paper and supplement assembly

- [ ] **T8.1 — Compile supplementary notes** (§12) from the task reports · M
- [ ] **T8.2 — Figures**: main and supplementary, from saved results only, with scripts in `figures/` · M
- [ ] **T8.3 — Data and code availability**: Zenodo DOI for code; derived data under CC BY-NC-SA · S
- [ ] **T8.4 — Reproducibility run**: from a clean environment, regenerate every table and figure; compare hashes · M
- [ ] **T8.5 — Claims audit**: claims ledger (§12.4); every claim maps to a result file and an evidence class. Optional: run the comp-bio-slop-check skill on the manuscript draft · S

---

## 12. Supplementary-material blueprint

### 12.1 Supplementary Notes

| Note | Title | Built from |
|---|---|---|
| S1 | Data sources, versions, licences and integrity | T0.2, T0.3 |
| S2 | Object inventory and metadata dictionary | T1.1 |
| S3 | Cohort, eligibility and technical confounders | T1.2–T1.6 |
| S4 | CD8 state reference: annotation audit, programs, scoring and calibration, functional axes | T2.1–T2.3, T2.5 |
| S5 | State graph, pseudotime and robustness | T2.4 |
| S6 | Environment features and dynamics-model specification (equations, priors, inference) | T2.6, T3.1 |
| S7 | Pre-specification, splits, lockbox, metrics and statistics | T0.4, T1.4, T3.4, T4.1 |
| S8 | Parameter recovery and identifiability | T3.2, T3.3 |
| S9 | Sensitivity analyses | T3.5 |
| S10 | Spatial agent-based model (ODD protocol) and mean-field consistency | T5.2–T5.6 |
| S11 | External and complementary validation | T4.2, T4.3 |
| S12 | Software, environments, compute, Jev provider settings and cost | T0.1, T6.x, T8.4 |
| S13 | Claim boundaries and limitations | §1.4, §2.1, T8.5 |

### 12.2 Supplementary Tables

| Table | Content |
|---|---|
| S1 | Data manifest (file, version, bytes, SHA-256, licence, URL) |
| S2 | Metadata dictionary of deposited objects |
| S3 | Per-biopsy table (patient, timepoint, arm, site, chemistry, cell counts; response flagged) |
| S4 | Program library with genes, sources, versions, coverage and overlaps |
| S5 | State-graph robustness (edges × methods × bootstrap frequency) |
| S6 | Environment feature values per biopsy (development; lockbox after T4.1) |
| S7 | Per-patient predictions and errors for every model (LOPO and lockbox) |
| S8 | Sensitivity-analysis results |
| S9 | Parameter estimates with intervals, per variant |
| S10 | Spatial run registry (bundle hashes, seeds, settings) |
| S11 | Jev run registry (prompt version, requests, cost, nondeterminism) |

### 12.3 Supplementary Figures

| Figure | Content |
|---|---|
| S1 | Flow diagram |
| S2 | CD8 state counts per biopsy |
| S3 | Technical confounders |
| S4 | State graph with edge confidence |
| S5 | Functional axes per state |
| S6 | Parameter recovery |
| S7 | Per-patient observed vs predicted compositions (all models) |
| S8 | Sensitivity forest plot |
| S9 | Spatial targets: observed vs simulated |
| S10 | Jev vs estimated policies |

Draft main-text figures:
1. Design and data.
2. State reference.
3. Held-out prediction (H1/H2).
4. Mechanism attribution and identifiability.
5. Spatial extension (if retained).

### 12.4 Claims ledger (template)

| Claim | Evidence class | Result file | Task | Status |
|---|---|---|---|---|
| e.g., "M3 reduces held-out Aitchison error vs M2 by X (95% CI)" | Estimated, held-out | `results/final/h2.parquet` | T4.1 | pending |

### 12.5 Draft availability statements

- **Data:** "Processed SPECIAL objects are available from KU Leuven RDR (doi:10.48804/GSAXBN, CC BY-NC-SA 4.0). Raw sequencing data are under controlled access at EGA (EGAS00001006488). Derived data are released under CC BY-NC-SA 4.0 at [Zenodo DOI]."
- **Code:** "Code is available at [repository] under the MIT licence, archived at [Zenodo DOI]."

### 12.6 Reporting standards to follow

- **ODD protocol** for the agent-based model.
- **TRIPOD+AI**: adapted items relevant to a prediction model evaluated on held-out patients (data, predictors, missingness, validation, performance).
- **FAIR** principles for released derived data.

---

## 13. Decision log and deviation log

### 13.1 Decisions (seeded from planning on 2026-10-27)

| ID | Decision | Rationale |
|---|---|---|
| D1 | SPECIAL (Poźniak/Marine) is the primary dataset | Paired pre/early on-treatment biopsies; public processed objects (immune, malignant, TME); spatial data in the same study |
| D2 | Non-spatial prediction first; space and Jev later | A predictive baseline is needed before mechanism layers can be judged |
| D3 | Pre-treatment is not an untreated control; no surrogate untreated arm | A biopsy at baseline does not tell what happens without therapy |
| D4 | Functional axes are separate (dysfunction ≠ non-proliferative) | Li et al., 2019 |
| D5 | Start from observed precursor states; no TCF7+ = Tpex assumption; priming not restricted to lymph nodes | Thompson et al., 2010; annotation audit |
| D6 | Avatars constrain state, not history; identity, clone and profile kept separate | Avoids conflating simulated jumps with observed lineage |
| D7 | Fit to observed change between timepoints, not to an assumed equilibrium | One biopsy does not show steady state |
| D8 | Jev only as a comparator after estimated probabilities | Jev choosing branches determines fate; it must be benchmarked |
| D9 | Code location: new analysis repository for Phases 1–4 (recommended), spatial runtime inside the VirtualTissue platform | Keeps the paper analysis independent of engine versioning (to be confirmed) |
| D10 | Anti-PD-1 monotherapy only for the primary analysis; combination analysed separately | Different intervention |

### 13.2 Deviation log (template)

| Date | Task | Planned | Actual | Reason | Impact on claims |
|---|---|---|---|---|---|

---

## 14. Glossary

- **Aitchison distance:** Euclidean distance between centred log-ratio transforms of two compositions.
- **Avatar:** a real training cell whose profile represents a simulated agent's current state.
- **Equivalence class (mechanisms):** model variants whose held-out predictions cannot be distinguished with the available data.
- **LOPO:** leave-one-patient-out cross-validation.
- **Lockbox:** patients withheld from all development and evaluated once after freezing.
- **Mean-field consistency:** the spatial model reduces to the non-spatial model when local inputs are replaced by sample means.
- **Virtual biopsy:** a stratified random subset of one patient's cells used to initialise a simulation.

---

## 15. References (verify each before citing)

**Cohort and related**
- Poźniak J, et al. A TCF4-dependent gene regulatory network confers resistance to immunotherapy in melanoma. *Cell* 2024. PMID 38181739.
- Roda N, et al. Cytotoxic NK cells impede response to checkpoint immunotherapy in melanoma with an immune-excluded phenotype. 2025. PMID 40530504; PMC12409281. Preprint: bioRxiv 10.1101/2023.12.14.571631.
- SPECIAL scRNA-seq dataset, EGA EGAD00001009291 (study EGAS00001006488).
- Poźniak J. Dissecting the melanoma ecosystem one cell at the time during immunotherapy. KU Leuven RDR, doi:10.48804/GSAXBN (CC BY-NC-SA 4.0).

**CD8 biology**
- Li H, et al. Dysfunctional CD8 T cells form a proliferative, dynamically regulated compartment within human melanoma. *Cell* 2019.
- Thompson ED, et al. Tumor masses support naive T cell infiltration, activation, and differentiation into effectors. *J Exp Med* 2010.
- Sade-Feldman M, et al. Defining T cell states associated with response to checkpoint immunotherapy in melanoma. *Cell* 2018.
- Yost KE, et al. Clonal replacement of tumor-specific T cells following PD-1 blockade. *Nat Med* 2019.
- Miller BC, et al. Subsets of exhausted CD8+ T cells differentially mediate tumor control and respond to checkpoint blockade. *Nat Immunol* 2019.
- Siddiqui I, et al. Intratumoral Tcf1+PD-1+CD8+ T cells with stem-like properties promote tumor control in response to vaccination and checkpoint blockade immunotherapy. *Immunity* 2019.
- Pauken KE, et al. Epigenetic stability of exhausted T cells limits durability of reinvigoration by PD-1 blockade. *Science* 2016.
- Simoni Y, et al. Bystander CD8+ T cells are abundant and phenotypically distinct in human tumour infiltrates. *Nature* 2018.
- Zheng L, et al. Pan-cancer single-cell landscape of tumor-infiltrating T cells. *Science* 2021.

**Methods**
- Aibar S, et al. SCENIC (AUCell). *Nat Methods* 2017.
- Badia-i-Mompel P, et al. decoupleR. *Bioinformatics Advances* 2022.
- Andreatta M, Carmona SJ. UCell. *Comput Struct Biotechnol J* 2021.
- Wolf FA, et al. PAGA. *Genome Biol* 2019.
- Haghverdi L, et al. Diffusion pseudotime. *Nat Methods* 2016.
- Setty M, et al. Palantir. *Nat Biotechnol* 2019.
- Street K, et al. Slingshot. *BMC Genomics* 2018.
- Weiler P, et al. CellRank 2. *Nat Methods* 2024.
- Schiebinger G, et al. Waddington-OT. *Cell* 2019.
- Klein D, et al. moscot. *Nature* 2025 (cross-timepoint optimal-transport baseline; verify).
- Büttner M, et al. scCODA. *Nat Commun* 2021.
- Dann E, et al. Milo. *Nat Biotechnol* 2022.
- Aitchison J. *The Statistical Analysis of Compositional Data.* 1986.
- van den Brink SC, et al. Single-cell sequencing reveals dissociation-induced gene expression in tissue subpopulations. *Nat Methods* 2017.
- Talts S, et al. Validating Bayesian inference algorithms with simulation-based calibration. arXiv 2018.

**Reporting**
- Grimm V, et al. The ODD protocol for describing agent-based and other simulation models: a second update. *JASSS* 2020.
- Collins GS, et al. TRIPOD+AI statement. *BMJ* 2024.

---

## Appendix A — Repository layout (proposed)

```text
VirtualTissue-SC/
├── config/study.yaml          # seeds, paths, thresholds, states, features, metrics
├── prereg/                    # FREEZE.json, LOCKBOX.json, SPATIAL_TARGETS.json
├── data/MANIFEST.json         # no data in git; manifest + hashes only
├── DATA_LICENSE.md
├── r/export_seurat.R          # T1.1
├── src/vtsc/
│   ├── io.py                  # neutral-format readers, lockbox guard
│   ├── cohort.py              # T1.2–T1.4
│   ├── scoring.py             # T2.3 (fold-aware)
│   ├── graph.py               # T2.4
│   ├── features.py            # T2.6
│   ├── models/                # M0–M3, variants, inference
│   ├── evaluate.py            # LOPO, lockbox, metrics, controls
│   └── report.py
├── resources/                 # programs_v1.json, state_graph_v1.json
├── results/                   # dev_lopo/, sensitivity/, final/, spatial/
├── figures/                   # scripts + outputs
├── reports/                   # one Markdown report per task
└── tests/
```

## Appendix B — `config/study.yaml` (initial values)

```yaml
master_seed: 20261027
dataset:
  doi: 10.48804/GSAXBN
  version: 8
  licence: CC-BY-NC-SA-4.0
cohort:
  arm_primary: nivolumab_monotherapy
  min_cd8_per_biopsy: 30
  sensitivity_min_cd8: [15, 50]
  lockbox_n: 4
  lockbox_strata: [response, site_pairing]
states:
  primary: [Memory, Activated, Cytotoxic, Dysfunctional]   # confirm in T2.1
  cycling_as_axis: true
scoring:
  method: aucell
  n_up: auto_from_depth
  calibration: within_lineage_by_depth_decile
  fold_aware: true
graph:
  bootstrap: 200
  edge_min_frequency: 0.8
  direction_rule: {min_spearman: 0.8, bootstrap_stable: true}
features:
  primary: [tumour_content, myeloid_ifng, site_lymph_node]
  exploratory: [nk_cytotoxic_fraction, malignant_state_mix, malignant_mhc1, cd8_fraction_immune, treg_fraction]
model:
  ladder: [M0, M1, M1b, M2, M3]
  variants: [V-T, V-P, V-I, V-TP, V-TPI]
  likelihood: dirichlet_multinomial
  pseudocount: 0.5
evaluation:
  primary_metric: aitchison
  key_secondary: log_predictive_density
  tests: fixed_sequence_wilcoxon_one_sided
  alpha: 0.05
  bootstrap: 2000
  interval: 0.8
jev:
  enabled: false
  pilot_cap: 100
```

## Appendix C — Session checklist (copy into each new conversation)

1. Read this plan and the latest task reports in the Project.
2. Show the task board status and propose the next task.
3. Confirm any decision marked "to be confirmed" before acting.
4. Run the task; save outputs with hashes; write `reports/<TASK>.md`.
5. Update this plan's ticks, Decision log and Deviation log; save to the Project.
6. Never: touch lockbox on-treatment data early, use response as input, spend Jev credits without approval, or edit the validated engine.
