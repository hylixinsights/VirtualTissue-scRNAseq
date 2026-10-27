# Cohort expansion review — 2026-10-27

The user selected revising/expanding the primary cohort after the SPECIAL
feasibility failure. This review preserves the original freeze. The user subsequently approved an
expanded **melanoma ICI audit with separate monotherapy and combination strata**.
No replacement cohort has yet passed eligibility. The CD8 threshold and early
biopsy target have not been relaxed.
This recommendation concerns feasibility auditing, not a claim that the expanded
model is adequately powered or identifiable.

## Result after auditing the deposited object

**The expanded cohort still fails the existing numerical prerequisite.**
With ≥30 author-labelled CD8 per biopsy, the available sources supply at most
**12 count-qualified patients: 4 monotherapy and 8 combination**. Reserving four
patients leaves at most eight for development, below the existing minimum of
ten. These are provisional candidates, not verified eligible patients. The
minimum is a feasibility floor, not a power calculation; stratification adds
requirements beyond that floor.

| Stratum | SPECIAL | New Lin patients | Combined count candidates |
|---|---:|---:|---:|
| Anti-PD-1 monotherapy | 3 | 1 (P7) | **4** |
| Anti-PD-1 + anti-CTLA-4 | 4 | 4 (P1, P2, P3, P5) | **8** |
| Total | 7 | 5 | **12** |

The primary output is `reports/evidence/cohort-expansion/lin-object-count-audit.json`;
sample counts, patient counts and clinical crosswalks are in
`tables/cohort_expansion_lin_object_{samples,patients}.tsv` and
`tables/cohort_expansion_lin_library_crosswalk.tsv`. At the pre-specified
sensitivity thresholds, totals are 17 at ≥15 (8 Mono, 9 Combo) and 11 at ≥50
(3 Mono, 8 Combo). The ≥15 result does not change the primary threshold or pass G1.

The integrated object has **25 in-house samples from 17 patients, only 8 with
both Pre and On**, whereas S1 lists 34 clinical rows and 16 clinically paired
patients. Eight patient/timepoint keys are absent from the deposit: P4 On,
P6 On, P8 On, P10 On, P11 On, P14 Pre, P15 Pre and P17 On. Absence is not a zero
cell count and its cause is unknown. P2 has one deposited On label but two S1
On biopsies (days 9 and 42); neither the sample ID nor the appropriate early
interval can be assigned. All clinical joins use provisional patient-number and
timepoint matching because the object lacks the E193xx/E233xx library IDs.

Lineage audit: coarse `celltype_main=CD8+T` also contains `gdT`, while some
`Cycling T/NK` cells have explicit CD8 labels in `celltype_r2`. Before recounting,
the audit declared `celltype_r2` beginning with `CD8_` as the specific author
lineage rule, including identified cycling CD8 and excluding `gdT`. The coarse
rule and their permissive union are saved as diagnostic comparisons. **All three
rules yield the same qualifying patients at 15, 30 and 50.** The original SPECIAL
rule is unchanged; its generic cycling population lacks CD8-specific identity.
Only total lineage counts and label vocabulary were inspected, not per-state
proportions or changes. These atlas annotations do not become a fitted state
reference for a future prediction study.

The atlas also includes **19 SPECIAL patients / 38 samples** under
`SKCM_Plozniak`. They are explicitly excluded from new-patient totals; the
original SPECIAL counts remain the primary source. Full cross-study patient
independence, exact drugs, histology, lesion matching and clinical library links
remain unverified. SPECIAL inter-biopsy dates are still missing. Lin chemistry
also has a source conflict: EGA names 5′ v2, while the publisher key-resource
table names 5′ v1.1; neither is silently selected.

## What the downloaded clinical sources establish

The initial clinical numbers below are generated in `reports/evidence/cohort-expansion/audit.json`
and the clinical `tables/cohort_expansion_*` tables. Sources, hashes, download
status and URLs are in `data/COHORT_EXPANSION_MANIFEST.json`. The clinical PDF was
rendered and its full table checked visually. Original files remain unmodified.

| Source | Verified scope | Contribution to the present primary question |
|---|---|---|
| SPECIAL | 3 monotherapy patients with ≥30 author-annotated CD8 in both biopsies | Count-only candidates; individual intervals still missing |
| Lin et al., Cancer Cell 2026, in-house cohort | 34 clinical rows, 17 patients, **16 clinical pairs**, of whom **6 are Mono and 10 Combo** | Clinical bound superseded by the deposited-object count audit above; specific drugs remain unresolved |
| Sade-Feldman et al., Cell 2018 | 32 clinical patients; 11 pre/post pairs; **5 pairs labelled PD1 without a CTLA4/sequential label** | Later-time candidate/transport analysis; cannot repair the early-biopsy cohort unchanged |

Lin S1 identifies P7–P12 as Mono, all with baseline day 0 and on-treatment day 7.
The source labels these pairs `Matched=Yes`; the exact operational meaning and
lesion identity still require reconciliation. S1 gives no nivolumab/pembrolizumab
names. P16 has baseline only; P2 has on-treatment days 9 and 42. Thus the headline
patient count must not be equated with the number of early paired scRNA patients.
Keep the day-42 sample in the inventory without selecting it for an early pair.
Vulvar/vaginal biopsy descriptions in the combination group also make histology
confirmation necessary; anatomical site alone does not establish melanoma subtype.

Lin's open [Zenodo v2 deposit](https://doi.org/10.5281/zenodo.19586163)
contains an integrated pan-cancer object. Downloadable Table S2A explicitly lists
Poźniak/SPECIAL among its reused studies. Count only independently identified
**in-house** patients as new; never add the atlas total to SPECIAL. The v1 record
is restricted, while v2 is open. The processed object has a deposited size of
4,021,369,094 bytes and MD5 `a3486788163889f8cc910ea0115d70ac`. Its download was
initially stopped after the monotherapy clinical upper bound failed. After the
user selected the expanded treatment audit, the full file was acquired using
validated HTTP ranges and **verified against that published MD5**. The source
SHA-256 is `43c9e52c5a7766fe903be29b8cbcc3aaee92c2e3a3b6801bc80a9f8b01875772`.
The earlier `.part` is retained as an invalid historical intermediate. No spatial
matrices were acquired. The successful audit extracts metadata only; full R
loading and official in-memory conversion exceeded the configured memory limits.

Sade-Feldman's five PD1-labelled pairs are P2, P3, P12, P15 and P20. First
on-treatment days are 35, 63, 77, 73 and 199; the corresponding inter-biopsy
intervals are 35, 63, **93**, 73 and 199 days because P12 baseline is day −16.
P1/P6 have sequential CTLA4/PD1 labels and are excluded from this simple
monotherapy bound. P23 has two on-treatment biopsies but no baseline. Do not
exclude patients for discordant lesion response. Smart-seq2 CD45-sorted immune
data also lack an unbiased whole-TME tumour-content denominator. Drug names,
CD8 totals and compatibility of state definitions need independent verification.

## Historical clinical bounds, superseded by object counts

These bounds motivated the audit and are preserved for provenance. They are
not the current candidate totals; the deposited-object result is 12 at ≥30.

| Hypothetical scope | Maximum candidate patients from the audited sources | Why this is not a G1 pass |
|---|---:|---|
| Early monotherapy: SPECIAL + new Lin patients | **9 = 3 + 6** | Below 10 development patients even without a lockbox; below the original 10 + 4 allocation |
| Monotherapy across early and later intervals, adding Sade-Feldman | **14 = 3 + 6 + 5** | Changes the time estimand; unverified CD8/drugs; heterogeneous assays; no sample-size margin; SPECIAL dates missing |
| Early ICI, with separate Mono/Combo strata: SPECIAL + Lin | **23 = (3 + 4) + (6 + 10)** | Changes the treatment estimand; 16 Lin pairs are clinical candidates only; timing/histology/source identity unresolved |

The four SPECIAL combination count candidates are GC28, GC32, GC38 and GC41.
These newly tabulated totals do not make them eligible for the frozen monotherapy
analysis. None of the optimistic sums establishes independence or eligibility;
all assume no duplicate patients and successful remaining checks. The original
minimum is a feasibility floor, not a power calculation.

## Other sources checked

The candidate register is `tables/cohort_expansion_candidates.tsv`.

- **GSE272993:** GEO explicitly describes PBMC sampling and a selected non-naive
  CD8 object. It cannot substitute for tumour biopsies or the full tumour CD8
  composition; retain only for complementary temporal/clonal work.
- **GSE164237:** GEO describes peripheral CD8 at baseline and 12 weeks. Wrong
  compartment for the present primary endpoint.
- **GSE123813 (Yost):** tumour single-cell/TCR data from BCC/SCC; useful for influx
  biology but requires a different cancer target, not automatic melanoma pooling.
- **Sun et al. 2024 (HRA005837):** the primary article describes 11 patients on
  toripalimab, 8 non-progressing patients with paired tissue and 3 progressing
  patients with pretreatment tissue only. The treatment is outside the approved
  drugs, and availability of the second biopsy depends on response. This is not
  an unbiased validation sample of all treated melanoma patients. Processed-data
  access was not established.
- **Wei et al. 2026 (PMC13130630):** longitudinal bulk/spatial multi-omics and
  single-cell spatial proteomics; the inspected data-availability section does
  not establish a new paired dissociated tumour scRNA cohort for this endpoint.
- **Liu et al. 2026 preprint (doi:10.64898/2026.09.12.751186):** relevant fixed-cell
  and spatial work, but page 13 says data/code will be released on peer-reviewed
  publication. No presently usable deposit was established. Treat as a lead,
  not as acquired data or additional eligible patients.

This is a targeted feasibility search, not a systematic review or proof that no
other qualifying melanoma cohort exists.

## Amendment status and next actions

The concrete draft is
`prereg/amendment-A02-cohort-expansion-draft-2026-10-27.md`. The user selected separate monotherapy/combination strata for the next audit.
The amended analysis is not frozen. Total-CD8 counts and explicit atlas reuse
are audited; raw-matrix alignment and final clinical identity remain open.
Because the necessary sample-size condition already fails, model fitting and
lockbox assignment remain stopped. The next evidence needed is:

1. Obtain the Lin library crosswalk and absent libraries or their documented
   exclusions; reconcile drugs, histology, chemistry, lesion matching and P2's
   On identity. `docs/data-request-draft.txt` is prepared, **not sent**.
2. Reconcile patient identities across source studies. Reused SPECIAL samples
   remain SPECIAL. Broad lineage labels may screen feasibility; inherited
   integrated states/embeddings may not supply a fitted reference or test labels.
3. Recompute eligibility separately by cohort and treatment using ≥30 CD8,
   retaining the 15/50 sensitivity thresholds. Resolve SPECIAL intervals before
   interpreting fitted parameters as rates per day. Do not infer dates from
   cohort-wide descriptions or nominal treatment cycles.
4. Pre-specify the time window, pair selection, treatment-aware minimal model,
   covariate budget and patient-level split before inspecting state outcomes.
   Keep all samples from the same patient in one partition. Require at least the
   original development/test floor and sufficient support in both treatment
   strata; test adequacy with simulation-based recovery rather than assuming
   a clinical sample total solves identifiability. Do not assume shared drug effects or claim
   causal comparisons between nonrandomized regimens.
5. Produce a new freeze with patient eligibility, estimand and split manifest;
   implement/test the lockbox access guard; then resume the original task ladder.
   Spatial work and Jev remain behind predictive gates.

Additional independent early paired tumour data are needed unless verified
missing libraries change these counts. A full-TME lineage reassessment of
SPECIAL would require its own fixed annotation/QC rules before recounting; it
must not tune labels to recover a desired sample size. Widening time or promoting
the ≥15 sensitivity threshold solely to clear the floor is not justified by this
audit. The broader question remains feasible to investigate only after a new
data/design assessment; this result does not establish that no other dataset exists.

## Exposure, reproducibility and licence notes

No expression-derived changes were computed, no model was fitted, and no
lockbox was assigned. An initial generic workbook inventory inadvertently
printed published state-composition/statistics rows from Sade-Feldman's S1.
That exposure is logged: Sade-Feldman may support development or labelled
transport checks, but **not an untouched final-test claim**. The reproducible
audit now whitelists only its clinical sheet and columns B/E/G/H/I.

Zenodo declares CC BY 4.0 for Lin's deposit. This does not settle the licence of
each upstream cohort inside the integrated object; SPECIAL-derived content
keeps its existing restrictions. Article supplements and GEO records retain
their own terms, with unresolved terms explicitly recorded. No data are
relicensed as MIT or released to a commercial product.

The 4 GB source expands to a 28,600,818,752-byte R serialization stream. A bounded
decoder verified the complete source's internal XXH3 checksum and output SHA-256;
an R-serialization parser retained metadata while explicitly discarding oversized
numeric arrays. It rejected any truncation within the metadata subtree. The
decoder matched the official R converter byte-for-byte on shuffled/unshuffled
synthetic fixtures and rejected corruption; metadata extraction matched R on
native/XDR Seurat fixtures including factors and missing logical values. This
validates the extraction route, **not an actual-source expression-matrix round
trip**. Matrix-to-metadata alignment remains unverified and no expression matrix
was exported for Lin. Reader extensions are isolated from the original freeze.

Reproduce the clinical-table extraction with the bundled document Python and
`scripts/extract_lin_clinical.py`, then run
`.venv/bin/python scripts/audit_cohort_expansion.py`. The extraction record
contains its PDF checksum, table coordinates and pdfplumber version. The audit
contains input hashes. See `docs/lin-metadata-reproduction.md` for the object
audit and validation commands; `reports/EVIDENCE_SHA256.json` indexes saved outputs.

## Primary sources

- Lin et al.: [article DOI](https://doi.org/10.1016/j.ccell.2026.08.007),
  [S1 PDF, page 11](https://ars.els-cdn.com/content/image/1-s2.0-S1535610826003843-mmc1.pdf),
  [S2A workbook](https://ars.els-cdn.com/content/image/1-s2.0-S1535610826003843-mmc2.xlsx).
- Sade-Feldman et al.: [GEO GSE120575](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575),
  [clinical workbook](https://ars.els-cdn.com/content/image/1-s2.0-S0092867418313941-mmc1.xlsx).
- [GSE272993](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272993),
  [GSE164237](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164237),
  [GSE123813](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813).
- [Sun et al.](https://www.thno.org/v14p2127.htm),
  [Wei et al.](https://doi.org/10.1016/j.xcrm.2026.102716),
  [Liu et al. preprint](https://doi.org/10.64898/2026.09.12.751186).
