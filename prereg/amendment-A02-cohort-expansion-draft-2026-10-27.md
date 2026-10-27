# A02 — Cohort expansion after feasibility failure

Status: clinical-source audit completed; **treatment/estimand change proposed,
not activated**. The user authorized revising/expanding the primary cohort.
The currently approved primary drugs remain nivolumab and pembrolizumab in
monotherapy. The frozen files referenced by `FREEZE.json` are unchanged.

## Trigger and completed work

SPECIAL supplied only three monotherapy patients above the original total-CD8
threshold in both biopsies. The expanded clinical audit found at most six new
early monotherapy pairs in Lin2026, still insufficient in aggregate. Exact
source rows, bounds and unresolved fields are in
`reports/evidence/cohort-expansion/audit.json` and
`reports/cohort_expansion_review.md`. No eligible-patient manifest or split has
been frozen for an amended study.

## Proposed route A: treatment-stratified early melanoma ICI

Audit SPECIAL and Lin's in-house patients as independent source cohorts. Keep
nivolumab/pembrolizumab monotherapy and anti-PD-1 + anti-CTLA-4 combination as
separate treatment strata, subject to exact drug/histology verification.
Do not include blood, BCC/SCC, toripalimab or reused atlas patients automatically.
The clinical/count-only upper bound is 23 patients, not a verified analysis n.

Proposed endpoint: early on-treatment tumour CD8 state composition conditional
on baseline composition/environment, known treatment regimen, and measured
inter-biopsy interval. Treatment and cohort must remain explicit in the minimal
model and baselines. Report monotherapy and combination performance separately.
Do not interpret a predictive treatment coefficient as a causal treatment effect.
The exact early time window and sample-selection rule must be justified and
frozen before state outcomes are opened; use all clinical samples for the audit.

Retain ≥30 total CD8 per biopsy, 15/50 sensitivity thresholds, patient as the
unit, no response-based eligibility, and no fitting from jointly integrated
embeddings. Keep at least 10 development patients and 4 held-out patients as
the existing numerical floor, then assess treatment-stratum coverage and model
identifiability. Passing the floor alone is not adequate statistical power.
A small final holdout remains estimation with uncertainty, not an exact
one-sided p<0.05 confirmation with four patients.

If eligible patients/strata cannot support this revised model, stop again; do
not add parameters or remove the lockbox to manufacture a pass. Whole-cohort
holdout is desirable only if enough development patients remain and treatment
is not confounded with cohort. Otherwise pre-specify a patient-level holdout
and label any leave-one-cohort-out evaluation as an additional transport test.

## Alternative route B: preserve the original monotherapy estimand

Acquire additional independent early paired tumour monotherapy data or audit
the complete SPECIAL TME with a fixed broad-lineage rule and technical QC.
The deposited immune counts remain the current evidence until this separate
audit demonstrates more eligible cells. Do not tune annotation/QC to recover
a desired sample size. Additional access requests or author correspondence
would require explicit instructions before sending messages.

## Why late-interval pooling is not selected by default

Sade-Feldman contributes five PD1-labelled clinical pairs, with 35–199 days
between first biopsies. Adding these would change the time target and assay
context while only reaching an optimistic 14 total candidates. CD8 totals,
drug names and some covariates remain unverified. Published state-summary
exposure during file inventory also rules out claiming this cohort as untouched
final validation. Keep it for labelled development/transport uses if appropriate.

## Task dependencies and supplementary destinations

| ID | Task / dependency | Deliverable and acceptance | Supplement |
|---|---|---|---|
| A02.1 | Clinical source screening; after failed feasibility precheck | Complete: candidate register, clinical rows, hashed bounds | S1–S2; S1 flow update after eligibility |
| A02.2 | Resolve new treatment scope; after A02.1 | Explicit user decision recorded; otherwise original scope remains | S1 protocol amendment |
| A02.3 | Lin in-house identity, raw-count availability, CD8 totals and overlap; after A02.2 for route A | Verified source object/subset and checksum; cell/sample alignment; cohort-qualified patient IDs; counts only | S1–S3; sample dictionary |
| A02.4 | Drug, timing, lesion, chemistry, histology and licence reconciliation; after A02.3 | Required fields verified or exclusion reasons; no guessed intervals | S1–S3; eligibility and licence tables |
| A02.5 | Freeze amended estimand and eligibility before states; after A02.4 | Feasible per-arm patient allocation, model-size justification, deterministic pair/split policy and access-guard tests | S1, statistical methods and claims ledger |
| T1.4 onward | Resume the original ladder; after A02.5 | Original acceptance criteria with explicit amendment mappings | Original S1–S13 plan |

The original FREEZE remains an immutable record. A02 must receive its own
freeze/version only after the scientific scope, actual eligibility and split
policy are finalized. No space or Jev gate is waived.
