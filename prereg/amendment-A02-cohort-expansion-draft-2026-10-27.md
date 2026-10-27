# A02 — Cohort expansion after feasibility failure

Status: **user approved auditing monotherapy and combination as separate strata**
on 2026-10-27. Clinical-source screening and source-metadata CD8 counting are
complete. The necessary count prerequisite fails again: 12 candidates at ≥30,
comprising 4 Mono and 8 Combo, before unresolved exclusions. Raw-matrix alignment
and clinical reconciliation remain incomplete. The final analysis estimand,
eligibility and split remain unfrozen.
The original monotherapy freeze referenced by `FREEZE.json` is unchanged;
`config/cohort-expansion-audit.yaml` records the separately approved audit scope.

## Trigger and completed work

SPECIAL supplied only three monotherapy patients above the original total-CD8
threshold in both biopsies. The expanded clinical audit found at most six new
early monotherapy pairs in Lin2026, still insufficient in aggregate. Exact
source rows, bounds and unresolved fields are in
`reports/evidence/cohort-expansion/audit.json` and
`reports/cohort_expansion_review.md`. No eligible-patient manifest or split has
been frozen for an amended study.

## Selected audit route A: treatment-stratified early melanoma ICI

Audit SPECIAL and Lin's in-house patients as independent source cohorts. Keep
nivolumab/pembrolizumab monotherapy and anti-PD-1 + anti-CTLA-4 combination as
separate treatment strata, subject to exact drug/histology verification.
Do not include blood, BCC/SCC, toripalimab or reused atlas patients automatically.
The initial clinical bound of 23 is superseded by the deposited-object result:
only 8 Lin patients have both Pre/On labels, and 5 meet ≥30 CD8 in both. Together
with 7 SPECIAL candidates this gives 12 (4 Mono, 8 Combo), not a verified analysis
n. Three lineage-label counting rules give the same qualifying patient sets;
the specific Lin rule includes explicitly labelled cycling CD8 and excludes gdT.
At ≥15 there are 17 candidates; this sensitivity result does not replace ≥30.
See `reports/evidence/cohort-expansion/lin-object-count-audit.json`.

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
| A02.2 | Resolve new treatment scope; after A02.1 | Complete: user selected monotherapy and combination in separate strata | S1 protocol amendment |
| A02.3 | Lin in-house identity, raw-count availability, CD8 totals and overlap; after A02.2 for route A | Partial: source checksum, metadata counts and explicit SPECIAL reuse audited; raw-count alignment and exact library mapping unresolved | S1–S3; sample dictionary |
| A02.4 | Drug, timing, lesion, chemistry, histology and licence reconciliation; after A02.3 | Partial: public clinical/EGA sources audited; missing library crosswalk, drug/interval/lesion/histology fields and chemistry conflict documented | S1–S3; eligibility and licence tables |
| A02.5 | Freeze amended estimand and eligibility before states; after A02.4 | Not started: numerical prerequisite fails; new source evidence required before model-size/split/access-guard work | S1, statistical methods and claims ledger |
| T1.4 onward | Resume the original ladder; after A02.5 | Original acceptance criteria with explicit amendment mappings | Original S1–S13 plan |

The original FREEZE remains an immutable record. A02 must receive its own
freeze/version only after the scientific scope, actual eligibility and split
policy are finalized. No space or Jev gate is waived.
