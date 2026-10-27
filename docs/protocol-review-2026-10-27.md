# Protocol review before state-outcome analysis

The original plan is preserved byte-for-byte in
`prereg/plan-v1-original-2026-10-27.md`. The working plan remains in the parent
`scRNAseq` directory. No state compositions or changes have been examined.

## Adopted amendment

The user explicitly selected both anti-PD-1 monotherapies on 2026-10-27.
The primary treatment definition is therefore nivolumab **or** pembrolizumab,
with combination therapy separate. `config/study.yaml` records this definition.
Both clinical tables agree on all shared clinical fields. The earlier
nivolumab-only wording was inconsistent with the published cohort.

## Timing and lesion identity

The final [2025 article](https://pmc.ncbi.nlm.nih.gov/articles/PMC12409281/)
describes early sampling as ±2–3 weeks; the supplied plan cites different
preprint wording. Neither clinical Table S1 contains biopsy dates, elapsed
days or patient-specific interval bounds. The immune object's only time
fields are categorical. Do not infer a patient's interval from these labels,
the drug name, or a dosing schedule. A nominal-time analysis would require
a documented amendment and would not recover measured rates.

The clinical tables describe tissue and anatomical site, not unique lesion
identifiers. Matching site text does not establish that two biopsies came
from the same lesion. Preserve an unknown lesion-match status.

## Lockbox and inferential resolution

`reports/evidence/lockbox-test-resolution.json` reproduces the design calculation:
with four nonzero paired differences, the smallest one-sided exact Wilcoxon
p-value is 1/16 = 0.0625. Zeros cannot improve that resolution. The lockbox
cannot confirm a test at alpha 0.05, even if every prediction improves.
An asymptotic approximation must not be substituted to evade this limit.
Estimation remains possible, but full LOPO would remain development-informed.
The lockbox size is unchanged; no patients have been assigned.

## Issues to resolve before any resumed development

- The count-only screening fails the necessary sample-size condition for G1.
  See `reports/T1.7_preliminary_feasibility.md`. State coarsening cannot create
  additional CD8 cells; lowering the threshold to the pre-specified alternative
  still does not meet the required development sample size.
- Multiple baseline and on-treatment biopsies for GC16 require an explicit
  pairing or aggregation rule before any future state-outcome inspection.
  Pooling different lesions must not be done silently.
- Freeze the feature/graph construction recipe, and refit learned graphs and
  representations inside every training fold. A graph fitted on all development
  patients cannot become a fixed input to their own held-out evaluations.
- The G4 alternative “informative” is not a numerical decision rule. If the
  study is redesigned and resumed, define that criterion before state outcomes
  are used for model selection or evaluation.

The local freeze records the current design and its limitations; it is not
a public registration or a statement that the study is feasible. Future
amendments need new versioned records rather than overwriting this freeze.
