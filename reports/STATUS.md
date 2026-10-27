# Execution status — 2026-10-27

The planned prediction study is stopped at a necessary feasibility condition.
Only three monotherapy patients meet the primary total-CD8 threshold in both
biopsies. The lower pre-specified threshold yields six, still below ten
development patients even without reserving a lockbox. These are upper bounds
because individual biopsy intervals are missing.

| Task | Status | Output |
|---|---|---|
| T0.1 | Complete | Isolated project; pinned and reproduced Python/R environments |
| T0.2 | Complete | Immune RDS and both clinical S1 workbooks; integrity manifest |
| T0.3 | Complete | Data-use policy and source-specific licence notes |
| T0.4 | Complete | Local plan/config freeze with hashes and Git source commit |
| T1.1 | Complete | Per-sample neutral export; exact count/metadata round trips |
| T1.2 | Complete with missing fields explicit | Full clinical sample table and reconciliation |
| T1.3 | Provisional | Count-only eligibility screening and preliminary flow figure |
| T1.4–T1.5 | Not started | Lockbox/state analyses blocked by feasibility |
| T1.6 | Partial | Site, chemistry and preservation audit; score analyses not run |
| T1.7 / G1 | Not completed | Necessary sample-size precheck fails |
| Later phases | Not started | Depend on a revised, feasible design |

Read `T1.7_preliminary_feasibility.md` for the finding and revision options.
The user approved including both nivolumab and pembrolizumab; no other scientific
threshold was relaxed. No patients were assigned to a lockbox, no SPECIAL state
changes were inspected, and no model/space/Jev work was started.

## Cohort expansion selected by the user

The clinical-source review is complete; read `cohort_expansion_review.md` and
the A02 draft in `prereg/`. Lin2026 adds at most six early monotherapy pairs,
giving an optimistic bound of nine with SPECIAL, still below the minimum.
These six have unverified paired scRNA presence, total CD8 and exact drug names.
The proposed Mono/Combo expansion has an optimistic bound of 23 candidates;
it requires an explicit new treatment-scope decision and further auditing.

Next action: resolve A02.2, then verify actual candidate data and eligibility.
The current freeze remains valid and the original primary scope remains
monotherapy. No gate has been passed by the expanded bounds.

No state changes were computed. Published Sade-Feldman state summaries were
inadvertently printed during workbook inventory; this exposure is recorded in
the review and rules out an untouched-final-test claim for that cohort.
