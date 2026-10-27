# Execution status — 2026-10-27

The planned prediction study remains stopped at a necessary feasibility
condition after the user-approved separate-stratum expansion audit.
SPECIAL + new Lin patients supply at most **12 count candidates at ≥30 CD8 per
biopsy: 4 monotherapy and 8 combination**. Four held out would leave only eight
for development, below the existing minimum of ten. The ≥15 sensitivity gives
17 candidates and ≥50 gives 11; the primary threshold remains ≥30.
Clinical identity, exact drugs, intervals and other checks can further reduce
these counts. No final eligibility or G1 pass is claimed.

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

The clinical review and deposited-object metadata counting are complete; read
`cohort_expansion_review.md` and the A02 draft. Lin has 25 deposited in-house
samples and 8 paired patients, versus 34 clinical rows and 16 clinical pairs.
Only five deposited pairs meet ≥30 total CD8. The historical clinical bound of
23 is superseded by the current count result of 12. Explicitly cycling CD8 were
included; alternative lineage rules give the same candidate sets.

Next dependency: resolve the missing libraries and patient/library crosswalk,
P2's On day, exact drug/histology/lesion fields, chemistry discrepancy and
SPECIAL intervals, or obtain additional qualifying data. An author data request
is drafted in `docs/data-request-draft.txt` and has not been sent.
A02.3 remains partial because raw-matrix alignment and exact library identity
are unverified; A02.4 is partial after public-source reconciliation. A02.5 is
not started because the numerical prerequisite fails.
The original monotherapy freeze remains valid. The approved separate-stratum
audit is recorded in `config/cohort-expansion-audit.yaml`; no amended model or
patient split is frozen yet. No gate has been passed by the expanded bounds.

No state changes were computed. Published Sade-Feldman state summaries were
inadvertently printed during workbook inventory; this exposure is recorded in
the review and rules out an untouched-final-test claim for that cohort.

The original freeze still verifies. Lin's full download matches the deposited
MD5; its internal checksum verifies. Bounded decoding and metadata extraction
passed independent synthetic R comparisons; an actual-source matrix round trip
has not been performed. Reader environments, provenance and reproduction steps
are saved separately from the original analysis freeze.
