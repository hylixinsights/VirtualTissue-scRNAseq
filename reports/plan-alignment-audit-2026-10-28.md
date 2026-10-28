# Original-plan alignment audit — 2026-10-28

## Verdict

**The v0.3 melanoma demonstration is not an implementation or completion of the
original SPECIAL predictive plan.** It is an independently usable, uncalibrated
software prototype. The change from PBMC to dissociated melanoma corrected the
example's tissue context, but did not implement the planned CD8 state model,
anti-PD-1 mechanism, parameter estimation or patient-level validation.

The research audit did follow the original ordering through the available
feasibility work and stopped before lockbox assignment. That stop remains valid
for the recorded design. Implementing a separate demonstration does not pass G1,
G4 or G5 and does not satisfy T3.1 or T5.3.

A03/A04 record the assistant's implementation choices in response to requests
for a testable pipeline and a melanoma viewer. They must not be read as user
approval to replace the scientific question, waive predictive gates or substitute
handwritten rules for inferred CD8 state dynamics. The previous delivery should
have made this distinction more explicit.

## Inputs inspected

- The original v1 plan and the v1.1 frozen plan, plus `config/study.yaml`.
- The user's pasted plan: byte-for-byte identical to the v1.2 working copy at
  audit start. It includes execution updates; it is not the pristine v1 file.
- The pasted cohort expansion review and its locally saved evidence.
- Current task board, source manifests, A03/A04, implementation modules,
  simulation validation report and README at source commit `ade322f`.

Source hashes, board counts and a task-by-task table are saved in the companion
`plan-alignment-audit-2026-10-28.json`. This review does not calculate new state
compositions, inspect a test set, change eligibility or refit a model.

## Requirement versus implementation

| Original requirement | What v0.3 actually does | Assessment |
|---|---|---|
| Predict early on-treatment **CD8 state composition** from the same patient's baseline biopsy, in unseen patients | Produces local actions, motion, damage, growth and death from one snapshot | Different output; primary question unanswered |
| SPECIAL paired biopsies; nivolumab/pembrolizumab; approved expansion audit keeps combination separate | Public 10x dissociated melanoma from one donor, no paired treatment endpoint | Appropriate software example, not a replacement predictive cohort |
| Memory, Activated, Cytotoxic and Dysfunctional CD8 states; cycling as a separate axis | Broad `CD8_T` lineage is fixed; no transitions among those four states | Central biological representation missing |
| Versioned programs, independent functional axes, AUCell and training-only calibration by depth | Small fixed marker panels; mean log-normalized expression transformed to [0,1] | Illustrative scoring, not T2.2/T2.3 |
| Robust state graph and directed/undirected transitions supported by patient resampling | No CD8 state graph in the tissue engine | T2.4 missing |
| Four axes: dysfunction, proliferation, cytotoxicity and persistence | Initial cytotoxicity, proliferation and interferon scores; fixed over time; interferon does not drive the engine | Does not satisfy T2.5 or planned dynamic agent state |
| M0/M1/M1b/M2/M3, latent baseline composition, Dirichlet-multinomial likelihood and estimated parameters | Legacy Markov simulator with supplied rates; tissue simulator with supplied local rules | Neither implements the fitted model ladder |
| Parameter recovery, identifiability, LOPO, lockbox and external direction validation | Software invariants, negative controls, seeded replay and notebook execution | Useful engineering verification; not predictive validation |
| Spatial hazards inherited from frozen nonspatial model; reproduce M3 at sample-mean local inputs | Separate attack/growth/migration rules; no M3 coupling | T5.3/T5.4 not satisfied; the old Markov expectation test does not count |
| Agent state/time-in-state, independent clone/patient/avatar fields, checkpoint term and anti-PD-1 schedule | Lineage/barcode/programs, fatigue/damage/cycle, parent ID; no checkpoint drug schedule | Partial reusable mechanics, not planned runtime |
| Stratified proportional baseline virtual biopsies; site templates and shuffled spatial controls | Uniform subset; assumed mixed/two-nest layouts; no held-out patient evaluation | Does not satisfy T5.5/G5 |
| Jev chooses allowed **state transitions**, compared with frozen hazards and controls | Deterministic priority chooses actions such as attack, migrate, grow or rest | No Jev implementation; no state-fate policy comparison |
| Full ODD protocol, S1–S13 and claim-to-result ledger | Technical model description and software validation; some early research evidence | Supplementary blueprint not assembled |
| Reproducibility, source terms, hashes, event replay and untouched reference engine | Many of these engineering components exist | Reusable work, not evidence that later scientific tasks are complete |

The critical code distinction is visible in `src/vtsc/tissue.py`: `decision_menu`
constructs actions; `simulate_tissue` chooses actions. The `catalog` retains a
fixed `cell_type` and initial `programs`. Fatigue or switching from migration to
attack is not a Memory→Activated→Cytotoxic/Dysfunctional state transition.
The inherited CD8 marker label is not a simulated evolving CD8 state.

`recognition-off` is a contact-killing control. It is **not** an untreated
anti-PD-1 arm or a validated no-treatment counterfactual. The prototype's
simulated deaths are not evidence of a patient's treatment effect.

## What has actually been completed in the study

The board has **46 entries: 6 recorded complete, 2 partial, 38 pending**. The
entries include explicit gates. These are execution statuses, not a percent
completion weighted by task difficulty.

Recorded complete: T0.1–T0.4 and T1.1–T1.2. T0.2 covered priority-1 acquisition;
T1.1 covered the immune object. This does not mean the whole-TME and malignant
objects were acquired/exported. Both are still marked missing in the manifest
and absent at the configured paths checked for this audit.

Partial: T1.3 count-based eligibility and T1.6 technical-confounder metadata.
The full G1 report, split/access guard, state reference, fitted models, final
validation and later scientific phases remain incomplete. No `LOCKBOX.json`
was assigned by the research workflow. Missing dates/site/chemistry information
remains explicit even where the clinical inventory task is recorded complete.

The existing 34 software tests and three executed notebooks demonstrate selected
implementation properties. They do not evaluate H1–H4 or complete the original
clean-environment reproduction of all paper figures/tables (T8.4).

## Feasibility evidence: use the latest count audit

The pasted SPECIAL-only result is consistent with the original count precheck:
3 monotherapy candidates at ≥30 CD8 per biopsy, 6 at ≥15, and 3 at ≥50.
The later deposited-object audit supersedes optimistic clinical-only totals:

| ≥30 CD8 in both biopsies | SPECIAL | New Lin patients | Total |
|---|---:|---:|---:|
| Anti-PD-1 monotherapy | 3 | 1 | 4 |
| Anti-PD-1 + anti-CTLA-4 | 4 | 4 | 8 |

These are upper bounds before unresolved eligibility checks. The total of 12
is not permission to pool treatment strata. Even that total would leave only
8 after reserving 4, below the existing 10-development floor; the separate-stratum
design additionally needs adequate support within strata. The earlier clinical
ceiling of 23 was not the number of usable deposited pairs.

Public availability does not fill absent libraries, identify an ambiguous
patient/library join or create inter-biopsy dates. Conversely, those omissions
do not prevent using the available public cells for a clearly labelled pilot
or engineering example. These are different claims.

A four-patient exact one-sided signed-rank test cannot reach p<0.05 (minimum
1/16 = 0.0625 without ties/zeros). That limitation was already recorded in the
protocol review; it is a design issue, not a software defect or access issue.

## Public-data-only route, without contacting authors

The user explicitly declined author email on 2026-10-28. **No author outreach,
data request, permission request or waiting for an author response is a project
dependency.** Existing draft correspondence is historical and inactive. No email
has been sent. Available public files are used under their recorded source terms.

Proceed using published articles, supplements, deposited files/versions and
metadata. If an attribute cannot be established there, preserve it as unknown
and apply the declared eligibility/sensitivity rule. Do not invent dates or count
an absent library as zero cells, and do not count reused SPECIAL patients as new.

The next protocol-aligned public-data work is:

1. Record a fixed, state-outcome-blind lineage/QC audit specification, then obtain
   and verify the already identified public `Entire_TME.rds` and
   `Malignant_cells.rds` files. Determine whether the full TME supplies additional
   valid CD8 observations and malignant/environment information. Recovery must
   be demonstrated; labels must not be tuned to clear G1.
2. Continue Lin expression-to-metadata/library alignment and source de-duplication
   using public files only. Unresolved matches remain unresolved/excluded as
   specified; an author reply is not required to finish the public-source audit.
3. Search additional independent, public **paired melanoma tumor** cohorts if
   needed. Blood-only, other-cancer and single-timepoint data do not repair the
   primary paired melanoma endpoint unchanged.
4. Recompute eligibility by cohort and treatment. If the frozen design still
   cannot be supported, present a concrete pilot/descriptive protocol revision
   with reduced claims before fitting a different scientific model. No silent
   threshold, time-window or endpoint changes.
5. Once the research gate is supportable, resume T1.4/T1.5/T1.7 and Phase 2, then
   the fitted nonspatial ladder. Scientific spatial/Jev work remains gated.

Steps above are the next work queue, not claims that those downloads/audits were
performed in this alignment review. The 10x prototype remains available for
software testing and can later reuse a properly fitted state model. It must
not be presented as the completion of the original plan.
