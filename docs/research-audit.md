# Separate SPECIAL/Lin research workflow

The new software-test purpose does not amend the frozen clinical eligibility
thresholds or grant a prediction/spatial gate pass. The user explicitly requested
an example-driven, bring-your-own-data prototype for group testing; see
[A03](../prereg/amendment-A03-software-testing-2026-10-27.md).

The research audit currently has at most 12 count candidates at ≥30 CD8 in both
biopsies, split into 4 monotherapy and 8 combination candidates. It remains
below the existing 10 development + 4 held-out floor before unresolved clinical
exclusions. No clinical model fitting or held-out evaluation was performed.

Read [STATUS](../reports/STATUS.md), [cohort expansion review](../reports/cohort_expansion_review.md),
[A02](../prereg/amendment-A02-cohort-expansion-draft-2026-10-27.md), and the
[working plan](../prereg/working-plan-v1.2-2026-10-27.md). `prereg/FREEZE.json`
still verifies the original frozen files. Never overwrite it to accommodate
the demonstration.

The historical Python environment is `requirements.lock.txt` and
`environment.yml`; R packages are in `renv.lock`, with an isolated reader
extension in `renv-expansion.lock`. Keep research and software environments
separate when reproducing the recorded audit. The software's `[research]` extra
adds workbook/PDF readers, but does not replace the full historical environment.

```bash
vtsc verify-freeze
vtsc status
```

Large source files remain outside Git. The Lin metadata-only fallback and its
limits are documented in [lin-metadata-reproduction.md](lin-metadata-reproduction.md).
The [data-request draft](data-request-draft.txt) has not been sent. No new SPECIAL
or Lin state compositions were opened to build the software example.
