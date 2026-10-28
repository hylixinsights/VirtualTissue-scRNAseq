# Group testing: melanoma from dissociated scRNA-seq

Run the README's real-melanoma workflow and its recognition-off control, then
repeat with your own reviewed biopsy. The synthetic melanoma fixture is useful
for offline tests and debugging; identify it explicitly in reports.

1. Confirm source checksums, sample identity and full/sampled lineage composition.
2. Open the viewer. Find interaction should jump to a real recorded attack and show its physical distance, allowed actions and target.
3. Click tumor, immune and unresolved cells. Inspect initial scores separately from simulated damage, fatigue and cycle state.
4. Play, pause, scrub backwards, zoom, drag, reset and toggle the field. Counts must match the selected frame.
5. The recognition-off control must have zero attacks and attack-mediated deaths; frozen must preserve all positions and counts.
6. Reopen the saved ZIP with `vtsc explore` and verify checksums. Same input/seed/config should reproduce the same trajectory.
7. Try `--layout mixed` and another seed. Geometry sensitivity is expected because positions are assumed.
8. Test your own H5AD/Cell Ranger H5/MTX and annotated metadata. Sample mixing and malformed annotations should produce clear errors.

Report the software commit/version, OS/Python, command, seed, input type, selected
sample, scenario/layout/config, error text and minimal reproducible example.
Attach `summary.json`, relevant decisions and a screenshot when appropriate.
Do not attach private expression data or identifying metadata to public issues.

A clean test run establishes software behavior, not biological truth. Please
review the fixed marker recipe, cell identity mapping and declared rules as
scientific assumptions separately from bugs in executing those rules.
