# Results, controls and inspection

A v2 run saves compact full records, JSONL decisions/events, trajectory TSV,
source provenance, model contract, summaries, SHA-256 manifest, offline HTML and
ZIP. The journal is hash chained to detect edits, not digitally signed identity.
`vtsc explore` verifies both archive hashes and journal before regenerating HTML.
Standalone HTML can be opened directly but does not authenticate its own source;
verify the ZIP before relying on a shared run.

Use `compare-melanoma A B OUT` for matched input, seed, layout and base contract.
Only arm mechanisms differ. The viewer replays a shared timeline and the most
recent saved frame at/before each time. Drug exposure is spatially uniform;
the generic cue overlay is a different field. Cell outlines/nuclei are schematic.
No contact or transition is invented by the renderer.

Select cells by clicking or keyboard identity selector. Find transition/attack
jumps to a recorded event. Read the prior-state fate signal separately from the
post-transition action signal. Both are recorded; they can differ when the state
changes its PD-1 activity. State age, four capacities, fatigue, permitted events,
hazards, probabilities and draw are inspectable. Full journal expands all fields.
RNA t0 and initial scores remain separate from assumed current capacity profiles.

Controls:

| Arm | Change |
|---|---|
| anti-pd1 | Configured normalized schedule |
| exposure-zero | Same schedule boundaries with zero exposure; simulated counterfactual |
| recognition-off | Reactivity/antigen recognition disabled; exposure retained |
| checkpoint-off | PD-1/PD-L1 term disabled; drug decoupled |
| transitions-off | CD8 states fixed; other mechanics/demography retained |
| spatial-shuffle | Identities permuted within equal-radius position groups |
| frozen | All evolution disabled |

Use ensembles, not a favorable single seed. The supplied sensitivity varies
reactive fraction, checkpoint strength and all transition rates by factors 0.5/2,
with the original value as baseline; fractions saturate at 1. Saved quantiles
represent Monte Carlo variation over the specified seeds, not clinical intervals.
The recording/population cap stops a run explicitly. Always inspect `complete`,
`stop_reason` and `observed_hours`; a capped plateau is not tumor control.

`recover-melanoma RUN NEW_OUTPUT [--hours H]` recomputes from the embedded input
and original seed/configuration, checks the prefix, then optionally extends.
It needs the same implementation and numerical environment; mismatches are
errors. Existing outputs are never overwritten. Use full time-step horizons
when extending; a different final partial step may make prefix equality fail.
