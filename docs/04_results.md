# Inspect a tissue run

Open `explorer.html` directly; no server or network is required. Find interaction
jumps to the next saved contact attack and zooms onto the actor. Cell selection
shows the barcode, current fate, latest allowed actions, decision, target, distance,
local cue, fatigue, damage, cycle and fixed initial RNA program scores. The timeline
is a replay of recorded frames. Position interpolation between frames is for display;
actions, counts and decisions use the selected recorded frame.

The signal overlay is a generic simulated cue, not measured expression. Cell
membranes and nuclei are schematic, and their visible shapes are not segmentations.
Drag to pan, use +/− to zoom and Fit to reset. A recognition-off run should report
no recorded attacks. At time zero there are no decisions yet.

`vtsc explore RUN_OR_ZIP` verifies the manifest. For a ZIP, add
`--html outputs/reopened.html --open` to regenerate a viewer. No simulation is rerun.

`run.json` contains the immutable lineage catalog, compact frame arrays (column
names in `frame_columns`), decision records and events. Agent indices refer to
`catalog`. Only present agents are in a frame; cleared cells remain in the catalog.
Daughters carry a parent index and birth time. Cell types and RNA programs are fixed.
`decisions.jsonl` and `events.jsonl` duplicate those records for line-oriented analysis.

Compare matched `trajectory.tsv` and `summary.json` files. Check initial bundle
hash, seed, scenario, layout, parameters and observed duration before interpreting
a difference. The legacy `vtsc compare` is for the older nonspatial engine and
rejects tissue recordings. This release offers descriptive comparisons, not a
calibrated statistical comparison or fitted treatment contrast.
