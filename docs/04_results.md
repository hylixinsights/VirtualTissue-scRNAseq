# Inspect, share and compare results

```bash
vtsc explore outputs/my_run --open
vtsc explore outputs/my_run.zip --html outputs/reopened.html --open
vtsc compare outputs/my_run outputs/my_control outputs/comparison
```

`explore` verifies recorded file hashes before reopening the run. A changed
event table, model or figure fails verification. A checksum establishes file
integrity relative to the manifest, not independent authenticity.

Open `explorer.html` in any modern browser. No server, internet connection or API
key is needed. Move the slider or press Play; click a live agent to inspect its
recorded history, including daughters. The visualization displays up to 2,000
live cells in an arbitrary grid; a position is never a tissue coordinate.
Cells that have died can be inspected in `events.tsv`; the current live-cell grid
does not provide a deceased-cell picker.

`trajectory.tsv` contains exact process states sampled at output times. Events
between samples may cancel; `events.tsv` contains the full-precision sequence.
`final_cells.tsv` includes live measured/supplied IDs and new `sim_...` IDs;
new identifiers are checked against all previous IDs to avoid collisions.
No input expression matrix is copied into the run.

The results ZIP contains the same files as the run directory. Share it only
when its IDs, state labels and source terms permit sharing. A ZIP remains
subject to the input's licence and your own data-governance rules.

Comparisons require the **same input bundle hash and state vocabulary**. The
plot ends at the shorter observed horizon. It is a descriptive comparison of
software assumptions; there is no automatic significance test or treatment-effect
estimate. Use multiple independent seeds for uncertainty assessment and report
any capped runs.

An incomplete run is not necessarily a software bug: the maximum population
or event count may have been reached. Inspect the summary before raising caps.
Large horizons and explosive birth rates can exceed memory/disk or overflow the
analytical expectation. Reduce them rather than interpreting numerical overflow
as biological saturation.
