# Group testing v0.4

1. Install in a new Python 3.12/3.13 environment. Run pytest and verify-freeze.
2. Run the 120-cell offline demo. Open comparison/explorer.html directly.
3. Find a CD8 transition. Verify state, prior/new capacities, fate/action signals,
   probabilities, draw and event time against decisions.jsonl and frames.
4. Play, pause, scrub; switch arm/color/field; select a cell by keyboard. Confirm
   both arms start at the same positions and exposure zero is explicitly labelled.
5. Download and prepare the 10x melanoma. Check Unknown/coverage and full versus
   sampled composition. Do not equate marker labels with validated malignancy.
6. Prepare a reviewed own-data biopsy. Try a normalized layer or mixed sample
   without selection: require a clear error. Preserve source permissions.
7. Reopen ZIP, recover into a new output and compare frames. Run multiple seeds
   and controls; report the seed/contract even when no transition or drug effect occurs.

Report software version, OS/Python, input format, command, complete/stop_reason,
minimal reproducible synthetic input and sanitized error log. Do not attach
private patient data to public issues. No remote CI/Colab success is implied by
local validation. Clinical interpretation requires the separate study gates.
