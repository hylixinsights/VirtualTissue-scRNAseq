# Melanoma tissue model and assumptions (v0.3)

This is a local rule agent prototype for software tests. It does not fit kinetics
or claim that dissociated RNA identifies spatial organization or treatment effects.
The previous exact nonspatial Gillespie engine remains in `model.py` for legacy
regression tests; the new tissue engine is `tissue.py`.

## Initial conditions

Each retained cell has its original barcode, a lineage label and three initial
RNA programs: cytotoxicity, proliferation and interferon. Fixed gene panels are
in `tissue_inputs.py`. Counts are normalized to 10,000 per cell and log1p-transformed.
For each panel, average over the full panel length (missing genes count as zero),
then map the mean `s` to `1-exp(-s/2)`. Scores stay fixed during a simulation.
The interferon score is exported for inspection; it does not currently drive a rule.
A lack of expression is not proof of absent protein or biological capability.

The optional marker recipe selects the highest mean lineage score only if it
is ≥0.5, has ≥2 detected panel genes, and exceeds the second score by ≥0.15.
Otherwise the cell is Unknown. A T-cell candidate with CD8A or CD8B detected
becomes CD8_T. These are declared heuristic thresholds, not a validated classifier.

Uniform subsampling uses a specified seed and no lineage enrichment. Cell
composition is reported before/after selection. Dissociation bias is not corrected.

## Geometry and local signal

Agents occupy non-overlapping circular exclusion discs. Radii are assumed:
melanoma 8 µm; CD8/other T/B 4.8; NK 5.3; myeloid 7; fibroblast 6.5;
endothelial/unknown 6. The square domain area equals total disc area divided
by packing fraction (default 0.32), with a minimum width of 100 µm.

`tumor-nests` places melanoma-labelled cells around two assumed centers before
placing other cells. `mixed` uses uniform placement for all. Rejection sampling
prevents initial overlaps. These positions are neither measured nor recovered
from RNA. Cell shapes/nuclei and extracellular fibers in the viewer are schematic;
fibers have no mechanical effect and exclusion uses the discs.

A generic cue is sourced by melanoma-labelled live cells on a 32×32 grid,
diffused by a Gaussian kernel with reflecting boundaries and exponentially
decayed. Its units are arbitrary. It is not identified as a measured cytokine.
Default diffusion is 90 µm²/h, decay 0.3/h and source 0.5 units/h/cell.
Mobile cells follow the local gradient plus seeded exploratory noise, with
speed 10 µm/h. Candidate moves that overlap any present cell or leave the
domain are rejected. There is no compression, force balance or tissue adhesion model.

## Time steps, legal actions and policy

The default time step is 0.25 h. Cells act in a seeded random order each step;
actions update the shared local environment immediately. At every decision,
the kernel constructs a menu. Priority is:

1. Myeloid clearance of a corpse in contact.
2. CD8/NK attack on a live melanoma target in contact, if recognition is explicitly on and fatigue <0.85.
3. Melanoma division if cycle is complete (placement can still fail).
4. Melanoma growth if cycle is incomplete (space can still block it).
5. Migration for mobile, non-fatigued lineages.
6. Rest.

Every menu, choice, reason, target, contact count, local cue and execution outcome
is recorded. This is a deterministic local decision policy with seeded exploratory
movement, not an LLM. Failed growth/movement/division has `blocked_by_space` outcome.

Contact means center distance ≤sum of current radii +2 µm. Attack damage is
`0.65 * (0.25 + 0.75 * initial_cytotoxicity) * dt`. The basal 0.25 is an explicit
scenario assumption, so a zero RNA score does not itself forbid an attack.
Damage reaching 1 causes death; dead agents cannot act or resurrect. Attacking
increases fatigue by 0.32/h, rest reduces it by 0.22/h, migration by half that.

Cycle progression is `(0.2 + 0.8 * initial_proliferation) * dt / 24 hours`.
Radius grows as `base_radius * sqrt(1+cycle)` only when space permits. At cycle 1,
two base-radius discs must fit for division; the parent retains its ID, the daughter
gets a new ID and parent reference, and both reset cycle to zero. Daughters inherit
initial programs and existing damage. Growth and division may not occur in a short
run, even though the corresponding tests exercise them with declared faster cycles.

Dead bodies persist until myeloid contact clearance or the 6 h timeout. Their
remaining exclusion disc blocks movement. The population cap counts all created identities (including dead cells) and
reports `population_cap`. A 100,000 recorded-cell budget limits portable file size
and reports `recording_cap` when reached; observed time must be checked rather than assumed complete.

## Controls and interpretation

- `recognition-off` (default): same agents, movement and growth; no attacks.
- `declared-recognition`: all eligible cytotoxic agents may recognize all melanoma-labelled targets; a deliberately broad scenario assumption.
- `frozen`: no evolution; positions/population are identical at all frames.

Matching seeds produces matching initial geometry between scenarios; subsequent
random draws and trajectories can diverge. Compare multiple seeds, both layouts
and time-step/parameter sensitivity. A single visual trajectory does not establish
a treatment effect, antigen specificity, recovered histology or clinical prediction.
There is no receptor-level communication model, vascular/oxygen/metabolic field,
anti-PD-1 mechanism, Jev adapter, fitted policy or external validation in this release.
