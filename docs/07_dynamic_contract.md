# Dynamic melanoma contract / ODD v2 (exploratory)

## Overview and entities

Purpose: reproducible experiments with melanoma and CD8 state dynamics, not
clinical response prediction. Input is one dissociated biopsy. RNA at t0 and
its checksum never change. Agents retain cell_id, sample_id, patient_id (nullable),
lineage, cd8_state, state age, parent_id, clone_id (nullable), molecular_source_id
and simulated reference_id. Daughters inherit state/clone/source, reset age to
zero, and have a simulated parent link. Influx has no observed RNA source.
Unknown CD8 stays Unknown; it is never silently mapped to Memory.

Four CD8 states have independent assumed capacity profiles. Dysfunctional can
proliferate. Fatigue is a reversible mechanical variable, not exhaustion.
The graph permits reversible and branching transitions. Every edge and profile
is an exploratory assumption; no graph was inferred from one sample.

## Shared dynamics

For state j, allowed rates q_jk are nonnegative per hour. Environment modifies
activation/effector edges by antigen signal and dysfunction edges by repeated
stimulation plus local tumor fraction. S = reactive * target * contact *
costimulation * (1 - checkpoint); checkpoint = enabled * strength * PD1 * PDL1 *
(1 - exposure). All receptor/ligand activities are assumed normalized values,
not RNA-derived protein or occupancy. Exposure is a piecewise constant [0,1]
schedule in hours; intervals are half-open. It is not a clinical dose or PK.
No edge directly maps Dysfunctional to Memory. Drug can increase stimulation
and dysfunction together; benefit is not guaranteed.

Competing transitions, birth and death have H=sum(h), stay probability exp(-H dt),
and p(event k)=(1-exp(-H dt))*h_k/H. One event per cell per step (tau-leap);
new daughters act next step. Influx is independent Poisson(rate*dt). The common
continuous-time expectation is dx/dt=(Q^T+diag(b-d))x+influx. The fixed-step
expectation uses the exact same one-event probabilities, including births;
convergence to the continuous-time solver is tested as dt decreases.

M0 persists baseline; M1 adds a training CLR shift; M1b scales it by elapsed time.
M2 uses state-only rates; M3 uses the shared environment-conditioned rates.
Penalized maximum likelihood uses pinned SciPy, Dirichlet-multinomial observation
noise and Monte Carlo marginalization over Dirichlet(n_pre+0.5) baseline draws.
Seeded fixed draws make the objective deterministic. Only identifiable subsets
are fitted; rank-deficient sensitivity/Jacobian results are reported explicitly.
No fitted synthetic model is used as a clinical model.

## Spatial process and observation

Reuse v0.3 non-overlap, bounded movement, cue diffusion, contact damage, tumor
cycle/division and corpse clearance. Cue is a generic assumed chemotactic field;
local tumor/myeloid fractions use live neighbors as denominator. Missing measured
TME features remain null. The drug field is uniform (no measured transport).
CD8 hazards execute before actions; changed capabilities affect that step's
movement, damage and proliferation. Birth can fail due to space and is logged.
Death is irreversible. Frames contain actual states, ages and capacities.
The journal records causal options, rates, probabilities, random draw, context,
prior/new states, provider, outcome and seed with a hash chain.

Matched arms share input, initial positions and keyed random streams. Arms:
anti-pd1, exposure-zero (simulated counterfactual), recognition-off,
checkpoint-off, transitions-off, spatial-shuffle and frozen. Shuffle preserves
positions/composition/radii by permuting identities within radius groups.
Mean-field mode fixes context and disables movement, attacks, clearance and
space constraints on CD8 birth; this is a software consistency experiment.

## Constraints and provenance

Input schema v2; run schema v2. Legacy v1 inputs require explicit migration and
cannot invent CD8 states. Legacy v0.3 run-tissue remains reproducible. Initial
2–2000 agents; max 5000 created; <=96 h, dt 0.05–0.5 h; recording budget and stop
reason explicit. Intended group examples: 100–400 cells, 24 h; runtime/memory
measured in release evidence. No network is required after example download.
No output overwrite; recovery recomputes deterministically into a new directory.
Every parameter has units and evidence class in the exported contract.

## Mechanism sources (qualitative only; no numerical rates derived)

- Li et al., Cell 2019, https://doi.org/10.1016/j.cell.2018.11.043 : dysfunctional
  melanoma CD8 can proliferate; this does not identify transition rates.
- Hui et al., Science 2017, https://pubmed.ncbi.nlm.nih.gov/28280247/ : PD-1/PD-L1
  inhibits costimulation; the normalized scalar model is our simplification.
- Pauken et al., Science 2016, https://doi.org/10.1126/science.aaf2807 : checkpoint
  blockade does not fully reset exhausted-cell epigenetic state (mouse context).
- Andreatta & Carmona, 2021, https://doi.org/10.1016/j.csbj.2021.06.043 : rank-based
  scoring motivation. The implemented capped-rank score is explicitly specified,
  not claimed to reproduce an unpinned external package version.

## Exact hazard factors and capacity semantics

For class activation/reinvigoration the rate multiplier is S; effector uses
0.1+0.9S; rest uses 1-S; dysfunction uses 0.15+S+local tumor fraction. Division
uses its state-specific rate times 0.2+0.8S; death uses its state-specific rate.
Influx is specified separately by state. These multipliers and numerical rates
are assumptions. Myeloid fraction is measured in context but has no fitted
suppressive coefficient in v0.4. Site/chemistry are metadata, not fitted effects.
M3 implements context-conditioned rates with selected base-rate estimation;
it is not the full prespecified hierarchical covariate/regimen model.

Profiles expose four separate axes. Cytotoxicity multiplies actual contact
damage; dysfunction reduces migration speed by a factor 1-0.4*dysfunction.
Proliferation/persistence are descriptive capacity axes; their kinetic
counterparts are the independently editable state-specific division/death
rates. They are not multiplied into these rates a second time. State changes
select new profiles **and** new division/death hazards. Persistence is not
measured survival, and state age is tracked but does not yet modify hazards.

Contact damage is attack_damage_per_hour * current cytotoxicity * recognition
* (1-checkpoint) * dt. NK uses its initial RNA-derived cytotoxic program and no
state-specific PD1 activity; this is a declared simplification. Antigen/TCR
reactivity is a seeded Bernoulli assignment at reactive_fraction, not inferred
from receptor RNA. Division inherits that assignment. No specificity assay is
implied. Exposure-zero retains the treatment window boundaries with zero
amplitude, preserving the paired discretization for non-grid schedules.
