# Model, assumptions and numerical behavior

The v0.2 software prototype is a **non-spatial continuous-time branching Markov
process**. Each live agent has a stable ID and one discrete state. Parameters
are read from the bundle's YAML or an explicit `--model` override. The default
model is illustrative and carries no treatment interpretation.

| Event | Hazard | Effect |
|---|---|---|
| State transition i → j | q_ij × number in i | One agent changes state; its ID is retained |
| Division in i | b_i × number in i | Net +1 birth; parent remains and new ID inherits state |
| Death in i | d_i × number in i | One live agent is removed |
| Influx into i | a_i | One new agent is added, without a measured parent |

Transition, division and death rates are **per cell per hour**; influx is
**cells per hour for the entire selected population**. Subsampling a population
does not automatically rescale influx. Choose that assumption explicitly.
Division has no explicit cell cycle, molecular inheritance or size constraint.

The engine uses the Gillespie direct method: draw the next waiting time from
an exponential distribution with the sum of all hazards; select an event
channel proportionally to its hazard; select a source agent uniformly within
the affected state. This is exact for the declared time-homogeneous linear
process, up to numerical precision. There is no tau-leaping or per-step
approximation. `--interval` changes observation snapshots, not event timing.

For the column-vector expectation m, the uncapped model obeys
`dm/dt = A m + a`. Off-diagonal `A[j,i] = q_ij`; diagonal
`A[i,i] = b_i - d_i - sum_j(q_ij)`. An augmented matrix exponential supplies
the exact expectation including constant influx. Dashed expectation curves are
neither confidence intervals nor observed data. They describe the uncapped
process; stopped runs are reported only through their observed interval.

## Editing assumptions

```bash
vtsc model-template my_model.yaml
vtsc check outputs/input.zip --model my_model.yaml --hours 24
vtsc run outputs/input.zip outputs/changed_assumptions --model my_model.yaml --hours 24
```

Start from [demo-model.yaml](../examples/demo-model.yaml). Every state is explicit;
an omitted state rate is zero. A transition is listed once, with different source
and target states. Negative/nonfinite rates, duplicate edges and unknown states
are rejected. `rate_provenance` must document the source or assumption.
The input labels must be compatible with the model vocabulary.

The `no-events` scenario sets every transition, division, death and influx rate
to zero. It is a software negative control, not an untreated biological control.
Baseline means the supplied model; it does not mean untreated patients.

## Seeds, caps and scope

The seeded NumPy generator determines event times and choices. Identical inputs,
configuration and seed in the tested numerical environment reproduce events.
Changing a model changes its event sequence even with the same seed; paired
seeds are not guaranteed variance reduction. Artifact ZIP timestamps and HDF5
container metadata need not be byte-identical across environments.

Maximum population/event limits stop a run before the next unrecorded event.
They are safeguards, not carrying capacities. Inspect `complete`, `stop_reason`,
`observed_hours` and the event journal. Zero hazards advance to the horizon with
unchanged state; extinction stays empty unless influx is present. There is no
checkpoint resume in this version; restart from the verified input and seed.

The simulation does not establish identifiable rates, biological lineage
relationships, clinical response, spatial structure or causal treatment effects.
The SPECIAL model ladder and held-out-patient gates are unchanged. Software
testing is a separately authorized purpose, recorded in A03.
