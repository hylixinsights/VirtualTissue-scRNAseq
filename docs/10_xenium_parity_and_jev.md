# Comparison with the latest published Xenium release

Reference: [VirtualTissue-Xenium v0.1.3](https://github.com/hylixinsights/VirtualTissue-Xenium/releases/tag/v0.1.3),
published 2026-10-28, commit `16d0a307d118a42143a7034c81b5dc366ca2d62c`.
Version 0.1.3 was published during this task; the final comparison was refreshed
against its tag. All bundled v2 runtime files are unchanged from v0.1.2. The
provider was originally adapted from v0.1.2; attribution is retained in NOTICE.md.
No Xenium runtime files were edited. Metadata/hashes are in
`reports/xenium-parity/reference.json`.

## Architectural comparison

| Capability | Xenium v0.1.3 | scRNAseq v0.5.0 |
|---|---|---|
| Starting geometry | Measured Xenium coordinates/polygons with fitted cells | Assumed positions and non-overlapping discs from one dissociated biopsy |
| Movement | Mechanical displacement, bounded programs, target approach and field gradients | Bounded cue/exploratory migration; target approach in provider modes; exclusion and blocked outcomes |
| Behavior decisions | Jev Choice ranks kernel-permitted menus | Optional real Jev Choice, explicit mock transport, or legacy local rules |
| Physical fate | Kernel controls damage, death, cycle, division and clearance | Kernel controls the same classes of outcomes plus stochastic four-state CD8 dynamics |
| CD8 state dynamics | Xenium-specific state/competence assumptions | Memory/Activated/Cytotoxic/Dysfunctional hazards, age, Unknown and independent capacity axes |
| Drug mechanism | No named-drug experiment in published tutorial | Explicit normalized anti-PD-1 schedule/checkpoint model; assumed, not clinical PK |
| Molecular evidence | 20 GO:BP AUCell scores conditioned on type/depth | Provisional rank-state axes and log-count program summaries; **not the same AUCell representation** |
| Fields | Ten modeled fields, including oxygen, metabolites and immune signals | One generic cue plus uniform exposure; no claim of full Xenium field parity |
| Decision scheduling | Type/event-dependent wake rules and persistent programs | First eligible checkpoint, then 1 h (0.5 h for contacted effectors), at most 24 requests/cell/day |
| Prompt context | User choice: none/custom/bundle, vtx-choice-3 | Fixed, explicit exploratory melanoma setting, vtsc-behavior-1; no organ inference from filenames |
| Rounding | Preserves reported values; compatibility copy ties a chosen option within 0.01 of the maximum | Same reported tolerance; direct validated choice needs no tie adjustment in the kernel |
| Diagnostics | AUCell coverage, wake eligibility and aggregate cycle-probability summaries | Exact first-round payload/count preview and individual behavior/fate journals; no aggregate AUCell/cycle report |
| Request controls | ≤20 questions / 40 KB, grouped effectors, fixed model, request ceiling | Same bounds and isolated CD8/NK groups; explicit paid opt-in; no retry or redirect |
| Recovery | Configuration-bound durable request journal, replay from t0 | Same architecture; cached answers replay without HTTP; explicit new run/ceiling required for further paid calls |
| Inspection | Recorded action probabilities, events, fields, contacts | Behavior probabilities and separate CD8 fate hazards/draws, actual movement/contact and synchronized arms |

The projects now share the decision/kernel separation and major interaction
classes. They are **not identical simulators**: dissociated RNA cannot recover
measured histology, and the current scRNAseq environment is simpler. Neither
software verification nor this comparison establishes clinical validity.

## Select the decision mode

`local` preserves the v0.4 programmed behavior trajectory. CD8 fates still use
the common stochastic hazard kernel. `mock-jev` exercises the actual request,
validation, batching, journaling and behavior execution code with deterministic
programmed answers. `jev` calls the pinned TypeSafe model. Mock success is not
live Jev or biological validation.

```bash
# Complete free run; no key or network required
vtsc demo outputs/jev-contract-test --cells 120 --hours 12 --decision-mode mock-jev
vtsc explore outputs/jev-contract-test/comparison --open

# Inspect exactly what the first decision round would send (no paid calls)
vtsc check-melanoma outputs/my-input.zip outputs/preflight

# Only the operator explicitly authorizes paid calls; key comes from the environment
# Set TYPESAFE_API_KEY securely in your shell, never in a notebook/config or Git.
vtsc run-melanoma outputs/my-input.zip outputs/jev-run --decision-mode jev \
  --allow-paid --max-http-attempts 100 --jev-model jev-1.14.0 --hours 12
```

The request limit is cumulative HTTP attempts including failures, not money or
a guarantee of finishing the requested horizon. Model pricing/account limits
are outside this software. Live mode requires a key plus **both** explicit flags;
preflight and result exploration do not read a key. No live Jev calls were made
during this implementation/publication task.

The three notebooks expose DECISION_MODE, ALLOW_PAID and MAX_HTTP_ATTEMPTS.
Their default is the free mock contract test. Keep the same mode/model/seed/base
configuration across arms; comparisons reject mixed behavior policies.

## What Jev receives and owns

The endpoint is `POST https://api.typesafe.ai/v1/systemone`. The default pinned
model is `jev-1.14.0`, prompt `vtsc-behavior-1`. See TypeSafe's
[Choice schema](https://docs.typesafe.ai/primitives/choice) and
[System One API](https://docs.typesafe.ai/api/system-one).
Requests include internal integer cell indices, lineage, limited initial score
summaries, simulated state/capacities, contacts, and permitted behavior descriptions.
Only CD8 observations include the simulated drug/checkpoint context. Effectors
are isolated from other lineages in batches. No original barcodes, patient/sample
IDs, raw counts, gene lists, coordinates, full contract, desired response or
kinetic rate is sent. Review the saved preflight before using private data.

The API must return the pinned model, exactly the offered options and finite
probabilities/confidence. Reported probabilities are preserved; normalization
within rounding tolerance is explicit. A chosen option may trail the rounded
maximum by at most 0.01. Choices are used directly, not sampled as kinetic rates.
The state-transition/division/death hazards remain separate and kernel-owned.

Behavior rounds read a pre-step snapshot and validate all answers before any
step is committed. Cells execute the selected program until another checkpoint;
current menus are checked again after state changes and earlier cell actions.
An unavailable choice becomes logged rest, never an invented target or an illegal
attack. Melanoma growth is a choice; completed-cycle division is physical and
space-checked. A provider cannot assign identity, create a death event or set a
rate. Fate changes can alter the capacities used by the chosen behavior.

## Journal, limits and resumption

A durable sibling `RUN.provider.jsonl` is written **before** each HTTP send and
fsynced. Valid responses are reused by a request hash bound to input, runtime,
model, prompt, configuration, seed and arm. Budget checks precede a whole round;
partial/failed rounds do not mutate tissue. The saved run has an explicit stop
reason and its complete prefix remains inspectable. The final ZIP also contains
the provider journal, including model/usage and a hash chain. Credentials are
never persisted. Unresolved or failed attempts are never automatically resent.

```bash
# Recreate the recorded prefix without sending new requests, even for a live run
vtsc recover-melanoma outputs/jev-run.zip outputs/replayed
# Continue a budget-limited run only through explicit paid authorization and a new output
vtsc run-melanoma outputs/my-input.zip outputs/jev-continued --decision-mode jev \
  --allow-paid --max-http-attempts 200 --replay-journal outputs/jev-run.provider.jsonl \
  --jev-model jev-1.14.0 --hours 12
```

Reuse exactly the source configuration/seed/arm/layout/population cap. The new
ceiling is lifetime attempts represented by the copied journal, not 200 extra
calls. A failed/unresolved journal stops rather than guessing whether a paid
request was billed. Recovery with a longer horizon stops when no recorded answer
exists; it never silently switches to live access. Earlier v0.4 recordings need
their original version for deterministic recovery and remain readable for replay.
