# v0.5 verification and Xenium comparison

Scope: exploratory software implementation, not scientific validation. The
owner requested comparison with the latest published Xenium release and then
publication to `hylixinsights/VirtualTissue-scRNAseq`. Paid Jev calls remain
unauthorized and **zero actual Jev HTTP requests were made**.

## Reference and resulting behavior

GitHub identified [Xenium v0.1.3](https://github.com/hylixinsights/VirtualTissue-Xenium/releases/tag/v0.1.3)
as the latest published release: commit
`16d0a307d118a42143a7034c81b5dc366ca2d62c`, published 2026-10-28.
Version 0.1.3 appeared during this task; the comparison was refreshed against
its tagged source. The bundled runtime is byte-identical to v0.1.2. Metadata/source hashes are in [reference.json](reference.json).
No Xenium files were edited in this task.

The main gap was the absence of a live behavior-provider adapter. Version 0.5
adds optional Jev Choice decisions through the published API contract, isolated
effector batches, bounded target approach, action revalidation, a pinned model,
request ceilings and a durable configuration/model-bound journal. Kernel-owned
state transitions, division, death and physical constraints remain separate.
The explorer shows reported behavior probabilities alongside the separate
stochastic fate hazards/draws. The default local policy preserves v0.4 motion.

[The comparison table and operating guide](../../docs/10_xenium_parity_and_jev.md)
state the important differences: assumed discs versus measured geometry,
one cue versus ten fields, rank-based RNA axes versus GO:BP AUCell, and periodic
versus type/event-dependent scheduling. No numerical or biological equivalence
with Xenium is claimed. Mock outputs are programmed contract tests, not Jev predictions.

## Evidence

| Check | Result | Evidence |
|---|---|---|
| Regression and new provider/mechanics tests | 61 passed | [tests.txt](tests.txt), [pytest.xml](pytest.xml) |
| Three notebooks through real local kernels | Passed; synthetic input and mock transport | [notebook-validation.json](notebook-validation.json) |
| Wheel installed into a clean venv, isolated imports | Demo, ZIP reopening, preflight, cached recovery and dependency check passed | [wheel-validation.json](wheel-validation.json) |
| Synthetic paired 120-cell, 12 h run | Both arms completed; actual state transition and contact damage recorded | [integration.json](integration.json) |
| Public 10x, proportional 400-cell, 12 h runs | Anti-PD-1, zero exposure and recognition-off completed; 75/75/0 attacks | [integration.json](integration.json) |
| Cached recovery | Reproduced frames and decision prefix without transport calls | Tests and integration report |
| Previous local policy | All 400-cell 24 h recorded frames equal the accepted v0.4 run exactly | Integration report |
| Browser | Behavior/fate distinction, journal expansion, synchronized final counts and synthetic transition | [browser-validation.json](browser-validation.json) |
| Frozen study and historical evidence | Hashes verified unchanged | [publication-audit.json](publication-audit.json) |

The real example has 27 provisionally annotated CD8 and eight melanoma-like
cells. Its 12 h mock runs had no CD8 transitions, tumor deaths or tumor response;
that is the recorded outcome, not a reason to enrich cells or change seeds.
The synthetic run demonstrates dynamic states. The real annotation audit and
its limitations remain in [v0.4 verification](../product-v04/VALIDATION.md).
`http_attempts` in mock usage records counts injected transport batches, **not
external HTTP requests or charges**. The 400-cell run used 286 such batches;
the first-round preflight has 363 questions in 20 batches.

Reproduction: `python -m pytest -q`, `python scripts/validate_notebooks.py`,
and `python -I scripts/validate_wheel.py` from the clean wheel environment.
`scripts/validate_jev_integration.py` uses the already verified local v0.4
public-input artifact; a fresh clone can regenerate that input using the README.
Generated run ZIPs and raw matrices are intentionally outside Git.

## Limits retained

The first remote run passed Linux checks but found a Windows checkout converting
frozen JSON resources to CRLF. Git attributes now preserve the required LF bytes
without changing the model or its expected hashes; archived reports/tables keep
their original bytes. See [the diagnosis](windows-checkout-fix.json) and the
subsequent Actions run for remote confirmation.

Live Jev behavior/latency/cost has not been tested. Hosted Colab and direct
file-protocol browser use have not been verified here. Remote OS checks are
reported by GitHub Actions, not inferred from macOS results. No clinical model
was fitted, no paired SPECIAL state changes were inspected, and no scientific
gate was advanced. The original research board, freeze, claim ledger and
earlier product evidence remain separate and unchanged.
