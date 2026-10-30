# Jev (TypeSafe) API: live test notes

Tested live on 2026-10-30 with `scripts/jev_smoke_test.py`. Raw parsed results
(no key, no headers, no cookies) are in
[`reports/jev-smoke/2026-10-30/smoke_results.json`](../reports/jev-smoke/2026-10-30/smoke_results.json).

## Endpoint and authentication

| Item | Value |
|---|---|
| Endpoint | `POST https://api.typesafe.ai/v1/systemone` |
| Auth | `Authorization: Bearer <key>` |
| Key source | environment variable `TYPESAFE_API_KEY` (or `JEV_API_KEY`); never in files, Git, logs or notebooks |
| Model | `jev-1.14.0` (pinned; an unknown version returns HTTP 400) |
| Content type | `application/json` in and out |
| Latency seen | 0.15 to 0.4 s per request with 1 to 2 questions |
| Network | the host `api.typesafe.ai` must be reachable (allow-listed in cloud environments) |

## Request format

```json
{
  "model": "jev-1.14.0",
  "state": {"setting": "...", "prompt_version": "...", "semantics": "..."},
  "questions": {
    "c1": {
      "type": "choice",
      "instructions": {"task": "Choose one permitted local behavior for this cell.",
                       "local_observation": {"cell_type": "basal keratinocyte", "exposure": {"IL17A": "high"}}},
      "criteria": {"proliferate": "description", "differentiate": "description", "rest": "description"}
    }
  }
}
```

- `state` is required (HTTP 422 when missing). Free-form keys inside it are accepted.
- `questions` is a map of question id to question; many cells can share one request
  (the repo batches up to 20 per request, see `src/vtsc/jev_payload.py`).
- `instructions` takes free-form JSON. `criteria` is the permitted action menu:
  option id to a plain-language description.

## Response format

```json
{
  "model": "jev-1.14.0",
  "answers": {
    "c1": {"type": "choice", "choice": "proliferate", "confidence": 0.8,
           "probabilities": {"proliferate": 0.85, "differentiate": 0.01,
                             "secrete_inflammatory": 0.14, "rest": 0.0}}
  },
  "usage": {"input_tokens": 859, "output_tokens": 117}
}
```

**Yes, Jev returns action probabilities.** Each answer has a probability for every
option in `criteria` (two-decimal values summing to 1), a `choice` equal to the
argmax, and a separate `confidence`. Usage is reported per request. This matches
what `validate_response` in `src/vtsc/jev_provider.py` already expects, so the
existing provider works unchanged against the live API.

Errors are JSON with a `detail` field:

| Case | HTTP | `detail` |
|---|---|---|
| Bad key | 401 | `error_type: authentication_error` |
| Unknown model | 400 | `error_type: api_usage_error`, `Unknown model: ...` |
| Missing field | 422 | list of validation errors; **echoes the request body**, so never log it raw |

## Behaviour observed in the psoriasis probes

Option menu: `proliferate`, `differentiate`, `secrete_inflammatory`, `rest`.

| Probe | Result |
|---|---|
| Basal keratinocyte with high IL-17A/IL-22, KRT16/KRT6A/MKI67 up | proliferate 0.85, secrete_inflammatory 0.14 |
| Spinous keratinocyte, no IL-17A, KRT10/KRT1 up | differentiate 0.90 |
| Same two cells, `state.setting` lesional vs non-lesional | almost identical (within 0.03) |
| Quiet basal keratinocyte, same request sent twice | 0.64 / 0.65 proliferate: not fully deterministic, small jitter |
| Quiet basal keratinocyte plus a `guidebook_prior` (secrete_inflammatory 0.55) in the observation | secrete_inflammatory 0.76 (0.91 in an earlier probe without `semantics`) |

What this means for the design:

1. **Context must be in the local observation, not only in `state.setting`.** The
   lesional/non-lesional label alone barely moved the answers; the cytokine
   exposure and marker genes drove them. The perturbation "prompt" therefore has
   to change what each cell observes (cytokine levels, neighbours, guidebook
   entry), which is how the plan's context block and guidebook are designed.
2. **Do not put the guidebook prior inside the prompt and also mix it in the
   kernel.** Jev treats a prior in the observation as strong evidence and sharpens
   it (0.55 became 0.76 to 0.91). Pick one route. The recommended route is the
   plan's kernel-side mixture `(1-w)·p_jev + w·p_prior` with the prior kept out of
   the prompt, so `w` keeps a clear meaning.
3. **Answers vary slightly between identical calls and with wording** (`semantics`
   text changed the quiet-basal answer from 0.69 to 0.64). Keep the response
   journal (`jev_provider.py`) for exact replay, pin `prompt_version`, and treat
   differences under about 0.05 as noise.
4. Probabilities are rounded to 0.01, so a rare action with a true probability
   below 0.005 is reported as 0. Use a small floor if a sampled action must stay possible.

## Running the smoke test

```bash
export TYPESAFE_API_KEY=...        # do not commit; do not put in notebooks
python scripts/jev_smoke_test.py reports/jev-smoke/$(date +%F)
```

Eight requests, five billed (about 3.5k input and 0.4k output tokens). The script
refuses to write its output if the key string appears in it.
