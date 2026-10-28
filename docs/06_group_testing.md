# Group testing and bug reports

Start with `vtsc demo outputs/test_01 --cells 300 --hours 48 --seed 21`.
Open the baseline explorer and the no-events control. Check that the control
has no events and unchanged states/population. Run `python -m pytest -q`.
Then try PBMC3k, followed by one small, reviewed sample from your own data.

Useful tests include reopening the results ZIP on another computer, changing
the random seed, using a custom state vocabulary/rate file, deliberately reaching
an event cap, and correcting a deliberately mismatched metadata barcode.
The last case must fail visibly instead of silently attaching the wrong labels.

For an actionable bug report include:

- Commit/version, Python version and operating system.
- Exact command, model YAML, seed and expected behavior.
- Complete error text, or `summary.json` plus the observed discrepancy.
- The smallest synthetic/permitted input that reproduces it.
- Whether the same problem appears in the default synthetic demo.

Do not attach patient expression or identifying metadata to a public issue.
The generated synthetic example is the preferred reproducer. A local sharing
intent does not make a public GitHub repository private; do not commit the
group's private inputs or output ZIPs.

CI runs core tests, the freeze check and a short offline demo. It does not fetch
external biological data or run the clinical research study. A green CI run
does not establish biological validation. Platform results are only confirmed
after the corresponding jobs actually complete.
