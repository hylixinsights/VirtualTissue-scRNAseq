# Local v0.4 handoff

The package is prepared locally. No push, remote release, PR or external message
was performed. The intended owner-controlled destination is
https://github.com/hylixinsights/VirtualTissue-scRNAseq.git .

## Install from source or wheel

Use Python 3.12/3.13. Source: `python -m pip install -e ".[test,notebooks]"`.
Wheel: `python -m pip install dist/virtualtissue_sc-0.4.0-py3-none-any.whl`.
`requirements-product-v04.lock.txt` records the tested macOS environment;
pyproject.toml pins the portable direct dependencies. R exports use the existing
renv environment. No SPECIAL/Lin data are included in the wheel or source archive.

Run `vtsc demo outputs/group-test --cells 120 --hours 24`, then
`vtsc explore outputs/group-test/comparison --open`. Own-data instructions are
in docs/02_your_own_data.md. A local HTTP server is an alternative if your viewer
host prohibits file URLs; serve only the generated output directory.

## Owner publication steps (not executed by the agent)

Review `git status`, `git log` and DATA_LICENSE.md. Confirm the destination
repository exists and the remote is correct before using these commands:

```bash
git remote add origin https://github.com/hylixinsights/VirtualTissue-scRNAseq.git
git push -u origin HEAD
```

If origin already exists, inspect it rather than adding/replacing it blindly.
Do not force-push. CI is configured for Linux/Windows Python 3.12/3.13, plus
notebook and clean-wheel smoke tests. Remote CI has not run in this local task.
After publication, users can clone that URL and follow README installation.
The local source archive and wheel can also be shared for group software tests,
retaining the code licence and source attributions.
