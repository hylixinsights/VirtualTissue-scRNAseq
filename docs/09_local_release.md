# Install and test v0.5

The destination is [VirtualTissue-scRNAseq](https://github.com/hylixinsights/VirtualTissue-scRNAseq).
The owner explicitly authorized pushing this implementation after comparison
with the latest published Xenium release. This supersedes the local-only
publication restriction in the historical handoff prompt; it does not authorize
paid Jev calls or change the scientific study gates.

Use Python 3.12/3.13. From the cloned repository:

```bash
python -m pip install -e ".[test,notebooks]"
vtsc demo outputs/group-test --cells 120 --hours 12 --decision-mode mock-jev
vtsc explore outputs/group-test/comparison --open
```

The demo uses synthetic cells and programmed mock answers. For real data, use
the pinned public melanoma downloader or [your own biopsy](02_your_own_data.md).
The [Jev guide](10_xenium_parity_and_jev.md) describes preflight, paid opt-in,
request limits, recovery and differences from Xenium v0.1.3.

To build and install a wheel:

```bash
python -m pip wheel --no-deps . --wheel-dir dist
python -m pip install dist/virtualtissue_sc-0.5.0-py3-none-any.whl
```

`requirements-product-v05.lock.txt` records the tested macOS dependency set;
portable direct dependencies are pinned in pyproject.toml. R exports keep the
existing renv environment. No SPECIAL/Lin expression matrices are in the wheel
or source archive. Preserve the code licence and source/data attributions.

Use a new output name for each execution. If your viewer host prohibits file
URLs, serve only the generated output directory on localhost. The browser reads
recorded events and does not call Jev. File-protocol browser execution and hosted
Colab remain separate from the tested localhost/local-kernel workflows.

CI covers Linux/Windows Python 3.12/3.13, notebooks and a wheel installed outside
the checkout. Inspect the [Actions runs](https://github.com/hylixinsights/VirtualTissue-scRNAseq/actions)
for remote results. Local evidence and its limits are in
[v0.5 verification](../reports/xenium-parity/VALIDATION.md); v0.4 evidence remains
unchanged in `reports/product-v04/`.
