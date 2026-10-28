# Notebooks: local and Google Colab

Use these in order for the first test:

1. `Part1_Prepare_scRNAseq.ipynb`: default synthetic data → verified `input.zip`.
2. `Part2-3_Run_Simulation.ipynb`: inspect model → baseline + no-events control.
3. `Part4_Explore_Results.ipynb`: verify and reopen saved outputs.

All notebooks use the installed `vtsc` implementation; they do not carry a
second simulation engine. Core software works without notebooks or internet
after installation. Run notebooks in a Jupyter-compatible editor with the
environment from `pip install -e ".[notebooks]"` selected as the kernel.

## Local paths

Setup finds the repository above the current working directory. Set
`INSTALL_PACKAGE=False` if the correct environment is already installed.
`WORK` defaults to `outputs/notebook_demo` inside the repository. Use the same
absolute `WORK` or exact bundle/run paths across notebooks. Change the output
name for a new run. An existing bundle/run is deliberately not overwritten.

## Colab paths and persistence

After the repository is pushed, open each notebook through Colab's GitHub
picker and save a copy in your own Drive. Use a Python 3 CPU runtime; Python
3.12/3.13 are supported. Setup clones the public code if needed and installs it
in the current session. It does not mount Drive automatically.

The `/content` runtime is temporary and is not shared automatically between
notebooks. Either download the Part 1 bundle and upload it to Part 2, or
explicitly mount **your own** Drive in each notebook:

```python
from google.colab import drive
drive.mount('/content/drive')
WORK = Path('/content/drive/MyDrive/VirtualTissue_scRNAseq/my_test')
WORK.mkdir(parents=True, exist_ok=True)
```

Run this after setup and before data/run cells, or change the setup's `WORK`
line to this path. All three notebooks must use the same folder, or explicit
paths to the shared ZIPs. Drive mounting is optional; no project service receives
your inputs. Data and outputs remain subject to your sharing permissions.

Colab may request a kernel restart after dependency installation. Restart, set
`INSTALL_PACKAGE=False`, then rerun setup and subsequent cells. An inline HTML
display may restrict JavaScript depending on notebook trust/browser policy;
download `explorer.html` and open it directly for the fully offline viewer.

## Own data

In Part 1 choose `DATA_SOURCE='own'`, then set `OWN_INPUT`, `OWN_METADATA` if
needed, `COUNT_LAYER`, `STATE_COLUMN`, `SAMPLE_COLUMN`, `SAMPLE_ID`, and
`MODEL_PATH`. Read `docs/02_your_own_data.md` first. Defaults are a synthetic
software example, not automatic annotation for arbitrary human data.
