# Notebook workflow

Install `pip install -e ".[notebooks]"` with Python 3.12/3.13.

1. Part 1 prepares the **public dissociated melanoma** (default), an offline
   synthetic melanoma fixture, or your own annotated biopsy.
2. Parts 2–3 run spatial agents with dynamic CD8 states, anti-PD-1, exposure-zero,
   recognition-off and transitions-off controls. Behavior choices default to mocked Jev responses. Live Jev requires an environment key, explicit paid opt-in and a request ceiling. State hazards remain kernel-owned.
3. Part 4 verifies saved results and exports the synchronized offline viewer, inspects actual decisions and tests replay recovery.

Use the same `WORK` directory across parts. Outputs are immutable by name: choose
a new directory for a new experiment. Real input supports Cell Ranger H5,
H5AD/count layers and 10x Matrix Market with barcode-aligned metadata.

The notebook defaults to the real melanoma download. The local validation script
executes all three notebooks with the synthetic melanoma path and installation
disabled, in fresh temporary workspaces and real kernels. Hosted Colab is not
yet verified, and GitHub clone instructions need this release to be uploaded.
