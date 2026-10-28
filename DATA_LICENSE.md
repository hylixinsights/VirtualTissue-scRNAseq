# Data licence and permitted project use

## Software-testing examples (v0.2)

Original simulator code and programmatically generated synthetic fixtures are
MIT-licensed. Synthetic fixtures contain no patient measurements. The optional
10x PBMC3k example is CC BY 4.0, explicitly declared on its original dataset
page; the source and attribution are recorded in `examples/datasets.json` and
the locally downloaded `SOURCE.json`. PBMC-derived bundles and results retain
the source terms and should identify the subset/annotation/simulation changes.
PBMC1k v3 is also explicitly CC BY 4.0 but has not yet been tested here.
See `docs/05_example_datasets.md` for URLs and attribution text.

Own-data inputs retain the user's applicable permissions/terms. Neither a
public software repository nor a `--source-license` field grants rights to
redistribute somebody else's data. The example pipeline does not read SPECIAL
or the Lin atlas, and does not depend on their licences or clinical eligibility.

## Historical SPECIAL research assets

The SPECIAL processed-data deposit (Joanna Poźniak, KU Leuven RDR,
[doi:10.48804/GSAXBN](https://doi.org/10.48804/GSAXBN), version 8) declares
[CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/).
The deposited metadata snapshot is `reports/evidence/rdr-v8-metadata.json`.

For SPECIAL-derived matrices, metadata exports, avatars and simulation-input
bundles, this project retains attribution, the licence link, and a description
of modifications; distributes adaptations under CC BY-NC-SA 4.0; and prohibits
commercial use. Keep these assets separate from Hylix commercial services.
Do not label them MIT. The licence does not require new, independent analysis
code to use the data licence or change the existing VirtualTissue code licence.

Suggested attribution: “Processed single-cell data from the SPECIAL cohort
(Poźniak et al., Cell 2024; Poźniak, Roda et al., Cancer Discovery 2025),
Joanna Poźniak, KU Leuven RDR, doi:10.48804/GSAXBN, version 8,
CC BY-NC-SA 4.0. Changes: [describe the export, filtering or transformation].”

Clinical article supplements retain their own source terms, recorded separately
in the manifest; the dataset licence is not automatically assigned to them.
Only public processed data are acquired here. No EGA controlled-access raw
sequencing files or direct patient identifiers are requested. Study pseudonyms
are retained. Original downloads and derived data live outside the Git repository.

This is the project data-use policy. Any release must include the applicable
source notices and transformation provenance for each asset.

## Cohort-expansion source screening

Lin2026 Zenodo v2 (`doi:10.5281/zenodo.19586163`) declares CC BY 4.0 in its
downloaded metadata. Its integrated object also incorporates upstream studies,
including SPECIAL according to Table S2A. Preserve source-specific terms;
do not treat the deposit-level licence as proof that SPECIAL-derived content
has lost its CC BY-NC-SA restrictions. Article supplements, GEO files and
preprints keep their own source terms. Unverified terms remain unresolved in
`data/COHORT_EXPANSION_MANIFEST.json`; this screening is not a redistribution
licence audit. No candidate data have been released or assigned the code's MIT
licence.
