"""T1.2/T1.3 preparation: biopsy totals and metadata, never state composition."""
import itertools
import json
from collections import Counter
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq

from vtsc.provenance import write_json

ROOT = Path(__file__).resolve().parents[1]
CD8_LABELS = {f"{state}_CD8_Tcells" for state in ("Memory", "Activated", "Cytotoxic", "Dysfunctional")}


def main():
    sources = json.loads((ROOT / "reports/evidence/clinical-source-reconciliation.json").read_text())
    assert not sources["conflicts"] and not sources["samples_only_2024"] and not sources["samples_only_2025"]
    clinical = {r["Sample ID"]: r for r in sources["sources"]["2024"]["records"]}
    counts = {}; vocabulary = set()
    for path in sorted((ROOT.parent / "data/SPECIAL/derived/immune").glob("*/obs.parquet")):
        cols = ["sample_ID", "patient_ID", "Timepoint", "Response", "Tissue", "high_resolution_clusters"]
        obs = pq.read_table(path, columns=cols).to_pandas()
        sample = path.parent.name
        for column in ["sample_ID", "patient_ID", "Timepoint", "Response", "Tissue"]:
            assert obs[column].nunique(dropna=False) == 1, f"Inconsistent {sample}: {column}"
        row = clinical[sample]
        assert obs["sample_ID"].iloc[0] == sample
        assert int(obs["patient_ID"].iloc[0]) == int(row["Patient ID"])
        assert obs["Timepoint"].iloc[0] == row["Timepoint"]
        assert obs["Response"].iloc[0] == ("R" if row["Responder?"] else "NR")
        vocabulary.update(obs["high_resolution_clusters"].unique())
        counts[sample] = {"cells_in_immune_object": len(obs),
                          "author_annotated_cd8_cells": int(obs["high_resolution_clusters"].isin(CD8_LABELS).sum()),
                          "object_tissue": str(obs["Tissue"].iloc[0])}
    assert CD8_LABELS <= vocabulary
    rows = []
    for sample, r in sorted(clinical.items()):
        observed = counts.get(sample, {})
        rows.append({"sample_id": sample, "patient_id": f"GC{int(r['Patient ID'])}",
                     "timepoint": {"BT": "pre", "OT": "on"}[r["Timepoint"]],
                     "treatment": r["Treatment"],
                     "arm": "combination" if "+" in r["Treatment"] else "anti_pd1_monotherapy",
                     "tissue": r["Tissue"], "anatomical_site": r["Biopsy site"],
                     "lesion_id": None, "same_lesion_verified": False,
                     "chemistry": r["Library type"], "preservation": r["Fresh/frozen sample"],
                     "response_for_stratification_only": r["Responder?"], "response_criterion": r["Criterion"],
                     "interval_days": None, "interval_lower_days": None, "interval_upper_days": None,
                     "total_tme_cells": None, "present_in_immune_object": sample in counts,
                     "cells_in_immune_object": observed.get("cells_in_immune_object", 0),
                     "author_annotated_cd8_cells": observed.get("author_annotated_cd8_cells", 0),
                     "object_tissue": observed.get("object_tissue"),
                     "clinical_source_row": r["source_row"], "clinical_source_sheet": "Sup_table_1"})
    samples = pd.DataFrame(rows)
    samples.to_parquet(ROOT / "tables/samples.parquet", index=False)
    samples.to_csv(ROOT / "tables/samples.tsv", sep="\t", index=False)
    pairs = []; patients = []
    for patient, group in samples.groupby("patient_id", sort=True):
        assert group["arm"].nunique() == 1 and group["treatment"].nunique() == 1
        assert group["response_for_stratification_only"].nunique() == 1
        pre = group[(group.timepoint == "pre") & group.present_in_immune_object]
        on = group[(group.timepoint == "on") & group.present_in_immune_object]
        ambiguous = len(pre) > 1 or len(on) > 1
        patient_pairs = []
        for a, b in itertools.product(pre.to_dict("records"), on.to_dict("records")):
            pair = {"patient_id": patient, "pre_sample": a["sample_id"], "on_sample": b["sample_id"],
                    "arm": a["arm"], "pre_cd8": a["author_annotated_cd8_cells"], "on_cd8": b["author_annotated_cd8_cells"],
                    "tissue_match": a["tissue"] == b["tissue"], "site_text_match": a["anatomical_site"] == b["anatomical_site"],
                    "chemistry_match": a["chemistry"] == b["chemistry"], "preservation_match": a["preservation"] == b["preservation"],
                    "same_lesion_verified": False, "multiple_possible_pairs": ambiguous, "interval_available": False}
            for threshold in [15, 30, 50]:
                pair[f"passes_cd8_{threshold}"] = min(pair["pre_cd8"], pair["on_cd8"]) >= threshold
            pairs.append(pair); patient_pairs.append(pair)
        reasons = ["interval_unavailable"]
        if group["arm"].iloc[0] != "anti_pd1_monotherapy": reasons.append("combination_therapy")
        if pre.empty or on.empty: reasons.append("missing_paired_immune_biopsies")
        if ambiguous: reasons.append("multiple_biopsies_need_prespecified_pair_rule")
        if patient_pairs and not any(p["passes_cd8_30"] for p in patient_pairs): reasons.append("cd8_below_30")
        record = {"patient_id": patient, "arm": group["arm"].iloc[0], "pre_biopsies_in_object": len(pre),
                  "on_biopsies_in_object": len(on), "status": "pending_interval" if reasons == ["interval_unavailable"] else "not_currently_eligible",
                  "eligible": False, "reasons": ";".join(reasons)}
        for threshold in [15, 30, 50]:
            record[f"potential_count_eligible_{threshold}"] = group["arm"].iloc[0] == "anti_pd1_monotherapy" and any(p[f"passes_cd8_{threshold}"] for p in patient_pairs)
        patients.append(record)
    eligibility = pd.DataFrame(patients)
    eligibility.to_parquet(ROOT / "tables/eligibility_pending.parquet", index=False)
    eligibility.to_csv(ROOT / "tables/eligibility_pending.tsv", sep="\t", index=False)
    pd.DataFrame(pairs).to_csv(ROOT / "tables/candidate_pairs.tsv", sep="\t", index=False)
    cross = samples.groupby(["chemistry", "timepoint", "arm", "tissue"], dropna=False).size().reset_index(name="biopsies")
    cross.to_csv(ROOT / "tables/chemistry_site_timepoint.tsv", sep="\t", index=False)
    summary = {"clinical_biopsies": len(samples), "clinical_patients": samples.patient_id.nunique(),
               "object_biopsies": len(counts), "object_patients": samples[samples.present_in_immune_object].patient_id.nunique(),
               "object_cells": int(samples.cells_in_immune_object.sum()),
               "missing_immune_samples": samples.loc[~samples.present_in_immune_object, ["sample_id", "patient_id", "timepoint"]].to_dict("records"),
               "patients_by_treatment": samples.drop_duplicates("patient_id").treatment.value_counts().to_dict(),
               "clinical_patients_with_pre_and_on": int((samples.groupby("patient_id").timepoint.nunique() == 2).sum()),
               "multiple_biopsy_patients": samples.groupby(["patient_id", "timepoint"]).size().loc[lambda s:s>1].reset_index(name="biopsies").to_dict("records"),
               "count_eligible_upper_bounds": {str(t):int(eligibility[f"potential_count_eligible_{t}"].sum()) for t in [15,30,50]},
               "development_upper_bounds_after_four_lockbox": {str(t):(int(eligibility[f"potential_count_eligible_{t}"].sum())-4 if eligibility[f"potential_count_eligible_{t}"].sum() >= 4 else None) for t in [15,30,50]},
               "lockbox_four_allocatable_by_counts_only": {str(t):bool(eligibility[f"potential_count_eligible_{t}"].sum() >= 4) for t in [15,30,50]},
               "intervals_available": 0, "formal_eligibility_finalized": False,
               "source_note": "Upper bounds ignore unavailable intervals and unresolved multiple-biopsy pairing; author CD8 labels await T2.1 audit.",
               "cd8_label_definition": sorted(CD8_LABELS), "global_annotation_vocabulary": sorted(vocabulary),
               "state_compositions_or_changes_computed": False}
    write_json(ROOT / "reports/evidence/cohort-audit.json", summary)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
