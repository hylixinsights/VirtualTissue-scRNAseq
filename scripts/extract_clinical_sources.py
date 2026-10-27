"""Read the unmodified clinical workbooks and retain cell/row provenance.

Use the bundled Python runtime for the initial workbook extraction. This script
does not modify workbooks or inspect single-cell state outcomes.
"""
import json
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT.parent / "data/SPECIAL/raw/clinical"
SOURCES = [
    ("2024", "Pozniak_Cell2024_TableS1.xlsx", "Sup_table_1"),
    ("2025", "cd-24-1208_supplemental_table_s1_suppst1.xlsx", "A scRNA-seq"),
]


def main():
    extracted = {}
    for year, filename, sheet_name in SOURCES:
        workbook = openpyxl.load_workbook(RAW / filename, read_only=True, data_only=True)
        sheet = workbook[sheet_name]
        rows = sheet.iter_rows(values_only=True)
        headers = list(next(rows))
        records = []
        for index, row in enumerate(rows, 2):
            if row[0] is None:
                continue
            records.append({"source_row": index, **dict(zip(headers, row))})
        extracted[year] = {"filename": filename, "sheet": sheet_name, "headers": headers, "records": records}
        workbook.close()
    common = set(extracted["2024"]["headers"]) & set(extracted["2025"]["headers"])
    common -= {None, "Abbreviations"}
    by_year = {year: {row["Sample ID"]: row for row in source["records"]} for year, source in extracted.items()}
    for year, source in extracted.items():
        assert len(by_year[year]) == len(source["records"]), "Duplicate sample IDs"
    conflicts = []
    for sample in sorted(by_year["2024"].keys() & by_year["2025"].keys()):
        for column in sorted(common):
            a, b = by_year["2024"][sample][column], by_year["2025"][sample][column]
            if a != b:
                conflicts.append({"sample_id": sample, "column": column, "value_2024": a, "value_2025": b})
    summary = {year: {"samples": len(source["records"]),
                      "patients": len({r["Patient ID"] for r in source["records"]}),
                      "treatments": sorted({r["Treatment"] for r in source["records"]}),
                      "date_or_interval_columns": [h for h in source["headers"] if h and any(s in h.lower() for s in ("date", "day", "interval"))]}
               for year, source in extracted.items()}
    output = {"sources": extracted, "summary": summary, "conflicts": conflicts,
              "samples_only_2024": sorted(by_year["2024"].keys() - by_year["2025"].keys()),
              "samples_only_2025": sorted(by_year["2025"].keys() - by_year["2024"].keys()),
              "response_policy": "Stratification only; never a predictor",
              "same_anatomical_site_does_not_verify_same_lesion": True}
    path = ROOT / "reports/evidence/clinical-source-reconciliation.json"
    path.write_text(json.dumps(output, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps({"summary": summary, "conflicts": conflicts}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
