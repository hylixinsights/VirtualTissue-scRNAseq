"""Compare public EGA sample metadata with clinical S1, without expression.

Ordinal patient/day matches are provisional: they do not establish the E193xx
library-to-sample link or prove that a record is scRNA rather than Visium.
"""
import csv
import hashlib
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / 'reports/evidence/cohort-expansion'
records = []
for suffix, accession in [('665', 'EGAD50000000665'), ('2902', 'EGAD50000002902')]:
    path = EVIDENCE / f'lin-ega-samples-{suffix}.json'
    source = json.loads(path.read_text())
    assert isinstance(source, list)
    for sample in source:
        title = sample['title']
        match = re.fullmatch(r'(?:MM_pat|Patient\s*)(\d+)_(mono|combi)-day(\d+)', title)
        assert match, title
        patient, treatment, day = match.groups()
        records.append({
            'dataset': accession, 'ega_sample_id': sample['accession_id'],
            'title': title, 'source_subject_id': sample['subject_id'],
            'description': sample['description'], 'phenotype': sample['phenotype'],
            'provisional_s1_patient': f'P{int(patient)}', 'day': int(day),
            'arm_as_reported': treatment, 'specific_drug_reported': False,
            'sample_id_mapping_verified': False,
            'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        })
with (ROOT/'tables/cohort_expansion_lin_ega_samples.tsv').open('w', newline='') as h:
    w = csv.DictWriter(h, fieldnames=list(records[0]), delimiter='\t', lineterminator='\n')
    w.writeheader(); w.writerows(records)
with (ROOT/'tables/cohort_expansion_lin_samples.tsv').open() as h:
    clinical = list(csv.DictReader(h, delimiter='\t'))
ega_keys = {(x['provisional_s1_patient'], x['day']) for x in records}
missing = [x['sample_id'] for x in clinical if (x['patient_id'], int(x['day_as_reported'])) not in ega_keys]
summary = {
    'ega_records': len(records), 'unique_accessions': len({x['ega_sample_id'] for x in records}),
    'provisional_unique_patient_day_keys': len(ega_keys),
    'clinical_sample_ids_without_ordinal_day_match': missing,
    'specific_drug_names_available': False,
    'histologic_subtype_available': False,
    'clinical_pair_cd8_eligibility_inferred': False,
    'caution': 'Provisional patient/day matching only. Neither absence nor presence in these public records establishes paired scRNA availability. The mixed scRNA/Visium deposit must be reconciled at library level.',
    'chemistry_source_discrepancy': 'EGA dataset descriptions name 5-prime v2; the publisher key-resource table names 5-prime v1.1. Preserve both pending sample-level reconciliation.',
}
(EVIDENCE/'lin-ega-clinical-audit.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps(summary, indent=2))
