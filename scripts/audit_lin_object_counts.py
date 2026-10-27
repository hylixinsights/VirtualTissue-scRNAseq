"""Count lineage-labelled cells only; never read response or state proportions.

Clinical linkage by PatientN -> PN is provisional: deposited metadata omit
the E193xx library identifiers and P2 has two clinical on-treatment biopsies.
This is a necessary feasibility screen, not a final eligibility manifest.
"""
import csv
import hashlib
import json
import re
from pathlib import Path

import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
EVIDENCE = ROOT / 'reports/evidence/cohort-expansion'
META = ROOT.parent / 'data/Lin2026/derived/audit-quarantine/source_metadata.parquet'
THRESHOLDS = (15, 30, 50)
ARMS = {'aPD1': 'anti_pd1_monotherapy', 'aPD1+CTLA4': 'anti_pd1_plus_anti_ctla4'}


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def save_table(name, rows):
    assert rows
    with (ROOT / 'tables' / name).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


config_path = ROOT / 'config/cohort-expansion-audit.yaml'
config = yaml.safe_load(config_path.read_text())
assert config['lin_lineage_count_rule']['primary'] == 'celltype_r2 starts with CD8_'
columns = ['cohort', 'patient', 'sample', 'time_point', 'treatment',
           'celltype_main', 'celltype_adj', 'celltype_r2']
frame = pd.read_parquet(META, columns=columns)
assert frame.index.is_unique and not frame[columns].isna().any().any()
lin = frame.loc[frame.cohort.eq('SKCM_this study')].copy()
assert not lin.empty
assert set(lin.treatment.unique()) <= set(ARMS)
assert set(lin.time_point.unique()) == {'Pre', 'On'}
# The vocabulary, without per-state counts, explains the lineage rule.
vocabulary = lin[['celltype_main', 'celltype_adj', 'celltype_r2']].drop_duplicates()
vocabulary = vocabulary.sort_values(list(vocabulary.columns))
vocabulary['included_by_cd8_prefix'] = vocabulary.celltype_r2.str.startswith('CD8_')
save_table('cohort_expansion_lin_lineage_vocabulary.tsv', vocabulary.to_dict('records'))
lin['cd8_author_lineage'] = lin.celltype_r2.str.startswith('CD8_')
lin['cd8_coarse_label'] = lin.celltype_main.eq('CD8+T')
lin['cd8_union_diagnostic'] = lin.cd8_author_lineage | lin.cd8_coarse_label

clinical = pd.read_csv(ROOT / 'tables/cohort_expansion_lin_samples.tsv', sep='\t')
samples, crosswalk = [], []
for sample, group in lin.groupby('sample', observed=True):
    for column in ['patient', 'time_point', 'treatment']:
        assert group[column].nunique() == 1, (sample, column)
    patient = str(group.patient.iloc[0])
    match = re.fullmatch(r'SKCM_this study_Patient(\d+)', patient)
    assert match, patient
    pid = 'P' + match[1]
    time = str(group.time_point.iloc[0]).lower()
    arm = ARMS[str(group.treatment.iloc[0])]
    candidates = clinical.loc[clinical.patient_id.eq(pid) & clinical.timepoint.eq(time)]
    assert len(candidates) >= 1, (pid, time)
    assert set(candidates.arm_as_reported) == ({'Mono'} if arm == ARMS['aPD1'] else {'Combo'})
    samples.append({
        'cohort': 'Lin2026_in_house', 'patient_id': pid,
        'qualified_patient_id': f'Lin2026_in_house::{pid}',
        'source_patient_id': patient, 'source_sample_id': str(sample),
        'timepoint': time, 'arm': arm, 'total_cells': len(group),
        'cd8_author_lineage': int(group.cd8_author_lineage.sum()),
        'cd8_coarse_label': int(group.cd8_coarse_label.sum()),
        'cd8_union_diagnostic': int(group.cd8_union_diagnostic.sum()),
        'clinical_sample_candidates': ';'.join(candidates.sample_id),
        'clinical_day_candidates': ';'.join(candidates.day_as_reported.astype(str)),
        'clinical_link_status': 'ambiguous_multiple_on_biopsies' if len(candidates) > 1 else 'ordinal_timepoint_provisional',
        'exact_library_link_verified': False,
    })

for row in clinical.to_dict('records'):
    matches = [s for s in samples if s['patient_id'] == row['patient_id'] and s['timepoint'] == row['timepoint']]
    assert len(matches) <= 1
    crosswalk.append({
        'cohort': row['cohort'], 'patient_id': row['patient_id'],
        'clinical_sample_id': row['sample_id'], 'timepoint': row['timepoint'],
        'clinical_day': row['day_as_reported'], 'arm_as_reported': row['arm_as_reported'],
        'source_sample_candidate': matches[0]['source_sample_id'] if matches else '',
        'link_status': matches[0]['clinical_link_status'] if matches else 'patient_timepoint_absent_in_deposit',
        'exact_library_link_verified': False,
    })

patients = []
for pid in sorted({s['patient_id'] for s in samples}, key=lambda s: int(s[1:])):
    rows = [s for s in samples if s['patient_id'] == pid]
    assert len({s['arm'] for s in rows}) == 1
    pre, on = [s for s in rows if s['timepoint'] == 'pre'], [s for s in rows if s['timepoint'] == 'on']
    assert len(pre) <= 1 and len(on) <= 1
    patient = {'cohort': 'Lin2026_in_house', 'patient_id': pid,
               'qualified_patient_id': f'Lin2026_in_house::{pid}', 'arm': rows[0]['arm'],
               'paired_in_deposit': bool(pre and on),
               'missing_timepoint': 'pre' if not pre else ('on' if not on else ''),
               'exact_library_link_verified': False, 'final_eligibility': 'unresolved'}
    for rule in ['author_lineage', 'coarse_label', 'union_diagnostic']:
        for time, values in [('pre', pre), ('on', on)]:
            patient[f'{time}_cd8_{rule}'] = values[0][f'cd8_{rule}'] if values else None
        for cutoff in THRESHOLDS:
            patient[f'count_candidate_{rule}_{cutoff}'] = bool(
                pre and on and min(pre[0][f'cd8_{rule}'], on[0][f'cd8_{rule}']) >= cutoff)
    patients.append(patient)

# Reused SPECIAL is explicitly excluded from new-patient totals.
reused = frame.loc[frame.cohort.eq('SKCM_Plozniak')].copy()
reused['cd8_coarse_label'] = reused.celltype_main.eq('CD8+T')
reused_rows = []
for (patient, sample, time), group in reused.groupby(['patient', 'sample', 'time_point'], observed=True):
    reused_rows.append({'atlas_cohort': 'SKCM_Plozniak', 'source_patient': str(patient),
                        'source_sample': str(sample), 'timepoint': str(time),
                        'metadata_cells': len(group), 'coarse_cd8': int(group.cd8_coarse_label.sum()),
                        'count_as_new_patient': False, 'reason': 'SPECIAL_reused_in_atlas'})

special = pd.read_csv(ROOT / 'tables/cohort_expansion_special_by_arm.tsv', sep='\t')
assert all(special[f'count_only_candidate_{t}'].dtype == bool for t in THRESHOLDS)
summaries = {}
for cutoff in THRESHOLDS:
    per_rule = {}
    for rule in ['author_lineage', 'coarse_label', 'union_diagnostic']:
        values = {}
        for source_arm, arm in [('anti_pd1_monotherapy', ARMS['aPD1']), ('combination', ARMS['aPD1+CTLA4'])]:
            special_ids = sorted(special.loc[special.arm.eq(source_arm) & special[f'count_only_candidate_{cutoff}'], 'patient_id'])
            lin_ids = [p['patient_id'] for p in patients if p['arm'] == arm and p[f'count_candidate_{rule}_{cutoff}']]
            values[arm] = {'SPECIAL_ids': special_ids, 'Lin2026_ids': lin_ids,
                           'SPECIAL_n': len(special_ids), 'Lin2026_n': len(lin_ids),
                           'combined_count_candidates': len(special_ids) + len(lin_ids)}
        values['combined_total'] = sum(v['combined_count_candidates'] for v in values.values())
        per_rule[rule] = values
    summaries[str(cutoff)] = per_rule

outputs = {
    'cohort_expansion_lin_object_samples.tsv': samples,
    'cohort_expansion_lin_object_patients.tsv': patients,
    'cohort_expansion_lin_library_crosswalk.tsv': crosswalk,
    'cohort_expansion_reused_special.tsv': reused_rows,
}
for name, rows in outputs.items():
    save_table(name, rows)
summary = {
    'status': 'count_only_audit_no_G1_pass', 'lin_primary_lineage_rule': config['lin_lineage_count_rule'],
    'metadata_rows_all_cohorts': len(frame), 'atlas_cohorts': sorted(frame.cohort.unique().tolist()),
    'lin_metadata_cells': len(lin), 'lin_patients_in_deposit': len(patients),
    'lin_samples_in_deposit': len(samples), 'lin_paired_patients_in_deposit': sum(p['paired_in_deposit'] for p in patients),
    'lin_clinical_rows': len(clinical),
    'clinical_rows_without_patient_timepoint_in_deposit': [r['clinical_sample_id'] for r in crosswalk if r['link_status'] == 'patient_timepoint_absent_in_deposit'],
    'ambiguous_clinical_rows': [r['clinical_sample_id'] for r in crosswalk if r['link_status'] == 'ambiguous_multiple_on_biopsies'],
    'reused_special_patients_excluded_as_new': int(reused.patient.nunique()),
    'reused_special_samples_excluded_as_new': len(reused_rows),
    'threshold_summaries': summaries, 'original_development_plus_lockbox_floor': 14,
    'primary_count_floor_met': summaries['30']['author_lineage']['combined_total'] >= 14,
    'input_hashes': {str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p): sha(p) for p in [
        META, config_path, ROOT / 'tables/cohort_expansion_lin_samples.tsv', ROOT / 'tables/cohort_expansion_special_by_arm.tsv']},
    'output_hashes': {name: sha(ROOT / 'tables' / name) for name in outputs},
    'read_columns': columns, 'state_proportions_or_changes_computed': False,
    'response_used': False, 'raw_matrix_alignment_verified': False,
    'exact_drugs_verified': False, 'exact_library_link_verified': False,
    'cross_study_patient_independence_fully_verified': False,
    'remaining_checks': ['exact_library_IDs', 'P2_on_day_9_vs_42', 'missing_libraries',
                         'exact_drugs', 'histology', 'lesion_matching', 'chemistry',
                         'SPECIAL_interbiopsy_intervals', 'raw_matrix_alignment', 'cross_study_patient_independence'],
    'lockbox_assigned': False, 'model_fitted': False,
}
(EVIDENCE / 'lin-object-count-audit.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
print(json.dumps({k: summary[k] for k in ['lin_patients_in_deposit', 'lin_samples_in_deposit',
                 'lin_paired_patients_in_deposit', 'primary_count_floor_met', 'threshold_summaries']}, indent=2))
