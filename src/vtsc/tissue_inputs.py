"""Dissociated melanoma inputs. No spatial coordinates are inferred from RNA."""
from __future__ import annotations
import hashlib
import json
import subprocess
import zipfile
from collections import Counter
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .inputs import json_bytes, load_input
from .provenance import file_hashes

MELANOMA = {
    'id': '10x_melanoma_gemx_3p', 'title': "10k Human DTC Melanoma, Chromium GEM-X Single Cell 3′",
    'source_url': 'https://www.10xgenomics.com/datasets/10k-human-dtc-melanoma-GEM-X',
    'download_url': 'https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_Human_DTC_Melanoma_3p_gemx_10k_Human_DTC_Melanoma_3p_gemx/10k_Human_DTC_Melanoma_3p_gemx_10k_Human_DTC_Melanoma_3p_gemx_count_sample_filtered_feature_bc_matrix.h5',
    'sha256': 'a2472fd15fbc708c50b95aff7f762872929a0365b9d06454cc558d29d4400091',
    'bytes': 56774638, 'cells': 10645, 'features': 38606,
    'license': 'CC-BY-4.0', 'license_url': 'https://creativecommons.org/licenses/by/4.0/',
    'attribution': '10x Genomics (2024), 10k Human DTC Melanoma, Chromium GEM-X Single Cell 3′, Cell Ranger 8.0.0. Subsampling, marker annotation and simulated geometry by VirtualTissue-scRNAseq.',
    'verified_date': '2026-10-27',
    'source_discrepancy': 'The descriptive text says dissociated melanoma tissue; the structured Anatomical Entity field says Blood. This discrepancy is retained, not silently corrected.',
    'use': 'Single donor software example; no paired treatment outcomes or measured positions.'}

PANELS = {
    'Melanoma': ['MLANA', 'PMEL', 'TYRP1', 'DCT', 'SOX10', 'MITF'],
    'T_cell': ['CD3D', 'CD3E', 'TRAC', 'CD247'],
    'NK': ['NKG7', 'GNLY', 'KLRD1', 'NCR1'],
    'Myeloid': ['LST1', 'LYZ', 'FCER1G', 'CTSS', 'AIF1'],
    'Fibroblast': ['COL1A1', 'COL1A2', 'DCN', 'LUM'],
    'Endothelial': ['PECAM1', 'VWF', 'EMCN', 'KDR'],
    'B_cell': ['MS4A1', 'CD79A', 'CD79B', 'CD19'],
}
PROGRAMS = {'cytotoxicity': ['PRF1', 'GZMB', 'GZMH', 'GNLY'],
            'proliferation': ['MKI67', 'TOP2A', 'PCNA', 'MCM2'],
            'interferon': ['ISG15', 'IFIT1', 'IFIT3', 'STAT1']}
TYPES = ['Melanoma', 'CD8_T', 'T_cell', 'NK', 'Myeloid', 'Fibroblast', 'Endothelial', 'B_cell', 'Unknown']


def download_melanoma(destination):
    """Retrieve only a pinned public matrix; curl is available on macOS/Linux/Windows 10+."""
    dest = Path(destination); dest.mkdir(parents=True, exist_ok=True)
    target = dest / 'filtered_feature_bc_matrix.h5'
    if not target.exists():
        partial = target.with_suffix('.partial')
        subprocess.run(['curl', '--fail', '--location', '--retry', '2', '--connect-timeout', '30',
                        '--max-time', '600', '--output', str(partial), MELANOMA['download_url']], check=True)
        if partial.stat().st_size != MELANOMA['bytes'] or file_hashes(partial)['sha256'] != MELANOMA['sha256']:
            raise ValueError('Melanoma source checksum mismatch; the partial file is not accepted.')
        partial.replace(target)
    if target.stat().st_size != MELANOMA['bytes'] or file_hashes(target)['sha256'] != MELANOMA['sha256']:
        raise ValueError('Cached melanoma matrix differs from the pinned public source.')
    (dest / 'SOURCE.json').write_bytes(json_bytes(MELANOMA))
    return {'input': str(target), **MELANOMA}


def expression_features(data):
    names = data.var.get('gene_symbols', pd.Series(data.var_names, index=data.var_names)).astype(str).to_numpy()
    x = sparse.csr_matrix(data.X, dtype=float)
    x = sparse.diags(10000 / np.asarray(x.sum(axis=1)).ravel()) @ x
    x.data = np.log1p(x.data)
    scores, coverage, detected = {}, {}, {}
    for key, genes in {**PANELS, **PROGRAMS, 'CD8': ['CD8A', 'CD8B']}.items():
        mask = np.isin(names, genes)
        coverage[key] = [g for g in genes if g in names]
        scores[key] = np.asarray(x[:, mask].sum(axis=1)).ravel() / len(genes)
        detected[key] = np.asarray((x[:, mask] > 0).sum(axis=1)).ravel()
    return scores, detected, coverage


def annotate_demo(scores, detected):
    """Fixed marker heuristic, not malignancy inference. Ambiguous cells remain Unknown."""
    keys = list(PANELS)
    values = np.column_stack([scores[k] for k in keys])
    order = np.argsort(values, axis=1, kind='stable')
    best, second = order[:, -1], order[:, -2]
    rows = np.arange(len(best))
    margin = values[rows, best] - values[rows, second]
    labels = np.asarray(keys, dtype=object)[best]
    support = np.column_stack([detected[k] for k in keys])[rows, best]
    labels[(values[rows, best] < .5) | (support < 2) | (margin < .15)] = 'Unknown'
    labels[(labels == 'T_cell') & (detected['CD8'] > 0)] = 'CD8_T'
    return labels, margin


def write_bundle(payload, output):
    output = Path(output)
    if output.exists(): raise FileExistsError(output)
    content = json_bytes(payload)
    manifest = {'kind': 'vtsc_tissue_input', 'schema_version': 1,
                'files': {'input.json': hashlib.sha256(content).hexdigest()}}
    output.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        for name, value in [('input.json', content), ('manifest.json', json_bytes(manifest))]:
            entry = zipfile.ZipInfo(name, (2026, 1, 1, 0, 0, 0)); entry.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(entry, value)


def read_tissue_bundle(path):
    with zipfile.ZipFile(path) as z:
        infos = z.infolist()
        if len(infos) != 2 or {i.filename for i in infos} != {'input.json', 'manifest.json'} or sum(i.file_size for i in infos) > 50_000_000:
            raise ValueError('Invalid tissue input bundle members or size.')
        content = z.read('input.json'); manifest = json.loads(z.read('manifest.json'))
    if manifest != {'kind': 'vtsc_tissue_input', 'schema_version': 1, 'files': {'input.json': hashlib.sha256(content).hexdigest()}}:
        raise ValueError('Tissue input checksum or manifest mismatch.')
    payload = json.loads(content); cells = payload['cells']
    if not cells or len({c['cell_id'] for c in cells}) != len(cells): raise ValueError('Invalid cell identifiers.')
    for c in cells:
        if c['cell_type'] not in TYPES or not c['cell_id'].strip() or len(c['cell_id']) > 200:
            raise ValueError('Invalid cell type or identifier.')
        if set(c['programs']) != set(PROGRAMS) or any(not np.isfinite(v) or not 0 <= v <= 1 for v in c['programs'].values()):
            raise ValueError('Invalid expression program values.')
    return payload


def prepare_tissue(input_path, output, *, metadata=None, layer=None, cell_type_column='cell_type',
                   sample_id=None, recipe=None, max_cells=400, seed=21, source_license='user supplied'):
    if not 2 <= max_cells <= 2000 or seed < 0: raise ValueError('Use 2–2000 cells and a nonnegative seed.')
    data, sources = load_input(input_path, metadata, layer)
    input_n = data.n_obs
    if 'sample_id' in data.obs:
        values = data.obs.sample_id
        if values.isna().any() or values.astype(str).str.strip().eq('').any(): raise ValueError('Missing sample_id.')
        samples = values.astype(str).unique()
        if sample_id is None and len(samples) != 1: raise ValueError('Multiple samples: choose --sample-id; never pool biopsies.')
        sample_id = sample_id or samples[0]
        if sample_id not in samples: raise ValueError('Requested sample is absent.')
        data = data[values.astype(str).eq(sample_id)].copy()
    else: sample_id = sample_id or 'declared_single_biopsy'
    for col in ['patient_id', 'timepoint']:
        if col in data.obs and data.obs[col].nunique(dropna=False) > 1: raise ValueError(f'Multiple {col} values in selected biopsy.')
    scores, detected, coverage = expression_features(data)
    if any(not coverage[k] for k in PROGRAMS):
        raise ValueError('No genes found for an RNA program; supply human gene symbols in var_names or var.gene_symbols.')
    if recipe == 'melanoma-markers':
        if any(len(coverage[k]) < 2 for k in PANELS):
            raise ValueError('Marker recipe needs human gene symbols and at least two genes per lineage panel.')
        labels, margins = annotate_demo(scores, detected)
        annotation = 'provisional fixed marker recipe; Unknown retained; no CNV/malignancy validation'
    elif recipe is None:
        if cell_type_column not in data.obs: raise ValueError(f'Supply reviewed {cell_type_column!r} or explicitly use --recipe melanoma-markers.')
        labels = data.obs[cell_type_column].astype(str).to_numpy(); margins = np.zeros(data.n_obs)
        annotation = 'user supplied cell types; not validated by this software'
    else: raise ValueError('Unknown annotation recipe.')
    if not set(labels) <= set(TYPES): raise ValueError(f'Map your annotations to {TYPES}; unsupported: {set(labels) - set(TYPES)}')
    if 'Melanoma' not in labels: raise ValueError('No melanoma-labelled cells; review annotations before making a melanoma simulation.')
    # Uniform selection: preserve the sampling design; never silently add/enrich immune cells.
    idx = np.sort(np.random.default_rng(seed).choice(data.n_obs, min(max_cells, data.n_obs), replace=False))
    if 'Melanoma' not in labels[idx]: raise ValueError('Subsample contains no melanoma; increase --max-cells or choose another seed.')
    if len(idx) < 2: raise ValueError('The tissue prototype requires at least two cells.')
    cells = []
    for i in idx:
        cells.append({'cell_id': str(data.obs_names[i]), 'cell_type': str(labels[i]),
                      'programs': {k: round(float(1 - np.exp(-scores[k][i] / 2)), 6) for k in PROGRAMS},
                      'marker_margin': round(float(margins[i]), 5) if recipe else None})
    known = any(s['sha256'] == MELANOMA['sha256'] for s in sources)
    provenance = {'sources': sources, 'dataset': MELANOMA if known else None,
                  'input_kind': 'public_melanoma' if known else ('synthetic' if str(data.uns.get('source', '')).startswith('Synthetic melanoma') else 'own_data'),
                  'input_cells': input_n, 'biopsy_cells': data.n_obs, 'exported_cells': len(cells),
                  'sample_id': str(sample_id), 'annotation': annotation, 'annotation_recipe': recipe,
                  'full_composition': dict(Counter(labels)), 'sampled_composition': dict(Counter(labels[idx])),
                  'sampling': 'uniform without replacement; no immune enrichment', 'subsample_seed': seed,
                  'count_layer': layer or 'raw counts', 'gene_coverage': coverage,
                  'program_transform': 'mean log1p(CP10K) over fixed panel, including missing genes as zero; 1-exp(-score/2)',
                  'coordinates': 'absent; simulated geometry is a declared hypothesis',
                  'rates_fitted': False, 'source_license': MELANOMA['license'] if known else source_license,
                  'expression_exported': False}
    write_bundle({'schema_version': 1, 'cells': cells, 'provenance': provenance}, output)
    return provenance


def synthetic_melanoma(path, cells=300, seed=21):
    """Offline mixed melanoma fixture for CI. Explicitly generated, not a patient."""
    path = Path(path)
    if path.exists(): raise FileExistsError(path)
    if not 10 <= cells <= 2000 or seed < 0: raise ValueError('Synthetic fixture requires 10–2000 cells and seed >= 0.')
    rng = np.random.default_rng(seed)
    genes = sorted(set(sum(PANELS.values(), []) + sum(PROGRAMS.values(), []) + ['CD8A', 'CD8B', 'ACTB']))
    types = rng.choice(TYPES[:-1], cells, p=[.5, .17, .07, .05, .1, .06, .025, .025])
    types[0] = 'Melanoma'
    x = rng.poisson(.15, (cells, len(genes)))
    for i, t in enumerate(types):
        panel = PANELS['T_cell' if t == 'CD8_T' else t] + (['CD8A','CD8B'] if t == 'CD8_T' else [])
        for gene in panel: x[i, genes.index(gene)] += rng.poisson(8)
        for gene in PROGRAMS['cytotoxicity'] if t in ['CD8_T', 'NK'] else []: x[i, genes.index(gene)] += rng.poisson(6)
        for gene in PROGRAMS['proliferation'] if t == 'Melanoma' and i % 3 == 0 else []: x[i, genes.index(gene)] += rng.poisson(7)
    x[:, genes.index('ACTB')] += 1
    a = ad.AnnData(sparse.csr_matrix(x), obs=pd.DataFrame({'cell_type': types, 'sample_id': 'synthetic_melanoma'}, index=[f'synthetic_{i:05d}' for i in range(cells)]), var=pd.DataFrame(index=genes))
    a.uns['source'] = 'Synthetic melanoma test fixture; MIT; no patient data.'
    path.parent.mkdir(parents=True, exist_ok=True); a.write_h5ad(path, compression='gzip')
    return {'path': str(path), 'cells': cells, 'synthetic': True}
