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
from .programs import LIBRARY, score_programs, classify, audit_programs, gene_names
from .dynamics import STATES

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
        scores[key] = np.asarray(x[:, mask].sum(axis=1)).ravel() / max(1,len(coverage[key]))
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
    manifest = {'kind': 'vtsc_tissue_input', 'schema_version': payload.get('schema_version', 1),
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
    if manifest.get('schema_version') not in [1,2] or manifest != {'kind': 'vtsc_tissue_input', 'schema_version': manifest.get('schema_version'), 'files': {'input.json': hashlib.sha256(content).hexdigest()}}:
        raise ValueError('Tissue input checksum or manifest mismatch.')
    payload = json.loads(content); cells = payload['cells']
    if payload.get('schema_version') != manifest['schema_version']: raise ValueError('Input schema differs from manifest.')
    if not cells or len({c['cell_id'] for c in cells}) != len(cells): raise ValueError('Invalid cell identifiers.')
    for c in cells:
        if c['cell_type'] not in TYPES or not c['cell_id'].strip() or len(c['cell_id']) > 200:
            raise ValueError('Invalid cell type or identifier.')
        if set(c['programs']) != set(PROGRAMS) or any(not np.isfinite(v) or not 0 <= v <= 1 for v in c['programs'].values()):
            raise ValueError('Invalid expression program values.')
    if payload['schema_version'] == 2:
        for c in cells:
            if c.get('lineage') != c['cell_type'] or c.get('cd8_state') not in (STATES+['Unknown'] if c['lineage']=='CD8_T' else [None]):
                raise ValueError('Invalid lineage/CD8 state combination.')
        if len({c['sample_id'] for c in cells}) != 1: raise ValueError('Bundle contains multiple biopsies.')
    return payload


def prepare_tissue(input_path, output, *, metadata=None, layer=None, cell_type_column='cell_type',
                   sample_id=None, recipe=None, max_cells=400, seed=21, source_license='user supplied',
                   state_column='cd8_state', state_recipe=None, gene_map=None, sampling='stratified'):
    if not 2 <= max_cells <= 2000 or seed < 0: raise ValueError('Use 2–2000 cells and a nonnegative seed.')
    data, sources = load_input(input_path, metadata, layer)
    if gene_map:
        mapping = pd.read_csv(gene_map, sep='\t', dtype=str, keep_default_na=False)
        if set(mapping.columns) != {'gene_id','gene_symbol'} or mapping.gene_id.duplicated().any():
            raise ValueError('Gene mapping needs unique gene_id and gene_symbol columns.')
        mapping=mapping.set_index('gene_id').gene_symbol
        if not set(data.var_names) <= set(mapping.index): raise ValueError('Gene mapping must cover every expression feature.')
        data.var['gene_symbols']=mapping.loc[data.var_names].to_numpy()
        sources.append({'name':Path(gene_map).name,'bytes':Path(gene_map).stat().st_size,'sha256':file_hashes(gene_map)['sha256']})
    gene_names(data)
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
    rank_scores, rank_detected, rank_coverage = score_programs(data)
    inferred_states, state_margins = classify(rank_scores, rank_detected)
    cd8 = labels == 'CD8_T'
    states = np.full(data.n_obs, None, dtype=object)
    if state_column in data.obs and recipe is None and state_recipe is None:
        given = data.obs[state_column]
        if given[cd8].isna().any() or not set(given[cd8].astype(str)) <= set(STATES+['Unknown']):
            raise ValueError('CD8 annotations must be Memory, Activated, Cytotoxic, Dysfunctional or Unknown.')
        noncd8 = given[~cd8].dropna().astype(str)
        if (~noncd8.isin(['','Unknown','None'])).any(): raise ValueError('CD8 states cannot be assigned to another lineage.')
        states[cd8]=given[cd8].astype(str)
        state_source='user supplied; not validated'
    elif state_recipe == 'rank-markers-v1' or recipe == 'melanoma-markers':
        states[cd8]=inferred_states[cd8]; state_source=LIBRARY['version']+' provisional inference'
    elif cd8.any():
        raise ValueError('Supply reviewed cd8_state or explicitly select --state-recipe rank-markers-v1.')
    else: state_source='no CD8 cells'
    if state_recipe not in [None,'rank-markers-v1']: raise ValueError('Unknown state recipe.')
    n=min(max_cells,data.n_obs); rng=np.random.default_rng(seed)
    if sampling == 'uniform': idx=np.sort(rng.choice(data.n_obs,n,replace=False))
    elif sampling == 'stratified':
        groups=sorted(set(labels)); sizes=np.array([sum(labels==g) for g in groups]); quota=n*sizes/sizes.sum()
        counts=np.floor(quota).astype(int)
        for j in np.argsort(-(quota-counts),kind='stable')[:n-counts.sum()]: counts[j]+=1
        idx=np.sort(np.concatenate([rng.choice(np.flatnonzero(labels==g),count,replace=False) for g,count in zip(groups,counts)]))
    else: raise ValueError('Sampling must be stratified or uniform.')
    if 'Melanoma' not in labels[idx]: raise ValueError('Subsample contains no melanoma; increase --max-cells.')
    if len(idx)<2: raise ValueError('The tissue prototype requires at least two cells.')
    names=gene_names(data); tracked=sorted(set(sum(LIBRARY['state'].values(),[])+sum(LIBRARY['axes'].values(),[])+['PDCD1','CD274']))
    x=sparse.csr_matrix(data.X)
    cells=[]
    def optional(i,k):
        v=data.obs.iloc[i].get(k)
        return None if pd.isna(v) or v=='' else str(v)
    for i in idx:
        observed={g:int(x[i,np.flatnonzero(names==g)].sum()) for g in tracked if g in names}
        cells.append({'cell_id':str(data.obs_names[i]),'cell_type':str(labels[i]),'lineage':str(labels[i]),
          'cd8_state':states[i], 'time_in_state':0., 'sample_id':str(sample_id), 'patient_id':optional(i,'patient_id'),
          'clone_id':optional(i,'clone_id'), 'parent_id':None, 'molecular_source_id':str(data.obs_names[i]),
          'reference_id':'assumed-capacity-v1:'+str(states[i]), 'state_source':state_source,
          'state_margin':float(state_margins[i]) if 'provisional' in state_source and cd8[i] else None,
          'state_confidence':'unknown' if states[i]=='Unknown' else 'uncalibrated',
          'observed_rna_t0':observed,
          'observed_axes_t0':{k:float(rank_scores['axes:'+k][i]) if np.isfinite(rank_scores['axes:'+k][i]) else None for k in LIBRARY['axes']},
          'programs':{k:round(float(1-np.exp(-scores[k][i]/2)),6) for k in PROGRAMS},
          'marker_margin':round(float(margins[i]),5) if recipe else None})
    program_audit=audit_programs(data,labels,rank_scores,rank_detected,rank_coverage,states)
    known = any(s['sha256'] == MELANOMA['sha256'] for s in sources)
    provenance = {'sources': sources, 'dataset': MELANOMA if known else None,
                  'input_kind': 'public_melanoma' if known else ('synthetic' if str(data.uns.get('source', '')).startswith('Synthetic melanoma') else 'own_data'),
                  'input_cells': input_n, 'biopsy_cells': data.n_obs, 'exported_cells': len(cells),
                  'sample_id': str(sample_id), 'annotation': annotation, 'annotation_recipe': recipe,
                  'full_composition': dict(Counter(labels)), 'sampled_composition': dict(Counter(labels[idx])),
                  'sampling': sampling+' without replacement; proportional largest-remainder quotas by lineage when stratified; no enrichment', 'subsample_seed': seed,
                  'count_layer': layer or 'raw counts', 'gene_coverage': coverage,
                  'program_transform': 'mean log1p(CP10K) over available panel genes; coverage explicit; 1-exp(-score/2); dynamic CD8 axes use separately documented capped ranks',
                  'coordinates': 'absent; simulated geometry is a declared hypothesis',
                  'program_audit':program_audit, 'state_source':state_source,
                  'baseline_environment':{'tumor_fraction':float(np.mean(labels=='Melanoma')), 'myeloid_fraction':float(np.mean(labels=='Myeloid')), 'denominator':'all input cells of selected biopsy; provisional lineage labels', 'site':None, 'chemistry':None},
                  'rates_fitted': False, 'source_license': MELANOMA['license'] if known else source_license,
                  'expression_exported': False}
    write_bundle({'schema_version': 2, 'cells': cells, 'provenance': provenance}, output)
    return provenance


def synthetic_melanoma(path, cells=300, seed=21):
    """Offline mixed melanoma fixture for CI. Explicitly generated, not a patient."""
    path = Path(path)
    if path.exists(): raise FileExistsError(path)
    if not 10 <= cells <= 2000 or seed < 0: raise ValueError('Synthetic fixture requires 10–2000 cells and seed >= 0.')
    rng = np.random.default_rng(seed)
    genes = sorted(set(sum(PANELS.values(), []) + sum(PROGRAMS.values(), []) + sum(LIBRARY['state'].values(),[]) + sum(LIBRARY['held_out'].values(),[]) + sum(LIBRARY['axes'].values(),[]) + sum(LIBRARY['controls'].values(),[]) + ['CD8A', 'CD8B', 'ACTB']))
    types = rng.choice(TYPES[:-1], cells, p=[.5, .17, .07, .05, .1, .06, .025, .025])
    types[0] = 'Melanoma'
    x = rng.poisson(.15, (cells, len(genes)))
    for i, t in enumerate(types):
        panel = PANELS['T_cell' if t == 'CD8_T' else t] + (['CD8A','CD8B'] if t == 'CD8_T' else [])
        for gene in panel: x[i, genes.index(gene)] += rng.poisson(8)
        for gene in PROGRAMS['cytotoxicity'] if t in ['CD8_T', 'NK'] else []: x[i, genes.index(gene)] += rng.poisson(6)
        for gene in PROGRAMS['proliferation'] if t == 'Melanoma' and i % 3 == 0 else []: x[i, genes.index(gene)] += rng.poisson(7)
    cd8_states=np.full(cells,'',dtype=object)
    for i,t in enumerate(types):
        if t=='CD8_T':
            state=STATES[i%4]; cd8_states[i]=state
            for gene in LIBRARY['state'][state]+LIBRARY['held_out'][state]: x[i,genes.index(gene)]+=rng.poisson(12)
    x[:, genes.index('ACTB')] += 1
    a = ad.AnnData(sparse.csr_matrix(x), obs=pd.DataFrame({'cell_type': types, 'cd8_state':cd8_states, 'sample_id': 'synthetic_melanoma'}, index=[f'synthetic_{i:05d}' for i in range(cells)]), var=pd.DataFrame(index=genes))
    a.uns['source'] = 'Synthetic melanoma test fixture; MIT; no patient data.'
    path.parent.mkdir(parents=True, exist_ok=True); a.write_h5ad(path, compression='gzip')
    return {'path': str(path), 'cells': cells, 'synthetic': True}
