"""Explicit input contracts, checksum bundles and independent test data."""
from __future__ import annotations

import csv
import gzip
import io
import json
import tarfile
import urllib.request
import zipfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import io as sio, sparse

from . import __version__
from .model import DEMO_MODEL, checked_model, load_model
from .provenance import file_hashes

PBMC = {
    'id': '10x_pbmc3k', 'title': '3k PBMCs from a Healthy Donor',
    'source_url': 'https://www.10xgenomics.com/datasets/3-k-pbm-cs-from-a-healthy-donor-1-standard-1-1-0',
    'download_url': 'https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz',
    'license': 'CC-BY-4.0', 'license_url': 'https://creativecommons.org/licenses/by/4.0/',
    'attribution': '10x Genomics, 3k PBMCs from a Healthy Donor, Cell Ranger 1.1.0, 2016.',
    'sha256': '847d6ebd9a1ec9a768f2be7e40ca42cbfe75ebeb6d76a4c24167041699dc28b5',
    'bytes': 7621991, 'license_verified_date': '2026-10-27',
    'use': 'Input/annotation/software demonstration; not tumour or treatment validation.',
}


def json_bytes(value):
    return (json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n').encode()


def table_bytes(rows, columns):
    stream = io.StringIO(newline='')
    writer = csv.DictWriter(stream, fieldnames=columns, delimiter='\t', lineterminator='\n')
    writer.writeheader(); writer.writerows(rows)
    return stream.getvalue().encode()


def matrix_directory(path):
    def find(*names):
        candidates = [path / (name + suffix) for name in names for suffix in ['', '.gz']]
        present = [p for p in candidates if p.is_file()]
        if len(present) != 1:
            raise ValueError(f'Expected exactly one of {names}, optionally .gz, in {path}.')
        return present[0]
    matrix, features, barcodes = find('matrix.mtx'), find('features.tsv', 'genes.tsv'), find('barcodes.tsv')
    def opening(p):
        return gzip.open(p, 'rb') if p.suffix == '.gz' else p.open('rb')
    with opening(matrix) as handle:
        x = sparse.csr_matrix(sio.mmread(handle, spmatrix=True).T)
    with opening(features) as handle:
        genes = pd.read_csv(handle, sep='\t', header=None, dtype=str, keep_default_na=False)
    with opening(barcodes) as handle:
        ids = pd.read_csv(handle, sep='\t', header=None, dtype=str, keep_default_na=False)[0]
    if x.shape != (len(ids), len(genes)):
        raise ValueError('Matrix dimensions do not match barcodes/features. Expected genes × cells on disk.')
    if genes.shape[1] >= 3:
        keep = genes[2].eq('Gene Expression').to_numpy()
        if not keep.any():
            raise ValueError('No Gene Expression features in this 10x matrix.')
        x, genes = x[:, keep], genes.loc[keep]
    if genes[0].duplicated().any():
        raise ValueError('Feature identifiers must be unique.')
    obs = pd.DataFrame(index=pd.Index(ids, name='cell_id'))
    var = pd.DataFrame({'gene_symbols': genes[1].to_numpy() if genes.shape[1] >= 2 else genes[0].to_numpy()},
                       index=pd.Index(genes[0].to_numpy(), name='gene_id'))
    return ad.AnnData(x, obs=obs, var=var), [matrix, features, barcodes]


def load_input(path, metadata=None, layer=None):
    path = Path(path)
    if path.is_dir():
        if layer:
            raise ValueError('--layer applies only to H5AD inputs.')
        data, files = matrix_directory(path)
    elif path.suffix == '.h5':
        if layer:
            raise ValueError('--layer applies only to H5AD inputs.')
        import h5py
        with h5py.File(path, 'r') as handle:
            if 'matrix' not in handle:
                raise ValueError('Expected a Cell Ranger v3+ feature-barcode H5 matrix.')
            m = handle['matrix']
            x = sparse.csc_matrix((m['data'][:], m['indices'][:], m['indptr'][:]), shape=tuple(m['shape'][:])).T.tocsr()
            decode = lambda a: [v.decode() if isinstance(v, bytes) else str(v) for v in a]
            features = m['features']
            var = pd.DataFrame({'gene_symbols': decode(features['name'][:]),
                                'feature_type': decode(features['feature_type'][:])},
                               index=decode(features['id'][:]))
            data = ad.AnnData(x, obs=pd.DataFrame(index=decode(m['barcodes'][:])), var=var)
            data = data[:, data.var.feature_type.eq('Gene Expression').to_numpy()].copy()
        files = [path]
    elif path.suffix == '.h5ad':
        data, files = ad.read_h5ad(path), [path]
        if layer:
            if layer not in data.layers:
                raise ValueError(f'Missing layer {layer!r}; available: {list(data.layers)}')
            data.X = data.layers[layer].copy()
    else:
        raise ValueError('Input must be an H5AD file, Cell Ranger H5 file or a 10x Matrix Market directory.')
    if not data.obs_names.is_unique or not data.var_names.is_unique:
        raise ValueError('Cell and feature IDs must be unique; do not silently rename duplicates.')
    if data.n_obs == 0 or data.n_vars == 0:
        raise ValueError('Empty expression matrix.')
    data.obs_names = data.obs_names.astype(str)
    if metadata:
        table = pd.read_csv(metadata, sep='\t', dtype=str, keep_default_na=False)
        if 'cell_id' not in table or table.cell_id.duplicated().any():
            raise ValueError('Metadata needs unique cell_id values.')
        table = table.set_index('cell_id')
        if set(table.index) != set(data.obs_names):
            raise ValueError('Metadata and expression must have exactly the same cell IDs; no positional/intersection join.')
        # Explicit external metadata replaces, rather than ambiguously merges, obs.
        data.obs = table.loc[data.obs_names].copy()
        files.append(Path(metadata))
    if data.X is None:
        raise ValueError('The selected count matrix is empty; choose a raw-count layer.')
    values = data.X.data if sparse.issparse(data.X) else np.asarray(data.X).ravel()
    if not np.isfinite(values).all() or (values < 0).any() or not np.allclose(values, np.rint(values), atol=1e-8, rtol=0):
        raise ValueError('Provide raw nonnegative integer counts; choose --layer counts for normalized H5AD objects.')
    totals = np.asarray(data.X.sum(axis=1)).ravel()
    if not np.isfinite(totals).all() or (totals <= 0).any():
        raise ValueError('Every input cell must have positive total counts; review upstream QC.')
    sources = [{'name': p.name, 'bytes': p.stat().st_size, 'sha256': file_hashes(p)['sha256']} for p in files]
    return data, sources


def demo_pbmc_labels(data):
    """Fixed transparent heuristic only; never a validated reference classifier."""
    names = data.var.get('gene_symbols', pd.Series(data.var_names, index=data.var_names)).astype(str).to_numpy()
    def counts(genes):
        mask = np.isin(names, genes)
        if not mask.any():
            raise ValueError(f'PBMC recipe needs gene symbols from {genes}.')
        return np.asarray(data.X[:, mask].sum(axis=1)).ravel()
    selected = (counts(['CD3D', 'CD3E']) > 0) & (counts(['CD8A', 'CD8B']) > 0)
    programs = [counts(['CCR7', 'IL7R', 'TCF7']) / 3,
                counts(['GZMK']), counts(['NKG7', 'GNLY', 'PRF1', 'GZMB']) / 4]
    labels = np.asarray(DEMO_MODEL['states'])[np.argmax(np.stack(programs, axis=1), axis=1)]
    return selected, labels


def prepare(input_path, output, *, metadata=None, layer=None, model_path=None,
            state_column='state', sample_column='sample_id', sample_id=None,
            recipe=None, max_cells=None, seed=21, source_license='user-supplied; redistribution rights not assessed'):
    output = Path(output)
    if output.exists():
        raise FileExistsError(f'Output exists: {output}. Choose a new bundle name.')
    model = load_model(model_path)
    data, sources = load_input(input_path, metadata, layer)
    original_n = data.n_obs
    if sample_column in data.obs:
        values = data.obs[sample_column]
        if values.isna().any() or values.astype(str).str.strip().eq('').any():
            raise ValueError('Sample identifiers cannot be missing.')
        available = values.astype(str).unique().tolist()
        if sample_id is None and len(available) != 1:
            raise ValueError('Multiple samples detected: choose --sample-id explicitly; never pool patients/timepoints.')
        sample_id = sample_id or available[0]
        if sample_id not in available:
            raise ValueError('Requested sample ID is absent.')
        data = data[values.astype(str).eq(sample_id).to_numpy()].copy()
    else:
        sample_id = sample_id or 'user_declared_single_sample'
    for identity in ['patient_id', 'timepoint']:
        if identity in data.obs and data.obs[identity].nunique(dropna=False) > 1:
            raise ValueError(f'Selected input contains multiple {identity} values; prepare one sample/timepoint at a time.')
    annotation = 'user_supplied_states'
    if recipe:
        if recipe != 'pbmc-cd8-demo' or model['states'] != DEMO_MODEL['states']:
            raise ValueError('The PBMC demo recipe requires the default demo states.')
        mask, labels = demo_pbmc_labels(data)
        data.obs['state'] = labels
        data = data[mask].copy()
        state_column = 'state'
        annotation = 'illustrative_marker_recipe_not_validated; CD3D/E>0 and CD8A/B>0; mean-count argmax with fixed tie order'
    if state_column not in data.obs:
        raise ValueError(f'Missing annotation column {state_column!r}. Supply reviewed states in H5AD obs or --metadata.')
    values = data.obs[state_column]
    if values.isna().any() or not set(values.astype(str)) <= set(model['states']):
        raise ValueError(f'All selected cells must have a state in {model["states"]}.')
    if data.n_obs == 0:
        raise ValueError('No cells selected; revise explicit sample/annotation choices.')
    available_n = data.n_obs
    if max_cells is not None:
        if max_cells < 1 or seed < 0:
            raise ValueError('max_cells must be positive and seed nonnegative.')
        if data.n_obs > max_cells:
            selection = np.sort(np.random.default_rng(seed).choice(data.n_obs, max_cells, replace=False))
            data = data[selection].copy()
    rows = [{'cell_id': str(i), 'state': str(row[state_column]), 'sample_id': str(sample_id)}
            for i, row in data.obs.iterrows()]
    if any(not r['cell_id'].strip() or len(r['cell_id']) > 200 for r in rows):
        raise ValueError('Cell IDs must be nonempty and at most 200 characters.')
    provenance = {'schema_version': 1, 'software_version': __version__, 'sources': sources,
                  'input_cells': original_n, 'selected_before_subsampling': available_n,
                  'exported_cells': len(rows), 'features': data.n_vars, 'sample_id': sample_id,
                  'count_layer': layer or 'X / Matrix Market raw counts', 'annotation': annotation,
                  'subsample_seed': seed if len(rows) < available_n else None,
                  'source_license': source_license, 'expression_exported': False,
                  'rates_inferred_from_expression': False, 'coordinates_inferred': False,
                  'purpose': 'Uncalibrated software experiment; no clinical prediction claim.'}
    payload = {'cells.tsv': table_bytes(rows, ['cell_id', 'state', 'sample_id']),
               'model.yaml': yaml.safe_dump(model, sort_keys=False).encode(),
               'provenance.json': json_bytes(provenance)}
    import hashlib
    manifest = {'schema_version': 1, 'kind': 'vtsc_input_bundle',
                'files': {k: hashlib.sha256(v).hexdigest() for k, v in payload.items()}}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + '.partial')
    with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for key, value in {**payload, 'manifest.json': json_bytes(manifest)}.items():
            entry = zipfile.ZipInfo(key, date_time=(2026, 1, 1, 0, 0, 0))
            entry.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(entry, value)
    temporary.replace(output)
    return provenance


def read_bundle(path):
    import hashlib
    expected = {'cells.tsv', 'model.yaml', 'provenance.json', 'manifest.json'}
    with zipfile.ZipFile(path) as archive:
        infos = archive.infolist()
        if len(infos) != 4 or {i.filename for i in infos} != expected or sum(i.file_size for i in infos) > 50_000_000:
            raise ValueError('Invalid bundle members or uncompressed bundle exceeds 50 MB.')
        payload = {i.filename: archive.read(i) for i in infos}
    manifest = json.loads(payload['manifest.json'])
    if manifest.get('schema_version') != 1 or manifest.get('kind') != 'vtsc_input_bundle' or set(manifest.get('files', {})) != expected - {'manifest.json'}:
        raise ValueError('Unsupported bundle manifest.')
    for name, checksum in manifest['files'].items():
        if hashlib.sha256(payload[name]).hexdigest() != checksum:
            raise ValueError(f'Bundle checksum mismatch: {name}')
    model = checked_model(yaml.safe_load(payload['model.yaml']))
    reader = csv.DictReader(io.StringIO(payload['cells.tsv'].decode()), delimiter='\t')
    if reader.fieldnames != ['cell_id', 'state', 'sample_id']:
        raise ValueError('Invalid cells table columns.')
    cells = list(reader)
    if not cells or len({c['cell_id'] for c in cells}) != len(cells) or any(not c['cell_id'] or c['state'] not in model['states'] for c in cells):
        raise ValueError('Invalid cell IDs or states in bundle.')
    if len({c['sample_id'] for c in cells}) != 1:
        raise ValueError('Bundle must contain exactly one sample.')
    return cells, model, json.loads(payload['provenance.json'])


def synthetic_example(output, *, cells=300, seed=7):
    path = Path(output)
    if path.exists():
        raise FileExistsError(f'Output exists: {path}')
    if not 3 <= cells <= 100000 or seed < 0:
        raise ValueError('Synthetic example requires 3–100000 cells and nonnegative seed.')
    rng = np.random.default_rng(seed)
    genes = ['CD3D', 'CD3E', 'CD8A', 'CD8B', 'CCR7', 'IL7R', 'TCF7', 'GZMK', 'NKG7', 'GNLY', 'PRF1', 'GZMB', 'ACTB', 'GAPDH']
    labels = rng.choice(3, cells, p=[0.5, 0.3, 0.2])
    rates = np.full((cells, len(genes)), 1.0)
    rates[:, :4] = 5
    for state, columns in enumerate([[4, 5, 6], [7], [8, 9, 10, 11]]):
        rates[np.ix_(labels == state, columns)] = 10
    x = rng.poisson(rates)
    x[:, 12] += 1
    obs = pd.DataFrame({'state': np.asarray(DEMO_MODEL['states'])[labels], 'sample_id': 'synthetic_01'},
                       index=[f'synthetic_cell_{i:05d}' for i in range(cells)])
    data = ad.AnnData(sparse.csr_matrix(x), obs=obs, var=pd.DataFrame(index=genes))
    data.uns['source'] = 'Generated by VirtualTissue-scRNAseq; no patient data; code and generated fixture under MIT.'
    path.parent.mkdir(parents=True, exist_ok=True)
    data.write_h5ad(path, compression='gzip')
    return {'path': str(path), 'cells': cells, 'genes': len(genes), 'seed': seed, 'license': 'MIT'}


def download_pbmc(destination):
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    archive_path = destination / 'pbmc3k.tar.gz'
    if not archive_path.exists():
        partial = archive_path.with_suffix('.partial')
        request = urllib.request.Request(PBMC['download_url'], headers={'User-Agent': 'VirtualTissue-scRNAseq/0.2'})
        with urllib.request.urlopen(request, timeout=60) as response, partial.open('wb') as handle:
            received = 0
            while block := response.read(1024 * 1024):
                received += len(block)
                if received > PBMC['bytes']:
                    raise ValueError('Source exceeds expected size; review version before using it.')
                handle.write(block)
        if partial.stat().st_size != PBMC['bytes'] or file_hashes(partial)['sha256'] != PBMC['sha256']:
            raise ValueError('PBMC download checksum mismatch; partial file is not accepted.')
        partial.replace(archive_path)
    if file_hashes(archive_path)['sha256'] != PBMC['sha256']:
        raise ValueError('Cached PBMC source is corrupt; remove it and download again.')
    folder = destination / 'matrix'
    folder.mkdir(exist_ok=True)
    with tarfile.open(archive_path) as archive:
        for name in ['matrix.mtx', 'genes.tsv', 'barcodes.tsv']:
            member = archive.getmember('filtered_gene_bc_matrices/hg19/' + name)
            if not member.isfile() or member.size > 40_000_000:
                raise ValueError('Unexpected PBMC archive member.')
            content = archive.extractfile(member).read()
            target = folder / name
            if target.exists() and target.read_bytes() != content:
                raise ValueError(f'Existing example file differs from verified source: {target}')
            target.write_bytes(content)
    (destination / 'SOURCE.json').write_bytes(json_bytes(PBMC))
    return {'matrix_directory': str(folder), **PBMC}
