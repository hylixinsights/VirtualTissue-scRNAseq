import io
import json
import subprocess
import sys
import zipfile

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import io as sio, sparse

from vtsc.inputs import load_input, prepare, read_bundle, synthetic_example
from vtsc.results import run_bundle, read_run, write_html, compare_runs


def fixture(tmp_path, n=40):
    source = tmp_path / 'source.h5ad'
    synthetic_example(source, cells=n)
    return source


def test_input_schema_rejects_normalized_counts_duplicate_ids_and_metadata_mismatch(tmp_path):
    source = fixture(tmp_path)
    x = ad.read_h5ad(source)
    x.X = x.X.astype(float); x.X.data[0] = .5
    x.write_h5ad(tmp_path / 'normalized.h5ad')
    with pytest.raises(ValueError, match='raw nonnegative integer'): load_input(tmp_path / 'normalized.h5ad')
    x = ad.read_h5ad(source); x.obs_names = ['duplicate'] * x.n_obs
    x.write_h5ad(tmp_path / 'duplicate.h5ad')
    with pytest.raises(ValueError, match='unique'): load_input(tmp_path / 'duplicate.h5ad')
    meta = tmp_path / 'meta.tsv'; meta.write_text('cell_id\tstate\nwrong\tResting-like\n')
    with pytest.raises(ValueError, match='same cell IDs'): load_input(source, meta)


def test_metadata_join_is_by_barcode_not_position_and_mtx_transposes(tmp_path):
    directory = tmp_path / 'mtx'; directory.mkdir()
    sio.mmwrite(directory / 'matrix.mtx', sparse.coo_matrix([[2, 7], [3, 9], [1, 0]]))
    (directory / 'features.tsv').write_text('g1\tA\tGene Expression\ng2\tB\tGene Expression\ng3\tAntibody\tAntibody Capture\n')
    (directory / 'barcodes.tsv').write_text('c1\nc2\n')
    meta = tmp_path / 'meta.tsv'; meta.write_text('cell_id\tstate\nc2\tGZMK-like\nc1\tResting-like\n')
    data, _ = load_input(directory, meta)
    assert data.shape == (2, 2)
    assert data.obs['state'].tolist() == ['Resting-like', 'GZMK-like']
    assert np.array_equal(data.X.toarray(), [[2, 3], [7, 9]])


def test_own_h5ad_counts_layer_and_custom_state_model(tmp_path):
    import yaml
    source = fixture(tmp_path)
    data = ad.read_h5ad(source)
    data.layers['counts'] = data.X.copy()
    data.X = data.X.astype(float) / 10
    data.obs['reviewed_state'] = ['Quiescent', 'Effector'] * 20
    data.write_h5ad(source)
    config = {'schema_version': 1, 'name': 'own-data-fixture', 'time_unit': 'hours',
              'rate_provenance': 'Declared test hypothesis', 'states': ['Quiescent', 'Effector'],
              'transitions': [{'source': 'Quiescent', 'target': 'Effector', 'rate': .2}]}
    model = tmp_path / 'model.yaml'; model.write_text(yaml.safe_dump(config))
    bundle = tmp_path / 'custom.zip'
    prepare(source, bundle, layer='counts', state_column='reviewed_state', model_path=model)
    rows, parsed, provenance = read_bundle(bundle)
    assert len(rows) == 40 and parsed['states'] == config['states']
    assert provenance['count_layer'] == 'counts'
    output = tmp_path / 'custom_run'; run_bundle(bundle, output, hours=3, seed=9)
    run = read_run(output)
    assert run['trajectory'][-1]['Quiescent'] < 20
    assert sum(run['trajectory'][-1][s] for s in config['states']) == 40


def test_multisample_selection_is_required(tmp_path):
    source = fixture(tmp_path)
    data = ad.read_h5ad(source); data.obs['sample_id'] = ['one', 'two'] * 20; data.write_h5ad(source)
    with pytest.raises(ValueError, match='Multiple samples'): prepare(source, tmp_path / 'bad.zip')
    prepare(source, tmp_path / 'good.zip', sample_id='two')
    assert len(read_bundle(tmp_path / 'good.zip')[0]) == 20


def test_bundle_is_repeatable_and_rejects_tampering_and_extra_members(tmp_path):
    source = fixture(tmp_path)
    one, two = tmp_path / 'a.zip', tmp_path / 'b.zip'
    prepare(source, one, max_cells=25, seed=3); prepare(source, two, max_cells=25, seed=3)
    assert one.read_bytes() == two.read_bytes()
    with zipfile.ZipFile(one) as archive: payload = {n: archive.read(n) for n in archive.namelist()}
    payload['cells.tsv'] = payload['cells.tsv'].replace(b'Resting-like', b'Changed', 1)
    with zipfile.ZipFile(tmp_path / 'tampered.zip', 'w') as archive:
        for name, content in payload.items(): archive.writestr(name, content)
    with pytest.raises(ValueError, match='checksum'): read_bundle(tmp_path / 'tampered.zip')
    with zipfile.ZipFile(tmp_path / 'extra.zip', 'w') as archive:
        for name, content in payload.items(): archive.writestr(name, content)
        archive.writestr('../escape', b'x')
    with pytest.raises(ValueError, match='members'): read_bundle(tmp_path / 'extra.zip')


def test_end_to_end_outputs_reopen_control_and_comparison(tmp_path):
    source = fixture(tmp_path)
    bundle = tmp_path / 'input.zip'; prepare(source, bundle)
    output = tmp_path / 'run'
    run_bundle(bundle, output, hours=5)
    result = read_run(output)
    assert result == read_run(tmp_path / 'run.zip')
    assert (output / 'overview.png').stat().st_size > 1000
    assert result['summary']['complete']
    control = tmp_path / 'control'; run_bundle(bundle, control, scenario='no-events', hours=3)
    assert read_run(control)['events'] == []
    assert compare_runs(output, control, tmp_path / 'comparison')['shared_hours'] == 3
    with pytest.raises(FileExistsError): run_bundle(bundle, output)
    (output / 'events.tsv').write_text('corrupt')
    with pytest.raises(ValueError, match='checksum'): read_run(output)


def test_explorer_escapes_source_strings(tmp_path):
    path = tmp_path / 'explorer.html'
    write_html({'label': '</script><script>alert(1)</script>'}, path)
    assert '</script><script>alert(1)' not in path.read_text()
    assert '\\u003c/script\\u003e' in path.read_text()


def test_cli_reports_invalid_input_without_traceback(tmp_path):
    result = subprocess.run([sys.executable, '-m', 'vtsc.cli', 'validate', str(tmp_path / 'missing.h5ad')], capture_output=True, text=True)
    assert result.returncode == 2 and 'Traceback' not in result.stderr
