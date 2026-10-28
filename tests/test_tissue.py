from copy import deepcopy
import json
import math
import zipfile

import anndata as ad
import numpy as np
import pytest
from scipy.spatial.distance import pdist, squareform

from vtsc.tissue import simulate_tissue, decision_menu, DEFAULT
from vtsc.tissue_inputs import synthetic_melanoma, prepare_tissue, read_tissue_bundle, write_bundle
from vtsc.tissue_results import run_tissue, write_tissue_html
from vtsc.results import read_run


def payload(n=40):
    return {'cells': [{'cell_id': f'cell_{i}', 'cell_type': 'Melanoma' if i % 2 else 'CD8_T',
                      'programs': {'cytotoxicity': 1., 'proliferation': 1., 'interferon': .2}}
                     for i in range(n)], 'provenance': {'sample_id': 'test', 'dataset': None}}


def test_recognition_contact_and_fatigue_gate_attack():
    from vtsc.tissue import initialize
    cells, _ = initialize(payload(2)['cells'], np.random.default_rng(1), DEFAULT)
    cells[0].update(x=40., y=40.); cells[1].update(x=60., y=40.)
    assert 'attack' not in decision_menu(cells[0], [1], cells, DEFAULT, True)[0]
    cells[1]['x'] = 53
    assert 'attack' in decision_menu(cells[0], [1], cells, DEFAULT, True)[0]
    assert 'attack' not in decision_menu(cells[0], [1], cells, DEFAULT, False)[0]
    cells[0]['fatigue'] = .9
    assert 'attack' not in decision_menu(cells[0], [1], cells, DEFAULT, True)[0]
    cells[0]['fatigue'] = 0.; cells[1]['alive'] = False
    assert 'attack' not in decision_menu(cells[0], [1], cells, DEFAULT, True)[0]


def test_interactions_decisions_accounting_replay_and_negative_control():
    args = dict(hours=12, seed=9, config={'attack_damage_per_hour': 2, 'cycle_hours': 4})
    a = simulate_tissue(payload(), scenario='declared-recognition', **args)
    b = simulate_tissue(payload(), scenario='recognition-off', **args)
    assert a == simulate_tissue(payload(), scenario='declared-recognition', **args)
    assert a['frames'][0] == b['frames'][0]
    assert a['summary']['events'].get('attack', 0) > 0
    assert a['summary']['events'].get('death', 0) > 0
    assert a['summary']['events'].get('division', 0) > 0
    assert not any(e['event'] in ['attack','death'] for e in b['events'])
    assert a['summary']['final_alive'] == 40 + a['summary']['events'].get('division',0) - a['summary']['events'].get('death',0)
    for e in a['events']:
        if e['event'] == 'attack': assert e['distance_um'] <= e['contact_limit_um'] + 1e-8
    for dec in a['decisions']: assert dec['choice'] in dec['options']
    died = {e['cell']:e['time_hours'] for e in a['events'] if e['event']=='death'}
    for f in a['frames']:
        for c in f['cells']:
            if c[0] in died and f['time_hours'] >= died[c[0]]: assert not c[4]
    for dec in a['decisions']:
        if dec['cell'] in died: assert dec['time_hours'] <= died[dec['cell']]


def test_geometry_excludes_overlap_and_migration_is_speed_bounded():
    r = simulate_tissue(payload(), hours=5, seed=3, scenario='recognition-off')
    side = r['geometry']['width_um']
    for f in r['frames']:
        rows = f['cells']; xy=np.array([[c[1],c[2]] for c in rows]); radii=np.array([c[3] for c in rows])
        distances=squareform(pdist(xy)); np.fill_diagonal(distances,np.inf)
        assert np.all(distances + .004 >= radii[:,None] + radii)
        assert np.all(xy + .002 >= radii[:,None]) and np.all(xy <= side-radii[:,None]+.002)
    for prev, curr in zip(r['frames'], r['frames'][1:]):
        old={c[0]:c for c in prev['cells']}
        for c in curr['cells']:
            if r['catalog'][c[0]]['cell_type']=='CD8_T':
                assert math.dist(c[1:3],old[c[0]][1:3]) <= DEFAULT['immune_speed_um_h']*DEFAULT['dt_hours']+.003
    assert max(r['frames'][-1]['cue']) > 0


def test_frozen_control_and_seed_layout():
    a=simulate_tissue(payload(),scenario='frozen',hours=2)
    assert all(f['cells']==a['frames'][0]['cells'] for f in a['frames'])
    assert not a['events'] and not a['decisions']
    b=simulate_tissue(payload(),scenario='frozen',hours=2,seed=22)
    assert a['frames'][0]['cells'] != b['frames'][0]['cells']
    with pytest.raises(ValueError): simulate_tissue(payload(),hours=float('nan'))
    with pytest.raises(ValueError): simulate_tissue(payload(),config={'dt_hours':0})
    with pytest.raises(ValueError): simulate_tissue(payload(),config={'fake':1})


def test_own_data_tissue_bundle_and_saved_viewer(tmp_path):
    source=tmp_path/'own.h5ad';synthetic_melanoma(source,60)
    data=ad.read_h5ad(source);data.layers['counts']=data.X.copy();data.X=data.X/2
    data.obs = data.obs.rename(columns={'cell_type':'reviewed_type'});data.write_h5ad(source)
    bundle=tmp_path/'input.zip'
    p=prepare_tissue(source,bundle,layer='counts',cell_type_column='reviewed_type',max_cells=40)
    assert p['exported_cells']==40 and p['rates_fitted'] is False
    cells=read_tissue_bundle(bundle)['cells'];assert len(cells)==40
    run_tissue(bundle,tmp_path/'run',hours=1)
    result=read_run(tmp_path/'run.zip'); assert result['kind']=='vtsc_tissue_run'
    assert result['summary']['geometry_measured'] is False
    result['provenance']['sample_id']='</script><script>alert(1)</script>'
    write_tissue_html(result,tmp_path/'test.html')
    html=(tmp_path/'test.html').read_text();assert '</script><script>alert(1)' not in html
    assert 'Agent inspector' in html and 'ASSUMED SPATIAL ARRANGEMENT' in html
    with zipfile.ZipFile(bundle) as z: parts={n:z.read(n) for n in z.namelist()}
    parts['input.json']=parts['input.json'].replace(b'Melanoma',b'Changed')
    with zipfile.ZipFile(tmp_path/'bad.zip','w') as z:
        for k,v in parts.items():z.writestr(k,v)
    with pytest.raises(ValueError,match='checksum'):read_tissue_bundle(tmp_path/'bad.zip')


def test_annotation_does_not_hide_non_tumor_cells_and_multisample_rejected(tmp_path):
    path=tmp_path/'source.h5ad';synthetic_melanoma(path,100)
    p=prepare_tissue(path,tmp_path/'demo.zip',recipe='melanoma-markers',max_cells=100)
    assert p['exported_cells']==100
    assert p['full_composition']==p['sampled_composition']
    assert p['full_composition'].get('Melanoma',0)>0 and p['full_composition'].get('CD8_T',0)>0
    data=ad.read_h5ad(path);data.obs['sample_id']=['a','b']*50;data.write_h5ad(path)
    with pytest.raises(ValueError,match='Multiple samples'): prepare_tissue(path,tmp_path/'bad.zip')
    prepare_tissue(path,tmp_path/'one.zip',sample_id='a',max_cells=40)


def test_recording_cap_reports_observed_duration(monkeypatch):
    import vtsc.tissue as tissue
    monkeypatch.setattr(tissue, 'RECORDING_CELL_LIMIT', 100)
    result=simulate_tissue(payload(40),hours=10,scenario='frozen')
    assert result['summary']['stop_reason']=='recording_cap'
    assert result['summary']['observed_hours']==.25
    assert sum(len(f['cells']) for f in result['frames'])<=100


def test_cellranger_h5_preserves_raw_counts_and_excludes_antibody(tmp_path):
    import h5py
    from scipy import sparse
    from vtsc.inputs import load_input
    path=tmp_path/'10x.h5'
    x=sparse.csc_matrix([[1,2],[3,4],[8,9]])
    with h5py.File(path,'w') as f:
        m=f.create_group('matrix')
        for k,v in {'data':x.data,'indices':x.indices,'indptr':x.indptr,'shape':x.shape,'barcodes':[b'a',b'b']}.items():m.create_dataset(k,data=v)
        genes=m.create_group('features')
        for k,v in {'id':[b'g1',b'g2',b'a1'],'name':[b'MLANA',b'CD3D',b'antibody'],'feature_type':[b'Gene Expression',b'Gene Expression',b'Antibody Capture']}.items():genes.create_dataset(k,data=v)
    a,_=load_input(path)
    assert a.obs_names.tolist()==['a','b']
    assert a.var.gene_symbols.tolist()==['MLANA','CD3D']
    assert np.array_equal(a.X.toarray(),[[1,3],[2,4]])


def test_births_cannot_overflow_recording_budget(monkeypatch):
    import vtsc.tissue as tissue
    monkeypatch.setattr(tissue, 'RECORDING_CELL_LIMIT', 125)
    result=simulate_tissue(payload(40),hours=5,scenario='recognition-off',config={'cycle_hours':.25})
    assert result['summary']['stop_reason']=='recording_cap'
    assert sum(len(f['cells']) for f in result['frames'])<=125
