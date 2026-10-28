from copy import deepcopy
import base64,gzip,json,re
from pathlib import Path
import numpy as np
import pytest
from vtsc.tissue_inputs import synthetic_melanoma,prepare_tissue,read_tissue_bundle
from vtsc.melanoma import DEFAULT_CONFIG,simulate_melanoma
from vtsc.melanoma_results import run_melanoma,compare_melanoma,recover_melanoma,verify_dynamic,write_melanoma_html
from vtsc.results import read_run
from vtsc.inputs import load_input


def test_packaged_assumptions_match_frozen_program_and_graph():
    import hashlib
    import vtsc
    from vtsc.programs import LIBRARY
    from vtsc.dynamics import DEFAULT_MODEL
    root=Path(vtsc.__file__).parent/'resources'
    freeze=json.loads((root/'FREEZE.json').read_text())
    for name,digest in freeze['files'].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest()==digest
    assert json.loads((root/'programs_cd8_v1.json').read_text())==LIBRARY
    assert json.loads((root/'assumed_state_graph_v1.json').read_text())==DEFAULT_MODEL


def test_portable_paired_workflow_recovery_and_html_escaping(tmp_path):
    raw=tmp_path/'source.h5ad';synthetic_melanoma(raw,30)
    bundle=tmp_path/'input.zip';prepare_tissue(raw,bundle,source_license='MIT')
    a=run_melanoma(bundle,tmp_path/'a',hours=2)
    b=run_melanoma(bundle,tmp_path/'b',hours=2,arm='exposure-zero')
    compare_melanoma(tmp_path/'a.zip',tmp_path/'b.zip',tmp_path/'comparison')
    merged=verify_dynamic(read_run(tmp_path/'comparison.zip'));assert len(merged['runs'])==2
    recover_melanoma(tmp_path/'a.zip',tmp_path/'recovered',hours=3)
    assert read_run(tmp_path/'recovered')['frames'][:9]==read_run(tmp_path/'a')['frames']
    with pytest.raises(FileExistsError):run_melanoma(bundle,tmp_path/'a')
    run=read_run(tmp_path/'a');run['provenance']['label']='</script><script>alert(1)</script>'
    write_melanoma_html(run,tmp_path/'escaped.html');html=(tmp_path/'escaped.html').read_text()
    assert 'alert(1)' not in html
    packed=re.search(r'<script id="payload" type="application/octet-stream">([^<]+)</script>',html)[1]
    assert json.loads(gzip.decompress(base64.b64decode(packed)))['provenance']['label']==run['provenance']['label']
    run_melanoma(bundle,tmp_path/'different-seed',hours=2,seed=22)
    with pytest.raises(ValueError,match='same input'):compare_melanoma(tmp_path/'a',tmp_path/'different-seed',tmp_path/'invalid')


def test_arbitrary_schedule_boundaries_do_not_break_checkpoint_negative_control(tmp_path):
    raw=tmp_path/'source.h5ad';synthetic_melanoma(raw,30);bundle=tmp_path/'input.zip';prepare_tissue(raw,bundle)
    cfg=deepcopy(DEFAULT_CONFIG);cfg['model']['checkpoint_enabled']=False
    cfg['model']['schedule']=[{'start_hour':.37,'end_hour':1.19,'exposure':1.}]
    p=read_tissue_bundle(bundle)
    a=simulate_melanoma(p,hours=2,config=cfg);b=simulate_melanoma(p,hours=2,config=cfg,arm='exposure-zero')
    assert [f['cells'] for f in a['frames']]==[f['cells'] for f in b['frames']]
    assert [{k:v for k,v in e.items() if k!='exposure'} for e in a['events']]==[{k:v for k,v in e.items() if k!='exposure'} for e in b['events']]


def test_neutral_adapter_preserves_order_and_whitelists_metadata(tmp_path):
    pytest.importorskip('pyarrow')
    import pandas as pd
    from scipy import io,sparse
    x=np.array([[1,2,0],[4,1,3]])
    io.mmwrite(tmp_path/'counts.mtx',sparse.csr_matrix(x))
    pd.DataFrame({'cell_id':['c','a','b'],'sample_ID':['s']*3,'Response':['secret']*3,'CD8_state':['future']*3}).to_parquet(tmp_path/'obs.parquet')
    pd.DataFrame({'gene_id':['G1','G2']}).to_csv(tmp_path/'var.tsv',sep='\t',index=False)
    data,_=load_input(tmp_path)
    assert data.obs_names.tolist()==['c','a','b'] and np.array_equal(data.X.toarray(),x.T)
    assert set(data.obs.columns)=={'sample_id'}


def test_invalid_own_gene_identifiers_require_mapping(tmp_path):
    import anndata as ad
    raw=tmp_path/'source.h5ad';synthetic_melanoma(raw,30);a=ad.read_h5ad(raw)
    symbols=a.var_names.tolist();a.var_names=[f'ENSG{i:09}' for i in range(a.n_vars)];a.write_h5ad(raw)
    with pytest.raises(ValueError,match='mapping'):prepare_tissue(raw,tmp_path/'bad.zip')
    import pandas as pd
    pd.DataFrame({'gene_id':a.var_names,'gene_symbol':symbols}).to_csv(tmp_path/'map.tsv',sep='\t',index=False)
    prepare_tissue(raw,tmp_path/'good.zip',gene_map=tmp_path/'map.tsv')
