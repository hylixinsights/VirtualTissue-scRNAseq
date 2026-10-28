"""Reproducible release integration (public 10x already downloaded; no network)."""
from pathlib import Path
import json,sys,time
import anndata as ad
import numpy as np
from vtsc.inputs import json_bytes
from vtsc.tissue_inputs import prepare_tissue,synthetic_melanoma,read_tissue_bundle
from vtsc.melanoma_results import melanoma_demo,run_melanoma,compare_melanoma,recover_melanoma,verify_dynamic,write_melanoma_html
from vtsc.results import read_run
from vtsc.provenance import file_hashes
root=Path(__file__).resolve().parents[1]
out=root/(sys.argv[1] if len(sys.argv)>1 else 'outputs/v04-release')
if out.exists():raise FileExistsError(out)
out.mkdir(parents=True)
start=time.monotonic();report={'scope':'Implementation integration; no clinical validation','output':str(out)}
report['synthetic']=melanoma_demo(out/'synthetic',cells=120,hours=24,seed=21)
print('Synthetic complete',flush=True)
source=root/'example_data/melanoma10x/filtered_feature_bc_matrix.h5'
bundle=out/'real-input.zip'
report['real_input']=prepare_tissue(source,bundle,recipe='melanoma-markers',max_cells=400,seed=21)
(root/'reports/product-v04/melanoma-input-audit.json').write_bytes(json_bytes(report['real_input']))
print('Real input prepared',flush=True)
report['real_runs']={arm:run_melanoma(bundle,out/('real-'+arm),arm=arm,hours=24,seed=21) for arm in ['anti-pd1','exposure-zero','recognition-off']}
report['real_comparison']=compare_melanoma(out/'real-anti-pd1',out/'real-exposure-zero',out/'real-comparison')
reopened=verify_dynamic(read_run(out/'real-anti-pd1.zip'));write_melanoma_html(reopened,out/'real-reopened.html')
report['reopened_zip_verified']=True
print('Real arms and ZIP complete',flush=True)
# Own-data route: counts layer, exact barcode join, reviewed lineages/states, two biopsies require explicit selection.
fixture=out/'own.h5ad';synthetic_melanoma(fixture,60,22);a=ad.read_h5ad(fixture)
a.layers['counts']=a.X.copy();a.X=a.X.astype(float);a.X.data=np.log1p(a.X.data)
a.obs['sample_id']=['own_a' if i%2 else 'own_b' for i in range(a.n_obs)]
a.obs['patient_id']='own_fixture';a.write_h5ad(fixture)
meta=a.obs.copy();meta.index.name='cell_id';meta.iloc[::-1].to_csv(out/'own-metadata.tsv',sep='\t')
report['own_input']=prepare_tissue(fixture,out/'own-input.zip',metadata=out/'own-metadata.tsv',layer='counts',sample_id='own_b',max_cells=30,source_license='MIT; synthetic own-data interface fixture')
report['own_run']=run_melanoma(out/'own-input.zip',out/'own-run',hours=4,seed=21)
report['recovery']=recover_melanoma(out/'own-run.zip',out/'own-recovered',hours=6)
errors={}
for label,kw in [('normalized',{'sample_id':'own_b'}),('mixed_samples',{'layer':'counts'})]:
 try:prepare_tissue(fixture,out/f'{label}.zip',**kw)
 except ValueError as e:errors[label]=str(e)
 else:raise AssertionError(label+' was accepted')
report['invalid_input_errors']=errors
report['seconds']=time.monotonic()-start
report['artifact_hashes']={str(p.relative_to(out)):file_hashes(p)['sha256'] for p in out.rglob('summary.json')}
(root/'reports/product-v04/integration.json').write_bytes(json_bytes(report))
print('All integration workflows complete',flush=True)
