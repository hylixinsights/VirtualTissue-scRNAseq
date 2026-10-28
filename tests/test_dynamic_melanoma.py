from copy import deepcopy
import json
import numpy as np
import pytest
from scipy.special import gammaln

from vtsc.dynamics import DEFAULT_MODEL, STATES, checked_model, hazards, probabilities, solve, draw_decision, exposure_at, verify_journal, MockProvider
from vtsc.melanoma import DEFAULT_CONFIG, simulate_melanoma
from vtsc.tissue_inputs import synthetic_melanoma,prepare_tissue,read_tissue_bundle
from vtsc.inference import dm_logpmf,fit_ladder,predict

CTX={'contact':1.,'target':1.,'reactive':1.,'tumor_fraction':.5,'myeloid_fraction':.1}


def payload(n=30):
    cells=[]
    for i in range(n):
        lineage='Melanoma' if i==0 else 'CD8_T'
        cells.append({'cell_id':f'c{i}','cell_type':lineage,'lineage':lineage,'cd8_state':None if i==0 else STATES[(i-1)%4],
         'time_in_state':0.,'sample_id':'test','patient_id':None,'clone_id':None,'parent_id':None,'molecular_source_id':f'c{i}',
         'state_source':'synthetic','state_confidence':'assumed','reference_id':'assumed',
         'observed_rna_t0':{'TCF7':i+1},'observed_axes_t0':{},'programs':dict(cytotoxicity=.8,proliferation=.8,interferon=.1)})
    return {'schema_version':2,'cells':cells,'provenance':{'input_kind':'synthetic','source_license':'MIT'}}


def mean_config():
    c=deepcopy(DEFAULT_CONFIG);c['mean_context']=CTX;c['model']['schedule']=[];return c


def test_generator_conserves_mass_and_matches_oneway_analytic():
    m=deepcopy(DEFAULT_MODEL);m['demography_enabled']=False;m['schedule']=[]
    m['transitions']=[dict(m['transitions'][0],rate_per_hour=.3)]
    q=hazards('Memory',CTX,0,m)[0]['transition:Activated']
    t=np.linspace(0,15,16);x=solve([100,0,0,0],t,m,CTX)
    assert np.allclose(x[:,0],100*np.exp(-q*t))
    assert np.allclose(x.sum(axis=1),100)
    coarse=solve([100,0,0,0],[15],m,CTX,dt=.5)
    fine=solve([100,0,0,0],[15],m,CTX,dt=.05)
    assert np.linalg.norm(fine-x[-1])<=np.linalg.norm(coarse-x[-1])+1e-10


def test_competing_probabilities_extremes_and_input_rejection():
    for rates in [{},{'a':0,'b':0},{'a':1e9,'b':2e9},{'a':.1,'b':.2}]:
        p=probabilities(rates,.25)
        assert sum(p.values())==pytest.approx(1)
        assert all(0<=v<=1 for v in p.values())
    with pytest.raises(ValueError): probabilities({'bad':float('nan')},.25)
    m=deepcopy(DEFAULT_MODEL);m['schedule']=[{'start_hour':3,'end_hour':2,'exposure':1}]
    with pytest.raises(ValueError):checked_model(m)
    assert exposure_at(2,DEFAULT_MODEL)==.8 and exposure_at(96,DEFAULT_MODEL)==0


def test_drug_is_decoupled_without_checkpoint_and_recognition_is_separate():
    m=deepcopy(DEFAULT_MODEL)
    h0,s0=hazards('Memory',CTX,0,m);h1,s1=hazards('Memory',CTX,1,m)
    assert s1['antigen_signal']>s0['antigen_signal']>0
    assert h1['transition:Activated']>h0['transition:Activated']
    m['checkpoint_enabled']=False
    assert hazards('Memory',CTX,0,m)[0]==hazards('Memory',CTX,1,m)[0]
    for key in ['contact','target','reactive']:
        ctx={**CTX,key:0}
        assert hazards('Memory',ctx,1,m)[0]['transition:Activated']==0
    assert DEFAULT_MODEL['profiles']['Dysfunctional']['proliferation']>0
    assert not any(e['source']=='Dysfunctional' and e['target']=='Memory' for e in DEFAULT_MODEL['transitions'])


def test_runtime_transitions_really_change_axes_preserve_rna_and_journal():
    p=payload(40);before=deepcopy(p);c=mean_config()
    for e in c['model']['transitions']:e['rate_per_hour']=1.
    r=simulate_melanoma(p,hours=4,config=c)
    assert p==before and r['summary']['events']['transition']>0
    assert verify_journal(r['decisions'])
    for e in r['events']:
        if e['event']=='transition':
            d=next(d for d in r['decisions'] if d['cell']==e['cell'] and d['time_hours']==e['time_hours'])
            assert d['axes']==c['model']['profiles'][e['to_state']]
            assert d['state_decision']['time_in_state_after']==0
            assert r['catalog'][e['cell']]['lineage']=='CD8_T'
    r['decisions'][0]['choice']='tampered';assert not verify_journal(r['decisions'])


def test_provider_only_chooses_legal_state_transitions():
    class Illegal:
        name='bad'
        def choose(self,r):return 'transition:Myeloid'
    m=deepcopy(DEFAULT_MODEL);m['demography_enabled']=False
    with pytest.raises(ValueError,match='illegal'):
        draw_decision('Memory',0,CTX,0,.25,m,21,'cell',1,Illegal())
    d=draw_decision('Memory',0,CTX,0,.25,m,21,'cell',1,MockProvider())
    assert d['choice']=='stay' and d['provider']=='mock-stay-v1'


def test_actual_state_change_increases_contact_damage_without_changing_rna(monkeypatch):
    import vtsc.melanoma as engine
    initialize=engine.initialize
    def contacting(*args):
        cells,side=initialize(*args)
        cells[0].update(x=40.,y=40.)
        cells[1].update(x=40.+cells[0]['radius']+cells[1]['radius']+.5,y=40.)
        return cells,side
    monkeypatch.setattr(engine,'initialize',contacting)
    class Activate:
        name='test-legal-activation'
        def choose(self,request):return 'transition:Activated'
    cfg=deepcopy(DEFAULT_CONFIG);cfg['model']['demography_enabled']=False
    cfg['model']['reactive_fraction']=1.;cfg['model']['checkpoint_enabled']=False
    p=payload(2);before=deepcopy(p)
    stay=simulate_melanoma(p,hours=.25,config=cfg,provider=MockProvider())
    changed=simulate_melanoma(p,hours=.25,config=cfg,provider=Activate())
    damage=lambda run:sum(e['damage'] for e in run['events'] if e['event']=='attack')
    assert damage(changed)>damage(stay)>0
    assert p==before and changed['frames'][-1]['cells'][1][10]=='Activated'
    assert changed['catalog'][1]['observed_rna_t0']==stay['catalog'][1]['observed_rna_t0']


def test_control_identity_zero_drug_frozen_transition_off_and_dead_irreversible():
    p=payload();cfg=mean_config()
    a=simulate_melanoma(p,hours=3,config=cfg,arm='anti-pd1')
    b=simulate_melanoma(p,hours=3,config=cfg,arm='exposure-zero')
    assert a['frames']==b['frames'] and a['decisions']==b['decisions']
    c=deepcopy(cfg);c['model']['checkpoint_enabled']=False
    c['model']['schedule']=[{'start_hour':0,'end_hour':96,'exposure':1}]
    a=simulate_melanoma(p,hours=3,config=c);b=simulate_melanoma(p,hours=3,config=c,arm='exposure-zero')
    assert [f['cells'] for f in a['frames']]==[f['cells'] for f in b['frames']]
    frozen=simulate_melanoma(p,hours=3,arm='frozen');assert all(f['cells']==frozen['frames'][0]['cells'] for f in frozen['frames'])
    off=simulate_melanoma(p,hours=6,arm='transitions-off',config=cfg)
    assert not any(e['event']=='transition' for e in off['events'])
    cfg['model']['death_per_hour']=dict.fromkeys(STATES,2.)
    r=simulate_melanoma(p,hours=3,config=cfg)
    assert r['summary']['events']['death']>0
    for event in r['events']:
        if event['event']=='death':
            assert all(not row[4] for f in r['frames'] if f['time_hours']>=event['time_hours'] for row in f['cells'] if row[0]==event['cell'])


def test_division_influx_accounting_and_unknown_not_invented():
    c=mean_config();c['model']['division_per_hour']=dict.fromkeys(STATES,1.);c['model']['influx_cells_per_hour']['Memory']=2
    p=payload(10);p['cells'][1]['cd8_state']='Unknown'
    r=simulate_melanoma(p,hours=3,config=c)
    ev=r['summary']['events'];assert ev['division']>0 and ev['influx']>0
    assert r['summary']['final_alive']==10+ev.get('division',0)+ev.get('influx',0)-ev.get('death',0)
    assert all(f['cells'][1][10]=='Unknown' for f in r['frames'])
    for cat in r['catalog'][10:]:
        if cat['parent_id']:
            parent=next(x for x in r['catalog'] if x['cell_id']==cat['parent_id'])
            assert cat['observed_rna_t0']==parent['observed_rna_t0'] and cat['molecular_source_id']==parent['molecular_source_id']
        else:assert cat['observed_rna_t0']=={} and cat['molecular_source_id'] is None


def test_spatial_limits_recognition_off_and_paired_initial_positions():
    p=payload(25)
    a=simulate_melanoma(p,hours=3);b=simulate_melanoma(p,hours=3,arm='recognition-off')
    assert a['initial_signature']==b['initial_signature']
    assert not any(e['event']=='attack' for e in b['events'])
    side=a['geometry']['width_um']
    for f in a['frames']:
        for row in f['cells']:assert row[3]-1e-4<=row[1]<=side-row[3]+1e-4 and row[3]-1e-4<=row[2]<=side-row[3]+1e-4
        for i,row in enumerate(f['cells']):
            for other in f['cells'][i+1:]:assert np.hypot(row[1]-other[1],row[2]-other[2])+2e-4>=row[3]+other[3]
    for d in a['decisions']:
        if d['choice']=='attack':assert d['target'] in d['context']['target_ids'] and d['context']['recognition']>0
    for e in a['events']:
        if e['event']=='attack':assert e['distance_um']<=e['contact_limit_um']+1e-9


def test_common_kernel_mean_ensemble_not_a_single_seed():
    p=payload(81);cfg=mean_config();cfg['model']['demography_enabled']=False
    results=[]
    for seed in range(16):
        r=simulate_melanoma(p,hours=5,seed=seed,config=cfg)
        results.append([r['summary']['final_cd8_counts'][s] for s in STATES])
    results=np.array(results);expected=solve([20]*4,[5],cfg['model'],CTX,dt=.25)[0]
    # Monte Carlo mean must lie within 4 SE + one cell; tolerance fixed before execution.
    assert np.all(abs(results.mean(axis=0)-expected)<=4*results.std(axis=0,ddof=1)/np.sqrt(len(results))+1)


def test_dm_m0_shift_estimation_and_leakage_guard():
    assert dm_logpmf([0,0,0,0],[.25]*4,10)==pytest.approx(0)
    n=np.array([[100,50,20,10],[80,50,20,20],[100,30,20,10]])
    assert np.allclose(predict('M0',n,[1,2,3]),n/n.sum(axis=1,keepdims=True))
    on=n[:,::-1].copy();fit=fit_ladder(n,on,[1,2,3],family='M1',baseline_draws=4)
    assert fit['success'] and np.isclose(sum(fit['predicted_compositions'][0]),1)
    assert fit['clinical_validation'] is False
    with pytest.raises(ValueError,match='Leakage'):
        fit_ladder(n,on,[1,2,3],patient_ids=['a','b','c'],forbidden_ids=['a'])
    with pytest.raises(ValueError,match='intervals'):fit_ladder(n,on,[1,np.nan,3])


def test_input_v2_raw_layer_and_bad_states(tmp_path):
    import anndata as ad
    f=tmp_path/'cells.h5ad';synthetic_melanoma(f,40)
    a=ad.read_h5ad(f);a.raw=a.copy();a.X=a.X.astype(float);a.X.data=np.log1p(a.X.data);a.write_h5ad(f)
    prepare_tissue(f,tmp_path/'good.zip',layer='raw');p=read_tissue_bundle(tmp_path/'good.zip');assert p['schema_version']==2
    a=ad.read_h5ad(f);a.obs['cd8_state']=a.obs['cd8_state'].astype(str);a.obs.loc[a.obs.cell_type=='CD8_T','cd8_state']='Myeloid';a.write_h5ad(f)
    with pytest.raises(ValueError,match='CD8 annotations'):prepare_tissue(f,tmp_path/'bad.zip',layer='raw')
