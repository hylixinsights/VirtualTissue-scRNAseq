"""A05 synthetic recovery, identifiability and spatial/compartment consistency."""
from copy import deepcopy
import json,time
from pathlib import Path
import numpy as np
from vtsc.dynamics import DEFAULT_MODEL,STATES,solve
from vtsc.melanoma import DEFAULT_CONFIG,simulate_melanoma
from vtsc.inference import fit_ladder,predict,parameter_value,composition,clr,MEAN_CONTEXT
from vtsc.inputs import json_bytes

root=Path(__file__).resolve().parents[1];out=root/'reports/product-v04';start=time.monotonic()
variants={'V-T':{'transition:Memory:Activated':.07},'V-P':{'division:Activated':.035},'V-I':{'influx:Memory':1.5},
          'V-TP':{'transition:Memory:Activated':.07,'division:Activated':.035},
          'V-TPI':{'transition:Memory:Activated':.07,'division:Activated':.035,'influx:Memory':1.5}}
recoveries=[];datasets=[]
for variant,truth in variants.items():
 for seed in [21,22,23,24]:
    rng=np.random.default_rng(seed);m=deepcopy(DEFAULT_MODEL);m['schedule']=[]
    for edge in m['transitions']:edge['rate_per_hour']=0
    m['division_per_hour']=dict.fromkeys(STATES,0.);m['death_per_hour']=dict.fromkeys(STATES,0.)
    for name,value in truth.items():parameter_value(m,name,value)
    pre=np.array([rng.multinomial(800,rng.dirichlet(np.full(4,2.))) for i in range(24)])
    times=np.tile([6.,12.,24.],8);prob=predict('M2',pre,times,model=m)
    on=np.array([rng.multinomial(800,rng.dirichlet(p*1000)) for p in prob])
    initial=deepcopy(m)
    for name,value in truth.items():parameter_value(initial,name,value*.7)
    fit=fit_ladder(pre,on,times,family='M2',model=initial,parameters=list(truth),variant=variant,
                   concentration=1000,baseline_draws=8,seed=seed)
    recoveries.append({'variant':variant,'seed':seed,'truth':truth,'fit':fit})
    datasets.append({'variant':variant,'seed':seed,'pre':pre.tolist(),'on':on.tolist(),'times':times.tolist()})
    print('Recovery',variant,seed,fit['success'],flush=True)
summary={}
for variant,truth in variants.items():
 summary[variant]={}
 rr=[r for r in recoveries if r['variant']==variant]
 for name,value in truth.items():
    estimates=[r['fit']['estimates'][name]['value'] for r in rr]
    intervals=[r['fit']['estimates'][name]['local_95_interval'] for r in rr]
    summary[variant][name]={'truth':value,'estimates':estimates,'mean_bias':float(np.mean(estimates)-value),
        'relative_mean_bias':float(np.mean(estimates)/value-1),'local_95_coverage':float(np.mean([lo<=value<=hi for lo,hi in intervals])),
        'note':'Four synthetic replicates; insufficient to establish calibrated coverage. Baseline is marginalized; concentration held known.'}
# A deliberately non-identifiable composition-only birth/loss pair.
m=deepcopy(DEFAULT_MODEL);m['schedule']=[];m['transitions_enabled']=False
m['division_per_hour']=dict.fromkeys(STATES,.01);m['death_per_hour']=dict.fromkeys(STATES,.005)
rng=np.random.default_rng(31);pre=np.array([rng.multinomial(1000,rng.dirichlet(np.ones(4))) for i in range(12)])
times=np.full(12,12.);prob=predict('M2',pre,times,model=m);on=np.round(prob*1000).astype(int)
confounded=fit_ladder(pre,on,times,family='M2',model=m,parameters=['division:Activated','death:Activated'],concentration=1000,baseline_draws=4,variant='V-P')
assert confounded['identifiability']['jacobian_rank']<2
ladder={family:fit_ladder(pre,on,times,family=family,model=m,baseline_draws=4) for family in ['M0','M1','M1b','M2','M3']}
# Spatial engine, fixed common context; matched continuous and fixed-step expectations.
p={'schema_version':2,'provenance':{'input_kind':'synthetic'},'cells':[]}
for i in range(81):
 lineage='Melanoma' if i==0 else 'CD8_T'
 p['cells'].append({'cell_id':str(i),'cell_type':lineage,'lineage':lineage,'cd8_state':None if i==0 else STATES[(i-1)%4],
  'sample_id':'synthetic','patient_id':None,'clone_id':None,'parent_id':None,'molecular_source_id':str(i),
  'observed_rna_t0':{},'observed_axes_t0':{},'state_source':'synthetic','state_confidence':'assumed','programs':dict.fromkeys(['cytotoxicity','proliferation','interferon'],.5)})
cfg=deepcopy(DEFAULT_CONFIG);cfg['mean_context']=MEAN_CONTEXT;cfg['model']['schedule']=[]
cfg['model']['influx_cells_per_hour']['Memory']=.5
values=[]
for seed in range(64):
 r=simulate_melanoma(p,hours=12,seed=seed,config=cfg)
 values.append([r['summary']['final_cd8_counts'][s] for s in STATES])
values=np.asarray(values);step_expect=solve([20]*4,[12],cfg['model'],MEAN_CONTEXT,dt=.25)[0];continuous=solve([20]*4,[12],cfg['model'],MEAN_CONTEXT)[0]
se=values.std(axis=0,ddof=1)/np.sqrt(len(values));delta=values.mean(axis=0)-step_expect
consistent=bool(np.all(np.abs(delta)<=4*se+1));assert consistent
convergence={str(dt):float(np.max(np.abs(solve([20]*4,[12],cfg['model'],MEAN_CONTEXT,dt=dt)[0]-continuous))) for dt in [.5,.25,.1,.05]}
assert convergence['0.05']<convergence['0.5']
# Context/geometry gap is reported separately, not called agreement with the homogeneous solver.
spatial=[];localcfg=deepcopy(cfg);localcfg['mean_context']=None
for seed in range(8):
 r=simulate_melanoma(p,hours=12,seed=seed,config=localcfg)
 spatial.append([r['summary']['final_cd8_counts'][s] for s in STATES])
report={'scope':'Synthetic engineering verification, not patient validation','recovery':summary,'fits':recoveries,'ladder_smoke':ladder,
 'nonidentifiable_example':confounded,'mean_field':{'seeds':list(range(64)),'replicate_counts':values.tolist(),
 'mean':values.mean(axis=0).tolist(),'standard_error':se.tolist(),'fixed_step_expectation':step_expect.tolist(),
 'continuous_expectation':continuous.tolist(),'difference':delta.tolist(),'criterion':'Absolute difference <= 4 Monte Carlo SE + 1 cell per state, fixed prospectively','passed':consistent,
 'time_step_convergence_max_count_error':convergence,'spatial_replicates':spatial,'spatial_minus_homogeneous':(np.mean(spatial,axis=0)-continuous).tolist(),
 'scope':'Mean-context mode disables geometry constraints and extra actions; spatial gap includes heterogeneous contact/reactivity and space-limited births, not a pure Jensen gap.'},'seconds':time.monotonic()-start}
(out/'synthetic-design.json').write_bytes(json_bytes(datasets));(out/'dynamics-validation.json').write_bytes(json_bytes(report))
print('Dynamics validation complete',flush=True)
