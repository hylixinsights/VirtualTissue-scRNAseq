"""Portable v2 runs, matched controls, deterministic recovery and ensembles."""
from copy import deepcopy
import hashlib
import gzip
import base64
import json
import os
from pathlib import Path
import platform
import shutil
import sys
import time

import numpy as np
import scipy
import yaml
from . import __version__
from .dynamics import STATES, verify_journal
from .inputs import json_bytes, table_bytes
from .melanoma import ARMS, simulate_melanoma, validate_config
from .jev_provider import JevProvider
from .jev_payload import PARAMETERS, PROMPT_VERSION, mock_transport
from .behavior import BehaviorController
from .provenance import file_hashes
from .results import read_run
from .tissue_inputs import prepare_tissue, read_tissue_bundle, synthetic_melanoma


def compact(value):return json.dumps(value,separators=(',',':'),allow_nan=False).encode()


def implementation_hashes():
    return {name:file_hashes(Path(__file__).parent/name)['sha256'] for name in
            ['dynamics.py','melanoma.py','contract.py','programs.py','tissue.py','behavior.py','jev_provider.py','jev_payload.py']}


def write_melanoma_html(result,path):
    template=(Path(__file__).parent/'web/melanoma.html').read_text()
    data=base64.b64encode(gzip.compress(compact(result),mtime=0)).decode('ascii')
    path=Path(path)
    if path.exists():raise FileExistsError(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(template.replace('__PAYLOAD__',data))


def verify_dynamic(result):
    runs=result.get('runs',[result])
    for r in runs:
        if r.get('kind')!='vtsc_melanoma_run' or r.get('schema_version')!=2:raise ValueError('Not a dynamic melanoma v2 recording.')
        if not verify_journal(r['decisions']):raise ValueError('Decision journal hash chain failed verification.')
        if r['decisions'] and r['decisions'][-1]['record_hash']!=r['summary']['journal_root']:raise ValueError('Journal root mismatch.')
        if r.get('provider_journal'):
            from .jev_payload import fingerprint
            previous='0'*64
            for index,entry in enumerate(r['provider_journal']):
                value={k:v for k,v in entry.items() if k!='digest'}
                if entry['sequence']!=index or entry['previous']!=previous or fingerprint(value)!=entry['digest']:
                    raise ValueError('Provider journal hash chain failed verification.')
                previous=entry['digest']
            if previous!=r['behavior_policy']['journal_root']:raise ValueError('Provider journal root mismatch.')
    return result


def save_result(result,path):
    path=Path(path)
    if path.exists() or Path(str(path)+'.zip').exists():raise FileExistsError('Use a new run name; accepted results are never overwritten.')
    verify_dynamic(result);path.mkdir(parents=True)
    (path/'run.json').write_bytes(compact(result))
    (path/'summary.json').write_bytes(json_bytes(result['summary']))
    (path/'provenance.json').write_bytes(json_bytes(result['provenance']))
    (path/'assumptions.yaml').write_text(yaml.safe_dump(result['config'],sort_keys=False))
    (path/'decisions.jsonl').write_bytes(b''.join(compact(d)+b'\n' for d in result['decisions']))
    (path/'events.jsonl').write_bytes(b''.join(compact(e)+b'\n' for e in result['events']))
    if 'provider_journal' in result:
        (path/'provider_journal.jsonl').write_bytes(b''.join(compact(e)+b'\n' for e in result['provider_journal']))
    rows=[]
    for f in result['frames']:
        rows.append({'time_hours':f['time_hours'],**f['counts'],**{f'CD8_{k}':v for k,v in f['cd8_counts'].items()},'exposure':f['exposure']})
    (path/'trajectory.tsv').write_bytes(table_bytes(rows,list(rows[0])))
    write_melanoma_html(result,path/'explorer.html')
    manifest={'schema_version':1,'kind':'vtsc_results','files':{p.name:file_hashes(p)['sha256'] for p in sorted(path.iterdir())}}
    (path/'manifest.json').write_bytes(json_bytes(manifest));shutil.make_archive(str(path),'zip',root_dir=path)
    return {'output':str(path),'explorer':str(path/'explorer.html'),'results_zip':str(path)+'.zip',**result['summary']}


def run_melanoma(bundle,output,*,config_path=None,decision_mode='local',max_http_attempts=None,
                 allow_paid=False,jev_model='jev-1.14.0',replay_journal=None,**options):
    path=Path(output)
    if path.exists() or Path(str(path)+'.zip').exists():raise FileExistsError('Use a new run name.')
    payload=read_tissue_bundle(bundle);cfg=validate_config(yaml.safe_load(Path(config_path).read_text()) if config_path else None)
    if decision_mode not in ['local','mock-jev','jev']:raise ValueError('Unknown behavior decision mode.')
    if decision_mode=='local' and replay_journal:raise ValueError('A provider journal requires mock-jev or jev mode.')
    behavior=None;binding=None
    if decision_mode!='local':
        if cfg['mean_context'] is not None or cfg['policy']!='hazards':raise ValueError('Use spatial mode and kernel stochastic state hazards with Jev behaviors.')
        if decision_mode=='jev' and (not allow_paid or max_http_attempts is None):
            raise ValueError('Live Jev requires --allow-paid and an explicit --max-http-attempts ceiling. Preflight/mock uses no paid calls.')
        key=(os.environ.get('TYPESAFE_API_KEY') or os.environ.get('JEV_API_KEY')) if decision_mode=='jev' else ''
        if decision_mode=='jev' and not key:raise ValueError('Set TYPESAFE_API_KEY or JEV_API_KEY in the environment; never put credentials in config.')
        binding={'bundle_sha256':file_hashes(bundle)['sha256'],'config':cfg,'seed':options.get('seed',21),
                 'arm':options.get('arm','anti-pd1'),'layout':options.get('layout','tumor-nests'),
                 'max_cells':options.get('max_cells',1500),'implementation_hashes':implementation_hashes(),'prompt_version':PROMPT_VERSION,'mode':decision_mode}
        journal=Path(str(path)+'.provider.jsonl')
        if journal.exists():raise FileExistsError('Provider journal already exists. Preserve it; explicitly replay it into a new run name.')
        journal.parent.mkdir(parents=True,exist_ok=True)
        if replay_journal:shutil.copyfile(replay_journal,journal)
        api=JevProvider(key,jev_model,PARAMETERS,max_http_attempts if max_http_attempts is not None else 500,journal,binding,
                        transport=mock_transport if decision_mode=='mock-jev' else None,pause_s=0 if decision_mode=='mock-jev' else .25)
        behavior=BehaviorController(api,decision_mode)
    start=time.monotonic();result=simulate_melanoma(payload,config=cfg,behavior=behavior,**options)
    if behavior is not None:
        result.update(provider_journal=[json.loads(line) for line in behavior.provider.journal.read_text().splitlines()] if behavior.provider.journal.exists() else [],behavior_context=binding)
        result['summary'].update(behavior_mode=decision_mode,behavior_provider=behavior.provider.name,http_attempts=behavior.provider.requests)
    result.update(implementation_hashes=implementation_hashes(),initial_input=payload,software_version=__version__,bundle_sha256=file_hashes(bundle)['sha256'],
                  environment={'python':sys.version.split()[0],'numpy':np.__version__,'scipy':scipy.__version__,'platform':platform.system()},
                  execution={'seconds':time.monotonic()-start,'options':options})
    return save_result(result,output)


def compare_melanoma(first,second,output):
    a,b=verify_dynamic(read_run(first)),verify_dynamic(read_run(second))
    if a['bundle_sha256']!=b['bundle_sha256'] or a['summary']['seed']!=b['summary']['seed'] or a['layout']!=b['layout'] or a['config']!=b['config']:
        raise ValueError('Matched comparison requires the same input, seed, layout and base contract.')
    for key in ['mode','model','prompt_version','schedule']:
        if a.get('behavior_policy',{}).get(key)!=b.get('behavior_policy',{}).get(key):raise ValueError('Matched comparison requires the same behavior policy/model/schedule.')
    if a['scenario']!='spatial-shuffle' and b['scenario']!='spatial-shuffle' and a['initial_signature']!=b['initial_signature']:
        raise ValueError('Initial positions/states/reactivity differ.')
    end=min(a['summary']['observed_hours'],b['summary']['observed_hours'])
    times_a=[f['time_hours'] for f in a['frames'] if f['time_hours']<=end];times_b=[f['time_hours'] for f in b['frames'] if f['time_hours']<=end]
    # Different schedules may split a step. The viewer uses the latest recorded frame <= common clock.
    diffs={k:[a['effective_model'][k],b['effective_model'][k]] for k in a['effective_model'] if a['effective_model'][k]!=b['effective_model'][k]}
    result={'kind':'vtsc_melanoma_comparison','schema_version':2,'runs':[a,b],
      'comparison':{'same_input':True,'same_seed':True,'same_initial_positions':a['initial_signature']==b['initial_signature'],
       'shared_hours':end,'configuration_differences':diffs,'arm_a':a['scenario'],'arm_b':b['scenario'],
       'interpretation':'Paired simulated scenarios. Exposure-zero is a hypothesis-dependent counterfactual, not an observed untreated biopsy.'}}
    path=Path(output)
    if path.exists() or Path(str(path)+'.zip').exists():raise FileExistsError(path)
    path.mkdir(parents=True);(path/'run.json').write_bytes(compact(result));(path/'comparison.json').write_bytes(json_bytes(result['comparison']))
    write_melanoma_html(result,path/'explorer.html')
    (path/'manifest.json').write_bytes(json_bytes({'schema_version':1,'kind':'vtsc_results','files':{p.name:file_hashes(p)['sha256'] for p in sorted(path.iterdir())}}))
    shutil.make_archive(str(path),'zip',root_dir=path)
    return {'output':str(path),'explorer':str(path/'explorer.html'),**result['comparison']}


def recover_melanoma(source,output,*,hours=None):
    """Deterministic recomputation into a fresh directory, then optional extension."""
    old=verify_dynamic(read_run(source))
    if old['kind']!='vtsc_melanoma_run':raise ValueError('Recover an individual run.')
    if old.get('implementation_hashes') and old['implementation_hashes']!=implementation_hashes():
        raise ValueError('Recovery requires the recorded implementation; source hashes differ.')
    requested=old['summary']['requested_hours'] if hours is None else hours
    if requested<old['summary']['observed_hours']:raise ValueError('Recovery cannot shorten an accepted recording.')
    opts=old['execution']['options'].copy();opts.update(hours=requested,seed=old['summary']['seed'],arm=old['scenario'],layout=old['layout'])
    behavior=None
    if old.get('behavior_policy',{}).get('mode','local')!='local':
        mode=old['behavior_policy']['mode'];journal=Path(str(output)+'.provider.jsonl')
        if journal.exists() or Path(output).exists():raise FileExistsError('Use a new output and provider journal name.')
        journal.parent.mkdir(parents=True,exist_ok=True)
        journal.write_bytes(b''.join(compact(e)+b'\n' for e in old['provider_journal']))
        api=JevProvider('',old['behavior_policy']['model'],PARAMETERS,max(1,old['behavior_policy']['usage']['http_attempts']),journal,
                        old['behavior_context'],transport=mock_transport if mode=='mock-jev' else None,pause_s=0,replay_only=True)
        behavior=BehaviorController(api,mode)
    result=simulate_melanoma(old['initial_input'],config=old['config'],behavior=behavior,**opts)
    if behavior is not None:
        result.update(provider_journal=old['provider_journal'],behavior_context=old['behavior_context'])
        result['summary'].update(behavior_mode=mode,behavior_provider=behavior.provider.name,http_attempts=behavior.provider.requests)
    prefix=[f for f in result['frames'] if f['time_hours']<=old['summary']['observed_hours']]
    if prefix!=old['frames']:raise ValueError('Recovery prefix differs. Match runtime, environment and full-step horizon; no old outputs were modified.')
    decisions=[d for d in result['decisions'] if d['time_hours']<=old['summary']['observed_hours']]
    if decisions!=old['decisions']:raise ValueError('Recovery decision prefix differs; no old outputs were modified.')
    result.update(implementation_hashes=implementation_hashes(),initial_input=old['initial_input'],software_version=__version__,bundle_sha256=old['bundle_sha256'],environment=old['environment'],execution={'options':opts,'recovered_from':str(source),'method':'deterministic full replay; no hidden checkpoint or provider requests'})
    return save_result(result,output)


def ensemble_melanoma(bundle,output,*,seeds=(21,22,23,24,25),hours=24,config_path=None,cells_cap=1500):
    path=Path(output)
    if path.exists():raise FileExistsError(path)
    seeds=list(seeds)
    if len(seeds)<2 or len(set(seeds))!=len(seeds):raise ValueError('Use at least two distinct paired seeds.')
    payload=read_tissue_bundle(bundle);cfg=validate_config(yaml.safe_load(Path(config_path).read_text()) if config_path else None)
    rows=[]
    experiments=[('controls',None,1.,arm) for arm in ARMS if arm!='frozen']
    for name in ['reactive_fraction','checkpoint_strength','transition_scale']:
        for factor in [.5,2.]:
            for arm in ['anti-pd1','exposure-zero']:experiments.append(('sensitivity',name,factor,arm))
    for kind,name,factor,arm in experiments:
        config=deepcopy(cfg)
        if name=='transition_scale':
            for edge in config['model']['transitions']:edge['rate_per_hour']*=factor
        elif name:config['model'][name]=min(1.,config['model'][name]*factor)
        for seed in seeds:
            run=simulate_melanoma(payload,hours=hours,seed=seed,arm=arm,config=config,max_cells=cells_cap)
            rows.append({'experiment':kind,'parameter':name,'factor':factor,'arm':arm,'seed':seed,
                         'initial_signature':run['initial_signature'],'final_tumor':run['frames'][-1]['counts']['Melanoma'],
                         'final_cd8':run['frames'][-1]['cd8_counts'],'summary':run['summary']})
    deltas=[]
    for kind,name,factor,arm in experiments:
        if arm!='anti-pd1':continue
        aa=[r for r in rows if (r['experiment'],r['parameter'],r['factor'],r['arm'])==(kind,name,factor,'anti-pd1')]
        bb=[r for r in rows if (r['experiment'],r['parameter'],r['factor'],r['arm'])==(kind,name,factor,'exposure-zero')]
        values=[a['final_tumor']-b['final_tumor'] for a,b in zip(aa,bb)]
        deltas.append({'parameter':name,'factor':factor,'paired_tumor_differences':values,'mean':float(np.mean(values)),
                       'quantiles_10_50_90':np.quantile(values,[.1,.5,.9]).tolist()})
    path.mkdir(parents=True);result={'bundle_sha256':file_hashes(bundle)['sha256'],'seeds':seeds,'rows':rows,'paired_differences':deltas,
      'scope':'Monte Carlo variability and assumed-parameter sensitivity; not patient confidence intervals or clinical benefit',
      'config':cfg,'saturation':'Fractions/strength capped at 1 for factor 2; inspect effective value', 'all_completed':all(r['summary']['complete'] for r in rows)}
    (path/'ensemble.json').write_bytes(json_bytes(result));return {'output':str(path),'runs':len(rows),'all_completed':result['all_completed']}


def melanoma_demo(output,*,cells=120,hours=24,seed=21,decision_mode='local'):
    path=Path(output)
    if path.exists():raise FileExistsError(path)
    synthetic_melanoma(path/'synthetic_melanoma.h5ad',cells,seed);bundle=path/'input.zip'
    prepare_tissue(path/'synthetic_melanoma.h5ad',bundle,max_cells=cells,seed=seed,source_license='MIT; generated synthetic fixture')
    runs={arm:run_melanoma(bundle,path/arm,arm=arm,hours=hours,seed=seed,decision_mode=decision_mode) for arm in ['anti-pd1','exposure-zero']}
    runs['comparison']=compare_melanoma(path/'anti-pd1',path/'exposure-zero',path/'comparison')
    return runs


def preflight_melanoma(bundle,output,*,config_path=None,seed=21,arm='anti-pd1',layout='tumor-nests',jev_model='jev-1.14.0'):
    """Exact first-round payload preview through the real scheduler; zero HTTP calls."""
    import tempfile
    out=Path(output)
    if out.exists():raise FileExistsError('Use a new preflight name.')
    p=read_tissue_bundle(bundle);cfg=validate_config(yaml.safe_load(Path(config_path).read_text()) if config_path else None)
    bodies=[]
    def capture(body):bodies.append(deepcopy(body));return mock_transport(body)
    with tempfile.TemporaryDirectory(prefix='vtsc-preflight-') as temporary:
        api=JevProvider('',jev_model,PARAMETERS,100000,Path(temporary)/'preview.jsonl',{'preflight':True},transport=capture,pause_s=0)
        controller=BehaviorController(api,'mock-jev')
        run=simulate_melanoma(p,hours=cfg['mechanics']['dt_hours'],seed=seed,arm=arm,layout=layout,config=cfg,
                              max_cells=max(1500,len(p['cells'])),behavior=controller)
    report={'kind':'zero-cost behavior preflight','actual_http_requests':0,'model':jev_model,'prompt_version':PROMPT_VERSION,
            'first_round_questions':sum(len(b['questions']) for b in bodies),'first_round_batches':len(bodies),
            'first_round_payloads':bodies,'initial_cells':len(p['cells']),'assumed_geometry':layout,
            'schedule':controller.metadata()['schedule'],'limits':PARAMETERS,
            'note':'Future rounds depend on state/contact/growth. A request cap is not a currency budget or a guarantee of reaching the horizon.',
            'observed_programs_are_not_aucell':'vtsc rank/log-count summaries; no equivalence to Xenium depth-conditioned GO AUCell claimed'}
    out.mkdir(parents=True);(out/'preflight.json').write_bytes(json_bytes(report))
    return {'output':str(out),'preflight':str(out/'preflight.json'),'first_round_questions':report['first_round_questions'],
            'first_round_batches':report['first_round_batches'],'actual_http_requests':0}
