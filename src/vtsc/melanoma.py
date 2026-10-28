"""Dynamic-state melanoma runtime v2. Shares CD8 hazards with the compartment solver."""
from collections import Counter
from copy import deepcopy
import hashlib
import json
import math

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

from .dynamics import (STATES, DEFAULT_MODEL, checked_model, exposure_at, signals, draw_decision,
                       keyed_rng, chain_journal, MockProvider, LocalPolicy, context_values)
from .tissue import DEFAULT, RADII, MOBILE, initialize, free_disc, checked_config
from .tissue_inputs import TYPES
from .contract import parameter_ledger
from .behavior import options_for
from .jev_provider import BudgetExhausted, ProviderError

ARMS=['anti-pd1','exposure-zero','recognition-off','checkpoint-off','transitions-off','spatial-shuffle','frozen']
FRAME_COLUMNS=['index','x','y','radius','alive','action','target','damage','fatigue','cycle',
               'cd8_state','time_in_state','dysfunction','proliferation','cytotoxicity','persistence']
DEFAULT_CONFIG={'schema_version':2,'mechanics':DEFAULT,'model':DEFAULT_MODEL,'mean_context':None,
                'policy':'hazards','recording_cell_limit':100000}


def validate_config(config=None):
    cfg=deepcopy(DEFAULT_CONFIG if config is None else config)
    if set(cfg)!=set(DEFAULT_CONFIG) or cfg['schema_version']!=2: raise ValueError('Use the complete v2 melanoma-template configuration.')
    cfg['mechanics']=checked_config(cfg['mechanics']);cfg['model']=checked_model(cfg['model'])
    if cfg['policy'] not in ['hazards','local','mock']:raise ValueError('Only offline hazards, local and mock policies are available.')
    if cfg['mean_context'] is not None:context_values(cfg['mean_context'])
    limit=cfg['recording_cell_limit']
    if not isinstance(limit,int) or isinstance(limit,bool) or not 100<=limit<=200000:raise ValueError('Recording limit must be 100–200000 cell frames.')
    return cfg


def arm_model(model,arm):
    m=deepcopy(model)
    if arm=='exposure-zero':
        for window in m['schedule']:window['exposure']=0.
    if arm=='recognition-off':m['recognition_enabled']=False
    if arm=='checkpoint-off':m['checkpoint_enabled']=False
    if arm=='transitions-off':m['transitions_enabled']=False
    return m


def simulate_melanoma(payload, *, hours=24., seed=21, layout='tumor-nests', arm='anti-pd1', config=None,max_cells=1500,provider=None,behavior=None):
    cfg=validate_config(config);p=cfg['mechanics'];m=arm_model(cfg['model'],arm)
    if arm not in ARMS:raise ValueError('Unknown melanoma arm.')
    if payload.get('schema_version')!=2:raise ValueError('Dynamic runtime requires input schema v2. Reprepare raw counts with CD8 states; legacy inputs are not silently migrated.')
    if not math.isfinite(hours) or not 0<hours<=96 or not isinstance(seed,int) or isinstance(seed,bool) or seed<0:raise ValueError('Use hours in (0,96] and a nonnegative integer seed.')
    if not 2<=len(payload['cells'])<=2000 or not len(payload['cells'])<=max_cells<=5000:raise ValueError('Initial 2–2000 cells; cap between initial size and 5000.')
    ids=[c['cell_id'] for c in payload['cells']]
    if len(set(ids))!=len(ids):raise ValueError('Duplicate cell identities.')
    for c in payload['cells']:
        if c['lineage']!=c['cell_type'] or c['cd8_state'] not in (STATES+['Unknown'] if c['lineage']=='CD8_T' else [None]):raise ValueError('Invalid lineage/CD8 state.')
    initial,side=initialize(payload['cells'],keyed_rng(seed,'geometry'),p,layout)
    cells=deepcopy(initial); used_ids=set(ids)
    if arm=='spatial-shuffle':
        for radius in sorted(set(c['radius'] for c in cells)):
            indices=[i for i,c in enumerate(cells) if c['radius']==radius]
            positions=[(cells[i]['x'],cells[i]['y']) for i in indices]
            order=keyed_rng(seed,'spatial-shuffle',radius).permutation(len(indices))
            for k,i in enumerate(indices):cells[i]['x'],cells[i]['y']=positions[order[k]]
    mean_mode=cfg['mean_context'] is not None
    if behavior is not None and (mean_mode or cfg['policy']!='hazards' or provider is not None):
        raise ValueError('Behavior providers require spatial mode and kernel-owned stochastic CD8 hazards.')
    policy=provider or {'hazards':None,'local':LocalPolicy(),'mock':MockProvider()}[cfg['policy']]
    for c in cells:
        c['reactive']=bool(keyed_rng(seed,c['cell_id'],'reactivity').random()<m['reactive_fraction'])
        c['time_in_state']=0.;c['axes']=deepcopy(m['profiles'].get(c['cd8_state'],m['profiles']['Unknown']))
        c['parent_id']=None
    init_signature=hashlib.sha256(json.dumps([[c['cell_id'],c['x'],c['y'],c['cd8_state'],c['reactive']] for c in cells],sort_keys=True).encode()).hexdigest()
    ngrid=32;cue=np.zeros((ngrid,ngrid));dx=side/ngrid
    frames=[];decisions=[];events=[];t=0.;step=0;recorded=0;stopped='completed'
    def event(kind,i,target=-1,**extra):
        events.append({'time_hours':round(t,8),'event':kind,'cell':i,'target':target,**extra})
    def snapshot():
        nonlocal recorded
        rows=[]
        for i,c in enumerate(cells):
            if not c['present']:continue
            rows.append([i,round(c['x'],5),round(c['y'],5),round(c['radius'],5),int(c['alive']),c['action'],c['target'],
                         round(c['damage'],6),round(c['fatigue'],6),round(c['cycle'],6),c['cd8_state'],round(c['time_in_state'],8),
                         *[c['axes'][key] for key in ['dysfunction','proliferation','cytotoxicity','persistence']]])
        recorded+=len(rows)
        frames.append({'time_hours':round(t,8),'cells':rows,'cue':np.round(cue,4).ravel().tolist(),
                       'exposure':exposure_at(t,m),'counts':{k:sum(c['alive'] and c['lineage']==k for c in cells) for k in TYPES},
                       'cd8_counts':{k:sum(c['alive'] and c['lineage']=='CD8_T' and c['cd8_state']==k for c in cells) for k in STATES+['Unknown']}})
    def unique_id(prefix,i):
        value=f'{prefix}_{step}_{i}_{len(cells)}'
        while value in used_ids:value+='x'
        used_ids.add(value);return value
    def make_child(c,i,kind='division'):
        nonlocal stopped
        if len(cells)>=max_cells:stopped='population_cap';return 'population_cap'
        if recorded+len(cells)+1>cfg['recording_cell_limit']:stopped='recording_cap';return 'recording_cap'
        r=c['base_radius']; rng=keyed_rng(seed,c['cell_id'],step,kind); success=False
        if mean_mode:a=b=np.array([c['x'],c['y']]);success=True
        else:
            for angle in rng.uniform(0,2*math.pi,24):
                offset=r*1.01*np.array([math.cos(angle),math.sin(angle)])
                a=np.array([c['x'],c['y']])-offset;b=np.array([c['x'],c['y']])+offset
                if free_disc(a,r,cells,side,(i,)) and free_disc(b,r,cells,side,(i,)):success=True;break
        if not success:return 'blocked_by_space'
        child=deepcopy(c);j=len(cells)
        child.update(cell_id=unique_id('daughter',i),index=j,x=float(b[0]),y=float(b[1]),radius=r,cycle=0.,
                     parent=i,parent_id=c['cell_id'],born_at=t,action='born',target=-1,time_in_state=0.)
        c.update(x=float(a[0]),y=float(a[1]),radius=r,cycle=0.)
        cells.append(child);event('division',i,j,state=c['cd8_state'],parent_id=c['cell_id'],child_id=child['cell_id'])
        return 'completed'
    def kill(c,cause,i):
        c.update(alive=False,death_time=t,action='dead',target=-1);event('death',i,cause=cause,state=c['cd8_state'])
    snapshot()
    while t<hours-1e-10:
        if recorded+len(cells)>cfg['recording_cell_limit']:stopped='recording_cap';break
        dt=min(p['dt_hours'],hours-t)
        boundaries=[w[k] for w in m['schedule'] for k in ['start_hour','end_hour'] if w[k]>t+1e-10]
        if boundaries:dt=min(dt,min(boundaries)-t)
        if behavior is not None and arm!='frozen':
            try:behavior.prepare(cells,p,m,exposure_at(t,m),t)
            except BudgetExhausted:stopped='provider_budget_or_replay_exhausted';break
            except ProviderError:stopped='provider_failed';break
        start=t;exposure=exposure_at(start,m);t+=dt;step+=1
        if arm=='frozen':snapshot();continue
        if not mean_mode:
            for c in cells:
                if c['alive'] and c['lineage']=='Melanoma':
                    y,x=np.clip(np.array([c['y'],c['x']])/dx,0,ngrid-1).astype(int);cue[y,x]+=p['cue_source_per_hour']*dt
                if c['present'] and not c['alive'] and t-c['death_time']>=p['corpse_hours']:
                    c['present']=False;event('clearance_timeout',c['index'])
            cue=gaussian_filter(cue,math.sqrt(2*p['cue_diffusion_um2_h']*dt)/dx,mode='reflect')*math.exp(-p['cue_decay_per_hour']*dt)
        tree=cKDTree([[c['x'],c['y']] for c in cells]);step_n=len(cells)
        for i0 in keyed_rng(seed,step,'order').permutation(step_n):
            i=int(i0);c=cells[i]
            if not c['alive']:continue
            neighbors=[j for j in tree.query_ball_point([c['x'],c['y']],p['perception_um']+p['immune_speed_um_h']*dt) if j!=i and cells[j]['present'] and math.hypot(c['x']-cells[j]['x'],c['y']-cells[j]['y'])<=p['perception_um']]
            # Daughters born during this step also participate as physical obstacles/targets.
            neighbors += [j for j in range(step_n,len(cells)) if math.hypot(c['x']-cells[j]['x'],c['y']-cells[j]['y'])<=p['perception_um']]
            contacts=[j for j in neighbors if math.hypot(c['x']-cells[j]['x'],c['y']-cells[j]['y'])<=c['radius']+cells[j]['radius']+p['contact_gap_um']]
            targets=[j for j in contacts if cells[j]['alive'] and cells[j]['lineage']=='Melanoma']
            corpses=[j for j in contacts if not cells[j]['alive']]
            live=[j for j in neighbors if cells[j]['alive']]
            context=deepcopy(cfg['mean_context']) if mean_mode else {'contact':float(bool(targets)), 'target':float(bool(targets)),
               'reactive':float(c['reactive']), 'tumor_fraction':sum(cells[j]['lineage']=='Melanoma' for j in live)/max(1,len(live)),
               'myeloid_fraction':sum(cells[j]['lineage']=='Myeloid' for j in live)/max(1,len(live))}
            fate=None;fate_outcome='not_cd8'
            if c['lineage']=='CD8_T':
                fate=draw_decision(c['cd8_state'],c['time_in_state'],context,exposure,dt,m,seed,c['cell_id'],step,policy)
                fate_outcome='completed';choice=fate['choice'];c['time_in_state']+=dt
                if choice.startswith('transition:'):
                    previous=c['cd8_state'];c['cd8_state']=fate['new_state'];c['time_in_state']=0.
                    c['axes']=deepcopy(m['profiles'][c['cd8_state']]);c['reference_id']='assumed-capacity-v1:'+c['cd8_state']
                    event('transition',i,from_state=previous,to_state=c['cd8_state'])
                elif choice=='death':kill(c,'state_hazard',i)
                elif choice=='division':fate_outcome=make_child(c,i)
                fate['outcome']=fate_outcome;fate['time_in_state_after']=c['time_in_state'];fate['axes_after']=deepcopy(c['axes'])
            sig=signals(c['cd8_state'],context,exposure,m)
            behavior_decision=None
            action='rest';outcome='completed';target=-1;distance=None;options=['rest'];reason='stationary lineage or fatigue recovery'
            if not c['alive']:action='dead';options=[];reason='irreversible death'
            elif mean_mode:action='observe';options=['observe'];reason='mean-context consistency mode: extra mechanics disabled'
            else:
                if c['lineage'] in MOBILE and c['fatigue']<.85:options.append('migrate')
                # Unknown-state CD8 cannot acquire a made-up effector capacity.
                effector=c['axes']['cytotoxicity'] if c['lineage']=='CD8_T' else c['programs']['cytotoxicity']
                if c['lineage'] in ['CD8_T','NK'] and targets and sig['recognition']>0 and effector>0 and c['fatigue']<.85:options.append('attack')
                if c['lineage']=='Myeloid' and corpses:options.append('clear')
                if c['lineage']=='Melanoma':options.append('grow' if c['cycle']<1 else 'divide')
                if 'clear' in options:action='clear';reason='corpse within contact range'
                elif 'attack' in options:action='attack';reason='reactive agent + assumed antigen + contact + current effector capacity + checkpoint modulation'
                elif 'divide' in options:action='divide';reason='tumor cycle complete; placement checked'
                elif 'grow' in options:action='grow';reason='assumed tumor cycle rate modulated by initial proliferation program'
                elif 'migrate' in options:action='migrate';reason='cue gradient and exploration, modulated by current CD8 capacities'
                if behavior is not None and action!='divide':
                    options=options_for(c,cells,neighbors,contacts,targets,sig)
                    action,behavior_decision=behavior.choose(i,options)
                    reason='Recorded behavior choice; kernel revalidates current menu, contact and space'
                if action=='attack':
                    j=min(targets,key=lambda j:(math.hypot(c['x']-cells[j]['x'],c['y']-cells[j]['y']),j));target=j
                    distance=math.hypot(c['x']-cells[j]['x'],c['y']-cells[j]['y'])
                    damage=p['attack_damage_per_hour']*effector*sig['recognition']*(1-sig['checkpoint'])*dt
                    cells[j]['damage']=min(1.,cells[j]['damage']+damage);c['fatigue']=min(1.,c['fatigue']+p['fatigue_per_hour']*dt)
                    event('attack',i,j,damage=damage,distance_um=distance,contact_limit_um=c['radius']+cells[j]['radius']+p['contact_gap_um'],state=c['cd8_state'],checkpoint=sig['checkpoint'],exposure=exposure)
                    if cells[j]['damage']>=1:kill(cells[j],'contact_attack',j)
                elif action=='clear':
                    target=min(corpses);cells[target]['present']=False;event('phagocytosis',i,target)
                elif action in ['migrate','approach']:
                    y,x=np.clip(np.array([c['y'],c['x']])/dx,1,ngrid-2).astype(int)
                    grad=np.array([cue[y,x+1]-cue[y,x-1],cue[y+1,x]-cue[y-1,x]])
                    if np.linalg.norm(grad)>1e-9:grad/=np.linalg.norm(grad)
                    direction=grad+keyed_rng(seed,c['cell_id'],step,'move').normal(0,.65,2);direction/=max(np.linalg.norm(direction),1e-9)
                    if action=='approach':
                        target=min((j for j in neighbors if cells[j]['alive'] and cells[j]['lineage']=='Melanoma'),
                                   key=lambda j:math.hypot(c['x']-cells[j]['x'],c['y']-cells[j]['y']))
                        direction=np.array([cells[target]['x']-c['x'],cells[target]['y']-c['y']]);direction/=max(np.linalg.norm(direction),1e-9)
                    speed=p['immune_speed_um_h']*(1-.4*c['axes']['dysfunction'] if c['lineage']=='CD8_T' else 1)
                    outcome='blocked_by_space'
                    for angle in [0,.8,-.8,1.6,-1.6]:
                        cs,sn=math.cos(angle),math.sin(angle);d=np.array([cs*direction[0]-sn*direction[1],sn*direction[0]+cs*direction[1]])
                        xy=np.array([c['x'],c['y']])+d*speed*dt
                        if free_disc(xy,c['radius'],cells,side,(i,)):c['x'],c['y']=map(float,xy);outcome='completed';break
                    c['fatigue']=max(0.,c['fatigue']-p['recovery_per_hour']*dt*.5)
                elif action=='grow':
                    cycle=min(1.,c['cycle']+(.2+.8*c['programs']['proliferation'])*dt/p['cycle_hours']);radius=c['base_radius']*math.sqrt(1+cycle)
                    if free_disc([c['x'],c['y']],radius,cells,side,(i,)):c['cycle'],c['radius']=cycle,radius
                    else:outcome='blocked_by_space'
                elif action=='divide':outcome=make_child(c,i)
                elif action=='rest':c['fatigue']=max(0.,c['fatigue']-p['recovery_per_hour']*dt)
            c['action']=action if outcome=='completed' else 'blocked';c['target']=target
            y,x=np.clip(np.array([c['y'],c['x']])/dx,0,ngrid-1).astype(int)
            decisions.append({'time_hours':round(t,8),'step':step,'cell':i,'cell_id':c['cell_id'],'options':options,'choice':action,'reason':reason,
                'outcome':outcome,'target':target,'state_decision':fate,'context':{'contacts':len(contacts),'contact_ids':contacts,'target_ids':targets,
                'cue':round(float(cue[y,x]),6),'contact_distance_um':distance,**context,**sig},'seed':seed,
                'rule':'melanoma-runtime-v2','axes':deepcopy(c['axes']),'behavior_decision':behavior_decision})
        if m['demography_enabled'] and stopped=='completed':
            for state,rate in m['influx_cells_per_hour'].items():
                count=int(keyed_rng(seed,step,state,'influx').poisson(rate*dt))
                for serial in range(count):
                    if len(cells)>=max_cells:stopped='population_cap';break
                    if recorded+len(cells)+1>cfg['recording_cell_limit']:stopped='recording_cap';break
                    r=RADII['CD8_T'];rng=keyed_rng(seed,step,state,serial,'entry');position=None
                    for _ in range(200):
                        xy=rng.uniform(r,side-r,2)
                        if mean_mode or free_disc(xy,r,cells,side):position=xy;break
                    if position is None:event('influx_blocked',-1,state=state);continue
                    j=len(cells);new=deepcopy(cells[0]);ident=unique_id('influx',j)
                    new.update(cell_id=ident,index=j,lineage='CD8_T',cell_type='CD8_T',cd8_state=state,time_in_state=0.,
                      axes=deepcopy(m['profiles'][state]),observed_rna_t0={},observed_axes_t0={},programs=dict.fromkeys(['cytotoxicity','proliferation','interferon'],0.),
                      molecular_source_id=None,reference_id='assumed-capacity-v1:'+state,patient_id=None,clone_id=None,parent=None,parent_id=None,
                      born_at=t,x=float(position[0]),y=float(position[1]),radius=r,base_radius=r,alive=True,present=True,damage=0.,fatigue=0.,cycle=0.,
                      death_time=None,action='influx',target=-1,reactive=bool(rng.random()<m['reactive_fraction']),state_source='assumed influx',state_confidence='assumed')
                    cells.append(new);event('influx',j,state=state)
        snapshot()
        if stopped!='completed':break
    root=chain_journal(decisions);counts=Counter(e['event'] for e in events)
    catalog_keys=['cell_id','cell_type','lineage','sample_id','patient_id','clone_id','parent_id','molecular_source_id','state_source',
                  'state_confidence','state_margin','observed_rna_t0','observed_axes_t0','programs','parent','born_at','reactive']
    catalog=[{k:c.get(k) for k in catalog_keys} for c in cells]
    for i,c in enumerate(catalog):c['initial_cd8_state']=payload['cells'][i]['cd8_state'] if i<len(payload['cells']) else next((e.get('state') for e in events if e['target']==i and e['event']=='division' or e['cell']==i and e['event']=='influx'),None)
    return {'schema_version':2,'kind':'vtsc_melanoma_run','runtime':'dynamic-cd8-v2','config':cfg,'effective_model':m,'scenario':arm,
       'behavior_policy':behavior.metadata() if behavior else {'mode':'local','name':'local-priority-v1'},
       'parameter_ledger':parameter_ledger(cfg),'layout':layout,'initial_signature':init_signature,'geometry':{'width_um':side,'height_um':side,'measured':False,'field_grid':ngrid,'mean_context_mode':mean_mode},
       'frame_columns':FRAME_COLUMNS,'frames':frames,'decisions':decisions,'events':events,'catalog':catalog,'provenance':deepcopy(payload['provenance']),
       'summary':{'seed':seed,'initial_cells':len(initial),'final_alive':sum(c['alive'] for c in cells),'created_cells':len(cells),
         'observed_hours':round(t,8),'requested_hours':hours,'stop_reason':stopped,'complete':stopped=='completed','events':dict(counts),'decisions':len(decisions),
         'final_cd8_counts':frames[-1]['cd8_counts'],'policy':policy.name if policy else 'stochastic-hazards-v1','journal_root':root,
         'rates_fitted':m['evidence_class']=='estimated','evidence_class':'exploratory software experiment', 'clinical_validation':False,
         'counterfactual':arm=='exposure-zero','geometry_measured':False,'drug_units':'normalized assumed exposure, not clinical dose'}}
