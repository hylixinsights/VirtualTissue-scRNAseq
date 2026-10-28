"""Versioned shared CD8 hazards. Assumptions, not fitted SPECIAL biology."""
from __future__ import annotations
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
from typing import Protocol

import numpy as np
from scipy.linalg import expm

STATES = ['Memory', 'Activated', 'Cytotoxic', 'Dysfunctional']
AXES = ['dysfunction', 'proliferation', 'cytotoxicity', 'persistence']
SOURCE = 'A05 exploratory assumption; numerical values are not literature estimates'
PROFILES = {
 'Memory': dict(zip(AXES, [.05, .10, .08, .95])),
 'Activated': dict(zip(AXES, [.15, .65, .45, .65])),
 'Cytotoxic': dict(zip(AXES, [.25, .25, .95, .50])),
 'Dysfunctional': dict(zip(AXES, [.90, .70, .30, .45])),
 'Unknown': dict(zip(AXES, [0., 0., 0., 0.])),
}
DEFAULT_MODEL = {
 'schema_version': 2, 'name': 'assumed-cd8-anti-pd1-v1', 'evidence_class': 'assumed',
 'source': SOURCE, 'states': STATES,
 'transitions': [
  {'source': a, 'target': b, 'rate_per_hour': r, 'class': cls, 'evidence_class': 'assumed', 'provenance': SOURCE}
  for a,b,r,cls in [('Memory','Activated',.045,'activation'),('Activated','Memory',.012,'rest'),
   ('Activated','Cytotoxic',.035,'effector'),('Activated','Dysfunctional',.018,'dysfunction'),
   ('Cytotoxic','Dysfunctional',.022,'dysfunction'),('Cytotoxic','Memory',.006,'rest'),
   ('Dysfunctional','Activated',.004,'reinvigoration')]],
 'division_per_hour': dict(zip(STATES,[.001,.012,.003,.015])),
 'death_per_hour': dict(zip(STATES,[.001,.002,.003,.004])),
 'influx_cells_per_hour': dict.fromkeys(STATES,0.),
 'profiles': PROFILES,
 'pd1_activity': dict(zip(STATES,[.2,.6,.55,.9])),
 'pdl1_activity': .7, 'checkpoint_strength': .85, 'reactive_fraction': .65,
 'target_antigen': 1., 'costimulation': .8,
 'schedule': [{'start_hour':2.,'end_hour':96.,'exposure':.8}],
 'checkpoint_enabled': True, 'recognition_enabled': True, 'transitions_enabled': True,
 'demography_enabled': True,
 'parameter_provenance': {
  'transitions.*.rate_per_hour': {'unit':'1/hour','evidence_class':'assumed','source':SOURCE},
  'division_per_hour.*': {'unit':'1/hour','evidence_class':'assumed','source':SOURCE},
  'death_per_hour.*': {'unit':'1/hour','evidence_class':'assumed','source':SOURCE},
  'influx_cells_per_hour.*': {'unit':'cells/hour','evidence_class':'assumed','source':SOURCE},
  'profiles.*.*': {'unit':'normalized capacity [0,1]','evidence_class':'assumed','source':SOURCE},
  'pd1_activity.*,pdl1_activity,checkpoint_strength,reactive_fraction,target_antigen,costimulation':
   {'unit':'fraction [0,1]','evidence_class':'assumed','source':SOURCE},
  'schedule': {'unit':'hours; normalized exposure [0,1]','evidence_class':'assumed','source':'Piecewise constant exposure; not dose, pharmacokinetics or measured occupancy'},
 },
}


def checked_model(model=None):
    m = deepcopy(DEFAULT_MODEL if model is None else model)
    if set(m) != set(DEFAULT_MODEL) or m['schema_version'] != 2 or m['states'] != STATES:
        raise ValueError('Dynamic model requires the complete v2 template and four CD8 states.')
    if m['evidence_class'] not in ['assumed','estimated'] or not m['source'] or not m['parameter_provenance']:
        raise ValueError('Model needs evidence class, source and parameter provenance.')
    def number(x, limit, name):
        if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or not 0<=x<=limit:
            raise ValueError(f'Invalid {name}: finite value in [0, {limit}] required.')
    for k in ['pdl1_activity','checkpoint_strength','reactive_fraction','target_antigen','costimulation']:
        number(m[k],1,k)
    for k in ['checkpoint_enabled','recognition_enabled','transitions_enabled','demography_enabled']:
        if not isinstance(m[k],bool): raise ValueError(f'{k} must be boolean.')
    for k,lim in [('division_per_hour',2),('death_per_hour',2),('influx_cells_per_hour',100),('pd1_activity',1)]:
        if set(m[k])!=set(STATES): raise ValueError(f'{k} must name all four states.')
        for v in m[k].values(): number(v,lim,k)
    if set(m['profiles']) != set(STATES+['Unknown']): raise ValueError('Missing capacity profile.')
    for profile in m['profiles'].values():
        if set(profile)!=set(AXES): raise ValueError('Each profile requires four independent axes.')
        for v in profile.values(): number(v,1,'profile')
    edges=set()
    for e in m['transitions']:
        if set(e)!=set(DEFAULT_MODEL['transitions'][0]): raise ValueError('Invalid transition fields.')
        pair=(e['source'],e['target'])
        if pair[0] not in STATES or pair[1] not in STATES or pair[0]==pair[1] or pair in edges:
            raise ValueError('Invalid or duplicate transition edge.')
        edges.add(pair); number(e['rate_per_hour'],2,'rate_per_hour')
        if e['class'] not in ['activation','effector','rest','dysfunction','reinvigoration'] or e['evidence_class'] not in ['assumed','estimated'] or not e['provenance']:
            raise ValueError('Invalid transition class/provenance.')
    previous=0.
    for w in m['schedule']:
        if set(w)!= {'start_hour','end_hour','exposure'}: raise ValueError('Invalid exposure schedule.')
        for k in ['start_hour','end_hour']: number(w[k],10000,k)
        number(w['exposure'],1,'exposure')
        if w['start_hour']<previous or w['end_hour']<=w['start_hour']: raise ValueError('Schedule must be sorted and non-overlapping.')
        previous=w['end_hour']
    return m


def exposure_at(time, model):
    return next((w['exposure'] for w in model['schedule'] if w['start_hour']<=time<w['end_hour']),0.)


def context_values(context):
    keys = {'contact','target','reactive','tumor_fraction','myeloid_fraction'}
    if set(context)!=keys or any(not math.isfinite(v) or not 0<=v<=1 for v in context.values()):
        raise ValueError('Context requires contact, target, reactive, tumor_fraction, myeloid_fraction in [0,1].')
    return context


def signals(state, context, exposure, model):
    context_values(context)
    checkpoint = (model['checkpoint_strength']*model['pd1_activity'].get(state,0.)*model['pdl1_activity']*(1-exposure)
                  if model['checkpoint_enabled'] else 0.)
    recognition = context['reactive']*context['target']*model['target_antigen']*float(model['recognition_enabled'])
    stimulus = recognition*context['contact']*model['costimulation']*(1-checkpoint)
    return {'exposure':exposure,'checkpoint':checkpoint,'recognition':recognition,'antigen_signal':stimulus,
            'protein_activity_source':'assumed; not inferred from PDCD1/CD274 RNA'}


def hazards(state, context, exposure, model):
    """Single authority for spatial and compartment dynamics; rates per hour."""
    sig=signals(state,context,exposure,model); s=sig['antigen_signal']
    if state not in STATES: return {},sig
    factors={'activation':s,'effector':.1+.9*s,'rest':1-s,
             'dysfunction':.15+s+context['tumor_fraction'], 'reinvigoration':s}
    rates={f'transition:{e["target"]}':e['rate_per_hour']*factors[e['class']]
           for e in model['transitions'] if e['source']==state and model['transitions_enabled']}
    if model['demography_enabled']:
        rates['division']=model['division_per_hour'][state]*(.2+.8*s)
        rates['death']=model['death_per_hour'][state]
    return rates,sig


def probabilities(rates, dt):
    if not math.isfinite(dt) or dt<=0 or any(not math.isfinite(v) or v<0 for v in rates.values()):
        raise ValueError('Finite nonnegative hazards and positive dt required.')
    total=sum(rates.values()); event=-math.expm1(-total*dt)
    return {'stay':1-event, **{k:event*v/total if total else 0. for k,v in rates.items()}}


def keyed_rng(seed, *keys):
    digest=hashlib.sha256(json.dumps([seed,*keys],separators=(',',':')).encode()).digest()
    return np.random.default_rng(int.from_bytes(digest[:16],'little'))


class TransitionProvider(Protocol):
    name: str
    def choose(self, request: dict) -> str: ...


class MockProvider:
    """Offline adapter fixture. No network or paid calls."""
    name='mock-stay-v1'
    def choose(self, request): return 'stay'


class LocalPolicy:
    name='local-largest-transition-probability-v1'
    def choose(self, request):
        return max(request['options'],key=lambda k: request['probabilities'][k])


def draw_decision(state, age, context, exposure, dt, model, seed, cell_id, step, provider=None):
    rates,sig=hazards(state,context,exposure,model); probs=probabilities(rates,dt)
    draw=float(keyed_rng(seed,cell_id,step,'fate').random()); cumulative=0.; selected='stay'
    for option,p in probs.items():
        cumulative+=p
        if draw<cumulative: selected=option; break
    sampled=selected
    legal=['stay']+[k for k,v in probs.items() if k.startswith('transition:') and v>0]
    if provider is not None and selected not in ['death','division']:
        request={'state':state,'time_in_state':age,'context':deepcopy(context),'signals':sig,
                 'options':legal,'probabilities':{k:probs[k] for k in legal}}
        selected=provider.choose(deepcopy(request))
        if selected not in legal: raise ValueError(f'Provider selected illegal transition: {selected}')
    return {'prior_state':state,'new_state':selected.split(':',1)[1] if selected.startswith('transition:') else state,
            'time_in_state_before':age,'hazards_per_hour':rates,'probabilities':probs,'draw':draw,
            'sampled_choice':sampled,'choice':selected,'allowed_transitions':legal,'context':deepcopy(context),
            'signals':sig,'seed':seed,'dt_hours':dt,'provider':provider.name if provider else 'stochastic-hazards-v1',
            'probability_scope':'physical hazards; provider choice is separately recorded, never a physical rate'}


def system_matrix(model, context, exposure, dt=None):
    """Augmented generator, or exact expectation operator of the one-event step."""
    a=np.zeros((5,5)); a[4,4]=1 if dt is not None else 0
    for i,s in enumerate(STATES):
        rates,_=hazards(s,context,exposure,model)
        weights=probabilities(rates,dt) if dt is not None else rates
        a[i,i]=1 if dt is not None else 0
        for k,v in weights.items():
            if k.startswith('transition:'):
                a[i,i]-=v; a[STATES.index(k.split(':')[1]),i]+=v
            elif k=='division': a[i,i]+=v
            elif k=='death': a[i,i]-=v
        a[i,4]=model['influx_cells_per_hour'][s]*(dt if dt is not None else 1) if model['demography_enabled'] else 0
    return a


def solve(initial, times, model=None, context=None, dt=None):
    """Shared continuous expectation; optional fixed-step expectation used by ABM."""
    m=checked_model(model); context=context or {'contact':1.,'target':1.,'reactive':1.,'tumor_fraction':.5,'myeloid_fraction':.1}
    context_values(context); times=np.asarray(times,dtype=float); initial=np.asarray(initial,dtype=float)
    if initial.shape!=(4,) or not np.isfinite(initial).all() or (initial<0).any() or not np.isfinite(times).all() or (times<0).any() or (np.diff(times)<0).any():
        raise ValueError('Four nonnegative initial counts and sorted nonnegative times required.')
    if dt is not None and (not math.isfinite(dt) or dt<=0): raise ValueError('dt must be positive.')
    x=np.r_[initial,1.]; t=0.; out=[]
    edges=sorted({w[k] for w in m['schedule'] for k in ['start_hour','end_hour']})
    for end in times:
        while t<end-1e-10:
            stop=min([end]+[v for v in edges if v>t+1e-10]+([t+dt] if dt else [])); h=stop-t
            a=system_matrix(m,context,exposure_at(t,m),h if dt else None)
            x=(a if dt else expm(a*h))@x; t=stop
            if not np.isfinite(x).all(): raise ValueError('Dynamics overflow; reduce rates or horizon.')
        out.append(x[:4].copy())
    return np.asarray(out)


def chain_journal(records):
    previous='0'*64
    for i,row in enumerate(records):
        row['journal_index']=i; row['previous_hash']=previous
        previous=hashlib.sha256(json.dumps(row,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
        row['record_hash']=previous
    return previous


def verify_journal(records):
    previous='0'*64
    for i,row in enumerate(records):
        copy={k:v for k,v in row.items() if k!='record_hash'}
        if copy.get('previous_hash')!=previous or copy.get('journal_index')!=i: return False
        previous=hashlib.sha256(json.dumps(copy,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
        if row.get('record_hash')!=previous:return False
    return True
