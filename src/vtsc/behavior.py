"""Pre-step, atomic behavior decisions. The kernel retains ownership of cell fate."""
from copy import deepcopy
import math
from scipy.spatial import cKDTree
from .dynamics import signals
from .tissue import MOBILE
from .jev_payload import PROMPT_VERSION, payload, build_batches

SCHEDULE={'free_interval_hours':1.,'contact_interval_hours':.5,'max_per_cell_per_day':24,
          'semantics':'Ask at first eligible tick, then at periodic program checkpoints; contacted effectors use the shorter interval. No call for a single legal action.'}


def local_context(c,cells,p,indices,model,exposure):
    neighbors=[j for j in indices if j!=c['index'] and cells[j]['present'] and math.hypot(c['x']-cells[j]['x'],c['y']-cells[j]['y'])<=p['perception_um']]
    contacts=[j for j in neighbors if math.hypot(c['x']-cells[j]['x'],c['y']-cells[j]['y'])<=c['radius']+cells[j]['radius']+p['contact_gap_um']]
    targets=[j for j in contacts if cells[j]['alive'] and cells[j]['lineage']=='Melanoma']
    live=[j for j in neighbors if cells[j]['alive']]
    ctx={'contact':float(bool(targets)),'target':float(bool(targets)),'reactive':float(c['reactive']),
         'tumor_fraction':sum(cells[j]['lineage']=='Melanoma' for j in live)/max(1,len(live)),
         'myeloid_fraction':sum(cells[j]['lineage']=='Myeloid' for j in live)/max(1,len(live))}
    return neighbors,contacts,targets,ctx,signals(c['cd8_state'],ctx,exposure,model)


def options_for(c,cells,neighbors,contacts,targets,sig):
    options=['rest']
    if c['lineage'] in MOBILE and c['fatigue']<.85:options.append('migrate')
    effector=c['axes']['cytotoxicity'] if c['lineage']=='CD8_T' else c['programs']['cytotoxicity']
    if c['lineage'] in ['CD8_T','NK'] and targets and sig['recognition']>0 and effector>0 and c['fatigue']<.85:options.append('attack')
    if c['lineage'] in ['CD8_T','NK'] and not targets and c['reactive'] and c['fatigue']<.85 and any(cells[j]['alive'] and cells[j]['lineage']=='Melanoma' for j in neighbors):options.append('approach')
    if c['lineage']=='Myeloid' and any(not cells[j]['alive'] for j in contacts):options.append('clear')
    if c['lineage']=='Melanoma' and c['cycle']<1:options.append('grow')
    return options


def request_for(c,cells,p,indices,model,exposure,time):
    neighbors,contacts,targets,ctx,sig=local_context(c,cells,p,indices,model,exposure)
    options=options_for(c,cells,neighbors,contacts,targets,sig)
    if len(options)==1 or not c['alive'] or not c['present']:return None
    obs={'time_hours':round(time,8),'lineage':c['lineage'],
         't0_program_scores':{k:c['programs'][k] for k in ['cytotoxicity','proliferation','interferon']},
         'simulated':{'fatigue':c['fatigue'],'damage':c['damage'],'cycle':c['cycle']},
         'contacts':[{'index':j,'lineage':cells[j]['lineage'],'alive':cells[j]['alive']} for j in sorted(contacts)],
         'nearby_tumor_count':sum(cells[j]['alive'] and cells[j]['lineage']=='Melanoma' for j in neighbors)}
    if c['lineage']=='CD8_T':
        obs.update(cd8_state=c['cd8_state'],state_age_hours=c['time_in_state'],assumed_capacities=deepcopy(c['axes']),
                   inferred_axes_t0={k:c.get('observed_axes_t0',{}).get(k) for k in ['dysfunction','proliferation','cytotoxicity','persistence']},
                   simulated_exposure=exposure,simulated_checkpoint=sig['checkpoint'])
    return {'cell_id':c['index'],'options':options,'observation':obs}


class BehaviorController:
    def __init__(self,provider,mode):
        self.provider=provider;self.mode=mode;self.answers={};self.asked={};self.per_day={};self.rounds=[]
    def prepare(self,cells,p,model,exposure,time):
        tree=cKDTree([[c['x'],c['y']] for c in cells]);pending=[]
        for c in cells:
            row=request_for(c,cells,p,tree.query_ball_point([c['x'],c['y']],p['perception_um']),model,exposure,time)
            if row is None:continue
            i=c['index'];interval=.5 if row['observation']['contacts'] and c['lineage'] in ['CD8_T','NK'] else 1.
            if time-self.asked.get(i,-1e9)<interval-1e-9 or self.per_day.get((i,int(time//24)),0)>=24:continue
            pending.append(row)
        # All answers validate before this round can mutate the simulation or scheduler.
        responses=self.provider.choose(pending)
        for row in pending:
            i=row['cell_id'];answer=responses[i]
            if answer['choice'] not in row['options']:raise ValueError('Provider returned an illegal action.')
            self.answers[i]={**deepcopy(answer),'decision_time_hours':round(time,8),'offered_options':row['options'],
                             'provider':self.provider.name,'model':self.provider.model,'prompt_version':PROMPT_VERSION}
            self.asked[i]=time;day=(i,int(time//24));self.per_day[day]=self.per_day.get(day,0)+1
        self.rounds.append({'time_hours':round(time,8),'questions':len(pending),'http_batches':len(build_batches(pending,self.provider.model,self.provider.p))})
    def choose(self,index,legal):
        answer=deepcopy(self.answers.get(index));choice='rest' if answer is None else answer['choice']
        if answer is None:return 'rest',{'provider':self.provider.name,'status':'no_current_choice','current_legal_options':legal}
        answer['current_legal_options']=list(legal)
        if choice not in legal:
            answer.update(status='cancelled_by_kernel',executed_choice='rest')
            return 'rest',answer
        answer.update(status='applied',executed_choice=choice)
        return choice,answer
    def metadata(self):
        return {'mode':self.mode,'name':self.provider.name,'model':self.provider.model,'prompt_version':PROMPT_VERSION,
                'schedule':SCHEDULE,'usage':self.provider.usage(),'rounds':self.rounds,
                'journal_root':self.provider.chain,'request_context':self.provider.context,
                'live_validation':False,'probability_scope':'Behavior ranking only; not state-transition, division or death hazards'}
