"""Minimal typed observations for Jev behavior choice; never rates or raw RNA."""
import hashlib
import json

PROMPT_VERSION='vtsc-behavior-1'
PARAMETERS={'max_questions':20,'max_payload_bytes':40000,'max_response_bytes':1000000}
DESCRIPTIONS={
 'rest':'Rest and recover simulated fatigue until the next decision checkpoint.',
 'migrate':'Move along the local generic cue with bounded exploratory motion.',
 'approach':'Move toward a nearby live tumor target, subject to space and speed limits.',
 'attack':'Deliver contact-mediated cytotoxic damage to a permitted live target.',
 'clear':'Remove a contacted dead cell through the simulated clearance mechanism.',
 'grow':'Advance the assumed tumor growth clock; the kernel checks space and division.'}


def canonical(value):return json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
def fingerprint(value):return hashlib.sha256(canonical(value)).hexdigest()


def payload(rows,model):
    return {'model':model,'state':{'setting':'Independent local cells in an exploratory melanoma tissue reconstructed from dissociated scRNA-seq; positions are assumed.',
      'prompt_version':PROMPT_VERSION,'semantics':'Initial rank scores are inferred RNA summaries, not protein activity. Simulated state, capacities, contacts and exposure are assumptions. Choose a permitted behavior; probabilities rank choices, not kinetic rates. No clinical response prediction.'},
      'questions':{f'c{r["cell_id"]}':{'type':'choice','instructions':{'task':'Choose one permitted local behavior for this cell.',
        'local_observation':r['observation']},'criteria':{k:DESCRIPTIONS[k] for k in r['options']}} for r in rows}}


def build_batches(requests,model,p):
    batches=[]
    # Recognition-dependent effectors never share context with other lineages.
    for group in ['CD8_T','NK','other']:
        rows=sorted((r for r in requests if (r['observation']['lineage'] if r['observation']['lineage'] in ['CD8_T','NK'] else 'other')==group),key=lambda r:r['cell_id'])
        batch=[]
        for row in rows:
            candidate=batch+[row]
            if len(candidate)>p['max_questions'] or len(canonical(payload(candidate,model)))>p['max_payload_bytes']:
                if not batch:raise ValueError('A single choice exceeds the payload limit.')
                batches.append(batch);batch=[row]
                if len(canonical(payload(batch,model)))>p['max_payload_bytes']:raise ValueError('A single choice exceeds the payload limit.')
            else:batch=candidate
        if batch:batches.append(batch)
    return batches


def mock_transport(body):
    """Programmed offline answers exercising the real parser/journal, not Jev evidence."""
    answers={}
    for qid,q in body['questions'].items():
        choice=next(k for k in ['clear','attack','approach','grow','migrate','rest'] if k in q['criteria'])
        answers[qid]={'type':'choice','choice':choice,'probabilities':{k:float(k==choice) for k in q['criteria']},'confidence':1.}
    return {'model':body['model'],'answers':answers,'usage':{'input_tokens':0,'output_tokens':0}}
