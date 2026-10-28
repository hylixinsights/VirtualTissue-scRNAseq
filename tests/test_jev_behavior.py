from copy import deepcopy
from pathlib import Path
import json
import numpy as np
import pytest
from vtsc.jev_provider import JevProvider,ProviderError,BudgetExhausted,validate_response,NoRedirect
from vtsc.jev_payload import PARAMETERS,build_batches,payload,mock_transport,canonical
from vtsc.behavior import BehaviorController
from vtsc.melanoma import simulate_melanoma,DEFAULT_CONFIG
from vtsc.melanoma_results import run_melanoma,recover_melanoma,preflight_melanoma,verify_dynamic,compare_melanoma
from vtsc.tissue_inputs import synthetic_melanoma,prepare_tissue,read_tissue_bundle
from vtsc.results import read_run


def row(i=1,label='CD8_T'):
    return {'cell_id':i,'options':['rest','migrate'],'observation':{'lineage':label,'time_hours':0}}


def provider(tmp_path,transport=mock_transport,limit=500,context=None,name='journal.jsonl'):
    return JevProvider('FAKE_SECRET_FOR_TEST','jev-1.14.0',PARAMETERS,limit,tmp_path/name,context or {'run':'test'},transport=transport,pause_s=0)


def fixture(tmp_path,n=40):
    raw=tmp_path/'source.h5ad';synthetic_melanoma(raw,n);bundle=tmp_path/'input.zip';prepare_tissue(raw,bundle,source_license='MIT')
    return bundle


def test_batches_pin_model_isolate_effectors_and_enforce_bytes():
    rows=[row(i,'CD8_T' if i%3==0 else 'NK' if i%3==1 else 'Melanoma') for i in range(90)]
    batches=build_batches(rows,'jev-1.14.0',PARAMETERS)
    assert sum(map(len,batches))==90
    for batch in batches:
        assert len(batch)<=20 and len(canonical(payload(batch,'jev-1.14.0')))<=40000
        assert len({r['observation']['lineage'] for r in batch})==1
    with pytest.raises(ValueError):build_batches([{'cell_id':0,'options':['rest'],'observation':{'lineage':'CD8_T','extra':'x'*40001}}],'jev-1.14.0',PARAMETERS)


def test_initial_non_effector_payloads_do_not_receive_arm_flags(tmp_path):
    bundle=fixture(tmp_path)
    preflight_melanoma(bundle,tmp_path/'a',arm='anti-pd1')
    preflight_melanoma(bundle,tmp_path/'b',arm='recognition-off')
    def ordinary(path):
        d=json.loads(path.read_text())
        return [body for body in d['first_round_payloads'] if all(q['instructions']['local_observation']['lineage'] not in ['CD8_T','NK'] for q in body['questions'].values())]
    assert ordinary(tmp_path/'a/preflight.json')==ordinary(tmp_path/'b/preflight.json')


def test_comparison_rejects_mixed_behavior_policies(tmp_path):
    bundle=fixture(tmp_path)
    run_melanoma(bundle,tmp_path/'local',hours=1)
    run_melanoma(bundle,tmp_path/'mock',hours=1,decision_mode='mock-jev')
    with pytest.raises(ValueError,match='behavior policy'):
        compare_melanoma(tmp_path/'local',tmp_path/'mock',tmp_path/'bad')


def test_response_rounding_and_invalid_choices_are_audited():
    body=payload([row()],'jev-1.14.0');r=mock_transport(body)
    r['answers']['c1'].update(choice='rest',probabilities={'rest':.49,'migrate':.50})
    answer,_=validate_response(r,body)
    assert answer['c1']['reported_probabilities']['rest']==.49
    assert sum(answer['c1']['probabilities'].values())==pytest.approx(1)
    for field,value in [('choice','divide'),('probabilities',{'rest':0.,'migrate':float('nan')}),('probabilities',{'rest':.1,'migrate':.9})]:
        bad=deepcopy(r);bad['answers']['c1'][field]=value
        with pytest.raises(ProviderError):validate_response(bad,body)
    r['model']='jev-latest'
    with pytest.raises(ProviderError):validate_response(r,body)
    assert NoRedirect().redirect_request(None,None,302,'',{},'https://other.invalid') is None


def test_journal_secret_filter_replay_and_binding(tmp_path):
    calls=[]
    def transport(body):
        calls.append(body);r=mock_transport(body);r['secret']='FAKE_SECRET_FOR_TEST';r['usage']['extra']='FAKE_SECRET_FOR_TEST';return r
    api=provider(tmp_path,transport);answer=api.choose([row()]);assert len(calls)==1
    assert 'FAKE_SECRET_FOR_TEST' not in api.journal.read_text()
    replay=provider(tmp_path,lambda b:pytest.fail('Replay contacted provider'))
    assert replay.choose([row()])==answer and replay.usage()['reused_batches']==1
    with pytest.raises(ProviderError):provider(tmp_path,context={'run':'changed'})
    for model,parameters in [('jev-1.15.0',PARAMETERS),('jev-1.14.0',{**PARAMETERS,'max_questions':10})]:
        with pytest.raises(ProviderError):
            JevProvider('',model,parameters,500,api.journal,{'run':'test'},transport=mock_transport)
    content=api.journal.read_text().replace('"sequence":0','"sequence":99');api.journal.write_text(content)
    with pytest.raises(ProviderError):provider(tmp_path)


def test_failed_attempt_is_never_resent_and_cap_precedes_round(tmp_path):
    calls=[]
    def fail(body):calls.append(body);raise RuntimeError('FAKE_SECRET_FOR_TEST')
    api=provider(tmp_path,fail)
    with pytest.raises(ProviderError):api.choose([row()])
    assert len(calls)==1 and 'FAKE_SECRET_FOR_TEST' not in api.journal.read_text()
    retry=provider(tmp_path,fail)
    with pytest.raises(ProviderError):retry.choose([row()])
    assert len(calls)==1
    cap=provider(tmp_path,limit=1,name='cap.jsonl')
    with pytest.raises(BudgetExhausted):cap.choose([row(i) for i in range(21)])
    assert cap.requests==0 and not cap.journal.exists()


def test_mock_run_recovery_preflight_and_live_opt_in(tmp_path,monkeypatch):
    bundle=fixture(tmp_path)
    with pytest.raises(ValueError,match='allow-paid'):run_melanoma(bundle,tmp_path/'live',decision_mode='jev')
    preview=preflight_melanoma(bundle,tmp_path/'preview')
    assert preview['actual_http_requests']==0 and preview['first_round_questions']>0
    document=(tmp_path/'preview/preflight.json').read_text()
    for forbidden in ['observed_rna_t0','cell_id','sample_id','patient_id','molecular_source_id','rate_per_hour']:assert forbidden not in document
    run_melanoma(bundle,tmp_path/'run',hours=3,decision_mode='mock-jev')
    r=verify_dynamic(read_run(tmp_path/'run.zip'));assert r['behavior_policy']['usage']['http_attempts']>0
    assert any(d['behavior_decision'] and d['behavior_decision'].get('probabilities') for d in r['decisions'])
    assert all(d['state_decision']['provider']=='stochastic-hazards-v1' for d in r['decisions'] if d['state_decision'])
    import vtsc.melanoma_results as results
    monkeypatch.setattr(results,'mock_transport',lambda b:pytest.fail('Recovery sent a request'))
    recover_melanoma(tmp_path/'run.zip',tmp_path/'recovered')
    recovered=read_run(tmp_path/'recovered');assert recovered['frames']==r['frames'] and recovered['decisions']==r['decisions']
    recover_melanoma(tmp_path/'run.zip',tmp_path/'extended',hours=5)
    assert read_run(tmp_path/'extended')['summary']['stop_reason']=='provider_budget_or_replay_exhausted'


def test_provider_budget_returns_complete_prefix_and_failed_round_is_atomic(tmp_path):
    bundle=fixture(tmp_path,80);p=read_tissue_bundle(bundle)
    api=provider(tmp_path,limit=1);controller=BehaviorController(api,'mock-jev')
    r=simulate_melanoma(p,hours=2,behavior=controller)
    assert not r['summary']['complete'] and r['summary']['observed_hours']==0
    assert not r['decisions'] and len(r['frames'])==1 and api.requests==0
    def fail_after_first(body):
        fail_after_first.n+=1
        if fail_after_first.n==2:raise ProviderError('Second batch rejected')
        return mock_transport(body)
    fail_after_first.n=0
    api=provider(tmp_path,fail_after_first,name='failed.jsonl');controller=BehaviorController(api,'mock-jev')
    r=simulate_melanoma(p,hours=2,behavior=controller)
    assert r['summary']['stop_reason']=='provider_failed' and r['summary']['observed_hours']==0
    assert not r['decisions'] and not controller.answers and api.requests==2


def test_mock_approach_moves_toward_target_with_exclusion_and_kernel_hazards(tmp_path,monkeypatch):
    import vtsc.melanoma as engine
    from test_dynamic_melanoma import payload as make_payload
    init=engine.initialize
    def controlled(*args):
        cells,side=init(*args);cells[0].update(x=65.,y=50.);cells[1].update(x=30.,y=50.);return cells,max(100.,side)
    monkeypatch.setattr(engine,'initialize',controlled)
    cfg=deepcopy(DEFAULT_CONFIG);cfg['model']['reactive_fraction']=1.;cfg['model']['demography_enabled']=False;cfg['model']['transitions_enabled']=False
    controller=BehaviorController(provider(tmp_path),'mock-jev')
    r=simulate_melanoma(make_payload(2),hours=1,config=cfg,behavior=controller)
    positions=[next(row for row in f['cells'] if row[0]==1) for f in r['frames']]
    assert positions[-1][1]>positions[0][1]
    assert any(d['choice']=='approach' for d in r['decisions'])
    for a,b in zip(positions,positions[1:]):assert np.hypot(b[1]-a[1],b[2]-a[2])<=cfg['mechanics']['immune_speed_um_h']*.25+1e-4
    for f in r['frames']:
        a,b=f['cells'];assert np.hypot(a[1]-b[1],a[2]-b[2])+1e-4>=a[3]+b[3]


def test_dead_or_changed_target_cancels_cached_action_and_daily_ceiling(tmp_path):
    api=provider(tmp_path);controller=BehaviorController(api,'mock-jev')
    controller.answers[0]={'choice':'attack','probabilities':{'rest':0,'attack':1},'confidence':1}
    choice,record=controller.choose(0,['rest','migrate'])
    assert choice=='rest' and record['status']=='cancelled_by_kernel'
    bundle=fixture(tmp_path,30);p=read_tissue_bundle(bundle)
    r=simulate_melanoma(p,hours=26,behavior=controller)
    counts={}
    for line in api.journal.read_text().splitlines():
        entry=json.loads(line)
        if entry['phase']!='intent':continue
        for qid,q in entry['request']['questions'].items():
            day=int(q['instructions']['local_observation']['time_hours']//24);key=(qid,day);counts[key]=counts.get(key,0)+1
    assert max(counts.values())<=24
