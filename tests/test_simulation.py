import copy
import math

import numpy as np
import pytest

from vtsc.model import checked_model, load_model, simulate, mean_field


def simple_model(**overrides):
    config = {'schema_version': 1, 'time_unit': 'hours', 'name': 'analytic_fixture',
              'rate_provenance': 'Synthetic analytic test', 'states': ['A', 'B'],
              'transitions': [{'source': 'A', 'target': 'B', 'rate': .3}],
              'division_rates': {}, 'death_rates': {}, 'influx_rates': {}}
    config.update(overrides)
    return checked_model(config)


def cells(n=50):
    return [{'cell_id': str(i), 'state': 'A', 'sample_id': 'test'} for i in range(n)]


def test_no_events_control_is_identical_at_every_time():
    model = simple_model(transitions=[])
    result = simulate(cells(), model, hours=5, interval=2)
    assert [r['time_hours'] for r in result['trajectory']] == [0, 2, 4, 5]
    assert all(r['A'] == 50 and r['B'] == 0 for r in result['trajectory'])
    assert result['events'] == [] and result['summary']['complete']
    assert result['final_cells'] == [{'cell_id': str(i), 'state': 'A'} for i in sorted(range(50), key=str)]


def test_seed_conservation_event_replay_and_input_immutability():
    model = simple_model(division_rates={'A': .1}, death_rates={'B': .1}, influx_rates={'B': 2})
    original, original_model = cells(), copy.deepcopy(model)
    a = simulate(original, model, hours=5, seed=8)
    b = simulate(original, model, hours=5, seed=8)
    assert a == b and original == cells() and model == original_model
    live = {c['cell_id']: c['state'] for c in original}
    seen = set(live)
    times = [e['time_hours'] for e in a['events']]
    assert times == sorted(times)
    for event in a['events']:
        key, kind = event['cell_id'], event['event']
        if kind in ['division', 'influx']:
            assert key not in seen
            if kind == 'division': assert event['parent_id'] in live
            seen.add(key); live[key] = event['to_state']
        else:
            assert live[key] == event['from_state']
            if kind == 'death': del live[key]
            else: live[key] = event['to_state']
    assert live == {c['cell_id']: c['state'] for c in a['final_cells']}
    births = sum(e['event'] in ['division', 'influx'] for e in a['events'])
    deaths = sum(e['event'] == 'death' for e in a['events'])
    assert len(live) == len(original) + births - deaths


def test_analytic_transition_and_monte_carlo_agree():
    model = simple_model()
    expectation = mean_field([100, 0], [2], model)[0]
    p = math.exp(-.3 * 2)
    assert np.allclose(expectation, [100 * p, 100 * (1 - p)])
    observed = [simulate(cells(100), model, hours=2, interval=2, seed=i)['trajectory'][-1]['A'] for i in range(180)]
    standard_error = math.sqrt(100 * p * (1-p) / len(observed))
    assert abs(np.mean(observed) - 100 * p) < 5 * standard_error


def test_linear_birth_death_influx_expectation():
    model = simple_model(transitions=[], division_rates={'A': .2}, death_rates={'A': .1}, influx_rates={'A': 3})
    expected = 10 * math.exp(.1 * 4) + 3 * math.expm1(.1 * 4) / .1
    assert mean_field([10, 0], [4], model)[0, 0] == pytest.approx(expected)


def test_caps_are_reported_as_incomplete_and_do_not_invent_followup():
    result = simulate(cells(1), simple_model(transitions=[], division_rates={'A': 1}), hours=100, max_cells=1)
    assert result['summary']['stop_reason'] == 'max_cells_reached'
    assert not result['summary']['complete'] and result['summary']['observed_hours'] == 0
    assert result['events'] == []
    result = simulate(cells(100), simple_model(), hours=100, max_events=2)
    assert result['summary']['stop_reason'] == 'max_events_reached'
    assert len(result['events']) == 2
    assert result['trajectory'][-1]['time_hours'] == result['events'][-1]['time_hours'] < 100


def test_extinction_remains_zero_without_influx():
    result = simulate(cells(1), simple_model(transitions=[], death_rates={'A': 10}), hours=10)
    assert result['final_cells'] == [] and result['trajectory'][-1]['A'] == 0
    assert len(result['events']) == 1 and result['events'][0]['event'] == 'death'


@pytest.mark.parametrize('rate', [-1, float('nan'), float('inf'), True, '0.1'])
def test_bad_rates_rejected(rate):
    with pytest.raises(ValueError): simple_model(death_rates={'A': rate})


def test_bad_states_duplicate_edges_and_initial_cap_rejected():
    with pytest.raises(ValueError): simple_model(death_rates={'Typo': .1})
    edge = {'source': 'A', 'target': 'B', 'rate': .1}
    with pytest.raises(ValueError): simple_model(transitions=[edge, edge])
    with pytest.raises(ValueError): simulate(cells(), simple_model(), max_cells=1)
    with pytest.raises(ValueError): simulate(cells(), simple_model(), interval=0)
    with pytest.raises(ValueError): simple_model(states=['time_hours', 'B'])
    with pytest.raises(ValueError): simple_model(states=[{}, 'B'])
    with pytest.raises(ValueError): simple_model(transitions={'A': 'B'})
