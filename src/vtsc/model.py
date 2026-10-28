"""Exact continuous-time branching Markov process with declared, unfitted rates."""
from __future__ import annotations

import copy
import math
from pathlib import Path

import numpy as np
import yaml
from scipy.linalg import expm


DEMO_MODEL = {
    'schema_version': 1,
    'name': 'illustrative_cd8_states',
    'rate_provenance': 'User-declared illustrative assumptions; not fitted to expression or patients.',
    'time_unit': 'hours',
    'states': ['Resting-like', 'GZMK-like', 'Cytotoxic-like'],
    'transitions': [
        {'source': 'Resting-like', 'target': 'GZMK-like', 'rate': 0.025},
        {'source': 'GZMK-like', 'target': 'Resting-like', 'rate': 0.01},
        {'source': 'GZMK-like', 'target': 'Cytotoxic-like', 'rate': 0.02},
        {'source': 'Cytotoxic-like', 'target': 'GZMK-like', 'rate': 0.008},
    ],
    'division_rates': {'Resting-like': 0.001, 'GZMK-like': 0.004, 'Cytotoxic-like': 0.002},
    'death_rates': {'Resting-like': 0.001, 'GZMK-like': 0.002, 'Cytotoxic-like': 0.003},
    'influx_rates': {'Resting-like': 0.3},
}


def checked_model(config):
    c = copy.deepcopy(config)
    allowed = {'schema_version', 'name', 'rate_provenance', 'time_unit', 'states',
               'transitions', 'division_rates', 'death_rates', 'influx_rates'}
    if not isinstance(c, dict) or set(c) - allowed:
        raise ValueError('Model must be a mapping with documented keys only.')
    if c.get('schema_version') != 1 or c.get('time_unit') != 'hours':
        raise ValueError('Model schema_version must be 1 and time_unit must be hours.')
    states = c.get('states')
    if not isinstance(states, list) or not 1 <= len(states) <= 20:
        raise ValueError('Supply 1–20 unique states.')
    if any(not isinstance(s, str) or not s.strip() or len(s) > 100 for s in states):
        raise ValueError('State names must be nonempty strings of at most 100 characters.')
    if len(set(states)) != len(states) or 'time_hours' in states:
        raise ValueError('State names must be unique and cannot use the reserved name time_hours.')
    if not c.get('rate_provenance'):
        raise ValueError('Document the source/assumption of rates in rate_provenance.')

    def rate(value):
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise ValueError('Every rate must be a finite nonnegative number.')
        return float(value)

    edges = set()
    transitions = c.setdefault('transitions', [])
    if not isinstance(transitions, list):
        raise ValueError('transitions must be a list.')
    for edge in transitions:
        if not isinstance(edge, dict) or set(edge) != {'source', 'target', 'rate'}:
            raise ValueError('Each transition needs source, target and rate.')
        pair = edge['source'], edge['target']
        if any(s not in states for s in pair) or pair[0] == pair[1] or pair in edges:
            raise ValueError('Transition states must exist, differ and have no duplicate edge.')
        edges.add(pair)
        edge['rate'] = rate(edge['rate'])
    for key in ['division_rates', 'death_rates', 'influx_rates']:
        mapping = c.setdefault(key, {})
        if not isinstance(mapping, dict) or set(mapping) - set(states):
            raise ValueError(f'{key} contains an unknown state.')
        c[key] = {s: rate(mapping.get(s, 0)) for s in states}
    return c


def load_model(path=None, scenario='baseline'):
    c = checked_model(yaml.safe_load(Path(path).read_text()) if path else DEMO_MODEL)
    if scenario == 'no-events':
        for edge in c['transitions']:
            edge['rate'] = 0.0
        for key in ['division_rates', 'death_rates', 'influx_rates']:
            c[key] = dict.fromkeys(c['states'], 0.0)
    elif scenario != 'baseline':
        raise ValueError('Unknown scenario.')
    return c


def mean_field(initial, times, model):
    """Exact expectation of the uncapped linear process, including influx."""
    states = model['states']
    n = len(states)
    a = np.zeros((n + 1, n + 1))
    for edge in model['transitions']:
        i, j = states.index(edge['source']), states.index(edge['target'])
        a[j, i] += edge['rate']
        a[i, i] -= edge['rate']
    for i, state in enumerate(states):
        a[i, i] += model['division_rates'][state] - model['death_rates'][state]
        a[i, n] = model['influx_rates'][state]
    vector = np.r_[initial, 1.0]
    with np.errstate(over='raise', invalid='raise'):
        values = np.array([(expm(a * t) @ vector)[:n] for t in times])
    if not np.isfinite(values).all():
        raise ValueError('Rates/horizon produce an overflowing mean-field expectation; reduce them.')
    return values


def simulate(cells, model, *, hours=48.0, interval=1.0, seed=21, max_cells=5000, max_events=50000):
    """Gillespie events; each live agent has one state and a stable identifier.

    Division is an abstract net +1 birth: parent remains, daughter inherits state.
    Caps stop the run explicitly; they are not density-dependent biology.
    """
    model = checked_model(model)
    if not math.isfinite(hours) or not 0 < hours <= 10000:
        raise ValueError('hours must lie in (0, 10000].')
    if not math.isfinite(interval) or interval <= 0 or hours / interval > 5000:
        raise ValueError('Choose a positive interval producing at most 5000 snapshots.')
    if seed < 0 or max_cells < len(cells) or max_cells > 100000 or not 1 <= max_events <= 1000000:
        raise ValueError('Invalid seed or caps: initial cells <= max_cells <= 100000; 1 <= max_events <= 1000000.')
    states = model['states']
    ids = [str(c['cell_id']) for c in cells]
    if not ids or len(ids) != len(set(ids)) or any(c['state'] not in states for c in cells):
        raise ValueError('Initial cells need unique IDs and declared states.')
    live = {c['cell_id']: c['state'] for c in cells}
    groups = {s: [c['cell_id'] for c in cells if c['state'] == s] for s in states}
    positions = {cell: i for group in groups.values() for i, cell in enumerate(group)}
    all_ids = set(ids)
    rng = np.random.default_rng(seed)
    channels = [('transition', e['source'], e['target'], e['rate']) for e in model['transitions']]
    channels += [(kind, s, s, model[key][s]) for kind, key in
                 [('division', 'division_rates'), ('death', 'death_rates'), ('influx', 'influx_rates')]
                 for s in states]
    times = np.unique(np.r_[np.arange(0, hours, interval), hours])
    initial = [len(groups[s]) for s in states]
    expected = mean_field(initial, times, model)  # fail before starting an unstable run
    rows, events, index, clock, serial = [], [], 0, 0.0, 0
    stop_reason = 'horizon_reached'

    def snapshot(t):
        return {'time_hours': float(t), **{s: len(groups[s]) for s in states}}

    def remove(cell, state):
        group = groups[state]
        pos = positions.pop(cell)
        last = group.pop()
        if pos < len(group):
            group[pos] = last
            positions[last] = pos
        del live[cell]

    def add(cell, state):
        positions[cell] = len(groups[state])
        groups[state].append(cell)
        live[cell] = state

    def new_id():
        nonlocal serial
        while True:
            serial += 1
            value = f'sim_{serial:08d}'
            if value not in all_ids:
                all_ids.add(value)
                return value

    while clock < hours:
        hazards = np.array([rate * (1 if kind == 'influx' else len(groups[source]))
                            for kind, source, _, rate in channels])
        total = hazards.sum()
        if not math.isfinite(total):
            raise ValueError('Total event hazard overflow.')
        if total == 0:
            clock = hours
            break
        next_time = clock + float(rng.exponential(1 / total))
        if next_time <= clock:
            raise ValueError('Rates exceed floating-point time resolution.')
        if next_time > hours:
            clock = hours
            break
        if len(events) >= max_events:
            stop_reason = 'max_events_reached'
            break
        choice = int(rng.choice(len(channels), p=hazards / total))
        kind, source, target, _ = channels[choice]
        if kind in ['division', 'influx'] and len(live) >= max_cells:
            stop_reason = 'max_cells_reached'
            break
        while index < len(times) and times[index] < next_time:
            rows.append(snapshot(times[index])); index += 1
        clock = next_time
        parent = ''
        if kind == 'influx':
            cell = new_id(); add(cell, target)
        else:
            cell = groups[source][int(rng.integers(len(groups[source])))]
            if kind == 'division':
                parent = cell; cell = new_id(); add(cell, target)
            elif kind == 'death':
                remove(cell, source)
            else:
                remove(cell, source); add(cell, target)
        events.append({'event_id': len(events) + 1, 'time_hours': clock, 'event': kind,
                       'cell_id': cell, 'parent_id': parent,
                       'from_state': '' if kind in ['division', 'influx'] else source,
                       'to_state': '' if kind == 'death' else target})
    while index < len(times) and times[index] <= clock:
        rows.append(snapshot(times[index])); index += 1
    if not rows or rows[-1]['time_hours'] != clock:
        rows.append(snapshot(clock))
    observed_times = [r['time_hours'] for r in rows]
    expected = mean_field(initial, observed_times, model)
    expected_rows = [{'time_hours': t, **dict(zip(states, map(float, vector)))}
                     for t, vector in zip(observed_times, expected)]
    return {
        'trajectory': rows, 'mean_field': expected_rows, 'events': events,
        'final_cells': [{'cell_id': key, 'state': value} for key, value in sorted(live.items())],
        'summary': {'schema_version': 1, 'mode': 'uncalibrated_software_experiment', 'seed': seed,
                    'requested_hours': hours, 'observed_hours': clock, 'stop_reason': stop_reason,
                    'complete': stop_reason == 'horizon_reached', 'initial_cells': len(cells),
                    'final_cells': len(live), 'events': len(events), 'max_cells': max_cells,
                    'max_events': max_events, 'rates_fitted': False, 'spatial_coordinates_used': False,
                    'expectation_scope': 'Uncapped linear model; not a confidence interval or fitted prediction.'},
    }
