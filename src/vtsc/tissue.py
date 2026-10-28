"""Transparent 2D melanoma agent prototype; all kinetics and spatial priors are declared.

Decisions depend on local contacts, a diffusing generic cue, fatigue and free space.
No LLM, drug effect, antigen recognition or measured histology is inferred from RNA.
"""
from __future__ import annotations
from collections import Counter
from copy import deepcopy
import math

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

from .tissue_inputs import TYPES

DEFAULT = {
    'name': 'melanoma-local-policy-v1', 'dt_hours': .25,
    'packing_fraction': .32, 'immune_speed_um_h': 10.,
    'contact_gap_um': 2., 'perception_um': 55.,
    'cycle_hours': 24., 'attack_damage_per_hour': .65,
    'fatigue_per_hour': .32, 'recovery_per_hour': .22,
    'corpse_hours': 6., 'cue_diffusion_um2_h': 90., 'cue_decay_per_hour': .3,
    'cue_source_per_hour': .5,
    'provenance': 'Declared software hypotheses, not fitted biological rates. Generic cue has arbitrary units.'}
RADII = {'Melanoma': 8., 'CD8_T': 4.8, 'T_cell': 4.8, 'NK': 5.3, 'Myeloid': 7.,
         'Fibroblast': 6.5, 'Endothelial': 6., 'B_cell': 4.8, 'Unknown': 6.}
RECORDING_CELL_LIMIT = 100_000
MOBILE = {'CD8_T', 'T_cell', 'NK', 'Myeloid', 'B_cell'}
FRAME_COLUMNS = ['index', 'x', 'y', 'radius', 'alive', 'action', 'target', 'damage', 'fatigue', 'cycle']


def checked_config(config=None):
    p = deepcopy(DEFAULT)
    if config:
        if set(config) - set(p): raise ValueError(f'Unknown tissue parameters: {set(config) - set(p)}')
        p.update(config)
    for key, value in p.items():
        if key in ['name', 'provenance']: continue
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
            raise ValueError(f'Invalid nonnegative parameter {key}.')
    if not .05 <= p['dt_hours'] <= .5 or not .1 <= p['packing_fraction'] <= .45 or p['cycle_hours'] <= 0:
        raise ValueError('Use dt 0.05–0.5 h, packing fraction 0.1–0.45, and cycle_hours > 0.')
    if p['contact_gap_um'] > 5 or p['perception_um'] < p['contact_gap_um'] or p['corpse_hours'] <= 0:
        raise ValueError('Invalid contact/perception/clearance parameters.')
    limits = {'immune_speed_um_h': 100, 'perception_um': 500, 'cue_diffusion_um2_h': 10000,
              'cue_decay_per_hour': 100, 'cue_source_per_hour': 1000, 'attack_damage_per_hour': 100,
              'fatigue_per_hour': 100, 'recovery_per_hour': 100}
    if any(p[k] > limit for k, limit in limits.items()):
        raise ValueError('Parameter exceeds the numerical limits of this prototype.')
    return p


def initialize(cells, rng, config, layout='tumor-nests'):
    if layout not in ['tumor-nests', 'mixed']: raise ValueError('Unknown synthetic layout.')
    n = len(cells)
    side = max(100., math.sqrt(sum(math.pi * RADII[c['cell_type']] ** 2 for c in cells) / config['packing_fraction']))
    placed = []; out = [None] * n
    # Tumor placement first produces two assumed nests, with no claims about original histology.
    order = sorted(range(n), key=lambda i: cells[i]['cell_type'] != 'Melanoma')
    for i in order:
        c = cells[i]; r = RADII[c['cell_type']]
        for trial in range(20000):
            if layout == 'tumor-nests' and c['cell_type'] == 'Melanoma' and trial < 600:
                center = np.array([.36 if i % 2 else .64, .46 if i % 2 else .54]) * side
                xy = rng.normal(center, side * (.1 + trial / 600 * .25), 2)
            else: xy = rng.uniform(r, side - r, 2)
            if np.any(xy < r) or np.any(xy > side - r): continue
            if placed:
                a = np.asarray(placed)
                if np.any(np.linalg.norm(a[:, :2] - xy, axis=1) < a[:, 2] + r + .2): continue
            break
        else: raise ValueError('Could not place non-overlapping cells; reduce packing_fraction.')
        out[i] = {**deepcopy(c), 'index': i, 'x': float(xy[0]), 'y': float(xy[1]), 'radius': r,
                  'base_radius': r, 'alive': True, 'present': True, 'action': 'observe', 'target': -1,
                  'damage': 0., 'fatigue': 0., 'cycle': 0., 'death_time': None,
                  'parent': None, 'born_at': 0.}
        placed.append([*xy, r])
    return out, side


def free_disc(xy, radius, cells, side, exclude=()):
    if min(xy) < radius or max(xy) > side - radius: return False
    return all(not c['present'] or i in exclude or math.hypot(c['x']-xy[0], c['y']-xy[1]) >= c['radius'] + radius - 1e-7
               for i, c in enumerate(cells))


def decision_menu(c, neighbors, cells, config, recognition):
    """Kernel constructs legal actions; local policy ranks only these actions."""
    options = ['rest']
    contact = [j for j in neighbors if cells[j]['present'] and
               math.hypot(c['x']-cells[j]['x'], c['y']-cells[j]['y']) <= c['radius'] + cells[j]['radius'] + config['contact_gap_um']]
    targets = [j for j in contact if cells[j]['alive'] and cells[j]['cell_type'] == 'Melanoma']
    corpses = [j for j in contact if not cells[j]['alive']]
    if c['cell_type'] in MOBILE and c['fatigue'] < .85: options += ['migrate']
    if recognition and c['cell_type'] in ['CD8_T', 'NK'] and c['fatigue'] < .85 and targets:
        options += ['attack']
    if c['cell_type'] == 'Myeloid' and corpses: options += ['clear']
    if c['cell_type'] == 'Melanoma': options += ['grow' if c['cycle'] < 1 else 'divide']
    return options, contact, targets, corpses


def simulate_tissue(payload, *, hours=24, seed=21, layout='tumor-nests', scenario='recognition-off', config=None, max_cells=1500):
    config = checked_config(config)
    if scenario not in ['recognition-off', 'declared-recognition', 'frozen']: raise ValueError('Unknown scenario.')
    if not math.isfinite(hours) or not 0 < hours <= 96 or not isinstance(seed, int) or seed < 0:
        raise ValueError('Use hours in (0, 96] and a nonnegative integer seed.')
    if not 2 <= len(payload['cells']) <= 2000 or not len(payload['cells']) <= max_cells <= 5000:
        raise ValueError('Use 2–2000 initial cells and a population cap between initial size and 5000.')
    rng = np.random.default_rng(seed)
    cells, side = initialize(payload['cells'], rng, config, layout)
    field_n = 32; cue = np.zeros((field_n, field_n)); dx = side / field_n
    frames, decisions, events = [], [], []
    stopped = 'completed'; t = 0.; recorded_cells = 0
    def event(kind, i, target=-1, **extra):
        events.append({'time_hours': round(t, 6), 'event': kind, 'cell': i, 'target': target, **extra})
    def snapshot():
        nonlocal recorded_cells
        rows = [[i, round(c['x'], 3), round(c['y'], 3), round(c['radius'], 3), int(c['alive']),
                 c['action'], c['target'], round(c['damage'], 4), round(c['fatigue'], 4), round(c['cycle'], 4)]
                for i, c in enumerate(cells) if c['present']]
        recorded_cells += len(rows)
        frames.append({'time_hours': round(t, 6), 'cells': rows,
                       'cue': np.round(cue, 3).ravel().tolist(),
                       'counts': {k: sum(c['alive'] and c['cell_type'] == k for c in cells) for k in TYPES}})
    snapshot()
    while t < hours - 1e-9:
        if recorded_cells + len(cells) > RECORDING_CELL_LIMIT:
            stopped = 'recording_cap'; break
        dt = min(config['dt_hours'], hours - t); t += dt
        if scenario == 'frozen': snapshot(); continue
        # Diffusion/decay with reflecting boundaries; source from melanoma cells.
        for c in cells:
            if c['alive'] and c['cell_type'] == 'Melanoma':
                y, x = np.clip(np.array([c['y'], c['x']]) / dx, 0, field_n-1).astype(int)
                cue[y, x] += config['cue_source_per_hour'] * dt
        cue = gaussian_filter(cue, math.sqrt(2 * config['cue_diffusion_um2_h'] * dt) / dx, mode='reflect')
        cue *= math.exp(-config['cue_decay_per_hour'] * dt)
        for c in cells:
            if c['present'] and not c['alive'] and t - c['death_time'] >= config['corpse_hours']:
                c['present'] = False; event('clearance_timeout', c['index'])
        # Tree finds candidates at step start; recheck distances against current positions.
        xy = np.array([[c['x'], c['y']] for c in cells])
        tree = cKDTree(xy)
        order = rng.permutation(len(cells))
        for i0 in order:
            i = int(i0); c = cells[i]
            if not c['alive']: continue
            neighbors = [j for j in tree.query_ball_point([c['x'], c['y']], config['perception_um'] + config['immune_speed_um_h'] * dt) if j != i and cells[j]['present'] and math.hypot(c['x']-cells[j]['x'], c['y']-cells[j]['y']) <= config['perception_um']]
            options, contacts, targets, corpses = decision_menu(c, neighbors, cells, config, scenario == 'declared-recognition')
            c['target'] = -1
            action, reason = 'rest', 'stationary lineage or fatigue recovery'
            outcome = 'completed'
            if 'clear' in options: action, reason = 'clear', 'corpse within contact range'
            elif 'attack' in options: action, reason = 'attack', 'declared recognition + live melanoma contact + available energy'
            elif 'divide' in options: action, reason = 'divide', 'cycle complete; daughter placement must pass exclusion test'
            elif 'grow' in options: action, reason = 'grow', 'cycle progresses according to initial RNA program and declared cycle time'
            elif 'migrate' in options: action, reason = 'migrate', 'local cue gradient plus exploratory motion'
            distance = None; start_xy = [c['x'], c['y']]
            if action == 'attack':
                j = min(targets, key=lambda j: (math.hypot(c['x']-cells[j]['x'], c['y']-cells[j]['y']), j))
                target = cells[j]; c['target'] = j
                distance = math.hypot(c['x']-target['x'], c['y']-target['y'])
                damage = config['attack_damage_per_hour'] * (.25 + .75*c['programs']['cytotoxicity']) * dt
                target['damage'] = min(1., target['damage'] + damage)
                c['fatigue'] = min(1., c['fatigue'] + config['fatigue_per_hour'] * dt)
                event('attack', i, j, distance_um=distance, contact_limit_um=c['radius']+target['radius']+config['contact_gap_um'], damage=damage)
                if target['damage'] >= 1:
                    target.update(alive=False, death_time=t, action='dead', target=-1)
                    event('death', j, i, cause='contact_attack')
            elif action == 'clear':
                j = min(corpses); c['target'] = j; cells[j]['present'] = False; event('phagocytosis', i, j)
            elif action == 'migrate':
                y, x = np.clip(np.array([c['y'], c['x']]) / dx, 1, field_n-2).astype(int)
                grad = np.array([cue[y, x+1]-cue[y, x-1], cue[y+1, x]-cue[y-1, x]])
                if np.linalg.norm(grad) > 1e-9: grad /= np.linalg.norm(grad)
                direction = grad + rng.normal(0, .65, 2); direction /= max(np.linalg.norm(direction), 1e-9)
                moved = False
                for turn in [0., .8, -.8, 1.6, -1.6]:
                    cs, sn = math.cos(turn), math.sin(turn)
                    d = np.array([cs*direction[0]-sn*direction[1], sn*direction[0]+cs*direction[1]])
                    candidate = np.array(start_xy) + d * config['immune_speed_um_h'] * dt
                    if free_disc(candidate, c['radius'], cells, side, (i,)):
                        c['x'], c['y'] = map(float, candidate); moved = True; break
                if not moved: outcome = 'blocked_by_space'
                c['fatigue'] = max(0., c['fatigue'] - config['recovery_per_hour'] * dt * .5)
            elif action == 'grow':
                cycle = min(1., c['cycle'] + (.2 + .8*c['programs']['proliferation']) * dt / config['cycle_hours'])
                radius = c['base_radius'] * math.sqrt(1 + cycle)
                if free_disc(start_xy, radius, cells, side, (i,)): c['radius'], c['cycle'] = radius, cycle
                else: outcome = 'blocked_by_space'
            elif action == 'divide':
                if recorded_cells + len(cells) + 1 > RECORDING_CELL_LIMIT:
                    outcome = 'recording_cap'; stopped = 'recording_cap'
                elif len(cells) >= max_cells:
                    outcome = 'population_cap'; stopped = 'population_cap'
                else:
                    success = False; r = c['base_radius']
                    for angle in rng.uniform(0, 2*math.pi, 12):
                        offset = r * 1.01 * np.array([math.cos(angle), math.sin(angle)])
                        a, b = np.array(start_xy) - offset, np.array(start_xy) + offset
                        if free_disc(a, r, cells, side, (i,)) and free_disc(b, r, cells, side, (i,)):
                            child = deepcopy(c); j = len(cells)
                            child.update(cell_id=f'born_{j:06d}', index=j, x=float(b[0]), y=float(b[1]), radius=r,
                                         cycle=0., parent=i, born_at=t, action='born', target=-1)
                            c.update(x=float(a[0]), y=float(a[1]), radius=r, cycle=0.)
                            cells.append(child); event('division', i, j); success = True; break
                    if not success: outcome = 'blocked_by_space'
            else: c['fatigue'] = max(0., c['fatigue'] - config['recovery_per_hour'] * dt)
            c['action'] = action if outcome == 'completed' else 'blocked'
            y, x = np.clip(np.array(start_xy[::-1]) / dx, 0, field_n-1).astype(int)
            decisions.append({'time_hours': round(t, 6), 'cell': i, 'options': options, 'choice': action,
                              'reason': reason, 'outcome': outcome, 'target': c['target'],
                              'context': {'contacts': len(contacts), 'cue': round(float(cue[y,x]), 4),
                                          'contact_distance_um': distance}})
        snapshot()
        if stopped != 'completed': break
    counts = Counter(e['event'] for e in events)
    return {'schema_version': 1, 'kind': 'vtsc_tissue_run', 'config': config, 'scenario': scenario,
            'layout': layout, 'geometry': {'width_um': side, 'height_um': side, 'measured': False, 'field_grid': field_n},
            'frame_columns': FRAME_COLUMNS, 'frames': frames, 'decisions': decisions, 'events': events,
            'catalog': [{k: c[k] for k in ['cell_id','cell_type','programs','parent','born_at']} for c in cells],
            'provenance': payload['provenance'],
            'summary': {'seed': seed, 'initial_cells': len(payload['cells']), 'final_alive': sum(c['alive'] for c in cells),
                        'observed_hours': round(t,6), 'requested_hours': hours, 'stop_reason': stopped,
                        'decisions': len(decisions), 'events': dict(counts), 'geometry_measured': False,
                        'policy': 'deterministic local priority with seeded exploratory movement', 'rates_fitted': False}}
