"""Offline parity integration. Requires the accepted local v0.4 real input/run."""
from pathlib import Path
import json
import time
from vtsc import __version__
from vtsc.melanoma_results import (melanoma_demo, preflight_melanoma, run_melanoma,
                                 compare_melanoma, recover_melanoma, implementation_hashes)
from vtsc.melanoma import simulate_melanoma
from vtsc.results import read_run
from vtsc.tissue_inputs import read_tissue_bundle

root = Path(__file__).resolve().parents[1]
out = root / 'outputs/v05-final'
if out.exists():
    raise FileExistsError('Preserve accepted outputs; choose a new output name in this script.')
start = time.monotonic()
bundle = root / 'outputs/v04-final/real-input.zip'
report = {'actual_jev_http_requests': 0, 'software_version': __version__,
          'implementation_hashes': implementation_hashes()}
report['synthetic'] = melanoma_demo(out / 'synthetic', cells=120, hours=12, decision_mode='mock-jev')
report['preflight'] = preflight_melanoma(bundle, out / 'real-preflight')
report['real'] = {arm: run_melanoma(bundle, out / arm, hours=12, arm=arm, decision_mode='mock-jev')
                  for arm in ['anti-pd1', 'exposure-zero', 'recognition-off']}
report['comparison'] = compare_melanoma(out / 'anti-pd1', out / 'exposure-zero', out / 'comparison')
report['recovery'] = recover_melanoma(out / 'anti-pd1.zip', out / 'recovered')
old = read_run(root / 'outputs/v04-final/real-anti-pd1')
current = simulate_melanoma(read_tissue_bundle(bundle), hours=24, seed=21, arm='anti-pd1')
assert current['frames'] == old['frames']
report['v04_local_frames_preserved_exactly'] = True
assert all(r['complete'] for r in report['real'].values())
assert report['real']['recognition-off']['events'].get('attack', 0) == 0
report['seconds'] = time.monotonic() - start
(root / 'reports/xenium-parity/integration.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps({'completed': True, 'seconds': report['seconds'], 'actual_jev_http_requests': 0}))
