"""Run with a clean wheel environment and -I, outside editable source imports."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import tempfile

import vtsc
from vtsc.melanoma_results import verify_dynamic
from vtsc.results import read_run

root = Path(__file__).resolve().parents[1]
package = Path(vtsc.__file__).resolve().parent
assert 'site-packages' in str(package) and '/src/' not in str(package)
for name in ['web/melanoma.html', 'resources/programs_cd8_v1.json',
             'resources/assumed_state_graph_v1.json', 'resources/input_v2.schema.json',
             'resources/run_v2.schema.json', 'resources/config_v2.schema.json',
             'resources/FREEZE.json', 'jev_provider.py', 'jev_payload.py', 'behavior.py']:
    assert (package / name).is_file(), name
work = Path(tempfile.mkdtemp(prefix='vtsc-wheel-'))
cli = str(Path(sys.executable).with_name('vtsc'))

def run(*args):
    subprocess.run([cli, *map(str, args)], cwd=work, check=True, stdout=subprocess.DEVNULL)

run('demo', work / 'demo', '--cells', 30, '--hours', 3, '--decision-mode', 'mock-jev')
run('explore', work / 'demo/anti-pd1.zip', '--html', work / 'reopened.html')
run('check-melanoma', work / 'demo/input.zip', work / 'preflight')
run('recover-melanoma', work / 'demo/anti-pd1.zip', work / 'recovered')
original = read_run(work / 'demo/anti-pd1.zip')
recovered = read_run(work / 'recovered.zip')
assert original['frames'] == recovered['frames']
verify_dynamic(read_run(work / 'demo/comparison.zip'))
assert 'DecompressionStream' in (work / 'reopened.html').read_text()
subprocess.run([sys.executable, '-m', 'pip', 'check'], cwd=work, check=True)
wheel = root / f'dist/virtualtissue_sc-{vtsc.__version__}-py3-none-any.whl'
report = {
    'version': vtsc.__version__, 'python': sys.version.split()[0],
    'installed_module': str(package), 'isolated_import': True,
    'working_directory': str(work), 'all_assets_present': True,
    'mock_paired_demo_completed': True, 'zip_reopened_verified': True,
    'preflight_completed': True, 'recovered_frames_identical': True,
    'actual_jev_http_requests': 0, 'wheel_bytes': wheel.stat().st_size,
    'wheel_sha256': hashlib.sha256(wheel.read_bytes()).hexdigest(),
    'scope': 'Clean venv without system site packages; installed dependencies; imports isolated from checkout',
}
(root / 'reports/xenium-parity/wheel-validation.json').write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(report, indent=2))
