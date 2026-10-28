"""Portable tissue recordings and offline visual inspection."""
import hashlib
import json
import shutil
from pathlib import Path
import yaml

from . import __version__
from .inputs import json_bytes, table_bytes
from .provenance import file_hashes
from .results import read_run
from .tissue import simulate_tissue
from .tissue_inputs import read_tissue_bundle


def write_tissue_html(result, path):
    template = (Path(__file__).parent / 'web/tissue.html').read_text()
    data = json.dumps(result, separators=(',', ':'), allow_nan=False).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')
    Path(path).write_text(template.replace('__PAYLOAD__', data), encoding='utf-8')


def run_tissue(bundle, output, *, config_path=None, **options):
    path = Path(output)
    if path.exists() or Path(str(path) + '.zip').exists(): raise FileExistsError('Use a new tissue run output name.')
    payload = read_tissue_bundle(bundle)
    config = yaml.safe_load(Path(config_path).read_text()) if config_path else None
    result = simulate_tissue(payload, config=config, **options)
    result.update(software_version=__version__, bundle_sha256=file_hashes(bundle)['sha256'])
    path.mkdir(parents=True)
    (path / 'run.json').write_bytes(json_bytes(result))
    (path / 'summary.json').write_bytes(json_bytes(result['summary']))
    (path / 'assumptions.yaml').write_text(yaml.safe_dump(result['config'], sort_keys=False))
    (path / 'provenance.json').write_bytes(json_bytes(result['provenance']))
    (path / 'decisions.jsonl').write_text(''.join(json.dumps(d, separators=(',',':'))+'\n' for d in result['decisions']))
    (path / 'events.jsonl').write_text(''.join(json.dumps(e, separators=(',',':'))+'\n' for e in result['events']))
    rows = [{'time_hours': f['time_hours'], **f['counts']} for f in result['frames']]
    (path / 'trajectory.tsv').write_bytes(table_bytes(rows, list(rows[0])))
    write_tissue_html(result, path / 'explorer.html')
    manifest = {'schema_version': 1, 'kind': 'vtsc_results',
                'files': {p.name: file_hashes(p)['sha256'] for p in sorted(path.iterdir())}}
    (path / 'manifest.json').write_bytes(json_bytes(manifest))
    shutil.make_archive(str(path), 'zip', root_dir=path)
    return {'output': str(path), 'explorer': str(path / 'explorer.html'), **result['summary']}


def tissue_demo(output, *, cells=300, hours=24, seed=21):
    from .tissue_inputs import synthetic_melanoma, prepare_tissue
    path = Path(output)
    if path.exists(): raise FileExistsError(path)
    synthetic_melanoma(path / 'synthetic_melanoma.h5ad', cells, seed)
    bundle = path / 'input.zip'
    prepare_tissue(path / 'synthetic_melanoma.h5ad', bundle, max_cells=cells, seed=seed, source_license='MIT; synthetic fixture')
    runs = {}
    for scenario in ['declared-recognition', 'recognition-off']:
        runs[scenario] = run_tissue(bundle, path / scenario, scenario=scenario, hours=hours, seed=seed)
    return runs
