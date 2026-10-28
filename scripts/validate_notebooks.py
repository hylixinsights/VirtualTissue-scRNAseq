"""Execute the three notebooks through real kernels in a fresh temporary workspace."""
import json
import os
import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import nbformat
from nbclient import NotebookClient

ROOT = Path(__file__).resolve().parents[1]
results = []
with tempfile.TemporaryDirectory(prefix='vtsc-notebooks-') as temporary:
    base = Path(temporary)
    kernel = base / 'jupyter/kernels/vtsc-test'
    kernel.mkdir(parents=True)
    (kernel / 'kernel.json').write_text(json.dumps({
        'argv': [sys.executable, '-m', 'ipykernel_launcher', '-f', '{connection_file}'],
        'display_name': 'VTSC validation', 'language': 'python'}))
    environment = {'JUPYTER_PATH': str(base / 'jupyter'), 'JUPYTER_RUNTIME_DIR': str(base / 'runtime'),
                   'IPYTHONDIR': str(base / 'ipython'), 'MPLCONFIGDIR': str(base / 'matplotlib'),
                   'VTSC_NOTEBOOK_WORK': str(base / 'work')}
    with patch.dict(os.environ, environment):
        for path in sorted((ROOT / 'notebooks').glob('*.ipynb')):
            nb = nbformat.read(path, as_version=4)
            for cell in nb.cells:
                if cell.cell_type == 'code':
                    cell.source = cell.source.replace('INSTALL_PACKAGE = True', 'INSTALL_PACKAGE = False')
                    cell.source = cell.source.replace('DATA_SOURCE = "melanoma"', 'DATA_SOURCE = "synthetic"')
            NotebookClient(nb, timeout=180, kernel_name='vtsc-test',
                           resources={'metadata': {'path': str(ROOT)}}).execute()
            results.append({'notebook': str(path.relative_to(ROOT)),
                            'code_cells': sum(c.cell_type == 'code' for c in nb.cells),
                            'executed_without_errors': True})
            print(path.name, 'PASS', flush=True)
(ROOT / 'reports/xenium-parity/notebook-validation.json').write_text(json.dumps({
    'python': sys.version.split()[0], 'kernel': 'isolated software environment',
    'scope': 'Real local kernels; synthetic melanoma agent path; installation tested separately; hosted Colab not tested',
    'notebooks': results}, indent=2) + '\n')
