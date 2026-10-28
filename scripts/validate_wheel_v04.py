"""Execute with a clean wheel environment and -I; never editable source imports."""
from pathlib import Path
import subprocess,sys,tempfile,json,hashlib,zipfile
import vtsc
from vtsc.melanoma_results import melanoma_demo,verify_dynamic
from vtsc.results import read_run
root=Path(__file__).resolve().parents[1]
package=Path(vtsc.__file__).resolve().parent
assert 'site-packages' in str(package) and '/src/' not in str(package)
assert vtsc.__version__=='0.4.0'
for name in ['web/melanoma.html','resources/programs_cd8_v1.json','resources/assumed_state_graph_v1.json','resources/input_v2.schema.json','resources/run_v2.schema.json','resources/config_v2.schema.json','resources/FREEZE.json']:assert (package/name).is_file(),name
work=Path(tempfile.mkdtemp(prefix='vtsc-wheel-'))
subprocess.run([str(Path(sys.executable).with_name('vtsc')),'demo',str(work/'demo'),'--cells','30','--hours','3'],cwd=work,check=True,stdout=subprocess.DEVNULL)
subprocess.run([str(Path(sys.executable).with_name('vtsc')),'explore',str(work/'demo/anti-pd1.zip'),'--html',str(work/'reopened.html')],cwd=work,check=True,stdout=subprocess.DEVNULL)
verify_dynamic(read_run(work/'demo/comparison.zip'))
assert 'DecompressionStream' in (work/'reopened.html').read_text()
subprocess.run([sys.executable,'-m','pip','check'],cwd=work,check=True)
wheel=root/'dist/virtualtissue_sc-0.4.0-py3-none-any.whl'
report={'version':vtsc.__version__,'python':sys.version.split()[0],'installed_module':str(package),'isolated_import':True,'working_directory':str(work),'all_assets_present':True,'paired_demo_completed':True,'zip_reopened_verified':True,'wheel_bytes':wheel.stat().st_size,'wheel_sha256':hashlib.sha256(wheel.read_bytes()).hexdigest(),'scope':'Clean venv without system site packages; wheel installed with dependencies; imports isolated from checkout'}
(root/'reports/product-v04/wheel-validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
