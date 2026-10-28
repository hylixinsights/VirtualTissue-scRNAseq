"""A05: pinned public acquisition only. No metadata/state analysis."""
import json,subprocess,hashlib,time
from pathlib import Path
root=Path(__file__).resolve().parents[1]
m=json.loads((root/'data/MANIFEST.json').read_text()); results=[]
for rec in m['files']:
 if rec['filename'] not in ['Entire_TME.rds','Malignant_cells.rds']: continue
 target=root.parent/'data'/rec['relative_path']; partial=target.with_suffix('.rds.partial')
 url='https://rdr.kuleuven.be/api/access/datafile/'+('210951' if rec['filename']=='Entire_TME.rds' else '28868')
 start=time.monotonic()
 if not target.exists():
  r=subprocess.run(['curl','--fail','--location','--retry','2','--connect-timeout','30','--max-time','1800','--output',str(partial),url],capture_output=True,text=True)
  (root/'reports/product-v04'/f"{target.stem}-download.log").write_text(r.stderr)
  if r.returncode: results.append({'file':target.name,'status':'download_failed','returncode':r.returncode}); continue
 else: partial=target
 hs={a:hashlib.new(a) for a in ['md5','sha256']}
 with partial.open('rb') as f:
  while b:=f.read(2**20):
   for h in hs.values(): h.update(b)
 good=partial.stat().st_size==rec['expected_bytes'] and hs['md5'].hexdigest()==rec['source_checksum']['value']
 row={'file':target.name,'url':url,'bytes':partial.stat().st_size,**{a:h.hexdigest() for a,h in hs.items()},'status':'verified' if good else 'quarantined_checksum_failure','seconds':time.monotonic()-start,'license':'CC-BY-NC-SA-4.0'};results.append(row)
 if good:
  if partial!=target: partial.rename(target)
  rec['local']={k:v for k,v in row.items() if k in ['bytes','sha256','md5','status']}
  (root/'data/MANIFEST.json').write_text(json.dumps(m,indent=2)+'\n')
 (root/'reports/product-v04/tme-acquisition.json').write_text(json.dumps(results,indent=2)+'\n')
 print(target.name,row['status'],flush=True)
(root/'reports/product-v04/tme-acquisition.json').write_text(json.dumps(results,indent=2)+'\n')
