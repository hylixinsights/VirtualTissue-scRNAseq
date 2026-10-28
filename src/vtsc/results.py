"""Portable run artifacts, static scientific plots and an offline event explorer."""
import hashlib
import html
import json
import platform
import shutil
import sys
import zipfile
from pathlib import Path

import numpy as np

from . import __version__
from .inputs import json_bytes, table_bytes, read_bundle
from .model import load_model, simulate
from .provenance import file_hashes

COLORS = ['#3870ca', '#dd8b25', '#169880', '#9b60b0', '#d04a59', '#647480', '#b07d50', '#2d9cb6']


def plot_run(result, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator
    states = result['model']['states']
    rows, expected = result['trajectory'], result['mean_field']
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), layout='constrained')
    times = [r['time_hours'] for r in rows]
    totals = np.array([sum(r[s] for s in states) for r in rows])
    for i, state in enumerate(states):
        color = COLORS[i % len(COLORS)]
        axes[0].step(times, [r[state] for r in rows], where='post', color=color, label=state)
        axes[0].plot(times, [r[state] for r in expected], '--', color=color, alpha=.65)
    proportions = np.array([[r[s] for r in rows] for s in states], dtype=float)
    np.divide(proportions, totals, out=proportions, where=totals != 0)
    axes[1].stackplot(times, proportions, colors=[COLORS[i % len(COLORS)] for i in range(len(states))])
    axes[0].set(ylabel='Simulated cells', title='Sampled trajectory and uncapped expectation')
    axes[1].set(ylabel='State fraction (zero if extinct)', ylim=(0, 1), title='Composition of the simulated population')
    for ax in axes:
        ax.set_xlabel('Simulated time (hours)'); ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', alpha=.2); ax.xaxis.set_major_locator(MaxNLocator(5))
    axes[0].legend(fontsize=8, loc='best')
    fig.suptitle('Software experiment · declared rates · no clinical calibration', fontsize=13)
    fig.savefig(path, dpi=150); plt.close(fig)


HTML = '''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>VirtualTissue-scRNAseq · Run explorer</title>
<style>
:root{font:15px system-ui;color:#172b40;background:#f3f6fa}body{margin:0}main{max-width:1160px;margin:auto;padding:32px 24px}h1{font-size:30px;margin:6px 0}h2{font-size:17px;margin-top:0}p{line-height:1.5}.kicker{letter-spacing:.12em;font-size:11px;font-weight:700;color:#327c76}.note{border-left:4px solid #d99d39;background:#fff5df;padding:12px 16px}.card{background:white;border:1px solid #dce3eb;border-radius:12px;padding:20px;margin:18px 0}.grid{display:grid;grid-template-columns:2fr 1fr;gap:18px}.metric{font-size:26px;font-weight:700}.muted{color:#596b7e;font-size:13px}button,select{padding:9px 14px;border:1px solid #b8c7d7;border-radius:6px;background:white;color:#172b40;cursor:pointer}button:hover{background:#eaf1f9}input[type=range]{width:75%;vertical-align:middle}.legend{display:flex;gap:14px;flex-wrap:wrap;font-size:12px;margin:14px 0}.swatch{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:5px}canvas{width:100%;height:290px;display:block}#agents{height:340px;cursor:crosshair}.metrics{display:flex;gap:35px;flex-wrap:wrap}table{width:100%;border-collapse:collapse;font-size:12px}td,th{text-align:left;padding:8px 5px;border-bottom:1px solid #e9eef3}.scroll{max-height:290px;overflow:auto}code{word-break:break-all}a{color:#236aab}@media(max-width:800px){.grid{display:block}main{padding:20px 12px}}
</style><main><div class="kicker">VIRTUALTISSUE · SINGLE CELL</div><h1>Cell-state experiment</h1>
<p class="muted">Inspect the starting cells, declared assumptions and every simulated event. This file works offline.</p>
<div class="note"><strong>Uncalibrated software experiment.</strong> Expression initializes labels only. Rates are supplied assumptions. The cell grid is an index layout, not tissue geometry. No treatment efficacy is inferred.</div>
<section class="card"><div class="metrics"><div><div class="muted">Observed time</div><div class="metric" id="clock"></div></div><div><div class="muted">Live agents</div><div class="metric" id="population"></div></div><div><div class="muted">Run status</div><div id="status"></div></div><div><div class="muted">Seed</div><div id="seed"></div></div></div>
<p><button id="play">Play</button> <input aria-label="Simulation timeline" id="time" type="range" min="0" step="1" value="0"></p><div class="legend" id="legend"></div></section>
<div class="grid"><section class="card"><h2>Population trajectory</h2><canvas id="curve" width="1000" height="380"></canvas><p class="muted">Recorded state counts at output times. Full-precision event times are available in events.tsv.</p></section>
<section class="card"><h2>Run provenance</h2><div id="details"></div></section></div>
<div class="grid"><section class="card"><h2>Live agents · click to inspect</h2><canvas id="agents" width="1000" height="450"></canvas><p class="muted" id="gridnote"></p></section>
<section class="card"><h2>Cell history</h2><p id="selection">Choose an agent from the grid.</p><div class="scroll"><table><thead><tr><th>Hours</th><th>Event</th><th>State</th></tr></thead><tbody id="history"></tbody></table></div></section></div>
<section class="card"><h2>Model assumptions</h2><p id="assumptions"></p><div class="scroll"><pre id="model"></pre></div></section>
<p class="muted">Generated by VirtualTissue-scRNAseq. Compare seeds and the no-events control before interpreting a run.</p>
</main><script id="payload" type="application/json">__PAYLOAD__</script><script>
const d=JSON.parse(document.getElementById('payload').textContent), states=d.model.states, colors=__COLORS__;
const el=id=>document.getElementById(id), color=s=>colors[states.indexOf(s)%colors.length];
const timeline=el('time');timeline.max=d.trajectory.length-1;
el('status').textContent=d.summary.stop_reason;el('seed').textContent=d.summary.seed;
el('assumptions').textContent=d.model.rate_provenance;el('model').textContent=JSON.stringify(d.model,null,2);
for(const state of states){const span=document.createElement('span'),dot=document.createElement('i');dot.className='swatch';dot.style.background=color(state);span.append(dot,document.createTextNode(state));el('legend').append(span);}
for(const [key,value] of Object.entries({Sample:d.provenance.sample_id,Annotation:d.provenance.annotation,Cells:d.initial_cells.length,Events:d.events.length,Version:d.software_version,License:d.provenance.source_license})){const p=document.createElement('p'),b=document.createElement('strong');b.textContent=key+': ';p.append(b,document.createTextNode(String(value)));el('details').append(p);}
let shown=[],selected=null, timer=null;
function aliveAt(t){const map=new Map(d.initial_cells.map(c=>[c.cell_id,c.state]));for(const e of d.events){if(e.time_hours>t)break;if(e.event==='death')map.delete(e.cell_id);else map.set(e.cell_id,e.to_state);}return map;}
function history(id,t){el('selection').textContent=id;el('history').replaceChildren();const initial=d.initial_cells.find(c=>c.cell_id===id);const entries=[];if(initial)entries.push({time_hours:0,event:'measured / supplied start',to_state:initial.state});entries.push(...d.events.filter(e=>(e.cell_id===id||e.parent_id===id)&&e.time_hours<=t));for(const e of entries){const row=document.createElement('tr');for(const text of [e.time_hours.toFixed(3),e.parent_id===id?'daughter created':e.event,e.to_state||'removed']){const td=document.createElement('td');td.textContent=text;row.append(td);}el('history').append(row);}}
function draw(){const row=d.trajectory[Number(timeline.value)],t=row.time_hours,map=aliveAt(t);el('clock').textContent=t.toFixed(2)+' h';el('population').textContent=map.size;
 const cv=el('curve'),c=cv.getContext('2d'),w=cv.width,h=cv.height,left=65,top=25,right=w-25,bottom=h-45,maxT=d.summary.observed_hours||1,maxY=Math.max(1,...d.trajectory.flatMap(r=>states.map(s=>r[s])));
 c.clearRect(0,0,w,h);c.font='18px system-ui';c.fillStyle='#596b7e';c.strokeStyle='#dbe3ec';
 for(let k=0;k<=4;k++){let y=bottom-(bottom-top)*k/4;c.beginPath();c.moveTo(left,y);c.lineTo(right,y);c.stroke();c.fillText(Math.round(maxY*k/4),5,y+5);c.fillText((maxT*k/4).toFixed(1),left+(right-left)*k/4-12,h-10);}
 for(const s of states){c.beginPath();c.strokeStyle=color(s);c.lineWidth=3;d.trajectory.forEach((r,i)=>{const x=left+r.time_hours/maxT*(right-left),y=bottom-r[s]/maxY*(bottom-top);if(i===0)c.moveTo(x,y);else{const previous=d.trajectory[i-1];c.lineTo(x,bottom-previous[s]/maxY*(bottom-top));c.lineTo(x,y);}});c.stroke();}
 c.strokeStyle='#172b40';c.lineWidth=1;c.setLineDash([5,5]);const marker=left+t/maxT*(right-left);c.beginPath();c.moveTo(marker,top);c.lineTo(marker,bottom);c.stroke();c.setLineDash([]);
 const av=el('agents'),a=av.getContext('2d');a.clearRect(0,0,av.width,av.height);const entries=Array.from(map).slice(0,2000),columns=Math.max(1,Math.ceil(Math.sqrt(entries.length*av.width/av.height))),size=Math.min(30,av.width/columns,av.height/Math.max(1,Math.ceil(entries.length/columns)));
 shown=entries.map(([id,state],i)=>({id,state,x:(i%columns+.5)*size,y:(Math.floor(i/columns)+.5)*size,r:Math.max(1,size*.34)}));for(const dot of shown){a.fillStyle=color(dot.state);a.beginPath();a.arc(dot.x,dot.y,dot.r,0,Math.PI*2);a.fill();if(dot.id===selected){a.strokeStyle='#172b40';a.lineWidth=2;a.stroke();}}
 el('gridnote').textContent='Showing '+shown.length+' of '+map.size+' live agents. Grid position has no biological meaning.';if(selected)history(selected,t);
}
timeline.addEventListener('input',draw);el('agents').addEventListener('click',e=>{const rect=e.target.getBoundingClientRect(),x=(e.clientX-rect.left)*e.target.width/rect.width,y=(e.clientY-rect.top)*e.target.height/rect.height;const dot=shown.find(p=>(p.x-x)**2+(p.y-y)**2<=(p.r+4)**2);if(dot){selected=dot.id;draw();}});
el('play').addEventListener('click',()=>{if(timer){clearInterval(timer);timer=null;el('play').textContent='Play';return;}el('play').textContent='Pause';timer=setInterval(()=>{if(Number(timeline.value)>=Number(timeline.max)){clearInterval(timer);timer=null;el('play').textContent='Play';return;}timeline.value=Number(timeline.value)+1;draw();},180);});draw();
</script></html>'''


def write_html(result, path):
    # Prevent source labels containing HTML/script terminators from escaping JSON.
    payload = json.dumps(result, ensure_ascii=True, allow_nan=False).replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')
    Path(path).write_text(HTML.replace('__PAYLOAD__', payload).replace('__COLORS__', json.dumps(COLORS)), encoding='utf-8')


def run_bundle(bundle, output, *, model_path=None, scenario='baseline', **options):
    path = Path(output)
    zip_path = path.parent / (path.name + '.zip')
    if path.exists() or zip_path.exists():
        raise FileExistsError('Run directory or results ZIP already exists. Choose a new run name.')
    cells, embedded, provenance = read_bundle(bundle)
    model = load_model(model_path, scenario) if model_path else embedded
    if scenario == 'no-events':
        for edge in model['transitions']:
            edge['rate'] = 0.0
        for key in ['division_rates', 'death_rates', 'influx_rates']:
            model[key] = dict.fromkeys(model['states'], 0.0)
    result = simulate(cells, model, **options)
    result.update({'model': model, 'provenance': provenance, 'initial_cells': cells,
                   'software_version': __version__, 'scenario': scenario,
                   'bundle_sha256': file_hashes(bundle)['sha256'],
                   'environment': {'python': sys.version.split()[0], 'platform': platform.system(), 'numpy': np.__version__}})
    path.mkdir(parents=True)
    (path / 'run.json').write_bytes(json_bytes(result))
    for name in ['trajectory', 'mean_field']:
        (path / (name + '.tsv')).write_bytes(table_bytes(result[name], ['time_hours', *model['states']]))
    (path / 'events.tsv').write_bytes(table_bytes(result['events'], ['event_id', 'time_hours', 'event', 'cell_id', 'parent_id', 'from_state', 'to_state']))
    (path / 'final_cells.tsv').write_bytes(table_bytes(result['final_cells'], ['cell_id', 'state']))
    (path / 'summary.json').write_bytes(json_bytes(result['summary']))
    (path / 'model.yaml').write_text(__import__('yaml').safe_dump(model, sort_keys=False))
    plot_run(result, path / 'overview.png')
    write_html(result, path / 'explorer.html')
    manifest = {'schema_version': 1, 'kind': 'vtsc_results',
                'files': {p.name: file_hashes(p)['sha256'] for p in sorted(path.iterdir()) if p.is_file()}}
    (path / 'manifest.json').write_bytes(json_bytes(manifest))
    shutil.make_archive(str(path), 'zip', root_dir=path)
    return {'output': str(path), 'results_zip': str(zip_path), **result['summary']}


def read_run(path):
    path = Path(path)
    if path.is_dir():
        manifest = json.loads((path / 'manifest.json').read_text())
        def read(name):
            if Path(name).name != name or name in ['.', '..']:
                raise ValueError('Invalid result path.')
            return (path / name).read_bytes()
    else:
        with zipfile.ZipFile(path) as archive:
            infos = archive.infolist()
            if len({i.filename for i in infos}) != len(infos) or sum(i.file_size for i in infos) > 400_000_000:
                raise ValueError('Duplicate result members or results ZIP exceeds 400 MB.')
            payload = {i.filename: archive.read(i) for i in infos}
        manifest = json.loads(payload['manifest.json'])
        read = lambda name: payload[name]
    if manifest.get('schema_version') != 1 or manifest.get('kind') != 'vtsc_results' or 'run.json' not in manifest.get('files', {}):
        raise ValueError('Unsupported results manifest.')
    for name, checksum in manifest['files'].items():
        if hashlib.sha256(read(name)).hexdigest() != checksum:
            raise ValueError(f'Result checksum mismatch: {name}')
    return json.loads(read('run.json'))


def compare_runs(first, second, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    a, b = read_run(first), read_run(second)
    if a.get('kind') in ['vtsc_tissue_run','vtsc_melanoma_run','vtsc_melanoma_comparison'] or b.get('kind') in ['vtsc_tissue_run','vtsc_melanoma_run','vtsc_melanoma_comparison']:
        raise ValueError('compare is for legacy nonspatial runs; compare tissue trajectory.tsv and summary.json or use the notebooks.')
    if a['bundle_sha256'] != b['bundle_sha256'] or a['model']['states'] != b['model']['states']:
        raise ValueError('Comparison requires the same initial bundle and state vocabulary.')
    end = min(a['summary']['observed_hours'], b['summary']['observed_hours'])
    path = Path(output)
    if path.exists():
        raise FileExistsError(path)
    path.mkdir(parents=True)
    fig, ax = plt.subplots(figsize=(8, 4), layout='constrained')
    for run, style, label in [(a, '-', 'A'), (b, '--', 'B')]:
        rows = [r for r in run['trajectory'] if r['time_hours'] <= end]
        for i, state in enumerate(run['model']['states']):
            ax.step([r['time_hours'] for r in rows], [r[state] for r in rows], where='post',
                    linestyle=style, color=COLORS[i % len(COLORS)], label=f'{state} · {label}')
    ax.set(xlabel='Simulated hours', ylabel='Simulated cells', title='Shared observed interval · illustrative rates')
    ax.legend(fontsize=8); fig.savefig(path / 'comparison.png', dpi=150); plt.close(fig)
    result = {'shared_hours': end, 'same_initial_bundle': True, 'seed_a': a['summary']['seed'],
              'seed_b': b['summary']['seed'], 'scenario_a': a['scenario'], 'scenario_b': b['scenario'],
              'interpretation': 'Descriptive software comparison, not a biological treatment effect or statistical test.'}
    (path / 'comparison.json').write_bytes(json_bytes(result))
    return result
