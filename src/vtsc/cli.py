"""Project status, integrity registration and immutable pre-analysis freeze."""

import argparse
import subprocess
import json
import sys
from zipfile import BadZipFile
import yaml
from pathlib import Path

from .provenance import file_hashes, inspect_download, utc_now, verify_freeze, write_json


def _main():
    parser = argparse.ArgumentParser(description='Melanoma simulation from dissociated scRNA-seq and independent research audits')
    parser.add_argument("--root", type=Path, default=Path.cwd())
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("status")
    sub.add_parser("verify-freeze")
    register = sub.add_parser("register")
    register.add_argument("filename")
    register.add_argument("--data-root", type=Path, required=True)
    freeze = sub.add_parser("freeze")
    freeze.add_argument("--file", action="append", required=True)
    freeze.add_argument("--commit", required=True)
    demo = sub.add_parser('demo', help='Offline dynamic CD8 melanoma with paired anti-PD-1/exposure-zero arms')
    demo.add_argument('output', type=Path)
    demo.add_argument('--decision-mode', choices=['local','mock-jev'], default='local')
    demo.add_argument('--cells', type=int, default=120)
    demo.add_argument('--seed', type=int, default=21)
    demo.add_argument('--hours', type=float, default=24)
    synth = sub.add_parser('synthetic', help='Generate a raw-count annotated H5AD fixture')
    synth.add_argument('output', type=Path)
    synth.add_argument('--cells', type=int, default=300)
    synth.add_argument('--seed', type=int, default=7)
    download = sub.add_parser('download-example', help='Download checksum-pinned CC BY 4.0 dissociated melanoma')
    download.add_argument('output', type=Path)
    download.add_argument('--dataset', choices=['melanoma', 'pbmc3k-legacy'], default='melanoma')
    tissue_template = sub.add_parser('tissue-template', help='Write declared tissue assumptions')
    tissue_template.add_argument('output', type=Path)
    tissue_prepare = sub.add_parser('prepare-tissue', help='Prepare one melanoma biopsy with its microenvironment')
    tissue_prepare.add_argument('input', type=Path)
    tissue_prepare.add_argument('output', type=Path)
    tissue_prepare.add_argument('--metadata', type=Path)
    tissue_prepare.add_argument('--layer')
    tissue_prepare.add_argument('--state-column', default='cd8_state')
    tissue_prepare.add_argument('--state-recipe', choices=['rank-markers-v1'])
    tissue_prepare.add_argument('--gene-map', type=Path)
    tissue_prepare.add_argument('--sampling', choices=['stratified','uniform'], default='stratified')
    tissue_prepare.add_argument('--cell-type-column', default='cell_type')
    tissue_prepare.add_argument('--sample-id')
    tissue_prepare.add_argument('--recipe', choices=['melanoma-markers'])
    tissue_prepare.add_argument('--max-cells', type=int, default=400)
    tissue_prepare.add_argument('--seed', type=int, default=21)
    tissue_prepare.add_argument('--source-license', default='user supplied')
    tissue_run = sub.add_parser('run-tissue', help='Run interacting agents on assumed 2D geometry')
    tissue_run.add_argument('bundle', type=Path)
    tissue_run.add_argument('output', type=Path)
    tissue_run.add_argument('--scenario', choices=['recognition-off', 'declared-recognition', 'frozen'], default='recognition-off')
    tissue_run.add_argument('--layout', choices=['tumor-nests', 'mixed'], default='tumor-nests')
    tissue_run.add_argument('--config', type=Path)
    tissue_run.add_argument('--hours', type=float, default=24)
    tissue_run.add_argument('--seed', type=int, default=21)
    tissue_run.add_argument('--max-cells', type=int, default=1500)
    dynamic_template = sub.add_parser('melanoma-template', help='Write complete versioned CD8/drug/mechanics contract')
    dynamic_template.add_argument('output', type=Path)
    dynamic = sub.add_parser('run-melanoma', help='Dynamic CD8 states and explicit normalized anti-PD-1 exposure')
    dynamic.add_argument('bundle', type=Path)
    dynamic.add_argument('output', type=Path)
    dynamic.add_argument('--arm', choices=['anti-pd1','exposure-zero','recognition-off','checkpoint-off','transitions-off','spatial-shuffle','frozen'], default='anti-pd1')
    dynamic.add_argument('--config', type=Path)
    dynamic.add_argument('--decision-mode', choices=['local','mock-jev','jev'], default='local')
    dynamic.add_argument('--max-http-attempts', type=int)
    dynamic.add_argument('--allow-paid', action='store_true')
    dynamic.add_argument('--jev-model', default='jev-1.14.0')
    dynamic.add_argument('--replay-journal', type=Path)
    check_dynamic=sub.add_parser('check-melanoma', help='Zero-cost first-round Jev payload and activity preview')
    check_dynamic.add_argument('bundle', type=Path)
    check_dynamic.add_argument('output', type=Path)
    check_dynamic.add_argument('--config', type=Path)
    check_dynamic.add_argument('--seed', type=int, default=21)
    check_dynamic.add_argument('--arm', choices=['anti-pd1','exposure-zero','recognition-off'], default='anti-pd1')
    check_dynamic.add_argument('--layout', choices=['mixed','tumor-nests'], default='tumor-nests')
    check_dynamic.add_argument('--jev-model', default='jev-1.14.0')
    dynamic.add_argument('--hours', type=float, default=24)
    dynamic.add_argument('--seed', type=int, default=21)
    dynamic.add_argument('--max-cells', type=int, default=1500)
    dynamic.add_argument('--layout', choices=['mixed','tumor-nests'], default='tumor-nests')
    matched = sub.add_parser('compare-melanoma', help='Synchronized offline inspection of matched arms')
    matched.add_argument('first', type=Path)
    matched.add_argument('second', type=Path)
    matched.add_argument('output', type=Path)
    recovery = sub.add_parser('recover-melanoma', help='Recompute a verified recording deterministically into a new output; optionally extend')
    recovery.add_argument('input', type=Path)
    recovery.add_argument('output', type=Path)
    recovery.add_argument('--hours', type=float)
    ens = sub.add_parser('ensemble-melanoma', help='Paired controls and predeclared parameter sensitivity')
    ens.add_argument('bundle', type=Path)
    ens.add_argument('output', type=Path)
    ens.add_argument('--seeds', type=int, nargs='+', default=[21,22,23,24,25])
    ens.add_argument('--hours', type=float, default=24)
    ens.add_argument('--config', type=Path)
    template = sub.add_parser('model-template', help='Write an editable assumption file')
    template.add_argument('output', type=Path)
    for name in ['validate', 'prepare']:
        command = sub.add_parser(name, help='Validate raw counts' if name == 'validate' else 'Prepare a verified input bundle')
        command.add_argument('input', type=Path)
        command.add_argument('--metadata', type=Path)
        command.add_argument('--layer')
        if name == 'prepare':
            command.add_argument('output', type=Path)
            command.add_argument('--model', type=Path)
            command.add_argument('--state-column', default='state')
            command.add_argument('--sample-column', default='sample_id')
            command.add_argument('--sample-id')
            command.add_argument('--recipe', choices=['pbmc-cd8-demo'])
            command.add_argument('--max-cells', type=int)
            command.add_argument('--seed', type=int, default=21)
            command.add_argument('--source-license', default='user-supplied; redistribution rights not assessed')
    for name in ['check', 'run']:
        command = sub.add_parser(name, help='Inspect an input bundle' if name == 'check' else 'Run the exact non-spatial cell-state process')
        command.add_argument('bundle', type=Path)
        command.add_argument('--model', type=Path)
        command.add_argument('--scenario', choices=['baseline', 'no-events'], default='baseline')
        command.add_argument('--hours', type=float, default=48)
        if name == 'run':
            command.add_argument('output', type=Path)
            command.add_argument('--seed', type=int, default=21)
            command.add_argument('--interval', type=float, default=1)
            command.add_argument('--max-cells', type=int, default=5000)
            command.add_argument('--max-events', type=int, default=50000)
    explore = sub.add_parser('explore', help='Verify saved results; optionally regenerate the offline viewer')
    explore.add_argument('input', type=Path)
    explore.add_argument('--html', type=Path)
    explore.add_argument('--open', action='store_true')
    compare = sub.add_parser('compare', help='Compare runs from the same bundle on their shared interval')
    compare.add_argument('first', type=Path)
    compare.add_argument('second', type=Path)
    compare.add_argument('output', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    if args.command in ['check-melanoma','run-melanoma','melanoma-template','compare-melanoma','recover-melanoma','ensemble-melanoma','prepare-tissue', 'run-tissue', 'tissue-template', 'demo', 'synthetic', 'download-example', 'model-template', 'validate', 'prepare', 'check', 'run', 'explore', 'compare']:
        result = software_command(args)
        print(json.dumps(result, indent=2, allow_nan=False))
        return
    if args.command == "status":
        print((root / "reports/STATUS.md").read_text())
        return
    if args.command == "verify-freeze":
        payload = json.loads((root / "prereg/FREEZE.json").read_text())
        failures = verify_freeze(root, payload)
        print(json.dumps({"valid": not failures, "mismatches": failures}))
        raise SystemExit(bool(failures))
    if args.command == "freeze":
        target = root / "prereg/FREEZE.json"
        if target.exists():
            parser.error("Freeze already exists. Record an amendment; never overwrite it.")
        records = []
        for name in args.file:
            path = (root / name).resolve()
            if not path.is_relative_to(root):
                parser.error("Freeze inputs must be inside the repository")
            records.append({"path": str(path.relative_to(root)), "sha256": file_hashes(path)["sha256"]})
        write_json(target, {"schema_version": 1, "stage": "T0.4", "created_at_utc": utc_now(),
                            "source_commit": args.commit, "public_registration": False, "files": records})
        print("Created prereg/FREEZE.json (local freeze, not public registration).")
        return
    manifest_path = root / "data/MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    records = [r for r in manifest["files"] if r["filename"] == args.filename]
    if len(records) != 1:
        parser.error("Filename must identify one manifest record")
    record = records[0]
    data_root = args.data_root.resolve()
    path = (data_root / record["relative_path"]).resolve()
    if not path.is_relative_to(data_root):
        parser.error("Data path escapes data root")
    result = inspect_download(path, record.get("expected_bytes"), record.get("source_checksum"))
    # Never silently replace a previously recorded local hash after modification.
    previous_hash = record.get("local", {}).get("sha256")
    if previous_hash and previous_hash != result.get("sha256"):
        parser.error("Registered file has changed; investigate and log a deviation")
    record["local"] = result
    write_json(manifest_path, manifest)
    print(json.dumps({"filename": args.filename, **result}))
    raise SystemExit(result["status"] not in {"verified", "hashed_only"})


def software_command(args):
    from .inputs import synthetic_example, download_pbmc, load_input, prepare, read_bundle
    from .model import load_model, mean_field
    from .results import run_bundle, read_run, write_html, compare_runs
    import yaml
    if args.command == 'prepare-tissue':
        from .tissue_inputs import prepare_tissue
        return prepare_tissue(args.input, args.output, metadata=args.metadata, layer=args.layer,
                              cell_type_column=args.cell_type_column, sample_id=args.sample_id,
                              recipe=args.recipe, max_cells=args.max_cells, seed=args.seed, source_license=args.source_license,
                              state_column=args.state_column, state_recipe=args.state_recipe, gene_map=args.gene_map, sampling=args.sampling)
    if args.command == 'run-tissue':
        from .tissue_results import run_tissue
        return run_tissue(args.bundle, args.output, scenario=args.scenario, layout=args.layout,
                          config_path=args.config, hours=args.hours, seed=args.seed, max_cells=args.max_cells)
    if args.command == 'tissue-template':
        from .tissue import DEFAULT
        if args.output.exists(): raise FileExistsError(args.output)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(yaml.safe_dump(DEFAULT, sort_keys=False))
        return {'config': str(args.output), 'rates_fitted': False}
    if args.command == 'demo':
        from .melanoma_results import melanoma_demo
        return melanoma_demo(args.output, cells=args.cells, hours=args.hours, seed=args.seed,decision_mode=args.decision_mode)
    if args.command == 'melanoma-template':
        from .melanoma import DEFAULT_CONFIG
        if args.output.exists(): raise FileExistsError(args.output)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(yaml.safe_dump(DEFAULT_CONFIG, sort_keys=False))
        return {'config':str(args.output),'evidence_class':'assumed','clinical_validation':False}
    if args.command == 'check-melanoma':
        from .melanoma_results import preflight_melanoma
        return preflight_melanoma(args.bundle,args.output,config_path=args.config,seed=args.seed,arm=args.arm,layout=args.layout,jev_model=args.jev_model)
    if args.command == 'run-melanoma':
        from .melanoma_results import run_melanoma
        return run_melanoma(args.bundle,args.output,config_path=args.config,arm=args.arm,layout=args.layout,hours=args.hours,seed=args.seed,max_cells=args.max_cells,decision_mode=args.decision_mode,max_http_attempts=args.max_http_attempts,allow_paid=args.allow_paid,jev_model=args.jev_model,replay_journal=args.replay_journal)
    if args.command == 'compare-melanoma':
        from .melanoma_results import compare_melanoma
        return compare_melanoma(args.first,args.second,args.output)
    if args.command == 'recover-melanoma':
        from .melanoma_results import recover_melanoma
        return recover_melanoma(args.input,args.output,hours=args.hours)
    if args.command == 'ensemble-melanoma':
        from .melanoma_results import ensemble_melanoma
        return ensemble_melanoma(args.bundle,args.output,seeds=args.seeds,hours=args.hours,config_path=args.config)
    if args.command == 'model-template':
        if args.output.exists():
            raise FileExistsError(args.output)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(yaml.safe_dump(load_model(), sort_keys=False))
        return {'model': str(args.output), 'rates_fitted': False}
    if args.command == 'synthetic':
        return synthetic_example(args.output, cells=args.cells, seed=args.seed)
    if args.command == 'download-example':
        from .tissue_inputs import download_melanoma
        return download_pbmc(args.output) if args.dataset == 'pbmc3k-legacy' else download_melanoma(args.output)
    if args.command == 'validate':
        data, sources = load_input(args.input, args.metadata, args.layer)
        return {'valid_counts': True, 'cells': data.n_obs, 'features': data.n_vars,
                'metadata_columns': list(data.obs.columns), 'sources': sources}
    if args.command == 'prepare':
        return prepare(args.input, args.output, metadata=args.metadata, layer=args.layer, model_path=args.model,
                       state_column=args.state_column, sample_column=args.sample_column, sample_id=args.sample_id,
                       recipe=args.recipe, max_cells=args.max_cells, seed=args.seed, source_license=args.source_license)
    if args.command == 'check':
        import math
        cells, model, provenance = read_bundle(args.bundle)
        if not math.isfinite(args.hours) or not 0 < args.hours <= 10000:
            raise ValueError('hours must lie in (0, 10000].')
        if args.model:
            model = load_model(args.model, args.scenario)
        if args.scenario == 'no-events':
            for e in model['transitions']: e['rate'] = 0
            for key in ['division_rates', 'death_rates', 'influx_rates']:
                model[key] = dict.fromkeys(model['states'], 0.0)
        if any(c['state'] not in model['states'] for c in cells):
            raise ValueError('Bundle states are incompatible with this model.')
        counts = {s: sum(c['state'] == s for c in cells) for s in model['states']}
        expected = mean_field(list(counts.values()), [args.hours], model)[0]
        return {'valid_bundle': True, 'initial_cells': len(cells), 'initial_states': counts,
                'uncapped_expected_final_cells': float(expected.sum()), 'model': model, 'provenance': provenance}
    if args.command == 'run':
        return run_bundle(args.bundle, args.output, model_path=args.model, scenario=args.scenario, hours=args.hours,
                          interval=args.interval, seed=args.seed, max_cells=args.max_cells, max_events=args.max_events)
    if args.command == 'compare':
        return compare_runs(args.first, args.second, args.output)
    result = read_run(args.input)
    if result.get('kind') in ['vtsc_melanoma_run','vtsc_melanoma_comparison']:
        from .melanoma_results import verify_dynamic
        verify_dynamic(result)
    target = args.html
    if target:
        if target.exists():
            raise FileExistsError(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        if result.get('kind') in ['vtsc_melanoma_run','vtsc_melanoma_comparison']:
            from .melanoma_results import write_melanoma_html
            write_melanoma_html(result,target)
        elif result.get('kind') == 'vtsc_tissue_run':
            from .tissue_results import write_tissue_html
            write_tissue_html(result, target)
        else:
            write_html(result, target)
    elif args.input.is_dir():
        target = args.input / 'explorer.html'
    if args.open:
        if target is None:
            raise ValueError('For a results ZIP, specify --html before --open.')
        import webbrowser
        webbrowser.open(target.resolve().as_uri())
    return {'verified': True, 'summary': result.get('summary',result.get('comparison')), 'explorer': str(target) if target else None}


def main():
    try:
        _main()
    except (ValueError, OSError, KeyError, BadZipFile, yaml.YAMLError, subprocess.CalledProcessError) as error:
        print(f'vtsc: {error}', file=sys.stderr)
        raise SystemExit(2) from None


if __name__ == '__main__':
    main()
