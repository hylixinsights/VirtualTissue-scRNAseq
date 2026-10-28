"""Project status, integrity registration and immutable pre-analysis freeze."""

import argparse
import json
import sys
from zipfile import BadZipFile
import yaml
from pathlib import Path

from .provenance import file_hashes, inspect_download, utc_now, verify_freeze, write_json


def _main():
    parser = argparse.ArgumentParser(description=__doc__)
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
    demo = sub.add_parser('demo', help='Offline synthetic end-to-end software test')
    demo.add_argument('output', type=Path)
    demo.add_argument('--cells', type=int, default=300)
    demo.add_argument('--seed', type=int, default=21)
    demo.add_argument('--hours', type=float, default=48)
    synth = sub.add_parser('synthetic', help='Generate a raw-count annotated H5AD fixture')
    synth.add_argument('output', type=Path)
    synth.add_argument('--cells', type=int, default=300)
    synth.add_argument('--seed', type=int, default=7)
    download = sub.add_parser('download-example', help='Download checksum-pinned CC BY 4.0 PBMC3k')
    download.add_argument('output', type=Path)
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
    if args.command in ['demo', 'synthetic', 'download-example', 'model-template', 'validate', 'prepare', 'check', 'run', 'explore', 'compare']:
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
    if args.command == 'model-template':
        if args.output.exists():
            raise FileExistsError(args.output)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(yaml.safe_dump(load_model(), sort_keys=False))
        return {'model': str(args.output), 'rates_fitted': False}
    if args.command == 'synthetic':
        return synthetic_example(args.output, cells=args.cells, seed=args.seed)
    if args.command == 'download-example':
        return download_pbmc(args.output)
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
    if args.command == 'demo':
        if args.output.exists():
            raise FileExistsError('Demo output already exists. Use a new directory.')
        synthetic_example(args.output / 'synthetic.h5ad', cells=args.cells, seed=args.seed)
        bundle = args.output / 'input.zip'
        prepare(args.output / 'synthetic.h5ad', bundle, source_license='MIT; generated synthetic fixture, no patient data')
        baseline = run_bundle(bundle, args.output / 'baseline', hours=args.hours, seed=args.seed)
        control = run_bundle(bundle, args.output / 'no_events', hours=args.hours, seed=args.seed, scenario='no-events')
        compare_runs(args.output / 'baseline', args.output / 'no_events', args.output / 'comparison')
        return {'baseline': baseline, 'control': control, 'explorer': str(args.output / 'baseline/explorer.html')}
    if args.command == 'compare':
        return compare_runs(args.first, args.second, args.output)
    result = read_run(args.input)
    target = args.html
    if target:
        if target.exists():
            raise FileExistsError(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        write_html(result, target)
    elif args.input.is_dir():
        target = args.input / 'explorer.html'
    if args.open:
        if target is None:
            raise ValueError('For a results ZIP, specify --html before --open.')
        import webbrowser
        webbrowser.open(target.resolve().as_uri())
    return {'verified': True, 'summary': result['summary'], 'explorer': str(target) if target else None}


def main():
    try:
        _main()
    except (ValueError, OSError, KeyError, BadZipFile, yaml.YAMLError) as error:
        print(f'vtsc: {error}', file=sys.stderr)
        raise SystemExit(2) from None


if __name__ == '__main__':
    main()
