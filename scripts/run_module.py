#!/usr/bin/env python3
"""Run explicit RESIST steps with preflight, logs, and output provenance."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import uuid

ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('module', choices=list('ABCD'))
    parser.add_argument('--steps', help='Comma-separated step IDs; use --list for choices')
    parser.add_argument('--config', type=Path, default=ROOT / 'config/config.yaml')
    parser.add_argument('--apa-config', type=Path, default=ROOT / 'config/apa_config.yaml')
    parser.add_argument('--check', action='store_true', help='Validate only; do not run analyses')
    parser.add_argument('--list', action='store_true', help='Show steps and expected products')
    args = parser.parse_args()
    manifest = json.loads((ROOT / 'config/steps.json').read_text())
    if args.list:
        print('Default steps:', ', '.join(manifest['defaults'][args.module]))
        for key, spec in manifest['steps'].items():
            if key.startswith(args.module):
                print(f"{key}: {spec['title']}\n  Requires: {spec['requires']}\n  Outputs: {spec['outputs']}")
        return 0
    steps = [v.strip() for v in args.steps.split(',')] if args.steps else manifest['defaults'][args.module]
    if not steps or len(set(steps)) != len(steps):
        parser.error('Choose at least one step and do not repeat IDs.')
    for step in steps:
        if step not in manifest['steps'] or not step.startswith(args.module):
            parser.error(f'{step!r} is not a runnable step in Module {args.module}; use --list.')
    # Respect the user's order, but reject reversed producer/consumer steps.
    for child, parent in [('B2','B1'), ('B3','B2'), ('B6','B1'), ('B7','B6'), ('C8','C7'), ('C10','C9')]:
        if child in steps and parent in steps and steps.index(child) < steps.index(parent):
            parser.error(f'{parent} must precede {child}.')
    if shutil.which('Rscript') is None:
        parser.error('Rscript is missing. See config/setup/ and TUTORIAL.md, Step 1.')
    config = args.config.expanduser().resolve()
    apa = args.apa_config.expanduser().resolve()
    if not config.is_file():
        parser.error(f'Configuration does not exist: {config}')
    env = dict(os.environ, RESIST_HOME=str(ROOT), RESIST_CONFIG=str(config), RESIST_APA_CONFIG=str(apa))
    env['RESIST_STEPS'] = ','.join(steps)
    print(f"RESIST v3_02 | Module {args.module} | Steps: {', '.join(steps)}", flush=True)
    preflight = subprocess.run(['Rscript', str(ROOT/'scripts/preflight.R')], env=env,
                               cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    # The final machine-readable path line is emitted only after validation passes.
    lines = preflight.stdout.splitlines()
    output_line = next((s for s in lines if s.startswith('RESIST_OUTPUT_ROOT=')), None)
    print('\n'.join(s for s in lines if not s.startswith('RESIST_OUTPUT_ROOT=')), flush=True)
    if preflight.returncode:
        return preflight.returncode
    if args.check:
        print('CHECK PASSED. No analyses were run.')
        return 0
    if not output_line:
        print('Preflight did not return an output directory.', file=sys.stderr)
        return 1
    output = Path(output_line.split('=', 1)[1])
    run_id = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + args.module + '-' + uuid.uuid4().hex[:6]
    logdir = output/'logs'/run_id
    logdir.mkdir(parents=True)
    report = {'version':'v3_02', 'run_id':run_id, 'steps':steps, 'config':str(config),
              'config_sha256':digest(config), 'started_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
              'status':'running', 'step_results':[]}
    shutil.copy2(config, logdir/'config.yaml')
    if args.module == 'D':
        shutil.copy2(apa, logdir/'apa_config.yaml')
        report['apa_config_sha256'] = digest(apa)
    (logdir/'preflight.txt').write_text(preflight.stdout)
    report_path = logdir/'run.json'
    report_path.write_text(json.dumps(report, indent=2))
    rc = 0
    for step in steps:
        spec = manifest['steps'][step]
        script = ROOT/spec['script']
        outdir = output/args.module
        before = {str(p):(p.stat().st_mtime_ns,p.stat().st_size) for p in outdir.rglob('*') if p.is_file()} if outdir.exists() else {}
        print(f"[{step}] {spec['title']}", flush=True)
        log = logdir/(step+'.log')
        with log.open('w') as stream:
            process = subprocess.Popen(['Rscript',str(script)], cwd=ROOT, env=env, stdout=subprocess.PIPE,
                                       stderr=subprocess.STDOUT, text=True)
            for line in process.stdout:
                print(line, end='', flush=True)
                stream.write(line)
            rc = process.wait()
        produced = []
        for p in sorted(outdir.rglob('*')):
            if p.is_file() and before.get(str(p)) != (p.stat().st_mtime_ns,p.stat().st_size):
                produced.append({'path':str(p.relative_to(output)), 'bytes':p.stat().st_size, 'sha256':digest(p)})
        status = 'failed' if rc else ('completed' if produced else 'no_new_outputs')
        report['step_results'].append({'step':step,'script_sha256':digest(script),'exit_code':rc,
                                       'status':status,'outputs':produced})
        report_path.write_text(json.dumps(report, indent=2))
        print(f"[{step}] {status}: {len(produced)} files written; log {log}", flush=True)
        if not rc and not produced:
            print('  This step wrote no new files. Check the log for eligibility or significance filters.')
        if rc:
            break
    report['status'] = 'failed' if rc else ('completed_with_no_output_steps' if any(s['status']=='no_new_outputs' for s in report['step_results']) else 'completed')
    report['finished_utc'] = dt.datetime.now(dt.timezone.utc).isoformat()
    report_path.write_text(json.dumps(report, indent=2))
    print(f"Run status: {report['status']}\nResults: {output/args.module}\nRecord: {report_path}")
    return rc

if __name__ == '__main__':
    sys.exit(main())
