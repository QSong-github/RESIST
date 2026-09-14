#!/usr/bin/env python3
"""Run one RESIST module with step logs and a record of the installed software."""
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
    version = (ROOT / "VERSION").read_text(encoding="utf-8").strip()
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=f"RESIST {version}")
    parser.add_argument('module', choices=list('ABCD'))
    parser.add_argument('--steps', help='Comma-separated step IDs; use --list for choices')
    parser.add_argument('--config', type=Path, default=ROOT / 'config/config.yaml')
    parser.add_argument('--apa-config', type=Path, default=ROOT / 'config/apa_config.yaml')
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
    print(f"RESIST {version} | Module {args.module} | Steps: {', '.join(steps)}", flush=True)
    # Read the result path and software metadata only; do not load study data.
    # R reads YAML itself so path handling agrees with the analysis configuration.
    setup = subprocess.run(['Rscript', '--vanilla', '-e', r"""
        prefix <- Sys.getenv('CONDA_PREFIX')
        if (nzchar(prefix)) {
            if (!identical(normalizePath(R.home()),
                normalizePath(file.path(prefix, 'lib/R'), mustWork=FALSE)))
                stop('Rscript is outside the active conda environment. Activate resist.')
            .libPaths(.Library, include.site=FALSE)
        }
        stopifnot(requireNamespace('yaml', quietly=TRUE),
                  requireNamespace('jsonlite', quietly=TRUE))
        cfg <- yaml::read_yaml(Sys.getenv('RESIST_CONFIG'))
        path <- cfg$paths$results
        if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path))
            stop('Configuration must define one nonempty paths.results value.')
        if (grepl('^~', path)) stop('Use an absolute path instead of ~ in paths.results.')
        pkgs <- as.data.frame(installed.packages()[, c('Package', 'Version', 'LibPath')],
                              stringsAsFactors=FALSE)
        cc <- find.package('CellChat', quiet=TRUE)
        sha <- NA_character_
        if (length(cc)) {
            desc <- read.dcf(file.path(cc, 'DESCRIPTION'))
            if ('RemoteSha' %in% colnames(desc)) sha <- desc[1, 'RemoteSha']
        }
        cat(jsonlite::toJSON(list(results=path, CellChat_source=sha,
            session=capture.output(sessionInfo()), packages=pkgs,
            R_home=R.home(), libraries=.libPaths()), auto_unbox=TRUE))
    """], env=env, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if setup.returncode:
        print(setup.stderr, file=sys.stderr)
        return setup.returncode
    try:
        software = json.loads(setup.stdout)
        output = Path(software['results'])
    except (ValueError, KeyError, TypeError) as exc:
        print(f'Cannot read configuration/software metadata: {exc}', file=sys.stderr)
        return 1
    if env.get('CONDA_PREFIX'):
        env['R_LIBS_USER'] = str(Path(software['R_home']) / 'library')
        env['R_LIBS_SITE'] = env['R_LIBS_USER']
    if not output.is_absolute():
        output = ROOT / output
    output = output.resolve()
    if args.module == 'D' and not apa.is_file():
        parser.error(f'APA configuration does not exist: {apa}')
    run_id = dt.datetime.now(dt.timezone.utc).strftime('%Y%m%dT%H%M%SZ') + '-' + args.module + '-' + uuid.uuid4().hex[:6]
    logdir = output/'logs'/run_id
    logdir.mkdir(parents=True)
    report = {'version':version, 'run_id':run_id, 'steps':steps, 'config':str(config),
              'config_sha256':digest(config), 'started_utc':dt.datetime.now(dt.timezone.utc).isoformat(),
              'status':'running', 'step_results':[]}
    shutil.copy2(config, logdir/'config.yaml')
    if args.module == 'D':
        shutil.copy2(apa, logdir/'apa_config.yaml')
        report['apa_config_sha256'] = digest(apa)
    (logdir/'software.json').write_text(json.dumps(software, indent=2))
    report['software_record'] = 'software.json'
    report['dispatcher_sha256'] = digest(Path(__file__).resolve())
    report['manifest_sha256'] = digest(ROOT/'config/steps.json')
    report['shared_script_sha256'] = {str(p.relative_to(ROOT)): digest(p)
        for p in sorted((ROOT/'scripts/lib').glob('*.R'))}
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
            process = subprocess.Popen(['Rscript','--vanilla',str(script)], cwd=ROOT, env=env, stdout=subprocess.PIPE,
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
