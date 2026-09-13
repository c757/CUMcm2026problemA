from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

ROOT=Path(__file__).resolve().parent
SOURCES=('solve.py','run_all.py','make_figures.py','verify_moving_domain.py',
         'export_result2_full.py','test_sync.py','reproduce.py','requirements.txt',
         'A题.pdf','题目分析报告.md','术语表格.md')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--smoke',action='store_true')
    parser.add_argument('--compute-only',action='store_true')
    args=parser.parse_args()
    start=time.perf_counter()
    out=Path(tempfile.mkdtemp(prefix='reproduced_',dir=ROOT))
    for name in SOURCES:
        shutil.copy2(ROOT/name,out/name)
    for name in ('附件','utils','assets'):
        shutil.copytree(ROOT/name,out/name,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    if env.get('PYTHONPATH'):
        env['PYTHONPATH']=os.pathsep.join(str(Path(p or '.').resolve())
                                        for p in env['PYTHONPATH'].split(os.pathsep))
    print(f'Reproduction directory: {out}',flush=True)
    command=[sys.executable,'-B','test_sync.py'] if args.smoke else [sys.executable,'-B','run_all.py','--q2-horizon','full']
    if args.compute_only and not args.smoke:
        command.append('--compute-only')
    subprocess.run(command,cwd=out,env=env,check=True)
    from utils.repro_manifest import build_manifest
    inputs=([ROOT/n for n in SOURCES]+sorted((ROOT/'附件').rglob('*.xlsx'))
            +sorted((ROOT/'utils').glob('*.py'))+[p for p in sorted((ROOT/'assets').rglob('*')) if p.is_file()])
    parameters=(dict(n=80,rtol=1e-8,horizon_s=600,smoke=True)
                if args.smoke else dict(n=json.loads((out/'results/summary.json').read_text())['grid_n'],
                                        rtol=1e-8,q2_horizon='full',smoke=False,figures=not args.compute_only))
    suffix=' --smoke' if args.smoke else ' --compute-only' if args.compute_only else ''
    manifest=build_manifest(inputs,20260910,parameters,
                            'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B reproduce.py'+suffix,
                            ['numpy','scipy','pandas','openpyxl','matplotlib'])
    manifest['runtime_s']=time.perf_counter()-start
    manifest['outputs']=[dict(path=str(p.relative_to(out)),bytes=p.stat().st_size,
                              sha256=hashlib.sha256(p.read_bytes()).hexdigest())
                         for folder in ('results','figures') if (out/folder).exists()
                         for p in sorted((out/folder).rglob('*'))
                         if p.is_file() and p.name != '复现清单.json']
    (out/'results').mkdir(exist_ok=True)
    (out/'results/复现清单.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(dict(status='PASS',output_directory=str(out),runtime_s=manifest['runtime_s']),ensure_ascii=False),flush=True)


if __name__=='__main__':
    main()
