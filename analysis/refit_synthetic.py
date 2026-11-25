"""Optional one-case refitting; never invoked by plot/verification commands."""
from pathlib import Path
import argparse,json,subprocess,sys

ROOT=Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case',required=True)
    ap.add_argument('--method',choices=['MCMC','free_LA'],required=True)
    args=ap.parse_args()
    manifest=json.loads((ROOT/'data/synthetic/benchmark_manifest.json').read_text())
    cases={c['case_id']:c for c in manifest['cases']}
    if args.case not in cases:ap.error('Unknown case; see data/synthetic/benchmark_manifest.json')
    dest=ROOT/'results/synthetic_refits';dest.mkdir(parents=True,exist_ok=True)
    c=cases[args.case]
    if args.method=='MCMC':
        cmd=[sys.executable,str(ROOT/'ce3d_xps_fitter_core.py'),'--input',str(ROOT/c['input']),
             '--outdir',str(dest/'MCMC'/args.case),'--background-mode','provided','--model','LA',
             '--walkers','80','--steps','8000','--burn','4000','--seed',str(c['seed'])]
        subprocess.run(cmd,check=True,cwd=ROOT)
    else:
        import free_LA_fit as model
        for r in manifest['cases']:r['input']=str(ROOT/r['input'])
        (dest/'benchmark_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
        model.WORK=dest;model.DATA=dest
        model.fit_case(args.case)

if __name__=='__main__':main()
