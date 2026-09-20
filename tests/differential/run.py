"""Public VS capture and per-plane comparison. Calibration never constitutes acceptance."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tests/vapoursynth'))
from fixtures import cases

def metrics(a,b,mask):
    yy,xx=np.nonzero(mask)
    if not len(yy): return None
    av,bv=a[mask].astype(np.float64),b[mask].astype(np.float64)
    error=np.abs(av-bv)
    idx=int(error.argmax())
    return dict(max=float(error[idx]),fraction=float(np.count_nonzero(error)/len(error)),
        mae=float(error.mean()),rmse=float(np.sqrt(np.mean(error*error))),
        location=[int(yy[idx]),int(xx[idx])],reference=float(av[idx]),candidate=float(bv[idx]))

def compare(root,catalog,budgets):
    rows=[]; failed=[]; maxima={}; exceptions=[]
    records={}
    for algorithm in ('FFT3D','DFTTest'):
        for label in ('old1','new1','old0','new0'):
            manifest=json.loads((root/algorithm/label/'manifest.json').read_text())
            if label=='new0' and manifest['kernel_info']['target'] in ('scalar','SCALAR','EMU128','scalar (Highway disabled)'):
                raise RuntimeError('required SIMD capture used a scalar target')
            records[algorithm,label]={r['case']['id']:r for r in manifest['cases']}
    for case in catalog:
        a=case['algorithm']; cid=case['id']; bits=case['bits']
        edge=max(case['params'].get('bw',32),case['params'].get('bh',32)) if a=='FFT3D' else case['params'].get('sbsize',16)
        for label,left,right in (('A','old1','new1'),('B','old1','new0'),('C','old0','new0'),('D','new1','new0')):
            lr=records[a,left].get(cid); rr=records[a,right].get(cid)
            if not lr or not rr or lr['status'] not in ('ok','interface_exception') or rr['status']!='ok':
                failed.append(dict(case=cid,pair=label,error='capture failed or manifest entry missing'))
                continue
            if lr['input_sha256']!=rr['input_sha256'] or lr['case']!=rr['case']:
                raise RuntimeError('capture inputs or case parameters differ')
            if a=='DFTTest' and case['params'].get('planes')==[] and label!='D':
                record=lr
                if record['status']!='interface_exception': raise RuntimeError('unverified reference interface exception')
                exceptions.append(dict(case=cid,pair=label,reason=record['interface_exception']))
                continue
            lp,rp=root/a/left/(cid+'.npz'),root/a/right/(cid+'.npz')
            if not lp.exists() or not rp.exists(): failed.append(dict(case=cid,pair=label,error='missing output')); continue
            with np.load(lp) as lhs,np.load(rp) as rhs:
                if set(lhs.files)!=set(rhs.files): raise RuntimeError('captured frame/plane sets differ')
                for plane in lhs.files:
                    av,bv=lhs[plane],rhs[plane]
                    if av.shape!=bv.shape or av.dtype!=bv.dtype: raise RuntimeError('inconsistent captured format')
                    h,w=av.shape; yy,xx=np.indices((h,w))
                    ey=(yy<edge)|(yy>=h-edge); ex=(xx<edge)|(xx>=w-edge)
                    regions={'full':np.ones((h,w),bool),'edge':ey|ex,'interior':~(ey|ex),'corner':ey&ex}
                    result={name:metrics(av,bv,mask) for name,mask in regions.items()}
                    key=f'{a}/{bits}/pocketfft/{label}'
                    maxima[key]=max(maxima.get(key,0),result['full']['max'])
                    row=dict(case=cid,plane=plane,pair=label,budget_key=key,regions=result)
                    if budgets is not None:
                        budget=budgets['limits'].get(key)
                        if budget is None: raise RuntimeError('missing frozen budget '+key)
                        error=np.abs(av.astype(np.float64)-bv.astype(np.float64))
                        limit=budget['atol']+budget.get('rtol',0)*np.maximum(np.abs(av.astype(float)),np.abs(bv.astype(float)))
                        row['pass']=bool(np.all(error<=limit))
                        if not row['pass']: failed.append(row)
                    rows.append(row)
    (root/'metrics.json').write_text(json.dumps(rows,indent=2))
    summary=dict(comparisons=len(rows),failures=failed,maxima=maxima,interface_exceptions=exceptions,budget_sha256=None)
    if budgets is not None: summary['budget_sha256']=hashlib.sha256(json.dumps(budgets,sort_keys=True).encode()).hexdigest()
    (root/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(dict(comparisons=len(rows),failures=len(failed),maxima=maxima),indent=2))
    return bool(failed)

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--plugin',type=Path,required=True)
    p.add_argument('--fft3d-reference',type=Path,required=True);p.add_argument('--dfttest-reference',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--seed',type=int,required=True)
    p.add_argument('--budgets',type=Path);p.add_argument('--compare-only',action='store_true')
    args=p.parse_args()
    if args.seed!=1 and not args.budgets: p.error('holdout requires a frozen budget file')
    budget=json.loads(args.budgets.read_text()) if args.budgets else None
    catalog=cases();args.output.mkdir(parents=True,exist_ok=True)
    casefile=args.output/'cases.json'
    if not args.compare_only:
        casefile.write_text(json.dumps(catalog,indent=2))
        for a,ref in [('FFT3D',args.fft3d_reference),('DFTTest',args.dfttest_reference)]:
            for label,opt,plugin,new in [('old1',1,ref,False),('new1',1,args.plugin,True),('old0',0,ref,False),('new0',0,args.plugin,True)]:
                cmd=[sys.executable,str(ROOT/'tests/vapoursynth/worker.py'),'--plugin',str(plugin),'--algorithm',a,
                     '--cases',str(casefile),'--output',str(args.output/a/label),'--seed',str(args.seed),'--opt',str(opt)]
                if new: cmd.append('--new')
                result=subprocess.run(cmd,text=True,capture_output=True)
                dest=args.output/a/label;dest.mkdir(parents=True,exist_ok=True)
                (dest/'process.log').write_text(result.stdout+result.stderr)
                print(result.stdout,flush=True)
                if result.returncode:
                    print(result.stderr,flush=True)
                    raise RuntimeError(f'capture {a}/{label} failed; see {dest / "process.log"}')
    return compare(args.output,catalog,budget)

if __name__=='__main__': raise SystemExit(main())
