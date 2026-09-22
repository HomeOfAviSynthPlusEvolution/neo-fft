"""Original AVS binary versus Neo-FFT VS, without reading algorithm sources."""
import argparse,json,sys,hashlib
from pathlib import Path
import numpy as np
import vapoursynth as vs
from avs_reference import AviSynth
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'vapoursynth'))
from fixtures import environment

def main():
    ap=argparse.ArgumentParser()
    for name in ('plugin','reference','avisynth','fftw','report'):ap.add_argument('--'+name,type=Path,required=True)
    ap.add_argument('--seed',type=int,default=1);ap.add_argument('--opt',type=int,default=1);ap.add_argument('--quick',action='store_true');args=ap.parse_args()
    for path in (args.plugin,args.reference,args.avisynth,args.fftw):
        if not path.is_file():raise FileNotFoundError(path)
    policy=environment(vs);c=vs.core;c.num_threads=1;c.std.LoadPlugin(path=str(args.plugin.resolve()))
    avs=AviSynth(args.avisynth.resolve(),args.fftw.resolve());results=[];uncovered=[]
    # Fixed before holdouts; normalized float and integer output limits are separate.
    budgets={8:1,16:4,32:2e-5}
    try:
        avs.evaluate('SetMemoryMax(64)\nLoadPlugin("'+args.reference.resolve().as_posix()+'")\nreturn BlankClip(length=1)')
        formats=[('gray',8),('420',8),('rgb',8),('gray',16),('420',16),('rgb',16),('gray',32),('420',32),('rgb',32)]
        for family,bits in formats:
            convert={'gray':'.ConvertToY8()','420':'.ConvertToYUV420()','rgb':'.ConvertToPlanarRGB()'}[family]
            pieces=['base.Trim(0,0).Levels(0,1,255,'+str((n*n*3+args.seed*7)%32)+','+str(255-(n*7+args.seed*11)%48)+',coring=false)' for n in range(9)]
            script='base=ColorBars(width=32,height=24,pixel_type="YV24")\nglobal p5src=('+ '+'.join(pieces)+')'+convert+('.ConvertBits('+str(bits)+')' if bits!=8 else '')+'\nreturn p5src'
            source=avs.evaluate(script);samples=[avs.frame(source,n,bits,family) for n in range(9)]
            color=vs.GRAY if family=='gray' else (vs.RGB if family=='rgb' else vs.YUV)
            fmt=c.query_video_format(color,vs.FLOAT if bits==32 else vs.INTEGER,bits,1 if family=='420' else 0,1 if family=='420' else 0)
            base=c.std.BlankClip(width=32,height=24,format=fmt.id,length=9)
            def fill(n,f):
                out=f.copy()
                for p in range(out.format.num_planes):np.asarray(out[p])[:]=samples[n][p]
                return out
            clip=c.std.ModifyFrame(base,clips=base,selector=fill)
            variants=[]
            for T,O in ((2,0),(2,1),(3,1),(4,2),(5,2),(6,4)):
                for mode in (0,1):
                    for ft in range(5):variants.append((T,O,mode,ft,7,True,"scalar"))
            for window in range(12):
                for T,O in ((4,0),(5,2),(6,4)):variants.append((T,O,1,0,window,False,"scalar"))
            for model in ('curve','axis','sample'):
                for T,O in ((4,0),(4,2),(5,2),(6,4)):
                    variants.append((T,O,1,0,7,True,model))
            if args.quick:variants=variants[:2]
            for T,O,mode,ft,window,mean,model in variants:
                kw=dict(tbsize=T,tmode=1,tosize=O,sbsize=5 if mode==0 else 8,smode=mode,sosize=4,swin=0,twin=window,sigma=8 if ft<2 else .7,sigma2=.3,pmin=0,pmax=500,ftype=ft,zmean=mean,dither=0,threads=1)
                if model=='curve':kw['slocation']=[0,4,1,12]
                if model=='axis':kw['sst']=[0,4,1,12]
                if model=='sample':kw.update(nlocation=[1,0,2,2,1,0,2,2],alpha=4)
                refkw=kw.copy()
                for key,old in (('slocation','sstring'),('sst','sst')):
                    if key in refkw:
                        values=refkw.pop(key);refkw[old]=' '.join(format(float(np.float32(values[i])),'.9g')+':'+format(float(np.float32(values[i+1])),'.9g') for i in range(0,len(values),2))
                if 'nlocation' in refkw:
                    values=refkw.pop('nlocation');alpha=refkw.pop('alpha')
                    refkw['nstring']='a:'+str(alpha)+' '+' '.join(','.join(str(x) for x in values[i:i+4]) for i in range(0,len(values),4))
                def val(v):return ('"'+v+'"') if isinstance(v,str) else (str(v).lower() if isinstance(v,bool) else str(v))
                text=','.join(k+'='+val(v) for k,v in refkw.items())
                try:ref=avs.evaluate('return p5src.DFTTest('+text+',opt=1,lsb=false,lsb_in=false)')
                except RuntimeError as error:
                    known_errors={
                        'Filter error: GetPlaneHeightSubsampling called with unsupported plane.',
                        'Filter error: GetPlaneHeightSubsampling not available on greyscale pixel type.',
                    }
                    if model!='sample' or str(error) not in known_errors:raise
                    uncovered.append(dict(family=family,bits=bits,parameters=kw,reason=str(error)));continue
                own=c.neo_fft.DFTTest(clip,opt=args.opt,**kw);expected=[];worst=0
                for n in range(9):
                    original=avs.frame(ref,n,bits,family);expected.append(original);frame=own.get_frame(n)
                    for p,b in enumerate(original):
                        delta=np.abs(np.asarray(frame[p]).astype(float)-b.astype(float))
                        assert np.isfinite(delta).all(),('nonfinite',family,bits,kw,n,p)
                        worst=max(worst,float(np.max(delta)))
                assert worst<=budgets[bits],(family,bits,kw,worst)
                avs.drop(ref)
                # Fresh reference instance: characterize ring reuse independently.
                ref=avs.evaluate('return p5src.DFTTest('+text+',opt=1,lsb=false,lsb_in=false)')
                for n in (8,1,5,0,7,2):
                    actual=avs.frame(ref,n,bits,family)
                    for p,b in enumerate(actual):assert np.array_equal(b,expected[n][p]),('reference request-order',family,bits,kw,n)
                avs.drop(ref);results.append(dict(family=family,bits=bits,parameters=kw,max_error=worst))
            avs.drop(source);print('reference',family,bits,'passed',sum(r['family']==family and r['bits']==bits for r in results),'uncovered',sum(r['family']==family and r['bits']==bits for r in uncovered),flush=True)
    finally:
        avs.close()
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(json.dumps(dict(seed=args.seed,budgets=budgets,reference_sha256=hashlib.sha256(args.reference.read_bytes()).hexdigest(),fftw_sha256=hashlib.sha256(args.fftw.read_bytes()).hexdigest(),results=results,uncovered=uncovered),indent=2))
    print('Original AVS differential cases:',len(results),'uncovered sampled combinations:',len(uncovered))
if __name__=='__main__':main()
