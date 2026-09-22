import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import vapoursynth as vs
from fixtures import environment,source

def starts(n,T,O):
    return [s for s in range(-max(T-O,O),n+1,T-O) if s+T>n]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plugin',type=Path,required=True);ap.add_argument('--opt',type=int,default=1);args=ap.parse_args()
    policy=environment(vs);c=vs.core;c.num_threads=4;c.std.LoadPlugin(path=str(args.plugin.resolve()))
    src=c.std.BlankClip(width=5,height=3,format=vs.GRAYS,length=31)
    def fill(n,f):
        out=f.copy();y,x=np.indices((3,5));np.asarray(out[0])[:]=((n*n*7+y*13+x*19)%113)/128.;out.props['Index']=n;return out
    src=c.std.ModifyFrame(src,clips=src,selector=fill)
    samples=[np.asarray(src.get_frame(n)[0]).copy() for n in range(31)]
    cases=0
    for T in range(1,16):
        for O in range(T):
            if O>T//2 and T%(T-O):continue
            h=np.array([1/np.sqrt(len(range(z%(T-O),T,T-O)))/np.sqrt(T) for z in range(T)],np.float32)
            for mean in (False,True):
                kw=dict(tmode=1,tbsize=T,tosize=O,smode=0,sbsize=1,swin=7,twin=7,ftype=2,sigma=0 if mean else .375,zmean=mean,opt=args.opt)
                out=c.neo_fft.DFTTest(src,**kw)
                for n in sorted(set((0,1,T-1,T,29,30))):
                    expected=np.zeros((3,5),np.float64)
                    for s in starts(n,T,O):
                        value=sum(samples[max(0,min(30,s+z))]*float(h[z]) for z in range(T))/float(sum(h)) if mean else samples[n]*.375
                        expected+=value*(T*float(h[n-s])**2)
                    frame=out.get_frame(n)
                    assert np.max(np.abs(np.asarray(frame[0])-expected))<3e-6,(T,O,mean,n)
                    assert frame.props['Index']==n
                cases+=1
    # Independent sampled-noise calibration, including even T and non-lattice fn.
    # B=1 leaves a full complex temporal axis; use numpy's full DFT, not rfft.
    for T,O in ((2,0),(4,0),(4,2),(5,2),(6,4)):
        h=np.array([1/np.sqrt(len(range(z%(T-O),T,T-O)))/np.sqrt(T) for z in range(T)],np.float32)
        h2=np.full(T,1/np.sqrt(T),np.float32)
        E=np.float32(0);E2=np.float32(0)
        for z in range(T):E=np.float32(E+np.float32(h[z]*h[z]));E2=np.float32(E2+np.float32(h2[z]*h2[z]))
        sample=np.array([samples[1+z][0,0]*255*h2[z] for z in range(T)])
        power=np.abs(np.fft.fft(sample))**2*(float(E)/float(E2))*4
        out=c.neo_fft.DFTTest(src,tmode=1,tbsize=T,tosize=O,smode=0,sbsize=1,swin=7,twin=7,zmean=False,ftype=0,nlocation=[1,0,0,0,1,0,0,0],alpha=4,opt=args.opt)
        for n in (0,2,3,30):
            expected=np.zeros((3,5),np.float64)
            for s in starts(n,T,O):
                block=np.stack([samples[max(0,min(30,s+z))]*255*float(h[z]) for z in range(T)])
                spectrum=np.fft.fft(block,axis=0);psd=np.abs(spectrum)**2
                spectrum*=np.maximum((psd-power[:,None,None])/(psd+1e-15),0)
                inverse=np.fft.ifft(spectrum,axis=0).real
                expected+=inverse[n-s]*T*float(h[n-s])/255
            actual=np.asarray(out.get_frame(n)[0]);assert np.max(np.abs(actual-expected))<1e-5,(T,O,n,'sample oracle')
    # Matrix with real spatial overlap, models and final diffusion. Fresh requests
    # ensure neither model readiness nor frame cache can hide ordering differences.
    for family,bits in (('gray',8),('420',16),('rgb',32)):
        clip,_,_=source(vs,dict(format=family,bits=bits,width=40,height=32,frames=9),17)
        for T,O in ((2,0),(3,1),(4,2),(5,2),(6,4)):
            for mode in (0,1):
                for ftype in range(5):
                    kw=dict(tmode=1,tbsize=T,tosize=O,sbsize=5 if mode==0 else 4,smode=mode,sosize=2,ftype=ftype,sigma=8 if ftype<2 else .7,zmean=True,opt=args.opt,dither=4,dither_seed=17)
                    if ftype==0:kw['sst']=[0,5,1,9]
                    if ftype==1:kw['nlocation']=[1,0,0,0,1,0,0,0]
                    serial=c.neo_fft.DFTTest(clip,threads=1,**kw)
                    expected={n:[np.asarray(serial.get_frame(n)[p]).tobytes() for p in range(clip.format.num_planes)] for n in range(9)}
                    parallel=c.neo_fft.DFTTest(clip,threads=3,**kw)
                    for n in (8,1,6,0,5,2,8):assert [np.asarray(parallel.get_frame(n)[p]).tobytes() for p in range(clip.format.num_planes)]==expected[n]
                    concurrent=c.neo_fft.DFTTest(clip,threads=1,**kw)
                    def check(n):
                        f=concurrent.get_frame(n);assert [np.asarray(f[p]).tobytes() for p in range(f.format.num_planes)]==expected[n]
                    with ThreadPoolExecutor(4) as pool:list(pool.map(check,(8,1,8,4,0,6)))
                    cases+=1
    clip,original,_=source(vs,dict(format='420',bits=8,width=40,height=32,frames=9),101)
    for planes in ([0],[0,0],[]):
        result=c.neo_fft.DFTTest(clip,tmode=1,tbsize=4,tosize=2,sbsize=4,sosize=2,planes=planes,opt=args.opt)
        for n in (8,1,5):
            f=result.get_frame(n)
            for p in (range(3) if not planes else (1,2)):assert np.asarray(f[p]).tobytes()==original[n][p].tobytes()
    for T in range(1,16):
        for O in (-1,T,T+1):
            try:c.neo_fft.DFTTest(src,tmode=1,tbsize=T,tosize=O,planes=[])
            except vs.Error:pass
            else:raise AssertionError((T,O))
    # No-window pass-through still checks raw temporal controls, and requests n only.
    def guarded(n,f):
        if n!=8:raise vs.Error('unexpected source')
        return f
    safe=c.std.ModifyFrame(src,clips=src,selector=guarded)
    c.neo_fft.DFTTest(safe,tmode=1,tbsize=6,tosize=4,planes=[]).get_frame(8)
    print('Phase 5 temporal mean/identity and deterministic integration configurations:',cases)
if __name__=='__main__':main()
