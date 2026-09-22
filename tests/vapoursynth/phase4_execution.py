import argparse
from pathlib import Path
import numpy as np
import vapoursynth as vs
from fixtures import environment,source

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plugin',type=Path,required=True);args=ap.parse_args()
    policy=environment(vs);c=vs.core;c.num_threads=1;c.std.LoadPlugin(path=str(args.plugin.resolve()))
    src,_,_=source(vs,dict(format='rgb',bits=8,width=128,height=96,frames=7),41)
    for bt in (-1,0,1,3,5):
        kw=dict(bt=bt,bw=8,bh=8,l=2,t=2,r=2,b=2,interlaced=True,opt=1)
        a=c.neo_fft.FFT3D(src,mt=False,ncpu=1,**kw);b=c.neo_fft.FFT3D(src,mt=True,ncpu=2147483647,**kw)
        for n in (6,0,2,4):
            fa=a.get_frame(n);fb=b.get_frame(n)
            for p in range(3):assert np.asarray(fa[p]).tobytes()==np.asarray(fb[p]).tobytes(),(bt,n,p)
    for T in (1,3):
        for mode in (0,1):
            for dither in (0,1,4):
                kw=dict(tbsize=T,sbsize=9 if mode==0 else 8,smode=mode,sosize=4,dither=dither,dither_seed=17,opt=1)
                a=c.neo_fft.DFTTest(src,threads=1,fft_threads=1,**kw)
                b=c.neo_fft.DFTTest(src,threads=3,fft_threads=2147483647,**kw)
                for n in (6,2,0):
                    fa=a.get_frame(n);fb=b.get_frame(n)
                    for p in range(3):assert np.asarray(fa[p]).tobytes()==np.asarray(fb[p]).tobytes(),(T,mode,dither,n,p)
    for name,control,aliases in (('FFT3D','bt',(-2147483648,-1,0,2,3,8,2147483647)),('DFTTest','tbsize',(0,2,3,8))):
        call=getattr(c.neo_fft,name);kw={control:1};expected=call(src,opt=0,**kw).get_frame(0)
        for opt in aliases:
            actual=call(src,opt=opt,**kw).get_frame(0)
            for p in range(3):assert np.asarray(actual[p]).tobytes()==np.asarray(expected[p]).tobytes()
    for opt in (-1,4,7,9,2147483647):
        try:c.neo_fft.DFTTest(src,opt=opt)
        except vs.Error:pass
        else:raise AssertionError(opt)
    print('Phase 4 execution: mt/threads bitwise invariant; opt aliases; requested FFT workers INT_MAX, effective PocketFFT workers 1')
if __name__=='__main__':main()
