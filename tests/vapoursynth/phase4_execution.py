import argparse
from pathlib import Path
import numpy as np
import vapoursynth as vs
from fixtures import environment,source

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plugin',type=Path,required=True);ap.add_argument('--opt',type=int,default=1);args=ap.parse_args()
    policy=environment(vs);c=vs.core;c.num_threads=4;c.std.LoadPlugin(path=str(args.plugin.resolve()))
    src,_,_=source(vs,dict(format='rgb',bits=8,width=128,height=96,frames=7),41)
    for bt in (-1,0,1,3,5):
        kw=dict(bt=bt,bw=8,bh=8,l=2,t=2,r=2,b=2,interlaced=True,opt=args.opt)
        a=c.neo_fft.FFT3D(src,**kw);b=c.neo_fft.FFT3D(src,**kw)
        order=(6,0,2,4)
        expected=[a.get_frame(n) for n in order]
        pending=[b.get_frame_async(n) for n in order]
        for n,fa,future in zip(order,expected,pending):
            fb=future.result()
            for p in range(3):assert np.asarray(fa[p]).tobytes()==np.asarray(fb[p]).tobytes(),(bt,n,p)
    for removed_mt in (False,True):
        try:c.neo_fft.FFT3D(src,mt=removed_mt)
        except vs.Error:pass
        else:raise AssertionError('removed mt parameter accepted')
    for T in (1,3):
        for mode in (0,1):
            for dither in (0,1,4):
                kw=dict(tbsize=T,sbsize=9 if mode==0 else 8,smode=mode,sosize=4,dither=dither,dither_seed=17,opt=args.opt)
                a=c.neo_fft.DFTTest(src,threads=1,**kw)
                b=c.neo_fft.DFTTest(src,threads=16,**kw)
                order=(6,2,0,6)
                expected=[a.get_frame(n) for n in order]
                pending=[b.get_frame_async(n) for n in order]
                for n,fa,future in zip(order,expected,pending):
                    fb=future.result()
                    for p in range(3):assert np.asarray(fa[p]).tobytes()==np.asarray(fb[p]).tobytes(),(T,mode,dither,n,p)
    kw=dict(tbsize=1,sbsize=8,sosize=4,opt=args.opt)
    expected=c.neo_fft.DFTTest(src,**kw).get_frame(0)
    for threads in (-2147483648,-1,0,1,3,16,17,2147483647):
        actual=c.neo_fft.DFTTest(src,threads=threads,**kw).get_frame(0)
        for p in range(3):assert np.asarray(actual[p]).tobytes()==np.asarray(expected[p]).tobytes(),threads
    for threads in (-2147483649,2147483648):
        try:c.neo_fft.DFTTest(src,threads=threads,**kw)
        except vs.Error:pass
        else:raise AssertionError(('out-of-int32 threads accepted',threads))
    for name,control,aliases in (('FFT3D','bt',(-2147483648,-1,0,2,3,8,2147483647)),('DFTTest','tbsize',(0,2,3,8))):
        call=getattr(c.neo_fft,name);kw={control:1};expected=call(src,opt=0,**kw).get_frame(0)
        for opt in aliases:
            actual=call(src,opt=opt,**kw).get_frame(0)
            for p in range(3):assert np.asarray(actual[p]).tobytes()==np.asarray(expected[p]).tobytes()
    for opt in (-1,4,7,9,2147483647):
        try:c.neo_fft.DFTTest(src,opt=opt)
        except vs.Error:pass
        else:raise AssertionError(opt)
    print('Phase 4 execution: host concurrency invariant; mt removed; DFTTest reserved threads domain and output invariant; opt aliases; effective PocketFFT workers 1')
if __name__=='__main__':main()
