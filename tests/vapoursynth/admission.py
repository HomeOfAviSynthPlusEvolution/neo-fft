import argparse
from pathlib import Path
import numpy as np
import vapoursynth as vs
from fixtures import environment

def fails(fn, fragment=None):
    try: fn()
    except (vs.Error, ValueError, TypeError) as e:
        if fragment and fragment.lower() not in str(e).lower(): raise AssertionError(str(e))
    else: raise AssertionError('expected VS error')

def main():
    p=argparse.ArgumentParser();p.add_argument('--plugin',type=Path,required=True);args=p.parse_args()
    policy=environment(vs); c=vs.core;c.num_threads=4
    c.std.LoadPlugin(path=str(args.plugin.resolve()))
    src=c.std.BlankClip(width=128,height=96,format=vs.GRAY8,color=[128])
    for name,temporal in [('FFT3D','bt'),('DFTTest','tbsize')]:
        call=getattr(c.neo_fft,name)
        fails(lambda:call(src),'unsupported')
        for kwargs in [dict(opt=2),dict(fft_backend='fftw'),dict(fft_backend='gpu'),dict(planes=[-1]),dict(planes=[1]),
                       dict(sigma=float('nan')),dict(sigma=float('inf')),dict(sigma=1e100),dict(opt=1<<40),dict(sigma=-1)]:
            fails(lambda:call(src,**{temporal:1},**kwargs))
        fails(lambda:call(src,**{temporal:1},sigma='wrong'))
        fails(lambda:call(src,**{temporal:1},unknown=0))
    for kwargs in [dict(sharpen=.1),dict(sigma2=1),dict(mt=True),dict(interlaced=True),dict(l=-1),dict(kratio=3),
                   dict(wintype=3),dict(beta=0),dict(ncpu=0),dict(ow=17),dict(bw=1)]:
        fails(lambda:c.neo_fft.FFT3D(src,bt=1,**kwargs))
    for kwargs in [dict(dither=1),dict(threads=2),dict(fft_threads=2),dict(nlocation=[0]),dict(ssx=[1.]),
                   dict(ftype=5),dict(f0beta=0),dict(pmin=2,pmax=1),dict(sbsize=8,sosize=5),dict(swin=12),
                   dict(smode=0,sbsize=4),dict(dither_seed=-1),dict(alpha=9)]:
        fails(lambda:c.neo_fft.DFTTest(src,tbsize=1,planes=[],**kwargs))
    c.neo_fft.DFTTest(src,tbsize=1,smode=0,sbsize=3,sosize=-999,tosize=-999,threads=-1,fft_threads=-1).get_frame(0)
    small=c.std.BlankClip(width=1,height=1,format=vs.GRAY8,color=[12])
    c.neo_fft.DFTTest(small,tbsize=1,smode=0,sbsize=1).get_frame(0)
    fails(lambda:c.neo_fft.DFTTest(small,tbsize=1,smode=1,sbsize=1,sosize=0),'reflection')
    fails(lambda:c.neo_fft.FFT3D(c.std.BlankClip(width=8,height=8,format=vs.GRAY8),bt=1,bw=8,bh=8,ow=0,oh=0),'reflection')
    fails(lambda:c.neo_fft.DFTTest(src,tbsize=1,sbsize=4,sosize=0,swin=6,zmean=True),'DC')
    c.neo_fft.DFTTest(src,tbsize=1,sbsize=4,sosize=0,swin=6,zmean=False).get_frame(0)
    floatclip=c.std.BlankClip(width=128,height=96,format=vs.YUV444PS)
    def unusual(n,f):
        out=f.copy()
        for p in range(3):
            arr=np.asarray(out[p]);arr[:]=.25
            if p==1 and n==0: arr.view(np.uint32)[:]=0x7fc12345
            if p==2: arr.view(np.uint32)[:]=0x80000000
        out.props['Keep']=b'payload\x00test'
        return out
    unusualclip=c.std.ModifyFrame(floatclip,clips=floatclip,selector=unusual)
    for name,temporal in [('FFT3D','bt'),('DFTTest','tbsize')]:
        call=getattr(c.neo_fft,name)
        for opt in (0,1):
            kwargs={temporal:1,'opt':opt}
            bad=call(unusualclip,planes=[1],**kwargs)
            fails(lambda:bad.get_frame(0),'non-finite')
            bad.get_frame(1)  # Recover using the same filter instance after a failed request.
            fails(lambda:call(unusualclip,**kwargs).get_frame(0),'non-finite')
            good=call(unusualclip,planes=[0],**kwargs).get_frame(0)
            original=unusualclip.get_frame(0)
            assert dict(good.props)==dict(original.props)
            for p in (1,2): assert np.asarray(good[p]).tobytes()==np.asarray(original[p]).tobytes()
            # Failure in one instance must not poison a subsequent valid request.
            call(src,**kwargs).get_frame(0)
    info=c.neo_fft.KernelInfo()
    print('VS admission, negative inputs, selected NaN errors, bitwise copies, properties passed:',info)

if __name__=='__main__':main()
