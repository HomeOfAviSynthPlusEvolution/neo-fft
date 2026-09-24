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
    assert c.neo_fft.identifier == 'org.neofilters.neo_fft'
    src=c.std.BlankClip(width=128,height=96,format=vs.GRAY8,color=[128])
    # A finite spectrum can overflow inside inverse FFT before normalization.
    # Every backend must either produce finite values or reject the request.
    impulse=c.std.BlankClip(width=64,height=64,format=vs.GRAYS,length=1)
    def make_impulse(n,f):
        out=f.copy();np.asarray(out[0])[0,0]=1;return out
    impulse=c.std.ModifyFrame(impulse,clips=impulse,selector=make_impulse)
    for opt in (0,1):
        overflow=c.neo_fft.DFTTest(impulse,tbsize=1,sbsize=16,sosize=0,swin=7,twin=7,zmean=False,ftype=2,sigma=1e36,opt=opt)
        try:
            frame=overflow.get_frame(0)
        except vs.Error as e:
            assert 'finite' in str(e)
        else:
            assert np.all(np.isfinite(np.asarray(frame[0])))
    # Phase-2 temporal defaults bt=3, tbsize=3 succeed on clips with N >= 3
    c.neo_fft.FFT3D(src).get_frame(0)
    c.neo_fft.DFTTest(src).get_frame(0)

    # Invalid temporal parameters rejection
    c.neo_fft.FFT3D(src,bt=-1).get_frame(0)
    c.neo_fft.FFT3D(src,bt=0).get_frame(0)
    fails(lambda:c.neo_fft.FFT3D(src,bt=6),'outside -1..5')
    fails(lambda:c.neo_fft.DFTTest(src,tbsize=0))
    fails(lambda:c.neo_fft.DFTTest(src,tbsize=2))
    fails(lambda:c.neo_fft.DFTTest(src,tbsize=4))
    fails(lambda:c.neo_fft.DFTTest(src,tbsize=17))
    fails(lambda:c.neo_fft.DFTTest(src,tbsize=3,tmode=2))

    # Short clip rejection for DFTTest vs dynamic fallback for FFT3D
    short_clip = c.std.BlankClip(width=128, height=96, format=vs.GRAY8, length=2)
    fails(lambda:c.neo_fft.DFTTest(short_clip,tbsize=3),'less than or equal to the number of frames')
    c.neo_fft.FFT3D(short_clip,bt=5).get_frame(0) # Degrades to 2D without error

    for name,temporal in [('FFT3D','bt'),('DFTTest','tbsize')]:
        call=getattr(c.neo_fft,name)
        for kwargs in [dict(planes=[-1]),dict(planes=[1]),
                       dict(sigma=float('nan')),dict(sigma=float('inf')),dict(sigma=1e100),dict(opt=1<<40),dict(sigma=-1)]:
            fails(lambda:call(src,**{temporal:1},**kwargs))
        fails(lambda:call(src,**{temporal:1},sigma='wrong'))
        fails(lambda:call(src,**{temporal:1},unknown=0))
    for kwargs in [dict(sharpen=-.1),dict(sigma2=-1),dict(l=-1),dict(kratio=-1),
                   dict(wintype=3),dict(beta=0),dict(ow=17),dict(bw=1)]:
        fails(lambda:c.neo_fft.FFT3D(src,bt=1,**kwargs))
    for kwargs in [dict(dither=-1),dict(nlocation=[0]),dict(ssx=[1.]),
                   dict(ftype=5),dict(f0beta=0),dict(pmin=2,pmax=1),dict(sbsize=8,sosize=5),dict(swin=12),
                   dict(smode=0,sbsize=4),dict(dither_seed=-1),dict(alpha=0)]:
        fails(lambda:c.neo_fft.DFTTest(src,tbsize=1,planes=[],**kwargs))
    # Legacy controls are registered but never read by the plugin/DS2 parser.
    for name,ignored in [('FFT3D',dict(ncpu=(-(1<<40),0,1<<40),measure=(False,True),mt=(False,True))),
                         ('DFTTest',dict(fft_threads=(-(1<<40),0,1<<40),threads=(-(1<<40),0,1<<40)))]:
        call=getattr(c.neo_fft,name)
        ignored['fft_backend']=('pocketfft','fftw','gpu',b'unknown\x00backend')
        for parameter,values in ignored.items():
            for value in values:
                call(src,**{parameter:value}).get_frame(0)
    c.neo_fft.DFTTest(src,tbsize=1,smode=0,sbsize=3,sosize=-999,tosize=-999,threads=-1).get_frame(0)
    small=c.std.BlankClip(width=1,height=1,format=vs.GRAY8,color=[12])
    c.neo_fft.DFTTest(small,tbsize=1,smode=0,sbsize=1).get_frame(0)
    fails(lambda:c.neo_fft.DFTTest(small,tbsize=1,smode=1,sbsize=1,sosize=0),'reflection')
    fails(lambda:c.neo_fft.FFT3D(c.std.BlankClip(width=8,height=8,format=vs.GRAY8),bt=1,bw=8,bh=8,ow=0,oh=0),'reflection')
    fails(lambda:c.neo_fft.FFT3D(c.std.BlankClip(width=12,height=12,format=vs.GRAY8),bt=1,bw=8,bh=8,ow=0,oh=0),'reflection')
    fails(lambda:c.neo_fft.DFTTest(src,tbsize=1,sbsize=4,sosize=0,swin=6,zmean=True),'DC')
    c.neo_fft.DFTTest(src,tbsize=1,sbsize=4,sosize=0,swin=6,zmean=False).get_frame(0)

    # Temporal execution at boundaries and center
    for bt in (2, 3, 4, 5):
        c.neo_fft.FFT3D(src, bt=bt).get_frame(0)
        c.neo_fft.FFT3D(src, bt=bt).get_frame(1)
        c.neo_fft.FFT3D(src, bt=bt).get_frame(src.num_frames - 1)
    for tbsize in (3, 5):
        c.neo_fft.DFTTest(src, tbsize=tbsize).get_frame(0)
        c.neo_fft.DFTTest(src, tbsize=tbsize).get_frame(1)
        c.neo_fft.DFTTest(src, tbsize=tbsize).get_frame(src.num_frames - 1)

    # Empty planes on temporal filter must strictly isolate and only request target frame
    def only_frame_4(n, f):
        if n != 4:
            raise vs.Error("neighbor frame accessed")
        return f
    guarded_clip = c.std.ModifyFrame(src, clips=src, selector=only_frame_4)
    c.neo_fft.DFTTest(guarded_clip, tbsize=5, planes=[]).get_frame(4)
    fails(lambda: c.neo_fft.FFT3D(src, bt=5, bw=8, bh=8, sigma=1.1e18), 'finite')

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
            for value in (float('inf'), -float('inf')):
                def infinite(n, f):
                    out=f.copy(); np.asarray(out[0])[:]=value; return out
                infclip=c.std.ModifyFrame(floatclip,clips=floatclip,selector=infinite)
                fails(lambda:call(infclip,planes=[0],**kwargs).get_frame(0),'non-finite')
    info=c.neo_fft.KernelInfo()
    assert 'fft_threads' not in info
    assert info['fft_backend']
    print('VS admission, negative inputs, selected NaN errors, bitwise copies, properties passed:',info)

if __name__=='__main__':main()
