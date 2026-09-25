"""Host contract: bounded cold requests, checkpoint continuation and warmup API."""
import argparse
from pathlib import Path
import numpy as np
import vapoursynth as vs
from fixtures import environment


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--plugin',type=Path,required=True)
    ap.add_argument('--opt',type=int,default=1)
    args=ap.parse_args()
    policy=environment(vs)
    c=vs.core;c.num_threads=2
    c.std.LoadPlugin(path=str(args.plugin.resolve()))
    for bits,fmt in ((8,vs.GRAY8),(16,vs.GRAY16),(32,vs.GRAYS)):
        blank=c.std.BlankClip(width=32,height=24,format=fmt,length=10005)
        yy,xx=np.indices((24,32))
        def make_source(start=0,virtual=False):
            trace=[]
            allowed=set()
            def fill(n,f):
                assert not (virtual and n==0), 'virtual initial frame was fetched'
                original=start+n-1 if virtual else n
                assert original in allowed, ('out-of-budget source',original,sorted(allowed))
                trace.append(original)
                value=(64+(original*7+xx*11+yy*3)%17)/255
                out=f.copy()
                array=value.astype(np.float32) if bits==32 else np.rint(value*((1<<bits)-1)).astype(np.uint8 if bits==8 else np.uint16)
                np.copyto(np.asarray(out[0]),array)
                out.props['OriginalFrame']=original
                return out
            node=c.std.ModifyFrame(blank,clips=blank,selector=fill)
            c.std.SetVideoCache(node,mode=0)
            return node,trace,allowed
        kw=dict(bt=0,bw=8,bh=8,ow=4,oh=4,sigma=12,opt=args.opt)
        def output(node,n):
            f=node.get_frame(n)
            assert np.isfinite(np.asarray(f[0])).all()
            return np.asarray(f[0]).tobytes()
        for warmup in (0,4,8,16):
            src,trace,allowed=make_source()
            node=c.neo_fft.FFT3D(src,kalman_warmup=warmup,**kw)
            c.std.SetVideoCache(node,mode=0)
            start=10000-warmup
            allowed.update(range(start,10001))
            first=output(node,10000)
            assert trace==list(range(start,10001))
            ref_src,ref_trace,ref_allowed=make_source(start,True)
            ref_allowed.update(range(start,10001))
            ref=c.neo_fft.FFT3D(ref_src,kalman_warmup=2147483647,**kw)
            assert first==output(ref,warmup+1)
            for n in (10001,10002):
                trace.clear();allowed.clear();allowed.add(n);ref_allowed.add(n)
                assert output(node,n)==output(ref,n-start+1)
                assert trace==[n]
            trace.clear();allowed.clear();allowed.add(10000)
            assert output(node,10000)==first and trace==[10000]
            # Backward request cannot use a future state and restarts locally.
            trace.clear();allowed.clear();allowed.update(range(9900-warmup,9901))
            backward=output(node,9900)
            back_src,_,back_allowed=make_source(9900-warmup,True)
            back_allowed.update(allowed)
            assert backward==output(c.neo_fft.FFT3D(back_src,kalman_warmup=2147483647,**kw),warmup+1)
            assert trace==list(range(9900-warmup,9901))
        src,trace,allowed=make_source();allowed.update(range(9992,10001))
        default=output(c.neo_fft.FFT3D(src,**kw),10000)
        assert trace==list(range(9992,10001))
        assert default==output(c.neo_fft.FFT3D(src,kalman_warmup=8,**kw),10000)
        # Sample model acquisition is separate from the recurrence budget.
        src,trace,allowed=make_source();allowed.update([11,*range(9992,10001)])
        sampled=c.neo_fft.FFT3D(src,pfactor=1,pframe=11,px=2,py=2,**kw)
        output(sampled,10000);assert trace==[11,*range(9992,10001)]
        # Frame zero bypasses both warmup and sampling; other modes ignore W.
        src,trace,allowed=make_source();allowed.add(0)
        output(c.neo_fft.FFT3D(src,kalman_warmup=16,pfactor=1,pframe=11,px=2,py=2,**kw),0)
        assert trace==[0]
        for value in (-1,2147483648):
            try:c.neo_fft.FFT3D(src,bt=1,planes=[],kalman_warmup=value)
            except vs.Error:pass
            else:raise AssertionError(('accepted invalid warmup',value))
        trace.clear();allowed.clear();allowed.add(10000)
        a=output(c.neo_fft.FFT3D(src,bt=1,bw=8,bh=8,kalman_warmup=0,opt=args.opt),10000)
        b=output(c.neo_fft.FFT3D(src,bt=1,bw=8,bh=8,kalman_warmup=2147483647,opt=args.opt),10000)
        assert a==b and trace==[10000,10000]
    print('Kalman warmup: u8/u16/float cold 10000, continuation, backward/exact checkpoints, sampled model and API validation passed')


if __name__=='__main__':main()
