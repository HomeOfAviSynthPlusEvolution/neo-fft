import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import vapoursynth as vs
from fixtures import environment,source

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plugin',type=Path,required=True);ap.add_argument('--reference',type=Path);ap.add_argument('--opt',type=int,default=1);args=ap.parse_args()
    policy=environment(vs);c=vs.core;c.num_threads=4;c.std.LoadPlugin(path=str(args.plugin.resolve()))
    if args.reference: c.std.LoadPlugin(path=str(args.reference.resolve()))
    count=0
    for family in ('gray','420','rgb'):
        for bits in (8,16,32):
            src,original,_=source(vs,dict(format=family,bits=bits,width=64,height=48,frames=9),37)
            for model in ({},{'sigma2':9},{'pfactor':.7,'pframe':8,'px':2,'py':2}):
                for fields in (False,True):
                    kw=dict(bt=0,bw=8,bh=8,ow=4,oh=4,sigma=12,kratio=2,opt=args.opt,l=4,t=4,r=4,b=4,interlaced=fields,**model)
                    baseline=c.neo_fft.FFT3D(src,**kw)
                    expected=[]
                    for n in range(9):
                        f=baseline.get_frame(n);expected.append([np.asarray(f[p]).tobytes() for p in range(f.format.num_planes)])
                    assert expected[0]==[a.tobytes() for a in original[0]]
                    result=c.neo_fft.FFT3D(src,**kw)
                    def check(n):
                        f=result.get_frame(n)
                        assert [np.asarray(f[p]).tobytes() for p in range(f.format.num_planes)]==expected[n],(family,bits,model,fields,n)
                        assert dict(f.props)==dict(src.get_frame(n).props)
                    for n in (8,2,7,1,5,3,8):check(n)
                    result=c.neo_fft.FFT3D(src,**kw)
                    with ThreadPoolExecutor(4) as pool:list(pool.map(check,(8,2,8,6,1,4)))
                    count+=1
    # Frame zero bypasses both nonfinite pixels and a forbidden pattern source.
    raw=c.std.BlankClip(width=64,height=48,format=vs.GRAYS,length=9)
    def zero(n,f):
        if n!=0: raise vs.Error('unexpected model/replay source')
        out=f.copy();np.asarray(out[0]).view(np.uint32)[:]=0x7fc12345;return out
    raw=c.std.ModifyFrame(raw,clips=raw,selector=zero)
    f=c.neo_fft.FFT3D(raw,bt=0,pfactor=1,pframe=8,px=2,py=2,bw=8,bh=8).get_frame(0)
    assert np.all(np.asarray(f[0]).view(np.uint32)==0x7fc12345)
    # Uniform zero is identity; pattern models retain sigma1 initial covariance.
    src,original,_=source(vs,dict(format='gray',bits=8,width=64,height=48,frames=9),19)
    identity=c.neo_fft.FFT3D(src,bt=0,sigma=0,bw=8,bh=8,opt=args.opt)
    for n in (8,1,0):assert np.asarray(identity.get_frame(n)[0]).tobytes()==original[n][0].tobytes()
    base=dict(bt=0,bw=8,bh=8,ow=4,oh=4,sigma=12,kratio=2,pframe=8,px=2,py=2,opt=args.opt)
    sampled=c.neo_fft.FFT3D(src,pfactor=1,**base)
    for strength in (.2,.7,2,1e38):
        other=c.neo_fft.FFT3D(src,pfactor=strength,**base)
        for n in (1,8):assert np.asarray(other.get_frame(n)[0]).tobytes()==np.asarray(sampled.get_frame(n)[0]).tobytes()
    # Shared 720p YUV420 state fits 64 MiB; 1080p still exceeds it. Both must
    # reuse the preceding checkpoint, while backward seeks stay canonical.
    for width,height in ((1280,720),(1920,1080)):
        hd,_,_=source(vs,dict(format='420',bits=8,width=width,height=height,frames=4),53)
        requested=[]
        def record(n,f):
            requested.append(n);return f
        hd=c.std.ModifyFrame(hd,clips=hd,selector=record)
        c.std.SetVideoCache(hd,mode=0)
        large=c.neo_fft.FFT3D(hd,bt=0,bw=32,bh=32,ow=16,oh=16,sigma=2,opt=args.opt,ncpu=1)
        c.std.SetVideoCache(large,mode=0)
        def pixels(n):
            frame=large.get_frame(n)
            return [np.asarray(frame[p]).tobytes() for p in range(frame.format.num_planes)]
        first=pixels(1);requested.clear();pixels(2)
        assert requested==[2],('checkpoint replayed history',width,height,requested)
        assert pixels(1)==first
    if args.reference:
        for model in ({},{'sigma2':9},{'pfactor':.7,'pframe':8,'px':2,'py':2}):
            kw=dict(bt=0,bw=8,bh=8,ow=4,oh=4,sigma=12,kratio=2,**model)
            own=c.neo_fft.FFT3D(src,opt=args.opt,**kw);ref=c.neo_fft3d.FFT3D(src,opt=0,**kw)
            for n in range(9):
                a=np.asarray(own.get_frame(n)[0]).copy();b=np.asarray(ref.get_frame(n)[0]).copy()
                d=np.abs(a.astype(int)-b.astype(int));print('reference',model,n,int(np.count_nonzero(d)),int(d.max()))
                assert d.max()<=1,(model,n,int(d.max()))
    print('Phase 4 Kalman:',count,'format/model/ROI/field configurations; canonical requests and frame-zero bypass passed')
if __name__=='__main__':main()
