import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import vapoursynth as vs
from fixtures import environment
F=np.float32

def mix(v):
    v^=v>>16; v=(v*0x7feb352d)&0xffffffff; v^=v>>15; v=(v*0x846ca68b)&0xffffffff
    return v^(v>>16)
def hash32(seed,n,p,y,x):
    h=mix(seed^0xa511e9b3)
    for v in (n,p,y,x): h=mix(h^v)
    return h
def oracle(a,mode,seed,n,p):
    h,w=a.shape; out=np.empty((h,w),np.uint8); current=np.zeros(w,np.float32)
    scale=F(F(mode-1)+F(.5)); off=F(scale*F(.5))
    for y in range(h):
        nxt=np.zeros(w,np.float32)
        for x in range(w):
            e=a[y,x]
            if mode==1: v=F(F(e+current[x])+F(.5))
            else:
                u=F(F(hash32(seed,n,p,y,x)>>8)*F(2**-24))
                v=F(F(F(F(e+F(u*scale))-off)+current[x])+F(.5))
            d=int(np.clip(v,0,255)); out[y,x]=d; err=F(e-F(d))
            if x: nxt[x-1]=F(nxt[x-1]+F(err*F(.1875)))
            nxt[x]=F(nxt[x]+F(err*F(.3125)))
            if x+1<w:
                current[x+1]=F(current[x+1]+F(err*F(.4375)))
                nxt[x+1]=F(nxt[x+1]+F(err*F(.0625)))
        current=nxt
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plugin',type=Path,required=True);ap.add_argument('--opt',type=int,default=1);args=ap.parse_args()
    policy=environment(vs);c=vs.core;c.num_threads=4;c.std.LoadPlugin(path=str(args.plugin.resolve()))
    for coords,h in (((0,0,0,0,0),0x642a1d14),((0,1,0,0,0),0xd32e9dca),((1,0,0,0,0),0xf383f4c5),((17,23,2,5,7),0x9eb639f9),((2147483647,2147483646,2,107,255),0x5064b117)):
        assert hash32(*coords)==h
    count=0
    for w in (1,2,7,17,33):
        for fmt in (vs.GRAY8,vs.YUV444P8):
            src=c.std.BlankClip(width=w,height=7,format=fmt,length=5)
            def fill(n,f):
                out=f.copy()
                for p in range(out.format.num_planes):
                    y,x=np.indices((7,w)); np.asarray(out[p])[:]=(x*31+y*17+n*23+p*47)%256
                return out
            src=c.std.ModifyFrame(src,clips=src,selector=fill)
            for mode in (1,2,4,2147483647):
                for seed in (0,17,2147483647):
                    kw=dict(tbsize=1,sbsize=1,smode=0,swin=7,twin=7,zmean=False,ftype=2,sigma=.375,opt=args.opt,dither=mode,dither_seed=seed)
                    result=c.neo_fft.DFTTest(src,**kw)
                    def check(n):
                        actual=result.get_frame(n); original=src.get_frame(n)
                        for p in range(actual.format.num_planes):
                            e=F(np.asarray(original[p]))*F(.375)
                            assert np.array_equal(np.asarray(actual[p]),oracle(e,mode,seed,n,p)),(w,mode,seed,n,p)
                    for n in (4,1,3,0,2,4): check(n)
                    result=c.neo_fft.DFTTest(src,**kw)
                    with ThreadPoolExecutor(4) as pool: list(pool.map(check,(3,1,4,3,2)))
                    count+=1
    for fmt in (vs.GRAY16,vs.GRAYS):
        src=c.std.BlankClip(width=17,height=9,format=fmt,color=[.25 if fmt==vs.GRAYS else 10000],length=1)
        kw=dict(tbsize=1,sbsize=3,smode=0,opt=args.opt)
        a=c.neo_fft.DFTTest(src,dither=2147483647,dither_seed=17,**kw).get_frame(0)
        b=c.neo_fft.DFTTest(src,dither=0,**kw).get_frame(0)
        assert np.asarray(a[0]).tobytes()==np.asarray(b[0]).tobytes()
    print('Phase 4 dither independent binary32 oracle:',count,'configurations, shuffled/concurrent requests and ignored formats passed')
if __name__=='__main__':main()
