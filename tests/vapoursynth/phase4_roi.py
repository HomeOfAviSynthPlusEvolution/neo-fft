import argparse
from pathlib import Path
import numpy as np
import vapoursynth as vs
from fixtures import environment, source

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--plugin',type=Path,required=True); args=ap.parse_args()
    policy=environment(vs); c=vs.core; c.num_threads=4
    c.std.LoadPlugin(path=str(args.plugin.resolve()))
    count=0
    for family in ('gray','420','rgb'):
        for bits in (8,16,32):
            src,arrays,_=source(vs,dict(format=family,bits=bits,width=96,height=80,frames=5),17)
            for interlaced in (False,True):
                l,t,r,b=4,4,8,4
                cropped=c.std.Crop(src,left=l,top=t,right=r,bottom=b)
                def pack(n,f):
                    out=f.copy()
                    for p in range(out.format.num_planes):
                        a=np.asarray(f[p]); rows=list(range(0,a.shape[0],2))+list(range(a.shape[0]-1,0,-2))
                        np.copyto(np.asarray(out[p]),a[rows] if interlaced else a)
                    return out
                packed=c.std.ModifyFrame(cropped,clips=cropped,selector=pack)
                for bt in (-1,1,2,3,4,5):
                    for model in ({},{'sigma2':5},{'pfactor':1,'px':2,'py':2,'pframe':4},{'sigma2':5,'pshow':True}):
                        kw=dict(bt=bt,bw=8,bh=8,ow=4,oh=4,opt=1,**model)
                        n = (0,2,4)[count % 3]
                        actual=c.neo_fft.FFT3D(src,l=l,t=t,r=r,b=b,interlaced=interlaced,**kw).get_frame(n)
                        expected=c.neo_fft.FFT3D(packed,**kw).get_frame(n)
                        assert dict(actual.props)==dict(src.get_frame(n).props)
                        for p in range(actual.format.num_planes):
                            sw=1 if family=='420' and p else 0
                            sh=sw; lp,tp=l>>sw,t>>sh
                            e=np.asarray(expected[p]); h,w=e.shape
                            rows=list(range(0,h,2))+list(range(h-1,0,-2)) if interlaced else list(range(h))
                            reference=arrays[n][p].copy(); reference[tp+np.array(rows),lp:lp+w]=e
                            assert reference.tobytes()==np.asarray(actual[p]).tobytes(),(family,bits,bt,model,p)
                        count+=1
    # Odd absolute top is allowed; only selected chroma requires aligned margins.
    yuv=c.std.BlankClip(width=96,height=80,format=vs.YUV420P8)
    odd,data,_=source(vs,dict(format='420',bits=8,width=96,height=80,frames=1),29)
    result=c.neo_fft.FFT3D(odd,bt=1,planes=[0],l=1,t=1,r=1,b=1,interlaced=True,bw=8,bh=8,opt=1).get_frame(0)
    gray=c.std.BlankClip(width=94,height=78,format=vs.GRAY8,length=1)
    rows=list(range(0,78,2))+list(range(77,0,-2))
    def odd_pack(n,f):
        out=f.copy(); np.copyto(np.asarray(out[0]),data[0][0][1:-1,1:-1][rows]); return out
    gray=c.std.ModifyFrame(gray,clips=gray,selector=odd_pack)
    filtered=c.neo_fft.FFT3D(gray,bt=1,bw=8,bh=8,opt=1).get_frame(0)
    expected=data[0][0].copy(); expected[1+np.array(rows),1:-1]=np.asarray(filtered[0])
    assert expected.tobytes()==np.asarray(result[0]).tobytes()
    for p in (1,2): assert data[0][p].tobytes()==np.asarray(result[p]).tobytes()
    for kwargs in ({'l':1},{'t':1},{'interlaced':True,'b':2},{'l':2147483647,'r':2147483647}):
        try: c.neo_fft.FFT3D(yuv,bt=1,bw=8,bh=8,**kwargs)
        except vs.Error: pass
        else: raise AssertionError(kwargs)
    # Selected exterior NaNs must be copied without entering transforms.
    fl=c.std.BlankClip(width=96,height=80,format=vs.GRAYS)
    def exterior(n,f):
        out=f.copy(); a=np.asarray(out[0]); a.view(np.uint32)[:]=0x7fc12345; a[4:-4,4:-4]=.25; return out
    fl=c.std.ModifyFrame(fl,clips=fl,selector=exterior)
    result=c.neo_fft.FFT3D(fl,bt=1,l=4,t=4,r=4,b=4,interlaced=True,bw=8,bh=8).get_frame(0)
    assert np.asarray(result[0]).view(np.uint32)[0,0]==0x7fc12345
    print('Phase 4 ROI/field scalar oracle:',count,'cases; alignment, exterior NaNs and properties passed')
if __name__=='__main__': main()
