"""FFT3D phase-3 public matrix; limits fixed at 1 LSB / 4e-6 float."""
import argparse
import json
from pathlib import Path
import numpy as np
import vapoursynth as vs
from fixtures import environment, source
from admission import fails

def main():
 p=argparse.ArgumentParser();p.add_argument('--plugin',type=Path,required=True);p.add_argument('--reference',type=Path);p.add_argument('--report',type=Path);args=p.parse_args()
 policy=environment(vs);c=vs.core;c.num_threads=4
 c.std.LoadPlugin(path=str(args.plugin.resolve()));call=c.neo_fft.FFT3D
 ref=None
 if args.reference:
  c.std.LoadPlugin(path=str(args.reference.resolve()));ref=c.neo_fft3d.FFT3D
 base=c.std.BlankClip(width=64,height=48,format=vs.GRAY8,length=7,color=[64])
 for kw in (dict(smin=5,smax=4),dict(scutoff=0),dict(svr=-1),dict(hr=0),dict(ht=-1),dict(sigma3=-1)):
  fails(lambda:call(base,**kw))
 call(base,bt=-1,sigma=1e30,ht=1e30).get_frame(0)
 fails(lambda:call(base,bt=1,sigma=1e30),'finite')
 fails(lambda:call(base,bt=-1,ht=1e30,dehalo=1),'finite')
 def only(n,f):
  if n!=3: raise vs.Error('unexpected neighbor')
  return f
 guarded=c.std.ModifyFrame(base,clips=base,selector=only)
 call(guarded,bt=-1,sigma2=3).get_frame(3)
 records=[]
 for bits,family in ((8,'gray'),(10,'420'),(16,'444'),(32,'rgb'),(32,'420')):
  for W,H in ((8,8),(7,5)):
   for pattern in (1,2,3,5,6):
    src,inputs,_=source(vs,dict(format=family,bits=bits,width=64,height=48,frames=7,pattern=pattern),17)
    for bt in (-1,1,2,3,4,5):
     for vary in (False,True):
      for enhance in ({},dict(sharpen=.4),dict(dehalo=.2),dict(sharpen=.4,dehalo=.2)):
       kw=dict(bt=bt,bw=W,bh=H,ow=W//2,oh=H//2,sigma=2,degrid=1,**enhance)
       if vary: kw.update(sigma2=3,sigma3=4,sigma4=5)
       clips=[call(src,opt=opt,**kw) for opt in (1,0)]
       if ref: clips.append(ref(src,opt=0,fft_backend='pocketfft',**kw))
       maxerr=0
       for n in (0,3,6):
        frames=[clip.get_frame(n) for clip in clips]
        for plane in range(src.format.num_planes):
         a=np.asarray(frames[0][plane]).astype(float)
         for f in frames[1:]:
          error=float(np.max(np.abs(a-np.asarray(f[plane]).astype(float))))
          maxerr=max(maxerr,error)
          assert error <= (4e-6 if bits==32 else 1),(bits,family,W,H,pattern,kw,n,plane,error)
         if bt==-1 and not enhance:
          assert np.max(np.abs(a-(np.clip(inputs[n][plane].astype(float),0,1) if bits==32 else inputs[n][plane].astype(float)))) <= (4e-6 if bits==32 else 1)
        assert dict(frames[0].props)==dict(src.get_frame(n).props)
       records.append(dict(bits=bits,family=family,W=W,H=H,pattern=pattern,params=kw,max_error=maxerr))
 if args.report: args.report.write_text(json.dumps(records,indent=2))
 print('FFT3D models:',len(records),'cases passed')
if __name__=='__main__':main()
