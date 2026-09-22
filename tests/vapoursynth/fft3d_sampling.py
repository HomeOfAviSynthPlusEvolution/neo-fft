"""Sample publication/order, preview and declared dependency regressions."""
import argparse
from pathlib import Path
import numpy as np
import vapoursynth as vs
from fixtures import environment, source
from admission import fails

def main():
 p=argparse.ArgumentParser();p.add_argument('--plugin',type=Path,required=True);a=p.parse_args()
 environment(vs);c=vs.core;c.num_threads=4;c.std.LoadPlugin(path=str(a.plugin.resolve()));call=c.neo_fft.FFT3D
 count=0
 for bits,family in ((8,'gray'),(10,'420'),(16,'444'),(32,'rgb'),(32,'420')):
  src,inputs,_=source(vs,dict(format=family,bits=bits,width=64,height=48,frames=7,pattern=6),19)
  for W,H in ((8,8),(7,5)):
   for bt in (1,2,3,4,5):
    for auto in (False,True):
     kw=dict(bw=W,bh=H,ow=W//2,oh=H//2,bt=bt,pfactor=.7,pframe=6,px=0 if auto else 1,py=0 if auto else 1)
     nodes=[call(src,opt=o,**kw) for o in (0,1)]
     # Boundary-first, then interior; independent instances requested in reverse.
     expected={n:[np.array(x) for x in nodes[0].get_frame(n)] for n in (0,3,6)}
     for n in (6,3,0):
      frame=nodes[1].get_frame(n)
      for plane,x in enumerate(frame):
       assert np.max(np.abs(np.asarray(x).astype(float)-expected[n][plane].astype(float))) <= (4e-6 if bits==32 else 1)
     fresh=call(src,**kw)
     futures=[fresh.get_frame_async(n) for n in (3,0,6)]
     for n,f in zip((3,0,6),futures):
      frame=f.result()
      for plane,x in enumerate(frame): assert np.array_equal(x,expected[n][plane])
     count+=1
   # Rectangular preview oracle, manual grid 1,1 starts at source 0,0.
   for analytic in (False,True):
    kw=dict(bw=W,bh=H,ow=0,oh=0,px=1,py=1,pshow=True,bt=5,sharpen=1e30,ht=1e30,dehalo=1e30,pframe=6)
    kw.update(dict(sigma2=3) if analytic else dict(pfactor=1,sigma=1e30))
    def only(n,f):
     if n!=3: raise vs.Error('preview requested a different frame')
     return f
    guarded=c.std.ModifyFrame(src,clips=src,selector=only)
    frame=call(guarded,**kw).get_frame(3)
    for plane,x in enumerate(frame):
     midpoint=(1<<(bits-1)) if bits!=32 and family in ('420','444') and plane else 0
     want=np.full_like(inputs[3][plane],midpoint)
     want[:H,:W]=inputs[3][plane][:H,:W]
     if bits==32: want=np.clip(want,0,1)
     assert np.max(np.abs(np.asarray(x).astype(float)-want.astype(float))) <= (4e-6 if bits==32 else 1)
    assert dict(frame.props)==dict(src.get_frame(3).props)
    count+=1
 base=c.std.BlankClip(width=64,height=48,length=7,format=vs.GRAY8,color=[64])
 for kw in (dict(px=-1),dict(py=-1),dict(pcutoff=0),dict(pfactor=-1),dict(pfactor=1,px=100)):
  fails(lambda:call(base,**kw))
 # pframe is clamped, and bt=-1 discards dormant invalid grid coordinates.
 for pf in (-100,100): call(base,pfactor=1,pframe=pf,bw=8,bh=8).get_frame(3)
 call(base,bt=-1,pfactor=1,px=100,py=100).get_frame(3)
 print('FFT3D sampling:',count,'cases passed')
if __name__=='__main__':main()
