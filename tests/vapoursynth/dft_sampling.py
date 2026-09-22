"""DFT sampling public oracle and black-box matrix. Frozen 1 LSB / 4e-6."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import vapoursynth as vs
from fixtures import environment, source
from admission import fails

def main():
 p=argparse.ArgumentParser();p.add_argument('--plugin',type=Path,required=True);p.add_argument('--reference',type=Path);p.add_argument('--report',type=Path);a=p.parse_args()
 environment(vs);c=vs.core;c.num_threads=4;c.std.LoadPlugin(path=str(a.plugin.resolve()));call=c.neo_fft.DFTTest
 ref=None
 if a.reference: c.std.LoadPlugin(path=str(a.reference.resolve()));ref=c.neo_dfttest.DFTTest
 base=c.std.BlankClip(width=32,height=24,length=15,format=vs.YUV420P8)
 for loc in ([0],[-1,0,0,0],[13,0,0,0],[0,3,0,0],[0,1,12,0],[0,1,0,14],[0,0,0,0]*501):
  for extra in ({},dict(planes=[]),dict(ftype=2)):
   fails(lambda:call(base,sbsize=3,sosize=1,nlocation=loc,**extra))
 for alpha in (0,-1,float('nan'),float('inf')): fails(lambda:call(base,alpha=alpha,planes=[]))
 # A declared but output-unselected sample patch is read; the rest is not.
 raw=c.std.BlankClip(width=32,height=24,length=9,format=vs.YUV444PS,color=[.2,.1,.1])
 def damaged(n,f,inside):
  out=f.copy();np.asarray(out[1])[1 if inside else 20,1 if inside else 20]=np.nan;return out
 for inside in (False,True):
  clip=c.std.ModifyFrame(raw,clips=raw,selector=lambda n,f:damaged(n,f,inside))
  out=call(clip,tbsize=3,sbsize=3,sosize=1,nlocation=[0,1,0,0],planes=[0])
  if inside: fails(lambda:out.get_frame(5),'finite')
  else: out.get_frame(5)
 # Inactive samples/empty planes: n is the sole dependency with T=1.
 def only(n,f):
  if n!=5: raise vs.Error('inactive sample requested')
  return f
 guarded=c.std.ModifyFrame(raw,clips=raw,selector=only)
 for extra in (dict(planes=[]),dict(ftype=2,sigma=1)):
  call(guarded,tbsize=1,sbsize=3,sosize=1,nlocation=[0,1,0,0],**extra).get_frame(5)
 # Active sample replaces otherwise overflowing derived primary profile.
 call(raw,tbsize=3,sbsize=3,sosize=1,nlocation=[0,1,0,0],slocation=[0,3e38,1,3e38]).get_frame(5)
 # Independent complete central-pixel oracle with rectangular 3x3x3 windows.
 src,inputs,_=source(vs,dict(format='444',bits=32,width=24,height=20,frames=7,pattern=6),47)
 locations=[(0,1,2,3),(3,2,4,5)]
 for mean in (False,True):
  for alpha in (.5,2):
   powers=[]
   for fn,plane,y,x in locations:
    block=np.stack([inputs[n][plane][y:y+3,x:x+3] for n in range(fn,fn+3)]).astype(float)*255/np.sqrt(27)
    spectrum=np.fft.rfftn(block)
    if mean:spectrum[0,0,0]=0
    powers.append(np.abs(spectrum)**2)
   table=np.mean(powers,axis=0)*alpha
   for ftype in (0,1):
    kw=dict(tbsize=3,sbsize=3,smode=0,swin=7,twin=7,zmean=mean,alpha=alpha,ftype=ftype,nlocation=sum((list(t) for t in locations),[]))
    result=call(src,**kw).get_frame(3)
    for plane in range(3):
     for y,x in ((5,5),(10,12),(15,20)):
      block=np.stack([inputs[n][plane][y-1:y+2,x-1:x+2] for n in (2,3,4)]).astype(float)*255/np.sqrt(27)
      spectrum=np.fft.rfftn(block);dc=spectrum[0,0,0]
      if mean:spectrum[0,0,0]=0
      q=np.abs(spectrum)**2
      gain=np.maximum((q-table)/(q+1e-15),0) if ftype==0 else (q>=table)
      spectrum*=gain
      if mean:spectrum[0,0,0]=dc
      expected=np.fft.irfftn(spectrum,s=(3,3,3),axes=(0,1,2))[1,1,1]*np.sqrt(27)/255
      assert abs(np.asarray(result[plane])[y,x]-expected)<=4e-6,(mean,alpha,ftype,plane,expected,np.asarray(result[plane])[y,x])
 records=[]
 for bits,family in ((8,'gray'),(10,'420'),(16,'444'),(32,'rgb'),(32,'420')):
  for T,S in ((1,8),(3,5),(5,8),(15,3)):
   for pattern in range(7):
    src,inputs,_=source(vs,dict(format=family,bits=bits,width=32,height=24,frames=max(7,T+2),pattern=pattern),23)
    locations=[0,0,0,0,2,src.format.num_planes-1,1,1]
    for mean in (False,True):
     for ftype in (0,1):
      kw=dict(tbsize=T,sbsize=S,sosize=S//2,ftype=ftype,zmean=mean,nlocation=locations,alpha=2,fft_backend='pocketfft')
      nodes=[call(src,opt=o,**kw) for o in (1,0)]
      if ref:nodes.append(ref(src,opt=1,**kw))
      err=0
      for n in (0,src.num_frames//2,src.num_frames-1):
       frames=[node.get_frame(n) for node in nodes]
       for plane in range(src.format.num_planes):
        baseline=np.asarray(frames[0][plane]).astype(float)
        for frame in frames[1:]:
         delta=float(np.max(np.abs(baseline-np.asarray(frame[plane]).astype(float))))
         assert np.isfinite(delta),(bits,family,T,S,pattern,mean,ftype,n,plane,delta)
         err=max(err,delta)
       assert err<=(4e-6 if bits==32 else 1),(bits,family,T,S,pattern,mean,ftype,n,err)
       assert dict(frames[0].props)==dict(src.get_frame(n).props)
      # Repeated whole lists preserve averaging, and simultaneous first calls.
      dup=call(src,**dict(kw,nlocation=locations*2))
      for n,future in [(n,dup.get_frame_async(n)) for n in (src.num_frames-1,0,src.num_frames//2)]:
       frame=future.result();expected=nodes[0].get_frame(n)
       for plane in range(src.format.num_planes):assert np.max(np.abs(np.asarray(frame[plane]).astype(float)-np.asarray(expected[plane]).astype(float)))<=(4e-6 if bits==32 else 1)
      records.append(dict(bits=bits,family=family,T=T,S=S,pattern=pattern,zmean=mean,ftype=ftype,max_error=err))
 if a.report:a.report.write_text(json.dumps(dict(candidate_sha256=hashlib.sha256(a.plugin.read_bytes()).hexdigest(),reference_sha256=hashlib.sha256(a.reference.read_bytes()).hexdigest() if a.reference else None,cases=records),indent=2))
 print('DFT sampling:',len(records),'cases and independent oracle passed')
if __name__=='__main__':main()
