"""Final cross-model reference matrix, original image tolerances retained."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import vapoursynth as vs
from fixtures import environment, source

def fft_composition_oracle(c):
 # Rectangular analysis/synthesis removes overlap geometry from this oracle;
 # independent FFT/formulas verify fused versus two-stage gains and temporal order.
 src,inputs,_=source(vs,dict(format='gray',bits=32,width=64,height=48,frames=7,pattern=6),61)
 for W,H in ((4,4),(7,5)):
  ky,kx=np.indices((H,W//2+1));norm=1/(W*H)
  fy=(H-2*np.abs(ky-H//2))/H;fx=kx/(W//2+1);f=np.sqrt((fx*fx+fy*fy)/2)
  a=np.sqrt(.5)/4;b=2*a
  sig=np.where(f<a,5+(4-5)*f/a,np.where(f<b,4+(3-4)*(f-a)/(b-a),2+(3-2)*(1-f)/(1-b)))/255
  analytic=sig*sig/norm
  r=(2*np.minimum(ky,H-ky)/H)**2+(2*kx/W)**2;weight=r/(r+.1**2)
  dy=np.where(ky<H//2,ky,H-ky);d2=dy**2/(H//2)**2+kx**2/(W//2)**2
  ws=1-np.exp(-d2/(2*.3**2));wh=np.exp(-.7*d2*4)-np.exp(-d2*4);wh/=wh.max()
  A=(4/255)**2/norm;B=(20/255)**2/norm;C=.1**2/norm
  for degrid in (0,1):
   def residual(z):
    mean=np.zeros_like(z)
    if degrid:mean[0,0]=z[0,0].real
    return z-mean,mean
   def gain(q):return (1+.4*ws*np.sqrt(q*B/((q+A)*(q+B))))*(q+C)/(q+C+.2*wh*q)
   sampled,_=residual(np.fft.rfft2(inputs[0][0][:H,:W]));sampled=np.abs(sampled)**2*weight
   for model in ('uniform','analytic','sample'):
    for T in (-1,1,2,3,4,5):
     kw=dict(bw=W,bh=H,ow=0,oh=0,bt=T,sigma=2,degrid=degrid,sharpen=.4,dehalo=.2,ht=.1)
     if model=='analytic':kw.update(sigma2=3,sigma3=4,sigma4=5)
     if model=='sample':kw.update(pfactor=1,pframe=0,px=1,py=1)
     power=analytic if model=='analytic' else sampled if model=='sample' else (2/255)**2/norm
     spatial=np.fft.rfft2(inputs[3][0][:H,:W])
     if T<=1:
      R,M=residual(spatial);q=np.abs(R)**2+1e-15
      if T==1:
       wiener=np.maximum((q-power)/q,0)
       if model=='uniform':result=R*wiener*gain(q)+M
       else:
        R,M=residual(R*wiener+M);result=R*gain(np.abs(R)**2+1e-15)+M
      else:result=R*gain(q)+M
     else:
      spectra=np.stack([np.fft.rfft2(inputs[n][0][:H,:W]) for n in range(3-T//2,3+(T-1)//2+1)])
      temporal=np.fft.fft(spectra,axis=0);removed=np.zeros_like(temporal)
      if degrid:removed[0,0,0]=T*spatial[0,0].real
      R=temporal-removed;q=np.abs(R)**2+1e-15
      spatial=np.fft.ifft(R*np.maximum((q-T*power)/q,0)+removed,axis=0)[T//2]
      R,M=residual(spatial);result=R*gain(np.abs(R)**2+1e-15)+M
     want=np.clip(np.fft.irfft2(result,s=(H,W)),0,1)
     for opt in (0,1):
      actual=np.asarray(c.neo_fft.FFT3D(src,opt=opt,**kw).get_frame(3)[0])[:H,:W]
      delta=float(np.max(np.abs(actual-want)))
      assert delta<=4e-6,('composition',W,H,degrid,model,T,opt,delta)

def main():
 p=argparse.ArgumentParser();p.add_argument('--plugin',type=Path,required=True);p.add_argument('--fft-reference',type=Path);p.add_argument('--dft-reference',type=Path);p.add_argument('--report',type=Path);a=p.parse_args()
 environment(vs);c=vs.core;c.num_threads=4;c.std.LoadPlugin(path=str(a.plugin.resolve()))
 if a.fft_reference:c.std.LoadPlugin(path=str(a.fft_reference.resolve()))
 if a.dft_reference:c.std.LoadPlugin(path=str(a.dft_reference.resolve()))
 fft_composition_oracle(c)
 records=[]
 def compare(nodes,src,inputs,order,selected,case):
  error=0;limit=4e-6 if src.format.bits_per_sample==32 else 1
  for n in order:
   frames=[node.get_frame(n) for node in nodes]
   for plane in range(src.format.num_planes):
    baseline=np.asarray(frames[0][plane]).astype(float)
    for frame in frames[1:]:
     delta=float(np.max(np.abs(baseline-np.asarray(frame[plane]).astype(float))))
     assert delta<=limit,(case,n,plane,delta)
     error=max(error,delta)
    if plane not in selected:
     for frame in frames:assert np.asarray(frame[plane]).tobytes()==inputs[n][plane].tobytes()
   for frame in frames:assert dict(frame.props)==dict(src.get_frame(n).props)
  records.append(dict(case=case,max_error=error))
 for bits,family in ((8,'gray'),(10,'420'),(16,'444'),(32,'rgb'),(32,'420')):
  for W,H in ((8,8),(7,5)):
   for pattern in range(7):
    src,inputs,_=source(vs,dict(format=family,bits=bits,width=64,height=48,frames=7,pattern=pattern),37)
    for bt in (1,2,3,4,5):
     for degrid in (0,1):
      # Reference equality is specified at pfactor=1 only.
      kw=dict(bw=W,bh=H,ow=W//2,oh=H//2,bt=bt,pfactor=1,pframe=2,px=1,py=1,degrid=degrid,sharpen=.3,dehalo=.2,fft_backend='pocketfft')
      nodes=[c.neo_fft.FFT3D(src,opt=o,**kw) for o in (1,0)]
      if a.fft_reference:nodes.append(c.neo_fft3d.FFT3D(src,opt=0,**kw))
      compare(nodes,src,inputs,(3,0,6),range(src.format.num_planes),('FFT3D',bits,family,W,H,pattern,bt,degrid))
  # Curves on chroma as well as luma, both mean settings and plane selections.
  for T,S in ((3,5),(5,8),(15,3)):
   for mean in (False,True):
    src,inputs,_=source(vs,dict(format=family,bits=bits,width=32,height=24,frames=T+2,pattern=6),41)
    for ftype in range(5):
     for system in (0,1):
      selected=[src.format.num_planes-1] if system else list(range(src.format.num_planes))
      lo,hi=(2,12) if ftype<2 else (.25,.8)
      kw=dict(tbsize=T,sbsize=S,sosize=S//2,ftype=ftype,slocation=[0,lo,1,hi],ssystem=system,zmean=mean,planes=selected,fft_backend='pocketfft')
      nodes=[c.neo_fft.DFTTest(src,opt=o,**kw) for o in (1,0)]
      if a.dft_reference:nodes.append(c.neo_dfttest.DFTTest(src,opt=1,**kw))
      compare(nodes,src,inputs,(T//2,0,T+1),selected,('DFTTest',bits,family,T,S,mean,ftype,system))
 # Fresh reference preview instances avoid legacy mutable-window history.
 for bits in (16,32):
  src,inputs,_=source(vs,dict(format='rgb',bits=bits,width=64,height=48,frames=7,pattern=6),53)
  for auto in (False,True):
   for n in (0,3,6):
    kw=dict(bw=8,bh=8,ow=4,oh=4,bt=5,pfactor=1,pshow=True,px=0 if auto else 1,py=0 if auto else 1,fft_backend='pocketfft')
    nodes=[c.neo_fft.FFT3D(src,opt=o,**kw) for o in (1,0)]
    if a.fft_reference:nodes.append(c.neo_fft3d.FFT3D(src,opt=0,**kw))
    compare(nodes,src,inputs,(n,),range(3),('preview',bits,auto,n))
 if a.report:
  def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest() if path else None
  a.report.write_text(json.dumps(dict(candidate=digest(a.plugin),fft_reference=digest(a.fft_reference),dft_reference=digest(a.dft_reference),cases=records),indent=2))
 print('Phase 3 final matrix:',len(records),'cases passed')
if __name__=='__main__':main()
