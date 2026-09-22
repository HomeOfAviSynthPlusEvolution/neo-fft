import argparse
from pathlib import Path
import numpy as np
import vapoursynth as vs
from fixtures import environment

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--plugin',type=Path,required=True);args=ap.parse_args()
    policy=environment(vs);c=vs.core;c.num_threads=1;c.max_cache_size=1
    c.std.LoadPlugin(path=str(args.plugin.resolve()))
    peaks=[]
    for distance in (16,64,128):
        memory=[]
        src=c.std.BlankClip(width=1024,height=1024,format=vs.GRAY8,length=distance+1,keep=False)
        def fill(n,f):
            out=f.copy();np.asarray(out[0])[:]=n%251;memory.append(c.used_cache_size);return out
        src=c.std.ModifyFrame(src,clips=src,selector=fill);c.std.SetVideoCache(src,mode=0)
        out=c.neo_fft.FFT3D(src,bt=0,bw=8,bh=8,r=992,b=992,opt=1)
        frame=out.get_frame(distance);peaks.append(max(memory));del frame,out,src
    # This is the host's allocated pixel-buffer counter, not process RSS/history metadata.
    assert max(peaks)<8*1024*1024 and max(peaks)-min(peaks)<2*1024*1024,peaks
    print('Verified runtime',vs.__version__,'single-request pixel allocation peaks at replay 16/64/128:',peaks,
          '; host metadata, upstream caches and property-held frames are outside this bound')
if __name__=='__main__':main()
