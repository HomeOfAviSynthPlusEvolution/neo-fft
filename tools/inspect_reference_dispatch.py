"""Read Highway dispatch state in an unmodified, symbol-bearing GCC x64 reference DLL.

This is a version-specific diagnostic for Highway 1.3/1.4, not a plugin ABI.
Function RVAs come from this exact PE's COFF symbols. No target is forced.
"""
import argparse
import ctypes
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tests/vapoursynth'))
from fixtures import environment
import vapoursynth as vs

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--plugin', type=Path, required=True)
    p.add_argument('--nm', type=Path, required=True)
    p.add_argument('--algorithm', choices=('FFT3D', 'DFTTest'), required=True)
    p.add_argument('--output', type=Path, required=True)
    args=p.parse_args()
    data=args.plugin.read_bytes()
    pe=struct.unpack_from('<I',data,0x3c)[0]
    if data[pe:pe+4]!=b'PE\0\0' or struct.unpack_from('<H',data,pe+4)[0]!=0x8664:
        raise RuntimeError('diagnostic requires a Windows x64 PE image')
    optional=pe+24
    preferred=struct.unpack_from('<Q',data,optional+24)[0]
    image_size=struct.unpack_from('<I',data,optional+56)[0]
    symbols={}
    for line in subprocess.check_output([str(args.nm),str(args.plugin)],text=True).splitlines():
        fields=line.split()
        if len(fields)==3:
            symbols[fields[2]]=(int(fields[0],16),fields[1])
    policy=environment(vs)
    core=vs.core; core.num_threads=4
    before={p.identifier for p in core.plugins()}
    core.std.LoadPlugin(path=str(args.plugin.resolve()))
    plugins=[p for p in core.plugins() if p.identifier not in before]
    if len(plugins)!=1: raise RuntimeError('expected one explicitly loaded reference')
    plugin=plugins[0]
    names=[f.name for f in plugin.functions() if args.algorithm.lower() in f.name.lower()]
    if len(names)!=1: raise RuntimeError('ambiguous registration')
    library=ctypes.WinDLL(str(args.plugin.resolve()))
    base=library._handle
    def address(symbol):
        va,kind=symbols[symbol]; rva=va-preferred
        if not 0<=rva<image_size: raise RuntimeError('symbol is outside PE image')
        return base+rva
    def contains(pointer,size=8):
        if not base<=pointer<=base+image_size-size: raise RuntimeError('pointer outside loaded image')
        return pointer
    chosen=ctypes.CFUNCTYPE(ctypes.c_void_p)(address('_ZN3hwy15GetChosenTargetEv'))
    supported=ctypes.CFUNCTYPE(ctypes.c_int64)(address('_ZN3hwy16SupportedTargetsEv'))()
    def mask(): return ctypes.c_uint64.from_address(contains(chosen())).value
    initial=mask()
    src=core.std.BlankClip(width=128,height=96,format=vs.GRAYS,color=[.25])
    params=dict(fft_backend='pocketfft')
    params.update(dict(bt=1,ncpu=1,mt=False) if args.algorithm=='FFT3D' else dict(tbsize=1,threads=1,fft_threads=1))
    call=getattr(plugin,names[0])
    call(src,opt=1,**params).get_frame(0)
    scalar_mask=mask()
    call(src,opt=0,**params).get_frame(0)
    selected_mask=mask()
    if selected_mask==1: raise RuntimeError('public opt=0 did not initialize Highway dispatch')
    rows=[]
    for name,(va,kind) in symbols.items():
        if kind.lower()!='r' or not name.endswith('HighwayDispatchTableE'): continue
        if not ('Apply2D_' in name or 'ProcessSpatial_' in name): continue
        entries=(ctypes.c_uint64*17).from_address(contains(address(name),17*8))
        available=sum(1<<i for i,v in enumerate(entries) if v)
        bits=selected_mask & available
        if not bits: raise RuntimeError('no supported entry in spatial dispatch table')
        index=(bits & -bits).bit_length()-1
        selected_va=entries[index]-base+preferred
        matches=[s for s,(a,k) in symbols.items() if a==selected_va and k.lower()=='t' and not s.startswith('.')]
        if not matches: raise RuntimeError('selected function lacks a COFF symbol')
        rows.append(dict(table=name,index=index,function=matches[0]))
    if not rows: raise RuntimeError('no known spatial dispatch tables')
    targets={3:'AVX10_2',4:'AVX3_SPR',6:'AVX3_ZEN4',7:'AVX3_DL',8:'AVX3',9:'AVX2',11:'SSE4',12:'SSSE3',14:'SSE2'}
    report=dict(plugin=str(args.plugin.resolve()),sha256=hashlib.sha256(data).hexdigest(),
                method='read-only COFF symbol/RVA inspection of Highway 1.3/1.4 x64 dispatch after public VS requests',
                supported_mask=hex(supported),supported_targets=[n for b,n in targets.items() if supported&(1<<b)],
                initial_mask=initial,after_scalar_mask=scalar_mask,after_highway_mask=selected_mask,spatial_tables=rows)
    args.output.write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))

if __name__=='__main__': main()
